//! Optional post-search scoring. Search/CFD/representative policy stays upstream.
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write};
use std::path::PathBuf;
use std::process::Command;

use clap::{Args as ClapArgs, ValueEnum};

#[derive(Clone, Copy, Debug, PartialEq, ValueEnum)]
pub(crate) enum ScoringMode {
    Cfd,
    CrisprBulge,
    Both,
}
impl ScoringMode {
    fn name(self) -> &'static str {
        match self {
            Self::Cfd => "cfd",
            Self::CrisprBulge => "crispr-bulge",
            Self::Both => "both",
        }
    }
}

#[derive(ClapArgs, Debug)]
pub(crate) struct ScoringArgs {
    /// Scoring mode; neural modes are experimental and retain unsupported hits
    #[arg(long, value_enum, default_value = "cfd")]
    pub scoring: ScoringMode,
    /// Score an existing unannotated PAF instead of running any search
    #[arg(long)]
    pub score_paf: Option<PathBuf>,
    /// NEW exclusive result directory, required for neural modes
    #[arg(long)]
    pub scoring_output_dir: Option<PathBuf>,
    /// Absolute dedicated Python interpreter with optional TensorFlow dependencies
    #[arg(long)]
    pub crispr_bulge_python: Option<PathBuf>,
    /// Existing checksum-pinned CRISPR-Bulge model checkout (no runtime downloads)
    #[arg(long)]
    pub crispr_bulge_model_root: Option<PathBuf>,
    /// Scoring script in the matching source bundle; override when relocating
    #[arg(long, default_value = concat!(env!("CARGO_MANIFEST_DIR"), "/research/crispr_bulge/score_paf.py"))]
    pub scoring_script: PathBuf,
    /// Optional trusted reference SHA256 from the original search
    #[arg(long)]
    pub scoring_reference_sha256: Option<String>,
    /// Optional prepared trusted FASTA offset index (requires both SHA256 pins)
    #[arg(long)]
    pub scoring_reference_index: Option<PathBuf>,
    /// Trusted SHA256 of the prepared FASTA offset index
    #[arg(long)]
    pub scoring_reference_index_sha256: Option<String>,
    /// Number of eligible hits per neural inference chunk (models load once)
    #[arg(long, default_value_t = 512, value_parser = clap::value_parser!(u32).range(1..=4096))]
    pub scoring_batch_size: u32,
    /// Neural worker timeout in seconds, excluding reference validation
    #[arg(long, default_value_t = 300, value_parser = clap::value_parser!(u32).range(1..))]
    pub scoring_worker_timeout: u32,
}

impl ScoringArgs {
    pub fn validate(&self) -> Result<(), String> {
        if self.scoring == ScoringMode::Cfd {
            if self.scoring_output_dir.is_some()
                || self.crispr_bulge_python.is_some()
                || self.crispr_bulge_model_root.is_some()
                || self.scoring_reference_index.is_some()
                || self.scoring_reference_index_sha256.is_some()
                || self.scoring_reference_sha256.is_some()
            {
                return Err("neural scoring options require --scoring crispr-bulge or both".into());
            }
            return Ok(());
        }
        let directory = self
            .scoring_output_dir
            .as_ref()
            .ok_or("neural scoring requires --scoring-output-dir")?;
        if directory.exists() {
            return Err("--scoring-output-dir must not exist".into());
        }
        let python = self
            .crispr_bulge_python
            .as_ref()
            .ok_or("neural scoring requires --crispr-bulge-python")?;
        if !python.is_absolute() || !python.is_file() {
            return Err(
                "--crispr-bulge-python must be an existing absolute interpreter path".into(),
            );
        }
        if !self
            .crispr_bulge_model_root
            .as_ref()
            .is_some_and(|p| p.is_dir())
            || !self.scoring_script.is_file()
        {
            return Err(
                "neural scoring requires an existing model root and matching scoring script bundle"
                    .into(),
            );
        }
        if self.scoring_reference_index.is_some() != self.scoring_reference_index_sha256.is_some()
            || (self.scoring_reference_index.is_some() && self.scoring_reference_sha256.is_none())
        {
            return Err(
                "prepared reference index requires index path, index SHA256 and reference SHA256"
                    .into(),
            );
        }
        for sha in [
            &self.scoring_reference_sha256,
            &self.scoring_reference_index_sha256,
        ]
        .into_iter()
        .flatten()
        {
            if sha.len() != 64
                || !sha
                    .bytes()
                    .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
            {
                return Err("scoring SHA256 pins must be 64 lowercase hex characters".into());
            }
        }
        Ok(())
    }
}

pub(crate) struct Output {
    writer: Box<dyn Write>,
    source: Option<PathBuf>,
}
impl Output {
    pub fn new(args: &crate::Args) -> Result<Self, String> {
        let a = &args.scoring;
        if a.scoring == ScoringMode::Cfd {
            return Ok(Self {
                writer: Box::new(io::stdout()),
                source: a.score_paf.clone(),
            });
        }
        let directory = a.scoring_output_dir.as_ref().unwrap();
        fs::create_dir(directory).map_err(|e| format!("cannot create scoring directory: {e}"))?;
        fs::write(directory.join("status"), "INCOMPLETE\n").map_err(|e| e.to_string())?;
        if let Some(source) = &a.score_paf {
            return Ok(Self {
                writer: Box::new(io::sink()),
                source: Some(source.clone()),
            });
        }
        let source = directory.join("input.paf");
        let file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&source)
            .map_err(|e| e.to_string())?;
        Ok(Self {
            writer: Box::new(file),
            source: Some(source),
        })
    }

    pub fn finish(mut self, args: &crate::Args) -> Result<(), String> {
        self.writer.flush().map_err(|e| e.to_string())?;
        let a = &args.scoring;
        if a.scoring == ScoringMode::Cfd {
            if let Some(source) = &self.source {
                io::copy(
                    &mut File::open(source).map_err(|e| e.to_string())?,
                    &mut io::stdout(),
                )
                .map_err(|e| e.to_string())?;
            }
            return Ok(());
        }
        let root = a.scoring_output_dir.as_ref().unwrap();
        let result = root.join("scoring");
        let mut command = Command::new(a.crispr_bulge_python.as_ref().unwrap());
        command
            .args(["-I", "-B"])
            .arg(&a.scoring_script)
            .arg("--input")
            .arg(self.source.as_ref().unwrap())
            .arg("--out-dir")
            .arg(&result)
            .arg("--scoring")
            .arg(a.scoring.name())
            .arg("--annotations")
            .arg("--preserve-cfd")
            .arg("--reference")
            .arg(&args.reference)
            .arg("--model-root")
            .arg(a.crispr_bulge_model_root.as_ref().unwrap())
            .arg("--model-python")
            .arg(a.crispr_bulge_python.as_ref().unwrap())
            .arg("--batch-size")
            .arg(a.scoring_batch_size.to_string())
            .arg("--worker-timeout")
            .arg(a.scoring_worker_timeout.to_string());
        if let Some(guides) = &args.guides {
            command.arg("--guides").arg(guides);
        } else {
            command.arg("--guide").arg(args.guide.as_ref().unwrap());
        }
        if let Some(sha) = &a.scoring_reference_sha256 {
            command.arg("--expected-reference-sha256").arg(sha);
        }
        if let Some(index) = &a.scoring_reference_index {
            command
                .arg("--reference-index")
                .arg(index)
                .arg("--reference-index-sha256")
                .arg(a.scoring_reference_index_sha256.as_ref().unwrap());
        }
        // Child writes no PAF to stdout. Publish only a complete result bundle.
        command.stdout(std::process::Stdio::null());
        let status = command
            .status()
            .map_err(|e| format!("scoring process failed to start: {e}"))?;
        if !status.success()
            || !fs::read_to_string(result.join("status"))
                .unwrap_or_default()
                .starts_with("COMPLETE ")
            || !result.join("annotations.tsv").is_file()
            || !result.join("output.paf").is_file()
        {
            return Err(format!(
                "optional scoring failed; no PAF published; inspect {}",
                root.display()
            ));
        }
        fs::write(
            root.join("status"),
            "COMPLETE scoring; original representatives retained\n",
        )
        .map_err(|e| e.to_string())?;
        io::copy(
            &mut File::open(result.join("output.paf")).map_err(|e| e.to_string())?,
            &mut io::stdout(),
        )
        .map_err(|e| e.to_string())?;
        Ok(())
    }
}
impl Write for Output {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.writer.write(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}
