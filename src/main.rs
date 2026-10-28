use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::sync::Arc;

use bio::io::fasta;
use clap::{ArgGroup, Parser};
use flate2::read::MultiGzDecoder;
use lib_wfa2::affine_wavefront::AffineWavefronts;
use rayon::prelude::*;

mod cfd_score;
mod columba;
mod hit;
mod reporting;
mod verification;

use columba::{
    candidate_edit_distance_bound, load_guide_fasta, parse_columba_sam_file,
    run_columba_candidate_generation, ColumbaQuery, ColumbaRunConfig,
};
use hit::Hit;
use reporting::{report_filtered_hits, report_filtered_hits_for_query};
use verification::{
    build_verified_hit, load_reference_records, scan_window, uppercase_ascii_sequence,
    verify_columba_candidates,
};

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'G' => b'C',
            b'C' => b'G',
            b'N' => b'N',
            _ => b'N',
        })
        .collect()
}

#[derive(Parser)]
#[command(
    author,
    version,
    about = "CRISPR guide RNA off-target scanner",
    group(ArgGroup::new("guide_input").required(true).args(["guide", "guides"]))
)]
pub(crate) struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long)]
    guide: Option<String>,

    /// Multi-guide FASTA with unique, stable record IDs (-G)
    #[arg(short = 'G', long)]
    guides: Option<PathBuf>,

    /// Required final two bases of a complete three-base PAM (GG requires NGG)
    #[arg(short = 'p', long, default_value = "GG")]
    pam: String,

    /// Keep noncanonical or unavailable PAMs; score using observed PAM bases
    #[arg(long)]
    allow_noncanonical_pam: bool,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Minimum fraction of guide that must match (0.0-1.0)
    #[arg(short = 'f', long, default_value = "0.75")]
    min_match_fraction: f32,

    /// Path to mismatch scores file for CFD calculation
    #[arg(long, default_value = "mismatch_scores.txt")]
    mismatch_scores: PathBuf,

    /// Path to PAM scores file for CFD calculation
    #[arg(long, default_value = "pam_scores.txt")]
    pam_scores: PathBuf,

    /// Size of sequence window to scan (bp, default: 4x guide length)
    #[arg(short = 'w', long)]
    window_size: Option<usize>,

    /// Number of threads to use (default: number of logical CPUs)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Bypass sequence/PAM eligibility; overlap representative selection still applies
    #[arg(long)]
    no_filter: bool,

    /// Existing Columba SAM file to parse for candidate locations
    #[arg(long)]
    columba_sam: Option<PathBuf>,

    /// Path to Columba executable
    #[arg(long)]
    columba_bin: Option<PathBuf>,

    /// Existing Columba FM-index prefix
    #[arg(long)]
    columba_index: Option<PathBuf>,

    /// Preserve temporary Columba SAM output
    #[arg(long)]
    keep_columba_sam: bool,
}

fn open_reference_reader(path: &PathBuf) -> fasta::Reader<BufReader<Box<dyn BufRead>>> {
    let file = File::open(path).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if path.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    fasta::Reader::new(reader)
}

fn main() {
    let mut args = Args::parse();
    args.pam.make_ascii_uppercase();
    if args.pam.len() != 2
        || !args
            .pam
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
    {
        eprintln!("--pam must contain exactly two DNA bases (A/C/G/T), e.g. GG for NGG");
        std::process::exit(2);
    }
    let allow_noncanonical_pam = args.allow_noncanonical_pam || args.no_filter;

    if let Err(e) = cfd_score::init_score_matrices(
        args.mismatch_scores
            .to_str()
            .unwrap_or("mismatch_scores.txt"),
        args.pam_scores.to_str().unwrap_or("pam_scores.txt"),
    ) {
        eprintln!("Warning: CFD scoring disabled - {}", e);
    }

    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    let batch_mode = args.guides.is_some();
    let queries = if let Some(path) = &args.guides {
        load_guide_fasta(path).unwrap_or_else(|e| {
            eprintln!("{}", e);
            std::process::exit(1);
        })
    } else {
        vec![ColumbaQuery {
            query_name: "guide".to_string(),
            sequence: args.guide.clone().expect("clap requires --guide"),
        }]
    };

    if batch_mode
        && args.columba_sam.is_none()
        && args.columba_bin.is_none()
        && args.columba_index.is_none()
    {
        eprintln!("--guides requires --columba-bin and --columba-index (or --columba-sam)");
        std::process::exit(1);
    }

    if args.columba_sam.is_some() || args.columba_bin.is_some() || args.columba_index.is_some() {
        let reference_records = load_reference_records(&args.reference).unwrap_or_else(|e| {
            eprintln!("{}", e);
            std::process::exit(1);
        });

        let candidates = if let Some(path) = &args.columba_sam {
            parse_columba_sam_file(path).unwrap_or_else(|e| {
                eprintln!("Failed to parse Columba SAM file: {}", e);
                std::process::exit(1);
            })
        } else {
            let columba_bin = args.columba_bin.as_ref().unwrap_or_else(|| {
                eprintln!("--columba-bin is required when --columba-index is supplied");
                std::process::exit(1);
            });
            let columba_index = args.columba_index.as_ref().unwrap_or_else(|| {
                eprintln!("--columba-index is required when --columba-bin is supplied");
                std::process::exit(1);
            });
            let candidate_edit_distance = candidate_edit_distance_bound(
                args.max_mismatches,
                args.max_bulges,
                args.max_bulge_size,
            );
            run_columba_candidate_generation(&ColumbaRunConfig {
                columba_bin,
                index_prefix: columba_index,
                queries: &queries,
                candidate_edit_distance,
                threads: args.threads,
                keep_sam: args.keep_columba_sam,
            })
            .unwrap_or_else(|e| {
                eprintln!("{}", e);
                std::process::exit(1);
            })
        };

        if batch_mode {
            let query_names: HashSet<&str> = queries
                .iter()
                .map(|query| query.query_name.as_str())
                .collect();
            let mut candidates_by_query = HashMap::new();
            for candidate in candidates {
                if !query_names.contains(candidate.query_name.as_str()) {
                    eprintln!(
                        "Columba candidate has unknown QNAME '{}' not present in --guides",
                        candidate.query_name
                    );
                    std::process::exit(1);
                }
                candidates_by_query
                    .entry(candidate.query_name.clone())
                    .or_insert_with(Vec::new)
                    .push(candidate);
            }

            for query in &queries {
                let guide_fwd = Arc::new(query.sequence.as_bytes().to_vec());
                let query_candidates = candidates_by_query
                    .remove(&query.query_name)
                    .unwrap_or_default();
                let hits = verify_columba_candidates(
                    &query_candidates,
                    &reference_records,
                    &guide_fwd,
                    &args,
                );
                report_filtered_hits_for_query(
                    hits,
                    &args.pam,
                    &query.query_name,
                    allow_noncanonical_pam,
                );
            }
        } else {
            let guide_fwd = Arc::new(queries[0].sequence.as_bytes().to_vec());
            let hits =
                verify_columba_candidates(&candidates, &reference_records, &guide_fwd, &args);
            report_filtered_hits(hits, &args.pam, allow_noncanonical_pam);
        }
        return;
    }

    let guide_fwd = Arc::new(queries[0].sequence.as_bytes().to_vec());
    let guide_len = guide_fwd.len();
    let reader = open_reference_reader(&args.reference);
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let record_id = record.id().to_string();

        let window_size = args.window_size.unwrap_or(guide_len * 4);
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();

        let hits: Vec<Hit> = windows
            .into_par_iter()
            .map_init(
                || AffineWavefronts::with_penalties(0, 3, 5, 1),
                |aligner, (i, end)| {
                    let window = uppercase_ascii_sequence(&seq[i..end]);
                    let mut hits = Vec::new();
                    if window.len() < guide_len {
                        return hits;
                    }
                    for strand in ['+', '-'] {
                        let oriented = if strand == '+' {
                            window.clone()
                        } else {
                            reverse_complement(&window)
                        };
                        if let Some((score, cigar, _, _, _, offset)) = scan_window(
                            aligner,
                            &guide_fwd,
                            &oriented,
                            args.max_mismatches,
                            args.max_bulges,
                            args.max_bulge_size,
                            args.min_match_fraction,
                            args.no_filter,
                        ) {
                            let span = verification::alignment_reference_span(&cigar);
                            let pos = if strand == '+' {
                                i + offset
                            } else {
                                end - offset - span
                            };
                            hits.push(build_verified_hit(
                                record_id.clone(),
                                &seq,
                                pos,
                                strand,
                                score,
                                cigar,
                                Arc::clone(&guide_fwd),
                                &oriented,
                                offset,
                                &args,
                            ));
                        }
                    }
                    hits
                },
            )
            .flatten()
            .collect();

        report_filtered_hits(hits, &args.pam, allow_noncanonical_pam);
    }
}
