//! Experimental adapter for the separately built Columba prototype.
use crate::columba::ColumbaQuery;
use crate::{Args, Hit};
use std::collections::{HashMap, HashSet};
use std::io::Write;
use std::process::{Command, Stdio};
use std::sync::Arc;

pub(crate) fn search(
    args: &Args,
    queries: &[ColumbaQuery],
) -> Result<HashMap<String, Vec<Hit>>, String> {
    if args.max_bulges > 1 || args.max_bulge_size > 2 {
        return Err("experimental backend supports -b 0..1 and -z 0..2 only".into());
    }
    if !args.min_match_fraction.is_finite() || !(0.0..=1.0).contains(&args.min_match_fraction) {
        return Err("experimental backend requires -f between 0 and 1".into());
    }
    if args.allow_noncanonical_pam
        || args.no_filter
        || args.columba_sam.is_some()
        || args.columba_bin.is_some()
        || args.keep_columba_sam
        || args.window_size.is_some()
        || args.threads.is_some_and(|n| n != 1)
    {
        return Err("experimental backend does not support noncanonical PAMs, --no-filter, Columba SAM/bin options, --window-size, or threads other than 1".into());
    }
    let index = args
        .columba_index
        .as_ref()
        .ok_or("experimental backend requires --columba-index")?;
    let binary = args
        .pam_search_bin
        .as_ref()
        .ok_or("experimental backend requires --pam-search-bin")?;
    let mut fasta = String::new();
    let mut guides = HashMap::new();
    for q in queries {
        let seq = q.sequence.to_ascii_uppercase();
        if seq.is_empty()
            || !seq.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            || q.query_name.chars().any(char::is_whitespace)
        {
            return Err(
                "experimental backend requires nonempty A/C/G/T guides and whitespace-free IDs"
                    .into(),
            );
        }
        fasta.push_str(&format!(">{}\n{}\n", q.query_name, seq));
        if guides
            .insert(q.query_name.clone(), Arc::new(seq.into_bytes()))
            .is_some()
        {
            return Err("duplicate guide ID".into());
        }
    }
    let python = option_env!("CRISPR_PAM_PYTHON").unwrap_or("/usr/bin/python3");
    if !std::path::Path::new(python).is_absolute() {
        return Err("CRISPR_PAM_PYTHON must be an absolute build-time path".into());
    }
    let mut child = Command::new(python)
        .arg("-I")
        .arg(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/research/pam_anchored_poc/integration/run_snapshot.py"
        ))
        .args(["--binary"])
        .arg(binary)
        .arg("--reference")
        .arg(&args.reference)
        .arg("--index")
        .arg(index)
        .args([
            "--mismatches",
            &args.max_mismatches.to_string(),
            "--bulges",
            &args.max_bulges.to_string(),
            "--size",
            &args.max_bulge_size.to_string(),
            "--fraction",
            &args.min_match_fraction.to_string(),
            "--pam",
            &args.pam,
        ])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("cannot start experimental backend: {e}"))?;
    child
        .stdin
        .take()
        .unwrap()
        .write_all(fasta.as_bytes())
        .map_err(|e| format!("cannot write guides: {e}"))?;
    let output = child
        .wait_with_output()
        .map_err(|e| format!("experimental backend failed: {e}"))?;
    eprint!("{}", String::from_utf8_lossy(&output.stderr));
    if !output.status.success() {
        return Err("experimental backend failed; no results reported".into());
    }
    let text = String::from_utf8(output.stdout).map_err(|_| "invalid backend output encoding")?;
    let mut hits: HashMap<String, Vec<Hit>> = HashMap::new();
    let mut seen = HashSet::new();
    let mut lines = text.lines();
    if !lines
        .next()
        .is_some_and(|s| s.starts_with("guide\tcontig\tstrand\tstart\tend\tanchor\tpam\tcigar\t"))
    {
        return Err("invalid experimental backend header".into());
    }
    for line in lines {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() != 16 {
            return Err("invalid experimental backend row".into());
        }
        let num = |i: usize| {
            f[i].parse::<usize>()
                .map_err(|_| "invalid experimental backend number".to_string())
        };
        let guide = guides
            .get(f[0])
            .ok_or("unknown experimental backend guide")?
            .clone();
        let strand = match f[2] {
            "+" => '+',
            "-" => '-',
            _ => return Err("invalid strand".into()),
        };
        let (pos, end, anchor) = (num(3)?, num(4)?, num(5)?);
        if anchor != if strand == '+' { end } else { pos }
            || !seen.insert((f[0].to_owned(), f[1].to_owned(), strand, anchor))
        {
            return Err("invalid or duplicate PAM anchor".into());
        }
        let target_len = num(15)?;
        if end > target_len
            || end < pos
            || f[14].len() != end - pos
            || f[6].len() != 3
            || &f[6][1..] != args.pam
        {
            return Err("invalid backend target or PAM".into());
        }
        let ops = crate::verification::cigar_operations(f[7]).ok_or("invalid backend CIGAR")?;
        let mut cost = 0i32;
        let mut gap = false;
        for (n, op) in ops {
            match op {
                'X' => {
                    cost += 3 * n as i32;
                    gap = false;
                }
                '=' | 'M' => gap = false,
                'I' | 'D' => {
                    if !gap {
                        cost += 5;
                    }
                    cost += n as i32;
                    gap = true;
                }
                _ => return Err("invalid backend CIGAR operation".into()),
            }
        }
        hits.entry(f[0].into()).or_default().push(Hit {
            ref_id: f[1].into(),
            pos,
            strand,
            score: -cost,
            cigar: f[7].into(),
            guide,
            target_len,
            max_mismatches: args.max_mismatches,
            max_bulges: args.max_bulges,
            max_bulge_size: args.max_bulge_size,
            cfd_score: None,
            target_seq: f[14].as_bytes().to_vec(),
            pam_seq: Some(f[6][1..].into()),
        });
    }
    Ok(hits)
}
