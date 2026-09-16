#![cfg(unix)]

use std::fs;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const GUIDE_ALPHA: &str = "GAGTCCGAGCAGAAGAAGAA";
const GUIDE_BETA: &str = "ATCGATCGATCGATCGATCG";

fn unique_test_dir() -> PathBuf {
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|duration| duration.as_nanos())
        .unwrap_or(0);
    let path = std::env::temp_dir().join(format!(
        "crisprapido-batch-equivalence-{}-{}",
        std::process::id(),
        now
    ));
    fs::create_dir(&path).unwrap();
    path
}

fn write_mock_columba(path: &Path) {
    let script = r#"#!/bin/sh
set -eu
query=""
output=""
edit_distance=""
while [ "$#" -gt 0 ]; do
    case "$1" in
        -f)
            shift
            query="$1"
            ;;
        -o)
            shift
            output="$1"
            ;;
        -e)
            shift
            edit_distance="$1"
            ;;
    esac
    shift
done
printf '%s\n' "$edit_distance" >> "$MOCK_COLUMBA_LOG"
printf '@HD\tVN:1.6\n' > "$output"
awk '
function emit(    position) {
    if (sequence == "GAGTCCGAGCAGAAGAAGAA") {
        position = 6
    } else if (sequence == "ATCGATCGATCGATCGATCG") {
        position = 38
    } else {
        exit 9
    }
    printf "%s\t0\tchr1\t%d\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0\n", name, position
}
/^>/ {
    if (name != "") emit()
    name = substr($1, 2)
    sequence = ""
    next
}
{
    sequence = sequence $0
}
END {
    if (name != "") emit()
}
' "$query" >> "$output"
"#;
    fs::write(path, script).unwrap();
    let mut permissions = fs::metadata(path).unwrap().permissions();
    permissions.set_mode(0o755);
    fs::set_permissions(path, permissions).unwrap();
}

fn run_crisprapido(
    root: &Path,
    reference: &Path,
    columba: &Path,
    index_prefix: &Path,
    invocation_log: &Path,
    guide_args: &[&str],
) -> Output {
    Command::new(env!("CARGO_BIN_EXE_crisprapido"))
        .current_dir(root)
        .arg("--reference")
        .arg(reference)
        .args(guide_args)
        .args(["--pam", "GG"])
        .args(["--max-mismatches", "2"])
        .args(["--max-bulges", "1"])
        .args(["--max-bulge-size", "1"])
        .args(["--min-match-fraction", "0.75"])
        .arg("--mismatch-scores")
        .arg(root.join("mismatch_scores.txt"))
        .arg("--pam-scores")
        .arg(root.join("pam_scores.txt"))
        .arg("--columba-bin")
        .arg(columba)
        .arg("--columba-index")
        .arg(index_prefix)
        .env("MOCK_COLUMBA_LOG", invocation_log)
        .output()
        .unwrap()
}

fn successful_stdout(output: Output) -> String {
    assert!(
        output.status.success(),
        "CRISPRapido failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

fn normalize_single_guide_output(output: &str, query_name: &str) -> Vec<String> {
    output
        .lines()
        .map(|line| {
            let (_, remainder) = line.split_once('\t').unwrap();
            format!("{}\t{}", query_name, remainder)
        })
        .collect()
}

#[test]
fn batch_mode_matches_repeated_single_guide_mode() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let reference = root.join("tests/fixtures/batch_reference.fa");
    let guides = root.join("tests/fixtures/batch_guides.fa");
    let dir = unique_test_dir();
    let columba = dir.join("columba");
    let index_prefix = dir.join("index");
    let batch_log = dir.join("batch.log");
    let single_log = dir.join("single.log");
    write_mock_columba(&columba);
    fs::write(dir.join("index.meta"), b"mock index").unwrap();

    let batch_output = successful_stdout(run_crisprapido(
        &root,
        &reference,
        &columba,
        &index_prefix,
        &batch_log,
        &["--guides", guides.to_str().unwrap()],
    ));
    let alpha_output = successful_stdout(run_crisprapido(
        &root,
        &reference,
        &columba,
        &index_prefix,
        &single_log,
        &["--guide", GUIDE_ALPHA],
    ));
    let beta_output = successful_stdout(run_crisprapido(
        &root,
        &reference,
        &columba,
        &index_prefix,
        &single_log,
        &["--guide", GUIDE_BETA],
    ));

    let mut normalized_batch: Vec<_> = batch_output.lines().map(str::to_string).collect();
    let mut normalized_singles = normalize_single_guide_output(&alpha_output, "guide_alpha");
    normalized_singles.extend(normalize_single_guide_output(&beta_output, "guide_beta"));
    normalized_batch.sort();
    normalized_singles.sort();

    assert_eq!(normalized_batch, normalized_singles);
    assert_eq!(
        fs::read_to_string(&batch_log).unwrap(),
        "3\n0\n",
        "one complete batch must be passed to exactly one e=3 and one e=0 Columba run"
    );
    assert_eq!(
        fs::read_to_string(&single_log).unwrap(),
        "3\n0\n3\n0\n",
        "the control performs two Columba runs per single guide"
    );
    assert!(alpha_output.lines().all(|line| line.starts_with("Guide\t")));
    assert!(beta_output.lines().all(|line| line.starts_with("Guide\t")));
    assert!(batch_output
        .lines()
        .any(|line| line.starts_with("guide_alpha\t")));
    assert!(batch_output
        .lines()
        .any(|line| line.starts_with("guide_beta\t")));

    fs::remove_dir_all(dir).unwrap();
}
