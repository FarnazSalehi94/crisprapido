use std::{
    fs,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static NEXT: AtomicUsize = AtomicUsize::new(0);
const GUIDE: &str = "GAAACAGTCGATTTTATCAC";

#[test]
fn invalid_cfd_thresholds_fail_before_search() {
    for value in ["NaN", "inf", "-0.1", "1.1", "garbage"] {
        let out = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args([
                "-r",
                "not-opened.fa",
                "-g",
                GUIDE,
                &format!("--min-cfd={value}"),
            ])
            .output()
            .unwrap();
        assert!(!out.status.success());
        assert!(out.stdout.is_empty());
        assert!(String::from_utf8_lossy(&out.stderr).contains("min-cfd"));
    }
}

#[test]
fn cfd_thresholds_preserve_exact_sites_batch_ids_and_missing_scores() {
    let dir = std::env::temp_dir().join(format!(
        "crispr-min-cfd-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir(&dir).unwrap();
    fs::write(dir.join("ref.fa"), format!(">chr\n{GUIDE}AGG\n")).unwrap();
    fs::write(dir.join("guides.fa"), format!(">one\n{GUIDE}\n")).unwrap();
    fs::write(
        dir.join("in.sam"),
        "one\t0\tchr\t1\t255\t20M\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:0\n",
    )
    .unwrap();
    let run = |opts: &[&str]| {
        let out = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args([
                "-r",
                dir.join("ref.fa").to_str().unwrap(),
                "-G",
                dir.join("guides.fa").to_str().unwrap(),
                "--columba-sam",
                dir.join("in.sam").to_str().unwrap(),
                "-t",
                "1",
                "-m",
                "0",
                "-b",
                "0",
            ])
            .args(opts)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        String::from_utf8(out.stdout).unwrap()
    };
    let default = run(&[]);
    assert_eq!(default.lines().count(), 1);
    assert!(default.starts_with("one\t"));
    for tau in ["0", "0.1", "0.2", "0.5", "1"] {
        assert_eq!(default, run(&["--min-cfd", tau]));
    }
    let missing = dir.join("missing.txt");
    let unknown = run(&[
        "--min-cfd",
        "0.5",
        "--mismatch-scores",
        missing.to_str().unwrap(),
    ]);
    assert_eq!(unknown.lines().count(), 1);
    assert!(!unknown.contains("cf:f:"));
    fs::remove_dir_all(dir).unwrap();
}
