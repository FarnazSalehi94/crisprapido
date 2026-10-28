use std::fs;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

const GUIDE: &str = "ACGTCAGTACGATCGTACGA";
static NEXT: AtomicUsize = AtomicUsize::new(0);
struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "crisprapido-pam-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
    fn run(&self, fasta: &str, sam: &str, options: &[&str]) -> String {
        fs::write(self.0.join("reference.fa"), fasta).unwrap();
        fs::write(self.0.join("candidates.sam"), sam).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args([
                "-r",
                self.0.join("reference.fa").to_str().unwrap(),
                "-g",
                GUIDE,
                "--columba-sam",
                self.0.join("candidates.sam").to_str().unwrap(),
                "-p",
                "GG",
                "-t",
                "1",
            ])
            .args(options)
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout).unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
fn rc(s: &str) -> String {
    s.chars()
        .rev()
        .map(|b| match b {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect()
}
fn sam(name: &str, reverse: bool, pos: usize, cigar: &str) -> String {
    format!(
        "guide\t{}\t{}\t{}\t255\t{}\t*\t0\t0\t*\t*\tNM:i:0\tAS:i:0\n",
        if reverse { 16 } else { 0 },
        name,
        pos + 1,
        cigar
    )
}
fn tag<'a>(f: &'a [&str], prefix: &str) -> &'a str {
    f.iter().find_map(|s| s.strip_prefix(prefix)).unwrap()
}

#[test]
fn canonical_policy_and_actual_cfd_on_both_strands() {
    let fixture = Fixture::new();
    let mut fa = String::new();
    let mut candidates = String::new();
    for reverse in [false, true] {
        for pam in ["AGG", "TGG", "CGG", "GGG", "GGA", "ACC", "AAA"] {
            let name = format!("{}_{}", if reverse { "minus" } else { "plus" }, pam);
            let oriented = format!("TTTTT{}{}TTTTT", GUIDE, pam);
            fa.push_str(&format!(
                ">{}\n{}\n",
                name,
                if reverse { rc(&oriented) } else { oriented }
            ));
            candidates.push_str(&sam(&name, reverse, if reverse { 8 } else { 5 }, "20M"));
        }
    }
    let strict = fixture.run(&fa, &candidates, &["-b", "0", "-m", "0", "-z", "1"]);
    assert_eq!(strict.lines().count(), 8);
    for line in strict.lines() {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(tag(&f, "cf:f:"), "1.0000");
        assert_eq!(&f[2..4], &["0", "20"]);
        assert!(f[5].ends_with("GG"));
        assert_eq!(
            (f[7], f[8]),
            if f[4] == "+" {
                ("5", "25")
            } else {
                ("8", "28")
            }
        );
    }
    let permissive = fixture.run(
        &fa,
        &candidates,
        &["-b", "0", "-m", "0", "--allow-noncanonical-pam"],
    );
    assert_eq!(permissive.lines().count(), 14);
    for line in permissive.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f[5].ends_with("GGA") {
            assert_eq!(tag(&f, "cf:f:"), "0.0694");
        }
        if f[5].ends_with("ACC") {
            assert_eq!(tag(&f, "cf:f:"), "0.0000");
        }
    }
    let bypass = fixture.run(&fa, &candidates, &["-b", "0", "-m", "0", "--no-filter"]);
    assert_eq!(bypass.lines().count(), 14);
}

#[test]
fn incomplete_or_ambiguous_pams_are_not_scored() {
    let fixture = Fixture::new();
    let mut fa = String::new();
    let mut candidates = String::new();
    for reverse in [false, true] {
        for (i, pam) in ["", "G", "GG", "NGG", "ANG", "RAG"].iter().enumerate() {
            let name = format!("edge_{}_{}", reverse, i);
            let oriented = format!("{}{}", GUIDE, pam);
            fa.push_str(&format!(
                ">{}\n{}\n",
                name,
                if reverse { rc(&oriented) } else { oriented }
            ));
            candidates.push_str(&sam(
                &name,
                reverse,
                if reverse { pam.len() } else { 0 },
                "20M",
            ));
        }
    }
    assert!(fixture
        .run(&fa, &candidates, &["-b", "0", "-m", "0"])
        .is_empty());
    let output = fixture.run(
        &fa,
        &candidates,
        &["-b", "0", "-m", "0", "--allow-noncanonical-pam"],
    );
    assert_eq!(output.lines().count(), 12);
    assert!(!output.contains("cf:f:"));
}

fn consumption(cigar: &str) -> (usize, usize) {
    let (mut n, mut q, mut r) = (0, 0, 0);
    for c in cigar.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d as usize;
        } else {
            match c {
                '=' | 'X' | 'M' => {
                    q += n;
                    r += n;
                }
                'I' => q += n,
                'D' => r += n,
                _ => panic!("unexpected op"),
            };
            n = 0;
        }
    }
    (q, r)
}
#[test]
fn mismatch_and_indel_coordinates_remain_consistent_on_both_strands() {
    let fixture = Fixture::new();
    let mut fa = String::new();
    let mut candidates = String::new();
    let mut mismatch = GUIDE.to_string();
    mismatch.replace_range(9..10, "A");
    let mut multiple = mismatch.clone();
    multiple.replace_range(17..18, "T");
    let mut dna = GUIDE.to_string();
    dna.insert(10, 'T');
    let mut rna = GUIDE.to_string();
    rna.remove(10);
    let leading = format!("T{}", GUIDE);
    let trailing = format!("{}T", GUIDE);
    for reverse in [false, true] {
        for (name, target, cigar) in [
            ("mismatch", mismatch.as_str(), "20M"),
            ("multiple", multiple.as_str(), "20M"),
            ("dna", dna.as_str(), "10M1D10M"),
            ("rna", rna.as_str(), "10M1I9M"),
            ("leading", leading.as_str(), "1D20M"),
            ("trailing", trailing.as_str(), "20M1D"),
        ] {
            let name = format!("{}_{}", name, reverse);
            let oriented = format!("TTTTT{}AGGTTTTT", target);
            fa.push_str(&format!(
                ">{}\n{}\n",
                name,
                if reverse { rc(&oriented) } else { oriented }
            ));
            candidates.push_str(&sam(&name, reverse, if reverse { 8 } else { 5 }, cigar));
        }
    }
    let output = fixture.run(&fa, &candidates, &["-b", "1", "-z", "1", "-m", "2"]);
    assert_eq!(output.lines().count(), 12);
    for line in output.lines() {
        let f: Vec<_> = line.split('\t').collect();
        let (q, r) = consumption(tag(&f, "cg:Z:"));
        assert_eq!(q, 20);
        assert_eq!(&f[2..4], &["0", "20"]);
        let s: usize = f[7].parse().unwrap();
        let e: usize = f[8].parse().unwrap();
        assert_eq!(e - s, r);
        let reference = fa
            .lines()
            .skip_while(|x| *x != format!(">{}", f[5]))
            .nth(1)
            .unwrap();
        let pam = if f[4] == "+" {
            reference[e..e + 3].to_string()
        } else {
            rc(&reference[s - 3..s])
        };
        assert_eq!(&pam[1..], "GG");
        if f[5].starts_with("dna_") {
            assert!(tag(&f, "cg:Z:").contains('D'));
        }
        if f[5].starts_with("rna_") {
            assert!(tag(&f, "cg:Z:").contains('I'));
        }
        if f[5].starts_with("mismatch_") {
            assert!(tag(&f, "cg:Z:").contains('X'));
        }
    }
}

#[test]
fn direct_scan_checks_both_strands_in_one_window() {
    let fixture = Fixture::new();
    let sequence = format!("{}AGGTTTTTCCT{}", GUIDE, rc(GUIDE));
    fs::write(fixture.0.join("both.fa"), format!(">both\n{}\n", sequence)).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
        .args([
            "-r",
            fixture.0.join("both.fa").to_str().unwrap(),
            "-g",
            GUIDE,
            "-b",
            "0",
            "-m",
            "0",
            "-w",
            "80",
            "-t",
            "1",
        ])
        .output()
        .unwrap();
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert_eq!(text.lines().count(), 2, "{}", text);
    assert!(text.lines().any(|l| l.split('\t').nth(4) == Some("+")));
    assert!(text.lines().any(|l| l.split('\t').nth(4) == Some("-")));
    assert!(text.lines().all(|l| l.contains("cf:f:1.0000")));
}

#[test]
fn invalid_requested_pam_is_rejected_before_opening_reference() {
    for pam in ["NGG", "G", "GN", "", "GGG"] {
        let output = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args(["-r", "nonexistent.fa", "-g", GUIDE, "-p", pam])
            .output()
            .unwrap();
        assert_eq!(output.status.code(), Some(2));
        assert!(String::from_utf8_lossy(&output.stderr).contains("exactly two DNA bases"));
    }
}
