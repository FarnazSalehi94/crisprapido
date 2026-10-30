use std::{
    collections::{HashMap, HashSet},
    fs,
    path::Path,
    process::Command,
};

#[test]
fn recovered_pam_anchors_from_independent_indexed_panels() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/coverage_oracles_2026-10-29_v1");
    let panels = [
        (
            "indexed_synthetic",
            vec![
                ("g1", "g1_exact_plus", '+', 81),
                ("g1", "g1_mm1_minus", '-', 62),
                ("g1", "g1_mm1_plus", '+', 81),
                ("g2", "g2_mm1_minus", '-', 62),
                ("g2", "g2_mm1_plus", '+', 81),
            ],
        ),
        (
            "indexed_real_chm13",
            vec![("real1", "chm13_window_1", '+', 175)],
        ),
        (
            "indexed_12_windows",
            vec![
                ("real1", "CHM13_window_1", '+', 175),
                ("real1", "HG04184_mat_window_1", '+', 329),
            ],
        ),
    ];
    for (name, expected) in panels {
        let p = root.join(name);
        let output = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args([
                "-r",
                p.join("ref.fa").to_str().unwrap(),
                "-G",
                p.join("guides.fa").to_str().unwrap(),
                "--columba-sam",
                p.join("edit.sam").to_str().unwrap(),
                "-p",
                "GG",
                "-m",
                "4",
                "-b",
                "1",
                "-z",
                "1",
                "-f",
                "0.75",
                "-t",
                "1",
            ])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let reference = fs::read_to_string(p.join("ref.fa")).unwrap();
        let mut refs = HashMap::new();
        let mut current = "";
        for line in reference.lines() {
            if let Some(header) = line.strip_prefix('>') {
                current = header.split_whitespace().next().unwrap();
                refs.insert(current, String::new());
            } else {
                refs.get_mut(current).unwrap().push_str(line);
            }
        }
        let mut actual = HashSet::new();
        for line in String::from_utf8(output.stdout).unwrap().lines() {
            let f: Vec<_> = line.split('\t').collect();
            let start: usize = f[7].parse().unwrap();
            let end: usize = f[8].parse().unwrap();
            let strand = f[4].chars().next().unwrap();
            let anchor = if strand == '+' { end } else { start };
            let sequence = refs[f[5]].as_bytes();
            if strand == '+' {
                assert!(end + 3 <= sequence.len());
                assert!(matches!(sequence[end], b'A' | b'C' | b'G' | b'T'));
                assert_eq!(&sequence[end + 1..end + 3], b"GG");
            } else {
                assert!(start >= 3);
                assert_eq!(&sequence[start - 3..start - 1], b"CC");
                assert!(matches!(sequence[start - 1], b'A' | b'C' | b'G' | b'T'));
            }
            assert!(
                actual.insert((f[0].to_owned(), f[5].to_owned(), strand, anchor)),
                "duplicate anchor: {}",
                line
            );
        }
        for (guide, contig, strand, anchor) in expected {
            assert!(
                actual.contains(&(guide.to_owned(), contig.to_owned(), strand, anchor)),
                "{} missing {} {} {} {}",
                name,
                guide,
                contig,
                strand,
                anchor
            );
        }
    }
}
