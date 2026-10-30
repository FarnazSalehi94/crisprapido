use std::{collections::HashSet, fs, path::Path, process::Command};

#[test]
fn all_sixteen_broader_oracle_pam_sites_are_recovered() {
    let base =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/broader_coverage_1eb2f40_2026-10-30");
    let cases = [
        ("synthetic", "frozen2", "frozen2_I_10_minus", '-', 48),
        ("synthetic", "frozen2", "frozen2_I_10_plus", '+', 65),
        ("synthetic", "repeat", "repeat_D_19_minus", '-', 48),
        ("synthetic", "repeat", "repeat_D_5_minus", '-', 44),
        ("synthetic", "repeat", "repeat_D_5_plus", '+', 67),
        ("synthetic", "repeat", "repeat_I_19_minus", '-', 45),
        ("synthetic", "repeat", "repeat_I_19_minus", '-', 46),
        ("synthetic", "repeat", "repeat_I_5_minus", '-', 43),
        ("synthetic", "repeat", "repeat_I_5_minus", '-', 48),
        ("synthetic", "repeat", "repeat_I_5_plus", '+', 65),
        ("synthetic", "repeat", "repeat_exact_plus", '+', 70),
        ("synthetic", "repeat", "repeat_four_a_minus", '-', 46),
        ("synthetic", "repeat", "repeat_four_b_minus", '-', 48),
        ("synthetic", "repeat", "repeat_mm_0_minus", '-', 42),
        ("synthetic", "repeat", "repeat_mm_0_minus", '-', 43),
        ("real", "r3_I", "GRCh38_window_1", '+', 136),
    ];
    for panel in ["synthetic", "real"] {
        let dir = base.join(panel);
        let tmp = std::env::temp_dir().join(format!(
            "crispr-bounded-{}-{}.sam",
            std::process::id(),
            panel
        ));
        let combined = format!(
            "{}{}",
            fs::read_to_string(dir.join("edit.sam")).unwrap(),
            fs::read_to_string(dir.join("exact.sam")).unwrap()
        );
        fs::write(&tmp, combined).unwrap();
        let out = Command::new(env!("CARGO_BIN_EXE_crisprapido"))
            .args([
                "-r",
                dir.join("ref.fa").to_str().unwrap(),
                "-G",
                dir.join("guides.fa").to_str().unwrap(),
                "--columba-sam",
                tmp.to_str().unwrap(),
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
        let _ = fs::remove_file(&tmp);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let mut found = HashSet::new();
        for line in String::from_utf8(out.stdout).unwrap().lines() {
            let f: Vec<_> = line.split('\t').collect();
            let strand = f[4].chars().next().unwrap();
            let anchor: usize = f[if strand == '+' { 8 } else { 7 }].parse().unwrap();
            assert!(
                found.insert((f[0].to_owned(), f[5].to_owned(), strand, anchor)),
                "duplicate anchor: {}",
                line
            );
        }
        for &(case_panel, guide, contig, strand, anchor) in &cases {
            if case_panel == panel {
                assert!(
                    found.contains(&(guide.to_string(), contig.to_string(), strand, anchor)),
                    "missing {} {} {} {} {}",
                    panel,
                    guide,
                    contig,
                    strand,
                    anchor
                );
            }
        }
    }
}
