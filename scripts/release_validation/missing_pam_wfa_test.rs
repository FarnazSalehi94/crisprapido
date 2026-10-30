#[cfg(test)]
mod independent_missing_pam_stage_probe {
    use super::*;
    #[test]
    fn missing_sites_pass_anchored_verification_at_frozen_limits() {
        for (name, guide, target) in [
            ("g1_exact_plus", "ACGTCAGTACGATCGTACGA", "ACGTCAGTACGATCGTACGAA"),
            ("g1_mm1_plus", "ACGTCAGTACGATCGTACGA", "ACGACAGTACGATCGTACGAA"),
            ("g2_mm1_plus", "GAGTCCGAGCAGAAGAAGAA", "GAGACCGAGCAGAAGAAGAAA"),
            ("CHM13_plus", "GACTAAACCTCTCAACAAAT", "GACTAAACCTCTCAACAAATT"),
            ("HG04184_plus", "GACTAAACCTCTCAACAAAT", "GACAAAACCTCTCAACAAATT"),
        ] {
            let mut aligner = AffineWavefronts::with_penalties(0, 3, 5, 1);
            let result = anchored_candidate_alignment(&mut aligner, guide.as_bytes(), target.as_bytes(), 4, 1, 1, 0.75, false);
            println!("{}: {:?}", name, result);
            assert!(result.is_some(), "{} did not verify", name);
        }
    }
}
