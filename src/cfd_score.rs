use std::collections::HashMap;
use std::fs::read_to_string;
use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref MISMATCH_SCORES: Mutex<HashMap<String, f64>> = Mutex::new(HashMap::new());
    static ref PAM_SCORES: Mutex<HashMap<String, f64>> = Mutex::new(HashMap::new());
    static ref INITIALIZED: Mutex<bool> = Mutex::new(false);
}

fn parse_scores(content: &str) -> Result<HashMap<String, f64>, &'static str> {
    let mut map = HashMap::new();
    for line in content.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();
        if parts.len() != 2 {
            continue;
        }
        let score: f64 = parts[1].parse().map_err(|_| "invalid score")?;
        map.insert(parts[0].to_string(), score);
    }
    Ok(map)
}

pub fn init_score_matrices(mis_path: &str, pam_path: &str) -> Result<(), &'static str> {
    let mut initialized = INITIALIZED.lock().unwrap();
    if *initialized {
        return Ok(());
    }
    *initialized = true;

    let mis_content = read_to_string(mis_path).map_err(|_| "failed to read mismatch_scores.txt")?;
    let mut mis = MISMATCH_SCORES.lock().unwrap();
    *mis = parse_scores(&mis_content).map_err(|_| "failed to parse mismatch_scores.txt")?;

    let pam_content = read_to_string(pam_path).map_err(|_| "failed to read pam_scores.txt")?;
    let mut pam = PAM_SCORES.lock().unwrap();
    *pam = parse_scores(&pam_content).map_err(|_| "failed to parse pam_scores.txt")?;

    Ok(())
}

pub fn calculate_cfd(spacer: &str, protospacer: &str, pam: &str) -> Result<f64, &'static str> {
    let spacer_rna = spacer.replace('T', "U");
    if spacer_rna.len() != 20 || protospacer.len() != 20 {
        return Err("sequences must be exactly 20 nt long");
    }

    let mis = MISMATCH_SCORES.lock().unwrap();
    if mis.is_empty() {
        return Err("score matrices not initialized. Call init_score_matrices first.");
    }

    let pams = PAM_SCORES.lock().unwrap();
    let pam_score = *pams.get(pam).unwrap_or(&0.0);

    let mut total_score = pam_score;
    for (i, (spacer_nt, proto_nt)) in spacer_rna.chars().zip(protospacer.chars()).enumerate() {
        let position = i + 1;
        let spacer_equiv = if spacer_nt.to_ascii_lowercase() == 'u' { 'T' } else { spacer_nt.to_ascii_uppercase() };
        let proto_equiv = proto_nt.to_ascii_uppercase();
        if spacer_equiv == proto_equiv {
            total_score *= 1.0;
            continue;
        }
        let key = format!("r{}:d{},{}", spacer_nt.to_ascii_uppercase(), reverse_complement_nt(proto_nt), position);
        let mismatch_score = *mis.get(&key).unwrap_or(&0.0);
        total_score *= mismatch_score;
    }
    Ok(total_score)
}

pub fn get_cfd_score(guide: &[u8], target_seq: &[u8], _cigar: &str, pam: &str) -> Option<f64> {
    let guide_str = std::str::from_utf8(guide).ok()?;
    let target_str = std::str::from_utf8(target_seq).ok()?;
    calculate_cfd(guide_str, target_str, pam).ok()
}

fn reverse_complement_nt(nucleotide: char) -> char {
    match nucleotide {
        'A' => 'T',
        'C' => 'G',
        'T' | 'U' => 'A',
        'G' => 'C',
        '-' => '-',
        _ => nucleotide,
    }
}

    }
    *initialized = true;

    let mis_content = read_to_string(mis_path).map_err(|_| "failed to read mismatch_scores.txt")?;
    let mut mis = MISMATCH_SCORES.lock().unwrap();
    *mis = parse_scores(&mis_content).map_err(|_| "failed to parse mismatch_scores.txt")?;

    let pam_content = read_to_string(pam_path).map_err(|_| "failed to read pam_scores.txt")?;
    let mut pam = PAM_SCORES.lock().unwrap();
    *pam = parse_scores(&pam_content).map_err(|_| "failed to parse pam_scores.txt")?;

    Ok(())
}

pub fn calculate_cfd(spacer: &str, protospacer: &str, pam: &str) -> Result<f64, &'static str> {
    let spacer_rna = spacer.replace('T', "U");
    if spacer_rna.len() != 20 || protospacer.len() != 20 {
        return Err("sequences must be exactly 20 nt long");
    }

    let mis = MISMATCH_SCORES.lock().unwrap();
    if mis.is_empty() {
        return Err("score matrices not initialized. Call init_score_matrices first.");
    }

    let pams = PAM_SCORES.lock().unwrap();
    let pam_score = *pams.get(pam).unwrap_or(&0.0);

    let mut total_score = pam_score;
    for (i, (spacer_nt, proto_nt)) in spacer_rna.chars().zip(protospacer.chars()).enumerate() {
        let position = i + 1;
        let spacer_equiv = if spacer_nt.to_ascii_lowercase() == 'u' { 'T' } else { spacer_nt.to_ascii_uppercase() };
        let proto_equiv = proto_nt.to_ascii_uppercase();
        if spacer_equiv == proto_equiv {
            total_score *= 1.0;
            continue;
        }
        let key = format!("r{}:d{},{}", spacer_nt.to_ascii_uppercase(), proto_nt.to_ascii_uppercase(), position);
        let mismatch_score = *mis.get(&key).unwrap_or(&0.0);
        total_score *= mismatch_score;
    }
    Ok(total_score)
}

<<<<<<< HEAD
#[cfg(test)]
mod cfd_unit_tests {
    use super::*;

    /// Test perfect match with canonical PAM gives score 1.0
    #[test]
    fn test_perfect_match() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let spacer = "GAAACAGTCGATTTTATCAC";
        let score = calculate_cfd(spacer, spacer, "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Perfect match should be 1.0, got {}", score);
    }

    /// Test that non-canonical PAMs reduce score
    #[test]
    fn test_pam_effects() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let spacer = "GAAACAGTCGATTTTATCAC";

        // GG should be highest (1.0 for perfect match)
        let gg_score = calculate_cfd(spacer, spacer, "GG").unwrap();
        let ag_score = calculate_cfd(spacer, spacer, "AG").unwrap();
        let cg_score = calculate_cfd(spacer, spacer, "CG").unwrap();

        assert!(gg_score > ag_score, "GG should score higher than AG");
        assert!(ag_score > cg_score, "AG should score higher than CG");
    }

    /// Test gap at position 1 (PAM-distal) has no penalty
    #[test]
    fn test_position_1_gap_no_penalty() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        // Gap in spacer at position 1
        let score = calculate_cfd("-AAACAGTCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Gap at position 1 should have no penalty, got {}", score);

        // Gap in protospacer at position 1
        let score = calculate_cfd("GAAACAGTCGATTTTATCAC", "-AAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 1.0).abs() < 1e-9, "Gap at position 1 should have no penalty, got {}", score);
    }

    /// Test gaps at other positions have penalties
    #[test]
    fn test_gaps_at_other_positions_have_penalties() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        // Gap at position 7 should have a penalty (score 0.0 per scoring matrix)
        let score = calculate_cfd("GAAACA-TCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!((score - 0.0).abs() < 1e-9, "Gap at position 7 should be 0.0, got {}", score);

        // Gap at position 2 should reduce score
        let score = calculate_cfd("G-AACAGTCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG").unwrap();
        assert!(score < 1.0, "Gap at position 2 should reduce score, got {}", score);
    }

    /// Test case insensitivity
    #[test]
    fn test_case_insensitivity() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let upper = "GAAACAGTCGATTTTATCAC";
        let lower = "gaaacagtcgattttatcac";

        let score_upper = calculate_cfd(upper, upper, "GG").unwrap();
        let score_mixed = calculate_cfd(lower, upper, "gg").unwrap();

        assert!((score_upper - score_mixed).abs() < 1e-9,
            "Case should not matter: upper={}, mixed={}", score_upper, score_mixed);
    }

    /// Test T/U equivalence
    #[test]
    fn test_t_u_equivalence() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let with_t = "GAAACAGTCGATTTTATCAC";
        let with_u = "GAAACAGUCGAUUUUAUCAC";

        let score_t = calculate_cfd(with_t, with_t, "GG").unwrap();
        let score_u = calculate_cfd(with_u, with_t, "GG").unwrap();

        assert!((score_t - score_u).abs() < 1e-9,
            "T and U should be equivalent: T={}, U={}", score_t, score_u);
    }

    /// Test that scores are in valid range [0, 1]
    #[test]
    fn test_scores_in_valid_range() {
        init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
            .expect("Failed to initialize scoring matrices");

        let base = "GAAACAGTCGATTTTATCAC";
        let bases = ['A', 'C', 'G', 'T'];

        // Test a variety of mismatches
        for pos in [1, 5, 10, 15, 20] {
            for &new_base in &bases {
                let mut proto: Vec<char> = base.chars().collect();
                proto[pos - 1] = new_base;
                let proto_str: String = proto.into_iter().collect();

                let score = calculate_cfd(base, &proto_str, "GG").unwrap();
                assert!(score >= 0.0 && score <= 1.0,
                    "Score out of range [0,1]: {} for pos {} -> {}", score, pos, new_base);
            }
        }
    }

    /// Test CIGAR-based alignment preparation
    #[test]
    fn test_prepare_aligned_sequences() {
        let guide = b"ATCGATCGATCGATCGATCG";
        let target = b"ATCGATCGATCGATCGATCG";

        // Perfect match
        let (spacer, proto) = prepare_aligned_sequences(guide, target, "20=");
        assert_eq!(spacer.len(), 20);
        assert_eq!(proto.len(), 20);
        assert_eq!(spacer, proto);

        // Insertion (gap in target)
        let (spacer, proto) = prepare_aligned_sequences(guide, &target[..19], "10=1I9=");
        assert!(proto.contains('-'), "Insertion should create gap in protospacer");

        // Deletion (gap in query)
        let guide_short = b"ATCGATCGACGATCGATCG"; // 19bp
        let (spacer, proto) = prepare_aligned_sequences(guide_short, target, "10=1D9=");
        assert!(spacer.contains('-'), "Deletion should create gap in spacer");
    }
}

/// Get CFD score using CIGAR-based alignment
///
/// This function takes raw guide and target sequences along with a CIGAR string
/// and produces properly aligned sequences with gap characters for CFD scoring.
pub fn get_cfd_score(guide: &[u8], target_seq: &[u8], cigar: &str, pam: &str) -> Option<f64> {
    // Use CIGAR to create properly aligned sequences with gap characters
    let (spacer, protospacer) = prepare_aligned_sequences(guide, target_seq, cigar);

    match calculate_cfd(&spacer, &protospacer, pam) {
        Ok(score) => Some(score),
        Err(_) => None,
    }
}
=======
pub fn get_cfd_score(guide: &[u8], target_seq: &[u8], _cigar: &str, pam: &str) -> Option<f64> {
    let guide_str = std::str::from_utf8(guide).ok()?;
    let target_str = std::str::from_utf8(target_seq).ok()?;
    calculate_cfd(guide_str, target_str, pam).ok()
}
>>>>>>> 9c12ff6 (fix: CFD matches Python reference (T->U spacer, revcomp keys, gap handling, unignore tests))
