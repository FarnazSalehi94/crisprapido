extern crate crisprapido;

use crisprapido::cfd_score;
use std::fs;
use std::path::Path;

/// Ensure score files exist for testing
fn ensure_score_files() {
    // We'll skip the file creation since you already have the files
    // Make sure mismatch_scores.txt and pam_scores.txt are in the root directory
}

/// Test function to validate our CFD score calculation against Python implementation
#[test]
fn test_cfd_score_against_python() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Define test cases based on the Python implementation
    // Format: (spacer, protospacer, pam, expected_python_score)
    let test_cases = vec![
        // Perfect match with GG PAM
        ("ATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "GG", 1.0),
        
        // Single mismatch at position 1 (PAM-distal) with GG PAM
        // This tests the special case for position 1
        ("ATCGATCGATCGATCGATCG", "TTCGATCGATCGATCGATCG", "GG", 1.0),
        
        // Single mismatch at position 10 with GG PAM
        ("ATCGATCGATCGATCGATCG", "ATCGATCGAACGATCGATCG", "GG", 0.857142857),
        
        // Single mismatch at position 20 (PAM-proximal) with GG PAM
        ("ATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCT", "GG", 0.7),
        
        // Multiple mismatches with GG PAM
        ("ATCGATCGATCGATCGATCG", "TTCGATCGAACGATCGATCT", "GG", 0.5999999999),
        
        // Perfect match with non-canonical PAM (AG)
        ("ATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "AG", 0.25925925899999996),
        
        // Perfect match with non-canonical PAM (TG)
        ("ATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "TG", 0.038961038999999996),
        
        // Test with gap/bulge at position 1 (PAM-distal)
        // This tests the special case for position 1 gap
        ("-TCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "GG", 1.0),
        
        // Test with gap/bulge at other positions
        ("ATCG-TCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "GG", 0.0),
        
        // Real example from paper
        ("GAAACAGTCGATTTTATCAC", "GAAACAGTCGATTTTATCAC", "GG", 1.0),
        ("GAAACAGTCGATTTTATCAC", "GAAACAGGCGATTTTATCAC", "GG", 0.733333333),
    ];
    
    // Run test cases
    for (i, (spacer, protospacer, pam, expected_python_score)) in test_cases.iter().enumerate() {
        // Calculate CIGAR string based on the alignment
        let mut cigar = String::with_capacity(spacer.len());
        if spacer.len() == protospacer.len() {
            for (s, p) in spacer.chars().zip(protospacer.chars()) {
                if s == '-' || p == '-' {
                    if s == '-' {
                        cigar.push('D'); // Genomic target base aligned to a spacer gap
                    } else {
                        cigar.push('I'); // Guide base aligned to a genomic target gap
                    }
                } else if s == p {
                    cigar.push('M');
                } else {
                    cigar.push('X');
                }
            }
        } else {
            panic!("Test case {}: Spacer and protospacer must have the same length", i+1);
        }
        
        println!("Test case {}: spacer={}, protospacer={}, cigar={}, pam={}",
                 i+1, spacer, protospacer, cigar, pam);
        
        // Test approach 1: Use direct calculation
        match cfd_score::calculate_cfd(spacer, protospacer, pam) {
            Ok(score) => {
                println!("  Direct calculation: CFD score = {:.6} (expected {:.6})",
                         score, expected_python_score);
                
                // Allow some floating point tolerance
                let tolerance = 0.0001;
                assert!((score - expected_python_score).abs() < tolerance,
                        "Test case {} direct calculation failed: got {:.6} but expected {:.6}",
                        i+1, score, expected_python_score);
            },
            Err(e) => {
                println!("  Direct calculation failed: {}", e);
                // If we expect a score of 0.0, it's okay if the calculation fails
                if *expected_python_score > 0.0 {
                    panic!("Test case {} direct calculation failed unexpectedly: {}", i+1, e);
                }
            }
        }
        
        // get_cfd_score accepts ungapped guide/target sequences and applies the
        // CIGAR itself. Pre-aligned gap cases above belong to calculate_cfd.
        if spacer.contains('-') || protospacer.contains('-') {
            continue;
        }

        // Test approach 2: reconstruct an ungapped alignment via CIGAR.
        let spacer_bytes = spacer.as_bytes();
        let protospacer_bytes = protospacer.as_bytes();
        
        match cfd_score::get_cfd_score(spacer_bytes, protospacer_bytes, &cigar, pam) {
            Some(score) => {
                println!("  CIGAR calculation: CFD score = {:.6} (expected {:.6})",
                         score, expected_python_score);
                
                // Allow some floating point tolerance
                let tolerance = 0.0001;
                assert!((score - expected_python_score).abs() < tolerance,
                        "Test case {} CIGAR calculation failed: got {:.6} but expected {:.6}",
                        i+1, score, expected_python_score);
            },
            None => {
                println!("  CIGAR calculation failed");
                // If we expect a score of 0.0, it's okay if the calculation fails
                if *expected_python_score > 0.0 {
                    panic!("Test case {} CIGAR calculation failed unexpectedly", i+1);
                }
            }
        }
        
        println!("");
    }
}

/// Test with varied number of mismatches and their positions
#[test]
fn test_positional_effects() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Base perfect-match sequence
    let base_spacer = "ATCGATCGATCGATCGATCG";
    
    // Test mismatches at each position (1-based)
    for pos in 1..=20 {
        // Create a sequence with a single mismatch at position pos
        let mut protospacer = base_spacer.to_string();
        let base_at_pos = protospacer.chars().nth(pos-1).unwrap();
        let mismatch_base = match base_at_pos {
            'A' => 'T',
            'T' => 'G',
            'G' => 'C',
            'C' => 'A',
            _ => panic!("Unexpected base: {}", base_at_pos),
        };
        
        // Replace the base at position pos-1 (0-indexed)
        let mut chars: Vec<char> = protospacer.chars().collect();
        chars[pos-1] = mismatch_base;
        protospacer = chars.into_iter().collect();
        
        // Calculate CIGAR string
        let mut cigar = "M".repeat(20);
        let mut cigar_chars: Vec<char> = cigar.chars().collect();
        cigar_chars[pos-1] = 'X';
        cigar = cigar_chars.into_iter().collect();
        
        // Calculate CFD score
        if let Some(score) = cfd_score::get_cfd_score(
            base_spacer.as_bytes(), 
            protospacer.as_bytes(), 
            &cigar, 
            "GG"
        ) {
            println!("Position {} mismatch: CFD score = {:.6}", pos, score);
        } else {
            println!("CFD score calculation failed for position {}", pos);
        }
    }
}

/// Test with different PAM sequences
#[test]
fn test_pam_effects() {
    // Initialize score matrices
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");
    
    // Base perfect-match sequence
    let spacer = "ATCGATCGATCGATCGATCG";
    let protospacer = "ATCGATCGATCGATCGATCG";
    let cigar = "M".repeat(20);
    
    // Test different PAM sequences
    let pams = vec![
        ("GG", 1.0),               // Canonical PAM - highest score
        ("AG", 0.25925925899999996), // Non-canonical but functional
        ("TG", 0.038961038999999996), // Non-canonical but somewhat functional
        ("CG", 0.107142857),       // Non-canonical
        ("AT", 0.0),               // Non-functional
    ];
    
    for (pam, expected_score) in pams {
        if let Some(score) = cfd_score::get_cfd_score(
            spacer.as_bytes(), 
            protospacer.as_bytes(), 
            &cigar, 
            pam
        ) {
            println!("PAM {}: CFD score = {:.6} (expected {:.6})", 
                     pam, score, expected_score);
            
            // Allow some floating point tolerance
            let tolerance = 0.0001;
            assert!((score - expected_score).abs() < tolerance,
                    "PAM {} test failed: got {:.6} but expected {:.6}",
                    pam, score, expected_score);
        } else {
            println!("CFD score calculation failed for PAM {}", pam);
        }
    }
}

/// Test with bulges (insertions/deletions)
#[test]
fn test_bulge_effects() {
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("Failed to initialize scoring matrices");

    let guide = b"ATCGATCGATCGATCGATCG";

    for pos in 0..20 {
        let mut target_with_guide_insertion = guide.to_vec();
        target_with_guide_insertion.remove(pos);
        let guide_insertion_cigar = format!("{}I{}", "=".repeat(pos), "=".repeat(19 - pos));
        assert!(
            cfd_score::get_cfd_score(
                guide,
                &target_with_guide_insertion,
                &guide_insertion_cigar,
                "GG",
            )
            .is_some(),
            "guide-only I at aligned position {} should have a 20-column CFD representation",
            pos + 1
        );

        let mut target_with_genomic_insertion = guide.to_vec();
        target_with_genomic_insertion.insert(pos, b'A');
        let genomic_insertion_cigar = format!("{}D{}", "=".repeat(pos), "=".repeat(20 - pos));
        assert_eq!(
            cfd_score::get_cfd_score(
                guide,
                &target_with_genomic_insertion,
                &genomic_insertion_cigar,
                "GG",
            ),
            None,
            "genomic-only D at aligned position {} produces 21 columns",
            pos + 1
        );
    }
}
