//! CRISPR guide RNA scoring algorithms
//! 
//! This module implements various scoring algorithms for CRISPR guide-target pairs:
//! - CFD (Cutting Frequency Determination) score
//! - Elevation score (Doench lab)

use std::collections::HashMap;

/// Calculate CFD score for a guide-target pair
/// 
/// The CFD score estimates the cutting efficiency based on mismatch patterns
/// between the guide RNA and target DNA.
/// 
/// # Arguments
/// * `guide` - The guide RNA sequence (without PAM)
/// * `target` - The target DNA sequence
/// * `cigar` - The CIGAR string representing the alignment
/// 
/// # Returns
/// A score between 0.0 and 1.0, where higher values indicate better predicted cutting efficiency
pub fn calculate_cfd_score(guide: &[u8], target: &[u8], cigar: &str) -> f32 {
    // CFD score parameters from Doench et al. 2016
    // These are the penalty weights for different types of mismatches at different positions
    let mut mm_scores: HashMap<(usize, char, char), f32> = HashMap::new();
    
    // Position-specific mismatch penalties (simplified version)
    // Format: (position, guide_base, target_base) -> penalty
    mm_scores.insert((0, 'A', 'C'), 0.25);
    mm_scores.insert((0, 'A', 'G'), 0.25);
    mm_scores.insert((0, 'A', 'T'), 0.25);
    mm_scores.insert((0, 'C', 'A'), 0.25);
    mm_scores.insert((0, 'C', 'G'), 0.25);
    mm_scores.insert((0, 'C', 'T'), 0.25);
    mm_scores.insert((0, 'G', 'A'), 0.25);
    mm_scores.insert((0, 'G', 'C'), 0.25);
    mm_scores.insert((0, 'G', 'T'), 0.25);
    mm_scores.insert((0, 'T', 'A'), 0.25);
    mm_scores.insert((0, 'T', 'C'), 0.25);
    mm_scores.insert((0, 'T', 'G'), 0.25);
    
    // Middle positions have different penalties
    for pos in 1..guide.len()-1 {
        mm_scores.insert((pos, 'A', 'C'), 0.1);
        mm_scores.insert((pos, 'A', 'G'), 0.1);
        mm_scores.insert((pos, 'A', 'T'), 0.1);
        mm_scores.insert((pos, 'C', 'A'), 0.1);
        mm_scores.insert((pos, 'C', 'G'), 0.1);
        mm_scores.insert((pos, 'C', 'T'), 0.1);
        mm_scores.insert((pos, 'G', 'A'), 0.1);
        mm_scores.insert((pos, 'G', 'C'), 0.1);
        mm_scores.insert((pos, 'G', 'T'), 0.1);
        mm_scores.insert((pos, 'T', 'A'), 0.1);
        mm_scores.insert((pos, 'T', 'C'), 0.1);
        mm_scores.insert((pos, 'T', 'G'), 0.1);
    }
    
    // Last position has different penalties
    let last_pos = guide.len() - 1;
    mm_scores.insert((last_pos, 'A', 'C'), 0.15);
    mm_scores.insert((last_pos, 'A', 'G'), 0.15);
    mm_scores.insert((last_pos, 'A', 'T'), 0.15);
    mm_scores.insert((last_pos, 'C', 'A'), 0.15);
    mm_scores.insert((last_pos, 'C', 'G'), 0.15);
    mm_scores.insert((last_pos, 'C', 'T'), 0.15);
    mm_scores.insert((last_pos, 'G', 'A'), 0.15);
    mm_scores.insert((last_pos, 'G', 'C'), 0.15);
    mm_scores.insert((last_pos, 'G', 'T'), 0.15);
    mm_scores.insert((last_pos, 'T', 'A'), 0.15);
    mm_scores.insert((last_pos, 'T', 'C'), 0.15);
    mm_scores.insert((last_pos, 'T', 'G'), 0.15);
    
    // Bulge penalties
    let bulge_penalty = 0.5;
    
    // Start with perfect score
    let mut score = 1.0;
    
    // Parse CIGAR string to identify mismatches and bulges
    let mut guide_pos = 0;
    let mut target_pos = 0;
    
    for c in cigar.chars() {
        match c {
            'M' | '=' => {
                // Match - no penalty
                guide_pos += 1;
                target_pos += 1;
            },
            'X' => {
                // Mismatch - apply position-specific penalty
                if guide_pos < guide.len() && target_pos < target.len() {
                    let guide_base = guide[guide_pos] as char;
                    let target_base = target[target_pos] as char;
                    
                    if let Some(penalty) = mm_scores.get(&(guide_pos, guide_base, target_base)) {
                        score *= 1.0 - penalty;
                    } else {
                        // Default penalty if not in the table
                        score *= 0.8;
                    }
                }
                guide_pos += 1;
                target_pos += 1;
            },
            'I' | 'D' => {
                // Bulge - apply bulge penalty
                score *= 1.0 - bulge_penalty;
                if c == 'I' {
                    guide_pos += 1;
                } else {
                    target_pos += 1;
                }
            },
            _ => {}
        }
    }
    
    score
}

/// Calculate Elevation score for a guide-target pair
/// 
/// The Elevation score is a machine learning-based model from the Doench lab
/// that predicts off-target effects by considering the pattern of mismatches
/// and their positions.
/// 
/// # Arguments
/// * `guide` - The guide RNA sequence (without PAM)
/// * `target` - The target DNA sequence
/// * `cigar` - The CIGAR string representing the alignment
/// 
/// # Returns
/// A score between 0.0 and 1.0, where lower values indicate less off-target activity
pub fn calculate_elevation_score(guide: &[u8], target: &[u8], cigar: &str) -> f32 {
    // Simplified Elevation score implementation
    // In a real implementation, this would use the full Elevation model
    
    // Count mismatches and their positions
    let mut mismatch_count = 0;
    let mut mismatch_positions = Vec::new();
    let mut bulge_count = 0;
    
    let mut guide_pos = 0;
    let mut target_pos = 0;
    
    for c in cigar.chars() {
        match c {
            'M' | '=' => {
                guide_pos += 1;
                target_pos += 1;
            },
            'X' => {
                mismatch_count += 1;
                mismatch_positions.push(guide_pos);
                guide_pos += 1;
                target_pos += 1;
            },
            'I' | 'D' => {
                bulge_count += 1;
                if c == 'I' {
                    guide_pos += 1;
                } else {
                    target_pos += 1;
                }
            },
            _ => {}
        }
    }
    
    // Position-dependent weighting (seed region is more important)
    let mut position_weight = 0.0;
    for pos in &mismatch_positions {
        if *pos < 8 {  // Seed region (positions 0-7)
            position_weight += 0.2;
        } else {
            position_weight += 0.1;
        }
    }
    
    // Bulges are heavily penalized
    let bulge_weight = bulge_count as f32 * 0.3;
    
    // Calculate final score (lower is better for off-targets)
    let raw_score = 1.0 - (position_weight + bulge_weight);
    
    // Ensure score is between 0 and 1
    raw_score.max(0.0).min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_cfd_perfect_match() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        let cigar = "MMMMMMMMMM";
        
        let score = calculate_cfd_score(guide, target, cigar);
        assert_eq!(score, 1.0);
    }
    
    #[test]
    fn test_cfd_with_mismatches() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Mismatch at position 4
        let cigar = "MMMMXMMMMM";
        
        let score = calculate_cfd_score(guide, target, cigar);
        assert!(score < 1.0);
        assert!(score > 0.0);
    }
    
    #[test]
    fn test_cfd_with_multiple_mismatches() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Mismatches at positions 4, 8
        let cigar = "MMMMXMMXMX";
        
        let score = calculate_cfd_score(guide, target, cigar);
        // Multiple mismatches should result in a lower score
        assert!(score < 0.8);
    }
    
    #[test]
    fn test_cfd_with_bulge() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Insertion after position 4
        let cigar = "MMMMIMMMMMM";
        
        let score = calculate_cfd_score(guide, target, cigar);
        // Bulges are heavily penalized
        assert!(score <= 0.5);
    }
    
    #[test]
    fn test_elevation_perfect_match() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        let cigar = "MMMMMMMMMM";
        
        let score = calculate_elevation_score(guide, target, cigar);
        assert_eq!(score, 1.0);
    }
    
    #[test]
    fn test_elevation_with_mismatches() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Mismatch at position 4
        let cigar = "MMMMXMMMMM";
        
        let score = calculate_elevation_score(guide, target, cigar);
        assert!(score < 1.0);
        assert!(score > 0.0);
    }
}
