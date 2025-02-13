use std::path::PathBuf;
use std::sync::Arc;
use clap::Parser;
use bio::io::fasta;
use libwfa2::affine_wavefront::AffineWavefronts;
use std::fmt::Write;
use rayon::prelude::*;

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',
        _ => b'N',  // Convert any unexpected bases to N
    }).collect()
}

fn report_hit(ref_id: &str, pos: usize, _len: usize, strand: char, 
              _score: i32, cigar: &str, guide: &[u8], target_len: usize) {
    // Parse CIGAR to handle leading indels and calculate positions
    let mut ref_pos = pos;
    let mut ref_consumed = 0;
    let mut leading_indels = true;
    let mut leading_dels = 0;
    
    // First pass: count leading deletions and find first match/mismatch
    for c in cigar.chars() {
        if leading_indels {
            match c {
                'D' => leading_dels += 1,
                'I' => (), // ignore leading insertions
                _ => leading_indels = false
            }
        }
    }
    
    // Adjust start position and trim leading/trailing deletions
    ref_pos += leading_dels;
    let trimmed_cigar: String = cigar.chars()
        .skip_while(|&c| c == 'D' || c == 'I')
        .collect::<String>()
        .trim_end_matches(|c| c == 'D' || c == 'I')
        .to_string();
    
    // Calculate alignment statistics, accounting for N positions
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    for c in trimmed_cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    mismatches += 1;
                }
                ref_consumed += 1;
                pos += 1;
            },
            'I' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
            },
            'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                ref_consumed += 1;
            },
            'M' | '=' => {
                current_gap_size = 0;
                ref_consumed += 1;
            },
            _ => ()
        }
    }

    // Recalculate score based on the trimmed alignment, accounting for N positions
    let mut adjusted_score = 0;
    let mut in_gap = false;
    let mut pos = 0;
    for c in trimmed_cigar.chars() {
        match c {
            'X' => {
                // Only count mismatch if this position in the guide isn't N
                if pos < guide.len() && guide[pos] != b'N' {
                    adjusted_score += 3;  // Mismatch penalty
                }
                pos += 1;
            },
            'I' | 'D' => {
                if !in_gap {
                    adjusted_score += 5;  // Gap opening penalty
                    in_gap = true;
                }
                adjusted_score += 1;  // Gap extension penalty
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                in_gap = false;
                pos += 1;
            },
            _ => ()
        }
    }

    // Count matches from CIGAR
    let matches = trimmed_cigar.chars()
        .filter(|&c| c == 'M' || c == '=')
        .count();
    
    // Calculate block length (matches + mismatches + indels)
    let block_len = trimmed_cigar.len();
    
    // Convert guide length to string once
    let guide_len = guide.len();
    
    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("Window scan debug:");
    debug!("  CIGAR before trim: {}", cigar);
    debug!("  CIGAR after trim: {}", trimmed_cigar);
    debug!("  N-adjusted mismatches: {} (max: 4)", mismatches);
    debug!("  Gaps: {} (max: 1)", gaps);
    debug!("  Max gap size: {} (max: 2)", max_gap_size);
    debug!("  Guide sequence: {}", String::from_utf8_lossy(guide));
    
    // Apply filters after trimming and N-adjustment
    if mismatches > 4 || gaps > 1 || max_gap_size > 2 {
        debug!("  Passes filters: false");
        debug!("");
        return;
    }
    debug!("  Passes filters: true");
    debug!("");

    println!("Guide\t{}\t0\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}", 
        guide_len,                        // Query length
        guide_len,                        // Query end
        strand,                           // Strand (+/-)
        ref_id,                           // Target sequence name
        target_len,                       // Full target sequence length
        ref_pos,                          // Target start
        ref_pos + ref_consumed,           // Target end
        matches,                          // Number of matches
        block_len,                        // Total alignment block length
        adjusted_score,                   // AS:i alignment score
        mismatches,                       // NM:i number of mismatches
        gaps,                             // NG:i number of gaps
        max_gap_size,                     // BS:i biggest gap size
        convert_to_minimap2_cigar(&trimmed_cigar) // cg:Z CIGAR string
    );
}
#[cfg(test)]
use rand::{SeedableRng, RngCore, rngs::SmallRng};

#[cfg(test)]
mod tests {
    use super::*;

    fn generate_random_seq(rng: &mut SmallRng, length: usize) -> Vec<u8> {
        let bases = b"ACGT";
        (0..length)
            .map(|_| bases[rng.next_u32() as usize % 4])
            .collect()
    }

    fn create_flanked_sequence(rng: &mut SmallRng, core: &[u8], flank_size: usize) -> Vec<u8> {
        let mut seq = generate_random_seq(rng, flank_size);
        seq.extend_from_slice(core);
        seq.extend(generate_random_seq(rng, flank_size));
        seq
    }

    fn setup_aligner() -> AffineWavefronts {
        AffineWavefronts::with_penalties(
            0,     // match score
            3,     // mismatch penalty
            5,     // gap opening penalty
            1      // gap extension penalty
        )
    }

    #[test]
    fn test_perfect_match() {
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        
        let result = scan_window(&mut aligner, guide, target, 1, 1, 1);
        assert!(result.is_some());
        let (_score, cigar) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5
        
        let result = scan_window(&mut aligner, guide, target);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4
        
        let result = scan_window(&mut aligner, guide, target);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences() {
        let mut aligner = setup_aligner();
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        
        let result = scan_window(&mut aligner, guide, target);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510]);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar) = result.unwrap();
        assert_eq!(cigar, "MMMMMMMMMM");
    }

    #[test]
    fn test_with_mismatches_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510]);
        assert!(result.is_some(), "Should accept a single mismatch with flanks");
        let (_score, cigar) = result.unwrap();
        assert_eq!(cigar, "MMMMXMMMMM");
    }

    #[test]
    fn test_with_bulge_and_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..511]);
        assert!(result.is_some(), "Should accept a single base bulge with flanks");
        let (_score, cigar) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks() {
        let mut rng = SmallRng::seed_from_u64(42);
        let mut aligner = setup_aligner();
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window(&mut aligner, guide, &target[500..510]);
        assert!(result.is_none(), "Should reject sequence with too many mismatches even with flanks");
    }
}

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file (-r)
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM) (-g)
    #[arg(short, long)]
    guide: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,

    /// Size of sequence window to scan (bp)
    #[arg(short = 'w', long, default_value = "1000")]
    window_size: usize,
}

fn convert_to_minimap2_cigar(cigar: &str) -> String {
    let mut result = String::new();
    let mut count = 0;
    let mut current_op = None;

    for c in cigar.chars() {
        let op = match c {
            'M' => '=',
            'X' | 'I' | 'D' => c,
            _ => continue,
        };

        if Some(op) == current_op {
            count += 1;
        } else {
            if count > 0 {
                write!(result, "{}{}", count, current_op.unwrap()).unwrap();
            }
            current_op = Some(op);
            count = 1;
        }
    }

    if count > 0 && current_op.is_some() {
        write!(result, "{}{}", count, current_op.unwrap()).unwrap();
    }

    result
}

fn scan_window(aligner: &mut AffineWavefronts, guide: &[u8], window: &[u8], 
               max_mismatches: u32, max_bulges: u32, max_bulge_size: u32) 
               -> Option<(i32, String, u32, u32, u32)> {
    aligner.align(window, guide);  // Target sequence first, then guide sequence
    let score = aligner.score();
    let cigar = String::from_utf8_lossy(aligner.cigar()).to_string();
    
    // Count mismatches ignoring N positions in guide
    let mut n_adjusted_mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut pos = 0;
    
    for c in cigar.chars() {
        match c {
            'X' => {
                if pos < guide.len() && guide[pos] != b'N' {
                    n_adjusted_mismatches += 1;
                }
                pos += 1;
            },
            'I' | 'D' => {
                current_gap_size += 1;
                if current_gap_size == 1 {
                    gaps += 1;
                }
                max_gap_size = max_gap_size.max(current_gap_size);
                if c == 'I' { pos += 1; }
            },
            'M' | '=' => {
                current_gap_size = 0;
                pos += 1;
            },
            _ => ()
        }
    }

    // Debug macro for development/testing
    macro_rules! debug {
        ($($arg:tt)*) => {
            #[cfg(feature = "debug")]
            eprintln!($($arg)*);
        }
    }

    debug!("CIGAR: {}, N-adjusted Mismatches: {}, Gaps: {}, Max gap size: {}", 
           cigar, n_adjusted_mismatches, gaps, max_gap_size);

    #[cfg(test)]
    let (_max_m, _max_b, _max_bs) = (1, 1, 1);  // Stricter thresholds for tests
    
    #[cfg(not(test))]
    let (_max_m, _max_b, _max_bs) = (max_mismatches, max_bulges, max_bulge_size);

    // Always return the alignment result - filtering happens in report_hit
    Some((score, cigar, n_adjusted_mismatches, gaps, max_gap_size))
}

fn main() {
    let args = Args::parse();
    
    // Print PAF header as comment (disabled)
    // println!("#Query\tQLen\tQStart\tQEnd\tStrand\tTarget\tTLen\tTStart\tTEnd\tMatches\tBlockLen\tMapQ\tTags");
    
    // Import required WFA2 types
    use libwfa2::affine_wavefront::{AlignmentSpan, AffineWavefronts};

    // Set up WFA parameters with CRISPR-specific penalties and end-free alignment
    let mut aligner = AffineWavefronts::with_penalties(
        0,     // match score
        3,     // mismatch penalty
        5,     // gap opening penalty
        1      // gap extension penalty
    );
    
    // Configure end-free alignment with single-gap allowance
    aligner.set_alignment_span(AlignmentSpan::EndsFree {
        pattern_begin_free: 1,  // Start of guide RNA
        pattern_end_free: 1,    // End of guide RNA
        text_begin_free: 1,     // Start of genomic sequence
        text_end_free: 1        // End of genomic sequence
    });
    
    // Prepare guide sequences (forward and reverse complement)
    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();

    // Process reference sequences
    let reader = fasta::Reader::from_file(args.reference).expect("Failed to read FASTA file");
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let seq_len = seq.len();
        let record_id = record.id().to_string();
        
        // Generate all window positions first
        let step_size = args.window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + args.window_size).min(seq.len())))
            .collect();

        // Process windows in parallel
        windows.into_par_iter().for_each(|(i, end)| {
            let window = &seq[i..end];
            if window.len() < guide_len { return; }

            // Create a new aligner for each thread
            let mut aligner = AffineWavefronts::with_penalties(0, 3, 5, 1);
            
            // Try forward orientation
            if let Some((score, cigar, _mismatches, _gaps, _max_gap_size)) = 
                scan_window(&mut aligner, &guide_fwd, window,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size) {
                report_hit(&record_id, i, guide_len, '+', score, &cigar, &guide_fwd, seq_len);
            }
            
            // Try reverse complement orientation
            if let Some((score, cigar, _mismatches, _gaps, _max_gap_size)) = 
                scan_window(&mut aligner, &guide_rc, window,
                          args.max_mismatches, args.max_bulges, args.max_bulge_size) {
                report_hit(&record_id, i, guide_len, '-', score, &cigar, &guide_rc, seq_len);
            }
        });
    }
}
