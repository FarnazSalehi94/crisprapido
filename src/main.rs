use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use flate2::read::MultiGzDecoder;
use std::sync::Arc;
use std::collections::HashMap;
use clap::Parser;
use bio::io::fasta;
use sassy::profiles::Dna;
use sassy::search::Searcher;
// Remove the broken imports for now - we'll add correct ones later
// use sassy::{search, Alphabet, SearchConfig};
use std::fmt::Write;
use rayon::prelude::*;

mod cfd_score;

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

#[derive(Clone)]
struct Hit {
    ref_id: String,
    pos: usize,
    strand: char,
    score: i32,
    cigar: String,
    guide: Arc<Vec<u8>>,
    target_len: usize,
    max_mismatches: u32,
    max_bulges: u32,
    max_bulge_size: u32,
    cfd_score: Option<f64>,  // Add CFD score field
    target_seq: Vec<u8>,     // Add target sequence for CFD calculation
}

impl Hit {
    fn quality_score(&self) -> i32 {
        // Count matches and other stats
        let matches = self.cigar.chars()
            .filter(|&c| c == 'M' || c == '=')
            .count();
            
        let mismatches = self.cigar.chars()
            .filter(|&c| c == 'X')
            .count();
            
        let gaps = self.cigar.chars()
            .filter(|&c| c == 'I' || c == 'D')
            .count();
            
        // Higher score is better
        matches as i32 - (mismatches as i32) - (gaps as i32 * 2) - self.score
    }
    
    fn end_pos(&self) -> usize {
        // Calculate reference consumed bases
        let mut ref_consumed = 0;
        for c in self.cigar.chars() {
            match c {
                'M' | '=' | 'X' | 'D' => ref_consumed += 1,
                _ => {}
            }
        }
        self.pos + ref_consumed
    }
    
    fn overlaps_with(&self, other: &Hit) -> bool {
        self.strand == other.strand && 
        self.ref_id == other.ref_id &&
        self.pos < other.end_pos() && 
        other.pos < self.end_pos()
    }
}

fn report_hit(ref_id: &str, pos: usize, _len: usize, strand: char, 
              _score: i32, cigar: &str, guide: &[u8], target_len: usize,
              _max_mismatches: u32, _max_bulges: u32, _max_bulge_size: u32,
              target_seq: &[u8], pam: &str) {
    
    // Parse CIGAR to calculate positions and statistics
    let mut ref_pos = pos;
    let mut ref_consumed = 0;
    let mut query_pos = 0;
    let mut query_consumed = 0;
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap_size = 0;
    let mut max_gap_size = 0;
    let mut matches = 0;
    
    // Handle empty CIGAR (fallback to perfect match)
    let effective_cigar = if cigar.is_empty() {
        format!("{}=", guide.len())
    } else {
        cigar.to_string()
    };
    
    // Parse CIGAR string properly
    let mut chars = effective_cigar.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch.is_ascii_digit() {
            // Collect all digits
            let mut num_str = String::new();
            num_str.push(ch);
            while let Some(&next_ch) = chars.peek() {
                if next_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }
            
            // Get the operation
            if let Some(op) = chars.next() {
                if let Ok(count) = num_str.parse::<usize>() {
                    match op {
                        '=' | 'M' => {
                            matches += count;
                            ref_consumed += count;
                            query_consumed += count;
                            current_gap_size = 0;
                        },
                        'X' => {
                            // Only count mismatch if this position in the guide isn't N
                            for i in 0..count {
                                if query_pos + i < guide.len() && guide[query_pos + i] != b'N' {
                                    mismatches += 1;
                                }
                            }
                            ref_consumed += count;
                            query_consumed += count;
                            current_gap_size = 0;
                        },
                        'I' => {
                            if current_gap_size == 0 {
                                gaps += 1;
                            }
                            current_gap_size += count;
                            max_gap_size = max_gap_size.max(current_gap_size);
                            query_consumed += count;
                        },
                        'D' => {
                            if current_gap_size == 0 {
                                gaps += 1;
                            }
                            current_gap_size += count;
                            max_gap_size = max_gap_size.max(current_gap_size);
                            ref_consumed += count;
                        },
                        _ => {}
                    }
                }
            }
        }
    }
    
    // Calculate query start and end
    let query_start = 0; // Query always starts at 0 in local alignment
    let query_end = query_consumed;
    let query_length = guide.len(); // Total guide length
    
    // Calculate reference start and end
    let ref_start = ref_pos;
    let ref_end = ref_pos + ref_consumed;
    
    // Calculate adjusted score based on the alignment
    let mut adjusted_score = 0;
    let mut in_gap = false;
    for ch in effective_cigar.chars() {
        match ch {
            'X' => adjusted_score += 3,  // Mismatch penalty
            'I' | 'D' => {
                if !in_gap {
                    adjusted_score += 5;  // Gap opening penalty
                    in_gap = true;
                }
                adjusted_score += 1;  // Gap extension penalty
            },
            '=' | 'M' => in_gap = false,
            _ => {}
        }
    }
    
    // Calculate block length (total alignment length)
    let block_len = matches + mismatches + gaps;
    
    // Calculate CFD score
    let cfd_score = if !target_seq.is_empty() {
        cfd_score::get_cfd_score(guide, target_seq, &effective_cigar, pam)
    } else {
        None
    };

    // Add CFD score to output
    let cfd_tag = if let Some(score) = cfd_score {
        format!("\tcf:f:{:.4}", score)
    } else {
        String::new()
    };

    // Convert CIGAR to minimap2 format
    let minimap2_cigar = if effective_cigar.is_empty() {
        format!("{}=", guide.len())
    } else {
        convert_to_minimap2_cigar(&effective_cigar)
    };

    // Output in PAF format
    println!("Guide\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\tas:i:{}\tnm:i:{}\tng:i:{}\tbs:i:{}\tcg:Z:{}{}",
        query_length,                     // Query length (total guide length)
        query_start,                      // Query start (always 0 for local alignment)
        query_end,                        // Query end (bases consumed from query)
        strand,                           // Strand (+/-)
        ref_id,                           // Target sequence name
        target_len,                       // Full target sequence length
        ref_start,                        // Target start position
        ref_end,                          // Target end position
        matches,                          // Number of matches
        block_len,                        // Total alignment block length
        adjusted_score,                   // AS:i alignment score
        mismatches,                       // NM:i number of mismatches
        gaps,                             // NG:i number of gaps
        max_gap_size,                     // BS:i biggest gap size
        minimap2_cigar,                   // cg:Z CIGAR string
        cfd_tag                           // cf:f CFD score
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


    #[test]
    fn test_perfect_match_sassy() {
        let guide = b"ATCGATCGAT";
        let target = b"ATCGATCGAT";
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some());
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGAT";  // Single mismatch at position 5
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single mismatch");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "4=1X5=");
    }

    #[test]
    fn test_with_bulge_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGAATCGAT";  // Single base insertion after position 4
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single base bulge");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_sassy() {
        let guide =  b"ATCGATCGAT";
        let target = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        
        let result = scan_window_sassy(guide, target, 1, 1, 1, 0.75, false);
        assert!(result.is_none());
    }

    #[test]
    fn test_perfect_match_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let target = create_flanked_sequence(&mut rng, guide, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should match perfectly even with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "10=");
    }

    #[test]
    fn test_with_mismatches_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGAT";  // Single mismatch at position 5
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single mismatch with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert_eq!(cigar, "4=1X5=");
    }
    
    #[test]
    fn test_with_bulge_and_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGAATCGAT";  // Single base insertion after position 4
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..511], 1, 1, 1, 0.75, false);
        assert!(result.is_some(), "Should accept a single base bulge with flanks");
        let (_score, cigar, _mismatches, _gaps, _max_gap_size, _leading_dels) = result.unwrap();
        assert!(cigar.contains('I') || cigar.contains('D'), "Should contain an insertion or deletion");
    }

    #[test]
    fn test_too_many_differences_with_flanks_sassy() {
        let mut rng = SmallRng::seed_from_u64(42);
        let guide = b"ATCGATCGAT";
        let core = b"ATCGTTCGTT";  // Three mismatches at positions 5, 8, 9
        let target = create_flanked_sequence(&mut rng, core, 500);
        
        let result = scan_window_sassy(guide, &target[500..510], 1, 1, 1, 0.75, false);
        assert!(result.is_none(), "Should reject sequence with too many mismatches even with flanks");
    }
    
    #[test]
    fn test_hit_quality_scoring_and_filtering() {
        // Create Hit objects with different qualities
        let guide_seq = Arc::new(b"ATCGATCGAT".to_vec());
        
        // Perfect match hit
        let perfect_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 100,
            strand: '+',
            score: 0,
            cigar: "MMMMMMMMMM".to_string(),  // 10 perfect matches
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };
        
        // Hit with one mismatch
        let mismatch_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 105,  // Overlaps with perfect_hit
            strand: '+',
            score: 3,  // Higher score (worse)
            cigar: "MMMMXMMMMM".to_string(),  // 9 matches, 1 mismatch
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };
        
        // Hit with a bulge
        let bulge_hit = Hit {
            ref_id: "chr1".to_string(),
            pos: 110,  // Doesn't overlap with others
            strand: '+',
            score: 6,  // Even higher score (worse)
            cigar: "MMMDMMMMM".to_string(),  // Gap
            guide: Arc::clone(&guide_seq),
            target_len: 1000,
            max_mismatches: 4,
            max_bulges: 1,
            max_bulge_size: 2,
            cfd_score: None,
            target_seq: vec![],
        };
        
        // Verify overlapping detection
        assert!(perfect_hit.overlaps_with(&mismatch_hit), "Hits should overlap");
        assert!(mismatch_hit.overlaps_with(&perfect_hit), "Overlap should be symmetric");
        assert!(!perfect_hit.overlaps_with(&bulge_hit), "These hits shouldn't overlap");
        
        // Verify quality scoring
        assert!(perfect_hit.quality_score() > mismatch_hit.quality_score(),
                "Perfect match should have higher quality than mismatch");
        assert!(mismatch_hit.quality_score() > bulge_hit.quality_score(),
                "Mismatch should have higher quality than bulge");
                
        // Test end position calculation
        assert_eq!(perfect_hit.end_pos(), 110, "End position should be pos + matches");
        assert_eq!(mismatch_hit.end_pos(), 115, "End position includes mismatches");
        assert_eq!(bulge_hit.end_pos(), 119, "End position includes deletions");
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

    /// PAM sequence (to use for CFD scoring)
    #[arg(short = 'p', long, default_value = "GG")]
    pam: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short = 'b', long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,
    
    /// Minimum fraction of guide that must match (0.0-1.0)
    #[arg(short = 'f', long, default_value = "0.75")]
    min_match_fraction: f32,

    /// Path to mismatch scores file for CFD calculation
    #[arg(long, default_value = "mismatch_scores.txt")]
    mismatch_scores: PathBuf,

    /// Path to PAM scores file for CFD calculation
    #[arg(long, default_value = "pam_scores.txt")]
    pam_scores: PathBuf,

    /// Size of sequence window to scan (bp, default: 4x guide length)
    #[arg(short = 'w', long)]
    window_size: Option<usize>,

    /// Number of threads to use (default: number of logical CPUs)
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Disable all filtering (report every alignment)
    #[arg(long)]
    no_filter: bool,
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


// Convert SASSY's debug CIGAR format to standard format
fn parse_sassy_cigar_debug(debug_str: &str) -> String {
    let mut result = String::new();
    
    // Try to find CigarElem patterns in the debug string
    let mut pos = 0;
    while let Some(start) = debug_str[pos..].find("CigarElem { op: ") {
        let start = pos + start;
        
        // Extract operation type
        if let Some(op_start) = debug_str[start..].find("op: ") {
            let op_start = start + op_start + 4;
            if let Some(op_end) = debug_str[op_start..].find(",") {
                let op_end = op_start + op_end;
                let op = &debug_str[op_start..op_end];
                
                // Extract count
                if let Some(cnt_start) = debug_str[op_end..].find("cnt: ") {
                    let cnt_start = op_end + cnt_start + 5;
                    if let Some(cnt_end) = debug_str[cnt_start..].find(" }") {
                        let cnt_end = cnt_start + cnt_end;
                        if let Ok(count) = debug_str[cnt_start..cnt_end].parse::<u32>() {
                            let op_char = match op {
                                "Match" => '=',
                                "Sub" => 'X',
                                "Ins" => 'I',
                                "Del" => 'D',
                                _ => '='
                            };
                            result.push_str(&format!("{}{}", count, op_char));
                        }
                    }
                }
            }
        }
        pos = start + 1;
    }
    
    result
}

fn scan_window_sassy(
    guide: &[u8], 
    window: &[u8], 
    max_mismatches: u32, 
    max_bulges: u32, 
    max_bulge_size: u32,
    min_match_fraction: f32, 
    no_filter: bool
) -> Option<(i32, String, u32, u32, u32, usize)> {
    
    // Calculate maximum allowed errors
    let max_errors = (max_mismatches + max_bulges) as usize;
    
    // Create SASSY searcher with DNA profile
    let mut searcher: Searcher<Dna> = Searcher::new(false, None);
    
    // Convert window to a Vec so it implements SearchAble
    let window_vec = window.to_vec();
    
    // Search for matches using real SASSY
    let matches = searcher.search(guide, &window_vec, max_errors);

    if matches.is_empty() {
        return None;
    }
    
    // Take the best match (lowest cost)
    let best_match = matches.into_iter().min_by_key(|m| m.cost)?;
    
    // Convert SASSY results to CRISPRapido format
    let score = best_match.cost as i32;
    
    // Convert SASSY CIGAR to standard format
    let cigar_debug = format!("{:?}", best_match.cigar);
    let mut cigar_str = parse_sassy_cigar_debug(&cigar_debug);
    
    // If CIGAR parsing failed, create a fallback based on cost
    if cigar_str.is_empty() {
        if best_match.cost == 0 {
            // Perfect match
            cigar_str = format!("{}=", guide.len());
        } else {
            // Approximation: assume all errors are mismatches
            let matches = guide.len() - best_match.cost as usize;
            if matches > 0 {
                cigar_str = format!("{}={}", matches, best_match.cost);
                // Add 'X' for each mismatch
                for _ in 0..best_match.cost {
                    cigar_str.push('X');
                }
            } else {
                cigar_str = format!("{}X", guide.len());
            }
        }
    }
    
    // Calculate statistics from CIGAR
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut max_gap_size = 0;
    let mut current_gap_size = 0;
    let mut matches_count = 0;
    
    // Parse the CIGAR string to count operations
    let mut chars = cigar_str.chars().peekable();
    while let Some(ch) = chars.next() {
        if ch.is_ascii_digit() {
            let mut num_str = String::new();
            num_str.push(ch);
            while let Some(&next_ch) = chars.peek() {
                if next_ch.is_ascii_digit() {
                    num_str.push(chars.next().unwrap());
                } else {
                    break;
                }
            }
            
            if let Some(op) = chars.next() {
                if let Ok(count) = num_str.parse::<u32>() {
                    match op {
                        '=' | 'M' => {
                            matches_count += count as usize;
                            current_gap_size = 0;
                        },
                        'X' => {
                            mismatches += count;
                            current_gap_size = 0;
                        },
                        'I' | 'D' => {
                            if current_gap_size == 0 {
                                gaps += 1;
                            }
                            current_gap_size += count;
                            max_gap_size = max_gap_size.max(current_gap_size);
                        },
                        _ => {}
                    }
                }
            }
        }
    }
    
    // Apply filtering
    let non_n_positions = guide.iter().filter(|&&b| b != b'N').count();
    let match_percentage = if non_n_positions > 0 {
        (matches_count as f32 / non_n_positions as f32) * 100.0
    } else {
        0.0
    };

    if no_filter || (
        matches_count >= 1 && 
        match_percentage >= min_match_fraction * 100.0 && 
        mismatches <= max_mismatches && 
        gaps <= max_bulges && 
        max_gap_size <= max_bulge_size
    ) {
        // Find the actual position of the match in the window
        let mut actual_match_pos = 0;
        
        // For perfect matches, do exact substring search
        if best_match.cost == 0 {
            for i in 0..=(window.len().saturating_sub(guide.len())) {
                if &window[i..i+guide.len()] == guide {
                    actual_match_pos = i;
                    break;
                }
            }
        } else {
            // For matches with mismatches, find the best alignment position
            let mut best_score = std::i32::MAX;
            
            for i in 0..=(window.len().saturating_sub(guide.len())) {
                let mut score = 0;
                for j in 0..guide.len() {
                    if window[i + j] != guide[j] {
                        score += 1;
                    }
                }
                if score < best_score {
                    best_score = score;
                    actual_match_pos = i;
                }
            }
        }
        
        Some((score, cigar_str, mismatches, gaps, max_gap_size, actual_match_pos))
    } else {
        None
    }
}


fn main() {
    let args = Args::parse();
    
    // Initialize CFD score matrices
    if let Err(e) = cfd_score::init_score_matrices(
        args.mismatch_scores.to_str().unwrap_or("mismatch_scores.txt"),
        args.pam_scores.to_str().unwrap_or("pam_scores.txt")
    ) {
        eprintln!("Warning: CFD scoring disabled - {}", e);
    }
    
    // Print PAF header as comment (disabled)
    // println!("#Query\tQLen\tQStart\tQEnd\tStrand\tTarget\tTLen\tTStart\tTEnd\tMatches\tBlockLen\tMapQ\tTags");
    
    
    // Prepare guide sequences (forward and reverse complement)
    let guide_fwd = Arc::new(args.guide.as_bytes().to_vec());
    let guide_rc = Arc::new(reverse_complement(&guide_fwd));
    let guide_len = guide_fwd.len();

    // Set thread pool size if specified
    if let Some(n) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .expect("Failed to initialize thread pool");
    }

    // Process reference sequences
    // Create transparent reader that handles both plain and gzipped files
    let file = File::open(&args.reference).expect("Failed to open reference file");
    let reader: Box<dyn BufRead> = if args.reference.extension().map_or(false, |ext| ext == "gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };
    let reader = fasta::Reader::new(reader);
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq().to_vec();
        let seq_len = seq.len();
        let record_id = record.id().to_string();
        
        // Calculate window size as 4x guide length if not specified
        let window_size = args.window_size.unwrap_or(guide_len * 4);
        let step_size = window_size / 2;
        let windows: Vec<_> = (0..seq.len())
            .step_by(step_size)
            .map(|i| (i, (i + window_size).min(seq.len())))
            .collect();


        // Process windows in parallel and collect all hits
        let hits: Vec<Hit> = windows.into_par_iter()
            .map_init(
                || (),
                |_unit, (window_start, end)| {
                    let window = &seq[window_start..end];
                    if window.len() < guide_len { return None; }
    
                    // Try forward orientation
                    if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, match_offset_in_window)) = 
                        scan_window_sassy(&guide_fwd, window,
                                args.max_mismatches, args.max_bulges, args.max_bulge_size,
                                args.min_match_fraction, args.no_filter) {
                
                        // Calculate actual position in full sequence
                        let actual_pos = window_start + match_offset_in_window;
                
                        return Some(Hit {
                            ref_id: record_id.clone(),
                            pos: actual_pos,  // Use calculated position
                            strand: '+',
                            score,
                            cigar: cigar.clone(),
                            guide: Arc::clone(&guide_fwd),
                            target_len: seq_len,
                            max_mismatches: args.max_mismatches,
                            max_bulges: args.max_bulges,
                            max_bulge_size: args.max_bulge_size,
                            cfd_score: None,
                            target_seq: {
                                // Extract the actual target sequence for CFD calculation
                                let start = actual_pos;
                                let end = (actual_pos + guide_len).min(seq_len);
                                seq[start..end].to_vec()
                            },
                        });
                    }
            
                    // Try reverse complement orientation
                    if let Some((score, cigar, _mismatches, _gaps, _max_gap_size, match_offset_in_window)) = 
                        scan_window_sassy(&guide_rc, window,
                                args.max_mismatches, args.max_bulges, args.max_bulge_size,
                                args.min_match_fraction, args.no_filter) {
                
                        // Calculate actual position in full sequence
                        let actual_pos = window_start + match_offset_in_window;
                
                        return Some(Hit {
                            ref_id: record_id.clone(),
                            pos: actual_pos,  // Use calculated position
                            strand: '-',
                            score,
                            cigar: cigar.clone(),
                            guide: Arc::clone(&guide_rc),
                            target_len: seq_len,
                            max_mismatches: args.max_mismatches,
                            max_bulges: args.max_bulges,
                            max_bulge_size: args.max_bulge_size,
                            cfd_score: None,
                            target_seq: {
                                // Extract the actual target sequence for CFD calculation
                                let start = actual_pos;
                                let end = (actual_pos + guide_len).min(seq_len);
                                seq[start..end].to_vec()
                            },
                        });
                    }

                    None
                })
            .filter_map(|x| x)
            .collect();

        // Group hits by chromosome and strand
        let mut hits_by_group: HashMap<(String, char), Vec<Hit>> = HashMap::new();
        for hit in hits {
            hits_by_group.entry((hit.ref_id.clone(), hit.strand))
                        .or_insert_with(Vec::new)
                        .push(hit);
        }

        // For each group, filter overlapping hits
        for (_, mut group_hits) in hits_by_group {
            // Sort by position
            group_hits.sort_by_key(|hit| hit.pos);
            
            // Filter overlapping hits
            let _filtered_hits: Vec<Hit> = Vec::new(); // Unused, but keeping for future expansion
            let mut i = 0;
            while i < group_hits.len() {
                // Find all hits that overlap with the current one
                let mut best_idx = i;
                let mut best_quality = group_hits[i].quality_score();
                let mut j = i + 1;
                
                while j < group_hits.len() && group_hits[j].pos < group_hits[i].end_pos() {
                    if group_hits[j].overlaps_with(&group_hits[i]) {
                        let quality = group_hits[j].quality_score();
                        if quality > best_quality {
                            best_quality = quality;
                            best_idx = j;
                        }
                    }
                    j += 1;
                }
                
                // Add the best hit to filtered results
                let best_hit = &group_hits[best_idx];
                report_hit(
                    &best_hit.ref_id, 
                    best_hit.pos, 
                    best_hit.guide.len(), 
                    best_hit.strand, 
                    best_hit.score, 
                    &best_hit.cigar, 
                    &best_hit.guide, 
                    best_hit.target_len,
                    best_hit.max_mismatches, 
                    best_hit.max_bulges, 
                    best_hit.max_bulge_size,
                    &best_hit.target_seq,
                    &args.pam
                );
                
                // Move to the next non-overlapping hit
                i = j;
            }
        }
    }
}
