use std::fs::File;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Once;
use lazy_static::lazy_static;
use std::sync::Mutex;

lazy_static! {
    static ref MISMATCH_SCORES: Mutex<Option<HashMap<String, f64>>> = Mutex::new(None);
    static ref PAM_SCORES: Mutex<Option<HashMap<String, f64>>> = Mutex::new(None);
    static ref INIT: Once = Once::new();
}

/// Initialize score matrices if not already loaded
pub fn init_score_matrices(mismatch_path: &str, pam_path: &str) -> Result<(), String> {
    INIT.call_once(|| {
        let mm_scores = parse_scoring_matrix(mismatch_path)
            .map_err(|e| format!("Failed to load mismatch scores: {}", e));
        let pam_scores = parse_scoring_matrix(pam_path)
            .map_err(|e| format!("Failed to load PAM scores: {}", e));

        if let (Ok(mm), Ok(pam)) = (mm_scores, pam_scores) {
            *MISMATCH_SCORES.lock().unwrap() = Some(mm);
            *PAM_SCORES.lock().unwrap() = Some(pam);
        }
    });

    // Verify the matrices are loaded
    let mm_loaded = MISMATCH_SCORES.lock().unwrap().is_some();
    let pam_loaded = PAM_SCORES.lock().unwrap().is_some();

    if mm_loaded && pam_loaded {
        Ok(())
    } else {
        Err("Failed to initialize scoring matrices".to_string())
    }
}

/// Parse scoring matrix from space-delimited file
fn parse_scoring_matrix(file_path: &str) -> Result<HashMap<String, f64>, String> {
    // Open file
    let file = File::open(file_path)
        .map_err(|e| format!("Cannot open {}: {}", file_path, e))?;
    
    // Read file
    let reader = BufReader::new(file);
    let mut matrix = HashMap::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Error reading line: {}", e))?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 2 {
            let score = parts[1].parse::<f64>()
                .map_err(|e| format!("Invalid score format: {}", e))?;
            matrix.insert(parts[0].to_string(), score);
        }
    }
    Ok(matrix)
}

/// Get reverse complement of a single nucleotide (supports bulges)
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

/// Align spacer and target sequence for CFD calculation
fn prepare_aligned_sequences(guide: &[u8], target: &[u8], cigar: &str) -> (String, String) {
    let mut spacer = String::with_capacity(20);
    let mut protospacer = String::with_capacity(20);
    
    let mut guide_pos = 0;
    let mut target_pos = 0;
    
    for c in cigar.chars() {
        match c {
            'M' | '=' => {
                if guide_pos < guide.len() && target_pos < target.len() {
                    spacer.push(char::from(guide[guide_pos]));
                    protospacer.push(char::from(target[target_pos]));
                    guide_pos += 1;
                    target_pos += 1;
                }
            },
            'X' => {
                if guide_pos < guide.len() && target_pos < target.len() {
                    spacer.push(char::from(guide[guide_pos]));
                    protospacer.push(char::from(target[target_pos]));
                    guide_pos += 1;
                    target_pos += 1;
                }
            },
            'I' => {
                if guide_pos < guide.len() {
                    spacer.push(char::from(guide[guide_pos]));
                    protospacer.push('-');
                    guide_pos += 1;
                }
            },
            'D' => {
                if target_pos < target.len() {
                    spacer.push('-');
                    protospacer.push(char::from(target[target_pos]));
                    target_pos += 1;
                }
            },
            _ => {}
        }
    }
    
    // Pad to 20nt if needed
    while spacer.len() < 20 {
        spacer.push('-');
    }
    while protospacer.len() < 20 {
        protospacer.push('-');
    }
    
    // Truncate to 20nt if longer
    let spacer = spacer[0..20].to_string();
    let protospacer = protospacer[0..20].to_string();
    
    (spacer, protospacer)
}

/// Calculate CFD score
pub fn calculate_cfd(spacer: &str, protospacer: &str, pam: &str) -> Result<f64, String> {
    // Check for expected input lengths
    if spacer.len() != 20 || protospacer.len() != 20 {
        return Err(format!("Incorrect input sequence length, expected 20nt for both spacer and protospacer"));
    }

    // Get locked references to scoring matrices
    let mm_scores_lock = MISMATCH_SCORES.lock().unwrap();
    let pam_scores_lock = PAM_SCORES.lock().unwrap();
    
    // Verify matrices are initialized
    let mm_scores = mm_scores_lock.as_ref()
        .ok_or_else(|| "Mismatch scores not initialized".to_string())?;
    let pam_scores = pam_scores_lock.as_ref()
        .ok_or_else(|| "PAM scores not initialized".to_string())?;
    
    // Pre-process sequences
    let spacer_list: Vec<char> = spacer.to_uppercase().replace('T', "U").chars().collect();
    let protospacer_list: Vec<char> = protospacer.to_uppercase().replace('T', "U")
                                               .chars().collect();
    
    // Calculate CFD score for alignment by nucleotide
    let mut score = 1.0;
    for (i, &nt) in protospacer_list.iter().enumerate() {
        if spacer_list[i] == nt {
            // No penalty for perfect match
            continue;
        } else if i == 0 && (spacer_list[i] == '-' || nt == '-'){
            // No penalty for gap at most PAM-distal nucleotide
            continue;
        } else {
            // Incorporate score for given RNA-DNA basepair at this position
            let key = format!("r{}:d{},{}", spacer_list[i], reverse_complement_nt(nt), i + 1);
            let penalty = mm_scores.get(&key)
                .ok_or_else(|| format!("Invalid basepair: {}", key))?;
            score *= penalty;
        }
    }
    
    // Incorporate PAM score
    let pam_upper = pam.to_uppercase();
    let pam_penalty = pam_scores.get(&pam_upper)
        .ok_or_else(|| format!("Invalid PAM: {}", pam_upper))?;
    score *= pam_penalty;
    
    Ok(score)
}

/// Get CFD score for a hit
pub fn get_cfd_score(guide: &[u8], target: &[u8], cigar: &str, pam: &str) -> Option<f64> {
    // Prepare aligned sequences for CFD calculation
    let (spacer, protospacer) = prepare_aligned_sequences(guide, target, cigar);
    
    // Calculate CFD score
    match calculate_cfd(&spacer, &protospacer, pam) {
        Ok(score) => Some(score),
        Err(e) => {
            eprintln!("CFD score calculation error: {}", e);
            None
        }
    }
}
