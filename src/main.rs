fn main() {
    println!("Hello, world!");
}
use std::path::PathBuf;
use clap::Parser;
use bio::io::fasta;
use libwfa2::{wavefront_aligner_attr_default, wavefront_aligner_new, MemoryMode, DistanceMetric};
use libwfa2::affine_wavefront::AffineWavefronts;

#[derive(Parser)]
#[command(author, version, about = "CRISPR guide RNA off-target scanner")]
struct Args {
    /// Input reference FASTA file
    #[arg(short, long)]
    reference: PathBuf,

    /// Guide RNA sequence (without PAM)
    #[arg(short, long)]
    guide: String,

    /// Maximum number of mismatches allowed
    #[arg(short, long, default_value = "4")]
    max_mismatches: u32,

    /// Maximum number of bulges allowed
    #[arg(short, long, default_value = "1")]
    max_bulges: u32,

    /// Maximum size of each bulge in bp
    #[arg(short = 'z', long, default_value = "2")]
    max_bulge_size: u32,
}

fn scan_window(aligner: &mut AffineWavefronts, guide: &[u8], window: &[u8]) -> Option<(i32, String)> {
    aligner.align(guide, window);
    let score = aligner.score();
    let cigar = String::from_utf8_lossy(aligner.cigar()).to_string();
    
    // Count mismatches and gaps from CIGAR
    let mut mismatches = 0;
    let mut gaps = 0;
    let mut current_gap = 0;
    
    for c in cigar.chars() {
        match c {
            'X' => mismatches += 1,
            'I' | 'D' => {
                current_gap += 1;
                if current_gap == 1 {
                    gaps += 1;
                }
            },
            'M' => current_gap = 0,
            _ => ()
        }
    }

    if mismatches <= 4 && gaps <= 1 {
        Some((score, cigar))
    } else {
        None
    }
}

fn main() {
    let args = Args::parse();
    
    // Set up WFA parameters
    let mut attributes = wavefront_aligner_attr_default();
    attributes.memory_mode = MemoryMode::Ultralow;
    attributes.distance_metric = DistanceMetric::GapAffine;
    attributes.affine_penalties.mismatch = 3;
    attributes.affine_penalties.gap_opening = 5;
    attributes.affine_penalties.gap_extension = 1;

    let mut aligner = AffineWavefronts::new_with_attributes(&attributes);
    
    // Prepare guide sequence
    let guide = args.guide.as_bytes();
    let guide_len = guide.len();
    
    // Process reference sequences
    let reader = fasta::Reader::from_file(args.reference).expect("Failed to read FASTA file");
    
    for result in reader.records() {
        let record = result.expect("Error during FASTA record parsing");
        let seq = record.seq();
        
        // Scan windows
        for (i, window) in seq.windows(guide_len + 3).enumerate() {
            // Check if window ends with NGG PAM
            if window.len() >= guide_len + 3 
               && window[guide_len + 1] == b'G' 
               && window[guide_len + 2] == b'G' {
                
                if let Some((score, cigar)) = scan_window(&mut aligner, guide, &window[..guide_len]) {
                    println!("Hit in {} at position {}:", record.id(), i);
                    println!("Guide:     {}", String::from_utf8_lossy(guide));
                    println!("Target:    {}", String::from_utf8_lossy(&window[..guide_len]));
                    println!("PAM:       {}", String::from_utf8_lossy(&window[guide_len..guide_len+3]));
                    println!("Score:     {}", score);
                    println!("CIGAR:     {}", cigar);
                    println!();
                }
            }
        }
    }
}
