#[cfg(test)]
mod broad_stats_diagnostics {
    use super::*;
    use crate::columba::{deduplicate_columba_candidates, load_guide_fasta, parse_columba_sam_file};
    use std::path::PathBuf;
    #[test]
    fn per_query_stats() {
        let root = PathBuf::from(std::env::var("BROAD_PANELS").unwrap());
        for panel in ["synthetic", "real"] {
            let dir = root.join(panel);
            let refs = load_reference_records(&dir.join("ref.fa")).unwrap();
            let queries = load_guide_fasta(&dir.join("guides.fa")).unwrap();
            let mut candidates = parse_columba_sam_file(&dir.join("edit.sam")).unwrap();
            candidates.extend(parse_columba_sam_file(&dir.join("exact.sam")).unwrap());
            let candidates = deduplicate_columba_candidates(candidates);
            for query in queries {
                let chosen: Vec<_> = candidates.iter().filter(|c| c.query_name == query.query_name).collect();
                let owned: Vec<_> = chosen.iter().map(|c| crate::columba::ColumbaCandidate {
                    query_name: c.query_name.clone(), reference_name: c.reference_name.clone(),
                    reference_start: c.reference_start, reverse: c.reverse, cigar: c.cigar.clone(),
                    edit_distance: c.edit_distance, alignment_score: c.alignment_score,
                }).collect();
                let args = Args { reference: dir.join("ref.fa"), guide: Some(query.sequence.clone()), guides: None,
                    pam: "GG".into(), allow_noncanonical_pam: false, max_mismatches: 4, max_bulges: 1,
                    max_bulge_size: 1, min_match_fraction: 0.75, mismatch_scores: PathBuf::new(),
                    pam_scores: PathBuf::new(), window_size: None, threads: Some(1), no_filter: false,
                    columba_sam: None, columba_bin: None, columba_index: None, keep_columba_sam: false };
                let (hits, stats) = verify_columba_candidates_with_stats(&owned, &refs, &Arc::new(query.sequence.as_bytes().to_vec()), &args);
                let eligible: Vec<_> = hits.iter().filter(|h| h.pam_seq.as_deref() == Some("GG")).collect();
                let anchors: std::collections::HashSet<_> = eligible.iter().map(|h| (&h.ref_id, h.strand, if h.strand == '+' { h.end_pos() } else { h.pos })).collect();
                println!("ANCHOR_STATS\t{}\t{}\t{}\t{}\t{}", panel, query.query_name, hits.len(), eligible.len(), anchors.len());
                println!("EXPANSION_STATS\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", panel, query.query_name,
                    stats.original_candidates, stats.primary_accepted, stats.primary_rejected, stats.expansion_triggers,
                    stats.raw_alternatives_generated, stats.unique_alternatives_verified, stats.fallback_hits_accepted,
                    stats.duplicates_removed, hits.len());
            }
        }
    }
}
