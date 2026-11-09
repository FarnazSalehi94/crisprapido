//! Isolated benchmark using the unchanged production CFD function on saved hits.
#[path = "../../../src/cfd_score.rs"] mod cfd_score;
use std::{env,fs,time::Instant};
fn main() {
 let args: Vec<_>=env::args().collect(); assert_eq!(args.len(),2);
 let input=fs::read_to_string(&args[1]).unwrap();
 let records: Vec<Vec<&str>>=input.lines().map(|s|s.split('\t').collect()).collect();
 assert!(records.iter().all(|r|r.len()==4));
 cfd_score::init_score_matrices("mismatch_scores.txt","pam_scores.txt").unwrap();
 let guide=b"GCCAGGTACAGGTTCCGTCC";
 let t=Instant::now();
 let scores: Vec<_>=records.iter().map(|r|cfd_score::get_cfd_score(guide,r[0].as_bytes(),r[1],r[2])).collect();
 let elapsed=t.elapsed().as_secs_f64();
 for (r,s) in records.iter().zip(&scores) {
  match s {Some(x)=>assert_eq!(format!("{:.4}",x),r[3]),None=>assert_eq!(r[3],"NA")}
 }
 println!("{{\"records\":{},\"scored\":{},\"production_cfd_seconds\":{:.9}}}",records.len(),scores.iter().filter(|x|x.is_some()).count(),elapsed);
}
