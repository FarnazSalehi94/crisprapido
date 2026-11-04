use std::io::{self, BufRead};
use crisprapido::cfd_score::{get_cfd_score,init_score_matrices};
fn main() {
 init_score_matrices("mismatch_scores.txt","pam_scores.txt").expect("matrices");
 for line in io::stdin().lock().lines() {let line=line.unwrap();let parts:Vec<_>=line.split('\t').collect();assert_eq!(parts.len(),4);
  match get_cfd_score(parts[0].as_bytes(),parts[1].as_bytes(),parts[2],parts[3]) {Some(x)=>println!("{:.12}",x),None=>println!("NA")}
 }
}
