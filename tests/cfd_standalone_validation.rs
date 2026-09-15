//! One-off numerical equivalence harness for CRISPRapido and the pinned
//! lindayqlin/CFD_score_calculator implementation.
//! Pinned commit: 6e6eb86ab2c8c323a229f67a81de3eb312b79ec4.
//!
//! Run with:
//! CFD_STANDALONE_DIR=/path/to/CFD_score_calculator/python \
//!   cargo test --test cfd_standalone_validation -- --ignored --nocapture

use crisprapido::cfd_score;
use std::path::PathBuf;
use std::process::Command;

const TOLERANCE: f64 = 1.0e-6;
const PINNED_COMMIT: &str = "6e6eb86ab2c8c323a229f67a81de3eb312b79ec4";
// Scores emitted by the pinned Python implementation for validation_cases(),
// in the same order. The ignored live-oracle test below verifies these values.
const PINNED_EXPECTED: [f64; 42] = [
    1.0,
    0.857142857,
    0.8,
    0.5625,
    0.099999999750,
    0.492063492063,
    1.0,
    0.0,
    0.029017857156,
    0.259259259,
    0.0857142856,
    0.010934744159,
    1.0,
    0.9,
    1.0,
    0.857142857,
    0.7,
    0.5999999999,
    1.0,
    0.733333333,
    0.538461538,
    0.5,
    0.394871794354,
    0.6599999997,
    0.45,
    0.8,
    0.0,
    0.529411764706,
    1.0,
    1.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.005772005719,
    0.0,
    0.416666666667,
    0.186666666739,
    0.477777777778,
];

struct Case {
    name: &'static str,
    spacer: &'static str,
    protospacer: &'static str,
    pam: &'static str,
}

fn case(
    name: &'static str,
    spacer: &'static str,
    protospacer: &'static str,
    pam: &'static str,
) -> Case {
    assert_eq!(spacer.len(), 20, "{name}: spacer must have 20 columns");
    assert_eq!(
        protospacer.len(),
        20,
        "{name}: protospacer must have 20 columns"
    );
    assert_eq!(pam.len(), 2, "{name}: PAM must have two bases");
    Case {
        name,
        spacer,
        protospacer,
        pam,
    }
}

fn validation_cases() -> Vec<Case> {
    vec![
        // Representative synthetic cases. These exercise calculate_cfd()
        // directly; no CIGAR reconstruction occurs in this harness.
        case(
            "synthetic_perfect_GG",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT",
            "GG",
        ),
        case(
            "synthetic_mismatch_pos1",
            "ACGTACGTACGTACGTACGT",
            "CCGTACGTACGTACGTACGT",
            "GG",
        ),
        case(
            "synthetic_mismatch_pos8",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGAACGTACGTACGT",
            "GG",
        ),
        case(
            "synthetic_mismatch_pos20",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGA",
            "GG",
        ),
        case(
            "synthetic_two_mismatches",
            "ACGTACGTACGTACGTACGT",
            "ACATACGTACGTACGTTCGT",
            "GG",
        ),
        case(
            "synthetic_internal_protospacer_gap",
            "ACGTACGTACGTACGTACGT",
            "ACGTAC-TACGTACGTACGT",
            "GG",
        ),
        case(
            "synthetic_leading_protospacer_gap",
            "ACGTACGTACGTACGTACGT",
            "-CGTACGTACGTACGTACGT",
            "GG",
        ),
        case(
            "synthetic_trailing_protospacer_gap",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACG-",
            "GG",
        ),
        case(
            "synthetic_mismatch_plus_protospacer_gap",
            "ACGTACGTACGTACGTACGT",
            "ACGTATGTACGT-CGTACGT",
            "GG",
        ),
        case(
            "synthetic_perfect_AG_PAM",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGTACGTACGTACGT",
            "AG",
        ),
        case(
            "synthetic_mismatch_CG_PAM",
            "ACGTACGTACGTACGTACGT",
            "ACGTACGAACGTACGTACGT",
            "CG",
        ),
        case(
            "synthetic_gap_GC_PAM",
            "ACGTACGTACGTACGTACGT",
            "ACGTAC-TACGTACGTACGT",
            "GC",
        ),
        // Every exact input formerly intercepted by calculate_cfd()'s 18
        // hardcoded compatibility branches. DNA T is intentional; both
        // calculators independently normalize T to U.
        case(
            "hardcoded_01",
            "CTAACAGTTGCTTTTATCAC",
            "TTAACAGTTGCTTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_02",
            "GAAACAGTCGATTTTATCAC",
            "AAAACAGTCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_03",
            "ATCGATCGATCGATCGATCG",
            "TTCGATCGATCGATCGATCG",
            "GG",
        ),
        case(
            "hardcoded_04",
            "ATCGATCGATCGATCGATCG",
            "ATCGATCGAACGATCGATCG",
            "GG",
        ),
        case(
            "hardcoded_05",
            "ATCGATCGATCGATCGATCG",
            "ATCGATCGATCGATCGATCT",
            "GG",
        ),
        case(
            "hardcoded_06",
            "ATCGATCGATCGATCGATCG",
            "TTCGATCGAACGATCGATCT",
            "GG",
        ),
        case(
            "hardcoded_07",
            "-AAACAGTCGATTTTATCAC",
            "GAAACAGTCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_08",
            "GAAACAGTCGATTTTATCAC",
            "GAAACAGGCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_09",
            "GAAACAGTCGATTTTATCAC",
            "GAAACAGTCGATTTTATAAC",
            "GG",
        ),
        case(
            "hardcoded_10",
            "GAAACAGTCGATTTTATCAC",
            "GAAACAGTCGATTTTATCAA",
            "GG",
        ),
        case(
            "hardcoded_11",
            "GAAACAGTCGATTTTATCAC",
            "GAAACAGGCGATTTTATAAC",
            "GG",
        ),
        case(
            "hardcoded_12",
            "GAAACAGTCGATTTTATCAC",
            "AAAACAGGCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_13",
            "GAAACAGTCGATTTTATCAC",
            "AAAACAGTCGATTTTATCAA",
            "GG",
        ),
        case(
            "hardcoded_14",
            "CTAACAGTTGCTTTTATCAC",
            "CTAACAGATGCTTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_15",
            "GAAACAG-CGATTTTATCAC",
            "GAAACAGTCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_16",
            "GAAACAGTCGATTTTATCA-",
            "GAAACAGTCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_17",
            "GAAACAGTCGATTTTATCAC",
            "TAAACAGTCGATTTTATCAC",
            "GG",
        ),
        case(
            "hardcoded_18",
            "-TCGATCGATCGATCGATCG",
            "ATCGATCGATCGATCGATCG",
            "GG",
        ),
        // Real corrected-pilot I records, reconstructed from the reported
        // CIGAR and the oriented HG00097 haplotype sequence. PAM is the two
        // actual oriented genomic bases adjacent to each locus.
        case(
            "pilot_01_plus_internal_I_GT",
            "CATCTTCTTTCACCTGAACG",
            "CATCTTCTCTCAC-TCAACG",
            "GT",
        ),
        case(
            "pilot_02_plus_leading_I_CA",
            "CATCTTCTTTCACCTGAACG",
            "-AACTTCTTTCACCTGAAAG",
            "CA",
        ),
        case(
            "pilot_03_minus_trailing_I_AC",
            "CATCTTCTTTCACCTGAACG",
            "CTTCTCCTTTCACCTGAAC-",
            "AC",
        ),
        case(
            "pilot_04_plus_trailing_I_CT",
            "CATCTTCTTTCACCTGAACG",
            "CACCTTCTTTCACCTGGAC-",
            "CT",
        ),
        case(
            "pilot_05_minus_leading_I_TT",
            "CATCTTCTTTCACCTGAACG",
            "-ATCTTCTTTCTCCTCAACG",
            "TT",
        ),
        case(
            "pilot_06_plus_internal_I_TG",
            "CTCCGGGGAGAACTCCGGCG",
            "CTCCCGGGAGAACT-TGGCG",
            "TG",
        ),
        case(
            "pilot_07_plus_internal_I_TC",
            "TGGAAGTCCACTCCACTCAG",
            "TGGAAGTCCAC-CCGCTCAG",
            "TC",
        ),
        case(
            "pilot_08_plus_leading_I_GC",
            "ATGTCTCATGAACTACCCTG",
            "-TGTCCCATGAACCACCCTG",
            "GC",
        ),
        case(
            "pilot_09_minus_internal_I_AG",
            "GAACCTTAACATCCATTGTG",
            "GAACATTAACATCCAT-GGG",
            "AG",
        ),
        case(
            "pilot_10_plus_internal_I_GG",
            "GAACCTTAACATCCATTGTG",
            "GA-CCTTAACAGCCATTGTA",
            "GG",
        ),
        case(
            "pilot_11_plus_internal_I_GG",
            "CATCTTCTTTCACCTGAACG",
            "TATC-ACTTTCACCTGAACG",
            "GG",
        ),
        case(
            "pilot_12_minus_internal_I_GG",
            "CATCTTCTTTCACCTGAACG",
            "CAT-TTCCTTCACCTGAACT",
            "GG",
        ),
    ]
}

fn standalone_score(python_dir: &PathBuf, case: &Case) -> f64 {
    let script = python_dir.join("cfd_score_calculator.py");
    let output = Command::new("python3")
        .arg(&script)
        .arg(case.spacer)
        .arg(case.protospacer)
        .arg(case.pam)
        .current_dir(python_dir)
        .output()
        .unwrap_or_else(|error| panic!("{}: could not run pinned Python: {error}", case.name));

    assert!(
        output.status.success(),
        "{}: pinned Python failed: {}",
        case.name,
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout)
        .expect("pinned Python output was not UTF-8")
        .trim()
        .parse::<f64>()
        .unwrap_or_else(|error| panic!("{}: invalid pinned Python score: {error}", case.name))
}

fn initialize_crisprapido_scores() {
    cfd_score::init_score_matrices("mismatch_scores.txt", "pam_scores.txt")
        .expect("failed to initialize CRISPRapido CFD matrices");
}

fn assert_score(spacer: &str, protospacer: &str, pam: &str, expected: f64) {
    let actual = cfd_score::calculate_cfd(spacer, protospacer, pam)
        .unwrap_or_else(|error| panic!("CRISPRapido failed: {error}"));
    assert!(
        (actual - expected).abs() <= TOLERANCE,
        "spacer={spacer}, protospacer={protospacer}, PAM={pam}: expected {expected:.12}, got {actual:.12}"
    );
}

#[test]
fn calculate_cfd_matches_pinned_regression_values() {
    initialize_crisprapido_scores();
    let cases = validation_cases();
    assert_eq!(cases.len(), PINNED_EXPECTED.len());

    let mut failures = Vec::new();
    let mut maximum_difference: f64 = 0.0;
    for (case, expected) in cases.iter().zip(PINNED_EXPECTED) {
        let actual = cfd_score::calculate_cfd(case.spacer, case.protospacer, case.pam)
            .unwrap_or_else(|error| panic!("{}: CRISPRapido failed: {error}", case.name));
        let difference = (actual - expected).abs();
        maximum_difference = maximum_difference.max(difference);
        if difference > TOLERANCE {
            failures.push(case.name);
        }
    }

    println!(
        "pinned regression: total={}, passed={}, failed={}, max_abs_diff={:.12}",
        cases.len(),
        cases.len() - failures.len(),
        failures.len(),
        maximum_difference
    );
    assert!(
        failures.is_empty(),
        "{} of {} cases exceeded tolerance {}: {}",
        failures.len(),
        cases.len(),
        TOLERANCE,
        failures.join(", ")
    );
}

#[test]
fn required_gap_pam_and_former_hardcode_regressions() {
    initialize_crisprapido_scores();

    // Leading gaps are neutral at position 1 on either aligned sequence.
    assert_score("ACGTACGTACGTACGTACGT", "-CGTACGTACGTACGTACGT", "GG", 1.0);
    assert_score("-TCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCG", "GG", 1.0);

    assert_score(
        "ACGTACGTACGTACGTACGT",
        "ACGTAC-TACGTACGTACGT",
        "GG",
        0.492063492063,
    );
    assert_score("ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACG-", "GG", 0.0);
    assert_score(
        "ACGTACGTACGTACGTACGT",
        "ACGTATGTACGT-CGTACGT",
        "GG",
        0.029017857156,
    );
    assert_score(
        "ACGTACGTACGTACGTACGT",
        "ACGTACGTACGTACGTACGT",
        "AG",
        0.259259259,
    );

    // Representative inputs that were previously intercepted by hardcoded
    // production returns now exercise the general matrix algorithm.
    assert_score(
        "ATCGATCGATCGATCGATCG",
        "ATCGATCGAACGATCGATCG",
        "GG",
        0.857142857,
    );
    assert_score(
        "GAAACAGTCGATTTTATCA-",
        "GAAACAGTCGATTTTATCAC",
        "GG",
        0.529411764706,
    );
    assert_score("CTAACAGTTGCTTTTATCAC", "CTAACAGATGCTTTTATCAC", "GG", 0.8);
}

#[test]
#[ignore = "requires the pinned standalone Python source and pickle files"]
fn calculate_cfd_matches_pinned_standalone() {
    let python_dir = PathBuf::from(
        std::env::var("CFD_STANDALONE_DIR")
            .expect("set CFD_STANDALONE_DIR to the pinned repository's python directory"),
    );
    initialize_crisprapido_scores();

    println!("standalone_commit\t{PINNED_COMMIT}");
    println!("case\tspacer\tprotospacer\tPAM\tstandalone\tCRISPRapido\tabs_diff\tresult");
    let mut failures = Vec::new();
    let cases = validation_cases();
    assert_eq!(cases.len(), PINNED_EXPECTED.len());
    let mut maximum_difference: f64 = 0.0;
    for (case, expected) in cases.iter().zip(PINNED_EXPECTED) {
        let standalone = standalone_score(&python_dir, case);
        assert!(
            (standalone - expected).abs() <= TOLERANCE,
            "{}: stored pinned score {:.12} differs from live pinned score {:.12}",
            case.name,
            expected,
            standalone
        );
        let crisprapido = cfd_score::calculate_cfd(case.spacer, case.protospacer, case.pam)
            .unwrap_or_else(|error| panic!("{}: CRISPRapido failed: {error}", case.name));
        let difference = (standalone - crisprapido).abs();
        maximum_difference = maximum_difference.max(difference);
        let result = if difference <= TOLERANCE {
            "PASS"
        } else {
            failures.push(case.name);
            "FAIL"
        };
        println!(
            "{}\t{}\t{}\t{}\t{:.12}\t{:.12}\t{:.12}\t{}",
            case.name,
            case.spacer,
            case.protospacer,
            case.pam,
            standalone,
            crisprapido,
            difference,
            result
        );
    }

    println!(
        "summary\ttotal={}\tpassed={}\tfailed={}\tmax_abs_diff={:.12}",
        cases.len(),
        cases.len() - failures.len(),
        failures.len(),
        maximum_difference
    );

    assert!(
        failures.is_empty(),
        "{} of {} cases exceeded tolerance {}: {}",
        failures.len(),
        cases.len(),
        TOLERANCE,
        failures.join(", ")
    );
}
