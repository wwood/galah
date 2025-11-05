use crate::TrnaFinder;
use std::path::{Path, PathBuf};
use std::process::Command;

pub struct TrnascanAnalyser;

impl TrnaFinder for TrnascanAnalyser {
    fn find_trnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> usize {
        get_trnascan_output(genome_path, tmp_path)
    }

    fn method_name(&self) -> &str {
        "tRNAscan-SE"
    }
}

/// Given a genome path and temp dir, run trnascan for both modes and return the best mode, hit count, and output path
pub fn get_trnascan_output(genome_path: &str, tmp_path: &std::path::Path) -> usize {
    let mut best_trnas = 0;
    for mode in ["B", "A"] {
        let out_path = run_trnascan(genome_path, mode, tmp_path);
        let trnas = count_unique_standard_trnas(out_path.to_str().unwrap());
        if trnas > best_trnas {
            best_trnas = trnas;
        }
    }
    best_trnas
}

pub fn run_trnascan(genome_path: &str, mode: &str, out_dir: &Path) -> PathBuf {
    let genome_name = Path::new(genome_path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let out_path = out_dir.join(format!("{genome_name}.{mode}.trna.out"));
    let output = Command::new("tRNAscan-SE")
        .args([
            &format!("-{mode}"),
            "-o",
            out_path.to_str().unwrap(),
            genome_path,
            "--thread",
            "1",
        ])
        .output()
        .expect("Failed to run tRNAscan-SE");

    if !output.status.success() {
        info!(
            "tRNAscan-SE run on {} failed with {}.\nstdout:\n{}\nstderr:\n{}",
            genome_path,
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        panic!("tRNAscan-SE did not run successfully");
    }

    out_path
}

/// Parse tRNAscan-SE output and count unique standard tRNA types
pub fn count_unique_standard_trnas(out_path: &str) -> usize {
    let common_trnas = [
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met",
        "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
    ];
    use std::collections::HashSet;
    let mut unique_trnas = HashSet::new();
    let content = std::fs::read_to_string(out_path).unwrap();
    for line in content.lines().skip(3) {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let trna_type = fields[4];
        if common_trnas.contains(&trna_type) {
            unique_trnas.insert(trna_type);
        }
    }
    unique_trnas.len()
}

pub fn parse_trnascan_hits(out_path: &str) -> usize {
    let content = std::fs::read_to_string(out_path).unwrap();
    content
        .lines()
        .skip(3)
        .filter(|l| !l.trim().is_empty())
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::Write;

    #[test]
    fn test_parse_trnascan_hits() {
        let trna_content = "Header1\nHeader2\nHeader3\nchr1\t1\t1\t100\tType1\tCodon1\t\t\t10.0\tNote1\nchr1\t2\t101\t200\tType2\tCodon2\t\t\t20.0\tNote2\n";
        let path = "test_trna.out";
        let mut file = File::create(path).unwrap();
        file.write_all(trna_content.as_bytes()).unwrap();
        let hits = parse_trnascan_hits(path);
        assert_eq!(hits, 2);
        fs::remove_file(path).unwrap();
    }
}
