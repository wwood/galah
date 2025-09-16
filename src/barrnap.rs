use crate::RrnaFinder;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

pub struct BarrnapAnalyser;

impl RrnaFinder for BarrnapAnalyser {
    fn find_rrnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> (usize, usize, usize) {
        get_barrnap_output(genome_path, tmp_path)
    }

    fn method_name(&self) -> &str {
        "Barrnap"
    }
}

/// Given a genome path and temp dir, run barrnap for both kingdoms and return the best kingdom, hit count, and GFF path
pub fn get_barrnap_output(genome_path: &str, tmp_path: &std::path::Path) -> (usize, usize, usize) {
    let mut best_r5s = 0;
    let mut best_r16s = 0;
    let mut best_r23s = 0;
    for kingdom in ["bac", "arc"] {
        let gff_path = run_barrnap(genome_path, kingdom, 1, tmp_path);
        let (r5s, r16s, r23s) = parse_rrna_types(gff_path.to_str().unwrap());
        let total = r5s + r16s + r23s;
        let best_total = best_r5s + best_r16s + best_r23s;
        if total > best_total {
            best_r5s = r5s;
            best_r16s = r16s;
            best_r23s = r23s;
        }
    }
    (best_r5s, best_r16s, best_r23s)
}

pub fn run_barrnap(genome_path: &str, kingdom: &str, threads: usize, out_dir: &Path) -> PathBuf {
    let genome_name = Path::new(genome_path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let gff_path = out_dir.join(format!("{genome_name}.{kingdom}.gff"));
    let output = Command::new("barrnap")
        .args([
            "--kingdom",
            kingdom,
            "--threads",
            &threads.to_string(),
            genome_path,
        ])
        .output()
        .expect("Failed to run barrnap");

    if !output.status.success() {
        info!(
            "Barrnap run on {} failed with {}.\nstdout:\n{}\nstderr:\n{}",
            genome_path,
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        panic!("Barrnap did not run successfully");
    }

    fs::write(&gff_path, &output.stdout).expect("Failed to write barrnap output");
    gff_path
}

/// Parse barrnap GFF file and count 5S, 16S, 23S rRNAs
pub fn parse_rrna_types(gff_path: &str) -> (usize, usize, usize) {
    let content = std::fs::read_to_string(gff_path).unwrap();
    let mut r5s = 0;
    let mut r16s = 0;
    let mut r23s = 0;
    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if let Some(name) = fields[8].split(';').find_map(|kv| kv.strip_prefix("Name=")) {
            match name {
                "5S_rRNA" => r5s += 1,
                "16S_rRNA" => r16s += 1,
                "23S_rRNA" => r23s += 1,
                _ => {}
            }
        }
    }
    (r5s, r16s, r23s)
}

pub fn parse_barrnap_hits(gff_path: &str) -> usize {
    let content = std::fs::read_to_string(gff_path).unwrap();
    content.lines().filter(|l| !l.starts_with('#')).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::{self, File};
    use std::io::Write;

    #[test]
    fn test_parse_barrnap_hits() {
        let gff_content = "##gff-version 3\nchr1\tbarrnap\trRNA\t1\t100\t.\t+\t.\tID=rrna1\nchr1\tbarrnap\trRNA\t200\t300\t.\t-\t.\tID=rrna2\n";
        let path = "test_barrnap.gff";
        let mut file = File::create(path).unwrap();
        file.write_all(gff_content.as_bytes()).unwrap();
        let hits = parse_barrnap_hits(path);
        assert_eq!(hits, 2);
        fs::remove_file(path).unwrap();
    }
}
