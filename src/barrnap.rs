use crate::Domain;
use crate::RrnaFinder;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

pub struct BarrnapAnalyser;

impl RrnaFinder for BarrnapAnalyser {
    /// Run barrnap for bac, arc, and euk kingdoms and return the best result.
    fn find_rrnas(
        &self,
        genome_path: &str,
        tmp_path: &std::path::Path,
    ) -> (usize, usize, usize, usize, usize, usize) {
        get_barrnap_output_for_domains(
            genome_path,
            &[Domain::Bacteria, Domain::Archaea, Domain::Eukaryota],
            tmp_path,
        )
    }

    fn method_name(&self) -> &str {
        "Barrnap"
    }
}

/// Run barrnap for the kingdoms corresponding to `domains` and return the best result.
///
/// Returns `(r5s, r16s, r23s, r18s, r28s, r58s)`.
pub fn get_barrnap_output_for_domains(
    genome_path: &str,
    domains: &[Domain],
    tmp_path: &std::path::Path,
) -> (usize, usize, usize, usize, usize, usize) {
    let mut kingdoms: Vec<&str> = domains.iter().map(|d| d.barrnap_kingdom()).collect();
    kingdoms.dedup();

    let mut best = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    for kingdom in kingdoms {
        let gff_path = run_barrnap(genome_path, kingdom, 1, tmp_path);
        let result = parse_rrna_types_all(gff_path.to_str().unwrap());
        let total = result.0 + result.1 + result.2 + result.3 + result.4 + result.5;
        let best_total = best.0 + best.1 + best.2 + best.3 + best.4 + best.5;
        if total > best_total {
            best = result;
        }
    }
    best
}

pub fn run_barrnap(genome_path: &str, kingdom: &str, threads: usize, out_dir: &Path) -> PathBuf {
    let genome_name = Path::new(genome_path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let gff_path = out_dir.join(format!("{genome_name}.{kingdom}.gff"));
    if gff_path.is_file() {
        info!("Using cached Barrnap output: {:?}", gff_path);
        return gff_path;
    }
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

/// Parse all rRNA types (prokaryotic and eukaryotic) from a barrnap GFF file.
///
/// Returns `(r5s, r16s, r23s, r18s, r28s, r58s)`.
pub fn parse_rrna_types_all(gff_path: &str) -> (usize, usize, usize, usize, usize, usize) {
    let content = std::fs::read_to_string(gff_path).unwrap();
    let mut r5s = 0;
    let mut r16s = 0;
    let mut r23s = 0;
    let mut r18s = 0;
    let mut r28s = 0;
    let mut r58s = 0;
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
                "18S_rRNA" => r18s += 1,
                "28S_rRNA" => r28s += 1,
                "5.8S_rRNA" => r58s += 1,
                _ => {}
            }
        }
    }
    (r5s, r16s, r23s, r18s, r28s, r58s)
}

/// Parse a barrnap GFF file and count only prokaryotic rRNA types (5S, 16S, 23S).
///
/// Used when reading pre-computed GFF files generated before eukaryotic support was added.
pub fn parse_rrna_types(gff_path: &str) -> (usize, usize, usize) {
    let (r5s, r16s, r23s, _, _, _) = parse_rrna_types_all(gff_path);
    (r5s, r16s, r23s)
}
