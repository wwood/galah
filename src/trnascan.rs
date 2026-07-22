use crate::Domain;
use crate::TrnaFinder;
use flate2::read::GzDecoder;
use std::io::copy as io_copy;
use std::path::{Path, PathBuf};
use std::process::Command;

pub struct TrnascanAnalyser;

impl TrnaFinder for TrnascanAnalyser {
    /// Run tRNAscan-SE in all three modes (B, A, E) and return the best count.
    fn find_trnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> usize {
        get_trnascan_output_for_domains(
            genome_path,
            &[Domain::Bacteria, Domain::Archaea, Domain::Eukaryota],
            tmp_path,
        )
    }

    fn method_name(&self) -> &str {
        "tRNAscan-SE"
    }
}

/// Run tRNAscan-SE in the modes corresponding to `domains` and return the highest tRNA count.
pub fn get_trnascan_output_for_domains(
    genome_path: &str,
    domains: &[Domain],
    tmp_path: &std::path::Path,
) -> usize {
    let mut modes: Vec<&str> = domains.iter().map(|d| d.trnascan_mode()).collect();
    modes.dedup();

    // tRNAscan-SE does not support gzipped input; decompress once up front if needed
    // (reusing the file if another step, e.g. EukCC, already decompressed this genome).
    let decompressed_path;
    let effective_path: &str = if genome_path.ends_with(".gz") {
        let stem1 = Path::new(genome_path).file_stem().unwrap();
        let genome_name = Path::new(stem1)
            .file_stem()
            .unwrap_or(stem1)
            .to_string_lossy()
            .into_owned();
        let dest = tmp_path.join(format!("{genome_name}.fna"));
        if !dest.is_file() {
            let mut decoder = GzDecoder::new(
                std::fs::File::open(genome_path)
                    .unwrap_or_else(|e| panic!("Failed to open {}: {}", genome_path, e)),
            );
            let mut out = std::fs::File::create(&dest)
                .unwrap_or_else(|e| panic!("Failed to create {:?}: {}", dest, e));
            io_copy(&mut decoder, &mut out)
                .unwrap_or_else(|e| panic!("Failed to decompress {}: {}", genome_path, e));
        }
        decompressed_path = dest.to_string_lossy().into_owned();
        &decompressed_path
    } else {
        genome_path
    };

    let mut best = 0;
    for mode in modes {
        let out_path = run_trnascan(effective_path, mode, tmp_path);
        let trnas = count_unique_standard_trnas(out_path.to_str().unwrap());
        if trnas > best {
            best = trnas;
        }
    }
    best
}

pub fn run_trnascan(genome_path: &str, mode: &str, out_dir: &Path) -> PathBuf {
    let genome_name = Path::new(genome_path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let out_path = out_dir.join(format!("{genome_name}.{mode}.trna.out"));
    if out_path.is_file() {
        info!("Using cached tRNAscan-SE output: {:?}", out_path);
        return out_path;
    }
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

/// Parse tRNAscan-SE output and count unique standard tRNA types.
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
