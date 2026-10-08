use crate::QualityFinder;
use checkm::GenomeQuality;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::copy as io_copy;
#[cfg(target_family = "unix")]
use std::os::unix::fs::symlink;
#[cfg(target_family = "windows")]
use std::os::windows::fs::symlink_file as symlink;
use std::path::Path;
use std::process::Command;

fn build_checkm2_command() -> Command {
    if let Ok(cmd_str) = std::env::var("GALAH_CHECKM2_CMD") {
        let mut parts = cmd_str.split_whitespace();
        let prog = parts.next().expect("GALAH_CHECKM2_CMD must not be empty");
        let mut cmd = Command::new(prog);
        cmd.args(parts);
        return cmd;
    }
    if crate::pixi_env::is_on_path("checkm2") {
        return Command::new("checkm2");
    }
    info!("checkm2 not found on PATH, falling back to bundled pixi manifest");
    let manifest = crate::pixi_env::galah_manifest_path();
    let mut cmd = Command::new("pixi");
    cmd.args([
        "run",
        "--manifest-path",
        manifest.to_str().unwrap(),
        "-e",
        "checkm2",
        "checkm2",
    ]);
    cmd
}

pub struct CheckM2Analyser {
    // Cache for completeness and contamination results
    pub comp_cont_cache: HashMap<String, (f64, f64)>,
    pub database_path: String,
    pub quality_report_source_path: Option<std::path::PathBuf>,
}

impl CheckM2Analyser {
    pub fn new(database_path: String) -> Self {
        Self {
            comp_cont_cache: HashMap::new(),
            database_path,
            quality_report_source_path: None,
        }
    }

    pub fn copy_quality_report(&self, dest_path: &str) -> Result<(), String> {
        if let Some(src_path) = &self.quality_report_source_path {
            std::fs::copy(src_path, dest_path).map_err(|e| {
                format!("Failed to copy quality report from {src_path:?} to {dest_path}: {e}")
            })?;
            Ok(())
        } else {
            Err("No quality report available to copy (CheckM2 may not have been run)".to_string())
        }
    }
}

impl QualityFinder for CheckM2Analyser {
    fn prepare_comp_cont(&mut self, genome_paths: &[String], threads: usize, tmp_path: &Path) {
        let (cache, quality_report_path) =
            get_comp_cont(genome_paths, threads, tmp_path, &self.database_path);
        self.comp_cont_cache = cache;
        self.quality_report_source_path = Some(quality_report_path);
    }

    fn find_comp_cont(&self, genome_path: &str) -> (f64, f64) {
        self.comp_cont_cache
            .get(genome_path)
            .copied()
            .expect("Genome path not found in CheckM2 results")
    }

    fn method_name(&self) -> &str {
        "CheckM2"
    }
}

fn get_comp_cont(
    genome_paths: &[String],
    threads: usize,
    tmp_path: &Path,
    database_path: &str,
) -> (HashMap<String, (f64, f64)>, std::path::PathBuf) {
    let mut comp_cont_cache = HashMap::new();
    let checkm2_path = tmp_path.join("checkm2");

    let genomes_dir = tmp_path.join("genomes");
    std::fs::create_dir_all(&genomes_dir).expect("Failed to create genomes directory for CheckM2");
    for fasta in genome_paths {
        // Derive a plain .fna name regardless of whether the source is .fna or .fna.gz
        let stem1 = Path::new(fasta).file_stem().unwrap(); // strips .gz (or .fna)
        let stem = if fasta.ends_with(".gz") {
            Path::new(stem1).file_stem().unwrap_or(stem1) // strips .fna from .fna.gz
        } else {
            stem1
        };
        let dest = genomes_dir.join(format!("{}.fna", stem.to_string_lossy()));

        if fasta.ends_with(".gz") {
            let mut decoder = MultiGzDecoder::new(
                std::fs::File::open(fasta)
                    .unwrap_or_else(|e| panic!("Failed to open {}: {}", fasta, e)),
            );
            let mut out = std::fs::File::create(&dest)
                .unwrap_or_else(|e| panic!("Failed to create {:?}: {}", dest, e));
            io_copy(&mut decoder, &mut out)
                .unwrap_or_else(|e| panic!("Failed to decompress {}: {}", fasta, e));
        } else {
            let abs_fasta = std::fs::canonicalize(fasta)
                .expect("Failed to canonicalize genome path for symlink");
            symlink(&abs_fasta, &dest).expect("Failed to create symlink for genome to run CheckM2");
        }
    }

    info!("Running CheckM2 on provided genomes...");
    let output = build_checkm2_command()
        .args([
            "predict",
            "-o",
            &checkm2_path.to_string_lossy(),
            "--threads",
            &threads.to_string(),
            "-i",
            &genomes_dir.to_string_lossy(),
            "--database_path",
            database_path,
        ])
        .output()
        .expect("Failed to run CheckM2");

    if !output.status.success() {
        info!(
            "CheckM2 failed with {}.\nstdout:\n{}\nstderr:\n{}",
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        panic!("CheckM2 did not run successfully");
    }

    // Parse the quality_report.tsv to get genome_names, completeness, and contamination
    // Then populate the comp_cont_cache using genome_paths
    let quality_report_path = checkm2_path.join("quality_report.tsv");
    assert!(
        quality_report_path.is_file(),
        "CheckM2 did not produce quality_report.tsv at expected location: {:?}",
        quality_report_path
    );

    let checkm2_result = checkm::CheckM2QualityReport::read_file_path(
        quality_report_path
            .to_str()
            .expect("Invalid UTF-8 in quality_report_path"),
    )
    .expect("Failed to parse CheckM2 quality report");

    for genome_path in genome_paths {
        let genome_stem = Path::new(genome_path)
            .file_stem()
            .unwrap()
            .to_string_lossy();
        if let Ok(q) = checkm2_result.retrieve_via_fasta_path(genome_path) {
            comp_cont_cache.insert(
                genome_path.clone(),
                (
                    q.completeness() as f64 * 100.0,
                    q.contamination() as f64 * 100.0,
                ),
            );
        } else if let Some((_, q)) = checkm2_result
            .genome_to_quality
            .iter()
            .find(|(k, _)| **k == genome_stem)
        {
            comp_cont_cache.insert(
                genome_path.clone(),
                (
                    q.completeness() as f64 * 100.0,
                    q.contamination() as f64 * 100.0,
                ),
            );
        } else {
            panic!(
                "No CheckM2 quality found for genome {} (stem {})",
                genome_path, genome_stem
            );
        }
    }
    (comp_cont_cache, quality_report_path)
}
