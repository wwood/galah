use crate::QualityFinder;
use checkm::GenomeQuality;
use std::collections::HashMap;
#[cfg(target_family = "unix")]
use std::os::unix::fs::symlink;
#[cfg(target_family = "windows")]
use std::os::windows::fs::symlink_file as symlink;
use std::path::Path;
use std::process::Command;

pub struct CheckM2Analyser {
    // Cache for completeness and contamination results
    pub comp_cont_cache: HashMap<String, (f64, f64)>,
    pub database_path: String,
}

impl CheckM2Analyser {
    pub fn new(database_path: String) -> Self {
        Self {
            comp_cont_cache: HashMap::new(),
            database_path,
        }
    }
}

impl QualityFinder for CheckM2Analyser {
    fn prepare_comp_cont(&mut self, genome_paths: &[String], threads: usize, tmp_path: &Path) {
        self.comp_cont_cache = get_comp_cont(genome_paths, threads, tmp_path, &self.database_path);
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
) -> HashMap<String, (f64, f64)> {
    let mut comp_cont_cache = HashMap::new();
    let checkm2_path = tmp_path.join("checkm2");

    let genomes_dir = tmp_path.join("genomes");
    std::fs::create_dir_all(&genomes_dir).expect("Failed to create genomes directory for CheckM2");
    for fasta in genome_paths {
        let abs_fasta =
            std::fs::canonicalize(fasta).expect("Failed to canonicalize genome path for symlink");
        let new_filename = format!(
            "{}.fna",
            Path::new(fasta).file_stem().unwrap().to_string_lossy()
        );
        symlink(&abs_fasta, genomes_dir.join(new_filename))
            .expect("Failed to create symlink for genome to run CheckM2");
    }

    info!("Running CheckM2 on provided genomes...");
    let output = Command::new("checkm2")
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
    comp_cont_cache
}
