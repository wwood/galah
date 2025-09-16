use crate::QualityFinder;
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
    // Name    Completeness    Contamination   Completeness_Model_Used Translation_Table_Used  Coding_Density  Contig_N50      Average_Gene_Length     Genome_Size     GC_Content      Total_Coding_Sequences  Total_Contigs   Max_Contig_Length Additional_Notes
    // binchicken_co412.1081   68.37   2.91    Gradient Boost (General Model)  11      0.885   5745    235.3609865470852       355151  0.33    446     75      24150   None
    let quality_report_path = checkm2_path.join("quality_report.tsv");
    assert!(
        quality_report_path.is_file(),
        "CheckM2 did not produce quality_report.tsv at expected location: {:?}",
        quality_report_path
    );
    let quality_report = std::fs::read_to_string(quality_report_path)
        .expect("Failed to read CheckM2 quality report");

    for line in quality_report.lines().skip(1) {
        debug!("CheckM2 report line: {}", line);
        let fields: Vec<&str> = line.split('\t').collect();
        let genome_name = fields[0].to_string();
        let genome_path = genome_paths
            .iter()
            .find(|g| Path::new(g).file_stem().unwrap().to_string_lossy() == genome_name)
            .expect("Failed to find genome path")
            .to_string();
        let completeness = fields[1].parse().unwrap();
        let contamination = fields[2].parse().unwrap();
        comp_cont_cache.insert(genome_path, (completeness, contamination));
    }
    comp_cont_cache
}
