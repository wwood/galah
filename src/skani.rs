use std;
use std::io::BufReader;
use std::io::Write;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use bird_tool_utils::command::finish_command_safely;
use tempfile;

pub struct SkaniPreclusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
    pub small_genomes: bool,
    pub threads: u16,
}

impl PreclusterDistanceFinder for SkaniPreclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        precluster_skani(
            genome_fasta_paths,
            self.threshold,
            self.min_aligned_threshold,
            self.small_genomes,
            self.threads,
        )
    }

    fn distances_contigs(
        &self,
        genome_fasta_paths: &[&str],
        contig_names: &[&str],
    ) -> SortedPairGenomeDistanceCache {
        precluster_skani_contigs(
            genome_fasta_paths,
            self.threshold,
            self.min_aligned_threshold,
            self.small_genomes,
            self.threads,
            contig_names,
        )
    }

    fn distances_with_references(
        &self,
        genome_fasta_paths: &[&str],
        reference_genomes: &[&str],
    ) -> SortedPairGenomeDistanceCache {
        precluster_skani_with_references(
            genome_fasta_paths,
            reference_genomes,
            self.threshold,
            self.min_aligned_threshold,
            self.small_genomes,
            self.threads,
        )
    }

    fn method_name(&self) -> &str {
        "skani"
    }
}

fn precluster_skani(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    // Create a tempfile to list all the fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for fasta in genome_fasta_paths {
        writeln!(tf, "{fasta}").expect("Failed to write genome fasta paths to tempfile for skani");
    }

    // --sparse only outputs non-zero entries in an edge-list output
    // Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
    info!("Running skani to get distances ..");
    let mut cmd = std::process::Command::new("skani");
    cmd.arg("triangle")
        .arg("-t")
        .arg(format!("{threads}"))
        .arg("--sparse")
        .arg("--min-af")
        .arg(format!("{}", min_aligned_threshold * 100.0));

    if small_genomes {
        cmd.arg("--small-genomes");
    }

    cmd.arg("-l")
        .arg(tf.path().to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running skani command: {:?}", &cmd);

    // Parse the distances
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(stdout_reader);

    let mut distances = SortedPairGenomeDistanceCache::new();

    // Order is likely not conserved, so need to keep track
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // Get index of genomes within genome_fasta_paths
                let genome_id1 = genome_fasta_paths
                    .iter()
                    .position(|&x| x == &record[0])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find genome fasta path in genome_fasta_paths: {}",
                            &record[0]
                        )
                    });
                let genome_id2 = genome_fasta_paths
                    .iter()
                    .position(|&x| x == &record[1])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find genome fasta path in genome_fasta_paths: {}",
                            &record[1]
                        )
                    });

                let ani: f32 = record[2].parse().unwrap_or_else(|_| {
                    panic!("Failed to convert skani ANI to float value: {}", &record[2])
                });
                trace!("Found ANI {}", ani);
                if ani >= threshold {
                    trace!("Accepting ANI since it passed threshold");
                    distances.insert((genome_id1, genome_id2), Some(ani))
                }
            }
            Err(e) => {
                error!("Error parsing skani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    finish_command_safely(process, "skani")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for skani");
    debug!("Found skani distances: {:#?}", distances);
    info!("Finished skani triangle.");

    distances
}

fn precluster_skani_contigs(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
    contig_names: &[&str],
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    // Create a tempfile to list all the fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for fasta in genome_fasta_paths {
        writeln!(tf, "{fasta}").expect("Failed to write genome fasta paths to tempfile for skani");
    }

    // --sparse only outputs non-zero entries in an edge-list output
    info!("Running skani to get distances ..");
    let mut cmd = std::process::Command::new("skani");
    cmd.arg("triangle")
        .arg("-i")
        .arg("-t")
        .arg(format!("{threads}"))
        .arg("--sparse")
        .arg("--min-af")
        .arg(format!("{}", min_aligned_threshold * 100.0));

    if small_genomes {
        cmd.arg("--small-genomes");
    }

    cmd.arg("-l")
        .arg(tf.path().to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running skani command: {:?}", &cmd);

    // Parse the distances
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(stdout_reader);

    let mut distances = SortedPairGenomeDistanceCache::new();

    // Returning contig names from Ref_name and Query_name
    // Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // Get index of contigs within contig_names
                let contig_id1 = contig_names
                    .iter()
                    .position(|&x| x == &record[5])
                    .unwrap_or_else(|| {
                        panic!("Failed to find contig name in contig_names: {}", &record[5])
                    });
                let contig_id2 = contig_names
                    .iter()
                    .position(|&x| x == &record[6])
                    .unwrap_or_else(|| {
                        panic!("Failed to find contig name in contig_names: {}", &record[6])
                    });

                let ani: f32 = record[2].parse().unwrap_or_else(|_| {
                    panic!("Failed to convert skani ANI to float value: {}", &record[2])
                });
                trace!("Found ANI {}", ani);
                if ani >= threshold {
                    trace!("Accepting ANI since it passed threshold");
                    distances.insert((contig_id1, contig_id2), Some(ani))
                }
            }
            Err(e) => {
                error!("Error parsing skani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    finish_command_safely(process, "skani")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for skani");
    debug!("Found skani distances: {:#?}", distances);
    info!("Finished skani triangle.");

    distances
}

/// Create preclusters based on reference genomes using skani
// Assumes that both genome sets are already independently dereplicated
fn precluster_skani_with_references(
    combined_genomes: &[&str],
    reference_genomes: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    if small_genomes {
        panic!("Error: skani does not support small genomes with reference genome preclustering");
    }

    // Create a tempfile to list all the reference file paths
    let mut tf_ref = tempfile::Builder::new()
        .prefix("galah-input-reference-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for fasta in reference_genomes {
        writeln!(tf_ref, "{fasta}")
            .expect("Failed to write reference genome fasta paths to tempfile for skani");
    }

    // Create a tempdir to store the reference genome sketches
    let ref_db = tempfile::TempDir::new()
        .expect("Failed to create temporary directory for skani reference genomes");

    info!("Running skani to sketch reference genomes ..");
    let mut cmd_sketch = std::process::Command::new("skani");
    cmd_sketch
        .arg("sketch")
        .arg("-t")
        .arg(format!("{threads}"))
        .arg("-l")
        .arg(tf_ref.path().to_str().unwrap())
        .arg("-o")
        .arg(ref_db.path().join("galah-skani").to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running skani command: {:?}", &cmd_sketch);

    let mut process_sketch = cmd_sketch
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));

    // Wait for sketching to complete and check the result
    process_sketch
        .wait()
        .expect("Failed to wait for skani sketch");

    // Create a tempfile to list all the genome file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for fasta in combined_genomes {
        if reference_genomes.contains(fasta) {
            continue;
        }
        writeln!(tf, "{fasta}").expect("Failed to write genome fasta paths to tempfile for skani");
    }

    info!("Running skani search to get distances ..");
    let mut cmd = std::process::Command::new("skani");
    cmd.arg("search")
        .arg("-t")
        .arg(format!("{threads}"))
        .arg("--min-af")
        .arg(format!("{}", min_aligned_threshold * 100.0))
        .arg("--ql")
        .arg(tf.path().to_str().unwrap())
        .arg("-d")
        .arg(ref_db.path().join("galah-skani").to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running skani command: {:?}", &cmd);

    // Parse the distances
    // Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(stdout_reader);

    let mut distances = SortedPairGenomeDistanceCache::new();

    // Order is likely not conserved, so need to keep track
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // Get index of genomes within combined_genomes
                let genome_id1 = combined_genomes
                    .iter()
                    .position(|&x| x == &record[0])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find genome fasta path in combined_genomes: {}",
                            &record[0]
                        )
                    });
                let genome_id2 = combined_genomes
                    .iter()
                    .position(|&x| x == &record[1])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find genome fasta path in combined_genomes: {}",
                            &record[1]
                        )
                    });

                let ani: f32 = record[2].parse().unwrap_or_else(|_| {
                    panic!("Failed to convert skani ANI to float value: {}", &record[2])
                });
                trace!("Found ANI {}", ani);
                if ani >= threshold {
                    trace!("Accepting ANI since it passed threshold");
                    distances.insert((genome_id1, genome_id2), Some(ani))
                }
            }
            Err(e) => {
                error!("Error parsing skani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    finish_command_safely(process, "skani")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for skani");
    debug!("Found skani distances: {:#?}", distances);
    info!("Finished skani dist to references.");

    distances
}

pub struct SkaniClusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
    pub small_genomes: bool,
}

impl ClusterDistanceFinder for SkaniClusterer {
    fn initialise(&self) {
        assert!(self.threshold > 1.0);
    }

    fn method_name(&self) -> &str {
        "skani"
    }

    fn get_ani_threshold(&self) -> f32 {
        self.threshold
    }

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32> {
        Some(calculate_skani(
            fasta1,
            fasta2,
            self.small_genomes,
            self.min_aligned_threshold,
        ))
    }
}

pub fn calculate_skani(
    fasta1: &str,
    fasta2: &str,
    small_genomes: bool,
    min_aligned_threshold: f32,
) -> f32 {
    // --sparse only outputs non-zero entries in an edge-list output
    // Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
    let mut cmd = std::process::Command::new("skani");
    cmd.arg("dist")
        .arg("--min-af")
        .arg(format!("{}", min_aligned_threshold * 100.0));

    if small_genomes {
        cmd.arg("--small-genomes");
    }

    cmd.arg("-q")
        .arg(fasta1)
        .arg("-r")
        .arg(fasta2)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running skani command: {:?}", &cmd);

    // Parse the distances
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(stdout_reader);

    let mut num_records = 0;
    let mut to_return = 0.0;

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                if num_records > 0 {
                    error!("Unexpectedly found >1 result from skani");
                    std::process::exit(1);
                }

                assert!(record.len() == 7);
                to_return = record[2].parse().unwrap_or_else(|_| {
                    panic!("Failed to convert skani ANI to float value: {}", &record[2])
                });
                num_records += 1;
            }
            Err(e) => {
                error!("Error parsing skani output: {}", e);
                std::process::exit(1);
            }
        }
    }

    debug!("skani of {} against {} was {:?}", fasta1, fasta2, to_return);
    finish_command_safely(process, "skani")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for skani");
    to_return
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    #[should_panic(
        expected = "Error: skani produces inaccurate results with ANI less than 85%. Provided: 80"
    )]
    fn test_precluster_skani_with_low_ani() {
        init();
        precluster_skani(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            80.0,
            0.2,
            false,
            1,
        );
    }

    #[test]
    fn test_precluster_skani_with_valid_ani() {
        init();
        precluster_skani(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            95.0,
            0.2,
            false,
            1,
        );
    }
}
