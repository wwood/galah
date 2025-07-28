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
) -> SortedPairGenomeDistanceCache {
    let mut distances = SortedPairGenomeDistanceCache::new();

    for reference_genome in reference_genomes.iter() {
        // Find idx of reference genome within combined genomes
        let ref_idx_in_combined = combined_genomes
            .iter()
            .position(|&x| x == *reference_genome);

        for (genome_idx, target_genome) in combined_genomes.iter().enumerate() {
            // Skip reference genomes
            if reference_genomes.contains(target_genome) {
                continue;
            }

            // Calculate distance between this genome and the reference genome
            trace!(
                "Calculating ANI between {} and {}",
                target_genome,
                reference_genome
            );
            let ani = calculate_skani(
                target_genome,
                reference_genome,
                small_genomes,
                min_aligned_threshold,
            );

            trace!("Found ANI {}", ani);
            if ani >= threshold {
                trace!("Accepting ANI since it passed threshold");
                distances.insert((genome_idx, ref_idx_in_combined.unwrap()), Some(ani));
            }
        }
    }

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
