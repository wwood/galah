use std;
use std::io::BufReader;
use std::io::Write;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::PreclusterDistanceFinder;

use bird_tool_utils::command::finish_command_safely;
use tempfile;

pub struct DashingPreclusterer {
    /// Fraction, not percentage
    pub min_ani: f32,
    pub threads: u16,
}

impl PreclusterDistanceFinder for DashingPreclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        distances(genome_fasta_paths, self.min_ani, self.threads)
    }

    fn distances_contigs(
        &self,
        _genome_fasta_paths: &[&str],
        _contig_names: &[&str],
    ) -> SortedPairGenomeDistanceCache {
        // Dashing doesn't offer high-quality ANI with self-self comparisons, so we can't use it for contig comparisons.
        SortedPairGenomeDistanceCache::new()
    }

    fn method_name(&self) -> &str {
        "dashing"
    }
}

pub fn distances(
    genome_fasta_paths: &[&str],
    min_ani: f32,
    threads: u16,
) -> SortedPairGenomeDistanceCache {
    // Create a tempfile to list all the fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run dashing");

    for fasta in genome_fasta_paths {
        writeln!(tf, "{}", fasta)
            .expect("Failed to write genome fasta paths to tempfile for dashing");
    }

    // Run dashing to get distances
    info!("Running dashing to get approximate distances ..");
    let mut cmd = std::process::Command::new("dashing");
    cmd.arg("cmp")
        .arg("-o")
        .arg("/dev/null")
        .arg("--nthreads")
        .arg(format!("{}", threads))
        .arg("-M")
        // so that the order of the output remains consistent. Could get around
        // this by finding the indices afterwards but eh.
        .arg("--avoid-sorting")
        .arg("-F")
        .arg(tf.path().to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running dashing command: {:?}", &cmd);

    // Parse the distances
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "dashing"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_reader(stdout_reader);

    let mut distances = SortedPairGenomeDistanceCache::new();

    for (genome_id1, record_res) in rdr.records().enumerate() {
        match record_res {
            Ok(record) => {
                assert_eq!(genome_fasta_paths[genome_id1], &record[0]);
                debug!("Found dashing record {:?}", record);
                for genome_id2 in genome_id1 + 1..genome_fasta_paths.len() {
                    let dist: f32 = record[genome_id2 + 1].parse().unwrap_or_else(|_| {
                        panic!(
                            "Failed to convert dashing dist to float value: {}",
                            &record[genome_id2]
                        )
                    });
                    trace!("Found dist {}", dist);
                    if 1.0 - dist >= min_ani {
                        trace!("Accepting ANI since it passed threshold");
                        distances.insert((genome_id1, genome_id2), Some(1.0 - dist))
                    }
                }
            }
            Err(e) => {
                error!("Error parsing dashing output: {}", e);
                std::process::exit(1);
            }
        }
    }
    finish_command_safely(process, "dashing")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for dashing");
    debug!("Found dashing distances: {:#?}", distances);
    info!("Finished dashing genomes against each other.");

    distances
}
