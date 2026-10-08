use std;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use bird_tool_utils::command::finish_command_safely;
use rayon::prelude::*;
use tempfile;

pub struct SkaniPreclusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
    pub small_genomes: bool,
    pub threads: u16,
    pub low_memory: bool,
    pub skip_sanitize_headers: bool,
}

impl PreclusterDistanceFinder for SkaniPreclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        if self.low_memory {
            precluster_skani_lowmem(
                genome_fasta_paths,
                self.threshold,
                self.min_aligned_threshold,
                self.small_genomes,
                self.threads,
                self.skip_sanitize_headers,
            )
        } else {
            precluster_skani(
                genome_fasta_paths,
                self.threshold,
                self.min_aligned_threshold,
                self.small_genomes,
                self.threads,
                self.skip_sanitize_headers,
            )
        }
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
            self.skip_sanitize_headers,
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
            self.skip_sanitize_headers,
        )
    }

    fn method_name(&self) -> &str {
        "skani"
    }
}

/// Copy a FASTA file to a new tempfile, replacing any tab characters in header
/// lines with spaces. Returns the tempfile (kept alive by the caller).
/// skani uses TSV output, so tabs in sequence headers corrupt its output format.
/// Uses needletail to handle both plain and gzip-compressed FASTA files transparently.
fn sanitize_fasta_headers(fasta_path: &str) -> tempfile::TempPath {
    let mut tf = tempfile::Builder::new()
        .prefix("galah-sanitized-fasta")
        .suffix(".fna")
        .tempfile()
        .expect("Failed to create temporary file for sanitized fasta");

    let mut reader = needletail::parse_fastx_file(fasta_path)
        .unwrap_or_else(|e| panic!("Failed to open fasta file {}: {}", fasta_path, e));

    // Buffer writes so each contig doesn't cost 3 raw write() syscalls - with
    // tens of thousands of (often contig-fragmented) genomes to sanitize, the
    // syscall overhead otherwise dominates.
    let mut writer = BufWriter::new(&mut tf);

    while let Some(record) = reader.next() {
        let record =
            record.unwrap_or_else(|e| panic!("Failed to parse record in {}: {}", fasta_path, e));
        let header = std::str::from_utf8(record.id())
            .unwrap_or_else(|e| panic!("Non-UTF8 header in {}: {}", fasta_path, e))
            .replace('\t', " ");
        writeln!(writer, ">{}", header)
            .expect("Failed to write header to sanitized fasta tempfile");
        writer
            .write_all(&record.seq())
            .expect("Failed to write sequence to sanitized fasta tempfile");
        writeln!(writer).expect("Failed to write newline to sanitized fasta tempfile");
    }

    writer
        .flush()
        .expect("Failed to flush sanitized fasta tempfile writer");
    drop(writer);

    // Close the file handle (but keep the file on disk) to avoid exhausting
    // the OS file descriptor limit when sanitizing tens of thousands of genomes.
    tf.flush()
        .expect("Failed to flush sanitized fasta tempfile");
    tf.into_temp_path()
}

/// Check whether any header line in a FASTA file contains a tab character,
/// without paying the cost of rewriting the file. Lets callers skip
/// `sanitize_fasta_headers` (a full parse + rewrite of every contig) for the
/// overwhelmingly common case where no sanitizing is actually needed.
fn fasta_headers_contain_tab(fasta_path: &str) -> bool {
    let mut reader = needletail::parse_fastx_file(fasta_path)
        .unwrap_or_else(|e| panic!("Failed to open fasta file {}: {}", fasta_path, e));

    while let Some(record) = reader.next() {
        let record =
            record.unwrap_or_else(|e| panic!("Failed to parse record in {}: {}", fasta_path, e));
        if record.id().contains(&b'\t') {
            return true;
        }
    }
    false
}

/// A genome fasta path that either needed no sanitizing (so the original path
/// is used directly, avoiding a redundant copy) or was rewritten to a
/// tempfile with sanitized headers.
enum SanitizedFasta<'a> {
    Original(&'a str),
    Rewritten(tempfile::TempPath),
}

impl<'a> SanitizedFasta<'a> {
    fn path(&self) -> &str {
        match self {
            SanitizedFasta::Original(p) => p,
            SanitizedFasta::Rewritten(tp) => tp.to_str().unwrap(),
        }
    }
}

/// Only pay for `sanitize_fasta_headers`'s full parse + rewrite when the
/// fasta file actually has a tab in a header line.
fn sanitize_fasta_headers_if_needed(fasta_path: &str) -> SanitizedFasta<'_> {
    if fasta_headers_contain_tab(fasta_path) {
        SanitizedFasta::Rewritten(sanitize_fasta_headers(fasta_path))
    } else {
        SanitizedFasta::Original(fasta_path)
    }
}

/// As `sanitize_fasta_headers_if_needed`, but when `skip_sanitize_headers` is
/// set, skip even the tab-detection scan and pass the genome path straight
/// through unchanged. This exists for benchmarking against tools (e.g.
/// skDER) which perform no such sanitizing themselves - if any input genome
/// actually has a tab character in a header line, skani's TSV output will be
/// silently corrupted, so this should not be used on untrusted/unchecked
/// genome sets.
fn maybe_sanitize_fasta_headers(
    fasta_path: &str,
    skip_sanitize_headers: bool,
) -> SanitizedFasta<'_> {
    if skip_sanitize_headers {
        SanitizedFasta::Original(fasta_path)
    } else {
        sanitize_fasta_headers_if_needed(fasta_path)
    }
}

fn precluster_skani(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
    skip_sanitize_headers: bool,
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    // Sanitize FASTA headers to remove tabs, which corrupt skani's TSV output.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sanitized: Vec<SanitizedFasta> = genome_fasta_paths
        .par_iter()
        .map(|p| maybe_sanitize_fasta_headers(p, skip_sanitize_headers))
        .collect();

    // Create a tempfile to list all the sanitized fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for sf in &sanitized {
        writeln!(tf, "{}", sf.path())
            .expect("Failed to write sanitized genome fasta paths to tempfile for skani");
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

    // Order is likely not conserved, so need to keep track.
    // Map sanitized tempfile paths back to original genome_fasta_paths indices.
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // Match sanitized path back to original index
                let genome_id1 = sanitized
                    .iter()
                    .position(|sf| sf.path() == &record[0])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized genome path in sanitized list: {}",
                            &record[0]
                        )
                    });
                let genome_id2 = sanitized
                    .iter()
                    .position(|sf| sf.path() == &record[1])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized genome path in sanitized list: {}",
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

/// Create preclusters using skani in low-memory mode by sketching all genomes
/// then searching all genomes against the sketch database.
fn precluster_skani_lowmem(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
    skip_sanitize_headers: bool,
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    if small_genomes {
        panic!("Error: skani does not support small genomes with low-memory preclustering");
    }

    // Sanitize FASTA headers to remove tabs, which corrupt skani's TSV output.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sanitized: Vec<SanitizedFasta> = genome_fasta_paths
        .par_iter()
        .map(|p| maybe_sanitize_fasta_headers(p, skip_sanitize_headers))
        .collect();

    // Create a tempfile to list all the sanitized fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for sf in &sanitized {
        writeln!(tf, "{}", sf.path())
            .expect("Failed to write sanitized genome fasta paths to tempfile for skani");
    }

    // Create a tempdir to store all genome sketches
    let db_dir =
        tempfile::TempDir::new().expect("Failed to create temporary directory for skani sketches");

    info!("Running skani to sketch genomes for low-memory mode ..");
    let mut cmd_sketch = std::process::Command::new("skani");
    cmd_sketch
        .arg("sketch")
        .arg("-t")
        .arg(format!("{threads}"))
        .arg("-l")
        .arg(tf.path().to_str().unwrap())
        .arg("-o")
        .arg(db_dir.path().join("galah-skani").to_str().unwrap())
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null());
    debug!("Running skani command: {:?}", &cmd_sketch);

    let mut process_sketch = cmd_sketch
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));

    // Wait for sketching to complete and check the result
    process_sketch
        .wait()
        .expect("Failed to wait for skani sketch");

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
        .arg(db_dir.path().join("galah-skani").to_str().unwrap())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null());
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
    // Map sanitized tempfile paths back to original genome_fasta_paths indices.
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                if record[0] == record[1] {
                    // Ignore self matches
                    continue;
                }

                // Match sanitized path back to original index
                let genome_id1 = sanitized
                    .iter()
                    .position(|sf| sf.path() == &record[0])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized genome path in sanitized list: {}",
                            &record[0]
                        )
                    });
                let genome_id2 = sanitized
                    .iter()
                    .position(|sf| sf.path() == &record[1])
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized genome path in sanitized list: {}",
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
    debug!("Found skani distances (low-memory): {:#?}", distances);
    info!("Finished skani low-memory search.");

    distances
}

fn precluster_skani_contigs(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
    small_genomes: bool,
    threads: u16,
    contig_names: &[&str],
    skip_sanitize_headers: bool,
) -> SortedPairGenomeDistanceCache {
    if threshold < 85.0 {
        panic!(
            "Error: skani produces inaccurate results with ANI less than 85%. Provided: {}",
            threshold
        );
    }

    // Sanitize FASTA headers to remove tabs, which corrupt skani's TSV output.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sanitized: Vec<SanitizedFasta> = genome_fasta_paths
        .par_iter()
        .map(|p| maybe_sanitize_fasta_headers(p, skip_sanitize_headers))
        .collect();

    // Create a tempfile to list all the sanitized fasta file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for sf in &sanitized {
        writeln!(tf, "{}", sf.path())
            .expect("Failed to write sanitized genome fasta paths to tempfile for skani");
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

    // Returning contig names from Ref_name and Query_name (columns 5 and 6).
    // These are sanitized (tabs replaced with spaces), matching what was written
    // into the sanitized FASTA headers.
    // Ref_file Query_file ANI Align_fraction_ref Align_fraction_query Ref_name Query_name
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // Get index of contigs within contig_names.
                // Note: contig names with tabs will have had them replaced with
                // spaces in the sanitized FASTA, so we match against the
                // sanitized form here.
                let sanitized_ref = record[5].replace('\t', " ");
                let sanitized_query = record[6].replace('\t', " ");

                let contig_id1 = contig_names
                    .iter()
                    .position(|&x| x.replace('\t', " ") == sanitized_ref)
                    .unwrap_or_else(|| {
                        panic!("Failed to find contig name in contig_names: {}", &record[5])
                    });
                let contig_id2 = contig_names
                    .iter()
                    .position(|&x| x.replace('\t', " ") == sanitized_query)
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
    skip_sanitize_headers: bool,
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

    // Sanitize reference FASTA headers to remove tabs, which corrupt skani's TSV output.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sanitized_refs: Vec<SanitizedFasta> = reference_genomes
        .par_iter()
        .map(|p| maybe_sanitize_fasta_headers(p, skip_sanitize_headers))
        .collect();

    // Create a tempfile to list all the reference file paths
    let mut tf_ref = tempfile::Builder::new()
        .prefix("galah-input-reference-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for sf in &sanitized_refs {
        writeln!(tf_ref, "{}", sf.path())
            .expect("Failed to write sanitized reference genome fasta paths to tempfile for skani");
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
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null());
    debug!("Running skani command: {:?}", &cmd_sketch);

    let mut process_sketch = cmd_sketch
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "skani"));

    // Wait for sketching to complete and check the result
    process_sketch
        .wait()
        .expect("Failed to wait for skani sketch");

    // Sanitize non-reference combined genome FASTA headers.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sanitized_combined: Vec<(&&str, SanitizedFasta)> = combined_genomes
        .par_iter()
        .filter(|fasta| !reference_genomes.contains(fasta))
        .map(|p| (p, maybe_sanitize_fasta_headers(p, skip_sanitize_headers)))
        .collect();

    // Create a tempfile to list all the non-reference genome file paths
    let mut tf = tempfile::Builder::new()
        .prefix("galah-input-genomes")
        .suffix(".txt")
        .tempfile()
        .expect("Failed to open temporary file to run skani");

    for (_, sf) in &sanitized_combined {
        writeln!(tf, "{}", sf.path())
            .expect("Failed to write sanitized genome fasta paths to tempfile for skani");
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
        .stderr(std::process::Stdio::null());
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

    // Order is likely not conserved, so need to keep track.
    // Map sanitized paths back to indices in combined_genomes.
    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                debug!("Found skani record {:?}", record);

                // record[0] is a reference (sanitized), record[1] is a query (sanitized)
                let genome_id1 = sanitized_refs
                    .iter()
                    .position(|sf| sf.path() == &record[0])
                    .map(|i| {
                        combined_genomes
                            .iter()
                            .position(|&g| g == reference_genomes[i])
                            .unwrap_or_else(|| {
                                panic!(
                                    "Failed to find reference genome in combined_genomes: {}",
                                    reference_genomes[i]
                                )
                            })
                    })
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized reference path in sanitized_refs: {}",
                            &record[0]
                        )
                    });

                let genome_id2 = sanitized_combined
                    .iter()
                    .position(|(_, sf)| sf.path() == &record[1])
                    .map(|i| {
                        combined_genomes
                            .iter()
                            .position(|&g| g == *sanitized_combined[i].0)
                            .unwrap_or_else(|| {
                                panic!(
                                    "Failed to find genome in combined_genomes: {}",
                                    sanitized_combined[i].0
                                )
                            })
                    })
                    .unwrap_or_else(|| {
                        panic!(
                            "Failed to find sanitized genome path in sanitized_combined: {}",
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
    pub skip_sanitize_headers: bool,
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
            self.skip_sanitize_headers,
        ))
    }
}

pub fn calculate_skani(
    fasta1: &str,
    fasta2: &str,
    small_genomes: bool,
    min_aligned_threshold: f32,
    skip_sanitize_headers: bool,
) -> f32 {
    // Sanitize FASTA headers to remove tabs, which corrupt skani's TSV output.
    // Genomes whose headers have no tabs are used as-is (no copy needed).
    let sf1 = maybe_sanitize_fasta_headers(fasta1, skip_sanitize_headers);
    let sf2 = maybe_sanitize_fasta_headers(fasta2, skip_sanitize_headers);

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
        .arg(sf1.path())
        .arg("-r")
        .arg(sf2.path())
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
            false,
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
            false,
        );
    }

    #[test]
    fn test_precluster_skani_with_tab_in_headers() {
        init();
        // Genomes with tab characters in FASTA headers should not cause parse errors
        precluster_skani(
            &[
                "tests/data/abisko_tabs/73.20120800_S1D.21.fna",
                "tests/data/abisko_tabs/73.20110800_S2M.16.fna",
            ],
            95.0,
            0.2,
            false,
            1,
            false,
        );
    }

    #[test]
    fn test_fasta_headers_contain_tab_false_when_no_tabs() {
        init();
        assert!(!fasta_headers_contain_tab(
            "tests/data/abisko4/73.20120800_S1X.13.fna"
        ));
    }

    #[test]
    fn test_fasta_headers_contain_tab_true_when_tabs_present() {
        init();
        assert!(fasta_headers_contain_tab(
            "tests/data/abisko_tabs/73.20120800_S1D.21.fna"
        ));
    }

    #[test]
    fn test_sanitize_fasta_headers_if_needed_reuses_original_path() {
        init();
        let fasta_path = "tests/data/abisko4/73.20120800_S1X.13.fna";
        match sanitize_fasta_headers_if_needed(fasta_path) {
            SanitizedFasta::Original(p) => assert_eq!(p, fasta_path),
            SanitizedFasta::Rewritten(_) => {
                panic!("Expected original path to be reused when no tabs are present")
            }
        }
    }

    #[test]
    fn test_sanitize_fasta_headers_if_needed_rewrites_when_tabs_present() {
        init();
        let fasta_path = "tests/data/abisko_tabs/73.20120800_S1D.21.fna";
        match sanitize_fasta_headers_if_needed(fasta_path) {
            SanitizedFasta::Original(_) => {
                panic!("Expected a rewritten tempfile when tabs are present in headers")
            }
            SanitizedFasta::Rewritten(_) => {}
        }
    }

    #[test]
    fn test_maybe_sanitize_fasta_headers_skip_flag_bypasses_tab_check() {
        init();
        // Even a genome with tabs in its headers is passed straight through
        // unchanged when skip_sanitize_headers is set.
        let fasta_path = "tests/data/abisko_tabs/73.20120800_S1D.21.fna";
        match maybe_sanitize_fasta_headers(fasta_path, true) {
            SanitizedFasta::Original(p) => assert_eq!(p, fasta_path),
            SanitizedFasta::Rewritten(_) => {
                panic!("Expected original path to be reused when skip_sanitize_headers is set")
            }
        }
    }

    #[test]
    fn test_precluster_skani_with_skip_sanitize_headers() {
        init();
        // With no tabs in headers, skipping sanitizing should give the same
        // result as the normal path.
        let distances = precluster_skani(
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
            true,
        );
        assert!(distances.contains_key(&(0, 1)));
    }
}
