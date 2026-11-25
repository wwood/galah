use std;
use std::io::Write;
#[cfg(target_family = "unix")]
use std::os::unix::fs::symlink;
#[cfg(target_family = "windows")]
use std::os::windows::fs::symlink_file as symlink;

use crate::fastani::FastaniClusterer;
use crate::finch::FinchPreclusterer;
use crate::skani::SkaniClusterer;
use crate::skani::SkaniPreclusterer;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;
use crate::QualityFinder;
use crate::SortedPairGenomeDistanceCache;
use bird_tool_utils::clap_utils::*;
use bird_tool_utils::clap_utils::{default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Flag, Manual, Opt, Section};
use clap::*;
use needletail::parse_fastx_file;
use rayon::prelude::*;

use crate::genome_info_file;
use crate::genome_stats;

pub enum Preclusterer {
    Finch(FinchPreclusterer),
    Skani(SkaniPreclusterer),
}

impl PreclusterDistanceFinder for Preclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        match self {
            Preclusterer::Finch(f) => f.distances(genome_fasta_paths),
            Preclusterer::Skani(s) => s.distances(genome_fasta_paths),
        }
    }

    fn distances_contigs(
        &self,
        genome_fasta_paths: &[&str],
        contig_names: &[&str],
    ) -> SortedPairGenomeDistanceCache {
        match self {
            Preclusterer::Finch(f) => f.distances_contigs(genome_fasta_paths, contig_names),
            Preclusterer::Skani(s) => s.distances_contigs(genome_fasta_paths, contig_names),
        }
    }

    fn distances_with_references(
        &self,
        genome_fasta_paths: &[&str],
        reference_genomes: &[&str],
    ) -> SortedPairGenomeDistanceCache {
        match self {
            Preclusterer::Finch(f) => {
                f.distances_with_references(genome_fasta_paths, reference_genomes)
            }
            Preclusterer::Skani(s) => {
                s.distances_with_references(genome_fasta_paths, reference_genomes)
            }
        }
    }

    fn method_name(&self) -> &str {
        match self {
            Preclusterer::Finch(f) => f.method_name(),
            Preclusterer::Skani(s) => s.method_name(),
        }
    }
}

pub enum Clusterer {
    Fastani(FastaniClusterer),
    Skani(SkaniClusterer),
}

impl ClusterDistanceFinder for Clusterer {
    fn initialise(&self) {
        match self {
            Clusterer::Fastani(f) => f.initialise(),
            Clusterer::Skani(s) => s.initialise(),
        }
    }

    fn method_name(&self) -> &str {
        match self {
            Clusterer::Fastani(f) => f.method_name(),
            Clusterer::Skani(s) => s.method_name(),
        }
    }

    fn get_ani_threshold(&self) -> f32 {
        match self {
            Clusterer::Fastani(f) => f.get_ani_threshold(),
            Clusterer::Skani(s) => s.get_ani_threshold(),
        }
    }

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32> {
        match self {
            Clusterer::Fastani(f) => f.calculate_ani(fasta1, fasta2),
            Clusterer::Skani(s) => s.calculate_ani(fasta1, fasta2),
        }
    }
}

pub struct GalahClusterer<'a> {
    pub genome_fasta_paths: Vec<&'a str>,
    pub preclusterer: Preclusterer,
    pub clusterer: Clusterer,
    pub cluster_contigs: bool,
    pub contig_names: &'a Option<Vec<&'a str>>,
    pub reference_genomes: Option<Vec<String>>,
    pub dereplicate_input_genomes_first: bool,
}

pub struct GalahClustererCommandDefinition {
    pub dereplication_ani_argument: String,
    pub dereplication_prethreshold_ani_argument: String,
    pub dereplication_quality_formula_argument: String,
    pub dereplication_run_checkm2_argument: String,
    pub dereplication_checkm2_db_path_argument: String,
    pub dereplication_precluster_method_argument: String,
    pub dereplication_cluster_method_argument: String,
    pub dereplication_aligned_fraction_argument: String,
    pub dereplication_small_genomes_argument: String,
    pub dereplication_small_contigs_argument: String,
    pub dereplication_large_contigs_argument: String,
    pub dereplication_fraglen_argument: String,
    pub dereplication_cluster_contigs_argument: String,
    pub dereplication_low_memory_argument: String,
    pub dereplication_skip_sanitize_headers_argument: String,
    pub dereplication_both_min_aligned_fraction_argument: String,
    pub dereplication_reference_genomes_argument: String,
    pub dereplication_reference_genomes_list_argument: String,
    pub dereplication_skip_input_dereplication_argument: String,
    pub dereplication_output_cluster_definition_file: String,
    pub dereplication_output_representative_fasta_directory: String,
    pub dereplication_output_representative_fasta_directory_copy: String,
    pub dereplication_output_representative_list: String,
}

impl GalahClustererCommandDefinition {
    /// Galah's own argument names, e.g. `--ani`.
    pub fn galah() -> Self {
        Self::with_prefix("")
    }

    /// Argument names for embedding galah's clustering arguments in another
    /// tool's CLI, with `prefix` prepended to every clustering and output
    /// argument, e.g. `with_prefix("dereplication-")` gives
    /// `--dereplication-ani`. Help text refers to the prefixed names.
    ///
    /// CheckM2 execution arguments (`--run-checkm2`, `--checkm2-db-path`) are
    /// left unprefixed, as they are not dereplication-specific. Fields are
    /// public so individual names can be overridden afterwards.
    pub fn with_prefix(prefix: &str) -> Self {
        let p = |name: &str| format!("{prefix}{name}");
        GalahClustererCommandDefinition {
            dereplication_ani_argument: p("ani"),
            dereplication_prethreshold_ani_argument: p("precluster-ani"),
            dereplication_quality_formula_argument: p("quality-formula"),
            dereplication_run_checkm2_argument: "run-checkm2".to_string(),
            dereplication_checkm2_db_path_argument: "checkm2-db-path".to_string(),
            dereplication_precluster_method_argument: p("precluster-method"),
            dereplication_cluster_method_argument: p("cluster-method"),
            dereplication_aligned_fraction_argument: p("min-aligned-fraction"),
            dereplication_small_genomes_argument: p("small-genomes"),
            dereplication_small_contigs_argument: p("small-contigs"),
            dereplication_large_contigs_argument: p("large-contigs"),
            dereplication_fraglen_argument: p("fragment-length"),
            dereplication_cluster_contigs_argument: p("cluster-contigs"),
            dereplication_low_memory_argument: p("low-memory"),
            dereplication_skip_sanitize_headers_argument: p("skip-sanitise-headers"),
            dereplication_both_min_aligned_fraction_argument: p("both-min-aligned-fraction"),
            dereplication_reference_genomes_argument: p("reference-genomes"),
            dereplication_reference_genomes_list_argument: p("reference-genomes-list"),
            dereplication_skip_input_dereplication_argument: p("skip-input-dereplication"),
            dereplication_output_cluster_definition_file: p("output-cluster-definition"),
            dereplication_output_representative_fasta_directory: p(
                "output-representative-fasta-directory",
            ),
            dereplication_output_representative_fasta_directory_copy: p(
                "output-representative-fasta-directory-copy",
            ),
            dereplication_output_representative_list: p("output-representative-list"),
        }
    }
}

impl Default for GalahClustererCommandDefinition {
    fn default() -> Self {
        Self::galah()
    }
}

lazy_static! {
    static ref GALAH_COMMAND_DEFINITION: GalahClustererCommandDefinition =
        GalahClustererCommandDefinition::galah();
}

lazy_static! {
    static ref PROGRAM_NAME: String = std::env::current_exe()
        .ok()
        .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
        .and_then(|s| s.into_string().ok())
        .expect("Failed to find running program basename");
}

lazy_static! {
    static ref CLUSTER_HELP: String = format!(
        "
                     {}
              {}

{}

  {} cluster --genome-fasta-directory input_genomes/
    --output-representative-fasta-directory output_directory/

{}

  {} cluster --ani 95 --precluster-ani 90 --precluster-method finch
    --genome-fasta-list genomes.txt
    --output-cluster-definition clusters.tsv

{}

  {} cluster --genome-fasta-list new_genomes.txt
    --reference-genomes-list reference_genomes.txt
    --output-cluster-definition clusters.tsv

{}

  {} cluster --cluster-contigs --small-contigs --genome-fasta-files contigs.fasta
    --output-cluster-definition contig_clusters.tsv

See {} cluster --full-help for further options and further detail.
",
        ansi_term::Colour::Green.paint(format!(
            "{} cluster",
            std::env::current_exe()
                .ok()
                .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
                .and_then(|s| s.into_string().ok())
                .expect("Failed to find running program basename")
        )),
        ansi_term::Colour::Green.paint("Cluster (dereplicate) genomes"),
        ansi_term::Colour::Purple.paint(format!(
            "Example: Dereplicate at {}% (after pre-clustering at {}%) a directory of .fna\n\
            FASTA files and create a new directory of symlinked FASTA files of\n\representatives:",
            crate::DEFAULT_ANI,
            crate::DEFAULT_PRETHRESHOLD_ANI
        )),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        ansi_term::Colour::Purple.paint(
            "Example: Dereplicate a set of genomes with paths specified in genomes.txt at\n\
            95% ANI, after a preclustering at 90% using the MinHash finch method, and\n\
            output the cluster definition to clusters.tsv:"
        ),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        ansi_term::Colour::Purple.paint(
            "Example: Dereplicate a set of new genomes against a set of reference genomes,\n\
            output the cluster definition to clusters.tsv. New genomes are dereplicated\n\
            amongst themselves first, then only the resulting representative(s) are compared\n\
            against the reference genomes, so new genomes do not need to be pre-dereplicated.\n\
            The reference genomes themselves must already be dereplicated, since\n\
            reference-vs-reference comparisons are never made. Uses less memory and time than\n\
            clustering everything together."
        ),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        ansi_term::Colour::Purple.paint(
            "Example: Dereplicate a set of contigs within FASTA files using small-genomes settings\n\
            (recommended for contigs < 20kb). Can be used for virus or plasmid clustering:"
        ),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
    );
}

pub fn add_dereplication_filtering_parameters_to_section(
    section: Section,
    definition: &GalahClustererCommandDefinition,
) -> Section {
    section
        .option(Opt::new("PATH").long("--checkm2-quality-report").help(&format!(
            "CheckM version 2 quality_report.tsv (i.e. the {} in the output directory output of {}) for defining genome quality, \
            which is used both for filtering and to rank genomes during clustering.",
            monospace_roff("quality_report.tsv"),
            monospace_roff("checkm2 predict ..")
        )))
        .option(Opt::new("PATH").long("--checkm-tab-table").help(&format!(
            "CheckM tab table (i.e. the output of {}). The information contained is used like {}.",
            monospace_roff("checkm .. --tab_table -f PATH .."),
            &monospace_roff("--checkm2-quality-report")
        )))
        .option(Opt::new("PATH").long("--genome-info").help(&format!(
            "dRep style genome info table for defining \
        quality. The information contained is used like {}.",
            &monospace_roff("--checkm2-quality-report")
        )))
        .option(Opt::new("FLOAT").long("--min-completeness").help(
            "Ignore genomes with less completeness than \
        this percentage. [default: not set]",
        ))
        .option(Opt::new("FLOAT").long("--max-contamination").help(
            "Ignore genomes with more contamination than \
        this percentage. [default: not set]",
        ))
        .flag(
            Flag::new()
                .long(&format!("--{}", definition.dereplication_run_checkm2_argument))
                .help(&format!(
                    "Run CheckM2 to generate quality scoring used for clustering. Requires \
                    --{} or CHECKM2DB env variable to be set.",
                    definition.dereplication_checkm2_db_path_argument
                )),
        )
        .option(
            Opt::new("DB_PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_checkm2_db_path_argument
                ))
                .help(
                    "Path to CheckM2 database (required for running CheckM2) \
                    [default: from CHECKM2DB environment variable]",
                ),
        )
}

pub fn add_dereplication_clustering_parameters_to_section(
    section: Section,
    definition: &GalahClustererCommandDefinition,
) -> Section {
    section
        .option(
            Opt::new("FLOAT")
                .long(&format!("--{}", definition.dereplication_ani_argument))
                .help(&format!("Overall ANI level to dereplicate at with the primary clusterer. {}",
                    &default_roff(crate::DEFAULT_ANI))),
        )
        .option(
            Opt::new("FLOAT")
                .long(&format!(
                    "--{}",
                    definition.dereplication_aligned_fraction_argument
                ))
                .help(&format!(
                    "Min aligned fraction of two genomes for \
                clustering. Specify as a percentage (1-100) or a fraction (0-1). {}",
                    default_roff(crate::DEFAULT_ALIGNED_FRACTION)
                )),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_small_genomes_argument
                ))
                .help("Use small-genomes settings in skani calculation. Recommended for sequences < 20kb."),
        )
        .option(
            Opt::new("FLOAT")
                .long(&format!("--{}", definition.dereplication_fraglen_argument))
                .help(&format!(
                    "Length of fragment used in FastANI calculation \
                (i.e. {}). {}",
                &monospace_roff("--fragLen"),
                    default_roff(crate::DEFAULT_FRAGMENT_LENGTH)
                )),
        )
        .option(
            Opt::new("FORMULA")
                .long(&format!(
                    "--{}",
                    definition.dereplication_quality_formula_argument
                ))
                .help(
                    &format!("Scoring function for genome quality {}. One of: {}",
                    default_roff(crate::DEFAULT_QUALITY_FORMULA),
                    bird_tool_utils::clap_utils::table_roff(&[
                        &["formula","description"],
                        &[&monospace_roff("Parks2020_reduced"),
                            &format!("(default) A quality formula described in Parks et. al. 2020 \
                            https://doi.org/10.1038/s41587-020-0501-8 (Supplementary Table 19) but only including those scoring criteria that can be calculated from the sequence without homology searching: {}", &monospace_roff("completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000"))],
                        &[&monospace_roff(crate::HQ_PREFERRED_QUALITY_FORMULA),
                            &format!("As {}, but genomes assessed as MIMAG High quality are always ranked above those that are not, so a Medium quality genome (e.g. one missing a 16S rRNA gene) is never preferred over a High quality one. Only available in {}, since it requires the MIMAG quality assessment from the analyse step.",
                                &monospace_roff("Parks2020_reduced"), &monospace_roff("galah process"))],
                        &[&monospace_roff("completeness-4contamination"),&monospace_roff("completeness-4*contamination")],
                        &[&monospace_roff("completeness-5contamination"),&monospace_roff("completeness-5*contamination")],
                        &[&monospace_roff("dRep"),&monospace_roff("completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log10(N50)")],
                    ]),
                )),
        )
        .option(
            Opt::new("FLOAT")
                .long(&format!(
                    "--{}",
                    definition.dereplication_prethreshold_ani_argument
                ))
                .help(
                    &format!("Require at least this precluster-derived ANI \
                for preclustering and to avoid primary clustering on \
                distant lineages within preclusters. {}",
                    default_roff(crate::DEFAULT_PRETHRESHOLD_ANI)
                )),
        )
        .option(
            Opt::new("NAME")
                .long(&format!(
                    "--{}",
                    definition.dereplication_precluster_method_argument
                ))
                .help(&format!(
                    "method of calculating rough ANI for \
                dereplication. '{}' for finch MinHash, '{}' for Skani. {}",
                monospace_roff("finch"),
                monospace_roff("skani"),
                default_roff(crate::DEFAULT_PRECLUSTER_METHOD)
                )),
        )
        .option(
            Opt::new("NAME")
                .long(&format!(
                    "--{}",
                    definition.dereplication_cluster_method_argument
                ))
                .help(&format!(
                    "method of calculating ANI. \
                '{}' for FastANI, \
                '{}' for Skani. {}",
                monospace_roff("fastani"),
                monospace_roff("skani"),
                default_roff(crate::DEFAULT_CLUSTER_METHOD)
                )),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_cluster_contigs_argument
                ))
                .help(&format!("Cluster contigs within a fasta file instead of genomes. \
                When used, either --{} or --{} must be specified.",
                    definition.dereplication_small_contigs_argument,
                    definition.dereplication_large_contigs_argument)),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_small_contigs_argument
                ))
                .help(&format!("Use small-genomes settings in skani when clustering contigs. \
                Recommended for contigs < 20kb. Mutually exclusive with --{}.",
                    definition.dereplication_large_contigs_argument)),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_large_contigs_argument
                ))
                .help(&format!("Do not use small-genomes settings in skani when clustering contigs. \
                Recommended for contigs >= 20kb. Mutually exclusive with --{}.",
                    definition.dereplication_small_contigs_argument)),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_low_memory_argument
                ))
                .help("Reduce memory use by sketching to file and searching it instead."),
        )
        .flag(
            Flag::new()
                .long(&format!("--{}", definition.dereplication_both_min_aligned_fraction_argument))
                .help(&format!(
                    "Require --{} to be met by both genomes (skani --both-min-af) rather than either. \
                    Useful for e.g. plasmids, to avoid clustering a small sequence with a much larger one that contains it. \
                    Requires the skani cluster method.",
                    definition.dereplication_aligned_fraction_argument
                )),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_skip_sanitize_headers_argument
                ))
                .help(
                    "Skip checking/rewriting FASTA headers that contain tab characters before \
                    running skani, passing genome paths straight through unchanged. Mainly \
                    useful for benchmarking against tools which do not perform \
                    this sanitizing step. If any input genome actually has a tab character in \
                    a header line, skani's TSV output will be silently corrupted, so only use \
                    this on genome sets already known not to have tabs in their headers.",
                ),
        )
        .option(
            Opt::new("PATH ...")
                .long(&format!(
                    "--{}",
                    definition.dereplication_reference_genomes_argument
                ))
                .help(&format!(
                    "Reference genomes to cluster against. {} Mutually exclusive with --{}.",
                    reference_genomes_common_help(),
                    definition.dereplication_reference_genomes_list_argument
                )),
        )
        .option(
            Opt::new("PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_reference_genomes_list_argument
                ))
                .help(&format!(
                    "File containing paths to reference genomes (one per line). {} Mutually exclusive with --{}.",
                    reference_genomes_common_help(),
                    definition.dereplication_reference_genomes_argument
                )),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_skip_input_dereplication_argument
                ))
                .help(&format!(
                    "When used with --{}/--{}, skip dereplicating \
                    input genomes amongst themselves first: compare every input genome directly against \
                    the reference set instead, so near-duplicate input genomes are not merged before \
                    matching (restores the behaviour prior to this flag's introduction). Has no effect \
                    unless reference genomes are given.",
                    definition.dereplication_reference_genomes_argument,
                    definition.dereplication_reference_genomes_list_argument
                )),
        )
}

/// Help text shared by the reference genome options.
fn reference_genomes_common_help() -> &'static str {
    "These should be pre-clustered at the chosen %ANI - \
    reference-vs-reference comparisons are never made, so this is not checked. Input genomes \
    (--genome-fasta-files etc.) are dereplicated amongst themselves first, as if no reference \
    genomes were given at all, and only the resulting representative(s) are then compared against \
    these reference genomes - so input genomes do not need to be pre-dereplicated beforehand. \
    If quality is provided for representative selection, values for these genomes must also be provided. \
    Genomes within the precluster ANI cutoff of each reference will be placed in the same precluster."
}

pub fn add_dereplication_output_parameters_to_section(
    section: Section,
    definition: &GalahClustererCommandDefinition,
) -> Section {
    section
        .option(
            Opt::new("PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_output_cluster_definition_file
                ))
                .help("Output a file of representative<TAB>member lines. Genomes that do not cluster with any other genome are listed as their own representative."),
        )
        .option(
            Opt::new("PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_output_representative_fasta_directory
                ))
                .help("Symlink representative genomes into this directory."),
        )
        .option(
            Opt::new("PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_output_representative_fasta_directory_copy
                ))
                .help("Copy representative genomes into this directory."),
        )
        .option(
            Opt::new("PATH")
                .long(&format!(
                    "--{}",
                    definition.dereplication_output_representative_list
                ))
                .help(
                    "Print newline separated list of paths to representatives \
                        into this file.",
                ),
        )
}

pub struct GalahOutput {
    pub output_clusters_file: Option<std::fs::File>,
    pub output_representative_fasta_directory: Option<std::path::PathBuf>,
    pub output_representative_fasta_directory_copy: Option<std::path::PathBuf>,
    pub output_representative_list: Option<std::fs::File>,
}

pub fn setup_galah_outputs(
    m: &clap::ArgMatches,
    command_definition: &GalahClustererCommandDefinition,
) -> GalahOutput {
    let output_clusters_file = m
        .get_one::<String>(&command_definition.dereplication_output_cluster_definition_file)
        .map(|o| std::fs::File::create(o).expect("Failed to open output cluster definition file"));

    let output_representative_fasta_directory = setup_representative_output_directory(
        m,
        &command_definition.dereplication_output_representative_fasta_directory,
    );
    let output_representative_fasta_directory_copy = setup_representative_output_directory(
        m,
        &command_definition.dereplication_output_representative_fasta_directory_copy,
    );

    let output_representative_list = m
        .get_one::<String>(&command_definition.dereplication_output_representative_list)
        .map(|o| std::fs::File::create(o).expect("Failed to open output representative list file"));

    GalahOutput {
        output_clusters_file,
        output_representative_fasta_directory,
        output_representative_fasta_directory_copy,
        output_representative_list,
    }
}

pub fn run_cluster_subcommand(
    matches: &clap::ArgMatches,
    program_basename: &str,
    program_version: &str,
) {
    let m = matches.subcommand_matches("cluster").unwrap();
    set_log_level(m, true, program_basename, program_version);
    bird_tool_utils::clap_utils::print_full_help_if_needed(
        m,
        cluster_full_help(program_basename, program_version),
    );

    let num_threads = *m.get_one::<u16>("threads").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build_global()
        .expect("Programming error: rayon initialised multiple times");

    let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true)
        .unwrap()
        .into_iter()
        .map(|s| s.split('\t').next().unwrap_or("").to_string())
        .collect();

    if m.get_one::<String>(&GALAH_COMMAND_DEFINITION.dereplication_quality_formula_argument)
        .map(|s| s.as_str())
        == Some(crate::HQ_PREFERRED_QUALITY_FORMULA)
    {
        eprintln!(
            "Error: {}",
            hq_preferred_unavailable_message(&GALAH_COMMAND_DEFINITION)
        );
        std::process::exit(1);
    }

    let cluster_contigs =
        m.get_flag(&GALAH_COMMAND_DEFINITION.dereplication_cluster_contigs_argument);

    // Validate that when --cluster-contigs is used, exactly one of --small-contigs or --large-contigs is specified
    if cluster_contigs {
        let small_contigs =
            m.get_flag(&GALAH_COMMAND_DEFINITION.dereplication_small_contigs_argument);
        let large_contigs =
            m.get_flag(&GALAH_COMMAND_DEFINITION.dereplication_large_contigs_argument);

        match (small_contigs, large_contigs) {
            (false, false) => {
                eprintln!("Error: When --cluster-contigs is used, either --small-contigs or --large-contigs must be specified.");
                eprintln!(
                    "Use --small-contigs for contigs < 20kb, --large-contigs for contigs >= 20kb."
                );
                std::process::exit(1);
            }
            (true, true) => {
                eprintln!("Error: Cannot specify both --small-contigs and --large-contigs.");
                std::process::exit(1);
            }
            _ => {} // Valid: exactly one is specified
        }
    }

    let contig_names_owned = if cluster_contigs {
        if m.contains_id(
            &GALAH_COMMAND_DEFINITION.dereplication_output_representative_fasta_directory,
        ) || m.contains_id(
            &GALAH_COMMAND_DEFINITION.dereplication_output_representative_fasta_directory_copy,
        ) {
            panic!("Cannot specify --cluster-contigs with --output-representative-fasta-directory or --output-representative-fasta-directory-copy");
        }

        // Get all contig names from the fasta files
        let mut contig_names2 = vec![];
        for genome_fasta_path in &genome_fasta_files {
            let mut reader = parse_fastx_file(genome_fasta_path)
                .unwrap_or_else(|_| panic!("Failed to read contig names from file '{}' as there was a problem opening the file",
                    genome_fasta_path));
            while let Some(seq1) = reader.next() {
                let seq = seq1.expect("invalid record");
                let seq_id = {
                    let raw = String::from_utf8_lossy(seq.id());
                    raw.split('\t').next().unwrap_or("").to_string()
                };
                if contig_names2.contains(&seq_id) {
                    panic!(
                        "Duplicate contig name found in file '{}': {}",
                        genome_fasta_path, seq_id
                    );
                }
                contig_names2.push(seq_id);
            }
        }

        Some(contig_names2)
    } else {
        None
    };

    let contig_names = contig_names_owned
        .as_ref()
        .map(|c| c.iter().map(String::as_str).collect::<Vec<&str>>());

    // Handle reference genomes if provided
    let reference_genomes_owned = if let Some(refs) =
        m.get_many::<String>(&GALAH_COMMAND_DEFINITION.dereplication_reference_genomes_argument)
    {
        Some(refs.cloned().collect::<Vec<String>>())
    } else if let Some(ref_file) =
        m.get_one::<String>(&GALAH_COMMAND_DEFINITION.dereplication_reference_genomes_list_argument)
    {
        let content = std::fs::read_to_string(ref_file)
            .unwrap_or_else(|_| panic!("Failed to read reference genomes list file: {}", ref_file));
        Some(
            content
                .lines()
                .filter(|line| !line.trim().is_empty())
                .map(|s| s.split('\t').next().unwrap_or("").to_string())
                .collect::<Vec<String>>(),
        )
    } else {
        None
    };

    let reference_genomes = reference_genomes_owned
        .as_ref()
        .map(|refs| refs.iter().map(|s| s.as_str()).collect::<Vec<&str>>());

    if reference_genomes.is_some() {
        let num_reference_genomes = reference_genomes.as_ref().map_or(0, |r| r.len());
        info!(
            "Clustering against {} reference genomes",
            num_reference_genomes
        );
    }

    // Validate that reference genomes are not used with contig clustering
    if reference_genomes.is_some() && cluster_contigs {
        eprintln!(
            "Error: Reference genome clustering is not currently supported with --cluster-contigs"
        );
        std::process::exit(1);
    }

    // Combine input genomes with reference genomes (when available) for quality filtering
    let (combined_genomes, ref_genomes_for_clusterer) =
        if let Some(ref_genomes) = &reference_genomes {
            let mut combined = ref_genomes
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<String>>();
            combined.extend(genome_fasta_files.iter().cloned());
            (combined, Some(ref_genomes.clone()))
        } else {
            (genome_fasta_files, None)
        };

    let galah = generate_galah_clusterer(
        &combined_genomes,
        &contig_names,
        cluster_contigs,
        m,
        &GALAH_COMMAND_DEFINITION,
        ref_genomes_for_clusterer.as_deref(),
        None,
    )
    .expect("Failed to parse galah clustering arguments correctly");

    // Open file handles here so errors are caught before CPU-heavy clustering
    let output_definitions = setup_galah_outputs(m, &GALAH_COMMAND_DEFINITION);

    let passed_genomes = &galah.genome_fasta_paths;
    info!("Clustering {} genomes ..", passed_genomes.len());
    let clusters = galah.cluster();

    info!("Found {} genome clusters", clusters.len());

    write_galah_outputs(
        output_definitions,
        &clusters,
        passed_genomes,
        contig_names.as_ref(),
    );
    info!("Finished printing genome clusters");
}

pub fn write_galah_outputs(
    output_definitions: GalahOutput,
    clusters: &Vec<Vec<usize>>,
    passed_genomes: &[&str],
    contig_names: Option<&Vec<&str>>,
) {
    let references = match contig_names {
        Some(c) => c,
        None => passed_genomes,
    };
    if let Some(mut f) = output_definitions.output_clusters_file {
        for cluster in clusters {
            let rep_index = cluster[0];
            for genome_index in cluster {
                writeln!(
                    f,
                    "{}\t{}",
                    references[rep_index], references[*genome_index]
                )
                .expect("Failed to write to output clusters file");
            }
        }
    }

    write_cluster_reps_to_directory(
        clusters,
        references,
        &output_definitions.output_representative_fasta_directory,
        &|link, current_stab, rep| {
            symlink(link, std::path::Path::new(&current_stab)).unwrap_or_else(|_| {
                panic!(
                    "Failed to create symbolic link to representative genome {}",
                    rep
                )
            });
        },
    );
    write_cluster_reps_to_directory(
        clusters,
        references,
        &output_definitions.output_representative_fasta_directory_copy,
        &|link, current_stab, rep| {
            std::fs::copy(link, std::path::Path::new(&current_stab)).unwrap_or_else(|_| {
                panic!(
                    "Failed to create symbolic link to representative genome {}",
                    rep
                )
            });
        },
    );

    if let Some(mut f) = output_definitions.output_representative_list {
        for cluster in clusters {
            let rep_index = cluster[0];
            writeln!(f, "{}", references[rep_index])
                .expect("Failed to write to output representative list file");
        }
    }
}

fn setup_representative_output_directory(
    m: &clap::ArgMatches,
    argument: &str,
) -> Option<std::path::PathBuf> {
    match m.get_one::<String>(argument) {
        Some(ref d) => {
            let path = std::path::PathBuf::from(&d);
            if path.exists() {
                if path.is_dir() {
                    if std::fs::read_dir(&path)
                        .unwrap_or_else(|_| panic!("Error opening existing output directory {}", d))
                        .collect::<Vec<_>>()
                        .is_empty()
                    {
                        info!("Using pre-existing but empty {}", argument)
                    } else {
                        error!("The {} specified ({}) exists and is not empty", argument, d);
                        std::process::exit(1)
                    }
                } else {
                    error!(
                        "The {} path specified ({}) exists but is not a directory",
                        argument, d
                    );
                    std::process::exit(1)
                }
            } else {
                info!("Creating {} ..", argument);
                std::fs::create_dir_all(&path)
                    .unwrap_or_else(|_| panic!("Error creating {} ({})", argument, d))
            }
            Some(path)
        }
        None => None,
    }
}

fn write_cluster_reps_to_directory(
    clusters: &Vec<Vec<usize>>,
    passed_genomes: &[&str],
    output_representative_fasta_directory: &Option<std::path::PathBuf>,
    file_creation_fn: &dyn Fn(&std::path::Path, String, &str),
) {
    if let Some(dir) = output_representative_fasta_directory {
        let mut some_names_clashed = false;
        for cluster in clusters {
            let rep = passed_genomes[cluster[0]];
            let link = std::fs::canonicalize(rep).unwrap_or_else(|_| {
                panic!(
                    "Failed to convert representative path into an absolute path: {}",
                    rep
                )
            });

            let basename = std::path::Path::new(rep)
                .file_name()
                .unwrap_or_else(|| panic!("Failed to find file_name from {}", rep));
            // Check that no output files have the same name. If so, add .1.fna, .2.fna etc.
            let mut current_stab = dir.join(basename).to_str().unwrap().to_string();
            let mut counter = 0usize;
            while std::path::Path::new(&current_stab).exists() {
                if !some_names_clashed {
                    warn!("One or more sequence files have the same file name (e.g. ). Renaming clashes by adding .1.fna, .2.fna etc.");
                    some_names_clashed = true;
                }
                counter += 1;
                current_stab = format!("{}.{}.fna", dir.join(basename).to_str().unwrap(), counter);
            }
            file_creation_fn(&link, current_stab, rep);
        }
    };
}

enum CheckMResultEnum {
    GenomeInfoGenomeQuality {
        result: checkm::CheckMResult<genome_info_file::GenomeInfoGenomeQuality>,
    },
    CheckM1Result {
        result: checkm::CheckMResult<checkm::CheckM1GenomeQuality>,
    },
    CheckM2Result {
        result: checkm::CheckMResult<checkm::CheckM2GenomeQuality>,
    },
}

/// A quality report injected by a caller (e.g. `process`) rather than named directly by the
/// user on the command line, distinguished by format since the two are not interchangeable.
pub enum InjectedQualityReport {
    /// Path to a CheckM2 `quality_report.tsv`-format file.
    CheckM2(String),
    /// Path to a genome-info CSV (`genome,completeness,contamination`), as read by
    /// `genome_info_file::read_genome_info_file`.
    GenomeInfo(String),
    /// As `GenomeInfo`, plus the genome paths (as given in the input) assessed as MIMAG High
    /// quality, as required by the `Parks2020_reduced_hqpreferred` quality formula.
    GenomeInfoWithMimag {
        path: String,
        high_quality_genomes: std::collections::HashSet<String>,
    },
}

/// Error message for the `Parks2020_reduced_hqpreferred` formula being used without MIMAG
/// quality available (i.e. outside of `galah process`, or with reference genomes).
pub fn hq_preferred_unavailable_message(
    argument_definition: &GalahClustererCommandDefinition,
) -> String {
    format!(
        "--{} {} requires MIMAG quality from the analyse step, so it can only be used with \
        `galah process` (without reference genomes)",
        argument_definition.dereplication_quality_formula_argument,
        crate::HQ_PREFERRED_QUALITY_FORMULA,
    )
}

pub fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
    injected_quality_report: Option<InjectedQualityReport>,
) -> std::result::Result<Vec<&'a str>, String> {
    let quality_formula = clap_matches
        .get_one::<String>(&argument_definition.dereplication_quality_formula_argument)
        .unwrap()
        .as_str();
    let high_quality_genomes = match &injected_quality_report {
        Some(InjectedQualityReport::GenomeInfoWithMimag {
            high_quality_genomes,
            ..
        }) => Some(high_quality_genomes),
        _ => None,
    };
    if quality_formula == crate::HQ_PREFERRED_QUALITY_FORMULA && high_quality_genomes.is_none() {
        return Err(hq_preferred_unavailable_message(argument_definition));
    }

    if clap_matches.get_flag(&argument_definition.dereplication_cluster_contigs_argument) {
        return Ok(genome_fasta_files.iter().map(|s| &**s).collect());
    }

    match clap_matches.contains_id("checkm-tab-table")
        || clap_matches.contains_id("genome-info")
        || clap_matches.contains_id("checkm2-quality-report")
        || injected_quality_report.is_some()
        || (clap_matches.contains_id(&argument_definition.dereplication_run_checkm2_argument)
            && clap_matches.get_flag(&argument_definition.dereplication_run_checkm2_argument))
    {
        false => {
            warn!("Since CheckM input has not been provided and CheckM2 has been disabled, genomes are not being ordered by quality. Instead the order of their input is being used");
            Ok(genome_fasta_files.iter().map(|s| &**s).collect())
        }
        true => {
            let checkm = if clap_matches.contains_id("checkm-tab-table") {
                info!("Reading CheckM tab table ..");
                CheckMResultEnum::CheckM1Result {
                    result: checkm::CheckM1TabTable::read_file_path(
                        clap_matches.get_one::<String>("checkm-tab-table").unwrap(),
                    )
                    .unwrap(),
                }
            } else if clap_matches.contains_id("checkm2-quality-report") {
                info!("Reading CheckM2 Quality report ..");
                CheckMResultEnum::CheckM2Result {
                    result: checkm::CheckM2QualityReport::read_file_path(
                        clap_matches
                            .get_one::<String>("checkm2-quality-report")
                            .unwrap(),
                    )
                    .unwrap(),
                }
            } else if clap_matches.contains_id("genome-info") {
                if clap_matches
                    .get_one::<String>(&argument_definition.dereplication_quality_formula_argument)
                    .unwrap()
                    == "dRep"
                {
                    return Err(
                        "The dRep quality formula cannot be used with --genome-info".to_string()
                    );
                }
                info!(
                    "Reading genome info file {}",
                    clap_matches.get_one::<String>("genome-info").unwrap()
                );
                CheckMResultEnum::GenomeInfoGenomeQuality {
                    result: genome_info_file::read_genome_info_file(
                        clap_matches.get_one::<String>("genome-info").unwrap(),
                    )
                    .expect("Error parsing genomeInfo file"),
                }
            } else if let Some(report) = &injected_quality_report {
                match report {
                    InjectedQualityReport::CheckM2(path) => {
                        info!("Reading injected CheckM2 Quality report ..");
                        CheckMResultEnum::CheckM2Result {
                            result: checkm::CheckM2QualityReport::read_file_path(path).unwrap(),
                        }
                    }
                    InjectedQualityReport::GenomeInfo(path)
                    | InjectedQualityReport::GenomeInfoWithMimag { path, .. } => {
                        info!("Reading injected combined quality report ..");
                        CheckMResultEnum::GenomeInfoGenomeQuality {
                            result: genome_info_file::read_genome_info_file(path)
                                .expect("Error parsing injected combined quality report"),
                        }
                    }
                }
            } else if clap_matches
                .contains_id(&argument_definition.dereplication_run_checkm2_argument)
                && clap_matches.get_flag(&argument_definition.dereplication_run_checkm2_argument)
            {
                // Run CheckM2 as in analyse
                let db_path = clap_matches.get_one::<String>(&argument_definition.dereplication_checkm2_db_path_argument)
                    .map(|s| s.to_string())
                    .or_else(|| std::env::var("CHECKM2DB").ok())
                    .expect("CheckM2 database path must be provided via --checkm2-db-path or CHECKM2DB env var");
                let threads = *clap_matches
                    .get_one::<u16>("threads")
                    .expect("Failed to parse --threads argument")
                    as usize;
                use crate::checkm2::CheckM2Analyser;
                let tmpdir = tempfile::tempdir().expect("Failed to create tempdir for CheckM2");
                let tmp_path = tmpdir.path();
                let mut analyser = CheckM2Analyser::new(db_path);
                analyser.prepare_comp_cont(genome_fasta_files, threads, tmp_path);
                CheckMResultEnum::CheckM2Result {
                    result: {
                        let quality_report_path =
                            tmp_path.join("checkm2").join("quality_report.tsv");
                        checkm::CheckM2QualityReport::read_file_path(
                            quality_report_path.to_str().unwrap(),
                        )
                        .unwrap()
                    },
                }
            } else {
                panic!("Programming error");
            };

            let max_contamination = parse_percentage(clap_matches, "max-contamination")?;
            let min_completeness = parse_percentage(clap_matches, "min-completeness")?;
            debug!("{:?}", clap_matches);
            debug!(
                "Using max contamination {:?} and min completeness {:?}",
                max_contamination, min_completeness
            );

            let sorted_thresholded_genomes = match quality_formula {
                "completeness-4contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 4*contamination");
                    match &checkm {
                        CheckMResultEnum::CheckM1Result { result } => result
                            .order_fasta_paths_by_completeness_minus_4contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                        CheckMResultEnum::CheckM2Result { result } => result
                            .order_fasta_paths_by_completeness_minus_4contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                        CheckMResultEnum::GenomeInfoGenomeQuality { result } => result
                            .order_fasta_paths_by_completeness_minus_4contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                    }
                }
                "completeness-5contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 5*contamination");
                    match &checkm {
                        CheckMResultEnum::CheckM1Result { result } => result
                            .order_fasta_paths_by_completeness_minus_5contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                        CheckMResultEnum::CheckM2Result { result } => result
                            .order_fasta_paths_by_completeness_minus_5contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                        CheckMResultEnum::GenomeInfoGenomeQuality { result } => result
                            .order_fasta_paths_by_completeness_minus_5contamination(
                                &genome_fasta_files
                                    .iter()
                                    .map(|s| &**s)
                                    .collect::<Vec<&str>>(),
                                min_completeness,
                                max_contamination,
                            )
                            .unwrap(),
                    }
                }
                "Parks2020_reduced" | crate::HQ_PREFERRED_QUALITY_FORMULA => {
                    if quality_formula == crate::HQ_PREFERRED_QUALITY_FORMULA {
                        info!(
                            "Ranking {} MIMAG High quality genome(s) above all others, then by Parks2020_reduced",
                            high_quality_genomes.map_or(0, |hq| hq.len())
                        );
                    }
                    info!("Calculating num_contigs etc. for genome quality assessment ..");
                    let genome_appraisals = match &checkm {
                        CheckMResultEnum::CheckM1Result { result } => {
                            filter_and_calculate_genome_stats(
                                genome_fasta_files,
                                result,
                                min_completeness,
                                max_contamination,
                            )
                        }
                        CheckMResultEnum::CheckM2Result { result } => {
                            filter_and_calculate_genome_stats(
                                genome_fasta_files,
                                result,
                                min_completeness,
                                max_contamination,
                            )
                        }
                        CheckMResultEnum::GenomeInfoGenomeQuality { result } => {
                            filter_and_calculate_genome_stats(
                                genome_fasta_files,
                                result,
                                min_completeness,
                                max_contamination,
                            )
                        }
                    };
                    let mut appraisal: Vec<_> = genome_appraisals
                        .iter()
                        // Calculate quality score
                        .map(|(fasta_file, checkm_quality, stats)| {
                            let score = checkm_quality.completeness as f64 *100.
                                - 5.*checkm_quality.contamination as f64*100.
                                - 5.*stats.num_contigs as f64/100.
                                - 5.*stats.num_ambiguous_bases as f64/100000.;
                            debug!("For genome {} found quality score {}, from checkm {:?} and stats {:?}",
                                fasta_file, score, &checkm_quality, &stats);
                            (fasta_file, score)
                        })
                        .collect();

                    // sort descending
                    appraisal.sort_by(|a, b| {
                        b.1.partial_cmp(&a.1)
                            .expect("Arithmetic error while calculating genome quality")
                    });
                    if quality_formula == crate::HQ_PREFERRED_QUALITY_FORMULA {
                        // Stable sort, so Parks2020_reduced order is kept within each tier.
                        let hq = high_quality_genomes.expect("Checked above");
                        appraisal.sort_by_key(|(f, _)| !hq.contains(f.as_str()));
                    }
                    appraisal
                        .into_iter()
                        .map(|(f, _)| f.as_str())
                        .collect::<Vec<_>>()
                }
                "dRep" => {
                    info!("Calculating num_contigs etc. for genome quality assessment ..");
                    let checkm_result = match &checkm {
                        CheckMResultEnum::CheckM1Result { result } => result,
                        _ => panic!("dRep quality formula only works with CheckM v1 quality scoring since it include strain heterogeneity"),
                    };
                    let genome_appraisals = filter_and_calculate_genome_stats(
                        genome_fasta_files,
                        checkm_result,
                        min_completeness,
                        max_contamination,
                    );
                    let mut appraisal: Vec<_> = genome_appraisals
                        .iter()
                        // Calculate quality score
                        // completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log(N50)
                        // It's log10 specifically, see https://github.com/MrOlm/drep/blob/3ca43f20ec2b43c2c826d2f50e5473d54f4a4510/drep/d_choose.py#L202
                        .map(|(fasta_file, checkm_quality, stats)| {
                            let original_checkm1_quality = checkm_result.retrieve_via_fasta_path(fasta_file).unwrap();
                            let score = checkm_quality.completeness as f64 *100.
                                - 5.*checkm_quality.contamination as f64*100.
                                + checkm_quality.contamination as f64 * original_checkm1_quality.strain_heterogeneity as f64
                                + 0.5*((stats.n50 as f64).log10());
                            debug!("For genome {} found quality score {}, from checkm {:?} and stats {:?}",
                                fasta_file, score, &checkm_quality, &stats);
                            (fasta_file, score)
                        })
                        .collect();

                    // sort descending
                    appraisal.sort_by(|a, b| {
                        b.1.partial_cmp(&a.1)
                            .expect("Arithmetic error while calculating genome quality")
                    });
                    appraisal
                        .into_iter()
                        .map(|(f, _)| f.as_str())
                        .collect::<Vec<_>>()
                }
                _ => panic!("Programming error"),
            };
            info!(
                "Read in genome qualities for {} genomes. {} passed quality thresholds",
                match &checkm {
                    CheckMResultEnum::CheckM1Result { result } => {
                        result.genome_to_quality.len()
                    }
                    CheckMResultEnum::CheckM2Result { result } => {
                        result.genome_to_quality.len()
                    }
                    CheckMResultEnum::GenomeInfoGenomeQuality { result } => {
                        result.genome_to_quality.len()
                    }
                },
                sorted_thresholded_genomes.len()
            );
            if sorted_thresholded_genomes.is_empty() {
                warn!("No genomes passed the quality thresholds, so no clusters will be found");
            }
            Ok(sorted_thresholded_genomes)
        }
    }
}

#[derive(Debug, Copy, Clone)]
struct TwoStatGenomeQuality {
    completeness: f32,
    contamination: f32,
}

fn filter_and_calculate_genome_stats<
    'a,
    T: checkm::GenomeQuality + Clone + Copy + std::fmt::Debug + Send + Sync,
>(
    genome_fasta_files: &'a Vec<String>,
    checkm_stats: &checkm::CheckMResult<T>,
    min_completeness: Option<f32>,
    max_contamination: Option<f32>,
) -> Vec<(
    &'a String,
    TwoStatGenomeQuality,
    genome_stats::GenomeAssemblyStats,
)> {
    genome_fasta_files
        .par_iter()
        // convert to checkm::GenomeQuality
        .map(|fasta_file| {
            (
                fasta_file,
                TwoStatGenomeQuality {
                    completeness: checkm_stats
                        .retrieve_via_fasta_path(fasta_file)
                        .unwrap_or_else(|_| {
                            panic!("Failed to find CheckM statistics for {}", fasta_file)
                        })
                        .completeness(),
                    contamination: checkm_stats
                        .retrieve_via_fasta_path(fasta_file)
                        .unwrap_or_else(|_| {
                            panic!("Failed to find CheckM statistics for {}", fasta_file)
                        })
                        .contamination(),
                },
            )
        })
        // filter out poor checkm quality genomes
        .filter(|(_fasta_file, checkm_quality)| {
            let ok1 = match min_completeness {
                Some(m) => checkm_quality.completeness >= m,
                None => true,
            };
            ok1 && match max_contamination {
                Some(m) => checkm_quality.contamination <= m,
                None => true,
            }
        })
        // calculate stats for good genomes
        .map(|(fasta_file, checkm_quality)| {
            (
                fasta_file,
                checkm_quality,
                crate::genome_stats::calculate_genome_stats(fasta_file),
            )
        })
        .collect()
}

pub fn generate_galah_clusterer<'a>(
    genome_fasta_paths: &'a Vec<String>,
    contig_names: &'a Option<Vec<&'a str>>,
    cluster_contigs: bool,
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
    reference_genomes: Option<&[&str]>,
    injected_quality_report: Option<InjectedQualityReport>,
) -> std::result::Result<GalahClusterer<'a>, String> {
    let precluster_method = clap_matches
        .get_one::<String>(&argument_definition.dereplication_precluster_method_argument)
        .unwrap()
        .as_str();
    let cluster_method = clap_matches
        .get_one::<String>(&argument_definition.dereplication_cluster_method_argument)
        .unwrap()
        .as_str();

    if clap_matches.get_flag(&argument_definition.dereplication_both_min_aligned_fraction_argument)
        && cluster_method != "skani"
    {
        return Err(format!(
            "--{} is only supported when the cluster method is skani (got --{} {})",
            argument_definition.dereplication_both_min_aligned_fraction_argument,
            argument_definition.dereplication_cluster_method_argument,
            cluster_method,
        ));
    }

    // Only require the tools which are actually going to be used e.g. fastANI
    // need not be installed when clustering with skani.
    crate::external_command_checker::check_for_clustering_dependencies(
        precluster_method,
        cluster_method,
    );

    let skip_clusterer = precluster_method == cluster_method;

    match filter_genomes_through_checkm(
        genome_fasta_paths,
        clap_matches,
        argument_definition,
        injected_quality_report,
    ) {
        Err(e) => std::result::Result::Err(e),

        Ok(v2) => {
            let threads = *clap_matches
                .get_one::<u16>("threads")
                .expect("Failed to parse --threads argument");

            let small_genomes = determine_small_genomes_setting(
                clap_matches,
                argument_definition,
                cluster_contigs,
            )?;

            let both_min_af = clap_matches
                .get_flag(&argument_definition.dereplication_both_min_aligned_fraction_argument);
            let skip_sanitize_headers = clap_matches
                .get_flag(&argument_definition.dereplication_skip_sanitize_headers_argument);

            // Filter reference genomes to only include those that passed quality filtering
            let reference_genomes = reference_genomes.map(|refs| {
                let passed_genomes: std::collections::HashSet<&str> = v2.iter().copied().collect();
                refs.iter()
                    .filter(|ref_genome| passed_genomes.contains(**ref_genome))
                    .map(|s| s.to_string())
                    .collect::<Vec<String>>()
            });

            let dereplicate_input_genomes_first = !clap_matches
                .get_flag(&argument_definition.dereplication_skip_input_dereplication_argument);

            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                preclusterer: match clap_matches
                    .get_one::<String>(
                        &argument_definition.dereplication_precluster_method_argument,
                    )
                    .unwrap()
                    .as_str()
                {
                    "finch" => Preclusterer::Finch(FinchPreclusterer {
                        min_ani: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_prethreshold_ani_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse precluster-ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_prethreshold_ani_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse precluster-ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_prethreshold_ani_argument
                                )
                            )
                        }),
                        num_kmers: 1000,
                        kmer_length: 21,
                        low_memory: clap_matches
                            .get_flag(&argument_definition.dereplication_low_memory_argument),
                    }),
                    "skani" => Preclusterer::Skani(SkaniPreclusterer {
                        threshold: {
                            if skip_clusterer {
                                parse_percentage(
                                    clap_matches,
                                    &argument_definition.dereplication_ani_argument,
                                )
                                .unwrap_or_else(|_| {
                                    panic!(
                                        "Failed to parse ani {:?}",
                                        clap_matches.get_one::<f32>(
                                            &argument_definition.dereplication_ani_argument
                                        )
                                    )
                                })
                                .unwrap_or_else(|| {
                                    panic!(
                                        "Failed to parse ani {:?}",
                                        clap_matches.get_one::<f32>(
                                            &argument_definition.dereplication_ani_argument
                                        )
                                    )
                                }) * 100.
                            } else {
                                parse_percentage(
                                    clap_matches,
                                    &argument_definition.dereplication_prethreshold_ani_argument,
                                )
                                .unwrap_or_else(|_| {
                                    panic!(
                                        "Failed to parse precluster-ani {:?}",
                                        clap_matches.get_one::<f32>(
                                            &argument_definition
                                                .dereplication_prethreshold_ani_argument
                                        )
                                    )
                                })
                                .unwrap_or_else(|| {
                                    panic!(
                                        "Failed to parse precluster-ani {:?}",
                                        clap_matches.get_one::<f32>(
                                            &argument_definition
                                                .dereplication_prethreshold_ani_argument
                                        )
                                    )
                                }) * 100.
                            }
                        },
                        min_aligned_threshold: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_aligned_fraction_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        }),
                        small_genomes,
                        threads,
                        low_memory: clap_matches
                            .get_flag(&argument_definition.dereplication_low_memory_argument),
                        both_min_af,
                        skip_sanitize_headers,
                    }),
                    _ => panic!("Programming error"),
                },
                clusterer: match clap_matches
                    .get_one::<String>(&argument_definition.dereplication_cluster_method_argument)
                    .unwrap()
                    .as_str()
                {
                    "fastani" => Clusterer::Fastani(FastaniClusterer {
                        threshold: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_ani_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_ani_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_ani_argument
                                )
                            )
                        }) * 100.,
                        min_aligned_threshold: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_aligned_fraction_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        }),
                        fraglen: *clap_matches
                            .get_one::<u32>(&argument_definition.dereplication_fraglen_argument)
                            .unwrap_or_else(|| {
                                panic!(
                                    "Failed to parse fragment length {:?}",
                                    clap_matches.get_one::<u32>(
                                        &argument_definition.dereplication_fraglen_argument
                                    )
                                )
                            }),
                    }),
                    "skani" => Clusterer::Skani(SkaniClusterer {
                        threshold: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_ani_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_ani_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse ani {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_ani_argument
                                )
                            )
                        }) * 100.,
                        min_aligned_threshold: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_aligned_fraction_argument,
                        )
                        .unwrap_or_else(|_| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        })
                        .unwrap_or_else(|| {
                            panic!(
                                "Failed to parse min-aligned-fraction {:?}",
                                clap_matches.get_one::<f32>(
                                    &argument_definition.dereplication_aligned_fraction_argument
                                )
                            )
                        }),
                        small_genomes,
                        both_min_af,
                        skip_sanitize_headers,
                    }),
                    _ => panic!("Programming error"),
                },
                cluster_contigs,
                contig_names,
                reference_genomes,
                dereplicate_input_genomes_first,
            })
        }
    }
}

pub fn parse_percentage(
    m: &clap::ArgMatches,
    parameter: &str,
) -> std::result::Result<Option<f32>, String> {
    match m.contains_id(parameter) {
        true => {
            let mut percentage: f32 = *m.get_one::<f32>(parameter).unwrap();
            if (1.0..=100.0).contains(&percentage) {
                percentage /= 100.0;
            } else if !(0.0..=100.0).contains(&percentage) {
                error!("Invalid alignment percentage: '{}'", percentage);
                let err = std::result::Result::Err(format!(
                    "Invalid percentage specified for --{parameter}: '{percentage}'",
                ));
                return err;
            }
            debug!("Using {} {}%", parameter, percentage * 100.0);
            Ok(Some(percentage))
        }
        false => Ok(None),
    }
}

impl GalahClusterer<'_> {
    pub fn cluster(&self) -> Vec<Vec<usize>> {
        let reference_genomes_strs: Option<Vec<&str>> = self
            .reference_genomes
            .as_ref()
            .map(|refs| refs.iter().map(|s| s.as_ref()).collect());

        crate::clusterer::cluster(
            &self.genome_fasta_paths,
            &self.preclusterer,
            &self.clusterer,
            self.cluster_contigs,
            self.contig_names.as_deref(),
            reference_genomes_strs.as_deref(),
            self.dereplicate_input_genomes_first,
        )
    }
}

pub fn cluster_full_help(program_basename: &str, program_version: &str) -> Manual {
    let mut manual = Manual::new(&format!("{program_basename} cluster"))
        .about(format!("Cluster genome FASTA files by average nucleotide identity (version {program_version})"))
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description("This cluster mode dereplicates genomes, choosing a subset of the input genomes as representatives. \
            Required inputs are (1) a genome definition, and (2) an output format definition.\n\n\
            The source code for this program can be found at https://github.com/wwood/galah or https://github.com/wwood/coverm")
        .custom_synopsis_expansion("<GENOME_INPUTS> <OUTPUT_ARGUMENTS>");

    // input
    manual = manual.custom(
        bird_tool_utils::clap_utils::add_genome_specification_to_section(Section::new(
            "Genome input",
        )),
    );

    // filtering
    manual = manual.custom(add_dereplication_filtering_parameters_to_section(
        Section::new("Filtering parameters"),
        &GALAH_COMMAND_DEFINITION,
    ));

    // clustering
    manual = manual.custom(add_dereplication_clustering_parameters_to_section(
        Section::new("Clustering parameters"),
        &GALAH_COMMAND_DEFINITION,
    ));

    // output
    manual = manual.custom(add_dereplication_output_parameters_to_section(
        Section::new("Output"),
        &GALAH_COMMAND_DEFINITION,
    ));

    // general
    manual = manual.custom(
        Section::new("General parameters")
            .option(
                Opt::new("INT")
                    .short("-t")
                    .long("--threads")
                    .help(&format!("Number of threads. {}", default_roff("1"))),
            )
            .flag(
                Flag::new()
                    .short("-v")
                    .long("--verbose")
                    .help("Print extra debugging information"),
            )
            .flag(Flag::new().short("-q").long("--quiet").help(
                "Unless there is an error, do not print \
                log messages",
            ))
            .flag(
                Flag::new()
                    .short("-h")
                    .long("--help")
                    .help("Output a short usage message."),
            )
            .flag(
                Flag::new()
                    .long("--full-help")
                    .help("Output a full help message and display in 'man'."),
            )
            .flag(Flag::new().long("--full-help-roff").help(
                "Output a full help message in raw ROFF format for \
            conversion to other formats.",
            )),
    );
    manual
}

/// Add galah's clustering arguments to a clap command, named according to
/// `definition` (use `GalahClustererCommandDefinition::with_prefix` to embed
/// them in another tool's CLI). Includes the CheckM2 execution arguments but
/// not genome quality file/filter arguments (`--checkm2-quality-report`,
/// `--min-completeness` etc.), which are not prefixed and which embedding tools
/// may want to constrain differently.
///
/// Argument ids and long names both come from `definition`, so
/// `ArgMatches` lookups via the definition work for any prefix. `definition`
/// must be `'static` (e.g. held in a `lazy_static`) as clap argument names are
/// `&'static str`.
pub fn add_dereplication_clustering_clap_args(
    command: clap::Command,
    definition: &'static GalahClustererCommandDefinition,
) -> clap::Command {
    let d = definition;
    let f = |name: &str| format!("--{name}");
    command
        .arg(Arg::new(d.dereplication_ani_argument.as_str())
            .long(d.dereplication_ani_argument.as_str())
            .help("Average nucleotide identity threshold for clustering")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ANI))
        .arg(Arg::new(d.dereplication_small_genomes_argument.as_str())
            .long(d.dereplication_small_genomes_argument.as_str())
            .help("Use small-genomes settings in skani calculation. Recommended for sequences < 20kb.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(d.dereplication_fraglen_argument.as_str())
            .long(d.dereplication_fraglen_argument.as_str())
            .help("Length of fragment used in FastANI calculation (i.e. --fragLen)")
            .value_parser(clap::value_parser!(u32))
            .default_value(crate::DEFAULT_FRAGMENT_LENGTH))
        .arg(Arg::new(d.dereplication_aligned_fraction_argument.as_str())
            .long(d.dereplication_aligned_fraction_argument.as_str())
            .help("Min aligned fraction of two genomes for clustering. Specify as a percentage (1-100) or a fraction (0-1)")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ALIGNED_FRACTION))
        .arg(Arg::new(d.dereplication_quality_formula_argument.as_str())
            .long(d.dereplication_quality_formula_argument.as_str())
            .value_parser([
                "completeness-4contamination",
                "completeness-5contamination",
                "Parks2020_reduced",
                crate::HQ_PREFERRED_QUALITY_FORMULA,
                "dRep"])
            .default_value(crate::DEFAULT_QUALITY_FORMULA))
        .arg(Arg::new(d.dereplication_run_checkm2_argument.as_str())
            .long(d.dereplication_run_checkm2_argument.as_str())
            .help("Run CheckM2 for genome quality scoring during clustering")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(d.dereplication_checkm2_db_path_argument.as_str())
            .long(d.dereplication_checkm2_db_path_argument.as_str())
            .help("Path to CheckM2 database. If not specified, will use $CHECKM2DB environment variable if set."))
        .arg(Arg::new(d.dereplication_prethreshold_ani_argument.as_str())
            .long(d.dereplication_prethreshold_ani_argument.as_str())
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_PRETHRESHOLD_ANI))
        .arg(Arg::new(d.dereplication_precluster_method_argument.as_str())
            .long(d.dereplication_precluster_method_argument.as_str())
            .help("method of calculating rough ANI. 'finch' for finch MinHash, 'skani' for skani")
            .value_parser(crate::PRECLUSTER_METHODS)
            .default_value(crate::DEFAULT_PRECLUSTER_METHOD))
        .arg(Arg::new(d.dereplication_cluster_method_argument.as_str())
            .long(d.dereplication_cluster_method_argument.as_str())
            .help("method of calculating ANI. 'fastani' for FastANI, 'skani' for skani")
            .value_parser(crate::CLUSTER_METHODS)
            .default_value(crate::DEFAULT_CLUSTER_METHOD))
        .arg(Arg::new(d.dereplication_cluster_contigs_argument.as_str())
            .long(d.dereplication_cluster_contigs_argument.as_str())
            .help(format!("Cluster contigs within a fasta file instead of genomes. When used, either {} or {} must be specified.",
                f(&d.dereplication_small_contigs_argument), f(&d.dereplication_large_contigs_argument)))
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(d.dereplication_small_contigs_argument.as_str())
            .long(d.dereplication_small_contigs_argument.as_str())
            .help("Use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs < 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(d.dereplication_cluster_contigs_argument.as_str()))
        .arg(Arg::new(d.dereplication_large_contigs_argument.as_str())
            .long(d.dereplication_large_contigs_argument.as_str())
            .help("Do not use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs >= 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(d.dereplication_cluster_contigs_argument.as_str())
            .conflicts_with(d.dereplication_small_contigs_argument.as_str()))
        .arg(Arg::new(d.dereplication_low_memory_argument.as_str())
            .long(d.dereplication_low_memory_argument.as_str())
            .help("Reduce memory by sketching all genomes and searching instead of triangle")
            .action(clap::ArgAction::SetTrue)
            .conflicts_with(d.dereplication_reference_genomes_argument.as_str())
            .conflicts_with(d.dereplication_reference_genomes_list_argument.as_str()))
        .arg(Arg::new(d.dereplication_both_min_aligned_fraction_argument.as_str())
            .long(d.dereplication_both_min_aligned_fraction_argument.as_str())
            .help(format!("Require {} to be met by both genomes (skani --both-min-af) rather than either. Useful for e.g. plasmids, to avoid clustering a small sequence with a much larger one that contains it. Requires the skani cluster method.",
                f(&d.dereplication_aligned_fraction_argument)))
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(d.dereplication_skip_sanitize_headers_argument.as_str())
            .long(d.dereplication_skip_sanitize_headers_argument.as_str())
            .help("Skip checking/rewriting FASTA headers containing tabs before running skani, using genome paths as-is. Mainly for benchmarking against tools (e.g. skDER) which do not sanitize headers. If any input genome has a tab in a header line, skani's TSV output will be silently corrupted, so only use this when input genomes are known not to have tabs in their headers.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(d.dereplication_reference_genomes_argument.as_str())
            .long(d.dereplication_reference_genomes_argument.as_str())
            .help("Reference genomes to cluster against. These should already be dereplicated amongst themselves. Input genomes are dereplicated amongst themselves first, then only the resulting representative(s) are compared against these reference genomes - so input genomes do not need to be pre-dereplicated. Uses less memory and time than clustering everything together.")
            .value_delimiter(' ')
            .num_args(1..)
            .conflicts_with(d.dereplication_low_memory_argument.as_str())
            .conflicts_with(d.dereplication_reference_genomes_list_argument.as_str()))
        .arg(Arg::new(d.dereplication_reference_genomes_list_argument.as_str())
            .long(d.dereplication_reference_genomes_list_argument.as_str())
            .help("File containing paths to reference genomes (one per line). These should already be dereplicated amongst themselves. Input genomes are dereplicated amongst themselves first, then only the resulting representative(s) are compared against these reference genomes - so input genomes do not need to be pre-dereplicated. Uses less memory and time than clustering everything together.")
            .conflicts_with(d.dereplication_low_memory_argument.as_str())
            .conflicts_with(d.dereplication_reference_genomes_argument.as_str()))
        .arg(Arg::new(d.dereplication_skip_input_dereplication_argument.as_str())
            .long(d.dereplication_skip_input_dereplication_argument.as_str())
            .help(format!("When used with {}/{}, skip dereplicating input genomes amongst themselves first: compare every input genome directly against the reference set instead, so near-duplicate input genomes are not merged before matching (restores the behaviour prior to this flag's introduction). Has no effect unless reference genomes are given.",
                f(&d.dereplication_reference_genomes_argument), f(&d.dereplication_reference_genomes_list_argument)))
            .action(clap::ArgAction::SetTrue))
}

/// Add galah's output arguments (cluster definition file, representative
/// directory/copy/list) to a clap command, named according to `definition`.
/// No `required_unless_present` constraints are applied, as the embedding
/// command decides when output is required.
pub fn add_dereplication_output_clap_args(
    command: clap::Command,
    definition: &'static GalahClustererCommandDefinition,
) -> clap::Command {
    let d = definition;
    command
        .arg(Arg::new(d.dereplication_output_cluster_definition_file.as_str())
            .long(d.dereplication_output_cluster_definition_file.as_str())
            .help("Output a file of representative<TAB>member lines. Genomes that do not cluster with any other genome are listed as their own representative"))
        .arg(Arg::new(d.dereplication_output_representative_fasta_directory.as_str())
            .long(d.dereplication_output_representative_fasta_directory.as_str())
            .help("Symlink representative genomes into this directory"))
        .arg(Arg::new(d.dereplication_output_representative_fasta_directory_copy.as_str())
            .long(d.dereplication_output_representative_fasta_directory_copy.as_str())
            .help("Copy representative genomes into this directory"))
        .arg(Arg::new(d.dereplication_output_representative_list.as_str())
            .long(d.dereplication_output_representative_list.as_str())
            .help("Print newline separated list of paths to representatives into this file"))
}

pub fn add_cluster_subcommand(app: clap::Command) -> clap::Command {
    let mut cluster_subcommand = add_clap_verbosity_flags(Command::new("cluster"))
        .about("Cluster FASTA files by average nucleotide identity")
        .override_help(CLUSTER_HELP.as_str())
        .arg(Arg::new("full-help").long("full-help").action(clap::ArgAction::SetTrue))
        .arg(Arg::new("full-help-roff").long("full-help-roff").action(clap::ArgAction::SetTrue))
        .arg(Arg::new("checkm-tab-table")
            .long("checkm-tab-table")
            .help("Output of CheckM lineage_wf/taxonomy_wf/qa with --tab_table specified"))
        .arg(Arg::new("checkm2-quality-report")
            .long("checkm2-quality-report")
            .help("Output of CheckM2 predict"))
        .arg(Arg::new("genome-info")
            .long("genome-info")
            .help("genomeInfo file in same format as dRep i.e. a CSV with three header columns, first line 'genome,completeness,contamination'."))
        .arg(Arg::new("min-completeness")
            .long("min-completeness")
            .help("Genomes with less than this percentage of completeness are excluded")
            .value_parser(clap::value_parser!(f32))
            .default_value("0"))
        .arg(Arg::new("max-contamination")
            .long("max-contamination")
            .help("Genomes with greater than this percentage of contamination are excluded")
            .value_parser(clap::value_parser!(f32))
            .default_value("100"))
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .help("Number of CPU threads to use")
            .default_value("1")
            .value_parser(clap::value_parser!(u16)))
        .arg(Arg::new("output-cluster-definition")
            .short('o')
            .long("output-cluster-definition")
            .help("Output a file of representative<TAB>member lines. Genomes that do not cluster with any other genome are listed as their own representative")
            .required_unless_present_any([
                "output-representative-fasta-directory",
                "output-representative-fasta-directory-copy",
                "output-representative-list",
                "full-help",
                "full-help-roff",]))
        .arg(Arg::new("output-representative-fasta-directory")
            .long("output-representative-fasta-directory")
            .help("Symlink representative genomes into this directory")
            .required_unless_present_any([
                "output-cluster-definition",
                "output-representative-list",
                "output-representative-fasta-directory-copy",
                "full-help",
                "full-help-roff",]))
        .arg(Arg::new("output-representative-fasta-directory-copy")
            .long("output-representative-fasta-directory-copy")
            .help("Copy representative genomes into this directory")
            .required_unless_present_any([
                "output-cluster-definition",
                "output-representative-fasta-directory",
                "output-representative-list",
                "full-help",
                "full-help-roff",]))
        .arg(Arg::new("output-representative-list")
            .long("output-representative-list")
            .help("Print newline separated list of paths to representatives into this directory")
            .required_unless_present_any([
                "output-representative-fasta-directory",
                "output-representative-fasta-directory-copy",
                "output-cluster-definition",
                "full-help",
                "full-help-roff",]));

    cluster_subcommand =
        add_dereplication_clustering_clap_args(cluster_subcommand, &GALAH_COMMAND_DEFINITION);
    cluster_subcommand =
        bird_tool_utils::clap_utils::add_genome_specification_arguments(cluster_subcommand);

    app.subcommand(cluster_subcommand)
}

/// Determine the small_genomes setting based on contig-specific flags
fn determine_small_genomes_setting(
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
    cluster_contigs: bool,
) -> std::result::Result<bool, String> {
    if cluster_contigs {
        // When clustering contigs, exactly one of --small-contigs or --large-contigs must be set
        let small_contigs =
            clap_matches.get_flag(&argument_definition.dereplication_small_contigs_argument);
        let large_contigs =
            clap_matches.get_flag(&argument_definition.dereplication_large_contigs_argument);

        match (small_contigs, large_contigs) {
            (true, false) => Ok(true),
            (false, true) => Ok(false),
            (false, false) => Err(format!(
                "When --{} is used, either --{} or --{} must be specified",
                argument_definition.dereplication_cluster_contigs_argument,
                argument_definition.dereplication_small_contigs_argument,
                argument_definition.dereplication_large_contigs_argument
            )),
            (true, true) => {
                unreachable!("clap should prevent both flags from being set due to conflicts_with")
            }
        }
    } else {
        // When not clustering contigs, use the regular --small-genomes flag
        Ok(clap_matches.get_flag(&argument_definition.dereplication_small_genomes_argument))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    lazy_static! {
        static ref PREFIXED: GalahClustererCommandDefinition =
            GalahClustererCommandDefinition::with_prefix("dereplication-");
    }

    /// Names that `with_prefix` changes, in their unprefixed (galah) form.
    fn prefixed_galah_names() -> Vec<String> {
        let g = GalahClustererCommandDefinition::galah();
        vec![
            g.dereplication_ani_argument,
            g.dereplication_prethreshold_ani_argument,
            g.dereplication_quality_formula_argument,
            g.dereplication_precluster_method_argument,
            g.dereplication_cluster_method_argument,
            g.dereplication_aligned_fraction_argument,
            g.dereplication_small_genomes_argument,
            g.dereplication_small_contigs_argument,
            g.dereplication_large_contigs_argument,
            g.dereplication_fraglen_argument,
            g.dereplication_cluster_contigs_argument,
            g.dereplication_low_memory_argument,
            g.dereplication_skip_sanitize_headers_argument,
            g.dereplication_both_min_aligned_fraction_argument,
            g.dereplication_reference_genomes_argument,
            g.dereplication_reference_genomes_list_argument,
            g.dereplication_skip_input_dereplication_argument,
            g.dereplication_output_cluster_definition_file,
            g.dereplication_output_representative_fasta_directory,
            g.dereplication_output_representative_fasta_directory_copy,
            g.dereplication_output_representative_list,
        ]
    }

    fn rendered_sections(definition: &GalahClustererCommandDefinition) -> String {
        let manual = Manual::new("test")
            .about("test")
            .author(Author::new("test"))
            .custom(add_dereplication_filtering_parameters_to_section(
                Section::new("Filtering"),
                definition,
            ))
            .custom(add_dereplication_clustering_parameters_to_section(
                Section::new("Clustering"),
                definition,
            ))
            .custom(add_dereplication_output_parameters_to_section(
                Section::new("Output"),
                definition,
            ));
        // ROFF escapes hyphens
        manual.render().replace("\\-", "-")
    }

    #[test]
    fn test_with_prefix_empty_is_galah() {
        let g = GalahClustererCommandDefinition::galah();
        assert_eq!(g.dereplication_ani_argument, "ani");
        assert_eq!(g.dereplication_prethreshold_ani_argument, "precluster-ani");
        assert_eq!(
            g.dereplication_skip_sanitize_headers_argument,
            "skip-sanitise-headers"
        );
    }

    #[test]
    fn test_with_prefix_names() {
        let d = &*PREFIXED;
        assert_eq!(d.dereplication_ani_argument, "dereplication-ani");
        assert_eq!(
            d.dereplication_prethreshold_ani_argument,
            "dereplication-precluster-ani"
        );
        // CheckM2 execution arguments are not prefixed
        assert_eq!(d.dereplication_run_checkm2_argument, "run-checkm2");
        assert_eq!(d.dereplication_checkm2_db_path_argument, "checkm2-db-path");
    }

    #[test]
    fn test_help_text_uses_prefixed_names() {
        let text = rendered_sections(&PREFIXED);
        for name in prefixed_galah_names() {
            let bare = format!("--{name}");
            assert!(
                !text.contains(&bare),
                "help text refers to unprefixed {}",
                bare
            );
        }
        // And the prefixed cross-references are present
        for expected in [
            "--dereplication-small-contigs or --dereplication-large-contigs",
            "--dereplication-reference-genomes-list",
            "--dereplication-min-aligned-fraction",
        ] {
            assert!(text.contains(expected), "missing {}", expected);
        }
    }

    #[test]
    fn test_help_text_galah_names_unchanged() {
        let text = rendered_sections(&GalahClustererCommandDefinition::galah());
        assert!(text.contains("--small-contigs or --large-contigs"));
        assert!(text.contains("--min-aligned-fraction"));
        assert!(!text.contains("--dereplication-"));
    }

    #[test]
    fn test_help_text_custom_checkm2_names() {
        let mut d = GalahClustererCommandDefinition::with_prefix("x-");
        d.dereplication_checkm2_db_path_argument = "my-db".to_string();
        let text = rendered_sections(&d);
        assert!(text.contains("--my-db or CHECKM2DB"));
    }

    #[test]
    fn test_small_genomes_error_uses_definition_names() {
        let cmd = add_dereplication_clustering_clap_args(Command::new("t"), &PREFIXED);
        let m = cmd
            .try_get_matches_from(["t", "--dereplication-cluster-contigs"])
            .unwrap();
        let err = determine_small_genomes_setting(&m, &PREFIXED, true).unwrap_err();
        assert!(err.contains("--dereplication-cluster-contigs"));
        assert!(err.contains("--dereplication-small-contigs"));
    }

    #[test]
    fn test_clap_args_prefixed() {
        let cmd = add_dereplication_output_clap_args(
            add_dereplication_clustering_clap_args(Command::new("t"), &PREFIXED),
            &PREFIXED,
        );
        cmd.clone().debug_assert();
        let m = cmd
            .try_get_matches_from([
                "t",
                "--dereplication-ani",
                "97",
                "--dereplication-reference-genomes",
                "a.fna",
                "b.fna",
                "--dereplication-output-cluster-definition",
                "out.tsv",
            ])
            .unwrap();
        assert_eq!(
            *m.get_one::<f32>(&PREFIXED.dereplication_ani_argument)
                .unwrap(),
            97.0
        );
        assert_eq!(
            m.get_many::<String>(&PREFIXED.dereplication_reference_genomes_argument)
                .unwrap()
                .count(),
            2
        );
        assert!(m.contains_id(&PREFIXED.dereplication_output_cluster_definition_file));
        // Unprefixed names must not be accepted
        assert!(
            add_dereplication_clustering_clap_args(Command::new("t"), &PREFIXED)
                .try_get_matches_from(["t", "--ani", "97"])
                .is_err()
        );
    }

    #[test]
    fn test_clap_arg_relations_prefixed() {
        let build = || add_dereplication_clustering_clap_args(Command::new("t"), &PREFIXED);
        // small-contigs requires cluster-contigs
        assert!(build()
            .try_get_matches_from(["t", "--dereplication-small-contigs"])
            .is_err());
        // small and large contigs conflict
        assert!(build()
            .try_get_matches_from([
                "t",
                "--dereplication-cluster-contigs",
                "--dereplication-small-contigs",
                "--dereplication-large-contigs"
            ])
            .is_err());
        // low-memory conflicts with reference genomes
        assert!(build()
            .try_get_matches_from([
                "t",
                "--dereplication-low-memory",
                "--dereplication-reference-genomes-list",
                "x"
            ])
            .is_err());
    }

    #[test]
    fn test_clap_args_galah_defaults() {
        let m =
            add_dereplication_clustering_clap_args(Command::new("t"), &GALAH_COMMAND_DEFINITION)
                .try_get_matches_from(["t"])
                .unwrap();
        assert_eq!(
            m.get_one::<String>("precluster-method").unwrap(),
            crate::DEFAULT_PRECLUSTER_METHOD
        );
    }
}
