use std;
use std::io::Write;
#[cfg(target_family = "unix")]
use std::os::unix::fs::symlink;
#[cfg(target_family = "windows")]
use std::os::windows::fs::symlink_file as symlink;

use crate::dashing::DashingPreclusterer;
use crate::fastani::FastaniClusterer;
use crate::finch::FinchPreclusterer;
use crate::skani::SkaniClusterer;
use crate::skani::SkaniPreclusterer;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;
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
    Dashing(DashingPreclusterer),
    Finch(FinchPreclusterer),
    Skani(SkaniPreclusterer),
}

impl PreclusterDistanceFinder for Preclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        match self {
            Preclusterer::Dashing(d) => d.distances(genome_fasta_paths),
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
            Preclusterer::Dashing(d) => d.distances_contigs(genome_fasta_paths, contig_names),
            Preclusterer::Finch(f) => f.distances_contigs(genome_fasta_paths, contig_names),
            Preclusterer::Skani(s) => s.distances_contigs(genome_fasta_paths, contig_names),
        }
    }

    fn method_name(&self) -> &str {
        match self {
            Preclusterer::Dashing(d) => d.method_name(),
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
}

pub struct GalahClustererCommandDefinition {
    pub dereplication_ani_argument: String,
    pub dereplication_prethreshold_ani_argument: String,
    pub dereplication_quality_formula_argument: String,
    pub dereplication_precluster_method_argument: String,
    pub dereplication_cluster_method_argument: String,
    pub dereplication_aligned_fraction_argument: String,
    pub dereplication_small_genomes_argument: String,
    pub dereplication_small_contigs_argument: String,
    pub dereplication_large_contigs_argument: String,
    pub dereplication_fraglen_argument: String,
    pub dereplication_cluster_contigs_argument: String,
    // pub dereplication_ani_method_argument: String,
    pub dereplication_output_cluster_definition_file: String,
    pub dereplication_output_representative_fasta_directory: String,
    pub dereplication_output_representative_fasta_directory_copy: String,
    pub dereplication_output_representative_list: String,
}

lazy_static! {
    static ref GALAH_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
        GalahClustererCommandDefinition {
            dereplication_ani_argument: "ani".to_string(),
            dereplication_prethreshold_ani_argument: "precluster-ani".to_string(),
            dereplication_quality_formula_argument: "quality-formula".to_string(),
            dereplication_precluster_method_argument: "precluster-method".to_string(),
            dereplication_cluster_method_argument: "cluster-method".to_string(),
            dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
            dereplication_small_genomes_argument: "small-genomes".to_string(),
            dereplication_small_contigs_argument: "small-contigs".to_string(),
            dereplication_large_contigs_argument: "large-contigs".to_string(),
            dereplication_fraglen_argument: "fragment-length".to_string(),
            dereplication_cluster_contigs_argument: "cluster-contigs".to_string(),
            // dereplication_ani_method_argument: "ani-method".to_string(),
            dereplication_output_cluster_definition_file: "output-cluster-definition".to_string(),
            dereplication_output_representative_fasta_directory:
                "output-representative-fasta-directory".to_string(),
            dereplication_output_representative_fasta_directory_copy:
                "output-representative-fasta-directory-copy".to_string(),
            dereplication_output_representative_list: "output-representative-list".to_string(),
        }
    };
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

pub fn add_dereplication_filtering_parameters_to_section(section: Section) -> Section {
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
                clustering. {}",
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
                dereplication. '{}' for HyperLogLog, \
                '{}' for finch MinHash, '{}' for Skani. {}",
                monospace_roff("dashing"),
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
                .help("Cluster contigs within a fasta file instead of genomes. When used, either --small-contigs or --large-contigs must be specified."),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_small_contigs_argument
                ))
                .help("Use small-genomes settings in skani when clustering contigs. Recommended for contigs < 20kb. Mutually exclusive with --large-contigs."),
        )
        .flag(
            Flag::new()
                .long(&format!(
                    "--{}",
                    definition.dereplication_large_contigs_argument
                ))
                .help("Do not use small-genomes settings in skani when clustering contigs. Recommended for contigs >= 20kb. Mutually exclusive with --small-contigs."),
        )
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
                .help("Output a file of representative<TAB>member lines."),
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

    let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true).unwrap();

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
        if m.contains_id("output-representative-fasta-directory")
            || m.contains_id("output-representative-fasta-directory-copy")
        {
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
                let seq_id = String::from_utf8_lossy(seq.id()).into_owned();
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

    let galah = generate_galah_clusterer(
        &genome_fasta_files,
        &contig_names,
        cluster_contigs,
        m,
        &GALAH_COMMAND_DEFINITION,
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

pub fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
) -> std::result::Result<Vec<&'a str>, String> {
    if clap_matches.get_flag(&GALAH_COMMAND_DEFINITION.dereplication_cluster_contigs_argument) {
        return Ok(genome_fasta_files.iter().map(|s| &**s).collect());
    }

    match clap_matches.contains_id("checkm-tab-table")
        || clap_matches.contains_id("genome-info")
        || clap_matches.contains_id("checkm2-quality-report")
    {
        false => {
            warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
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

            let sorted_thresholded_genomes = match clap_matches
                .get_one::<String>(&argument_definition.dereplication_quality_formula_argument)
                .unwrap()
                .as_str()
            {
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
                "Parks2020_reduced" => {
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
) -> std::result::Result<GalahClusterer<'a>, String> {
    crate::external_command_checker::check_for_fastani();

    let skip_clusterer = {
        clap_matches
            .get_one::<String>(&argument_definition.dereplication_precluster_method_argument)
            .unwrap()
            .as_str()
            == clap_matches
                .get_one::<String>(&argument_definition.dereplication_cluster_method_argument)
                .unwrap()
                .as_str()
    };

    match filter_genomes_through_checkm(genome_fasta_paths, clap_matches, argument_definition) {
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

            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                preclusterer: match clap_matches
                    .get_one::<String>(
                        &argument_definition.dereplication_precluster_method_argument,
                    )
                    .unwrap()
                    .as_str()
                {
                    "dashing" => {
                        crate::external_command_checker::check_for_dashing();
                        Preclusterer::Dashing(DashingPreclusterer {
                            min_ani: parse_percentage(
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
                            }),
                            threads,
                        })
                    }
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
                    }),
                    _ => panic!("Programming error"),
                },
                cluster_contigs,
                contig_names,
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
        crate::clusterer::cluster(
            &self.genome_fasta_paths,
            &self.preclusterer,
            &self.clusterer,
            self.cluster_contigs,
            self.contig_names.as_deref(),
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

pub fn add_cluster_subcommand(app: clap::Command) -> clap::Command {
    let mut cluster_subcommand = add_clap_verbosity_flags(Command::new("cluster"))
        .about("Cluster FASTA files by average nucleotide identity")
        .override_help(CLUSTER_HELP.as_str())
        .arg(Arg::new("full-help").long("full-help").action(clap::ArgAction::SetTrue))
        .arg(Arg::new("full-help-roff").long("full-help-roff").action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_ani_argument)
            .long("ani")
            .help("Average nucleotide identity threshold for clustering")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ANI))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_small_genomes_argument)
            .long("small-genomes")
            .help("Use small-genomes settings in skani calculation. Recommended for sequences < 20kb.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_fraglen_argument)
            .long(&*GALAH_COMMAND_DEFINITION.dereplication_fraglen_argument)
            .help("Length of fragment used in FastANI calculation (i.e. --fragLen)")
            .value_parser(clap::value_parser!(u32))
            .default_value(crate::DEFAULT_FRAGMENT_LENGTH))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_aligned_fraction_argument)
            .long("min-aligned-fraction")
            .help("Min aligned fraction of two genomes for clustering")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ALIGNED_FRACTION))
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
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_quality_formula_argument)
            .long("quality-formula")
            .value_parser([
                "completeness-4contamination",
                "completeness-5contamination",
                "Parks2020_reduced",
                "dRep"])
            .default_value(crate::DEFAULT_QUALITY_FORMULA))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_prethreshold_ani_argument)
            .long("precluster-ani")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_PRETHRESHOLD_ANI))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_precluster_method_argument)
            .long("precluster-method")
            .help("method of calculating rough ANI. 'dashing' for HyperLogLog, 'finch' for finch MinHash, 'skani' for skani")
            .value_parser(crate::PRECLUSTER_METHODS)
            .default_value(crate::DEFAULT_PRECLUSTER_METHOD))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_cluster_method_argument)
            .long("cluster-method")
            .help("method of calculating ANI. 'fastani' for FastANI, 'skani' for skani")
            .value_parser(crate::CLUSTER_METHODS)
            .default_value(crate::DEFAULT_CLUSTER_METHOD))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_cluster_contigs_argument)
            .long("cluster-contigs")
            .help("Cluster contigs within a fasta file instead of genomes. When used, either --small-contigs or --large-contigs must be specified.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_small_contigs_argument)
            .long("small-contigs")
            .help("Use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs < 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(&*GALAH_COMMAND_DEFINITION.dereplication_cluster_contigs_argument))
        .arg(Arg::new(&*GALAH_COMMAND_DEFINITION.dereplication_large_contigs_argument)
            .long("large-contigs")
            .help("Do not use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs >= 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(&*GALAH_COMMAND_DEFINITION.dereplication_cluster_contigs_argument)
            .conflicts_with(&*GALAH_COMMAND_DEFINITION.dereplication_small_contigs_argument))
        .arg(Arg::new("threads")
            .short('t')
            .long("threads")
            .help("Number of CPU threads to use")
            .default_value("1")
            .value_parser(clap::value_parser!(u16)))
        .arg(Arg::new("output-cluster-definition")
            .short('o')
            .long("output-cluster-definition")
            .help("Output a file of representative<TAB>member lines")
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
            (false, false) => Err("When --cluster-contigs is used, either --small-contigs or --large-contigs must be specified".to_string()),
            (true, true) => unreachable!("clap should prevent both flags from being set due to conflicts_with"),
        }
    } else {
        // When not clustering contigs, use the regular --small-genomes flag
        Ok(clap_matches.get_flag(&argument_definition.dereplication_small_genomes_argument))
    }
}
