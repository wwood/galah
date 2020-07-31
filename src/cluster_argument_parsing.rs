use std;
use std::io::Write;

use bird_tool_utils::clap_utils::*;
use clap::*;
use rayon::prelude::*;

use crate::genome_info_file;
use crate::genome_stats;

pub enum Preclusterer {
    Dashing {
        min_ani: f32,
        threads: usize,
    },
    Finch {
        min_ani: f32,
        num_kmers: usize,
        kmer_length: u8,
    },
}

pub struct GalahClusterer<'a> {
    pub genome_fasta_paths: Vec<&'a str>,
    pub ani: f32,
    pub min_aligned_fraction: f32,
    pub fraglen: u32,
    pub preclusterer: Preclusterer,
}

pub struct GalahClustererCommandDefinition {
    pub dereplication_ani_argument: String,
    pub dereplication_prethreshold_ani_argument: String,
    pub dereplication_quality_formula_argument: String,
    pub dereplication_precluster_method_argument: String,
    pub dereplication_aligned_fraction_argument: String,
    pub dereplication_fraglen_argument: String,
}

pub fn run_cluster_subcommand(matches: &clap::ArgMatches) {
    let m = matches.subcommand_matches("cluster").unwrap();
    set_log_level(m, true, "Galah", crate_version!());

    let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .expect("Programming error: rayon initialised multiple times");

    let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true).unwrap();

    let galah = generate_galah_clusterer(
        &genome_fasta_files,
        &m,
        &GalahClustererCommandDefinition {
            dereplication_ani_argument: "ani".to_string(),
            dereplication_prethreshold_ani_argument: "prethreshold-ani".to_string(),
            dereplication_quality_formula_argument: "quality-formula".to_string(),
            dereplication_precluster_method_argument: "precluster-method".to_string(),
            dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
            dereplication_fraglen_argument: "fraglen".to_string(),
        },
    )
    .expect("Failed to parse galah clustering arguments correctly");

    // Open file handles here so errors are caught before CPU-heavy clustering
    let output_clusters_file = match m.value_of("output-cluster-definition") {
        Some(o) => {
            Some(std::fs::File::create(o).expect("Failed to open output cluster definition file"))
        }
        None => None,
    };

    let output_representative_fasta_directory =
        setup_representative_output_directory(&m, "output-representative-fasta-directory");
    let output_representative_fasta_directory_copy =
        setup_representative_output_directory(&m, "output-representative-fasta-directory-copy");

    let output_representative_list = match m.value_of("output-representative-list") {
        Some(o) => {
            Some(std::fs::File::create(o).expect("Failed to open output representative list file"))
        }
        None => None,
    };

    let passed_genomes = &galah.genome_fasta_paths;
    info!("Clustering {} genomes ..", passed_genomes.len());
    let clusters = galah.cluster();

    info!("Found {} genome clusters", clusters.len());

    match output_clusters_file {
        Some(mut f) => {
            for cluster in &clusters {
                let rep_index = cluster[0];
                for genome_index in cluster {
                    writeln!(
                        f,
                        "{}\t{}",
                        passed_genomes[rep_index], passed_genomes[*genome_index]
                    )
                    .expect("Failed to write to output clusters file");
                }
            }
        }
        None => {}
    }

    write_cluster_reps_to_directory(
        &clusters,
        passed_genomes,
        output_representative_fasta_directory,
        &|link, current_stab, rep| {
            std::os::unix::fs::symlink(link, std::path::Path::new(&current_stab)).expect(&format!(
                "Failed to create symbolic link to representative genome {}",
                rep
            ));
        },
    );
    write_cluster_reps_to_directory(
        &clusters,
        passed_genomes,
        output_representative_fasta_directory_copy,
        &|link, current_stab, rep| {
            std::fs::copy(link, std::path::Path::new(&current_stab)).expect(&format!(
                "Failed to create symbolic link to representative genome {}",
                rep
            ));
        },
    );

    match output_representative_list {
        Some(mut f) => {
            for cluster in &clusters {
                let rep_index = cluster[0];
                writeln!(f, "{}", passed_genomes[rep_index])
                    .expect("Failed to write to output representative list file");
            }
        }
        None => {}
    }
    info!("Finished printing genome clusters");
}

fn setup_representative_output_directory<'a>(
    m: &'a clap::ArgMatches,
    argument: &str,
) -> Option<&'a std::path::Path> {
    match m.value_of(argument) {
        Some(ref d) => {
            let path = std::path::Path::new(d.clone());
            if path.exists() {
                if path.is_dir() {
                    if std::fs::read_dir(path)
                        .expect(&format!("Error opening existing output directory {}", d))
                        .collect::<Vec<_>>()
                        .len()
                        == 0
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
                std::fs::create_dir_all(path)
                    .expect(&format!("Error creating {} ({})", argument, d))
            }
            Some(path)
        }
        None => None,
    }
}

fn write_cluster_reps_to_directory(
    clusters: &Vec<Vec<usize>>,
    passed_genomes: &Vec<&str>,
    output_representative_fasta_directory: Option<&std::path::Path>,
    file_creation_fn: &dyn Fn(&std::path::Path, String, &str) -> (),
) {
    match output_representative_fasta_directory {
        Some(dir) => {
            let mut some_names_clashed = false;
            for cluster in clusters {
                let rep = passed_genomes[cluster[0]];
                let link = std::fs::canonicalize(rep).expect(&format!(
                    "Failed to convert representative path into an absolute path: {}",
                    rep
                ));

                let basename = std::path::Path::new(rep)
                    .file_name()
                    .expect(&format!("Failed to find file_name from {}", rep));
                // Check that no output files have the same name. If so, add .1.fna, .2.fna etc.
                let mut current_stab = dir.join(basename).to_str().unwrap().to_string();
                let mut counter = 0usize;
                while std::path::Path::new(&current_stab).exists() {
                    if some_names_clashed == false {
                        warn!("One or more sequence files have the same file name (e.g. ). Renaming clashes by adding .1.fna, .2.fna etc.");
                        some_names_clashed = true;
                    }
                    counter += 1;
                    current_stab =
                        format!("{}.{}.fna", dir.join(basename).to_str().unwrap(), counter);
                }
                file_creation_fn(&link, current_stab, rep);
            }
        }
        None => {}
    };
}

pub fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
) -> std::result::Result<Vec<&'a str>, String> {
    match clap_matches.is_present("checkm-tab-table") || clap_matches.is_present("genome-info") {
        false => {
            warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
            Ok(genome_fasta_files.iter().map(|s| &**s).collect())
        }
        true => {
            let checkm = if clap_matches.is_present("checkm-tab-table") {
                info!("Reading CheckM tab table ..");
                checkm::CheckMTabTable::read_file_path(
                    clap_matches.value_of("checkm-tab-table").unwrap(),
                )
            } else if clap_matches.is_present("genome-info") {
                if clap_matches
                    .value_of(&argument_definition.dereplication_quality_formula_argument)
                    .unwrap()
                    == "dRep"
                {
                    return Err(
                        "The dRep quality formula cannot be used with --genome-info".to_string()
                    );
                }
                info!(
                    "Reading genome info file {}",
                    clap_matches.value_of("genome-info").unwrap()
                );
                genome_info_file::read_genome_info_file(
                    clap_matches.value_of("genome-info").unwrap(),
                )
                .expect("Error parsing genomeInfo file")
            } else {
                panic!("Programming error");
            };

            let max_contamination = match parse_percentage(clap_matches, "max-contamination") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e),
            };
            let min_completeness = match parse_percentage(clap_matches, "min-completeness") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e),
            };
            debug!("{:?}", clap_matches);
            debug!(
                "Using max contamination {:?} and min completeness {:?}",
                max_contamination, min_completeness
            );

            let sorted_thresholded_genomes = match clap_matches
                .value_of(&argument_definition.dereplication_quality_formula_argument)
                .unwrap()
            {
                "completeness-4contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 4*contamination");
                    checkm
                        .order_fasta_paths_by_completeness_minus_4contamination(
                            &genome_fasta_files.iter().map(|s| &**s).collect(),
                            min_completeness,
                            max_contamination,
                        )
                        .unwrap()
                }
                "completeness-5contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 5*contamination");
                    checkm
                        .order_fasta_paths_by_completeness_minus_5contamination(
                            &genome_fasta_files.iter().map(|s| &**s).collect(),
                            min_completeness,
                            max_contamination,
                        )
                        .unwrap()
                }
                "Parks2020_reduced" => {
                    info!("Calculating num_contigs etc. for genome quality assessment ..");
                    let genome_appraisals = filter_and_calculate_genome_stats(
                        genome_fasta_files,
                        &checkm,
                        min_completeness,
                        max_contamination,
                    );
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
                    let genome_appraisals = filter_and_calculate_genome_stats(
                        genome_fasta_files,
                        &checkm,
                        min_completeness,
                        max_contamination,
                    );
                    let mut appraisal: Vec<_> = genome_appraisals
                        .iter()
                        // Calculate quality score
                        // completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log(N50)
                        // It's log10 specifically, see https://github.com/MrOlm/drep/blob/3ca43f20ec2b43c2c826d2f50e5473d54f4a4510/drep/d_choose.py#L202
                        .map(|(fasta_file, checkm_quality, stats)| {
                            let score = checkm_quality.completeness as f64 *100.
                                - 5.*checkm_quality.contamination as f64*100.
                                + checkm_quality.contamination as f64 * checkm_quality.strain_heterogeneity as f64
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
                checkm.genome_to_quality.len(),
                sorted_thresholded_genomes.len()
            );
            Ok(sorted_thresholded_genomes)
        }
    }
}

fn filter_and_calculate_genome_stats<'a>(
    genome_fasta_files: &'a Vec<String>,
    checkm_stats: &checkm::CheckMResult,
    min_completeness: Option<f32>,
    max_contamination: Option<f32>,
) -> Vec<(
    &'a String,
    checkm::GenomeQuality,
    genome_stats::GenomeAssemblyStats,
)> {
    genome_fasta_files
        .par_iter()
        // convert to checkm::GenomeQuality
        .map(|fasta_file| {
            (
                fasta_file,
                checkm_stats
                    .retrieve_via_fasta_path(fasta_file)
                    .expect(&format!(
                        "Failed to find CheckM statistics for {}",
                        fasta_file
                    )),
            )
        })
        // filter out poor checkm quality genomes
        .filter(|(_fasta_file, checkm_quality)| {
            let ok1 = match min_completeness {
                Some(m) => (checkm_quality.completeness >= m),
                None => true,
            };
            ok1 && match max_contamination {
                Some(m) => (checkm_quality.contamination <= m),
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
    clap_matches: &clap::ArgMatches,
    argument_definition: &GalahClustererCommandDefinition,
) -> std::result::Result<GalahClusterer<'a>, String> {
    crate::external_command_checker::check_for_fastani();

    match filter_genomes_through_checkm(&genome_fasta_paths, &clap_matches, argument_definition) {
        Err(e) => std::result::Result::Err(e),

        Ok(v2) => {
            let threads = value_t!(clap_matches.value_of("threads"), usize)
                .expect("Failed to parse --threads argument");
            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                ani: parse_percentage(
                    &clap_matches,
                    &argument_definition.dereplication_ani_argument,
                )
                .expect(&format!(
                    "Failed to parse ani {:?}",
                    clap_matches.value_of(&argument_definition.dereplication_ani_argument)
                ))
                .expect(&format!(
                    "Failed to parse ani {:?}",
                    clap_matches.value_of(&argument_definition.dereplication_ani_argument)
                )),
                min_aligned_fraction: parse_percentage(
                    &clap_matches,
                    &argument_definition.dereplication_aligned_fraction_argument,
                )
                .expect(&format!(
                    "Failed to parse min-aligned-fraction {:?}",
                    clap_matches
                        .value_of(&argument_definition.dereplication_aligned_fraction_argument)
                ))
                .expect(&format!(
                    "Failed to parse min-aligned-fraction {:?}",
                    clap_matches
                        .value_of(&argument_definition.dereplication_aligned_fraction_argument)
                )),
                fraglen: value_t!(
                    clap_matches.value_of(&argument_definition.dereplication_fraglen_argument),
                    u32
                )
                .expect(&format!(
                    "Failed to parse fragment length {:?}",
                    clap_matches.value_of(&argument_definition.dereplication_fraglen_argument)
                )),
                preclusterer: match clap_matches
                    .value_of(&argument_definition.dereplication_precluster_method_argument)
                    .unwrap()
                {
                    "dashing" => {
                        crate::external_command_checker::check_for_dashing();
                        Preclusterer::Dashing {
                            min_ani: parse_percentage(
                                clap_matches,
                                &argument_definition.dereplication_prethreshold_ani_argument,
                            )
                            .expect(&format!(
                                "Failed to parse prethreshold-ani {:?}",
                                clap_matches.value_of(
                                    &argument_definition.dereplication_prethreshold_ani_argument
                                )
                            ))
                            .expect(&format!(
                                "Failed to parse prethreshold-ani {:?}",
                                clap_matches.value_of(
                                    &argument_definition.dereplication_prethreshold_ani_argument
                                )
                            )),
                            threads: threads,
                        }
                    }
                    "finch" => Preclusterer::Finch {
                        min_ani: parse_percentage(
                            clap_matches,
                            &argument_definition.dereplication_prethreshold_ani_argument,
                        )
                        .expect(&format!(
                            "Failed to parse prethreshold-ani {:?}",
                            clap_matches.value_of(
                                &argument_definition.dereplication_prethreshold_ani_argument
                            )
                        ))
                        .expect(&format!(
                            "Failed to parse prethreshold-ani {:?}",
                            clap_matches.value_of(
                                &argument_definition.dereplication_prethreshold_ani_argument
                            )
                        )),
                        num_kmers: 1000,
                        kmer_length: 21,
                    },
                    _ => panic!("Programming error"),
                },
            })
        }
    }
}

pub fn parse_percentage(
    m: &clap::ArgMatches,
    parameter: &str,
) -> std::result::Result<Option<f32>, String> {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                let err = std::result::Result::Err(format!(
                    "Invalid percentage specified for --{}: '{}'",
                    parameter, percentage
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
        let genomes = &self.genome_fasta_paths;
        let ani = self.ani * 100.;
        let fastani_min_aligned_threshold = self.min_aligned_fraction;
        let fastani_fraglen = self.fraglen;
        match self.preclusterer {
            Preclusterer::Dashing { min_ani, threads } => crate::clusterer::cluster(
                genomes,
                &crate::dashing::DashingPreclusterer {
                    min_ani: min_ani,
                    threads: threads,
                },
                ani,
                fastani_min_aligned_threshold,
                fastani_fraglen,
            ),
            Preclusterer::Finch {
                min_ani,
                num_kmers,
                kmer_length,
            } => crate::clusterer::cluster(
                genomes,
                &crate::finch::FinchPreclusterer {
                    min_ani: min_ani,
                    num_kmers: num_kmers,
                    kmer_length: kmer_length,
                },
                ani,
                fastani_min_aligned_threshold,
                fastani_fraglen,
            ),
        }
    }
}

pub fn add_cluster_subcommand<'a>(app: clap::App<'a, 'a>) -> clap::App<'a, 'a> {
    let mut cluster_subcommand = SubCommand::with_name("cluster")
        .about("Cluster FASTA files by average nucleotide identity")
        .arg(Arg::with_name("ani")
            .long("ani")
            .help("Average nucleotide identity threshold for clustering")
            .takes_value(true)
            .default_value(crate::DEFAULT_ANI))
        .arg(Arg::with_name("fraglen")
            .long("fragment-length")
            .help("Length of fragment used in FastANI calculation (i.e. --fragLen)")
            .takes_value(true)
            .default_value(crate::DEFAULT_FRAGMENT_LENGTH))
        .arg(Arg::with_name("min-aligned-fraction")
            .long("min-aligned-fraction")
            .help("Min aligned fraction of two genomes for clustering")
            .takes_value(true)
            .default_value(crate::DEFAULT_ALIGNED_FRACTION))
        .arg(Arg::with_name("checkm-tab-table")
            .long("checkm-tab-table")
            .help("Output of CheckM lineage_wf/taxonomy_wf/qa with --tab_table specified")
            .takes_value(true))
        .arg(Arg::with_name("genome-info")
            .long("genome-info")
            .help("genomeInfo file in same format as dRep i.e. a CSV with three header columns, first line 'genome,completeness,contamination'.")
            .takes_value(true))
        .arg(Arg::with_name("min-completeness")
            .long("min-completeness")
            .help("Genomes with less than this percentage of completeness are excluded")
            .takes_value(true)
            .default_value("0"))
        .arg(Arg::with_name("max-contamination")
            .long("max-contamination")
            .help("Genomes with greater than this percentage of contamination are excluded")
            .default_value("100")
            .takes_value(true))
        .arg(Arg::with_name("quality-formula")
            .long("quality-formula")
            .help("Scoring function for genome quality. \
                'Parks2020_reduced' for 'completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000' \
                which is reduced from a quality formula described in Parks et. al. 2020\
                https://www.nature.com/articles/s41587-020-0501-8\
                'completeness-4contamination' for 'completeness-4*contamination', \
                'completeness-5contamination' for 'completeness-5*contamination', \
                'dRep' for 'completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log10(N50)'")
            .possible_values(&[
                "completeness-4contamination",
                "completeness-5contamination",
                "Parks2020_reduced",
                "dRep"])
            .default_value(crate::DEFAULT_QUALITY_FORMULA)
            .takes_value(true))
        .arg(Arg::with_name("prethreshold-ani")
            .long("prethreshold-ani")
            .help("Require at least this dashing-derived ANI for preclustering and to avoid FastANI on distant lineages within preclusters")
            .takes_value(true)
            .default_value(crate::DEFAULT_PRETHRESHOLD_ANI))
        .arg(Arg::with_name("precluster-method")
            .long("precluster-method")
            .help("method of calculating rough ANI. 'dashing' for HyperLogLog, 'finch' for finch MinHash")
            .possible_values(&["dashing","finch"])
            .default_value(crate::DEFAULT_PRECLUSTER_METHOD)
            .takes_value(true))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .help("Number of CPU threads to use")
            .default_value("1")
            .takes_value(true))
        .arg(Arg::with_name("output-cluster-definition")
            .short("o")
            .long("output-cluster-definition")
            .help("Output a file of representative<TAB>member lines")
            .required_unless_one(&[
                "output-representative-fasta-directory",
                "output-representative-fasta-directory-copy",
                "output-representative-list",])
            .takes_value(true))
        .arg(Arg::with_name("output-representative-fasta-directory")
            .long("output-representative-fasta-directory")
            .help("Symlink representative genomes into this directory")
            .required_unless_one(&[
                "output-cluster-definition",
                "output-representative-list",
                "output-representative-fasta-directory-copy",])
            .takes_value(true))
        .arg(Arg::with_name("output-representative-fasta-directory-copy")
            .long("output-representative-fasta-directory-copy")
            .help("Copy representative genomes into this directory")
            .required_unless_one(&[
                "output-cluster-definition",
                "output-representative-fasta-directory",
                "output-representative-list"])
            .takes_value(true))
        .arg(Arg::with_name("output-representative-list")
            .long("output-representative-list")
            .help("Print newline separated list of paths to representatives into this directory")
            .required_unless_one(&[
                "output-representative-fasta-directory",
                "output-representative-fasta-directory-copy",
                "output-cluster-definition"])
            .takes_value(true));

    cluster_subcommand =
        bird_tool_utils::clap_utils::add_genome_specification_arguments(cluster_subcommand);

    return app.subcommand(cluster_subcommand);
}
