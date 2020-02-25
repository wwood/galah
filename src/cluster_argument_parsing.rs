use std;
use clap::*;
use bird_tool_utils::clap_utils::*;
use rayon::prelude::*;

use crate::genome_info_file;

pub enum Preclusterer {
    Dashing {
        min_ani: f32,
        threads: usize,
    },
    Finch {
        min_ani: f32,
        num_kmers: usize,
        kmer_length: u8,
    }
}

pub struct GalahClusterer<'a> {
    pub genome_fasta_paths: Vec<&'a str>,
    pub ani: f32,
    pub preclusterer: Preclusterer,
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

    let galah = generate_galah_clusterer(&genome_fasta_files, &m)
        .expect("Failed to parse galah clustering arguments correctly");

    let passed_genomes = &galah.genome_fasta_paths;
    info!("Clustering {} genomes ..", passed_genomes.len());
    let clusters = galah.cluster();

    info!("Found {} genome clusters", clusters.len());


    for cluster in clusters {
        let rep_index = cluster[0];
        for genome_index in cluster {
            println!("{}\t{}", passed_genomes[rep_index], passed_genomes[genome_index]);
        }
    }
    info!("Finished printing genome clusters");
}

pub fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches)
    -> std::result::Result<Vec<&'a str>, String> {

    match clap_matches.is_present("checkm-tab-table") || clap_matches.is_present("genome-info") {
        false => {
            warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
            Ok(genome_fasta_files.iter().map(|s| &**s).collect())
        },
        true => {
            let checkm = if clap_matches.is_present("checkm-tab-table") {
                info!("Reading CheckM tab table ..");
                checkm::CheckMTabTable::read_file_path(
                    clap_matches.value_of("checkm-tab-table").unwrap())
            } else if clap_matches.is_present("genome-info") {
                info!("Reading genome info file {}", clap_matches.value_of("genome-info").unwrap());
                genome_info_file::read_genome_info_file(clap_matches.value_of("genome-info").unwrap())
                    .expect("Error parsing genomeInfo file")
            } else {
                panic!("Programming error");
            };

            let max_contamination = match parse_percentage(clap_matches, "max-contamination") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e)
            };
            let min_completeness = match parse_percentage(clap_matches, "min-completeness") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e)
            };

            let sorted_thresholded_genomes = match clap_matches.value_of("quality-formula").unwrap() {
                "completeness-4contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 4*contamination");
                    checkm.order_fasta_paths_by_completeness_minus_4contamination(
                        &genome_fasta_files.iter().map(|s| &**s).collect(),
                        min_completeness,
                        max_contamination)
                        .unwrap()
                },
                "completeness-5contamination" => {
                    info!("Ordering genomes by quality formula: completeness - 5*contamination");
                    checkm.order_fasta_paths_by_completeness_minus_5contamination(
                        &genome_fasta_files.iter().map(|s| &**s).collect(),
                        min_completeness,
                        max_contamination)
                        .unwrap()
                },
                "Parks2020_reduced" => {
                    info!("Calculating num_contigs etc. for genome quality assessment ..");
                    let mut appraisal: Vec<_> = genome_fasta_files
                        .par_iter()
                        // convert to checkm::GenomeQuality
                        .map(|fasta_file| {
                            (
                                fasta_file,
                                checkm.retrieve_via_fasta_path(fasta_file)
                                    .expect(&format!("Failed to find CheckM statistics for {}", fasta_file))
                            )
                        })
                        // filter out poor checkm quality genomes
                        .filter(|(_fasta_file, checkm_quality)| {
                            checkm_quality.completeness >= min_completeness.unwrap() &&
                            checkm_quality.contamination <= max_contamination.unwrap()
                        })
                        // calculate stats for good genomes
                        .map(|(fasta_file, checkm_quality)| {
                            (fasta_file, checkm_quality, crate::genome_stats::calculate_genome_stats(fasta_file))
                        })
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
                    appraisal.sort_by(|a,b| b.1.partial_cmp(&a.1).expect("Arithmetic error while calculating genome quality"));
                    appraisal.into_iter().map(|(f,_)| &**f).collect::<Vec<_>>()
                }
                _ => panic!("Programming error")
            };
            info!("Read in genome qualities for {} genomes. {} passed quality thresholds",
                checkm.genome_to_quality.len(),
                sorted_thresholded_genomes.len());
            Ok(sorted_thresholded_genomes)
        }
    }
}

pub fn generate_galah_clusterer<'a>(
    genome_fasta_paths: &'a Vec<String>,
    clap_matches: &clap::ArgMatches)
    -> std::result::Result<GalahClusterer<'a>,String> {

    crate::external_command_checker::check_for_fastani();

    match filter_genomes_through_checkm(
        &genome_fasta_paths, &clap_matches) {

        Err(e) => {
            std::result::Result::Err(e)
        },

        Ok(v2) => {
            let threads = value_t!(clap_matches.value_of("threads"), usize).expect("Failed to parse --threads argument");
            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                ani: parse_percentage(&clap_matches, "ani")
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani")))
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani"))),
                preclusterer: match clap_matches.value_of("precluster-method").unwrap() {
                    "dashing" => {
                        crate::external_command_checker::check_for_dashing();
                        Preclusterer::Dashing {
                            min_ani: parse_percentage(clap_matches, "prethreshold-ani")
                                .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani")))
                                .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani"))),
                            threads: threads
                        }
                    },
                    "finch" => Preclusterer::Finch {
                        min_ani: parse_percentage(clap_matches, "prethreshold-ani")
                            .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani")))
                            .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani"))),
                        num_kmers: 1000,
                        kmer_length: 21,
                    },
                    _ => panic!("Programming error")
                },
            })
        },
    }
}


pub fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> std::result::Result<Option<f32>,String> {
    match m.is_present(parameter) {
        true => {
            let mut percentage = value_t!(m.value_of(parameter), f32).unwrap();
            if percentage >= 1.0 && percentage <= 100.0 {
                percentage = percentage / 100.0;
            } else if percentage < 0.0 || percentage > 100.0 {
                error!("Invalid alignment percentage: '{}'", percentage);
                let err = std::result::Result::Err(format!("Invalid percentage specified for --{}: '{}'", parameter, percentage));
                return err
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
        let ani = self.ani*100.;
        match self.preclusterer {
            Preclusterer::Dashing {
                min_ani,
                threads,
            } => {
                crate::clusterer::cluster(
                    genomes,
                    &crate::dashing::DashingPreclusterer {
                        min_ani: min_ani,
                        threads: threads
                    },
                    ani
                )
            }
            Preclusterer::Finch {
                min_ani,
                num_kmers,
                kmer_length,
            } => {
                crate::clusterer::cluster(
                    genomes,
                    &crate::finch::FinchPreclusterer {
                        min_ani: min_ani,
                        num_kmers: num_kmers,
                        kmer_length: kmer_length,
                    },
                    ani
                )
            }
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
            .default_value("99"))
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
            .help("Genomes with less than this percentage of completeness are exluded")
            .takes_value(true)
            .default_value("0"))
        .arg(Arg::with_name("max-contamination")
            .long("max-contamination")
            .help("Genomes with greater than this percentage of contamination are exluded")
            .default_value("100")
            .takes_value(true))
        .arg(Arg::with_name("quality-formula")
            .long("quality-formula")
            .help("Scoring function for genome quality. \
                'completeness-4contamination' for 'completeness-4*contamination', \
                'completeness-5contamination' for 'completeness-5*contamination', \
                'Parks2020_reduced' for 'completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000' \
                which is reduced from a quality formula described in Parks et. al. 2020\
                https://www.biorxiv.org/content/10.1101/771964v2.abstract")
            .possible_values(&[
                "completeness-4contamination",
                "completeness-5contamination",
                "Parks2020_reduced"])
            .default_value("Parks2020_reduced")
            .takes_value(true))
        .arg(Arg::with_name("prethreshold-ani")
            .long("prethreshold-ani")
            .help("Require at least this dashing-derived ANI for preclustering and to avoid FastANI on distant lineages within preclusters")
            .takes_value(true)
            .default_value("95"))
        .arg(Arg::with_name("precluster-method")
            .long("precluster-method")
            .help("method of calculating rough ANI. 'dashing' for HyperLogLog, 'finch' for finch MinHash")
            .possible_values(&["dashing","finch"])
            .takes_value(true))
        .arg(Arg::with_name("threads")
            .short("t")
            .long("threads")
            .help("Number of CPU threads to use")
            .default_value("1")
            .takes_value(true));

    cluster_subcommand = bird_tool_utils::clap_utils::add_genome_specification_arguments(
        cluster_subcommand);

    return app.subcommand(cluster_subcommand);
}