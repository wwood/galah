use std;
use clap::*;

use super::minhash_clusterer;

pub struct GalahClusterer<'a> {
    pub genome_fasta_paths: Vec<&'a str>,
    ani: f32,
    distance_method: ClusteringDistanceMethod
}

pub enum ClusteringDistanceMethod {
    DashingFastani {
        prethreshold_ani: f32
    }
}

fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches)
    -> std::result::Result<Vec<&'a str>, String> {

    match clap_matches.is_present("checkm-tab-table") {
        false => {
            warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
            Ok(genome_fasta_files.iter().map(|s| &**s).collect())
        },
        true => {
            info!("Reading CheckM tab table ..");
            let checkm = checkm::CheckMTabTable::read_file_path(
                clap_matches.value_of("checkm-tab-table").unwrap());

            let max_contamination = match parse_percentage(clap_matches, "max-contamination") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e)
            };
            let min_completeness = match parse_percentage(clap_matches, "min-completeness") {
                Ok(fraction_opt) => fraction_opt,
                Err(e) => return Err(e)
            };

            info!("Ordering genomes by CheckM quality: completeness - 4*contamination");
            let v2 = checkm.order_fasta_paths_by_completeness_minus_4contamination(
                &genome_fasta_files.iter().map(|s| &**s).collect(),
                min_completeness,
                max_contamination)
                .unwrap();
            info!("Read in genome qualities for {} genomes. {} passed quality thresholds",
                checkm.genome_to_quality.len(),
                v2.len());
            Ok(v2)
        }
    }
}

pub fn generate_galah_clusterer<'a>(
    genome_fasta_paths: &'a Vec<String>,
    clap_matches: &clap::ArgMatches)
    -> std::result::Result<GalahClusterer<'a>,String> {


    match filter_genomes_through_checkm(
        &genome_fasta_paths, &clap_matches) {

        Err(e) => {
            std::result::Result::Err(e)
        },

        Ok(v2) => {
            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                ani: parse_percentage(&clap_matches, "ani")
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani")))
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani"))),
                distance_method: ClusteringDistanceMethod::DashingFastani  {
                    prethreshold_ani: parse_percentage(clap_matches, "prethreshold-ani")
                        .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani")))
                        .expect(&format!("Failed to parse prethreshold-ani {:?}", clap_matches.value_of("prethreshold-ani")))
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
        match self.distance_method {
            ClusteringDistanceMethod::DashingFastani{ prethreshold_ani } =>
                minhash_clusterer::minhash_clusters(
                    &self.genome_fasta_paths, prethreshold_ani*100., self.ani*100.)
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
            .required(true))
        .arg(Arg::with_name("checkm-tab-table")
            .long("checkm-tab-table")
            .help("Output of CheckM lineage_wf/taxonomy_wf/qa with --tab_table specified")
            .takes_value(true))
        .arg(Arg::with_name("min-completeness")
            .long("min-completeness")
            .help("Genomes with less than this percentage of completeness are exluded")
            .takes_value(true)
            .default_value("0"))
        .arg(Arg::with_name("max-contamination")
            .long("max-contamination")
            .help("Genomes with greater than this percentage of contamination are exluded")
            .takes_value(true))
        .arg(Arg::with_name("prethreshold-ani")
            .long("minhash-prethreshold")
            .help("Require at least this dashing-derived ANI for preclustering and to avoid FastANI on distant lineages within preclusters")
            .takes_value(true)
            .default_value("95"))
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