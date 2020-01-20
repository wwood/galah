use std;
use clap::*;

use super::minhash_clusterer;

pub struct GalahClusterer<'a> {
    pub genome_fasta_paths: Vec<&'a str>,
    ani: f32,
    n_hashes: usize,
    kmer_length: u8,
    distance_method: ClusteringDistanceMethod
}

pub enum ClusteringDistanceMethod {
    MinHash,
    MinhashFastani {
        minhash_prethreshold: f32
    }
}

pub fn filter_genomes_through_checkm<'a>(
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
            let distance_method = match clap_matches.value_of("method") {
                Some("minhash") => ClusteringDistanceMethod::MinHash,
                Some("minhash+fastani") => ClusteringDistanceMethod::MinhashFastani {
                    minhash_prethreshold: parse_percentage(clap_matches, "minhash-prethreshold")
                        .expect(&format!("Failed to parse minhash prethreshold {:?}", clap_matches.value_of("minhash-prethreshold")))
                        .expect(&format!("Failed to parse minhash prethreshold {:?}", clap_matches.value_of("minhash-prethreshold")))
                },
                _ => panic!("Unexpectedly found clustering method {}",
                    clap_matches.value_of("method")
                        .expect(&format!("Failed to parse method {:?}", clap_matches.value_of("method"))))
            };

            Ok(GalahClusterer {
                genome_fasta_paths: v2,
                ani: parse_percentage(&clap_matches, "ani")
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani")))
                    .expect(&format!("Failed to parse ani {:?}", clap_matches.value_of("ani"))),
                n_hashes: value_t!(clap_matches.value_of("num-hashes"), usize).unwrap(),
                kmer_length: value_t!(clap_matches.value_of("kmer-length"), u8)
                    .expect(&format!("Failed to parse kmer-length '{:?}'", clap_matches.value_of("kmer-length"))),
                distance_method: distance_method
            })
        },
    }
}


fn parse_percentage(m: &clap::ArgMatches, parameter: &str) -> std::result::Result<Option<f32>,String> {
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
            ClusteringDistanceMethod::MinHash => minhash_clusterer::minhash_clusters(
                &self.genome_fasta_paths, self.ani*100., self.n_hashes, self.kmer_length, None),
            ClusteringDistanceMethod::MinhashFastani{ minhash_prethreshold } =>
                minhash_clusterer::minhash_clusters(
                    &self.genome_fasta_paths, minhash_prethreshold*100., self.n_hashes, self.kmer_length, Some(self.ani*100.))
        }
    }
}