use clap::*;
use std;

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
