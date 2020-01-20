use clap::*;

pub fn filter_genomes_through_checkm<'a>(
    genome_fasta_files: &'a Vec<String>,
    clap_matches: &clap::ArgMatches)
    -> Vec<&'a str> {

    match clap_matches.is_present("checkm-tab-table") {
        false => {
            warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
            genome_fasta_files.iter().map(|s| &**s).collect()
        },
        true => {
            info!("Reading CheckM tab table ..");
            let checkm = checkm::CheckMTabTable::read_file_path(
                clap_matches.value_of("checkm-tab-table").unwrap());

            info!("Ordering genomes by CheckM quality: completeness - 4*contamination");
            let max_contamination = match clap_matches.is_present("max-contamination") {
                true => Some(value_t!(clap_matches.value_of("max-contamination"), f32)
                    .expect("Failed to parse max-contamination to float") / 100.0),
                false => None
            };
            let v2 = checkm.order_fasta_paths_by_completeness_minus_4contamination(
                &genome_fasta_files.iter().map(|s| &**s).collect(),
                Some(value_t!(clap_matches.value_of("min-completeness"), f32)
                    .expect("Failed to parse min-completeness to float") / 100.0),
                max_contamination)
                .unwrap();
            info!("Read in genome qualities for {} genomes. {} passed quality thresholds", 
                checkm.genome_to_quality.len(), 
                v2.len());
            v2
        }
    }
}