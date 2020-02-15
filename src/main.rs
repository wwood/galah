extern crate galah;

extern crate clap;
use clap::*;
use std::env;

#[macro_use]
extern crate log;

extern crate bird_tool_utils;
use bird_tool_utils::clap_utils::*;

static PROGRAM_NAME: &str = "Galah";

fn main(){
    let app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false, PROGRAM_NAME, crate_version!());

    match matches.subcommand_name() {
        Some("cluster") => {
            let m = matches.subcommand_matches("cluster").unwrap();
            set_log_level(m, true, PROGRAM_NAME, crate_version!());

            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true).unwrap();

            let galah = galah::cluster_argument_parsing::generate_galah_clusterer(&genome_fasta_files, &m)
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
        },
        Some("cluster-validate") => {
            let m = matches.subcommand_matches("cluster-validate").unwrap();
            set_log_level(m, true, PROGRAM_NAME, crate_version!());

            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            let ani = galah::cluster_argument_parsing::parse_percentage(&m, "ani");

            galah::cluster_validation::validate_clusters(m.value_of("cluster-file").unwrap(), ani.unwrap().unwrap());

        },
        Some("dist") => {
            let m = matches.subcommand_matches("dist").unwrap();
            set_log_level(m, true, PROGRAM_NAME, crate_version!());

            let n_hashes = value_t!(m.value_of("num-hashes"), usize).unwrap();
            let kmer_length = value_t!(m.value_of("kmer-length"), u8).unwrap();

            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();

            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            info!("Reading CheckM tab table ..");
            let checkm = checkm::CheckMTabTable::read_file_path(
                m.value_of("checkm-tab-table").unwrap()
            );

            let genome_fasta_files = parse_list_of_genome_fasta_files(&m, true).unwrap();

            let qualities = genome_fasta_files.iter().map(|fasta|
                checkm.retrieve_via_fasta_path(fasta)
                    .expect(&format!("Failed to link genome fasta file {} to a CheckM quality", fasta))
                )
                .collect::<Vec<_>>();
            info!("Linked {} genomes to their CheckM quality", qualities.len());

            info!("Printing distances ..");
            galah::ani_correction::print_metaani_distances(
                &genome_fasta_files.iter().map(|s| s.as_str()).collect::<Vec<_>>().as_slice(),
                qualities.as_slice(),
                n_hashes, kmer_length);
            info!("Finished");
        },
        _ => panic!("Programming error")
    }
}


fn build_cli() -> App<'static, 'static> {
    let mut app = App::new("galah")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Metagenome assembled genome (MAG) dereplicator / clusterer")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("cluster-validate")
                .about("Verify clustering results")

                .arg(Arg::with_name("cluster-file")
                    .long("cluster-file")
                    .required(true)
                    .help("Output of 'cluster' subcommand")
                    .takes_value(true))
                .arg(Arg::with_name("ani")
                    .long("ani")
                    .default_value("99")
                    .help("ANI to validate against")
                    .takes_value(true))
                .arg(Arg::with_name("threads")
                    .short("t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true)));

        // .subcommand(
        //     SubCommand::with_name("dist")
        //         .about("Calculate pairwise distances between a set of genomes")

        //         .arg(Arg::with_name("checkm-tab-table")
        //             .long("checkm-tab-table")
        //             .required(true)
        //             .help("Output of CheckM lineage_wf/taxonomy_wf/qa with --tab_table specified")
        //             .takes_value(true))
        //         .arg(Arg::with_name("genome-fasta-files")
        //             .long("genome-fasta-files")
        //             .multiple(true)
        //             .required(true)
        //             .takes_value(true))
        //         .arg(Arg::with_name("num-hashes")
        //             .long("num-hashes")
        //             .takes_value(true)
        //             .default_value("1000"))
        //         .arg(Arg::with_name("kmer-length")
        //             .long("kmer-length")
        //             .takes_value(true)
        //             .default_value("21"))
        //         .arg(Arg::with_name("threads")
        //             .short("-t")
        //             .long("threads")
        //             .default_value("1")
        //             .takes_value(true)))

    app = galah::cluster_argument_parsing::add_cluster_subcommand(app);
    return app;
}
