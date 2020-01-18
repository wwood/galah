extern crate galah;

extern crate clap;
use clap::*;
use std::env;
use std::process;

#[macro_use]
extern crate log;
extern crate env_logger;
use log::LevelFilter;
use env_logger::Builder;
extern crate rayon;

fn main(){
    let app = build_cli();
    let matches = app.clone().get_matches();
    set_log_level(&matches, false);

    match matches.subcommand_name() {
        Some("cluster") => {
            let m = matches.subcommand_matches("cluster").unwrap();
            set_log_level(m, true);

            let num_threads = value_t!(m.value_of("threads"), usize).unwrap();
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .expect("Programming error: rayon initialised multiple times");

            let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m);

            let v2: Vec<&str> = match m.is_present("checkm-tab-table") {
                false => {
                    warn!("Since CheckM input is missing, genomes are not being ordered by quality. Instead the order of their input is being used");
                    genome_fasta_files.iter().map(|s| &**s).collect()
                },
                true => {
                    info!("Reading CheckM tab table ..");
                    let checkm = checkm::CheckMTabTable::read_file_path(m.value_of("checkm-tab-table").unwrap());

                    info!("Ordering genomes by CheckM quality: completeness - 4*contamination");
                    let v2 = checkm.order_fasta_paths_by_completeness_minus_4contamination(
                        &genome_fasta_files.iter().map(|s| &**s).collect(),
                        Some(value_t!(m.value_of("min-completeness"), f32).expect("Failed to parse min-completeness to float") / 100.0),
                        Some(value_t!(m.value_of("max-contamination"), f32).expect("Failed to parse max-contamination to float") / 100.0))
                        .unwrap();
                    info!("Read in genome qualities for {} genomes. {} passed quality thresholds", 
                        checkm.genome_to_quality.len(), 
                        v2.len());
                    v2
                }
            };
            info!("Clustering {} genomes ..", v2.len());

            let ani = value_t!(m.value_of("ani"), f32).unwrap();
            let n_hashes = value_t!(m.value_of("num-hashes"), usize).unwrap();
            let kmer_length = value_t!(m.value_of("kmer-length"), u8).unwrap();
            let clusters = match m.value_of("method") {
                Some("minhash") => galah::minhash_clusterer::minhash_clusters(
                    &v2, ani, n_hashes, kmer_length, None),
                Some("minhash+fastani") => galah::minhash_clusterer::minhash_clusters(
                    &v2, 
                    value_t!(m.value_of("minhash-prethreshold"), f32).expect("Failed to parse --minhash-prethreshold parameter"),
                    n_hashes, kmer_length, Some(ani)),
                _ => unreachable!()
            };
            info!("Found {} genome clusters", clusters.len());

            for cluster in clusters {
                let rep_index = cluster[0];
                for genome_index in cluster {
                    println!("{}\t{}", v2[rep_index], v2[genome_index]);
                }
            }
            info!("Finished printing genome clusters");
        },
        Some("dist") => {
            let m = matches.subcommand_matches("dist").unwrap();
            set_log_level(m, true);

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

            let genome_fasta_files = parse_list_of_genome_fasta_files(&m);
            
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

fn set_log_level(matches: &clap::ArgMatches, is_last: bool) {
    let mut log_level = LevelFilter::Info;
    let mut specified = false;
    if matches.is_present("verbose") {
        specified = true;
        log_level = LevelFilter::Debug;
    }
    if matches.is_present("quiet") {
        specified = true;
        log_level = LevelFilter::Error;
    }
    if specified || is_last {
        let mut builder = Builder::new();
        builder.filter_level(log_level);
        if env::var("RUST_LOG").is_ok() {
            builder.parse_filters(&env::var("RUST_LOG").unwrap());
        }
        if builder.try_init().is_err() {
            panic!("Failed to set log level - has it been specified multiple times?")
        }
    }
    if is_last {
        info!("Cockatoo version {}", crate_version!());
    }
}

fn parse_list_of_genome_fasta_files(m: &clap::ArgMatches) -> Vec<String> {
    match m.is_present("genome-fasta-files") {
        true => {
            m.values_of("genome-fasta-files").unwrap().map(|s| s.to_string()).collect()
        },
        false => {
            if m.is_present("genome-fasta-directory") {
                let dir = m.value_of("genome-fasta-directory").unwrap();
                let paths = std::fs::read_dir(dir).unwrap();
                let mut genome_fasta_files: Vec<String> = vec!();
                let extension = m.value_of("genome-fasta-extension").unwrap();
                for path in paths {
                    let file = path.unwrap().path();
                    match file.extension() {
                        Some(ext) => {
                            if ext == extension {
                                let s = String::from(file.to_string_lossy());
                                genome_fasta_files.push(s);
                            } else {
                                info!(
                                    "Not using directory entry '{}' as a genome FASTA file, as \
                                     it does not end with the extension '{}'",
                                    file.to_str().expect("UTF8 error in filename"),
                                    extension);
                            }
                        },
                        None => {
                            info!("Not using directory entry '{}' as a genome FASTA file",
                                  file.to_str().expect("UTF8 error in filename"));
                        }
                    }
                }
                if genome_fasta_files.len() == 0 {
                    error!("Found 0 genomes from the genome-fasta-directory, cannot continue.");
                    process::exit(1);
                }
                genome_fasta_files // return
            } else {
                error!("Either a separator (-s) or path(s) to genome FASTA files \
                        (with -d or -f) must be given");
                process::exit(1);
            }
        }
    }
}


fn build_cli() -> App<'static, 'static> {
    return App::new("meta_ani")
        .version(crate_version!())
        .author("Ben J. Woodcroft <benjwoodcroft near gmail.com>")
        .about("Quality-corrected ANI calculator for metagenome assembled genomes (MAGs)")
        .args_from_usage("-v, --verbose       'Print extra debug logging information'
             -q, --quiet         'Unless there is an error, do not print logging information'")
        .global_setting(AppSettings::ArgRequiredElseHelp)
        .subcommand(
            SubCommand::with_name("dist")
                .about("Calculate pairwise distances between a set of genomes")

                .arg(Arg::with_name("checkm-tab-table")
                    .long("checkm-tab-table")
                    .required(true)
                    .takes_value(true))
                .arg(Arg::with_name("genome-fasta-files")
                    .long("genome-fasta-files")
                    .multiple(true)
                    .required(true)
                    .takes_value(true))
                .arg(Arg::with_name("num-hashes")
                    .long("num-hashes")
                    .takes_value(true)
                    .default_value("1000"))
                .arg(Arg::with_name("kmer-length")
                    .long("kmer-length")
                    .takes_value(true)
                    .default_value("21"))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true)))

        .subcommand(
            SubCommand::with_name("cluster")
                .about("Cluster FASTA files by average nucleotide identity")
                .arg(Arg::with_name("ani")
                    .long("ani")
                    .takes_value(true)
                    .required(true))
                .arg(Arg::with_name("checkm-tab-table")
                    .long("checkm-tab-table")
                    .takes_value(true))
                .arg(Arg::with_name("min-completeness")
                    .long("min-completeness")
                    .takes_value(true)
                    .default_value("0"))
                .arg(Arg::with_name("max-contamination")
                    .long("max-contamination")
                    .takes_value(true)
                    .default_value("0"))
                .arg(Arg::with_name("num-hashes")
                    .long("num-hashes")
                    .takes_value(true)
                    .default_value("1000"))
                .arg(Arg::with_name("kmer-length")
                    .long("kmer-length")
                    .takes_value(true)
                    .default_value("21"))
                .arg(Arg::with_name("minhash-prethreshold")
                    .long("minhash-prethreshold")
                    .takes_value(true)
                    .default_value("90"))
                .arg(Arg::with_name("genome-fasta-files")
                        .short("f")
                        .long("genome-fasta-files")
                        .multiple(true)
                        .conflicts_with("genome-fasta-directory")
                        .conflicts_with("single-genome")
                        .required_unless_one(
                            &["genome-fasta-directory"])
                        .takes_value(true))
                .arg(Arg::with_name("genome-fasta-directory")
                        .long("genome-fasta-directory")
                        .conflicts_with("separator")
                        .conflicts_with("genome-fasta-files")
                        .conflicts_with("single-genome")
                        .required_unless_one(
                            &["genome-fasta-files"])
                        .takes_value(true))
                .arg(Arg::with_name("genome-fasta-extension")
                        .short("x")
                        .long("genome-fasta-extension")
                        // Unsure why, but uncommenting causes test failure (in
                        // genome mode, not sure about here) - clap bug?
                        //.requires("genome-fasta-directory")
                        .default_value("fna")
                        .takes_value(true))

                .arg(Arg::with_name("method")
                    .long("method")
                    .possible_values(&["minhash+fastani","minhash"])
                    .default_value("minhash+fastani")
                    .takes_value(true))
                .arg(Arg::with_name("threads")
                    .short("-t")
                    .long("threads")
                    .default_value("1")
                    .takes_value(true)))
            
}
