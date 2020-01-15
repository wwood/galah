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
}
