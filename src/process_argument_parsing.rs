use crate::analyse_argument_parsing::{
    setup_analyse_outputs, write_analyse_outputs, GalahAnalyserCommandDefinition,
    ANALYSE_COMMAND_DEFINITION,
};
use crate::cluster_argument_parsing::{
    setup_galah_outputs, write_galah_outputs, GalahClustererCommandDefinition, GalahOutput,
};
use bird_tool_utils::clap_utils::default_roff;
use bird_tool_utils::clap_utils::*;
use bird_tool_utils_man::prelude::{Author, Flag, Manual, Opt, Section};
use clap::*;

lazy_static! {
    pub static ref PROCESS_CLUSTER_COMMAND_DEFINITION: GalahClustererCommandDefinition = {
        GalahClustererCommandDefinition {
            dereplication_ani_argument: "ani".to_string(),
            dereplication_prethreshold_ani_argument: "precluster-ani".to_string(),
            dereplication_quality_formula_argument: "quality-formula".to_string(),
            dereplication_run_checkm2_argument: "run-checkm2".to_string(),
            dereplication_checkm2_db_path_argument: "checkm2-db-path".to_string(),  // Shared with analyse
            dereplication_precluster_method_argument: "precluster-method".to_string(),
            dereplication_cluster_method_argument: "cluster-method".to_string(),
            dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
            dereplication_small_genomes_argument: "small-genomes".to_string(),
            dereplication_cluster_contigs_argument: "cluster-contigs".to_string(),
            dereplication_small_contigs_argument: "small-contigs".to_string(),
            dereplication_large_contigs_argument: "large-contigs".to_string(),
            dereplication_fraglen_argument: "fragment-length".to_string(),
            dereplication_low_memory_argument: "low-memory".to_string(),
            dereplication_reference_genomes_argument: "reference-genomes".to_string(),
            dereplication_reference_genomes_list_argument: "reference-genomes-list".to_string(),
            dereplication_output_cluster_definition_file: "output-cluster-definition"
                .to_string(),
            dereplication_output_representative_fasta_directory:
                "output-representative-fasta-directory".to_string(),
            dereplication_output_representative_fasta_directory_copy:
                "output-representative-fasta-directory-copy".to_string(),
            dereplication_output_representative_list: "output-representative-list"
                .to_string(),
        }
    };
}

lazy_static! {
    static ref PROCESS_ANALYSE_COMMAND_DEFINITION: GalahAnalyserCommandDefinition = {
        GalahAnalyserCommandDefinition {
            quality_method_argument: "quality-method".to_string(),
            rrna_method_argument: "rrna-method".to_string(),
            trna_method_argument: "trna-method".to_string(),
            output_mimag_summary_argument: "output-mimag-summary".to_string(),
            output_quality_report_argument: "output-quality-report".to_string(),
            checkm2_db_path_argument: "checkm2-db-path".to_string(),  // Shared with cluster
            checkm2_quality_report_argument: "checkm2-quality-report".to_string(),  // Shared with cluster
            checkm_tab_table_argument: "checkm-tab-table".to_string(),  // Shared with cluster
            barrnap_gff_list_argument: "barrnap-gff-list".to_string(),
            trnascan_out_list_argument: "trnascan-out-list".to_string(),
        }
    };
}

lazy_static! {
    static ref PROCESS_HELP: String = format!(
        "
                     {}
              {}

{}

  {} process --genome-fasta-directory input_genomes/
    --output-mimag-summary mimag_summary.tsv
    --checkm2-db-path /path/to/checkm2_db
    --output-representative-fasta-directory output_directory/

See {} process --full-help for further options and further detail.
",
        ansi_term::Colour::Green.paint(format!(
            "{} process",
            std::env::current_exe()
                .ok()
                .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
                .and_then(|s| s.into_string().ok())
                .expect("Failed to find running program basename")
        )),
        ansi_term::Colour::Green
            .paint("Analyse (determine MIMAG status of) and Cluster (dereplicate) genomes"),
        ansi_term::Colour::Purple.paint(format!(
            "Example: Analyse the rRNA/tRNA content of a directory of .fna FASTA files, using\n\
            CheckM2 database specified by argument, and output a summary of gene counts and\n\
            MIMAG status to mimag_summary.tsv. Then dereplicate at {}% (after pre-clustering\n\
            at {}%) a directory of .fna FASTA files and create a new directory of symlinked\n\
            FASTA files of representatives:",
            crate::DEFAULT_ANI,
            crate::DEFAULT_PRETHRESHOLD_ANI
        )),
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

pub fn add_process_subcommand(app: clap::Command) -> clap::Command {
    let mut process_subcommand = add_clap_verbosity_flags(Command::new("process"))
        .about("Process genomes: run analyse and cluster together")
        .override_help(PROCESS_HELP.as_str())
        .arg(
            Arg::new("full-help")
                .long("full-help")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("full-help-roff")
                .long("full-help-roff")
                .action(clap::ArgAction::SetTrue),
        );

    // Add cluster-related arguments (ids are prefixed with cluster- to avoid collisions)
    process_subcommand = process_subcommand
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_ani_argument)
            .long("ani")
            .help("Average nucleotide identity threshold for clustering")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ANI))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_small_genomes_argument)
            .long("small-genomes")
            .help("Use small-genomes settings in skani calculation. Recommended for sequences < 20kb.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_fraglen_argument)
            .long(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_fraglen_argument)
            .help("Length of fragment used in FastANI calculation (i.e. --fragLen)")
            .value_parser(clap::value_parser!(u32))
            .default_value(crate::DEFAULT_FRAGMENT_LENGTH))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_aligned_fraction_argument)
            .long("min-aligned-fraction")
            .help("Min aligned fraction of two genomes for clustering")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_ALIGNED_FRACTION))
        .arg(Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.checkm_tab_table_argument)
            .long("checkm-tab-table")
            .help("Output of CheckM lineage_wf/taxonomy_wf/qa with --tab_table specified"))
        .arg(Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.checkm2_quality_report_argument)
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
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_quality_formula_argument)
            .long("quality-formula")
            .value_parser([
                "completeness-4contamination",
                "completeness-5contamination",
                "Parks2020_reduced",
                "dRep"])
            .default_value(crate::DEFAULT_QUALITY_FORMULA))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_run_checkm2_argument)
            .long("run-checkm2")
            .help("Run CheckM2 for genome quality scoring during clustering")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.checkm2_db_path_argument)
            .long("checkm2-db-path")
            .help("Path to CheckM2 database. If not specified, will use $CHECKM2_DB_PATH environment variable if set."))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_prethreshold_ani_argument)
            .long("precluster-ani")
            .value_parser(clap::value_parser!(f32))
            .default_value(crate::DEFAULT_PRETHRESHOLD_ANI))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_precluster_method_argument)
            .long("precluster-method")
            .help("method of calculating rough ANI. 'finch' for finch MinHash, 'skani' for skani")
            .value_parser(crate::PRECLUSTER_METHODS)
            .default_value(crate::DEFAULT_PRECLUSTER_METHOD))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_cluster_method_argument)
            .long("cluster-method")
            .help("method of calculating ANI. 'fastani' for FastANI, 'skani' for skani")
            .value_parser(crate::CLUSTER_METHODS)
            .default_value(crate::DEFAULT_CLUSTER_METHOD))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_cluster_contigs_argument)
            .long("cluster-contigs")
            .help("Cluster contigs within a fasta file instead of genomes. When used, either --small-contigs or --large-contigs must be specified.")
            .action(clap::ArgAction::SetTrue))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_small_contigs_argument)
            .long("small-contigs")
            .help("Use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs < 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_cluster_contigs_argument))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_large_contigs_argument)
            .long("large-contigs")
            .help("Do not use small-genomes settings in skani calculation when clustering contigs. Recommended for contigs >= 20kb.")
            .action(clap::ArgAction::SetTrue)
            .requires(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_cluster_contigs_argument)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_small_contigs_argument))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_low_memory_argument)
            .long(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_low_memory_argument)
            .help("Reduce memory by sketching all genomes and searching instead of triangle")
            .action(clap::ArgAction::SetTrue)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_argument)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_list_argument))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_argument)
            .long("reference-genomes")
            .help("Reference genomes to cluster against. These should be representatives already clustered. Galah will only form clusters across the two groups, never within. Uses less memory than clustering together.")
            .value_delimiter(' ')
            .num_args(1..)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_low_memory_argument)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_list_argument))
        .arg(Arg::new(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_list_argument)
            .long("reference-genomes-list")
            .help("File containing paths to reference genomes (one per line). These should be representatives already clustered. Galah will only form clusters across the two groups, never within. Uses less memory than clustering together.")
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_low_memory_argument)
            .conflicts_with(&*PROCESS_CLUSTER_COMMAND_DEFINITION.dereplication_reference_genomes_argument))
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

    // Add analyse-related arguments (ids are prefixed with analyse- to avoid collisions)
    process_subcommand = process_subcommand
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.rrna_method_argument)
                .long("rrna-method")
                .value_parser(crate::RRNA_METHODS)
                .default_value(crate::DEFAULT_RRNA_METHOD)
                .help("Method for rRNA analysis"),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.trna_method_argument)
                .long("trna-method")
                .value_parser(crate::TRNA_METHODS)
                .default_value(crate::DEFAULT_TRNA_METHOD)
                .help("Method for tRNA analysis"),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.quality_method_argument)
                .long("quality-method")
                .value_parser(crate::QUALITY_METHODS)
                .default_value(crate::DEFAULT_QUALITY_METHOD)
                .help("Method for quality analysis"),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.output_mimag_summary_argument)
                .long("output-mimag-summary")
                .value_name("SUMMARY")
                .help("Path to output MIMAG summary file")
                .required_unless_present_any([
                    "output-quality-report",
                    "full-help",
                    "full-help-roff",]),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.output_quality_report_argument)
                .long("output-quality-report")
                .value_name("REPORT")
                .help("Path to output CheckM2-format quality report (copied from CheckM2 output)")
                .required_unless_present_any([
                    "output-mimag-summary",
                    "full-help",
                    "full-help-roff",]),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.barrnap_gff_list_argument)
                .long("barrnap-gff-list")
                .value_name("FILE")
                .help("Two-column TSV mapping genome paths to Barrnap GFF paths (no headers). Prevents rRNA method being run"),
        )
        .arg(
            Arg::new(&*PROCESS_ANALYSE_COMMAND_DEFINITION.trnascan_out_list_argument)
                .long("trnascan-out-list")
                .value_name("FILE")
                .help("Two-column TSV mapping genome paths to tRNAscan-SE outputs (no headers). Prevents tRNA method being run"),
        );

    process_subcommand =
        bird_tool_utils::clap_utils::add_genome_specification_arguments(process_subcommand);

    app.subcommand(process_subcommand)
}

pub fn process_full_help(program_basename: &str, program_version: &str) -> Manual {
    // Build a full manual similar to analyse and cluster, using the same long option names
    let mut manual = Manual::new(&format!("{program_basename} process"))
        .about(format!("Process genomes: analyse and cluster (version {program_version})"))
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description("This process mode runs analyse (quality + rRNA/tRNA) then clusters genomes.\n\n\
            Required inputs are (1) a genome definition, (2) an analyse output argument, and (3) a clustering output format definition.")
        .custom_synopsis_expansion("<GENOME_INPUTS> <OUTPUT_ARGUMENTS>");

    // Genome input
    manual = manual.custom(
        bird_tool_utils::clap_utils::add_genome_specification_to_section(Section::new(
            "Genome input",
        )),
    );

    // Analyse: quality parameters (using same structure as PROCESS_ANALYSE_COMMAND_DEFINITION for consistency)
    let analyse_def_for_manual = GalahAnalyserCommandDefinition {
        quality_method_argument: "quality-method".to_string(),
        rrna_method_argument: "rrna-method".to_string(),
        trna_method_argument: "trna-method".to_string(),
        output_mimag_summary_argument: "output-mimag-summary".to_string(),
        output_quality_report_argument: "output-quality-report".to_string(),
        checkm2_db_path_argument: "checkm2-db-path".to_string(),
        checkm2_quality_report_argument: "checkm2-quality-report".to_string(),
        checkm_tab_table_argument: "checkm-tab-table".to_string(),
        barrnap_gff_list_argument: "barrnap-gff-list".to_string(),
        trnascan_out_list_argument: "trnascan-out-list".to_string(),
    };
    manual = manual.custom(
        crate::analyse_argument_parsing::add_analyse_quality_parameters_to_section(
            Section::new("Quality parameters"),
            &analyse_def_for_manual,
        ),
    );

    // Analyse: RNA parameters
    manual = manual.custom(
        crate::analyse_argument_parsing::add_analyse_rna_parameters_to_section(
            Section::new("RNA parameters"),
            &analyse_def_for_manual,
        ),
    );

    // Cluster: filtering parameters (shared long names)
    manual = manual.custom(
        crate::cluster_argument_parsing::add_dereplication_filtering_parameters_to_section(
            Section::new("Filtering parameters"),
        ),
    );

    // Cluster: clustering parameters (using same structure as cluster for consistency)
    let cluster_def_for_manual = GalahClustererCommandDefinition {
        dereplication_ani_argument: "ani".to_string(),
        dereplication_prethreshold_ani_argument: "precluster-ani".to_string(),
        dereplication_quality_formula_argument: "quality-formula".to_string(),
        dereplication_run_checkm2_argument: "run-checkm2".to_string(),
        dereplication_checkm2_db_path_argument: "checkm2-db-path".to_string(),
        dereplication_precluster_method_argument: "precluster-method".to_string(),
        dereplication_cluster_method_argument: "cluster-method".to_string(),
        dereplication_aligned_fraction_argument: "min-aligned-fraction".to_string(),
        dereplication_small_genomes_argument: "small-genomes".to_string(),
        dereplication_small_contigs_argument: "small-contigs".to_string(),
        dereplication_large_contigs_argument: "large-contigs".to_string(),
        dereplication_fraglen_argument: "fragment-length".to_string(),
        dereplication_cluster_contigs_argument: "cluster-contigs".to_string(),
        dereplication_low_memory_argument: "low-memory".to_string(),
        dereplication_reference_genomes_argument: "reference-genomes".to_string(),
        dereplication_reference_genomes_list_argument: "reference-genomes-list".to_string(),
        dereplication_output_cluster_definition_file: "output-cluster-definition".to_string(),
        dereplication_output_representative_fasta_directory:
            "output-representative-fasta-directory".to_string(),
        dereplication_output_representative_fasta_directory_copy:
            "output-representative-fasta-directory-copy".to_string(),
        dereplication_output_representative_list: "output-representative-list".to_string(),
    };
    manual = manual.custom(
        crate::cluster_argument_parsing::add_dereplication_clustering_parameters_to_section(
            Section::new("Clustering parameters"),
            &cluster_def_for_manual,
        ),
    );

    // Output
    let output_section = Section::new("Output");
    let output_section = crate::analyse_argument_parsing::add_analyse_output_parameters_to_section(
        output_section,
        &analyse_def_for_manual,
    );
    let output_section =
        crate::cluster_argument_parsing::add_dereplication_output_parameters_to_section(
            output_section,
            &cluster_def_for_manual,
        );
    manual = manual.custom(output_section);

    // General parameters
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

pub fn run_process_subcommand(
    matches: &clap::ArgMatches,
    program_basename: &str,
    program_version: &str,
) {
    use bird_tool_utils::clap_utils::*;

    let m = matches.subcommand_matches("process").unwrap();
    set_log_level(m, true, program_basename, program_version);
    bird_tool_utils::clap_utils::print_full_help_if_needed(
        m,
        process_full_help(program_basename, program_version),
    );

    let num_threads = *m.get_one::<u16>("threads").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build_global()
        .expect("Programming error: rayon initialised multiple times");

    let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true).unwrap();

    info!("Processing {} genomes ..", genome_fasta_files.len());

    // Open file handles here so errors are caught before CPU-heavy commands
    let analyse_output_definitions = setup_analyse_outputs(m, &ANALYSE_COMMAND_DEFINITION);

    let output_definitions: GalahOutput =
        setup_galah_outputs(m, &PROCESS_CLUSTER_COMMAND_DEFINITION);

    let (analysis, clusters, passed_genomes) = crate::process::process_command(
        &genome_fasta_files,
        num_threads as usize,
        m,
        &PROCESS_CLUSTER_COMMAND_DEFINITION,
        analyse_output_definitions
            .output_quality_report_path
            .clone(),
    )
    .expect("Failed to process genomes");

    write_analyse_outputs(analyse_output_definitions, &analysis, &genome_fasta_files);

    // Convert Vec<String> to Vec<&str> for write_galah_outputs
    let passed_genomes_refs: Vec<&str> = passed_genomes.iter().map(|s| s.as_str()).collect();
    write_galah_outputs(output_definitions, &clusters, &passed_genomes_refs, None);

    info!("Finished processing genomes");
}
