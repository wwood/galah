use std::io::Write;

use crate::analyse::GenomeOutput;
use crate::barrnap::BarrnapAnalyser;
use crate::checkm2::CheckM2Analyser;
use crate::trnascan::TrnascanAnalyser;
use crate::QualityFinder;
use crate::RrnaFinder;
use crate::TrnaFinder;
use bird_tool_utils::clap_utils::*;
use bird_tool_utils::clap_utils::{default_roff, monospace_roff};
use bird_tool_utils_man::prelude::{Author, Flag, Manual, Opt, Section};
use clap::*;
use std::collections::HashMap;

pub enum QualityAnalyser {
    CheckM2(crate::checkm2::CheckM2Analyser),
}

impl QualityFinder for QualityAnalyser {
    fn prepare_comp_cont(
        &mut self,
        genome_paths: &[String],
        threads: usize,
        tmp_path: &std::path::Path,
    ) {
        match self {
            QualityAnalyser::CheckM2(a) => a.prepare_comp_cont(genome_paths, threads, tmp_path),
        }
    }
    fn find_comp_cont(&self, genome_path: &str) -> (f64, f64) {
        match self {
            QualityAnalyser::CheckM2(a) => a.find_comp_cont(genome_path),
        }
    }
    fn method_name(&self) -> &str {
        match self {
            QualityAnalyser::CheckM2(a) => a.method_name(),
        }
    }
}

pub enum RrnaAnalyser {
    Barrnap(BarrnapAnalyser),
}

impl RrnaFinder for RrnaAnalyser {
    fn find_rrnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> (usize, usize, usize) {
        match self {
            RrnaAnalyser::Barrnap(r) => r.find_rrnas(genome_path, tmp_path),
        }
    }
    fn method_name(&self) -> &str {
        match self {
            RrnaAnalyser::Barrnap(r) => r.method_name(),
        }
    }
}

pub enum TrnaAnalyser {
    Trnascan(TrnascanAnalyser),
}

impl TrnaFinder for TrnaAnalyser {
    fn find_trnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> usize {
        match self {
            TrnaAnalyser::Trnascan(t) => t.find_trnas(genome_path, tmp_path),
        }
    }

    fn method_name(&self) -> &str {
        match self {
            TrnaAnalyser::Trnascan(t) => t.method_name(),
        }
    }
}

pub struct GalahAnalyser<'a> {
    pub genome_fasta_files: &'a [std::string::String],
    pub threads: usize,
    pub quality_analyser: QualityAnalyser,
    pub rrna_analyser: RrnaAnalyser,
    pub trna_analyser: TrnaAnalyser,
    pub checkm2_quality_report: Option<String>,
    pub checkm_tab_table: Option<String>,
    pub barrnap_gff_list: Option<String>,
    pub trnascan_out_list: Option<String>,
}

impl GalahAnalyser<'_> {
    pub fn analyse(
        &mut self,
        output_quality_report_path: &Option<String>,
    ) -> Result<std::collections::HashMap<String, GenomeOutput>, String> {
        crate::analyse::analyse(
            self.genome_fasta_files,
            self.threads,
            &mut self.quality_analyser,
            &self.rrna_analyser,
            &self.trna_analyser,
            &self.checkm2_quality_report,
            output_quality_report_path,
            &self.checkm_tab_table,
            &self.barrnap_gff_list,
            &self.trnascan_out_list,
        )
    }
}

pub struct GalahAnalyserCommandDefinition {
    pub quality_method_argument: String,
    pub rrna_method_argument: String,
    pub trna_method_argument: String,
    pub output_mimag_summary_argument: String,
    pub output_quality_report_argument: String,
    pub checkm2_db_path_argument: String,
    pub checkm2_quality_report_argument: String,
    pub checkm_tab_table_argument: String,
    pub barrnap_gff_list_argument: String,
    pub trnascan_out_list_argument: String,
}

lazy_static! {
    pub static ref ANALYSE_COMMAND_DEFINITION: GalahAnalyserCommandDefinition = {
        GalahAnalyserCommandDefinition {
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
        }
    };
}

lazy_static! {
    static ref PROGRAM_NAME: String = std::env::current_exe()
        .ok()
        .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
        .and_then(|s| s.into_string().ok())
        .expect("Failed to find running program basename");
}

lazy_static! {
    static ref ANALYSE_HELP: String = format!(
        "
                     {}
              {}

{}

  {} analyse --genome-fasta-directory input_genomes/
    --output-mimag-summary mimag_summary.tsv
    --checkm2-db-path /path/to/checkm2_db

{}

  CHECKM2DB=/path/to/checkm2_db
  {} analyse --genome-fasta-list genomes.txt
    --output-mimag-summary mimag_summary.tsv

{}

  {} analyse --genome-fasta-list genomes.txt
    --checkm2-quality-report quality_report.tsv
    --barrnap-gff-list barrnap_gff_list.tsv
    --trnascan-out-list trnascan_out_list.tsv
    --output-mimag-summary mimag_summary.tsv

See {} analyse --full-help for further options and further detail.
",
        ansi_term::Colour::Green.paint(format!(
            "{} analyse",
            std::env::current_exe()
                .ok()
                .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
                .and_then(|s| s.into_string().ok())
                .expect("Failed to find running program basename")
        )),
        ansi_term::Colour::Green.paint("Analyse (determine MIMAG status of) genomes"),
        ansi_term::Colour::Purple.paint(
            "Example: Analyse the rRNA/tRNA content of a directory of .fna\n\
            FASTA files, using CheckM2 database specified by argument,\n\
            and output a summary of gene counts and MIMAG status\n\
            to mimag_summary.tsv:"
                .to_string(),
        ),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        ansi_term::Colour::Purple.paint(
            "Example: Analyse a set of genomes with paths specified in genomes.txt,\n\
            using CheckM2 database specified by environment variable, and\n\
            output the MIMAG summary to mimag_summary.tsv:"
        ),
        std::env::current_exe()
            .ok()
            .and_then(|pb| pb.file_name().map(|s| s.to_os_string()))
            .and_then(|s| s.into_string().ok())
            .expect("Failed to find running program basename"),
        ansi_term::Colour::Purple.paint(
            "Example: Analyse a set of genomes with paths specified in genomes.txt,\n\
            using precomputed CheckM2, Barrnap, and tRNASCAN-SE results, and\n\
            output the MIMAG summary to mimag_summary.tsv:"
        ),
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

pub fn analyse_full_help(program_basename: &str, program_version: &str) -> Manual {
    let mut manual = Manual::new(&format!("{program_basename} analyse"))
        .about(format!("Analyse genome FASTA files to determine MIMAG status (version {program_version})"))
        .author(Author::new(crate::AUTHOR).email("benjwoodcroft near gmail.com"))
        .description("This analyse mode measures rRNA/tRNA content to determine MIMAG status. \
            Required inputs are (1) a genome definition, and (2) an output format definition.\n\n\
            The source code for this program can be found at https://github.com/wwood/galah or https://github.com/wwood/coverm")
        .custom_synopsis_expansion("<GENOME_INPUTS> <OUTPUT_ARGUMENTS>");

    // input
    manual = manual.custom(
        bird_tool_utils::clap_utils::add_genome_specification_to_section(Section::new(
            "Genome input",
        )),
    );

    // quality
    manual = manual.custom(add_analyse_quality_parameters_to_section(
        Section::new("Quality parameters"),
        &ANALYSE_COMMAND_DEFINITION,
    ));

    // rna
    manual = manual.custom(add_analyse_rna_parameters_to_section(
        Section::new("RNA parameters"),
        &ANALYSE_COMMAND_DEFINITION,
    ));

    // output
    manual = manual.custom(add_analyse_output_parameters_to_section(
        Section::new("Output"),
        &ANALYSE_COMMAND_DEFINITION,
    ));

    // general
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

pub fn add_analyse_subcommand(app: clap::Command) -> clap::Command {
    let mut analyse_subcommand = add_clap_verbosity_flags(Command::new("analyse"))
        .about("Analyse rRNAs/tRNAs of FASTA files for MIMAG status")
        .override_help(ANALYSE_HELP.as_str())
        .arg(
            Arg::new("full-help")
                .long("full-help")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("full-help-roff")
                .long("full-help-roff")
                .action(clap::ArgAction::SetTrue),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .help("Number of CPU threads to use")
                .default_value("1")
                .value_parser(clap::value_parser!(u16)),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.output_mimag_summary_argument)
                .long("output-mimag-summary")
                .value_name("SUMMARY")
                .help("Path to output MIMAG summary file")
            .required_unless_present_any([
                "output-quality-report",
                "full-help",
                "full-help-roff",]),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.output_quality_report_argument)
                .long("output-quality-report")
                .value_name("REPORT")
                .help("Path to output CheckM2-format quality report")
            .required_unless_present_any([
                "output-mimag-summary",
                "full-help",
                "full-help-roff",]),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.rrna_method_argument)
                .long("rrna-method")
                .value_parser(crate::RRNA_METHODS)
                .default_value(crate::DEFAULT_RRNA_METHOD)
                .help("Method for rRNA analysis"),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.trna_method_argument)
                .long("trna-method")
                .value_parser(crate::TRNA_METHODS)
                .default_value(crate::DEFAULT_TRNA_METHOD)
                .help("Method for tRNA analysis"),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.quality_method_argument)
                .long("quality-method")
                .value_parser(crate::QUALITY_METHODS)
                .default_value(crate::DEFAULT_QUALITY_METHOD)
                .help("Method for quality analysis"),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.checkm2_db_path_argument)
                .long("checkm2-db-path")
                .value_name("CHECKM2DB")
                .help("Path to CheckM2 database (required for checkm2 quality method) [default: from CHECKM2DB environment variable]")
                .required(false),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.checkm2_quality_report_argument)
                .long("checkm2-quality-report")
                .value_name("FILE")
                .help("Path to CheckM2 quality_report.tsv file. Prevents quality method being run")
                .required(false),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.checkm_tab_table_argument)
                .long("checkm-tab-table")
                .value_name("FILE")
                .help("Path to CheckM tab table file. Prevents quality method being run")
                .required(false),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.barrnap_gff_list_argument)
                .long("barrnap-gff-list")
                .value_name("FILE")
                .help("Two-column TSV file mapping genome paths (as given in input) to Barrnap GFF paths (no headers). Prevents rRNA method being run")
                .required(false),
        )
        .arg(
            Arg::new(&*ANALYSE_COMMAND_DEFINITION.trnascan_out_list_argument)
                .long("trnascan-out-list")
                .value_name("FILE")
                .help("Two-column TSV file mapping genome paths (as given in input) to tRNAscan-SE output paths (no headers). Prevents tRNA method being run")
                .required(false),
        );

    analyse_subcommand =
        bird_tool_utils::clap_utils::add_genome_specification_arguments(analyse_subcommand);

    app.subcommand(analyse_subcommand)
}

pub fn add_analyse_quality_parameters_to_section(
    section: Section,
    definition: &GalahAnalyserCommandDefinition,
) -> Section {
    section
        .option(
            Opt::new("NAME")
                .long(&format!("--{}", definition.quality_method_argument))
                .help(&format!(
                    "method for finding genome quality. '{}' for CheckM2. {}",
                    monospace_roff("checkm2"),
                    default_roff(crate::DEFAULT_QUALITY_METHOD)
                )),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.checkm2_db_path_argument))
                .help(
                    "Path to CheckM2 database (required for CheckM2 quality method). \
                If not given, will use CHECKM2DB environment variable if set.",
                ),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.checkm2_quality_report_argument))
                .help(
                    "Path to pre-generated CheckM2 quality_report.tsv file. \
                If given, will use this file instead of running quality method.",
                ),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.checkm_tab_table_argument))
                .help(
                    "Path to pre-generated CheckM tab table file. \
                If given, will use this file instead of running quality method.",
                ),
        )
}

pub fn add_analyse_rna_parameters_to_section(
    section: Section,
    definition: &GalahAnalyserCommandDefinition,
) -> Section {
    section
        .option(
            Opt::new("NAME")
                .long(&format!("--{}", definition.rrna_method_argument))
                .help(&format!(
                    "method for finding rRNA genes. '{}' for Barrnap. {}",
                    monospace_roff("barrnap"),
                    default_roff(crate::DEFAULT_RRNA_METHOD)
                )),
        )
        .option(
            Opt::new("NAME")
                .long(&format!("--{}", definition.trna_method_argument))
                .help(&format!(
                    "method for finding tRNA genes. '{}' for tRNAscan-SE. {}",
                    monospace_roff("trnascan"),
                    default_roff(crate::DEFAULT_TRNA_METHOD)
                )),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.barrnap_gff_list_argument))
                .help("Path to two-column TSV file mapping genome paths (as given in input) to Barrnap GFF paths (no headers). \
                If given, will use these files instead of running rRNA method."),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.trnascan_out_list_argument))
                .help("Path to two-column TSV file mapping genome paths (as given in input) to tRNAscan-SE output paths (no headers). \
                If given, will use these files instead of running tRNA method."),
        )
}

pub fn add_analyse_output_parameters_to_section(
    section: Section,
    definition: &GalahAnalyserCommandDefinition,
) -> Section {
    section
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.output_mimag_summary_argument))
                .help("Output a tsv file summarising the MIMAG status for each genome."),
        )
        .option(
            Opt::new("PATH")
                .long(&format!("--{}", definition.output_quality_report_argument))
                .help("Output a CheckM2-format quality report TSV file."),
        )
}

pub struct AnalyseOutput {
    pub output_mimag_summary: Option<std::fs::File>,
    pub output_quality_report_path: Option<String>,
}

pub fn setup_analyse_outputs(
    m: &clap::ArgMatches,
    command_definition: &GalahAnalyserCommandDefinition,
) -> AnalyseOutput {
    let output_mimag_summary = m
        .get_one::<String>(&command_definition.output_mimag_summary_argument)
        .map(|o| std::fs::File::create(o).expect("Failed to open output MIMAG summary file"));

    let output_quality_report_path = m
        .get_one::<String>(&command_definition.output_quality_report_argument)
        .map(|s| s.to_string());

    AnalyseOutput {
        output_mimag_summary,
        output_quality_report_path,
    }
}

pub fn run_analyse_subcommand(
    matches: &clap::ArgMatches,
    program_basename: &str,
    program_version: &str,
) {
    let m = matches.subcommand_matches("analyse").unwrap();
    set_log_level(m, true, program_basename, program_version);
    bird_tool_utils::clap_utils::print_full_help_if_needed(
        m,
        analyse_full_help(program_basename, program_version),
    );

    let num_threads = *m.get_one::<u16>("threads").unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads as usize)
        .build_global()
        .expect("Programming error: rayon initialised multiple times");

    let genome_fasta_files: Vec<String> = parse_list_of_genome_fasta_files(m, true).unwrap();

    let mut galah = generate_galah_analyser(&genome_fasta_files, m, &ANALYSE_COMMAND_DEFINITION)
        .expect("Failed to parse galah analyse arguments correctly");

    // Open file handles here so errors are caught before CPU-heavy commands
    let output_definitions = setup_analyse_outputs(m, &ANALYSE_COMMAND_DEFINITION);

    info!("Analysing {} genomes ..", genome_fasta_files.len());
    let analysis = galah
        .analyse(&output_definitions.output_quality_report_path)
        .expect("Failed to analyse genomes");

    write_analyse_outputs(output_definitions, &analysis, &genome_fasta_files);
    info!("Finished printing genome analysis");
}

fn generate_galah_analyser<'a>(
    genome_fasta_files: &'a [String],
    m: &ArgMatches,
    command_definition: &GalahAnalyserCommandDefinition,
) -> Result<GalahAnalyser<'a>, String> {
    let threads = *m.get_one::<u16>("threads").unwrap() as usize;

    // Quality analyser (CheckM2) input directly or with DB path from arg or env
    let checkm2_quality_report = m
        .get_one::<String>(&command_definition.checkm2_quality_report_argument)
        .map(|s| s.to_string());
    let checkm_tab_table = m
        .get_one::<String>(&command_definition.checkm_tab_table_argument)
        .map(|s| s.to_string());

    let checkm2_db_path = if checkm2_quality_report.is_none() && checkm_tab_table.is_none() {
        m.get_one::<String>(&command_definition.checkm2_db_path_argument)
            .map(|s| s.to_string())
            .or_else(|| std::env::var("CHECKM2DB").ok())
            .expect(
                "CheckM2 database path must be provided via --checkm2-db-path or CHECKM2DB env var",
            )
    } else {
        String::new()
    };

    let quality_analyser = match m
        .get_one::<String>(&command_definition.quality_method_argument)
        .map(|s| s.as_str())
    {
        Some("checkm2") => QualityAnalyser::CheckM2(CheckM2Analyser::new(checkm2_db_path)),
        _ => return Err("Invalid quality method specified".to_string()),
    };

    let rrna_analyser = match m
        .get_one::<String>(&command_definition.rrna_method_argument)
        .map(|s| s.as_str())
    {
        Some("barrnap") => RrnaAnalyser::Barrnap(BarrnapAnalyser),
        _ => return Err("Invalid rRNA method specified".to_string()),
    };

    let trna_analyser = match m
        .get_one::<String>(&command_definition.trna_method_argument)
        .map(|s| s.as_str())
    {
        Some("trnascan") => TrnaAnalyser::Trnascan(TrnascanAnalyser),
        _ => return Err("Invalid tRNA method specified".to_string()),
    };

    // Extract file input arguments
    let barrnap_gff_list = m
        .get_one::<String>(&command_definition.barrnap_gff_list_argument)
        .map(|s| s.to_string());
    let trnascan_out_list = m
        .get_one::<String>(&command_definition.trnascan_out_list_argument)
        .map(|s| s.to_string());

    Ok(GalahAnalyser {
        genome_fasta_files,
        threads,
        quality_analyser,
        rrna_analyser,
        trna_analyser,
        checkm2_quality_report,
        checkm_tab_table,
        barrnap_gff_list,
        trnascan_out_list,
    })
}

pub fn write_analyse_outputs(
    output_definitions: AnalyseOutput,
    analysis: &HashMap<String, GenomeOutput>,
    genome_fasta_files: &Vec<String>,
) {
    if let Some(mut f) = output_definitions.output_mimag_summary {
        writeln!(
            f,
            "genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality",
        )
        .unwrap();
        for genome in genome_fasta_files {
            if let Some(output_data) = analysis.get(&**genome) {
                writeln!(
                    f,
                    "{genome}\t{completeness:.2}\t{contamination:.2}\t{r5s}\t{r16s}\t{r23s}\t{trnas}\t{mimag_quality}",
                    genome = genome,
                    completeness = output_data.completeness,
                    contamination = output_data.contamination,
                    r5s = output_data.r5s,
                    r16s = output_data.r16s,
                    r23s = output_data.r23s,
                    trnas = output_data.trnas,
                    mimag_quality = output_data.mimag_quality
                )
                .unwrap();
            } else {
                writeln!(f, "{genome}\t0.0\t0.0\t0\t0\t0\t0\tMedium quality").unwrap();
            }
        }
    }
}
