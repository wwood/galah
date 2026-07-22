use crate::analyse::GenomeOutput;
use crate::barrnap::BarrnapAnalyser;
use crate::checkm2::CheckM2Analyser;
use crate::cluster_argument_parsing;
use crate::trnascan::TrnascanAnalyser;
use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

/// Double-strip a genome path down to its bare name, matching the convention used by the
/// `checkm` crate's `retrieve_via_fasta_path` (strips `.gz` then the remaining extension).
fn genome_name_stem(genome_path: &str) -> String {
    let stem1 = Path::new(genome_path).file_stem().unwrap();
    if genome_path.ends_with(".gz") {
        Path::new(stem1)
            .file_stem()
            .unwrap_or(stem1)
            .to_string_lossy()
            .into_owned()
    } else {
        stem1.to_string_lossy().into_owned()
    }
}

type ProcessResult = Result<(HashMap<String, GenomeOutput>, Vec<Vec<usize>>, Vec<String>), String>;

pub fn process_command(
    genomes: &[String],
    threads: usize,
    cluster_args: &clap::ArgMatches,
    cluster_def: &cluster_argument_parsing::GalahClustererCommandDefinition,
    output_quality_report_path: Option<String>,
    process_analyse_def: &crate::process_argument_parsing::ProcessAnalyseCommandDefinition,
) -> ProcessResult {
    // Domain choice
    let domain_choice_str = cluster_args
        .get_one::<String>(&process_analyse_def.domain_choice_argument)
        .map(|s| s.as_str())
        .unwrap_or(crate::DEFAULT_DOMAIN_CHOICE);
    let domain_choice = domain_choice_str
        .parse::<crate::DomainChoice>()
        .unwrap_or(crate::DomainChoice::Isiteuk);

    // Quality analyser (CheckM2) input directly or with DB path from arg or env
    let checkm2_quality_report = cluster_args
        .get_one::<String>(&process_analyse_def.checkm2_quality_report_argument)
        .map(|s| s.to_string());
    let checkm_tab_table = cluster_args
        .get_one::<String>(&process_analyse_def.checkm_tab_table_argument)
        .map(|s| s.to_string());

    let checkm2_db_path = if checkm2_quality_report.is_none()
        && checkm_tab_table.is_none()
        && !matches!(domain_choice, crate::DomainChoice::Eukaryota)
    {
        cluster_args
            .get_one::<String>(&process_analyse_def.checkm2_db_path_argument)
            .map(|s| s.to_string())
            .or_else(|| std::env::var("CHECKM2DB").ok())
            .unwrap_or_default()
    } else {
        String::new()
    };
    let mut quality_finder = crate::analyse_argument_parsing::QualityAnalyser::CheckM2(
        CheckM2Analyser::new(checkm2_db_path),
    );

    // rRNA and tRNA analysers
    let rrna_finder = crate::analyse_argument_parsing::RrnaAnalyser::Barrnap(BarrnapAnalyser);
    let trna_finder = crate::analyse_argument_parsing::TrnaAnalyser::Trnascan(TrnascanAnalyser);

    // Input overrides for analyse (allow pre-generated files)
    let barrnap_gff_list = cluster_args
        .get_one::<String>(&process_analyse_def.barrnap_gff_list_argument)
        .map(|s| s.to_string());
    let trnascan_out_list = cluster_args
        .get_one::<String>(&process_analyse_def.trnascan_out_list_argument)
        .map(|s| s.to_string());

    // Domain / eukaryote args
    let isiteuk_output = cluster_args
        .get_one::<String>(&process_analyse_def.isiteuk_output_argument)
        .map(|s| s.to_string());
    let isiteuk_metapackage = cluster_args
        .get_one::<String>(&process_analyse_def.isiteuk_metapackage_argument)
        .map(|s| s.to_string());
    let bacteria_domain_cutoff = cluster_args
        .get_one::<f64>(&process_analyse_def.isiteuk_bacteria_cutoff_argument)
        .copied()
        .unwrap_or_else(|| crate::DEFAULT_ISITEUK_BACTERIA_CUTOFF.parse().unwrap());
    let archaea_domain_cutoff = cluster_args
        .get_one::<f64>(&process_analyse_def.isiteuk_archaea_cutoff_argument)
        .copied()
        .unwrap_or_else(|| crate::DEFAULT_ISITEUK_ARCHAEA_CUTOFF.parse().unwrap());
    let eukaryota_domain_cutoff = cluster_args
        .get_one::<f64>(&process_analyse_def.isiteuk_eukaryota_cutoff_argument)
        .copied()
        .unwrap_or_else(|| crate::DEFAULT_ISITEUK_EUKARYOTA_CUTOFF.parse().unwrap());
    let eukcc_db_path = cluster_args
        .get_one::<String>(&process_analyse_def.eukcc_db_path_argument)
        .map(|s| s.to_string());
    let eukcc_quality_report = cluster_args
        .get_one::<String>(&process_analyse_def.eukcc_quality_report_argument)
        .map(|s| s.to_string());
    let working_dir = cluster_args
        .get_one::<String>(&process_analyse_def.working_dir_argument)
        .map(|s| s.to_string());

    // Run analyse
    let analysis = crate::analyse::analyse(
        genomes,
        threads,
        &mut quality_finder,
        &rrna_finder,
        &trna_finder,
        &checkm2_quality_report,
        &output_quality_report_path,
        &checkm_tab_table,
        &barrnap_gff_list,
        &trnascan_out_list,
        &domain_choice,
        &isiteuk_output,
        isiteuk_metapackage,
        eukcc_db_path,
        &eukcc_quality_report,
        bacteria_domain_cutoff,
        archaea_domain_cutoff,
        eukaryota_domain_cutoff,
        working_dir.as_deref(),
    )?;

    // Set up clustering context similar to cluster subcommand
    let cluster_contigs =
        cluster_args.get_flag(&cluster_def.dereplication_cluster_contigs_argument);

    // Clustering contigs not yet implemented
    if cluster_contigs {
        panic!("Clustering contigs is not yet implemented in process command");
    }

    // Handle reference genomes if provided
    let reference_genomes_owned = if let Some(refs) =
        cluster_args.get_many::<String>(&cluster_def.dereplication_reference_genomes_argument)
    {
        Some(refs.cloned().collect::<Vec<String>>())
    } else if let Some(ref_file) =
        cluster_args.get_one::<String>(&cluster_def.dereplication_reference_genomes_list_argument)
    {
        let content = std::fs::read_to_string(ref_file)
            .unwrap_or_else(|_| panic!("Failed to read reference genomes list file: {}", ref_file));
        Some(
            content
                .lines()
                .filter(|line| !line.trim().is_empty())
                .map(|s| s.to_string())
                .collect::<Vec<String>>(),
        )
    } else {
        None
    };
    let reference_genomes = reference_genomes_owned
        .as_ref()
        .map(|refs| refs.iter().map(|s| s.as_str()).collect::<Vec<&str>>());

    if reference_genomes.is_some() {
        let num_reference_genomes = reference_genomes.as_ref().map_or(0, |r| r.len());
        info!(
            "Clustering against {} reference genomes",
            num_reference_genomes
        );
    }

    if reference_genomes.is_some() && cluster_contigs {
        eprintln!(
            "Error: Reference genome clustering is not currently supported with --cluster-contigs"
        );
        std::process::exit(1);
    }

    // Combine input genomes with reference genomes (when available) for quality filtering
    let (combined_genomes, ref_genomes_for_clusterer) =
        if let Some(ref_genomes) = &reference_genomes {
            let mut combined = ref_genomes
                .iter()
                .map(|s| s.to_string())
                .collect::<Vec<String>>();
            combined.extend(genomes.iter().cloned());
            (combined, Some(ref_genomes.clone()))
        } else {
            (genomes.to_vec(), None)
        };

    // Build clusterer, injecting quality for representative ranking.
    //
    // `analysis` already merges CheckM2 (Bacteria/Archaea) and EukCC (Eukaryota) completeness/
    // contamination for every genome in `genomes`, regardless of domain. Write that combined
    // table out as a genome-info-style report so representative selection can use it directly,
    // rather than relying on the CheckM2-only report (which has no entry for eukaryotes, and
    // would panic when one is looked up).
    //
    // Reference genomes are not covered by `analysis` (analyse() only runs on `genomes`), so
    // when reference genomes are in play, fall back to the pre-existing behaviour of only
    // injecting quality when the user explicitly requested a CheckM2 report be written out.
    // Kept alive until the end of this function so the path below stays valid for
    // `generate_galah_clusterer` to read; cleaned up automatically on drop.
    let mut _combined_quality_guard: Option<tempfile::NamedTempFile> = None;
    let combined_quality_report = if reference_genomes.is_none() {
        let mut combined_quality_file = tempfile::Builder::new()
            .prefix("galah-process-combined-quality")
            .suffix(".csv")
            .tempfile()
            .expect("Failed to create combined quality report tempfile");
        writeln!(combined_quality_file, "genome,completeness,contamination")
            .expect("Failed to write combined quality report header");
        for (genome_path, output) in &analysis {
            writeln!(
                combined_quality_file,
                "{},{},{}",
                genome_name_stem(genome_path),
                output.completeness,
                output.contamination
            )
            .expect("Failed to write combined quality report row");
        }
        combined_quality_file
            .flush()
            .expect("Failed to flush combined quality report");
        let path = combined_quality_file.path().to_string_lossy().into_owned();
        _combined_quality_guard = Some(combined_quality_file);
        Some(cluster_argument_parsing::InjectedQualityReport::GenomeInfo(
            path,
        ))
    } else {
        output_quality_report_path
            .clone()
            .map(cluster_argument_parsing::InjectedQualityReport::CheckM2)
    };

    let galah = cluster_argument_parsing::generate_galah_clusterer(
        &combined_genomes,
        &None,
        cluster_contigs,
        cluster_args,
        cluster_def,
        ref_genomes_for_clusterer.as_deref(),
        combined_quality_report,
    )
    .expect("Failed to parse galah clustering arguments correctly");

    let passed_genomes_owned: Vec<String> = galah
        .genome_fasta_paths
        .iter()
        .map(|s| s.to_string())
        .collect();
    info!("Clustering {} genomes ..", passed_genomes_owned.len());
    let clusters = galah.cluster();
    info!("Found {} genome clusters", clusters.len());

    Ok((analysis, clusters, passed_genomes_owned))
}
