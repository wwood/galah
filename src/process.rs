use crate::analyse::GenomeOutput;
use crate::barrnap::BarrnapAnalyser;
use crate::checkm2::CheckM2Analyser;
use crate::cluster_argument_parsing;
use crate::trnascan::TrnascanAnalyser;
use std::collections::HashMap;

type ProcessResult = Result<(HashMap<String, GenomeOutput>, Vec<Vec<usize>>, Vec<String>), String>;

pub fn process_command(
    genomes: &[String],
    threads: usize,
    cluster_args: &clap::ArgMatches,
    cluster_def: &cluster_argument_parsing::GalahClustererCommandDefinition,
    output_quality_report_path: Option<String>,
) -> ProcessResult {
    // Quality analyser (CheckM2) with DB path from arg or env
    let checkm2_db_path = cluster_args
        .get_one::<String>("checkm2-db-path")
        .map(|s| s.to_string())
        .or_else(|| std::env::var("CHECKM2DB").ok())
        .unwrap_or_default();
    let mut quality_finder = crate::analyse_argument_parsing::QualityAnalyser::CheckM2(
        CheckM2Analyser::new(checkm2_db_path),
    );

    // rRNA and tRNA analysers
    let rrna_finder = crate::analyse_argument_parsing::RrnaAnalyser::Barrnap(BarrnapAnalyser);
    let trna_finder = crate::analyse_argument_parsing::TrnaAnalyser::Trnascan(TrnascanAnalyser);

    // Input overrides for analyse (allow pre-generated files)
    let checkm2_quality_report = cluster_args
        .get_one::<String>("checkm2-quality-report")
        .map(|s| s.to_string());
    let checkm_tab_table = cluster_args
        .get_one::<String>("checkm-tab-table")
        .map(|s| s.to_string());
    let barrnap_gff_list = cluster_args
        .get_one::<String>("barrnap-gff-list")
        .map(|s| s.to_string());
    let trnascan_out_list = cluster_args
        .get_one::<String>("trnascan-out-list")
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

    // Build clusterer, injecting the CheckM2 quality report path from analyse (if produced)
    let galah = cluster_argument_parsing::generate_galah_clusterer(
        &combined_genomes,
        &None,
        cluster_contigs,
        cluster_args,
        cluster_def,
        ref_genomes_for_clusterer.as_deref(),
        output_quality_report_path.clone(),
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
