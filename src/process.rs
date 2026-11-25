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

type ProcessResult = Result<
    (
        HashMap<String, Vec<GenomeOutput>>,
        Vec<Vec<usize>>,
        Vec<String>,
    ),
    String,
>;

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

    // --domain-choice all produces multiple rows per genome (one per domain), which has no
    // single quality value to cluster on - reject rather than guess which row to rank by.
    if matches!(domain_choice, crate::DomainChoice::All) {
        eprintln!(
            "Error: --domain-choice all is not supported by `process` (produces multiple\n\
             rows per genome, but clustering needs one quality value per genome).\n\
             Use `analyse --domain-choice all`, or `process --domain-choice completeness`."
        );
        std::process::exit(1);
    }

    // Reference genomes are not run through analyse, so have no MIMAG tier to rank on. Check
    // before running analyse so the user isn't left waiting for an error.
    let quality_formula = cluster_args
        .get_one::<String>(&cluster_def.dereplication_quality_formula_argument)
        .map(|s| s.as_str())
        .unwrap_or(crate::DEFAULT_QUALITY_FORMULA);
    if quality_formula == crate::HQ_PREFERRED_QUALITY_FORMULA
        && (cluster_args.contains_id(&cluster_def.dereplication_reference_genomes_argument)
            || cluster_args.contains_id(&cluster_def.dereplication_reference_genomes_list_argument))
    {
        eprintln!(
            "Error: {}",
            cluster_argument_parsing::hq_preferred_unavailable_message(cluster_def)
        );
        std::process::exit(1);
    }

    // Quality analyser (CheckM2) input directly or with DB path from arg or env
    let checkm2_quality_report = cluster_args
        .get_one::<String>(&process_analyse_def.checkm2_quality_report_argument)
        .map(|s| s.to_string());
    let checkm_tab_table = cluster_args
        .get_one::<String>(&process_analyse_def.checkm_tab_table_argument)
        .map(|s| s.to_string());
    // When provided, galah already has completeness/contamination for every genome, so CheckM2
    // should not be run either here (quality for the MIMAG summary) or during clustering
    // (rep-selection already honours --genome-info directly, see filter_genomes_through_checkm).
    let genome_info = cluster_args
        .get_one::<String>(&process_analyse_def.genome_info_argument)
        .map(|s| s.to_string());

    let checkm2_db_path = if checkm2_quality_report.is_none()
        && checkm_tab_table.is_none()
        && genome_info.is_none()
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
    let mut analysis = crate::analyse::analyse(
        genomes,
        threads,
        &mut quality_finder,
        &rrna_finder,
        &trna_finder,
        &checkm2_quality_report,
        &output_quality_report_path,
        &checkm_tab_table,
        &genome_info,
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
        for (genome_path, outputs) in &analysis {
            // `--domain-choice all` is rejected above, so every genome has exactly one row here.
            let output = outputs
                .first()
                .expect("Analysis produced no rows for genome");
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
        // `--domain-choice all` is rejected above, so every genome has exactly one row here.
        let mimag_quality = analysis
            .iter()
            .filter_map(|(genome_path, outputs)| {
                outputs
                    .first()
                    .map(|o| (genome_path.clone(), o.mimag_quality.clone()))
            })
            .collect();
        Some(
            cluster_argument_parsing::InjectedQualityReport::GenomeInfoWithMimag {
                path,
                mimag_quality,
            },
        )
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

    flag_lower_quality_representatives(
        &mut analysis,
        &clusters,
        &passed_genomes_owned,
        quality_formula,
        &cluster_def.dereplication_quality_formula_argument,
    );

    Ok((analysis, clusters, passed_genomes_owned))
}

fn append_note(output: &mut GenomeOutput, note: String) {
    if output.notes.is_empty() {
        output.notes = note;
    } else {
        output.notes = format!("{}; {}", output.notes, note);
    }
}

/// MIMAG tier of a genome, or None if it has no analyse row (e.g. a reference genome).
fn mimag_tier<'a>(
    analysis: &'a HashMap<String, Vec<GenomeOutput>>,
    genome: &str,
) -> Option<&'a str> {
    analysis
        .get(genome)
        .and_then(|o| o.first())
        .map(|o| o.mimag_quality.as_str())
}

/// Find clusters whose representative has a lower MIMAG tier than one or more of its members
/// (e.g. a Medium quality representative of a cluster containing a High quality genome), record
/// this in the notes of the affected genomes, and warn. Genomes without an analyse row (e.g.
/// reference genomes) are skipped, since their MIMAG tier is unknown.
fn flag_lower_quality_representatives(
    analysis: &mut HashMap<String, Vec<GenomeOutput>>,
    clusters: &[Vec<usize>],
    passed_genomes: &[String],
    quality_formula: &str,
    quality_formula_argument: &str,
) -> usize {
    use cluster_argument_parsing::mimag_tier_rank;

    let hq_preferred = quality_formula == crate::HQ_PREFERRED_QUALITY_FORMULA;
    let mut num_affected_clusters = 0;
    // (representative, its tier, higher quality member, its tier)
    let mut example: Option<(String, String, String, String)> = None;

    for cluster in clusters {
        let rep = &passed_genomes[cluster[0]];
        let rep_tier = match mimag_tier(analysis, rep) {
            Some(t) => t.to_string(),
            None => continue,
        };
        let rep_rank = mimag_tier_rank(Some(&rep_tier));
        let better_members: Vec<(String, String)> = cluster[1..]
            .iter()
            .map(|i| &passed_genomes[*i])
            .filter_map(|g| {
                mimag_tier(analysis, g)
                    .filter(|t| mimag_tier_rank(Some(t)) < rep_rank)
                    .map(|t| (g.clone(), t.to_string()))
            })
            .collect();
        if better_members.is_empty() {
            continue;
        }

        num_affected_clusters += 1;
        if example.is_none() {
            let (member, member_tier) = better_members[0].clone();
            example = Some((rep.clone(), rep_tier.clone(), member, member_tier));
        }

        let suggestion = if hq_preferred {
            String::new()
        } else {
            format!(
                "; use --{} {} to prefer higher quality representatives",
                quality_formula_argument,
                crate::HQ_PREFERRED_QUALITY_FORMULA
            )
        };
        let rep_note = format!(
            "{} cluster representative chosen over higher quality genome(s) {}{}",
            rep_tier,
            better_members
                .iter()
                .map(|(g, t)| format!("{} ({})", g, t))
                .collect::<Vec<_>>()
                .join(", "),
            suggestion
        );
        for output in analysis.get_mut(rep).unwrap() {
            append_note(output, rep_note.clone());
        }
        for (member, member_tier) in better_members {
            let member_note = format!(
                "{}, but {} genome {} chosen as cluster representative",
                member_tier, rep_tier, rep
            );
            for output in analysis.get_mut(&member).unwrap() {
                append_note(output, member_note.clone());
            }
        }
    }

    if let Some((rep, rep_tier, member, member_tier)) = example {
        let explanation = if hq_preferred {
            format!(
                "Even with {}, this can happen: a genome becomes a representative only if it is \
                below the ANI threshold to every higher-ranked representative, but each remaining \
                genome is then assigned to the representative it has the highest ANI to, which \
                may be a lower quality representative than the one it was within the threshold of.",
                crate::HQ_PREFERRED_QUALITY_FORMULA
            )
        } else {
            format!(
                "The {} quality formula ranks genomes by completeness, contamination and assembly \
                statistics rather than MIMAG tier, so e.g. a slightly more complete genome that is \
                missing a 16S rRNA gene can be chosen over a High quality genome. If this is a \
                concern, re-run with --{} {}.",
                quality_formula,
                quality_formula_argument,
                crate::HQ_PREFERRED_QUALITY_FORMULA
            )
        };
        warn!(
            "{} cluster(s) have a representative of lower MIMAG quality than another genome in the \
            cluster (e.g. {} genome {} was chosen over {} genome {}). Affected genomes are listed in \
            the notes column of the MIMAG summary. {}",
            num_affected_clusters, rep_tier, rep, member_tier, member, explanation
        );
    }

    num_affected_clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    fn output(mimag_quality: &str, notes: &str) -> Vec<GenomeOutput> {
        vec![GenomeOutput {
            domain: "Bacteria".to_string(),
            completeness: 0.0,
            contamination: 0.0,
            r5s: 0,
            r16s: 0,
            r23s: 0,
            r18s: 0,
            r28s: 0,
            r58s: 0,
            trnas: 0,
            mimag_quality: mimag_quality.to_string(),
            notes: notes.to_string(),
        }]
    }

    fn notes(analysis: &HashMap<String, Vec<GenomeOutput>>, genome: &str) -> String {
        analysis[genome][0].notes.clone()
    }

    #[test]
    fn test_flag_lower_quality_representatives() {
        let genomes: Vec<String> = ["mq", "hq", "hq2", "lq", "hq_rep", "mq2", "lq_rep", "mq3"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut analysis: HashMap<String, Vec<GenomeOutput>> = HashMap::new();
        analysis.insert("mq".into(), output("Medium quality", "existing note"));
        analysis.insert("hq".into(), output("High quality", ""));
        analysis.insert("hq2".into(), output("High quality", ""));
        analysis.insert("lq".into(), output("Low quality", ""));
        analysis.insert("hq_rep".into(), output("High quality", ""));
        analysis.insert("mq2".into(), output("Medium quality", ""));
        analysis.insert("lq_rep".into(), output("Low quality", ""));
        analysis.insert("mq3".into(), output("Medium quality", ""));
        // Cluster 1: MQ rep over two HQ members (LQ member not flagged).
        // Cluster 2: HQ rep, not flagged. Cluster 3: LQ rep over an MQ member.
        let clusters = vec![vec![0, 1, 2, 3], vec![4, 5], vec![6, 7]];

        let n = flag_lower_quality_representatives(
            &mut analysis,
            &clusters,
            &genomes,
            "Parks2020_reduced",
            "quality-formula",
        );
        assert_eq!(n, 2);
        assert_eq!(
            notes(&analysis, "mq"),
            "existing note; Medium quality cluster representative chosen over higher quality \
            genome(s) hq (High quality), hq2 (High quality); use --quality-formula \
            Parks2020_reduced_hqpreferred to prefer higher quality representatives"
        );
        assert_eq!(
            notes(&analysis, "hq"),
            "High quality, but Medium quality genome mq chosen as cluster representative"
        );
        assert_eq!(notes(&analysis, "hq2"), notes(&analysis, "hq"));
        assert_eq!(notes(&analysis, "lq"), "");
        assert_eq!(notes(&analysis, "hq_rep"), "");
        assert_eq!(notes(&analysis, "mq2"), "");
        assert_eq!(
            notes(&analysis, "lq_rep"),
            "Low quality cluster representative chosen over higher quality genome(s) mq3 \
            (Medium quality); use --quality-formula Parks2020_reduced_hqpreferred to prefer \
            higher quality representatives"
        );
        assert_eq!(
            notes(&analysis, "mq3"),
            "Medium quality, but Low quality genome lq_rep chosen as cluster representative"
        );
    }

    #[test]
    fn test_flag_hqpreferred_omits_suggestion() {
        let genomes: Vec<String> = vec!["lq".into(), "hq".into()];
        let mut analysis: HashMap<String, Vec<GenomeOutput>> = HashMap::new();
        analysis.insert("lq".into(), output("Low quality", ""));
        analysis.insert("hq".into(), output("High quality", ""));

        let n = flag_lower_quality_representatives(
            &mut analysis,
            &[vec![0, 1]],
            &genomes,
            crate::HQ_PREFERRED_QUALITY_FORMULA,
            "quality-formula",
        );
        assert_eq!(n, 1);
        assert_eq!(
            notes(&analysis, "lq"),
            "Low quality cluster representative chosen over higher quality genome(s) hq (High quality)"
        );
    }

    #[test]
    fn test_flag_skips_genomes_without_analysis() {
        // e.g. a reference genome as representative, which has no MIMAG tier.
        let genomes: Vec<String> = vec!["ref".into(), "hq".into()];
        let mut analysis: HashMap<String, Vec<GenomeOutput>> = HashMap::new();
        analysis.insert("hq".into(), output("High quality", ""));

        let n = flag_lower_quality_representatives(
            &mut analysis,
            &[vec![0, 1]],
            &genomes,
            "Parks2020_reduced",
            "quality-formula",
        );
        assert_eq!(n, 0);
        assert_eq!(notes(&analysis, "hq"), "");
    }
}
