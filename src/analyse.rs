use crate::barrnap;
use crate::eukcc::{parse_eukcc_quality_report, EukccAnalyser};
use crate::isiteuk::{parse_isiteuk_tsv, IsiTeukAnalyser};
use crate::trnascan;
use crate::Domain;
use crate::DomainChoice;
use crate::QualityFinder;
use crate::RrnaFinder;
use crate::TrnaFinder;
use checkm::GenomeQuality;
use std::collections::HashMap;
use std::path::Path;
use tempfile::tempdir;

/// (r5s, r16s, r23s, r18s, r28s, r58s)
type RrnaCounts = (usize, usize, usize, usize, usize, usize);

#[derive(Debug, Clone)]
pub struct GenomeOutput {
    pub domain: String,
    pub completeness: f64,
    pub contamination: f64,
    pub r5s: usize,
    pub r16s: usize,
    pub r23s: usize,
    pub r18s: usize,
    pub r28s: usize,
    pub r58s: usize,
    pub trnas: usize,
    pub mimag_quality: String,
}

#[allow(clippy::too_many_arguments)]
pub fn analyse<Q: QualityFinder, R: RrnaFinder, T: TrnaFinder>(
    genomes: &[String],
    threads: usize,
    quality_finder: &mut Q,
    rrna_finder: &R,
    trna_finder: &T,
    checkm2_quality_report: &Option<String>,
    output_quality_report_path: &Option<String>,
    checkm_tab_table: &Option<String>,
    barrnap_gff_list: &Option<String>,
    trnascan_out_list: &Option<String>,
    domain_choice: &DomainChoice,
    isiteuk_output: &Option<String>,
    isiteuk_metapackage: Option<String>,
    eukcc_db_path: Option<String>,
    eukcc_quality_report: &Option<String>,
    bacteria_domain_cutoff: f64,
    archaea_domain_cutoff: f64,
    eukaryota_domain_cutoff: f64,
    working_dir: Option<&str>,
) -> Result<std::collections::HashMap<String, GenomeOutput>, String> {
    let quality_method = quality_finder.method_name();
    let rrna_method = rrna_finder.method_name();
    let trna_method = trna_finder.method_name();
    info!(
        "Running {}, {} and {} on provided genomes...",
        quality_method, rrna_method, trna_method
    );

    let _guard: Option<tempfile::TempDir>;
    let tmp_path_buf: std::path::PathBuf;
    match working_dir {
        Some(dir) => {
            std::fs::create_dir_all(dir).expect("Failed to create working directory");
            tmp_path_buf = std::path::PathBuf::from(dir);
            _guard = None;
        }
        None => {
            let td = tempdir().expect("Failed to create tempdir");
            tmp_path_buf = td.path().to_path_buf();
            _guard = Some(td);
        }
    }
    let tmp_path = tmp_path_buf.as_path();

    // If working_dir contains a cached CheckM2 quality report, use it automatically.
    let cached_checkm2_path = tmp_path.join("checkm2").join("quality_report.tsv");
    let effective_checkm2_report: Option<String> = checkm2_quality_report.clone().or_else(|| {
        if cached_checkm2_path.is_file() {
            info!(
                "Using cached CheckM2 quality report: {:?}",
                cached_checkm2_path
            );
            Some(cached_checkm2_path.to_string_lossy().into_owned())
        } else {
            None
        }
    });

    // ── Step 1: Determine domain assignments ─────────────────────────────────
    let domain_assignments: HashMap<String, Vec<Domain>> = match domain_choice.fixed_domains() {
        Some(domains) => genomes
            .iter()
            .map(|g| (g.clone(), domains.clone()))
            .collect(),
        None => {
            // DomainChoice::Isiteuk
            if let Some(isiteuk_path) = isiteuk_output {
                info!("Using pre-computed isiteuk output: {isiteuk_path}");
                parse_isiteuk_tsv(
                    isiteuk_path,
                    genomes,
                    bacteria_domain_cutoff,
                    archaea_domain_cutoff,
                    eukaryota_domain_cutoff,
                )
            } else {
                let metapackage = isiteuk_metapackage
                    .or_else(|| std::env::var("ISITEUK_METAPACKAGE_PATH").ok())
                    .unwrap_or_default();
                let analyser = IsiTeukAnalyser::new(
                    metapackage,
                    bacteria_domain_cutoff,
                    archaea_domain_cutoff,
                    eukaryota_domain_cutoff,
                );
                analyser.classify_genomes(genomes, threads, tmp_path)
            }
        }
    };

    // Genomes with no domain assigned fall back to all three domains.
    let all_domains = vec![Domain::Bacteria, Domain::Archaea, Domain::Eukaryota];
    let get_domains = |g: &String| -> &[Domain] {
        domain_assignments
            .get(g)
            .filter(|d| !d.is_empty())
            .map(|d| d.as_slice())
            .unwrap_or(all_domains.as_slice())
    };

    // Split genomes into prokaryote vs eukaryote-only groups for quality analysis.
    let (prok_genomes, euk_only_genomes): (Vec<&String>, Vec<&String>) = genomes
        .iter()
        .partition(|g| get_domains(g).iter().any(|d| d.is_prokaryote()));

    // ── Step 2: Quality analysis ──────────────────────────────────────────────
    let mut quality_cache: HashMap<String, (f64, f64)> = HashMap::new();

    // Prokaryotic quality (CheckM2 or pre-computed)
    if let Some(checkm2_report_path) = &effective_checkm2_report {
        info!("Using pre-generated CheckM2 quality report: {checkm2_report_path}");
        let checkm2_result = checkm::CheckM2QualityReport::read_file_path(checkm2_report_path)
            .map_err(|e| {
                format!("Failed to parse CheckM2 quality report {checkm2_report_path}: {e}")
            })?;
        for genome_path in genomes {
            let genome_stem = Path::new(genome_path)
                .file_stem()
                .unwrap()
                .to_string_lossy();
            if let Ok(q) = checkm2_result.retrieve_via_fasta_path(genome_path) {
                quality_cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            } else if let Some((_, q)) = checkm2_result
                .genome_to_quality
                .iter()
                .find(|(k, _)| **k == genome_stem)
            {
                quality_cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            }
        }
    } else if let Some(checkm_tab_path) = checkm_tab_table {
        info!("Using pre-generated CheckM tab table: {checkm_tab_path}");
        let checkm1_result = checkm::CheckM1TabTable::read_file_path(checkm_tab_path)
            .map_err(|e| format!("Failed to parse CheckM tab table {checkm_tab_path}: {e}"))?;
        for genome_path in genomes {
            let genome_stem = Path::new(genome_path)
                .file_stem()
                .unwrap()
                .to_string_lossy();
            if let Ok(q) = checkm1_result.retrieve_via_fasta_path(genome_path) {
                quality_cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            } else if let Some((_, q)) = checkm1_result
                .genome_to_quality
                .iter()
                .find(|(k, _)| **k == genome_stem)
            {
                quality_cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            }
        }
    } else if !prok_genomes.is_empty() {
        let prok_paths: Vec<String> = prok_genomes.iter().map(|g| (*g).clone()).collect();
        quality_finder.prepare_comp_cont(&prok_paths, threads, tmp_path);

        if let Some(dest) = output_quality_report_path {
            let src = tmp_path.join("checkm2").join("quality_report.tsv");
            if let Some(parent) = std::path::Path::new(dest).parent() {
                if !parent.as_os_str().is_empty() {
                    std::fs::create_dir_all(parent).map_err(|e| {
                        format!("Failed to create parent directory for quality report output: {e}")
                    })?;
                }
            }
            std::fs::copy(&src, dest).map_err(|e| {
                format!(
                    "Failed to copy CheckM2 quality report from {} to {}: {}",
                    src.display(),
                    dest,
                    e
                )
            })?;
        }

        for g in &prok_paths {
            quality_cache.insert(g.clone(), quality_finder.find_comp_cont(g));
        }
    }

    // Eukaryotic quality (EukCC or pre-computed)
    if !euk_only_genomes.is_empty() {
        let euk_paths: Vec<String> = euk_only_genomes.iter().map(|g| (*g).clone()).collect();

        if let Some(eukcc_report_path) = eukcc_quality_report {
            info!("Using pre-computed EukCC quality report: {eukcc_report_path}");
            let euk_cache = parse_eukcc_quality_report(eukcc_report_path, &euk_paths)?;
            quality_cache.extend(euk_cache);
        } else {
            let db_path = eukcc_db_path
                .or_else(|| std::env::var("EUKCC2_DB").ok())
                .unwrap_or_default();
            let mut eukcc = EukccAnalyser::new(db_path);
            eukcc.prepare_comp_cont(&euk_paths, threads, tmp_path);
            for g in &euk_paths {
                quality_cache.insert(g.clone(), eukcc.find_comp_cont(g));
            }
        }
    }

    // ── Step 3: rRNA analysis ─────────────────────────────────────────────────
    // Returns (r5s, r16s, r23s, r18s, r28s, r58s)
    let rrna_cache: HashMap<String, (usize, usize, usize, usize, usize, usize)> =
        if let Some(barrnap_list_path) = barrnap_gff_list {
            info!("Using pre-generated Barrnap GFF list: {barrnap_list_path}");
            parse_barrnap_gff_list(barrnap_list_path)?
        } else {
            genomes
                .iter()
                .map(|g| {
                    let domains = get_domains(g);
                    let result = barrnap::get_barrnap_output_for_domains(g, domains, tmp_path);
                    (g.clone(), result)
                })
                .collect()
        };

    // ── Step 4: tRNA analysis ─────────────────────────────────────────────────
    let trna_cache: HashMap<String, usize> = if let Some(trnascan_list_path) = trnascan_out_list {
        info!("Using pre-generated tRNAscan-SE output list: {trnascan_list_path}");
        parse_trnascan_out_list(trnascan_list_path)?
    } else {
        genomes
            .iter()
            .map(|g| {
                let domains = get_domains(g);
                let count = trnascan::get_trnascan_output_for_domains(g, domains, tmp_path);
                (g.clone(), count)
            })
            .collect()
    };

    // ── Step 5: Assemble GenomeOutput ─────────────────────────────────────────
    let mut genome_outputs: HashMap<String, GenomeOutput> = HashMap::new();
    for genome_path in genomes {
        let (completeness, contamination) = quality_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("Quality data not found for genome: {}", genome_path));
        let rrna = rrna_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("rRNA data not found for genome: {}", genome_path));
        let (r5s, r16s, r23s, r18s, r28s, r58s) = rrna;
        let trnas = trna_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("tRNA data not found for genome: {}", genome_path));

        let domains = get_domains(genome_path);
        let domain_str = domains
            .iter()
            .map(|d| d.display_name())
            .collect::<Vec<_>>()
            .join(",");

        // Use primary domain (first in list) to determine MIMAG criteria.
        let primary_domain = domains.first();
        let mimag_quality =
            compute_mimag_quality(primary_domain, completeness, contamination, rrna, trnas);

        genome_outputs.insert(
            genome_path.to_string(),
            GenomeOutput {
                domain: domain_str,
                completeness,
                contamination,
                r5s,
                r16s,
                r23s,
                r18s,
                r28s,
                r58s,
                trnas,
                mimag_quality: mimag_quality.to_string(),
            },
        );
    }
    Ok(genome_outputs)
}

fn compute_mimag_quality(
    primary_domain: Option<&Domain>,
    completeness: f64,
    contamination: f64,
    rrna: RrnaCounts,
    trnas: usize,
) -> &'static str {
    let (r5s, r16s, r23s, r18s, r28s, r58s) = rrna;
    match primary_domain {
        Some(Domain::Eukaryota) => {
            if completeness < 50.0 || contamination >= 10.0 {
                "Low quality"
            } else if completeness >= 90.0
                && contamination < 5.0
                && (r18s >= 1 || r28s >= 1 || r58s >= 1 || r5s >= 1)
                && trnas >= 18
            {
                "High quality"
            } else {
                "Medium quality"
            }
        }
        _ => {
            // Bacteria, Archaea, or unknown: use prokaryotic MIMAG criteria
            if completeness < 50.0 || contamination >= 10.0 {
                "Low quality"
            } else if completeness >= 90.0
                && contamination < 5.0
                && r5s >= 1
                && r16s >= 1
                && r23s >= 1
                && trnas >= 18
            {
                "High quality"
            } else {
                "Medium quality"
            }
        }
    }
}

/// Parse two-column TSV mapping genome names to Barrnap GFF files.
fn parse_barrnap_gff_list(list_path: &str) -> Result<HashMap<String, RrnaCounts>, String> {
    let mut rrna_cache = HashMap::new();
    let content = std::fs::read_to_string(list_path)
        .map_err(|e| format!("Failed to read Barrnap GFF list {list_path}: {e}"))?;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err(format!(
                "Invalid line in Barrnap GFF list (expected 2 columns): {line}"
            ));
        }
        let genome_path = fields[0].to_string();
        let gff_path = fields[1];

        let result = barrnap::parse_rrna_types_all(gff_path);
        rrna_cache.insert(genome_path, result);
    }
    Ok(rrna_cache)
}

/// Parse two-column TSV mapping genome names to tRNAscan-SE output files.
fn parse_trnascan_out_list(list_path: &str) -> Result<HashMap<String, usize>, String> {
    let mut trna_cache = HashMap::new();
    let content = std::fs::read_to_string(list_path)
        .map_err(|e| format!("Failed to read tRNAscan-SE output list {list_path}: {e}"))?;

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 2 {
            return Err(format!(
                "Invalid line in tRNAscan-SE output list (expected 2 columns): {line}"
            ));
        }
        let genome_path = fields[0].to_string();
        let out_path = fields[1];

        let trnas = crate::trnascan::count_unique_standard_trnas(out_path);
        trna_cache.insert(genome_path, trnas);
    }
    Ok(trna_cache)
}
