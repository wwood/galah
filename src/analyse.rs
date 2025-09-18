use crate::QualityFinder;
use crate::RrnaFinder;
use crate::TrnaFinder;
use checkm::GenomeQuality;
use std::collections::HashMap;
use std::path::Path;
use tempfile::tempdir;

#[derive(Debug, Clone)]
pub struct GenomeOutput {
    pub completeness: f64,
    pub contamination: f64,
    pub r5s: usize,
    pub r16s: usize,
    pub r23s: usize,
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
    checkm_tab_table: &Option<String>,
    barrnap_gff_list: &Option<String>,
    trnascan_out_list: &Option<String>,
) -> Result<std::collections::HashMap<String, GenomeOutput>, String> {
    let quality_method = quality_finder.method_name();
    let rrna_method = rrna_finder.method_name();
    let trna_method = trna_finder.method_name();
    info!(
        "Running {}, {} and {} on provided genomes...",
        quality_method, rrna_method, trna_method
    );

    let tmpdir = tempdir().expect("Failed to create tempdir");
    let tmp_path = tmpdir.path();

    // Handle quality analysis - use pre-generated files if available
    let quality_cache = if let Some(checkm2_report_path) = checkm2_quality_report {
        info!("Using pre-generated CheckM2 quality report: {checkm2_report_path}");
        let checkm2_result = checkm::CheckM2QualityReport::read_file_path(checkm2_report_path)
            .map_err(|e| {
                format!("Failed to parse CheckM2 quality report {checkm2_report_path}: {e}")
            })?;
        let mut cache = HashMap::new();
        for genome_path in genomes {
            let genome_stem = Path::new(genome_path)
                .file_stem()
                .unwrap()
                .to_string_lossy();
            if let Ok(q) = checkm2_result.retrieve_via_fasta_path(genome_path) {
                cache.insert(
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
                cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            } else {
                return Err(format!(
                    "No CheckM2 quality found for genome {genome_path} (stem {genome_stem})"
                ));
            }
        }
        cache
    } else if let Some(checkm_tab_path) = checkm_tab_table {
        info!("Using pre-generated CheckM tab table: {checkm_tab_path}");
        let checkm1_result = checkm::CheckM1TabTable::read_file_path(checkm_tab_path)
            .map_err(|e| format!("Failed to parse CheckM tab table {checkm_tab_path}: {e}"))?;
        let mut cache = HashMap::new();
        for genome_path in genomes {
            let genome_stem = Path::new(genome_path)
                .file_stem()
                .unwrap()
                .to_string_lossy();
            if let Ok(q) = checkm1_result.retrieve_via_fasta_path(genome_path) {
                cache.insert(
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
                cache.insert(
                    genome_path.clone(),
                    (
                        q.completeness() as f64 * 100.0,
                        q.contamination() as f64 * 100.0,
                    ),
                );
            } else {
                return Err(format!(
                    "No CheckM1 quality found for genome {genome_path} (stem {genome_stem})"
                ));
            }
        }
        cache
    } else {
        quality_finder.prepare_comp_cont(genomes, threads, tmp_path);
        genomes
            .iter()
            .map(|g| (g.clone(), quality_finder.find_comp_cont(g)))
            .collect()
    };

    // Handle rRNA analysis - use pre-generated files if available
    let rrna_cache = if let Some(barrnap_list_path) = barrnap_gff_list {
        info!("Using pre-generated Barrnap GFF list: {barrnap_list_path}");
        parse_barrnap_gff_list(barrnap_list_path)?
    } else {
        genomes
            .iter()
            .map(|g| (g.clone(), rrna_finder.find_rrnas(g, tmp_path)))
            .collect()
    };

    // Handle tRNA analysis - use pre-generated files if available
    let trna_cache = if let Some(trnascan_list_path) = trnascan_out_list {
        info!("Using pre-generated tRNAscan-SE output list: {trnascan_list_path}");
        parse_trnascan_out_list(trnascan_list_path)?
    } else {
        genomes
            .iter()
            .map(|g| (g.clone(), trna_finder.find_trnas(g, tmp_path)))
            .collect()
    };

    use std::collections::HashMap;
    let mut genome_outputs: HashMap<String, GenomeOutput> = HashMap::new();
    for genome_path in genomes {
        let (completeness, contamination) = quality_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("Quality data not found for genome: {}", genome_path));
        let (r5s, r16s, r23s) = rrna_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("rRNA data not found for genome: {}", genome_path));
        let trnas = trna_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("tRNA data not found for genome: {}", genome_path));
        let mimag_quality = if completeness < 50.0 || contamination >= 10.0 {
            "Low quality"
        } else if completeness <= 90.0
            || contamination >= 5.0
            || r5s == 0
            || r16s == 0
            || r23s == 0
            || trnas < 18
        {
            "Medium quality"
        } else {
            "High quality"
        };
        genome_outputs.insert(
            genome_path.to_string(),
            GenomeOutput {
                completeness,
                contamination,
                r5s,
                r16s,
                r23s,
                trnas,
                mimag_quality: mimag_quality.to_string(),
            },
        );
    }
    Ok(genome_outputs)
}

/// Parse two-column TSV mapping genome names to Barrnap GFF files
fn parse_barrnap_gff_list(
    list_path: &str,
) -> Result<HashMap<String, (usize, usize, usize)>, String> {
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

        let (r5s, r16s, r23s) = crate::barrnap::parse_rrna_types(gff_path);
        rrna_cache.insert(genome_path, (r5s, r16s, r23s));
    }
    Ok(rrna_cache)
}

/// Parse two-column TSV mapping genome names to tRNAscan-SE output files  
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
