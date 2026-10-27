use crate::QualityFinder;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::io::copy as io_copy;
use std::path::Path;
use std::process::Command;

fn build_eukcc_command() -> Command {
    if let Ok(cmd_str) = std::env::var("GALAH_EUKCC_CMD") {
        let mut parts = cmd_str.split_whitespace();
        let prog = parts.next().expect("GALAH_EUKCC_CMD must not be empty");
        let mut cmd = Command::new(prog);
        cmd.args(parts);
        return cmd;
    }
    if crate::pixi_env::is_on_path("eukcc") {
        return Command::new("eukcc");
    }
    info!("eukcc not found on PATH, falling back to bundled pixi manifest");
    let manifest = crate::pixi_env::galah_manifest_path();
    let mut cmd = Command::new("pixi");
    cmd.args([
        "run",
        "--manifest-path",
        manifest.to_str().unwrap(),
        "-e",
        "eukcc",
        "eukcc",
    ]);
    cmd
}

pub struct EukccAnalyser {
    pub comp_cont_cache: HashMap<String, (f64, f64)>,
    pub database_path: String,
}

impl EukccAnalyser {
    pub fn new(database_path: String) -> Self {
        Self {
            comp_cont_cache: HashMap::new(),
            database_path,
        }
    }
}

impl QualityFinder for EukccAnalyser {
    fn prepare_comp_cont(&mut self, genome_paths: &[String], threads: usize, tmp_path: &Path) {
        info!(
            "Running EukCC on {} eukaryotic genomes...",
            genome_paths.len()
        );
        for genome_path in genome_paths {
            let (comp, cont) =
                run_eukcc_single(genome_path, threads, tmp_path, &self.database_path);
            self.comp_cont_cache
                .insert(genome_path.clone(), (comp, cont));
        }
    }

    fn find_comp_cont(&self, genome_path: &str) -> (f64, f64) {
        self.comp_cont_cache
            .get(genome_path)
            .copied()
            .unwrap_or_else(|| panic!("Genome path not found in EukCC results: {}", genome_path))
    }

    fn method_name(&self) -> &str {
        "EukCC"
    }
}

fn run_eukcc_single(
    genome_path: &str,
    threads: usize,
    tmp_path: &Path,
    database_path: &str,
) -> (f64, f64) {
    let stem1 = Path::new(genome_path).file_stem().unwrap();
    let genome_name = if genome_path.ends_with(".gz") {
        Path::new(stem1)
            .file_stem()
            .unwrap_or(stem1)
            .to_string_lossy()
            .into_owned()
    } else {
        stem1.to_string_lossy().into_owned()
    };
    let out_dir = tmp_path.join(format!("eukcc_{}", genome_name));
    // Do not pre-create out_dir — EukCC creates it itself. If the directory already exists
    // EukCC treats it as a previous run and may resume/skip steps, causing missing output.
    // EukCC has written its per-genome quality table as either `eukcc.tsv` (older versions)
    // or `eukcc.csv` (EukCC >= 2.1); the content is tab-separated in both cases.
    let tsv_path = out_dir.join("eukcc.tsv");
    let csv_path = out_dir.join("eukcc.csv");
    if let Some(existing) = [&tsv_path, &csv_path].iter().find(|p| p.is_file()) {
        info!("Using cached EukCC output: {:?}", existing);
        return parse_eukcc_tsv(existing.to_str().unwrap(), genome_path);
    }

    // EukCC does not support gzipped input; decompress beside out_dir (not inside it)
    // so EukCC's output directory stays clean.
    let decompressed_path;
    let effective_path: &str = if genome_path.ends_with(".gz") {
        let dest = tmp_path.join(format!("{}.fna", genome_name));
        let mut decoder = MultiGzDecoder::new(
            std::fs::File::open(genome_path)
                .unwrap_or_else(|e| panic!("Failed to open {}: {}", genome_path, e)),
        );
        let mut out = std::fs::File::create(&dest)
            .unwrap_or_else(|e| panic!("Failed to create {:?}: {}", dest, e));
        io_copy(&mut decoder, &mut out)
            .unwrap_or_else(|e| panic!("Failed to decompress {}: {}", genome_path, e));
        decompressed_path = dest.to_string_lossy().into_owned();
        &decompressed_path
    } else {
        genome_path
    };

    let mut cmd = build_eukcc_command();
    cmd.args([
        "single",
        effective_path,
        "--out",
        out_dir.to_str().unwrap(),
        "--threads",
        &threads.to_string(),
    ]);
    if !database_path.is_empty() {
        cmd.args(["--db", database_path]);
    }

    let output = cmd.output().expect("Failed to run EukCC");
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        match no_estimate_result(&stderr, &tsv_path, &csv_path, genome_path) {
            Some(result) => return result,
            None => {
                info!(
                    "EukCC run on {} failed with {}.\nstdout:\n{}\nstderr:\n{}",
                    genome_path,
                    output.status,
                    String::from_utf8_lossy(&output.stdout),
                    stderr
                );
                panic!("EukCC did not run successfully");
            }
        }
    }

    match [&tsv_path, &csv_path].iter().find(|p| p.is_file()) {
        Some(p) => parse_eukcc_tsv(p.to_str().unwrap(), genome_path),
        None => {
            warn!(
                "EukCC did not produce output {} or {} for {}. Returning (0.0, 0.0).",
                tsv_path.display(),
                csv_path.display(),
                genome_path
            );
            (0.0, 0.0)
        }
    }
}

/// EukCC exits non-zero when it cannot estimate quality, e.g. when no placement marker
/// genes are found, or when some markers are found but not enough to identify a suitable
/// marker gene set (see https://github.com/wwood/galah/issues/90). Newer versions write an
/// NA result table first (exit 201), while EukCC 2.1.2 exits 1 without writing any output.
/// Treat both as (0.0, 0.0) rather than a failed run. Returns None if the failure looks like
/// a genuine error.
fn no_estimate_result(
    stderr: &str,
    tsv_path: &Path,
    csv_path: &Path,
    genome_path: &str,
) -> Option<(f64, f64)> {
    if let Some(existing) = [tsv_path, csv_path].iter().find(|p| p.is_file()) {
        warn!(
            "EukCC could not estimate quality for {}. Returning (0.0, 0.0).",
            genome_path
        );
        return Some(parse_eukcc_tsv(existing.to_str().unwrap(), genome_path));
    }
    if stderr.contains("No placement marker genes found")
        || stderr.contains("No marker gene set could be found")
    {
        warn!(
            "EukCC found no suitable marker gene set for {}. Returning (0.0, 0.0).",
            genome_path
        );
        return Some((0.0, 0.0));
    }
    None
}

fn parse_eukcc_tsv(tsv_path: &str, genome_path: &str) -> (f64, f64) {
    let content = match std::fs::read_to_string(tsv_path) {
        Ok(c) => c,
        Err(e) => {
            warn!(
                "EukCC did not produce output {} for {}: {}. Returning (0.0, 0.0).",
                tsv_path, genome_path, e
            );
            return (0.0, 0.0);
        }
    };

    let mut comp_idx: Option<usize> = None;
    let mut cont_idx: Option<usize> = None;

    for line in content.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if comp_idx.is_none() {
            comp_idx = cols.iter().position(|&c| c == "completeness");
            cont_idx = cols.iter().position(|&c| c == "contamination");
            continue;
        }
        if let (Some(ci), Some(coi)) = (comp_idx, cont_idx) {
            if cols.len() > ci.max(coi) {
                let comp: f64 = cols[ci].parse().unwrap_or(0.0);
                let cont: f64 = cols[coi].parse().unwrap_or(0.0);
                return (comp, cont);
            }
        }
    }

    warn!("Could not parse EukCC quality data for {}", genome_path);
    (0.0, 0.0)
}

/// Parse a pre-computed merged EukCC TSV (with `fasta` column) into a quality cache.
pub fn parse_eukcc_quality_report(
    report_path: &str,
    genome_paths: &[String],
) -> Result<HashMap<String, (f64, f64)>, String> {
    let content = std::fs::read_to_string(report_path)
        .map_err(|e| format!("Failed to read EukCC quality report {}: {}", report_path, e))?;

    let mut cache: HashMap<String, (f64, f64)> = HashMap::new();
    let mut fasta_idx: Option<usize> = None;
    let mut comp_idx: Option<usize> = None;
    let mut cont_idx: Option<usize> = None;

    for line in content.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if fasta_idx.is_none() {
            fasta_idx = cols.iter().position(|&c| c == "fasta");
            comp_idx = cols.iter().position(|&c| c == "completeness");
            cont_idx = cols.iter().position(|&c| c == "contamination");
            continue;
        }
        if let (Some(fi), Some(ci), Some(coi)) = (fasta_idx, comp_idx, cont_idx) {
            let max_idx = fi.max(ci).max(coi);
            if cols.len() > max_idx {
                let fasta = cols[fi];
                let comp: f64 = cols[ci].parse().unwrap_or(0.0);
                let cont: f64 = cols[coi].parse().unwrap_or(0.0);

                let matched = genome_paths.iter().find(|p| {
                    p.as_str() == fasta || Path::new(p).file_stem() == Path::new(fasta).file_stem()
                });
                if let Some(gp) = matched {
                    cache.insert(gp.clone(), (comp, cont));
                }
            }
        }
    }

    Ok(cache)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_estimate_result_with_na_table() {
        // EukCC >= 2.1.3 writes an NA table before exiting 201
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv_path = tmpdir.path().join("eukcc.tsv");
        let csv_path = tmpdir.path().join("eukcc.csv");
        std::fs::write(
            &tsv_path,
            "fasta\tcompleteness\tcontamination\tncbi_lng\ngenome.fasta\tNA\tNA\tNA\n",
        )
        .unwrap();
        assert_eq!(
            no_estimate_result("", &tsv_path, &csv_path, "genome.fasta"),
            Some((0.0, 0.0))
        );
    }

    #[test]
    fn test_no_estimate_result_no_placement_markers_without_output() {
        // EukCC 2.1.2 exits 1 without writing a table
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv_path = tmpdir.path().join("eukcc.tsv");
        let csv_path = tmpdir.path().join("eukcc.csv");
        let stderr = "15-10-2026 14:34:04:  Searching for marker genes in base database\n\
                      15-10-2026 14:34:04:  No placement marker genes found for genome_metaeuk_cleaned.faa.\n";
        assert_eq!(
            no_estimate_result(stderr, &tsv_path, &csv_path, "genome.fasta"),
            Some((0.0, 0.0))
        );
    }

    #[test]
    fn test_no_estimate_result_no_suitable_marker_set_without_output() {
        // https://github.com/wwood/galah/issues/90: EukCC 2.1.2 finds some marker genes
        // (unlike issue #83's zero-marker case) but can't assemble a suitable set from them,
        // and exits 1 without writing a table.
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv_path = tmpdir.path().join("eukcc.tsv");
        let csv_path = tmpdir.path().join("eukcc.csv");
        let stderr = "Found 1 marker genes, placing them in the tree using epa-ng\n\
                      23-10-2026 12:20:56:  Automatically locating best SCMG set\n\
                      23-10-2026 12:20:58:  Could not identify a single suitable marker gene set\n\
                      23-10-2026 12:20:58:  No marker gene set could be found with these settings\n";
        assert_eq!(
            no_estimate_result(stderr, &tsv_path, &csv_path, "genome.fasta"),
            Some((0.0, 0.0))
        );
    }

    #[test]
    fn test_no_estimate_result_genuine_failure() {
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv_path = tmpdir.path().join("eukcc.tsv");
        let csv_path = tmpdir.path().join("eukcc.csv");
        assert_eq!(
            no_estimate_result(
                "No database was provided via --db or EUKCC2_DB env variable",
                &tsv_path,
                &csv_path,
                "genome.fasta"
            ),
            None
        );
    }
}
