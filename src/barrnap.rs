use crate::Domain;
use crate::RrnaFinder;
use std::fs;
use std::path::{Path, PathBuf};

pub struct BarrnapAnalyser;

impl RrnaFinder for BarrnapAnalyser {
    /// Run barrnap for bac, arc, and euk kingdoms and return the best result.
    fn find_rrnas(
        &self,
        genome_path: &str,
        tmp_path: &std::path::Path,
    ) -> (usize, usize, usize, usize, usize, usize) {
        get_barrnap_output_for_domains(
            genome_path,
            &[Domain::Bacteria, Domain::Archaea, Domain::Eukaryota],
            tmp_path,
        )
    }

    fn method_name(&self) -> &str {
        "Barrnap"
    }
}

/// Run barrnap for the kingdoms corresponding to `domains` and return the best result.
///
/// Returns `(r5s, r16s, r23s, r18s, r28s, r58s)`.
pub fn get_barrnap_output_for_domains(
    genome_path: &str,
    domains: &[Domain],
    tmp_path: &std::path::Path,
) -> (usize, usize, usize, usize, usize, usize) {
    let mut kingdoms: Vec<&str> = domains.iter().map(|d| d.barrnap_kingdom()).collect();
    kingdoms.dedup();

    let mut best = (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    for kingdom in kingdoms {
        let gff_path = run_barrnap(genome_path, kingdom, 1, tmp_path);
        let result = parse_rrna_types_all(gff_path.to_str().unwrap());
        let total = result.0 + result.1 + result.2 + result.3 + result.4 + result.5;
        let best_total = best.0 + best.1 + best.2 + best.3 + best.4 + best.5;
        if total > best_total {
            best = result;
        }
    }
    best
}

pub fn run_barrnap(genome_path: &str, kingdom: &str, threads: usize, out_dir: &Path) -> PathBuf {
    let genome_name = Path::new(genome_path)
        .file_stem()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let gff_path = out_dir.join(format!("{genome_name}.{kingdom}.gff"));
    if gff_path.is_file() {
        info!("Using cached Barrnap output: {:?}", gff_path);
        return gff_path;
    }
    let output = crate::pixi_env::command_with_pixi_fallback("barrnap", "default")
        .args([
            "--kingdom",
            kingdom,
            "--threads",
            &threads.to_string(),
            genome_path,
        ])
        .output()
        .expect("Failed to run barrnap");

    if !output.status.success() {
        info!(
            "Barrnap run on {} failed with {}.\nstdout:\n{}\nstderr:\n{}",
            genome_path,
            output.status,
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        panic!("Barrnap did not run successfully");
    }

    fs::write(&gff_path, &output.stdout).expect("Failed to write barrnap output");
    gff_path
}

/// Minimum fraction of the expected (Rfam model) length a hit must span to count as a complete
/// rRNA gene. This is the old barrnap `--lencutoff` default.
const MIN_RRNA_LENGTH_FRACTION: f64 = 0.8;

/// Rfam model length (CLEN) for each rRNA model used by barrnap >= 1.0, keyed by Rfam accession.
fn rfam_model_length(accession: &str) -> Option<usize> {
    match accession {
        "RF00001" => Some(120),  // 5S (all kingdoms)
        "RF00002" => Some(154),  // 5.8S
        "RF00177" => Some(1533), // SSU bacteria (16S)
        "RF01959" => Some(1478), // SSU archaea (16S)
        "RF01960" => Some(1831), // SSU eukarya (18S)
        "RF02540" => Some(2987), // LSU archaea (23S)
        "RF02541" => Some(2925), // LSU bacteria (23S)
        "RF02543" => Some(3401), // LSU eukarya (28S)
        _ => None,
    }
}

/// True if the GFF record describes a truncated rRNA hit.
///
/// Barrnap >= 1.0 uses global Rfam models and no longer labels partial hits, so hits are judged by
/// length against the model length, taken from the `Dbxref=Rfam:RFxxxxx` attribute. Barrnap 0.9
/// labelled partial hits `(partial)` in `product=`. Records with neither (e.g. hand-written GFFs
/// without a `Dbxref`, or an unknown accession) are not filtered.
fn is_partial_rrna(start: &str, end: &str, attributes: &str) -> bool {
    if attributes.contains("(partial)") {
        return true;
    }
    let accession = attributes
        .split(';')
        .find_map(|kv| kv.strip_prefix("Dbxref=Rfam:"));
    let (Some(model_length), Ok(start), Ok(end)) = (
        accession.and_then(rfam_model_length),
        start.parse::<usize>(),
        end.parse::<usize>(),
    ) else {
        return false;
    };
    let hit_length = end.abs_diff(start) + 1;
    (hit_length as f64) < MIN_RRNA_LENGTH_FRACTION * model_length as f64
}

/// Parse all rRNA types (prokaryotic and eukaryotic) from a barrnap GFF file.
///
/// Partial hits are skipped, since MIMAG requires complete rRNA genes: see `is_partial_rrna`.
///
/// Returns `(r5s, r16s, r23s, r18s, r28s, r58s)`.
pub fn parse_rrna_types_all(gff_path: &str) -> (usize, usize, usize, usize, usize, usize) {
    let content = std::fs::read_to_string(gff_path).unwrap();
    let mut r5s = 0;
    let mut r16s = 0;
    let mut r23s = 0;
    let mut r18s = 0;
    let mut r28s = 0;
    let mut r58s = 0;
    for line in content.lines() {
        if line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            continue;
        }
        if is_partial_rrna(fields[3], fields[4], fields[8]) {
            continue;
        }
        if let Some(name) = fields[8].split(';').find_map(|kv| kv.strip_prefix("Name=")) {
            match name {
                "5S_rRNA" => r5s += 1,
                "16S_rRNA" => r16s += 1,
                "23S_rRNA" => r23s += 1,
                "18S_rRNA" => r18s += 1,
                "28S_rRNA" => r28s += 1,
                "5.8S_rRNA" => r58s += 1,
                _ => {}
            }
        }
    }
    (r5s, r16s, r23s, r18s, r28s, r58s)
}

/// Parse a barrnap GFF file and count only prokaryotic rRNA types (5S, 16S, 23S).
///
/// Used when reading pre-computed GFF files generated before eukaryotic support was added.
pub fn parse_rrna_types(gff_path: &str) -> (usize, usize, usize) {
    let (r5s, r16s, r23s, _, _, _) = parse_rrna_types_all(gff_path);
    (r5s, r16s, r23s)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// One barrnap 1.10-style rRNA GFF record spanning `start..=end`.
    fn record(name: &str, accession: &str, start: usize, end: usize) -> String {
        format!(
            "c1\tbarrnap:1.10.5\trRNA\t{start}\t{end}\t1e-10\t+\t.\t\
             Name={name};Alias=model;Dbxref=Rfam:{accession};product={name} ribosomal RNA\n"
        )
    }

    fn parse(gff: &str) -> (usize, usize, usize, usize, usize, usize) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("x.gff");
        std::fs::write(&path, format!("##gff-version 3\n{gff}")).unwrap();
        parse_rrna_types_all(path.to_str().unwrap())
    }

    #[test]
    fn test_complete_hits_of_every_model_are_counted() {
        let gff = [
            record("5S_rRNA", "RF00001", 1, 120),
            record("16S_rRNA", "RF00177", 1, 1533), // bacteria
            record("16S_rRNA", "RF01959", 1, 1478), // archaea
            record("23S_rRNA", "RF02541", 1, 2925), // bacteria
            record("23S_rRNA", "RF02540", 1, 2987), // archaea
            record("18S_rRNA", "RF01960", 1, 1831),
            record("28S_rRNA", "RF02543", 1, 3401),
            record("5.8S_rRNA", "RF00002", 1, 154),
        ]
        .concat();
        assert_eq!(parse(&gff), (1, 2, 2, 1, 1, 1));
    }

    #[test]
    fn test_truncated_hit_of_every_model_is_skipped() {
        let gff = [
            record("5S_rRNA", "RF00001", 1, 50),
            record("16S_rRNA", "RF00177", 1, 503),
            record("16S_rRNA", "RF01959", 1, 503),
            record("23S_rRNA", "RF02541", 1, 800),
            record("23S_rRNA", "RF02540", 1, 800),
            record("18S_rRNA", "RF01960", 1, 500),
            record("28S_rRNA", "RF02543", 1, 1000),
            record("5.8S_rRNA", "RF00002", 1, 60),
        ]
        .concat();
        assert_eq!(parse(&gff), (0, 0, 0, 0, 0, 0));
    }

    #[test]
    fn test_length_threshold_boundary() {
        // 80% of 1533 is 1226.4, so 1227 bp counts and 1226 bp does not.
        assert_eq!(parse(&record("16S_rRNA", "RF00177", 1, 1227)).1, 1);
        assert_eq!(parse(&record("16S_rRNA", "RF00177", 1, 1226)).1, 0);
        // 80% of 120 is exactly 96, so 96 bp counts and 95 bp does not.
        assert_eq!(parse(&record("5S_rRNA", "RF00001", 1, 96)).0, 1);
        assert_eq!(parse(&record("5S_rRNA", "RF00001", 1, 95)).0, 0);
    }

    #[test]
    fn test_hit_longer_than_model_is_counted() {
        assert_eq!(parse(&record("16S_rRNA", "RF00177", 1, 1700)).1, 1);
    }

    #[test]
    fn test_reverse_strand_coordinates() {
        // Coordinates given high-to-low are measured the same way.
        assert_eq!(parse(&record("16S_rRNA", "RF00177", 1500, 1)).1, 1);
        assert_eq!(parse(&record("16S_rRNA", "RF00177", 503, 1)).1, 0);
    }

    #[test]
    fn test_only_truncated_copies_are_skipped_when_mixed() {
        let gff = [
            record("16S_rRNA", "RF00177", 1, 503),
            record("16S_rRNA", "RF00177", 5000, 6520),
        ]
        .concat();
        assert_eq!(parse(&gff).1, 1);
    }

    #[test]
    fn test_partial_label_is_skipped_even_if_long() {
        // Barrnap 0.9 style: labelled partial, no Dbxref.
        let gff = "c1\tbarrnap:0.9\trRNA\t1\t1400\t1e-10\t+\t.\t\
                   Name=16S_rRNA;product=16S ribosomal RNA (partial);note=aligned only 45 percent\n";
        assert_eq!(parse(gff).1, 0);
    }

    #[test]
    fn test_records_without_model_length_are_not_filtered() {
        // No Dbxref (hand-written / mock GFFs): short hits still count.
        let no_dbxref =
            "c1\tbarrnap\trRNA\t200\t300\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n";
        assert_eq!(parse(no_dbxref).1, 1);
        // Unknown accession.
        assert_eq!(parse(&record("16S_rRNA", "RF99999", 1, 100)).1, 1);
        // Unparseable coordinates.
        let bad = "c1\tbarrnap\trRNA\tNA\tNA\t.\t+\t.\tName=16S_rRNA;Dbxref=Rfam:RF00177\n";
        assert_eq!(parse(bad).1, 1);
    }

    #[test]
    fn test_comments_short_lines_and_other_features_are_ignored() {
        let gff = "# comment\nshort\tline\n\
                   c1\tbarrnap\ttRNA\t1\t70\t.\t+\t.\tName=tRNA-Ala;Dbxref=Rfam:RF00177\n";
        assert_eq!(parse(gff), (0, 0, 0, 0, 0, 0));
    }

    #[test]
    fn test_real_barrnap_output_with_truncated_16s() {
        // Real barrnap 1.10.5 (--kingdom bac) output: the 16S is 503 bp at the end of a contig.
        let counts = parse_rrna_types_all("tests/data/partial_16s/DRP004690_co1_28.bac.gff");
        assert_eq!(counts, (1, 0, 0, 0, 0, 0));
    }
}
