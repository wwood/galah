use crate::Domain;
use flate2::read::MultiGzDecoder;
use std::collections::HashMap;

use std::io::Write as IoWrite;
use std::path::Path;
use std::process::Command;

fn build_isiteuk_command() -> Command {
    if let Ok(cmd_str) = std::env::var("GALAH_ISITEUK_CMD") {
        let mut parts = cmd_str.split_whitespace();
        let prog = parts.next().expect("GALAH_ISITEUK_CMD must not be empty");
        let mut cmd = Command::new(prog);
        cmd.args(parts);
        return cmd;
    }
    if crate::pixi_env::is_on_path("isiteuk") {
        return Command::new("isiteuk");
    }
    info!("isiteuk not found on PATH, falling back to bundled pixi manifest");
    let manifest = crate::pixi_env::galah_manifest_path();
    let mut cmd = Command::new("pixi");
    cmd.args([
        "run",
        "--manifest-path",
        manifest.to_str().unwrap(),
        "-e",
        "isiteuk",
        "isiteuk",
    ]);
    cmd
}

pub struct IsiTeukAnalyser {
    pub metapackage_path: String,
    pub bacteria_cutoff: f64,
    pub archaea_cutoff: f64,
    pub eukaryota_cutoff: f64,
}

impl IsiTeukAnalyser {
    pub fn new(
        metapackage_path: String,
        bacteria_cutoff: f64,
        archaea_cutoff: f64,
        eukaryota_cutoff: f64,
    ) -> Self {
        Self {
            metapackage_path,
            bacteria_cutoff,
            archaea_cutoff,
            eukaryota_cutoff,
        }
    }

    /// Run isiteuk on all genomes and return per-genome domain assignments.
    pub fn classify_genomes(
        &self,
        genome_paths: &[String],
        threads: usize,
        tmp_path: &Path,
    ) -> HashMap<String, Vec<Domain>> {
        let output_path = tmp_path.join("isiteuk_output.tsv");
        let genome_list_path = tmp_path.join("isiteuk_genomes.txt");
        let genomes_dir = tmp_path.join("isiteuk_genomes");

        // Build effective path → original mapping (compute decompressed paths without creating files).
        let mut effective_paths: Vec<String> = Vec::with_capacity(genome_paths.len());
        let mut effective_to_original: HashMap<String, String> = HashMap::new();
        for fasta in genome_paths {
            if fasta.ends_with(".gz") {
                let stem1 = Path::new(fasta).file_stem().unwrap();
                let stem2 = Path::new(stem1).file_stem().unwrap_or(stem1);
                let effective = genomes_dir
                    .join(format!("{}.fna", stem2.to_string_lossy()))
                    .to_string_lossy()
                    .to_string();
                effective_to_original.insert(effective.clone(), fasta.clone());
                effective_paths.push(effective);
            } else {
                effective_to_original.insert(fasta.clone(), fasta.clone());
                effective_paths.push(fasta.clone());
            }
        }

        // If a cached output exists, parse and return it without re-running isiteuk.
        if output_path.is_file() {
            info!("Using cached isiteuk output: {:?}", output_path);
            return parse_isiteuk_tsv(
                output_path.to_str().unwrap(),
                &effective_paths,
                self.bacteria_cutoff,
                self.archaea_cutoff,
                self.eukaryota_cutoff,
            )
            .into_iter()
            .map(|(k, v)| (effective_to_original.get(&k).cloned().unwrap_or(k), v))
            .collect();
        }

        // Decompress .gz genomes into genomes_dir for the isiteuk run.
        std::fs::create_dir_all(&genomes_dir).expect("Failed to create isiteuk genomes dir");
        for fasta in genome_paths {
            if fasta.ends_with(".gz") {
                let effective = effective_to_original
                    .iter()
                    .find(|(_, v)| *v == fasta)
                    .map(|(k, _)| k.clone())
                    .unwrap();
                let dest = Path::new(&effective);
                if !dest.is_file() {
                    let mut decoder = MultiGzDecoder::new(
                        std::fs::File::open(fasta)
                            .unwrap_or_else(|e| panic!("Failed to open {}: {}", fasta, e)),
                    );
                    let mut out = std::fs::File::create(dest)
                        .unwrap_or_else(|e| panic!("Failed to create {:?}: {}", dest, e));
                    std::io::copy(&mut decoder, &mut out)
                        .unwrap_or_else(|e| panic!("Failed to decompress {}: {}", fasta, e));
                }
            }
        }

        let mut f =
            std::fs::File::create(&genome_list_path).expect("Failed to create isiteuk genome list");
        for p in &effective_paths {
            writeln!(f, "{}", p).expect("Failed to write genome path to isiteuk list");
        }

        info!("Running isiteuk on {} genomes...", genome_paths.len());
        let mut cmd = build_isiteuk_command();
        cmd.args([
            "process",
            "--output",
            output_path.to_str().unwrap(),
            "--genome-list",
            genome_list_path.to_str().unwrap(),
            "--threads",
            &threads.to_string(),
        ]);
        if !self.metapackage_path.is_empty() {
            cmd.args(["--metapackage", &self.metapackage_path]);
        }

        let output = cmd.output().expect("Failed to run isiteuk");
        if !output.status.success() {
            info!(
                "isiteuk failed with {}.\nstdout:\n{}\nstderr:\n{}",
                output.status,
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            panic!("isiteuk did not run successfully");
        }

        parse_isiteuk_tsv(
            output_path.to_str().unwrap(),
            &effective_paths,
            self.bacteria_cutoff,
            self.archaea_cutoff,
            self.eukaryota_cutoff,
        )
        .into_iter()
        .map(|(k, v)| (effective_to_original.get(&k).cloned().unwrap_or(k), v))
        .collect()
    }
}

/// Derive the stem that isiteuk uses in its output from a genome file path.
///
/// isiteuk strips all extensions from the filename and emits just the base name.
/// For `.fna.gz` we strip two extensions; for `.fna` / `.fa` etc. we strip one.
/// We apply this only to the *path* side — never call this on the TSV genome column,
/// which already IS a bare stem (and may itself contain dots, e.g. "GCF_002008365.1_genomic").
fn genome_path_stem(path: &str) -> String {
    let p = Path::new(path);
    let fname = p.file_name().unwrap_or(p.as_os_str());
    let after_gz = if path.ends_with(".gz") {
        Path::new(fname).file_stem().unwrap_or(fname)
    } else {
        fname
    };
    Path::new(after_gz)
        .file_stem()
        .unwrap_or(after_gz)
        .to_string_lossy()
        .into_owned()
}

/// Parse an isiteuk output TSV and return per-genome domain assignments.
///
/// A domain is assigned when its `num_in_target_domain` meets the per-domain cutoff.
/// Genomes with no domain above any cutoff get an empty vec; callers fall back to all domains.
pub fn parse_isiteuk_tsv(
    path: &str,
    genome_paths: &[String],
    bacteria_cutoff: f64,
    archaea_cutoff: f64,
    eukaryota_cutoff: f64,
) -> HashMap<String, Vec<Domain>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read isiteuk output {}: {}", path, e));

    let mut result: HashMap<String, Vec<Domain>> = HashMap::new();

    for line in content.lines() {
        if line.starts_with("genome\t") {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let genome_col = cols[0];
        let domain_str = cols[1];
        let num_in: f64 = cols[2].parse().unwrap_or(0.0);

        let domain = match Domain::from_isiteuk_str(domain_str) {
            Some(d) => d,
            None => continue,
        };

        let cutoff = match &domain {
            Domain::Bacteria => bacteria_cutoff,
            Domain::Archaea => archaea_cutoff,
            Domain::Eukaryota => eukaryota_cutoff,
        };
        if num_in < cutoff {
            continue;
        }

        let matched = genome_paths
            .iter()
            .find(|p| p.as_str() == genome_col || genome_path_stem(p.as_str()) == genome_col);
        if let Some(gp) = matched {
            result.entry(gp.clone()).or_default().push(domain);
        }
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_isiteuk_tsv(dir: &std::path::Path, rows: &[(&str, &str, f64)]) -> String {
        let path = dir.join("isiteuk.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(
            f,
            "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain"
        )
        .unwrap();
        for (genome, domain, count) in rows {
            writeln!(f, "{genome}\t{domain}\t{count}\t0").unwrap();
        }
        path.to_string_lossy().to_string()
    }

    #[test]
    fn test_genome_path_stem() {
        assert_eq!(
            genome_path_stem("tests/data/domain_examples/GCF_002008365.1_genomic.fna.gz"),
            "GCF_002008365.1_genomic",
            ".fna.gz double extension"
        );
        assert_eq!(
            genome_path_stem("/tmp/isiteuk_genomes/GCF_002008365.1_genomic.fna"),
            "GCF_002008365.1_genomic",
            ".fna single extension"
        );
        assert_eq!(
            genome_path_stem("GCF_002008365.1_genomic.fna"),
            "GCF_002008365.1_genomic",
            "bare filename with extension"
        );
    }

    #[test]
    fn test_gzip_decompression_reads_all_concatenated_members() {
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Read;

        let mut encoder1 = GzEncoder::new(Vec::new(), Compression::default());
        encoder1.write_all(b">seq1\nACGT\n").unwrap();
        let member1 = encoder1.finish().unwrap();

        let mut encoder2 = GzEncoder::new(Vec::new(), Compression::default());
        encoder2.write_all(b">seq2\nTTTT\n").unwrap();
        let member2 = encoder2.finish().unwrap();

        let mut concatenated = member1;
        concatenated.extend(member2);

        let mut decoded = String::new();
        MultiGzDecoder::new(&concatenated[..])
            .read_to_string(&mut decoded)
            .unwrap();
        assert_eq!(decoded, ">seq1\nACGT\n>seq2\nTTTT\n");
    }

    #[test]
    fn test_parse_isiteuk_tsv_fna_gz_path() {
        // isiteuk outputs stem-only ("GCF_002008365.1_genomic"), but genome_paths
        // contains the original .fna.gz path. Matching must strip both extensions.
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv = write_isiteuk_tsv(
            tmpdir.path(),
            &[("GCF_002008365.1_genomic", "d__Bacteria", 34.0)],
        );
        let paths = vec!["tests/data/domain_examples/GCF_002008365.1_genomic.fna.gz".to_string()];
        let result = parse_isiteuk_tsv(&tsv, &paths, 10.0, 10.0, 20.0);
        assert_eq!(result.len(), 1);
        assert!(
            result.contains_key(&paths[0]),
            "expected key {:?}, got keys: {:?}",
            paths[0],
            result.keys().collect::<Vec<_>>()
        );
        assert_eq!(result[&paths[0]], vec![crate::Domain::Bacteria]);
    }

    #[test]
    fn test_parse_isiteuk_tsv_decompressed_fna_path() {
        // When classify_genomes decompresses to a .fna temp path, matching must
        // still find the genome via the single-extension file_stem.
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv = write_isiteuk_tsv(
            tmpdir.path(),
            &[("GCF_002008365.1_genomic", "d__Bacteria", 34.0)],
        );
        let decompressed = "/tmp/isiteuk_genomes/GCF_002008365.1_genomic.fna".to_string();
        let paths = vec![decompressed.clone()];
        let result = parse_isiteuk_tsv(&tsv, &paths, 10.0, 10.0, 20.0);
        assert_eq!(result.len(), 1);
        assert!(result.contains_key(&decompressed));
    }

    #[test]
    fn test_parse_isiteuk_tsv_cutoff_filtering() {
        let tmpdir = tempfile::tempdir().unwrap();
        let tsv = write_isiteuk_tsv(
            tmpdir.path(),
            &[
                ("genome_a", "d__Bacteria", 5.0),  // below cutoff of 10
                ("genome_b", "d__Bacteria", 15.0), // above cutoff
            ],
        );
        let paths = vec!["genome_a".to_string(), "genome_b".to_string()];
        let result = parse_isiteuk_tsv(&tsv, &paths, 10.0, 10.0, 20.0);
        assert!(
            !result.contains_key("genome_a"),
            "below-cutoff genome should not appear"
        );
        assert!(result.contains_key("genome_b"));
    }
}
