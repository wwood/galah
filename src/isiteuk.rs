use crate::Domain;
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

        let mut f =
            std::fs::File::create(&genome_list_path).expect("Failed to create isiteuk genome list");
        for p in genome_paths {
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
            genome_paths,
            self.bacteria_cutoff,
            self.archaea_cutoff,
            self.eukaryota_cutoff,
        )
    }
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

        let matched = genome_paths.iter().find(|p| {
            p.as_str() == genome_col
                || Path::new(p).file_stem() == Path::new(genome_col).file_stem()
        });
        if let Some(gp) = matched {
            result.entry(gp.clone()).or_default().push(domain);
        }
    }

    result
}
