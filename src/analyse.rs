use crate::RrnaFinder;
use crate::TrnaFinder;
use tempfile::tempdir;

#[derive(Debug, Clone)]
pub struct GenomeOutput {
    pub r5s: usize,
    pub r16s: usize,
    pub r23s: usize,
    pub trnas: usize,
    pub mimag_quality: String,
}

pub fn analyse<R: RrnaFinder, T: TrnaFinder>(
    genomes: &[String],
    rrna_finder: &R,
    trna_finder: &T,
) -> std::collections::HashMap<String, GenomeOutput> {
    let rrna_method = rrna_finder.method_name();
    let trna_method = trna_finder.method_name();

    info!(
        "Running {} and {} on provided genomes...",
        rrna_method, trna_method
    );
    let tmpdir = tempdir().expect("Failed to create tempdir");
    let tmp_path = tmpdir.path();

    use std::collections::HashMap;
    let mut genome_outputs: HashMap<String, GenomeOutput> = HashMap::new();
    for genome_path in genomes {
        let (r5s, r16s, r23s) = rrna_finder.find_rrnas(genome_path, tmp_path);
        let trnas = trna_finder.find_trnas(genome_path, tmp_path);
        let mimag_quality = if r5s == 0 || r16s == 0 || r23s == 0 || trnas < 18 {
            "Medium quality"
        } else {
            "High quality"
        };
        genome_outputs.insert(
            genome_path.to_string(),
            GenomeOutput {
                r5s,
                r16s,
                r23s,
                trnas,
                mimag_quality: mimag_quality.to_string(),
            },
        );
    }

    genome_outputs
}
