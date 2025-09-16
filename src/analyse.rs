use crate::QualityFinder;
use crate::RrnaFinder;
use crate::TrnaFinder;
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

pub fn analyse<Q: QualityFinder, R: RrnaFinder, T: TrnaFinder>(
    genomes: &[String],
    threads: usize,
    quality_finder: &mut Q,
    rrna_finder: &R,
    trna_finder: &T,
) -> std::collections::HashMap<String, GenomeOutput> {
    let quality_method = quality_finder.method_name();
    let rrna_method = rrna_finder.method_name();
    let trna_method = trna_finder.method_name();
    info!(
        "Running {}, {} and {} on provided genomes...",
        quality_method, rrna_method, trna_method
    );

    let tmpdir = tempdir().expect("Failed to create tempdir");
    let tmp_path = tmpdir.path();

    quality_finder.prepare_comp_cont(genomes, threads, tmp_path);

    use std::collections::HashMap;
    let mut genome_outputs: HashMap<String, GenomeOutput> = HashMap::new();
    for genome_path in genomes {
        let (completeness, contamination) = quality_finder.find_comp_cont(genome_path);
        let (r5s, r16s, r23s) = rrna_finder.find_rrnas(genome_path, tmp_path);
        let trnas = trna_finder.find_trnas(genome_path, tmp_path);
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
    genome_outputs
}
