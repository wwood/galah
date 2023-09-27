use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use skani::chain;
use skani::file_io;
use skani::params::*;

pub struct SkaniPreclusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
}

impl PreclusterDistanceFinder for SkaniPreclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        precluster_skani(
            genome_fasta_paths,
            self.threshold,
            self.min_aligned_threshold,
        )
    }
}

fn precluster_skani(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
) -> SortedPairGenomeDistanceCache {
    let mut to_return = SortedPairGenomeDistanceCache::new();
    for (i, fasta1) in genome_fasta_paths.iter().enumerate() {
        for (j, fasta2) in genome_fasta_paths[(i + 1)..genome_fasta_paths.len()]
            .iter()
            .enumerate()
        {
            let genome_index2 = i + j + 1;
            let ani = calculate_skani(fasta1, fasta2, min_aligned_threshold);
            if ani >= threshold {
                to_return.insert((i, genome_index2), Some(ani));
            }
        }
    }
    to_return
}

pub struct SkaniClusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
}

impl ClusterDistanceFinder for SkaniClusterer {
    fn initialise(&self) {
        assert!(self.threshold > 1.0);
    }

    fn method_name(&self) -> &str {
        "skani"
    }

    fn get_ani_threshold(&self) -> f32 {
        self.threshold
    }

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32> {
        Some(calculate_skani(fasta1, fasta2, self.min_aligned_threshold))
    }
}

fn default_params(mode: Mode, min_aligned_frac: f32) -> (CommandParams, SketchParams) {
    let cmd_params = CommandParams {
        screen: false,
        screen_val: 0.00,
        mode,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        max_results: 10000000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: min_aligned_frac as f64,
        keep_refs: false,
        est_ci: false,
        learned_ani: true,
        learned_ani_cmd: false,
        detailed_out: false,
    };

    let m = 1000;
    let c = 125;
    let k = 15;
    let sketch_params = SketchParams::new(m, c, k, false, false);
    (cmd_params, sketch_params)
}

pub fn calculate_skani(fasta1: &str, fasta2: &str, min_aligned_frac: f32) -> f32 {
    //Vector of Strings
    let refs = vec![fasta1.to_string()];
    let queries = vec![fasta2.to_string()];

    let (command_params, sketch_params) = default_params(Mode::Dist, min_aligned_frac);
    let ref_sketch = &file_io::fastx_to_sketches(&refs, &sketch_params, true)[0];
    let query_sketch = &file_io::fastx_to_sketches(&queries, &sketch_params, true)[0];
    let map_params = chain::map_params_from_sketch(ref_sketch, false, &command_params);
    let ani_result = chain::chain_seeds(ref_sketch, query_sketch, map_params);

    ani_result.ani * 100.0
}
