use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use rayon::prelude::*;

use skani::chain;
use skani::file_io;
use skani::params::*;
use skani::screen;

use concurrent_queue::ConcurrentQueue;

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

    fn method_name(&self) -> &str {
        "skani"
    }
}

fn precluster_skani(
    genome_fasta_paths: &[&str],
    threshold: f32,
    min_aligned_threshold: f32,
) -> SortedPairGenomeDistanceCache {
    let (command_params, sketch_params) = default_params(Mode::Dist, min_aligned_threshold);

    debug!("Sketching genomes with skani ..");
    let fasta_strings = genome_fasta_paths
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<String>>();
    // Note that sketches is now shuffled!
    let sketches = &file_io::fastx_to_sketches(&fasta_strings, &sketch_params, true);

    // Right now implemented by parallel collection into a queue, and then
    // reprocessed into a BTreeMap. Could potentially be made more efficient by
    // direct collection into a concurrent BTreeMap, but eh for now.
    let queue = ConcurrentQueue::unbounded();

    // Screen genomes before ANI calculation
    let kmer_to_sketch = screen::kmer_to_sketch_from_refs(sketches);

    info!("Calculating ANI from skani sketches ..");
    sketches.par_iter().enumerate().for_each(|(i, ref_sketch)| {
        let ref_sketch_i = &sketches[i];
        let screened_refs = screen::screen_refs(
            0.80,
            &kmer_to_sketch,
            ref_sketch_i,
            &sketch_params,
            sketches,
        );
        debug!(
            "{} has {} refs passing screening.",
            ref_sketch_i.file_name,
            screened_refs.len()
        );

        screened_refs.into_par_iter().for_each(|j| {
            if j > i {
                let map_params = chain::map_params_from_sketch(ref_sketch, false, &command_params);
                let ani_result = chain::chain_seeds(ref_sketch, &sketches[j], map_params);
                let ani = ani_result.ani * 100.;
                let ref_name = ref_sketch.file_name.clone();
                let query_name = sketches[j].file_name.clone();

                debug!("Pushing ANI result for {} and {}", ref_name, query_name);
                let ref_index = genome_fasta_paths
                    .iter()
                    .position(|&r| r == ref_name)
                    .unwrap();
                let query_index = genome_fasta_paths
                    .iter()
                    .position(|&r| r == query_name)
                    .unwrap();
                if ani >= threshold {
                    queue
                        .push((ref_index, query_index, ani))
                        .expect("Failed to push to queue during preclustering");
                }
            }
        });
    });
    debug!("Converting ANI results into sparse cache ..");

    let mut to_return = SortedPairGenomeDistanceCache::new();
    while let Ok((i, j, ani)) = queue.pop() {
        to_return.insert((i, j), Some(ani));
    }
    debug!("Finished skani.");

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

    fn calculate_ani_contigs(&self, fasta1: &str) -> Option<f32> {
        Some(calculate_skani_contigs(fasta1, self.min_aligned_threshold))
    }
}

fn default_params(mode: Mode, min_aligned_frac: f32, cluster_contigs: bool) -> (CommandParams, SketchParams) {
    let cmd_params = CommandParams {
        screen: true,
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
        max_results: 10000000, // for Triange usize::MAX,
        individual_contig_q: cluster_contigs,
        individual_contig_r: cluster_contigs,
        min_aligned_frac: min_aligned_frac as f64,
        keep_refs: false,
        est_ci: false,
        learned_ani: true,
        learned_ani_cmd: false,
        detailed_out: false,
        // distance: false,
        // rescue_small: true,
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

    let (command_params, sketch_params) = default_params(Mode::Dist, min_aligned_frac, false);
    let ref_sketch = &file_io::fastx_to_sketches(&refs, &sketch_params, true)[0];
    let query_sketch = &file_io::fastx_to_sketches(&queries, &sketch_params, true)[0];
    let map_params = chain::map_params_from_sketch(ref_sketch, false, &command_params);
    let ani_result = chain::chain_seeds(ref_sketch, query_sketch, map_params);

    ani_result.ani * 100.0
}

pub fn calculate_skani_contigs(fasta1: &str, min_aligned_frac: f32) -> f32 {
    //Vector of Strings
    let refs = vec![fasta1.to_string()];

    let (command_params, sketch_params) = default_params(Mode::Dist, min_aligned_frac, true);
    let ref_sketch = &file_io::fastx_to_sketches(&refs, &sketch_params, true)[0];
    let map_params = chain::map_params_from_sketch(ref_sketch, false, &command_params);
    let ani_result = chain::chain_seeds(ref_sketch, map_params);

    ani_result.ani * 100.0
}
