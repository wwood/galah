use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::PreclusterDistanceFinder;

pub struct FinchPreclusterer {
    /// Fraction, not percentage
    pub min_ani: f32,
    pub num_kmers: usize,
    pub kmer_length: u8,
}

impl PreclusterDistanceFinder for FinchPreclusterer {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache {
        distances(
            genome_fasta_paths,
            self.min_ani,
            self.num_kmers,
            self.kmer_length,
        )
    }

    fn method_name(&self) -> &str {
        "finch"
    }
}

pub fn distances(
    genome_fasta_paths: &[&str],
    min_ani: f32,
    num_kmers: usize,
    kmer_length: u8,
) -> SortedPairGenomeDistanceCache {
    // Hash all the files
    let sketch_params = finch::sketch_schemes::SketchParams::Mash {
        kmers_to_sketch: num_kmers,
        final_size: num_kmers,
        no_strict: true, // Possibly not right.
        kmer_length,
        hash_seed: 0,
    };
    let filters = finch::filtering::FilterParams {
        filter_on: Some(false),
        abun_filter: (None, None),
        err_filter: 0.,
        strand_filter: 0.,
    };
    info!("Sketching MinHash representations of each genome with finch ..");
    let sketches_result = finch::sketch_files(genome_fasta_paths, &sketch_params, &filters);
    info!("Finished sketching genomes");

    let sketches = sketches_result.expect("Failed to sketch genomes with finch");

    let mut to_return = SortedPairGenomeDistanceCache::new();
    for (i, sketch1) in sketches.iter().enumerate() {
        for (j, sketch2) in sketches[(i + 1)..sketches.len()].iter().enumerate() {
            let genome_index2 = i + j + 1;
            let distance = 1.0
                - finch::distance::distance(sketch1, sketch2, false)
                    .unwrap_or_else(|_| {
                        panic!(
                            "Failed to compare finch sketches {} and {}",
                            i, genome_index2
                        )
                    })
                    .mash_distance;
            debug!(
                "Comparing {} and {}, distance {}",
                genome_fasta_paths[i], genome_fasta_paths[genome_index2], distance
            );
            if distance >= min_ani as f64 {
                to_return.insert((i, genome_index2), Some(distance as f32));
            }
        }
    }
    to_return
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_hello_world() {
        init();

        let distances1 = distances(
            &["tests/data/set1/1mbp.fna", "tests/data/set1/500kb.fna"],
            0.9,
            1000,
            21,
        );
        let mut expected1 = SortedPairGenomeDistanceCache::new();
        expected1.insert((0, 1), Some(0.9808188));
        assert_eq!(expected1, distances1);

        let distances2 = distances(
            &["tests/data/set1/1mbp.fna", "tests/data/set1/500kb.fna"],
            0.99,
            1000,
            21,
        );
        let expected2 = SortedPairGenomeDistanceCache::new();
        assert_eq!(expected2, distances2);
    }
}
