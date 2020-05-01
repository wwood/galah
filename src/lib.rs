pub mod cluster_argument_parsing;
pub mod cluster_validation;
pub mod clusterer;
pub mod dashing;
pub mod external_command_checker;
pub mod finch;
pub mod genome_info_file;
pub mod genome_stats;
pub mod sorted_pair_genome_distance_cache;

#[macro_use]
extern crate log;
extern crate clap;
extern crate partitions;
extern crate rayon;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;

pub trait PreclusterDistanceFinder {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache;
}
