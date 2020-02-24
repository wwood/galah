pub mod clusterer;
pub mod cluster_argument_parsing;
pub mod cluster_validation;
pub mod sorted_pair_genome_distance_cache;
pub mod dashing;
pub mod finch;
pub mod external_command_checker;
pub mod genome_info_file;

#[macro_use]
extern crate log;
extern crate rayon;
extern crate partitions;
extern crate clap;



use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;

pub trait PreclusterDistanceFinder {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache;
}
