pub mod cluster_argument_parsing;
pub mod cluster_validation;
pub mod clusterer;
pub mod dashing;
pub mod external_command_checker;
pub mod fastani;
pub mod finch;
pub mod genome_info_file;
pub mod genome_stats;
pub mod skani;
pub mod sorted_pair_genome_distance_cache;

#[macro_use]
extern crate log;
extern crate clap;
extern crate disjoint;
extern crate rayon;
#[macro_use]
extern crate lazy_static;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;

pub trait PreclusterDistanceFinder {
    fn distances(&self, genome_fasta_paths: &[&str]) -> SortedPairGenomeDistanceCache;
}

pub trait ClusterDistanceFinder {
    fn initialise(&self);

    fn method_name(&self) -> &str;

    fn get_ani_threshold(&self) -> f32;

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32>;
}

pub const DEFAULT_ALIGNED_FRACTION: &str = "50";
pub const DEFAULT_FRAGMENT_LENGTH: &str = "3000";
pub const DEFAULT_ANI: &str = "99";
pub const DEFAULT_PRETHRESHOLD_ANI: &str = "95";
pub const DEFAULT_QUALITY_FORMULA: &str = "Parks2020_reduced";
pub const DEFAULT_PRECLUSTER_METHOD: &str = "dashing";
pub const DEFAULT_CLUSTER_METHOD: &str = "fastani";

pub const AUTHOR: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, Queensland University of Technology";
