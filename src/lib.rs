pub mod analyse;
pub mod analyse_argument_parsing;
pub mod barrnap;
pub mod checkm2;
pub mod cluster_argument_parsing;
pub mod cluster_validation;
pub mod clusterer;
pub mod external_command_checker;
pub mod fastani;
pub mod finch;
pub mod genome_info_file;
pub mod genome_stats;
pub mod process;
pub mod process_argument_parsing;
pub mod skani;
pub mod sorted_pair_genome_distance_cache;
pub mod trnascan;

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

    fn distances_contigs(
        &self,
        genome_fasta_paths: &[&str],
        contig_names: &[&str],
    ) -> SortedPairGenomeDistanceCache;

    fn distances_with_references(
        &self,
        genome_fasta_paths: &[&str],
        reference_genomes: &[&str],
    ) -> SortedPairGenomeDistanceCache;

    fn method_name(&self) -> &str;
}

pub trait ClusterDistanceFinder {
    fn initialise(&self);

    fn method_name(&self) -> &str;

    fn get_ani_threshold(&self) -> f32;

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32>;
}

pub trait QualityFinder {
    fn prepare_comp_cont(
        &mut self,
        genome_paths: &[String],
        threads: usize,
        tmp_path: &std::path::Path,
    );
    fn find_comp_cont(&self, genome_path: &str) -> (f64, f64);
    fn method_name(&self) -> &str;
}

pub trait TrnaFinder {
    fn find_trnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> usize;
    fn method_name(&self) -> &str;
}

pub trait RrnaFinder {
    fn find_rrnas(&self, genome_path: &str, tmp_path: &std::path::Path) -> (usize, usize, usize);
    fn method_name(&self) -> &str;
}

pub const DEFAULT_ALIGNED_FRACTION: &str = "15";
pub const DEFAULT_FRAGMENT_LENGTH: &str = "3000";
pub const DEFAULT_ANI: &str = "95";
pub const DEFAULT_PRETHRESHOLD_ANI: &str = "90";
pub const DEFAULT_QUALITY_FORMULA: &str = "Parks2020_reduced";
pub const DEFAULT_PRECLUSTER_METHOD: &str = "skani";
pub const PRECLUSTER_METHODS: [&str; 2] = ["skani", "finch"];
pub const DEFAULT_CLUSTER_METHOD: &str = "skani";
pub const CLUSTER_METHODS: [&str; 2] = ["skani", "fastani"];
pub const DEFAULT_QUALITY_METHOD: &str = "checkm2";
pub const QUALITY_METHODS: [&str; 1] = ["checkm2"];
pub const DEFAULT_RRNA_METHOD: &str = "barrnap";
pub const RRNA_METHODS: [&str; 1] = ["barrnap"];
pub const DEFAULT_TRNA_METHOD: &str = "trnascan";
pub const TRNA_METHODS: [&str; 1] = ["trnascan"];

pub const AUTHOR: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, Queensland University of Technology";
