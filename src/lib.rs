pub mod analyse;
pub mod analyse_argument_parsing;
pub mod barrnap;
pub mod checkm2;
pub mod cluster_argument_parsing;
pub mod cluster_validation;
pub mod clusterer;
pub mod eukcc;
pub mod external_command_checker;
pub mod fastani;
pub mod finch;
pub mod genome_info_file;
pub mod genome_stats;
pub mod isiteuk;
pub mod pixi_env;
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
    /// Returns (r5s, r16s, r23s, r18s, r28s, r58s)
    fn find_rrnas(
        &self,
        genome_path: &str,
        tmp_path: &std::path::Path,
    ) -> (usize, usize, usize, usize, usize, usize);
    fn method_name(&self) -> &str;
}

/// Biological domain of a genome as determined by marker gene analysis.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Domain {
    Bacteria,
    Archaea,
    Eukaryota,
}

impl Domain {
    /// Parse from isiteuk output column value (e.g. "d__Bacteria").
    pub fn from_isiteuk_str(s: &str) -> Option<Domain> {
        match s {
            "d__Bacteria" => Some(Domain::Bacteria),
            "d__Archaea" => Some(Domain::Archaea),
            "d__Eukaryota" => Some(Domain::Eukaryota),
            _ => None,
        }
    }

    pub fn display_name(&self) -> &str {
        match self {
            Domain::Bacteria => "Bacteria",
            Domain::Archaea => "Archaea",
            Domain::Eukaryota => "Eukaryota",
        }
    }

    pub fn barrnap_kingdom(&self) -> &str {
        match self {
            Domain::Bacteria => "bac",
            Domain::Archaea => "arc",
            Domain::Eukaryota => "fun", // barrnap 1.10+ uses "fun" for eukaryotes; "euk" is not supported
        }
    }

    pub fn trnascan_mode(&self) -> &str {
        match self {
            Domain::Bacteria => "B",
            Domain::Archaea => "A",
            Domain::Eukaryota => "E",
        }
    }

    pub fn is_prokaryote(&self) -> bool {
        matches!(self, Domain::Bacteria | Domain::Archaea)
    }
}

/// How to determine the domain of each genome.
#[derive(Debug, Clone, PartialEq)]
pub enum DomainChoice {
    /// Run isiteuk to classify each genome (default).
    Isiteuk,
    /// Treat all genomes as Bacteria.
    Bacteria,
    /// Treat all genomes as Archaea.
    Archaea,
    /// Treat all genomes as Eukaryota.
    Eukaryota,
    /// Run all domain tools for every genome and take the best result.
    All,
}

impl std::str::FromStr for DomainChoice {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "isiteuk" => Ok(DomainChoice::Isiteuk),
            "bac" => Ok(DomainChoice::Bacteria),
            "arc" => Ok(DomainChoice::Archaea),
            "euk" => Ok(DomainChoice::Eukaryota),
            "all" => Ok(DomainChoice::All),
            _ => Err(format!("Unknown domain choice: {s}")),
        }
    }
}

impl DomainChoice {
    /// Returns fixed domain list, or None when domains must be determined at runtime.
    pub fn fixed_domains(&self) -> Option<Vec<Domain>> {
        match self {
            DomainChoice::Bacteria => Some(vec![Domain::Bacteria]),
            DomainChoice::Archaea => Some(vec![Domain::Archaea]),
            DomainChoice::Eukaryota => Some(vec![Domain::Eukaryota]),
            DomainChoice::All => Some(vec![Domain::Bacteria, Domain::Archaea, Domain::Eukaryota]),
            DomainChoice::Isiteuk => None,
        }
    }
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
pub const DEFAULT_DOMAIN_CHOICE: &str = "isiteuk";
pub const DOMAIN_CHOICES: [&str; 5] = ["isiteuk", "bac", "arc", "euk", "all"];
pub const DEFAULT_ISITEUK_BACTERIA_CUTOFF: &str = "10";
pub const DEFAULT_ISITEUK_ARCHAEA_CUTOFF: &str = "10";
pub const DEFAULT_ISITEUK_EUKARYOTA_CUTOFF: &str = "20";

pub const AUTHOR: &str =
    "Ben J. Woodcroft, Centre for Microbiome Research, Queensland University of Technology";
