use std;

/// Like a BTreeMap except the keys are sorted before insertion / get etc.
#[derive(Debug)]
pub struct SortedPairGenomeDistanceCache {
    internal: std::collections::BTreeMap<(usize, usize), Option<f32>>
}

impl SortedPairGenomeDistanceCache {
    pub fn new() -> SortedPairGenomeDistanceCache {
        SortedPairGenomeDistanceCache {
            internal: std::collections::BTreeMap::new()
        }
    }

    pub fn insert(
        &mut self,
        genome_ids: (usize, usize),
        distance: Option<f32>) {

        if genome_ids.0 < genome_ids.1 {
            self.internal.insert((genome_ids.0,genome_ids.1), distance);
        } else {
            self.internal.insert((genome_ids.1,genome_ids.0), distance);
        }
    }

    pub fn get(
        &self,
        genome_ids: &(usize, usize))
    -> Option<&Option<f32>> {
        if genome_ids.0 < genome_ids.1 {
            self.internal.get(&(genome_ids.0,genome_ids.1))
        } else {
            self.internal.get(&(genome_ids.1,genome_ids.0))
        }
    }

    pub fn contains_key(
        &self,
        genome_ids: &(usize, usize))
    -> bool {
        if genome_ids.0 < genome_ids.1 {
            self.internal.contains_key(&(genome_ids.0,genome_ids.1))
        } else {
            self.internal.contains_key(&(genome_ids.1,genome_ids.0))
        }
    }
}