use std;

/// Like a BTreeMap except the keys are sorted before insertion / get etc.
#[derive(Debug, PartialEq, Clone)]
pub struct SortedPairGenomeDistanceCache {
    internal: std::collections::BTreeMap<(usize, usize), Option<f32>>,
}

impl Default for SortedPairGenomeDistanceCache {
    fn default() -> Self {
        Self::new()
    }
}

impl SortedPairGenomeDistanceCache {
    pub fn new() -> SortedPairGenomeDistanceCache {
        SortedPairGenomeDistanceCache {
            internal: std::collections::BTreeMap::new(),
        }
    }

    pub fn insert(&mut self, genome_ids: (usize, usize), distance: Option<f32>) {
        if genome_ids.0 < genome_ids.1 {
            self.internal.insert((genome_ids.0, genome_ids.1), distance);
        } else {
            self.internal.insert((genome_ids.1, genome_ids.0), distance);
        }
    }

    pub fn get(&self, genome_ids: &(usize, usize)) -> Option<&Option<f32>> {
        if genome_ids.0 < genome_ids.1 {
            self.internal.get(&(genome_ids.0, genome_ids.1))
        } else {
            self.internal.get(&(genome_ids.1, genome_ids.0))
        }
    }

    pub fn contains_key(&self, genome_ids: &(usize, usize)) -> bool {
        if genome_ids.0 < genome_ids.1 {
            self.internal.contains_key(&(genome_ids.0, genome_ids.1))
        } else {
            self.internal.contains_key(&(genome_ids.1, genome_ids.0))
        }
    }

    /// Create a new cache containing a subset of distances.
    pub fn transform_ids(&self, input_ids: &[usize]) -> SortedPairGenomeDistanceCache {
        let mut to_return = SortedPairGenomeDistanceCache::new();
        for (i, genome_id1) in input_ids.iter().enumerate() {
            for (j, genome_id2) in input_ids.iter().enumerate().skip(i + 1) {
                trace!("Testing {} / {} i.e. {} / {}", i, j, genome_id1, genome_id2);
                if let Some(ani_opt) = self.get(&(*genome_id1, *genome_id2)) {
                    to_return.insert((i, j), *ani_opt)
                }
            }
        }
        to_return
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_transform_hello_world() {
        init();

        let mut cache = SortedPairGenomeDistanceCache::new();
        cache.insert((1, 2), Some(0.99));

        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {} }",
            format!("{:?}", cache.transform_ids(&[0, 3]))
        );
        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {(0, 1): Some(0.99)} }",
            format!("{:?}", cache.transform_ids(&[1, 2]))
        );
        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {} }",
            format!("{:?}", cache.transform_ids(&[1, 3]))
        );
    }

    #[test]
    fn test_transform_multiple() {
        init();

        let mut cache = SortedPairGenomeDistanceCache::new();
        cache.insert((1, 2), Some(0.99));
        cache.insert((1, 4), Some(0.98));

        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {} }",
            format!("{:?}", cache.transform_ids(&[0, 3]))
        );
        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {(0, 1): Some(0.99)} }",
            format!("{:?}", cache.transform_ids(&[1, 2]))
        );
        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {(0, 1): Some(0.98)} }",
            format!("{:?}", cache.transform_ids(&[1, 4]))
        );
        assert_eq!(
            "SortedPairGenomeDistanceCache { internal: {(0, 1): Some(0.99), (0, 2): Some(0.98)} }",
            format!("{:?}", cache.transform_ids(&[1, 2, 4]))
        );
    }
}
