use std;

/// Like a HashMap except the keys are sorted before insertion / get etc.
#[derive(Debug, PartialEq, Clone)]
pub struct SortedPairGenomeDistanceCache {
    internal: std::collections::HashMap<(usize, usize), Option<f32>>,
}

impl Default for SortedPairGenomeDistanceCache {
    fn default() -> Self {
        Self::new()
    }
}

impl SortedPairGenomeDistanceCache {
    pub fn new() -> SortedPairGenomeDistanceCache {
        SortedPairGenomeDistanceCache {
            internal: std::collections::HashMap::new(),
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

    /// All stored genome pairs, each as (lower index, higher index).
    pub fn keys(&self) -> impl Iterator<Item = &(usize, usize)> {
        self.internal.keys()
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

    /// HashMap iteration order is not deterministic, so compare contents via a
    /// sorted Vec rather than relying on Debug output order.
    fn sorted_entries(cache: &SortedPairGenomeDistanceCache) -> Vec<((usize, usize), Option<f32>)> {
        let mut entries: Vec<_> = cache.internal.iter().map(|(k, v)| (*k, *v)).collect();
        entries.sort_unstable_by_key(|(k, _)| *k);
        entries
    }

    #[test]
    fn test_transform_hello_world() {
        init();

        let mut cache = SortedPairGenomeDistanceCache::new();
        cache.insert((1, 2), Some(0.99));

        assert_eq!(
            Vec::<((usize, usize), Option<f32>)>::new(),
            sorted_entries(&cache.transform_ids(&[0, 3]))
        );
        assert_eq!(
            vec![((0, 1), Some(0.99))],
            sorted_entries(&cache.transform_ids(&[1, 2]))
        );
        assert_eq!(
            Vec::<((usize, usize), Option<f32>)>::new(),
            sorted_entries(&cache.transform_ids(&[1, 3]))
        );
    }

    #[test]
    fn test_transform_multiple() {
        init();

        let mut cache = SortedPairGenomeDistanceCache::new();
        cache.insert((1, 2), Some(0.99));
        cache.insert((1, 4), Some(0.98));

        assert_eq!(
            Vec::<((usize, usize), Option<f32>)>::new(),
            sorted_entries(&cache.transform_ids(&[0, 3]))
        );
        assert_eq!(
            vec![((0, 1), Some(0.99))],
            sorted_entries(&cache.transform_ids(&[1, 2]))
        );
        assert_eq!(
            vec![((0, 1), Some(0.98))],
            sorted_entries(&cache.transform_ids(&[1, 4]))
        );
        assert_eq!(
            vec![((0, 1), Some(0.99)), ((0, 2), Some(0.98))],
            sorted_entries(&cache.transform_ids(&[1, 2, 4]))
        );
    }
}
