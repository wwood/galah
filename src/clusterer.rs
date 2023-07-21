use std;
use std::sync::Mutex;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use partitions::partition_vec::PartitionVec;
use rayon::prelude::*;

/// Given a list of genomes, return them clustered. Use dashing for first pass
/// analysis, then fastani as the actual threshold.
pub fn cluster<P: PreclusterDistanceFinder, C: ClusterDistanceFinder + std::marker::Sync>(
    genomes: &[&str],
    preclusterer: &P,
    clusterer: &C,
) -> Vec<Vec<usize>> {
    clusterer.initialise();

    // Dashing all the genomes together
    let dashing_cache = preclusterer.distances(genomes);

    info!("Preclustering ..");
    let minhash_preclusters = partition_sketches(genomes, &dashing_cache);
    trace!("Found preclusters: {:?}", minhash_preclusters);

    let all_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![]);

    // Convert single linkage data structure into just a list of list of indices
    let mut preclusters: Vec<Vec<usize>> = minhash_preclusters
        .all_sets()
        .map(|cluster| {
            let mut indices: Vec<_> = cluster.map(|cluster_genome| *cluster_genome.1).collect();
            indices.sort_unstable();
            indices
        })
        .collect();

    // Sort preclusters so bigger clusters are started before smaller
    preclusters.sort_unstable_by_key(|c2| std::cmp::Reverse(c2.len()));
    debug!("After sorting, found preclusters {:?}", preclusters);
    info!(
        "Found {} preclusters. The largest contained {} genomes",
        preclusters.len(),
        preclusters[0].len()
    );

    info!("Finding representative genomes and assigning all genomes to these ..");
    preclusters
        .par_iter()
        .enumerate()
        .for_each(|(precluster_id, original_genome_indices)| {
            let precluster_dashing_cache = dashing_cache.transform_ids(original_genome_indices);
            let mut precluster_genomes = vec![];
            for original_genome_index in original_genome_indices {
                precluster_genomes.push(genomes[*original_genome_index]);
            }
            debug!(
                "Clustering pre-cluster {}, with genome indices {:?}",
                precluster_id, original_genome_indices
            );
            trace!(
                "Found precluster dashing cache {:#?}",
                precluster_dashing_cache
            );

            debug!(
                "Calculating genome representatives by dashing+fastani in precluster {} ..",
                precluster_id
            );
            let (clusters, calculated_fastanis) = clusterer
                .find_representatives(&precluster_dashing_cache, precluster_genomes.as_slice());
            debug!(
                "In precluster {}, found {} genome representatives",
                precluster_id,
                clusters.len()
            );

            debug!(
                "Assigning genomes to representatives by dashing+fastani in precluster {}..",
                precluster_id
            );
            let clusters = clusterer.find_memberships(
                &clusters,
                &precluster_dashing_cache,
                precluster_genomes.as_slice(),
                calculated_fastanis,
            );
            // Indices here are within this single linkage cluster only, so
            // here we map them back to their original indices, and then
            // insert back into the original array.
            for cluster in clusters {
                all_clusters.lock().unwrap().push(
                    cluster
                        .iter()
                        .map(|within_index| original_genome_indices[*within_index])
                        .collect(),
                );
            }
            debug!("Finished proccessing pre-cluster {}", precluster_id);
        });
    all_clusters.into_inner().unwrap()
}

/// Create sub-sets by single linkage clustering
fn partition_sketches(
    genomes: &[&str],
    dashing_cache: &SortedPairGenomeDistanceCache,
) -> PartitionVec<usize> {
    let mut to_return: PartitionVec<usize> = PartitionVec::with_capacity(genomes.len());
    for (i, _) in genomes.iter().enumerate() {
        to_return.push(i);
    }

    genomes.iter().enumerate().for_each(|(i, _)| {
        genomes[0..i].iter().enumerate().for_each(|(j, _)| {
            trace!("Testing precluster between {} and {}", i, j);
            if dashing_cache.contains_key(&(i, j)) {
                debug!(
                    "During preclustering, found a match between genomes {} and {}",
                    i, j
                );
                to_return.union(i, j)
            }
        });
    });
    to_return
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    // #[test]
    // fn test_minhash_hello_world() {
    //     init();
    //     let clusters = minhash_clusters(
    //         &["tests/data/abisko4/73.20120800_S1X.13.fna",
    //           "tests/data/abisko4/73.20120600_S2D.19.fna",
    //           "tests/data/abisko4/73.20120700_S3X.12.fna",
    //           "tests/data/abisko4/73.20110800_S2D.13.fna",
    //         ],
    //         95.0,
    //         1000,
    //         21,
    //         None,
    //     );
    //     assert_eq!(
    //         vec![vec![0,1,2,3]],
    //         clusters
    //     )
    // }

    // #[test]
    // fn test_minhash_two_clusters() {
    //     init();
    //     let clusters = minhash_clusters(
    //         &["tests/data/abisko4/73.20120800_S1X.13.fna",
    //           "tests/data/abisko4/73.20120600_S2D.19.fna",
    //           "tests/data/abisko4/73.20120700_S3X.12.fna",
    //           "tests/data/abisko4/73.20110800_S2D.13.fna",
    //         ],
    //         98.0,
    //         1000,
    //         21,
    //         None,
    //     );
    //     assert_eq!(
    //         vec![vec![0,1,3],vec![2]],
    //         clusters
    //     )
    // }

    #[test]
    fn test_minhash_fastani_hello_world() {
        init();
        let mut clusters = cluster(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21,
            },
            &crate::fastani::FastaniClusterer {
                threshold: 95.0,
                min_aligned_threshold: 0.2,
                fraglen: 3000,
            },
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 2, 3]], clusters)
    }

    #[test]
    fn test_minhash_fastani_two_clusters_same_ani() {
        init();
        let mut clusters = cluster(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21,
            },
            &crate::fastani::FastaniClusterer {
                threshold: 98.0,
                min_aligned_threshold: 0.2,
                fraglen: 3000,
            },
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 3], vec![2]], clusters)
    }

    #[test]
    fn test_minhash_fastani_two_clusters_low_minhash_ani() {
        init();
        let mut clusters = cluster(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21,
            },
            &crate::fastani::FastaniClusterer {
                threshold: 98.0,
                min_aligned_threshold: 0.2,
                fraglen: 3000,
            },
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 3], vec![2]], clusters)
    }

    // #[test]
    // fn test_minhash_skani_hello_world() {
    //     init();
    //     let mut clusters = cluster(
    //         &[
    //             "tests/data/abisko4/73.20120800_S1X.13.fna",
    //             "tests/data/abisko4/73.20120600_S2D.19.fna",
    //             "tests/data/abisko4/73.20120700_S3X.12.fna",
    //             "tests/data/abisko4/73.20110800_S2D.13.fna",
    //         ],
    //         &crate::finch::FinchPreclusterer {
    //             min_ani: 0.9,
    //             num_kmers: 1000,
    //             kmer_length: 21,
    //         },
    //         &crate::SkaniClusterer {},
    //         95.0,
    //         0.2,
    //         3000,
    //         "skani"
    //     );
    //     for cluster in clusters.iter_mut() {
    //         cluster.sort_unstable();
    //     }
    //     assert_eq!(vec![vec![0, 1, 2, 3]], clusters)
    // }
}
