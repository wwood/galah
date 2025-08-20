use std;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Mutex;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;
use crate::PreclusterDistanceFinder;

use disjoint::DisjointSetVec;
use rayon::prelude::*;

/// Given a list of genomes, return them clustered. Use preclusterer for first pass
/// analysis, then clusterer as the actual threshold.
pub fn cluster<P: PreclusterDistanceFinder, C: ClusterDistanceFinder + std::marker::Sync>(
    genomes: &[&str],
    preclusterer: &P,
    clusterer: &C,
    cluster_contigs: bool,
    contig_names: Option<&[&str]>,
    reference_genomes: Option<&[&str]>,
) -> Vec<Vec<usize>> {
    clusterer.initialise();

    let preclusterer_name = preclusterer.method_name();
    let clusterer_name = clusterer.method_name();

    info!(
        "Preclustering with {} and clustering with {}",
        preclusterer_name, clusterer_name
    );

    let mut skip_clusterer = false;
    if clusterer_name == preclusterer_name {
        info!("Preclustering and clustering methods are the same, so reusing ANI values");
        skip_clusterer = true;
    }

    if cluster_contigs {
        if ["finch", "dashing"].contains(&preclusterer_name) {
            panic!("{} does not support contig comparisons.", preclusterer_name);
        }
        info!("Clustering contigs using {} ..", preclusterer_name);
        skip_clusterer = true;
    }

    // Preclusterer all the genomes together
    let preclusterer_cache = if let Some(ref_genomes) = reference_genomes {
        // For reference-based clustering, we need to compare genomes with ref_genomes
        preclusterer.distances_with_references(genomes, ref_genomes)
    } else if cluster_contigs {
        preclusterer.distances_contigs(genomes, contig_names.unwrap())
    } else {
        preclusterer.distances(genomes)
    };

    info!("Preclustering ..");
    let single_linkage_preclusters = if cluster_contigs {
        partition_sketches(contig_names.unwrap(), &preclusterer_cache)
    } else {
        partition_sketches(genomes, &preclusterer_cache)
    };
    trace!("Found preclusters: {:?}", single_linkage_preclusters);

    let all_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![]);

    // Convert single linkage data structure into just a list of list of indices
    let mut preclusters: Vec<Vec<usize>> = single_linkage_preclusters
        .indices()
        .sets()
        .iter()
        .map(|cluster| {
            let mut indices: Vec<_> = cluster.to_vec();
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
            let precluster_preclusterer_cache =
                preclusterer_cache.transform_ids(original_genome_indices);
            let mut precluster_genomes = vec![];
            for original_genome_index in original_genome_indices {
                if !cluster_contigs {
                    precluster_genomes.push(genomes[*original_genome_index]);
                } else {
                    precluster_genomes.push(contig_names.unwrap()[*original_genome_index]);
                }
            }
            debug!(
                "Clustering pre-cluster {}, with genome indices {:?}",
                precluster_id, original_genome_indices
            );
            trace!(
                "Found precluster preclusterer cache {:#?}",
                precluster_preclusterer_cache
            );

            debug!(
                "Calculating genome representatives by {}+{} in precluster {} ..",
                preclusterer_name, clusterer_name, precluster_id
            );
            let (clusters, calculated_clusterer_anis) = find_precluster_cluster_representatives(
                clusterer,
                &precluster_preclusterer_cache,
                precluster_genomes.as_slice(),
                skip_clusterer,
            );
            debug!(
                "In precluster {}, found {} genome representatives",
                precluster_id,
                clusters.len()
            );

            debug!(
                "Assigning genomes to representatives by {}+{} in precluster {}..",
                preclusterer_name, clusterer_name, precluster_id
            );
            let clusters = find_precluster_cluster_memberships(
                clusterer,
                &clusters,
                &precluster_preclusterer_cache,
                precluster_genomes.as_slice(),
                calculated_clusterer_anis,
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

// /// Choose representatives, greedily assigning based on the min_ani threshold.
// fn find_minhash_representatives(
//     sketches: &[Sketch],
//     ani_threshold: f64)
//     -> BTreeSet<usize> {

//     let mut to_return: BTreeSet<usize> = BTreeSet::new();

//     for (i, sketch1) in sketches.iter().enumerate() {
//         let mut is_rep = true;
//         for j in &to_return {
//             let sketch2: &Sketch = &sketches[*j];
//             if distance(&sketch1.hashes, &sketch2.hashes, "", "", true)
//                 .expect("Failed to calculate distance by sketch comparison")
//                 .mashDistance
//                 <= ani_threshold {

//                 is_rep = false;
//                 break;
//             }
//         }
//         if is_rep {
//             to_return.insert(i);
//         }
//     }
//     return to_return;
// }

fn find_precluster_cluster_representatives(
    clusterer: &(impl ClusterDistanceFinder + std::marker::Sync),
    preclusterer_cache: &SortedPairGenomeDistanceCache,
    genomes: &[&str],
    skip_clusterer: bool,
) -> (BTreeSet<usize>, SortedPairGenomeDistanceCache) {
    let mut clusters_to_return: BTreeSet<usize> = BTreeSet::new();
    let mut clusterer_cache = SortedPairGenomeDistanceCache::new();

    for (i, _) in genomes.iter().enumerate() {
        // Gather a list of all genomes which pass the minhash threshold, sorted
        // so highest ANIs are first
        let mut minhash_indices_and_distances: Vec<_> = clusters_to_return
            .iter()
            .map(|j| (j, preclusterer_cache.get(&(i, *j))))
            .filter(|(_, ani_opt)| ani_opt.is_some())
            .map(|(j, ani_opt)| (j, ani_opt.unwrap()))
            .collect();
        minhash_indices_and_distances.sort_unstable_by(|a, b| a.1.partial_cmp(b.1).unwrap());
        let potential_refs: Vec<_> = minhash_indices_and_distances
            .into_iter()
            .map(|(rep_index, _)| *rep_index)
            .collect();

        // Clusterer all potential reps against the current genome
        let clusterer_anis = if skip_clusterer {
            compute_ani_from_preclusterer(
                preclusterer_cache,
                potential_refs.iter().collect::<Vec<_>>().as_slice(),
                &i,
            )
        } else {
            calculate_clusterer_many_to_one_pairwise_stop_early(
                clusterer,
                potential_refs
                    .iter()
                    .map(|ref_index| genomes[*ref_index])
                    .collect::<Vec<_>>()
                    .as_slice(),
                genomes[i],
            )
        };
        let mut is_rep = true;
        for (potential_ref, clusterer_ani) in potential_refs.into_iter().zip(clusterer_anis.iter())
        {
            debug!(
                "Read in clusterer {:?} from ref {} {} and genome {} {}",
                clusterer_ani, potential_ref, genomes[potential_ref], i, genomes[i]
            );
            // The reps always have a lower index that the genome under
            // consideration, and the cache key is in ascending order.
            if let Some(ani) = clusterer_ani {
                debug!(
                    "Inserting into cache {}/{} {:?}",
                    potential_ref, i, *clusterer_ani
                );
                if !skip_clusterer {
                    clusterer_cache.insert((potential_ref, i), *clusterer_ani);
                }
                if *ani >= clusterer.get_ani_threshold() {
                    is_rep = false
                }
            }
        }
        if is_rep {
            debug!("Genome designated representative: {} {}", i, genomes[i]);
            clusters_to_return.insert(i);
        }
    }

    // When skip_clusterer is true, return all ANI in preclusterer_cache
    // Fixes bug when transitive property is not satisfied (see test_contig_cluster_rep_bug_small)
    if skip_clusterer {
        (clusters_to_return, preclusterer_cache.clone())
    } else {
        (clusters_to_return, clusterer_cache)
    }
}

/// Calculate Clusterer ANI values, submitting each genome pair in parallel.
fn calculate_clusterer_many_to_one_pairwise(
    clusterer: &(impl ClusterDistanceFinder + std::marker::Sync),
    query_genome_paths: &[&str],
    ref_genome_path: &str,
) -> Vec<Option<f32>> {
    query_genome_paths
        .par_iter()
        .map(|query_genome| clusterer.calculate_ani(query_genome, ref_genome_path))
        .collect()
}

/// Calculate Clusterer values, submitting each genome pair in parallel until one
/// passes the ani_threshold. Then stop. Return a list of Clusterer ANI values.
/// These need to be screened to determine if any passed the threshold.
fn calculate_clusterer_many_to_one_pairwise_stop_early(
    clusterer: &(impl ClusterDistanceFinder + std::marker::Sync),
    query_genome_paths: &[&str],
    ref_genome_path: &str,
) -> Vec<Option<f32>> {
    let to_return: Mutex<Vec<Option<f32>>> = Mutex::new(vec![None; query_genome_paths.len()]);

    query_genome_paths
        .par_iter()
        .enumerate()
        .find_any(|(i, query_genome)| {
            let ani = clusterer.calculate_ani(query_genome, ref_genome_path);
            to_return.lock().unwrap()[*i] = ani;
            match ani {
                Some(real_ani) => real_ani >= clusterer.get_ani_threshold(),
                None => false,
            }
        });

    to_return.into_inner().unwrap()
}

fn compute_ani_from_preclusterer(
    preclusterer_cache: &SortedPairGenomeDistanceCache,
    query_genome_indexes: &[&usize],
    ref_genome_index: &usize,
) -> Vec<Option<f32>> {
    query_genome_indexes
        .iter()
        .map(|query_genome| {
            let ani = preclusterer_cache.get(&(**query_genome, *ref_genome_index));
            ani.copied()
        })
        .collect::<Vec<_>>()
        .into_iter()
        .flatten()
        .collect()
}

// /// For each genome (sketch) assign it to the closest representative genome:
// fn find_minhash_memberships(
//     representatives: &BTreeSet<usize>,
//     sketches: &[Sketch],
// ) -> Vec<Vec<usize>> {

//     let mut rep_to_index = BTreeMap::new();
//     for (i, rep) in representatives.iter().enumerate() {
//         rep_to_index.insert(rep, i);
//     }

//     let mut to_return: Vec<Vec<usize>> = vec![vec![]; representatives.len()];
//     for (i, sketch1) in sketches.iter().enumerate() {
//         if representatives.contains(&i) {
//             to_return[rep_to_index[&i]].push(i);
//         } else {
//             let mut best_rep_min_ani = None;
//             let mut best_rep = None;
//             for rep in representatives.iter() {
//                 let dist = distance(&sketch1.hashes, &sketches[*rep].hashes, "", "", true)
//                     .expect("Failed to calculate distance by sketch comparison")
//                     .mashDistance;
//                 if best_rep_min_ani.is_none() || dist < best_rep_min_ani.unwrap() {
//                     best_rep = Some(rep);
//                     best_rep_min_ani = Some(dist);
//                 }
//             }
//             to_return[rep_to_index[best_rep.unwrap()]].push(i);
//         }
//     }
//     return to_return;
// }

/// Given a list of representative genomes and a list of genomes, assign each
/// genome to the closest representative.
fn find_precluster_cluster_memberships(
    clusterer: &(impl ClusterDistanceFinder + std::marker::Sync),
    representatives: &BTreeSet<usize>,
    preclusterer_cache: &SortedPairGenomeDistanceCache,
    genomes: &[&str],
    calculated_clusterer_anis: SortedPairGenomeDistanceCache,
) -> Vec<Vec<usize>> {
    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    debug!(
        "Before re-assignment, clusterer cache is {:#?}",
        calculated_clusterer_anis
    );

    let to_return: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![vec![]; representatives.len()]);
    let calculated_clusterer_anis_lock = Mutex::new(calculated_clusterer_anis);

    // Push all reps first so they are at the beginning of the list.
    for i in representatives.iter() {
        to_return.lock().unwrap()[rep_to_index[&i]].push(*i)
    }

    genomes.par_iter().enumerate().for_each(|(i, _)| {
        if !representatives.contains(&i) {
            let potential_refs: Vec<&usize> = representatives
                .iter()
                .filter(|rep| {
                    if calculated_clusterer_anis_lock
                        .lock()
                        .unwrap()
                        .contains_key(&(i, **rep))
                    {
                        false // Clusterer not needed since already cached
                    } else {
                        preclusterer_cache.contains_key(&(i, **rep))
                    }
                })
                .collect();

            let new_clusterers = calculate_clusterer_many_to_one_pairwise(
                clusterer,
                &potential_refs
                    .iter()
                    .map(|ref_i| genomes[**ref_i])
                    .collect::<Vec<_>>(),
                genomes[i],
            );
            for (ref_i, ani_opt) in potential_refs.iter().zip(new_clusterers.iter()) {
                calculated_clusterer_anis_lock
                    .lock()
                    .unwrap()
                    .insert((i, **ref_i), *ani_opt);
            }

            // TODO: Is there Clusterer bugs here? Donovan/Pierre seem to think so
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let clusterer: Option<f32> = match i < *rep {
                    true => match calculated_clusterer_anis_lock
                        .lock()
                        .unwrap()
                        .get(&(i, *rep))
                    {
                        // ani could be known as f32, None for calculated
                        // but below threshold, of not calculated. If not calculated, it isn't nearby
                        Some(ani_opt) => ani_opt.as_ref().copied(),
                        None => None,
                    },
                    false => match calculated_clusterer_anis_lock
                        .lock()
                        .unwrap()
                        .get(&(*rep, i))
                    {
                        Some(ani_opt) => ani_opt.as_ref().copied(),
                        None => None,
                    },
                };
                debug!(
                    "Between rep {} {} and genome {} {}, after decaching found clusterer {:?}",
                    rep, genomes[*rep], i, genomes[i], clusterer
                );

                if clusterer.is_some()
                    && (best_rep_min_ani.is_none()
                        || clusterer.unwrap() > best_rep_min_ani.unwrap())
                {
                    best_rep = Some(rep);
                    best_rep_min_ani = clusterer;
                }
            }
            to_return.lock().unwrap()[rep_to_index[best_rep.unwrap()]].push(i);
        }
    });

    to_return.into_inner().unwrap()
}

/// Create sub-sets by single linkage clustering
fn partition_sketches(
    genomes: &[&str],
    preclusterer_cache: &SortedPairGenomeDistanceCache,
) -> DisjointSetVec<usize> {
    let mut to_return: DisjointSetVec<usize> = DisjointSetVec::with_capacity(genomes.len());
    for (i, _) in genomes.iter().enumerate() {
        to_return.push(i);
    }

    // Collect all pairs that need to be joined in parallel
    let pairs_to_join: Vec<(usize, usize)> = genomes
        .par_iter()
        .enumerate()
        .flat_map(|(i, _)| {
            (0..i).into_par_iter().filter_map(move |j| {
                trace!("Testing precluster between {} and {}", i, j);
                if preclusterer_cache.contains_key(&(i, j)) {
                    debug!(
                        "During preclustering, found a match between genomes {} and {}",
                        i, j
                    );
                    Some((i, j))
                } else {
                    None
                }
            })
        })
        .collect();

    // Join all pairs sequentially
    for (i, j) in pairs_to_join {
        to_return.join(i, j);
    }

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
            false,
            None,
            None,
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
            false,
            None,
            None,
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
            false,
            None,
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 3], vec![2]], clusters)
    }

    #[test]
    fn test_minhash_skani_hello_world() {
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
            &crate::skani::SkaniClusterer {
                threshold: 95.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
            },
            false,
            None,
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 2, 3]], clusters)
    }

    #[test]
    fn test_minhash_skani_two_clusters_same_ani() {
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
            &crate::skani::SkaniClusterer {
                threshold: 99.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
            },
            false,
            None,
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 3], vec![2]], clusters)
    }

    #[test]
    fn test_skani_skani_two_clusters_same_ani() {
        init();
        let mut clusters = cluster(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::skani::SkaniPreclusterer {
                threshold: 90.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
                threads: 1,
            },
            &crate::skani::SkaniClusterer {
                threshold: 99.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
            },
            false,
            None,
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        assert_eq!(vec![vec![0, 1, 3], vec![2]], clusters)
    }

    #[test]
    fn test_skani_skani_two_preclusters() {
        init();
        let mut clusters = cluster(
            &[
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
                "tests/data/antonio_mags/BE_RX_R2_MAG52.fna",
            ],
            &crate::skani::SkaniPreclusterer {
                threshold: 90.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
                threads: 1,
            },
            &crate::skani::SkaniClusterer {
                threshold: 99.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
            },
            false,
            None,
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        clusters.sort_unstable();
        assert_eq!(vec![vec![0, 1, 3], vec![2], vec![4]], clusters)
    }

    #[test]
    fn test_contig_cluster() {
        init();
        let mut clusters = cluster(
            &["tests/data/contigs/contigs.fna"],
            &crate::skani::SkaniPreclusterer {
                threshold: 90.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
                threads: 1,
            },
            &crate::skani::SkaniClusterer {
                threshold: 99.0,
                min_aligned_threshold: 0.2,
                small_genomes: false,
            },
            true,
            Some(&[
                "73.20110600_S2D.10_contig_13024",
                "73.20110600_S2D.10_contig_13024_2",
                "73.20110600_S2D.10_contig_50844",
                "73.20110600_S2D.10_contig_37820",
            ]),
            None,
        );
        for cluster in clusters.iter_mut() {
            cluster.sort_unstable();
        }
        clusters.sort_unstable();
        assert_eq!(vec![vec![0, 1], vec![2], vec![3]], clusters)
    }
}
