use std;
use std::collections::{BTreeSet,BTreeMap};
use std::io::BufReader;
use std::sync::Mutex;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::PreclusterDistanceFinder;

use rayon::prelude::*;
use bird_tool_utils::command::finish_command_safely;
use partitions::partition_vec::PartitionVec;

/// Given a list of genomes, return them clustered. Use dashing for first pass
/// analysis, then fastani as the actual threshold.
pub fn cluster<P: PreclusterDistanceFinder>(
    genomes: &[&str],
    preclusterer: &P,
    fastani_threshold: f32,
    fastani_min_aligned_threshold: f32,
) -> Vec<Vec<usize>> {

    assert!(fastani_threshold > 1.0);

    // Dashing all the genomes together
    let dashing_cache = preclusterer.distances(genomes);

    info!("Preclustering ..");
    let minhash_preclusters = partition_sketches(genomes, &dashing_cache);
    trace!("Found preclusters: {:?}", minhash_preclusters);

    let all_clusters: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![]);

    // Convert single linkage data structure into just a list of list of indices
    let mut preclusters: Vec<Vec<usize>> = minhash_preclusters.all_sets().map(|cluster| {
        let mut indices: Vec<_> = cluster.map(|cluster_genome| *cluster_genome.1).collect();
        indices.sort_unstable();
        indices
    }).collect();

    // Sort preclusters so bigger clusters are started before smaller
    preclusters.sort_unstable_by( |c1, c2| c2.len().cmp(&c1.len()));
    debug!("After sorting, found preclusters {:?}", preclusters);
    info!("Found {} preclusters. The largest contained {} genomes",
        preclusters.len(), preclusters[0].len());

    info!("Finding representative genomes and assigning all genomes to these ..");
    preclusters.par_iter().enumerate().for_each( |(precluster_id, original_genome_indices)| {
        let precluster_dashing_cache = dashing_cache.transform_ids(original_genome_indices);
        let mut precluster_genomes = vec![];
        for original_genome_index in original_genome_indices {
            precluster_genomes.push(genomes[*original_genome_index]);
        }
        debug!("Clustering pre-cluster {}, with genome indices {:?}", precluster_id, original_genome_indices);
        trace!("Found precluster dashing cache {:#?}", precluster_dashing_cache);

        debug!("Calculating genome representatives by dashing+fastani in precluster {} ..",
            precluster_id);
        let (clusters, calculated_fastanis) = find_dashing_fastani_representatives(
            &precluster_dashing_cache, precluster_genomes.as_slice(), fastani_threshold, fastani_min_aligned_threshold
        );
        debug!("In precluster {}, found {} genome representatives", precluster_id, clusters.len());

        debug!("Assigning genomes to representatives by dashing+fastani in precluster {}..",
            precluster_id);
        let clusters = find_dashing_fastani_memberships(
            &clusters, &precluster_dashing_cache, precluster_genomes.as_slice(), calculated_fastanis, fastani_min_aligned_threshold
        );
        // Indices here are within this single linkage cluster only, so
        // here we map them back to their original indices, and then
        // insert back into the original array.
        for cluster in clusters {
            all_clusters.lock().unwrap().push(
                cluster.iter().map(|within_index| original_genome_indices[*within_index]).collect()
            );
        }
        debug!("Finished proccessing pre-cluster {}", precluster_id);
    });
    return all_clusters.into_inner().unwrap();
}

/// Create sub-sets by single linkage clustering
fn partition_sketches(
    genomes: &[&str],
    dashing_cache: &SortedPairGenomeDistanceCache)
    -> PartitionVec<usize> {

    let mut to_return: PartitionVec<usize> = PartitionVec::with_capacity(genomes.len());
    for (i, _) in genomes.iter().enumerate() {
        to_return.push(i);
    }

    genomes.iter().enumerate().for_each(|(i, _)| {
        genomes[0..i].iter().enumerate().for_each(|(j, _)| {
            trace!("Testing precluster between {} and {}", i, j);
            if dashing_cache.contains_key(&(i,j)) {
                debug!("During preclustering, found a match between genomes {} and {}", i, j);
                to_return.union(i,j)
            }
        });
    });
    return to_return;
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

fn find_dashing_fastani_representatives(
    dashing_cache: &SortedPairGenomeDistanceCache,
    genomes: &[&str],
    fastani_threshold: f32,
    fastani_min_aligned_threshold: f32)
    -> (BTreeSet<usize>, SortedPairGenomeDistanceCache) {

    let mut clusters_to_return: BTreeSet<usize> = BTreeSet::new();
    let mut fastani_cache = SortedPairGenomeDistanceCache::new();

    for (i, _) in genomes.iter().enumerate() {
        // Gather a list of all genomes which pass the minhash threshold, sorted
        // so highest ANIs are first
        let mut minhash_indices_and_distances: Vec<_> = clusters_to_return
            .iter()
            .map( |j| (j, dashing_cache.get(&(i,*j))))
            .filter(|(_, ani_opt)| ani_opt.is_some())
            .map(|(j, ani_opt)| (j,ani_opt.unwrap()))
            .collect();
        minhash_indices_and_distances.sort_unstable_by(|a,b| a.1.partial_cmp(&b.1).unwrap());
        let potential_refs: Vec<_> = minhash_indices_and_distances
            .into_iter()
            .map( |(rep_index,_)| *rep_index )
            .collect();

        // FastANI all potential reps against the current genome
        let fastanis = calculate_fastani_many_to_one_pairwise_stop_early(
            potential_refs.iter().map(|ref_index| genomes[*ref_index]).collect::<Vec<_>>().as_slice(),
            genomes[i],
            fastani_threshold,
            fastani_min_aligned_threshold);
        let mut is_rep = true;
        for (potential_ref, fastani) in potential_refs.into_iter().zip(fastanis.iter()) {
            debug!("Read in fastani {:?} from ref {} {} and genome {} {}", 
                fastani, potential_ref, genomes[potential_ref], i, genomes[i]);
            // The reps always have a lower index that the genome under
            // consideration, and the cache key is in ascending order.
            match fastani {
                Some(ani) => {
                    debug!("Inserting into cache {}/{} {:?}", potential_ref, i, *fastani);
                    fastani_cache.insert((potential_ref,i), *fastani);
                    if *ani >= fastani_threshold {
                        is_rep = false
                    }
                },
                None => {}
            }
        }
        if is_rep {
            debug!("Genome designated representative: {} {}", i, genomes[i]);
            clusters_to_return.insert(i);
        }
    }
    return (clusters_to_return, fastani_cache);
}

/// Calculate FastANI values, submitting each genome pair in parallel.
fn calculate_fastani_many_to_one_pairwise(
    query_genome_paths: &[&str],
    ref_genome_path: &str,
    fastani_min_aligned_threshold: f32)
    -> Vec<Option<f32>> {

    query_genome_paths.par_iter().map(|query_genome| {
        calculate_fastani(query_genome, ref_genome_path, fastani_min_aligned_threshold)
    }).collect()
}


/// Calculate FastANI values, submitting each genome pair in parallel until one
/// passes the ani_threshold. Then stop. Return a list of FastANI values.
/// These need to be screened to determine if any passed the threshold.
fn calculate_fastani_many_to_one_pairwise_stop_early(
    query_genome_paths: &[&str],
    ref_genome_path: &str,
    ani_threshold: f32,
    fastani_min_aligned_threshold: f32)
    -> Vec<Option<f32>> {

    let to_return: Mutex<Vec<Option<f32>>> = Mutex::new(vec![None;query_genome_paths.len()]);

    query_genome_paths.par_iter().enumerate().find_any(|(i,query_genome)| {
        let ani = calculate_fastani(query_genome, ref_genome_path, fastani_min_aligned_threshold);
        to_return.lock().unwrap()[*i] = ani;
        match ani {
            Some(real_ani) => real_ani >= ani_threshold,
            None => false
        }
    });

    return to_return.into_inner().unwrap()
}

pub fn calculate_fastani(fasta1: &str, fasta2: &str, fastani_min_aligned_threshold: f32) -> Option<f32> {
    let one = calculate_fastani_one_way(fasta1,fasta2,fastani_min_aligned_threshold);
    match one {
        None => return None,
        Some(first) => {
            let two = calculate_fastani_one_way(fasta2,fasta1,fastani_min_aligned_threshold);
            match two {
                None => return None,
                Some(second) => return Some(
                    if first > second {
                        first
                    } else {
                        second
                    })
            }
        }
    }
}

fn calculate_fastani_one_way(fasta1: &str, fasta2: &str, fastani_min_aligned_threshold: f32) -> Option<f32> {
    let mut cmd = std::process::Command::new("fastANI");
    cmd
        .arg("-o")
        .arg("/dev/stdout")
        .arg("--minFraction")
        .arg(&format!("{}",fastani_min_aligned_threshold))
        .arg("--query")
        .arg(&fasta1)
        .arg("--ref")
        .arg(&fasta2)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running fastANI command: {:?}", &cmd);
    let mut process = cmd.spawn().expect(&format!("Failed to spawn {}", "fastANI"));
    let stdout = process.stdout.as_mut().unwrap();
    let stdout_reader = BufReader::new(stdout);

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(stdout_reader);

    let mut to_return = None;

    for record_res in rdr.records() {
        match record_res {
            Ok(record) => {
                assert!(record.len() == 5);
                let ani: f32 = record[2].parse().expect("Failed to convert fastani ANI to float value");
                if to_return.is_some() {
                    error!("Unexpectedly found >1 result from fastANI");
                    std::process::exit(1);

                }
                to_return = Some(ani);
            },
            Err(e) => {
                error!("Error parsing fastani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    debug!("FastANI of {} against {} was {:?}", fasta1, fasta2, to_return);
    finish_command_safely(process, "FastANI");
    return to_return;
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
fn find_dashing_fastani_memberships(
    representatives: &BTreeSet<usize>,
    dashing_cache: &SortedPairGenomeDistanceCache,
    genomes: &[&str],
    calculated_fastanis: SortedPairGenomeDistanceCache,
    fastani_min_aligned_threshold: f32,
) -> Vec<Vec<usize>> {

    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    debug!("Before re-assignment, fastani cache is {:#?}", calculated_fastanis);

    let to_return: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![vec![]; representatives.len()]);
    let calculated_fastanis_lock = Mutex::new(calculated_fastanis);

    // Push all reps first so they are at the beginning of the list.
    for i in representatives.iter() {
        to_return.lock().unwrap()[rep_to_index[&i]].push(*i)
    };

    genomes.par_iter().enumerate().for_each( |(i, _)| {
        if !representatives.contains(&i) {
            let potential_refs: Vec<&usize> = representatives.into_iter().filter(|rep| {
                if calculated_fastanis_lock.lock().unwrap().contains_key(&(i,**rep)) {
                    false // FastANI not needed since already cached
                } else {
                    dashing_cache.contains_key(&(i,**rep))
                }
            }).collect();

            let new_fastanis = calculate_fastani_many_to_one_pairwise(
                &potential_refs.iter().map(|ref_i| genomes[**ref_i]).collect::<Vec<_>>(),
                genomes[i],
                fastani_min_aligned_threshold,
            );
            for (ref_i, ani_opt) in potential_refs.iter().zip(new_fastanis.iter()) {
                calculated_fastanis_lock.lock().unwrap().insert((i,**ref_i), *ani_opt);
            };

            // TODO: Is there FastANI bugs here? Donovan/Pierre seem to think so
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let fastani: Option<f32> = match i < *rep {
                    true => match calculated_fastanis_lock.lock().unwrap().get(&(i, *rep)) {
                        // ani could be known as f32, None for calculated
                        // but below threshold, of not calculated. If not calculated, it isn't nearby
                        Some(ani_opt) => match ani_opt {
                            Some(ani) => Some(*ani),
                            None => None,
                        },
                        None => None
                    },
                    false => match calculated_fastanis_lock.lock().unwrap().get(&(*rep, i)) {
                        Some(ani_opt) => match ani_opt {
                            Some(ani) => Some(*ani),
                            None => None,
                        },
                        None => None
                    }
                };
                debug!("Between rep {} {} and genome {} {}, after decaching found fastani {:?}", rep, genomes[*rep], i, genomes[i], fastani);

                if fastani.is_some() && (best_rep_min_ani.is_none() || fastani.unwrap() > best_rep_min_ani.unwrap()) {
                    best_rep = Some(rep);
                    best_rep_min_ani = fastani;
                }
            }
            to_return.lock().unwrap()[rep_to_index[best_rep.unwrap()]].push(i);
        }
    });

    return to_return.into_inner().unwrap();
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
            &["tests/data/abisko4/73.20120800_S1X.13.fna",
              "tests/data/abisko4/73.20120600_S2D.19.fna",
              "tests/data/abisko4/73.20120700_S3X.12.fna",
              "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21
            },
            95.0,
            0.2,
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,2,3]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_same_ani() {
        init();
        let mut clusters = cluster(
            &["tests/data/abisko4/73.20120800_S1X.13.fna",
              "tests/data/abisko4/73.20120600_S2D.19.fna",
              "tests/data/abisko4/73.20120700_S3X.12.fna",
              "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21
            },
            98.0,
            0.2,
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }

    #[test]
    fn test_minhash_fastani_two_clusters_low_minhash_ani() {
        init();
        let mut clusters = cluster(
            &["tests/data/abisko4/73.20120800_S1X.13.fna",
              "tests/data/abisko4/73.20120600_S2D.19.fna",
              "tests/data/abisko4/73.20120700_S3X.12.fna",
              "tests/data/abisko4/73.20110800_S2D.13.fna",
            ],
            &crate::finch::FinchPreclusterer {
                min_ani: 0.9,
                num_kmers: 1000,
                kmer_length: 21
            },
            98.0,
            0.2,
        );
        for cluster in clusters.iter_mut() { cluster.sort_unstable(); }
        assert_eq!(
            vec![vec![0,1,3],vec![2]],
            clusters
        )
    }
}
