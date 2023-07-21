use std::collections::{BTreeMap, BTreeSet};
use std::io::BufReader;
use std::sync::Mutex;

use crate::sorted_pair_genome_distance_cache::SortedPairGenomeDistanceCache;
use crate::ClusterDistanceFinder;

use bird_tool_utils::command::finish_command_safely;
use rayon::prelude::*;

pub struct FastaniClusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
    pub fraglen: u32,
}

impl ClusterDistanceFinder for FastaniClusterer {
    fn initialise(&self) {
        assert!(self.threshold > 1.0);
    }

    fn find_representatives(
        &self,
        precluster_cache: &SortedPairGenomeDistanceCache,
        genomes: &[&str],
    ) -> (BTreeSet<usize>, SortedPairGenomeDistanceCache) {
        find_dashing_fastani_representatives(
            precluster_cache,
            genomes,
            self.threshold,
            self.min_aligned_threshold,
            self.fraglen,
        )
    }

    fn find_memberships(
        &self,
        representatives: &BTreeSet<usize>,
        precluster_cache: &SortedPairGenomeDistanceCache,
        genomes: &[&str],
        calculated_fastanis: SortedPairGenomeDistanceCache,
    ) -> Vec<Vec<usize>> {
        find_dashing_fastani_memberships(
            representatives,
            precluster_cache,
            genomes,
            calculated_fastanis,
            self.min_aligned_threshold,
            self.fraglen,
        )
    }
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
    fastani_min_aligned_threshold: f32,
    fastani_fraglen: u32,
) -> (BTreeSet<usize>, SortedPairGenomeDistanceCache) {
    let mut clusters_to_return: BTreeSet<usize> = BTreeSet::new();
    let mut fastani_cache = SortedPairGenomeDistanceCache::new();

    for (i, _) in genomes.iter().enumerate() {
        // Gather a list of all genomes which pass the minhash threshold, sorted
        // so highest ANIs are first
        let mut minhash_indices_and_distances: Vec<_> = clusters_to_return
            .iter()
            .map(|j| (j, dashing_cache.get(&(i, *j))))
            .filter(|(_, ani_opt)| ani_opt.is_some())
            .map(|(j, ani_opt)| (j, ani_opt.unwrap()))
            .collect();
        minhash_indices_and_distances.sort_unstable_by(|a, b| a.1.partial_cmp(b.1).unwrap());
        let potential_refs: Vec<_> = minhash_indices_and_distances
            .into_iter()
            .map(|(rep_index, _)| *rep_index)
            .collect();

        // FastANI all potential reps against the current genome
        let fastanis = calculate_fastani_many_to_one_pairwise_stop_early(
            potential_refs
                .iter()
                .map(|ref_index| genomes[*ref_index])
                .collect::<Vec<_>>()
                .as_slice(),
            genomes[i],
            fastani_threshold,
            fastani_min_aligned_threshold,
            fastani_fraglen,
        );
        let mut is_rep = true;
        for (potential_ref, fastani) in potential_refs.into_iter().zip(fastanis.iter()) {
            debug!(
                "Read in fastani {:?} from ref {} {} and genome {} {}",
                fastani, potential_ref, genomes[potential_ref], i, genomes[i]
            );
            // The reps always have a lower index that the genome under
            // consideration, and the cache key is in ascending order.
            match fastani {
                Some(ani) => {
                    debug!(
                        "Inserting into cache {}/{} {:?}",
                        potential_ref, i, *fastani
                    );
                    fastani_cache.insert((potential_ref, i), *fastani);
                    if *ani >= fastani_threshold {
                        is_rep = false
                    }
                }
                None => {}
            }
        }
        if is_rep {
            debug!("Genome designated representative: {} {}", i, genomes[i]);
            clusters_to_return.insert(i);
        }
    }
    (clusters_to_return, fastani_cache)
}

/// Calculate FastANI values, submitting each genome pair in parallel.
fn calculate_fastani_many_to_one_pairwise(
    query_genome_paths: &[&str],
    ref_genome_path: &str,
    fastani_min_aligned_threshold: f32,
    fastani_fraglen: u32,
) -> Vec<Option<f32>> {
    query_genome_paths
        .par_iter()
        .map(|query_genome| {
            calculate_fastani(
                query_genome,
                ref_genome_path,
                fastani_min_aligned_threshold,
                fastani_fraglen,
            )
        })
        .collect()
}

/// Calculate FastANI values, submitting each genome pair in parallel until one
/// passes the ani_threshold. Then stop. Return a list of FastANI values.
/// These need to be screened to determine if any passed the threshold.
fn calculate_fastani_many_to_one_pairwise_stop_early(
    query_genome_paths: &[&str],
    ref_genome_path: &str,
    ani_threshold: f32,
    fastani_min_aligned_threshold: f32,
    fastani_fraglen: u32,
) -> Vec<Option<f32>> {
    let to_return: Mutex<Vec<Option<f32>>> = Mutex::new(vec![None; query_genome_paths.len()]);

    query_genome_paths
        .par_iter()
        .enumerate()
        .find_any(|(i, query_genome)| {
            let ani = calculate_fastani(
                query_genome,
                ref_genome_path,
                fastani_min_aligned_threshold,
                fastani_fraglen,
            );
            to_return.lock().unwrap()[*i] = ani;
            match ani {
                Some(real_ani) => real_ani >= ani_threshold,
                None => false,
            }
        });

    to_return.into_inner().unwrap()
}

pub fn calculate_fastani(
    fasta1: &str,
    fasta2: &str,
    fastani_min_aligned_threshold: f32,
    fastani_fraglen: u32,
) -> Option<f32> {
    let one = calculate_fastani_one_way(
        fasta1,
        fasta2,
        //fastani_min_aligned_threshold,
        fastani_fraglen,
    );
    match one {
        None => None,
        Some(first) => {
            let two = calculate_fastani_one_way(
                fasta2,
                fasta1,
                //fastani_min_aligned_threshold,
                fastani_fraglen,
            );
            match two {
                None => None,
                Some(second) => {
                    // Calculate min aligned based on the fragment counts, not the genome length counts as fastANI currently does. See https://github.com/wwood/galah/issues/7
                    if first.fragments_matching as f32 / first.fragments_total as f32
                        >= fastani_min_aligned_threshold
                        || second.fragments_matching as f32 / second.fragments_total as f32
                            >= fastani_min_aligned_threshold
                    {
                        Some(if first.ani > second.ani {
                            first.ani
                        } else {
                            second.ani
                        })
                    } else {
                        None
                    }
                }
            }
        }
    }
}

#[derive(Debug)]
struct FastaniMatch {
    ani: f32,
    fragments_matching: u32,
    fragments_total: u32,
}

fn calculate_fastani_one_way(
    fasta1: &str,
    fasta2: &str,
    //fastani_min_aligned_threshold: f32,
    fastani_fraglen: u32,
) -> Option<FastaniMatch> {
    let mut cmd = std::process::Command::new("fastANI");
    cmd.arg("-o")
        .arg("/dev/stdout")
        // .arg("--minFraction")
        // .arg(&format!("{}", fastani_min_aligned_threshold))
        .arg("--fragLen")
        .arg(&format!("{}", fastani_fraglen))
        .arg("--query")
        .arg(fasta1)
        .arg("--ref")
        .arg(fasta2)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped());
    debug!("Running fastANI command: {:?}", &cmd);
    let mut process = cmd
        .spawn()
        .unwrap_or_else(|_| panic!("Failed to spawn {}", "fastANI"));
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
                let ani: f32 = record[2]
                    .parse()
                    .expect("Failed to convert fastani ANI to float value");
                let fragments_matching: u32 = record[3]
                    .parse()
                    .expect("Failed to convert fastani fragment count 1 to integer");
                let fragments_total: u32 = record[4]
                    .parse()
                    .expect("Failed to convert fastani fragment count 2 to integer");
                if to_return.is_some() {
                    error!("Unexpectedly found >1 result from fastANI");
                    std::process::exit(1);
                }
                to_return = Some(FastaniMatch {
                    ani,
                    fragments_matching,
                    fragments_total,
                });
            }
            Err(e) => {
                error!("Error parsing fastani output: {}", e);
                std::process::exit(1);
            }
        }
    }
    debug!(
        "FastANI of {} against {} was {:?}",
        fasta1, fasta2, to_return
    );
    finish_command_safely(process, "FastANI");
    to_return
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
    fastani_fraglen: u32,
) -> Vec<Vec<usize>> {
    let mut rep_to_index = BTreeMap::new();
    for (i, rep) in representatives.iter().enumerate() {
        rep_to_index.insert(rep, i);
    }

    debug!(
        "Before re-assignment, fastani cache is {:#?}",
        calculated_fastanis
    );

    let to_return: Mutex<Vec<Vec<usize>>> = Mutex::new(vec![vec![]; representatives.len()]);
    let calculated_fastanis_lock = Mutex::new(calculated_fastanis);

    // Push all reps first so they are at the beginning of the list.
    for i in representatives.iter() {
        to_return.lock().unwrap()[rep_to_index[&i]].push(*i)
    }

    genomes.par_iter().enumerate().for_each(|(i, _)| {
        if !representatives.contains(&i) {
            let potential_refs: Vec<&usize> = representatives
                .iter()
                .filter(|rep| {
                    if calculated_fastanis_lock
                        .lock()
                        .unwrap()
                        .contains_key(&(i, **rep))
                    {
                        false // FastANI not needed since already cached
                    } else {
                        dashing_cache.contains_key(&(i, **rep))
                    }
                })
                .collect();

            let new_fastanis = calculate_fastani_many_to_one_pairwise(
                &potential_refs
                    .iter()
                    .map(|ref_i| genomes[**ref_i])
                    .collect::<Vec<_>>(),
                genomes[i],
                fastani_min_aligned_threshold,
                fastani_fraglen,
            );
            for (ref_i, ani_opt) in potential_refs.iter().zip(new_fastanis.iter()) {
                calculated_fastanis_lock
                    .lock()
                    .unwrap()
                    .insert((i, **ref_i), *ani_opt);
            }

            // TODO: Is there FastANI bugs here? Donovan/Pierre seem to think so
            let mut best_rep_min_ani = None;
            let mut best_rep = None;
            for rep in representatives.iter() {
                let fastani: Option<f32> = match i < *rep {
                    true => match calculated_fastanis_lock.lock().unwrap().get(&(i, *rep)) {
                        // ani could be known as f32, None for calculated
                        // but below threshold, of not calculated. If not calculated, it isn't nearby
                        Some(ani_opt) => ani_opt.as_ref().copied(),
                        None => None,
                    },
                    false => match calculated_fastanis_lock.lock().unwrap().get(&(*rep, i)) {
                        Some(ani_opt) => ani_opt.as_ref().copied(),
                        None => None,
                    },
                };
                debug!(
                    "Between rep {} {} and genome {} {}, after decaching found fastani {:?}",
                    rep, genomes[*rep], i, genomes[i], fastani
                );

                if fastani.is_some()
                    && (best_rep_min_ani.is_none() || fastani.unwrap() > best_rep_min_ani.unwrap())
                {
                    best_rep = Some(rep);
                    best_rep_min_ani = fastani;
                }
            }
            to_return.lock().unwrap()[rep_to_index[best_rep.unwrap()]].push(i);
        }
    });

    to_return.into_inner().unwrap()
}
