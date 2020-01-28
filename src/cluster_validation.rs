
//use rayon::prelude::*;
use std;
use std::process;
use crate::minhash_clusterer::calculate_fastani;

pub fn validate_clusters(clustering_file: &str, ani_threshold: f32) {
    let ani = ani_threshold*100.;

    // Read cluster file
    let clusters = read_clustering_file(clustering_file);
    info!("Read in {} clusters", clusters.len());
    debug!("Clusters were: {:#?}", clusters);

    // FastANI within each cluster - each should be within the threshold
    for cluster in &clusters {
        let rep = &cluster[0];
        for genome in cluster {
            let fastani_res = calculate_fastani(&rep, &genome);
            match fastani_res {
                Some(fastani) => {
                    if fastani >= ani {
                        debug!("FastANI between {} and {} is ok: {}", rep, genome, fastani);
                    } else {
                        error!("FastANI between {} and {} is not ok: {}", rep, genome, fastani);
                    }
                },
                None => {
                    error!("FastANI between {} and {} is not ok: fastani was too divergent", rep, genome);
                }
            }
        }
    }
    
    // FastANI between each representative - each should be further than the ANI
    let reps: Vec<&str> = clusters.iter().map(|c| c[0].as_str()).collect();
    for (i, rep1) in reps.iter().enumerate() {
        for (j, rep2) in reps.iter().enumerate() {
            if i<j {
                let fastani_res = calculate_fastani(rep1, rep2);
                match fastani_res {
                    Some(fastani) => {
                        if fastani < ani {
                            debug!("FastANI between reps {} and {} is ok: {}", rep1, rep2, fastani);
                        } else {
                            error!("FastANI between reps {} and {} is not ok: {}", rep1, rep2, fastani);
                        }
                    },
                    None => {
                        error!("FastANI between reps {} and {} is ok: fastani was too divergent", rep1, rep2);
                    }
                }
            }
        }
    }

}

fn read_clustering_file(clustering_file: &str) 
-> Vec<Vec<String>> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_reader(std::fs::File::open(clustering_file).expect("Failed to open clustering file"));

    let mut current_cluster_rep: Option<String> = None;
    let mut all_clusters = vec![];
    let mut current_cluster = vec![];

    for record_res in rdr.records() {
        let record = record_res.expect("Failed to parse clustering file.");
        if record.len() != 2 {
            error!("Unexpectedly didn't find exactly 2 fields in clustering file: {:?}", record);
            process::exit(1);
        }

        if record[0] == record[1] {
            if current_cluster_rep.is_some() {
                all_clusters.push(current_cluster)
            }
            current_cluster = vec![];
            current_cluster_rep = Some(record[0].to_string());
        }
        current_cluster.push(record[1].to_string());
    }
    if current_cluster_rep.is_some() {
        all_clusters.push(current_cluster)
    }
    return all_clusters;
}