use std::io::BufReader;

use crate::ClusterDistanceFinder;

use bird_tool_utils::command::finish_command_safely;

pub struct FastaniClusterer {
    pub threshold: f32,
    pub min_aligned_threshold: f32,
    pub fraglen: u32,
}

impl ClusterDistanceFinder for FastaniClusterer {
    fn initialise(&self) {
        assert!(self.threshold > 1.0);
    }

    fn method_name(&self) -> &str {
        "FastANI"
    }

    fn get_ani_threshold(&self) -> f32 {
        self.threshold
    }

    fn calculate_ani(&self, fasta1: &str, fasta2: &str) -> Option<f32> {
        calculate_fastani(fasta1, fasta2, self.min_aligned_threshold, self.fraglen)
    }
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
        .arg(format!("{fastani_fraglen}"))
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
    finish_command_safely(process, "FastANI")
        .wait()
        .expect("Unexpected wait failure outside bird_tool_utils for FastANI");
    to_return
}
