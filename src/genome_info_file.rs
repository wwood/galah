use checkm;
use std;
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct GenomeInfoGenomeQuality {
    pub completeness: f32,
    pub contamination: f32,
}
impl checkm::GenomeQuality for GenomeInfoGenomeQuality {
    fn completeness(&self) -> f32 {
        self.completeness
    }
    fn contamination(&self) -> f32 {
        self.contamination
    }
}

/// Read a genome Info file defined by dRep i.e. ["genome"(basename of .fasta
/// file of that genome), "completeness"(0-100 value for completeness of the
/// genome), "contamination"(0-100 value of the contamination of the genome)].
pub fn read_genome_info_file(
    file_path: &str,
) -> Result<checkm::CheckMResult<GenomeInfoGenomeQuality>, String> {
    let mut qualities = std::collections::BTreeMap::new();
    let rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(std::path::Path::new(file_path));
    let mut total_seen = 0usize;

    let mut parse_result =
        rdr.unwrap_or_else(|_| panic!("Failed to parse genomeInfo file {}", file_path));
    if parse_result
        .headers()
        .expect("Failed to find headers in genomeInfo file")
        != vec!["genome", "completeness", "contamination"]
    {
        return Err("Incorrect headers found in genomeInfo file".to_string());
    }

    for result in parse_result.records() {
        let res = result.expect("Parsing error in genomeInfo file");
        if res.len() != 3 {
            return Err(format!(
                "Parsing error in genomeInfo file - didn't find 3 columns in line {res:?}",
            ));
        }
        let completeness: f32 = res[1]
            .parse::<f32>()
            .expect("Error parsing completeness in checkm tab table");
        let contamination: f32 = res[2]
            .parse::<f32>()
            .expect("Error parsing contamination in checkm tab table");
        trace!(
            "For {}, found completeness {} and contamination {}",
            &res[0],
            completeness,
            contamination
        );
        match qualities.insert(
            res[0].to_string(),
            GenomeInfoGenomeQuality {
                completeness: completeness / 100.,
                contamination: contamination / 100.,
            },
        ) {
            None => {}
            Some(_) => {
                return Err(format!(
                    "The genome {} was found multiple times in the checkm file {}",
                    &res[0], file_path
                ));
            }
        };
        total_seen += 1;
    }
    debug!("Read in {} genomes from {}", total_seen, file_path);
    Ok(checkm::CheckMResult {
        genome_to_quality: qualities,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn init() {
        let _ = env_logger::builder().is_test(true).try_init();
    }

    #[test]
    fn test_hello_world() {
        init();

        let genome_info = read_genome_info_file("tests/data/set1/genomeInfo.csv").unwrap();
        let mut map = std::collections::BTreeMap::new();
        map.insert(
            "500kb".to_string(),
            GenomeInfoGenomeQuality {
                completeness: 0.5,
                contamination: 0.01,
            },
        );
        map.insert(
            "1mbp".to_string(),
            GenomeInfoGenomeQuality {
                completeness: 1.0,
                contamination: 0.,
            },
        );
        assert_eq!(genome_info.genome_to_quality, map);
    }

    #[test]
    fn test_fail_on_checkm_tab_table() {
        init();

        assert!(read_genome_info_file("tests/data/set1/checkm.tsv").is_err());
    }
}
