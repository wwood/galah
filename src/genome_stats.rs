use needletail::parse_sequence_path;
use needletail::sequence::Sequence;

#[derive(Debug,PartialEq)]
pub struct GenomeAssemblyStats {
    pub num_contigs: usize,
    pub num_ambiguous_bases: usize,
}

pub fn calculate_genome_stats(fasta_path: &str) -> GenomeAssemblyStats {
    let mut num_contigs = 0;
    let mut num_ambiguous = 0;

    parse_sequence_path(
        fasta_path,
        |_| {},
        |seq| {
            num_contigs += 1;
            for base in seq.sequence() {
                if base == &b'N' || base == &b'n' {
                    num_ambiguous += 1
                }
            }
        }
    ).expect(
        &format!("Failed to calculate genome statistics for file {}", fasta_path));

    GenomeAssemblyStats {
        num_contigs: num_contigs,
        num_ambiguous_bases: num_ambiguous
    }
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

        assert_eq!(
            GenomeAssemblyStats {
                num_contigs: 161,
                num_ambiguous_bases: 6506
            },
            calculate_genome_stats(&"tests/data/abisko4/73.20110600_S2D.10.fna")
        )
    }
}