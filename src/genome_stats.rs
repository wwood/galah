use needletail::parse_sequence_path;
use needletail::sequence::Sequence;

#[derive(Debug, PartialEq)]
pub struct GenomeAssemblyStats {
    pub num_contigs: usize,
    pub num_ambiguous_bases: usize,
    pub n50: usize,
}

pub fn calculate_genome_stats(fasta_path: &str) -> GenomeAssemblyStats {
    let mut num_contigs = 0;
    let mut num_ambiguous = 0;
    let mut contig_lengths = vec![];
    let mut total_length = 0usize;

    parse_sequence_path(
        fasta_path,
        |_| {},
        |seq| {
            num_contigs += 1;
            let s = seq.sequence();
            contig_lengths.push(s.len());
            total_length += s.len();
            for base in s {
                if base == &b'N' || base == &b'n' {
                    num_ambiguous += 1
                }
            }
        },
    )
    .expect(&format!(
        "Failed to calculate genome statistics for file {}",
        fasta_path
    ));

    // Calculate n50
    contig_lengths.sort();
    let n50_cutoff = total_length / 2;
    let mut n50 = None;
    let mut n50_sum = 0usize;
    for length in contig_lengths {
        n50_sum += length;
        if n50_sum >= n50_cutoff {
            n50 = Some(length);
            break;
        }
    }

    GenomeAssemblyStats {
        num_contigs: num_contigs,
        num_ambiguous_bases: num_ambiguous,
        n50: n50.expect(&format!("Failed to calculate n50 from {}", fasta_path)),
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
                num_ambiguous_bases: 6506,
                n50: 8289
            },
            calculate_genome_stats(&"tests/data/abisko4/73.20110600_S2D.10.fna")
        )
    }

    #[test]
    fn test_one_contig_n50() {
        init();

        assert_eq!(
            GenomeAssemblyStats {
                num_contigs: 1,
                num_ambiguous_bases: 0,
                n50: 1_000_000
            },
            calculate_genome_stats(&"tests/data/set1/1mbp.fna")
        )
    }
}
