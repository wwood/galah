extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    use std::env;
    use std::fs;
    use std::path::Path;
    use tempfile::tempdir;

    fn setup_mock_bin(dir: &Path, rrna_5s: usize, rrna_16s: usize, rrna_23s: usize, trnas: usize) {
        // Barrnap outputs to stdout
        let mut barrnap_script = String::from("#!/bin/bash\n");
        if rrna_5s > 0 {
            barrnap_script.push_str("echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t1\t100\t.\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA'\n");
        }
        if rrna_16s > 0 {
            barrnap_script.push_str("echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t200\t300\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA'\n");
        }
        if rrna_23s > 0 {
            barrnap_script.push_str("echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t400\t500\t.\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA'\n");
        }
        let barrnap = dir.join("barrnap");
        fs::write(&barrnap, barrnap_script).unwrap();

        // tRNAscan-SE writes to a file specified by -o
        let common_trnas = [
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
            "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Fake",
        ];
        let mut trnascan_script = String::from("#!/bin/bash\n");
        trnascan_script.push_str("out=\"\"\n");
        trnascan_script.push_str("while [[ $# -gt 0 ]]; do\n");
        trnascan_script.push_str("  case $1 in\n");
        trnascan_script.push_str("    -o) out=$2; shift 2;;\n");
        trnascan_script.push_str("    *) shift;;\n");
        trnascan_script.push_str("  esac\n");
        trnascan_script.push_str("done\n");

        trnascan_script.push_str("echo -e 'Sequence                      \t\ttRNA \tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      ' > \"$out\"\n");
        trnascan_script.push_str("echo -e 'Name                          \ttRNA #\tBegin\tEnd  \tType\tCodon\tBegin\tEnd\tScore\tNote' >> \"$out\"\n");
        trnascan_script.push_str("echo -e '--------                      \t------\t-----\t------\t----\t-----\t-----\t----\t------\t------' >> \"$out\"\n");
        for trna in common_trnas.iter().take(trnas) {
            trnascan_script.push_str(&format!(
                "echo -e 'mock_contig\t1\t101\t200\t{trna}\tGCC\t0\t0\t20.0\tNote' >> \"$out\"\n"
            ));
        }
        let trnascan = dir.join("tRNAscan-SE");
        fs::write(&trnascan, trnascan_script).unwrap();

        for script in [&barrnap, &trnascan] {
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(script)
                .status();
        }
    }

    #[test]
    fn test_analyse_real() {
        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/set1/set1.1.fna",
                "tests/data/set1/set1.2.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/set1.1.fna\t0\t0\t0\t0\tMedium quality\n\
            tests/data/set1/set1.2.fna\t0\t0\t0\t0\tMedium quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t1\t1\t19\tHigh quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t1\t1\t1\t19\tHigh quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(tmpdir.path(), 1, 1, 1, 20);
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t1\t1\t20\tHigh quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_fake() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(tmpdir.path(), 1, 1, 1, 21);
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t1\t1\t20\tHigh quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_lower() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(tmpdir.path(), 1, 0, 0, 15);
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t0\t0\t15\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_no_16s() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(tmpdir.path(), 1, 0, 1, 20);
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t0\t1\t20\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_insufficient_trnas() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(tmpdir.path(), 1, 1, 1, 16);
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t1\t1\t1\t16\tMedium quality\n")
            .unwrap();
    }
}
