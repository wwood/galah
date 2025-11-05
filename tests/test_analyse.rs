extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    use std::env;
    use std::fs;
    use std::path::Path;
    use tempfile::tempdir;

    fn setup_mock_bin(
        dir: &Path,
        genome: String,
        completeness: f64,
        contamination: f64,
        rrna_5s: usize,
        rrna_16s: usize,
        rrna_23s: usize,
        trnas: usize,
    ) {
        // CheckM2 writes to quality_report.tsv within the directory specified by -o
        let mut checkm2_script = String::from("#!/bin/bash\n");
        checkm2_script.push_str("out=\"\"\n");
        checkm2_script.push_str("while [[ $# -gt 0 ]]; do\n");
        checkm2_script.push_str("  case $1 in\n");
        checkm2_script.push_str("    -o) out=$2; shift 2;;\n");
        checkm2_script.push_str("    *) shift;;\n");
        checkm2_script.push_str("  esac\n");
        checkm2_script.push_str("done\n");

        checkm2_script.push_str("mkdir -p \"$out\"\n");

        checkm2_script.push_str("echo -e 'Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes' > \"$out/quality_report.tsv\"\n");
        checkm2_script.push_str(&format!(
            "echo -e '{genome}\t{completeness}\t{contamination}\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3609865470852\t355151\t0.33\t446\t75\t24150\tNone' >> \"$out/quality_report.tsv\"\n"
        ));
        let checkm2 = dir.join("checkm2");
        fs::write(&checkm2, checkm2_script).unwrap();

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

        for script in [&checkm2, &barrnap, &trnascan] {
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(script)
                .status();
        }
    }

    #[test]
    #[ignore]
    fn test_analyse_real() {
        let checkm2_db_path = std::env::var("CHECKM2DB")
            .expect("CHECKM2DB environment variable must be set to run this test");
        println!("Using CheckM2 database at {}", checkm2_db_path);

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t6.35\t0.67\t0\t0\t0\t0\tLow quality\n\
            tests/data/set1/500kb.fna\t4.08\t0.02\t0\t0\t0\t0\tLow quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t82.17\t0.00\t1\t1\t1\t19\tMedium quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t84.95\t0.03\t1\t1\t1\t20\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            1,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t1\t1\t20\tHigh quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_fake() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            1,
            1,
            21,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t1\t1\t20\tHigh quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_lower() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            0,
            0,
            15,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t0\t0\t15\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_no_16s() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            0,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t0\t1\t20\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_insufficient_trnas() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            1,
            1,
            16,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t1\t1\t16\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_insufficient_completeness() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            89.9,
            2.0,
            1,
            1,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t89.90\t2.00\t1\t1\t1\t20\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_over_contamination() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            5.1,
            1,
            1,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t5.10\t1\t1\t1\t20\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_low_completeness() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            49.0,
            2.0,
            1,
            1,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t49.00\t2.00\t1\t1\t1\t20\tLow quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_high_contamination() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            11.0,
            1,
            1,
            1,
            20,
        );
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
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t11.00\t1\t1\t1\t20\tLow quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_with_checkm2_quality_report() {
        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--checkm2-quality-report",
                "tests/data/analyse_file_inputs/checkm2_quality_report.tsv",
                "--barrnap-gff-list",
                "tests/data/analyse_file_inputs/barrnap_gff_list.tsv",
                "--trnascan-out-list",
                "tests/data/analyse_file_inputs/trnascan_out_list.tsv",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t95.50\t1.20\t1\t1\t1\t19\tHigh quality\n\
            tests/data/set1/500kb.fna\t68.37\t2.91\t0\t1\t1\t10\tMedium quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.50\t1.20\t0\t0\t1\t1\tMedium quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t95.37\t2.91\t0\t0\t0\t0\tMedium quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_with_checkm_tab_table() {
        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--checkm-tab-table",
                "tests/data/analyse_file_inputs/checkm_tab_table.tsv",
                "--barrnap-gff-list",
                "tests/data/analyse_file_inputs/barrnap_gff_list.tsv",
                "--trnascan-out-list",
                "tests/data/analyse_file_inputs/trnascan_out_list.tsv",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t95.50\t1.20\t1\t1\t1\t19\tHigh quality\n\
            tests/data/set1/500kb.fna\t58.37\t12.91\t0\t1\t1\t10\tLow quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t48.37\t1.20\t0\t0\t1\t1\tLow quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t38.37\t2.91\t0\t0\t0\t0\tLow quality\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_with_quality_report_output() {
        let tmpdir = tempdir().unwrap();
        setup_mock_bin(
            tmpdir.path(),
            String::from("73.20120800_S1D.21"),
            95.0,
            2.0,
            1,
            1,
            1,
            20,
        );
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        let output_dir = tempdir().unwrap();
        let quality_report_path = output_dir.path().join("quality_report.tsv");

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--output-quality-report",
                quality_report_path.to_str().unwrap(),
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t1\t1\t20\tHigh quality\n")
            .unwrap();

        // Check that the quality report was written in CheckM2 format
        let quality_report_contents = std::fs::read_to_string(&quality_report_path).unwrap();
        assert!(quality_report_contents.contains("Name\tCompleteness\tContamination"));
        assert!(quality_report_contents.contains("73.20120800_S1D.21\t95\t2\t"));
    }
}
