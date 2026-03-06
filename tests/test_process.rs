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
        // genome, completeness, contamination, rrna_5s, rrna_16s, rrna_23s, trnas
        genomes: &[(String, f64, f64, usize, usize, usize, usize)],
    ) {
        // CheckM2 mock: write quality_report.tsv for all provided genomes
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

        for (name, comp, cont, _, _, _, _) in genomes.iter() {
            checkm2_script.push_str(&format!(
                "echo -e '{name}\t{comp}\t{cont}\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3609865470852\t355151\t0.33\t446\t75\t24150\tNone' >> \"$out/quality_report.tsv\"\n",
                name = name,
                comp = comp,
                cont = cont
            ));
        }
        let checkm2 = dir.join("checkm2");
        fs::write(&checkm2, checkm2_script).unwrap();

        // Barrnap mock: inspect input fasta basename and emit rRNA lines according to the matching genome config
        let mut barrnap_script = String::from("#!/bin/bash\n");
        barrnap_script.push_str("infile=\"\"\n");
        barrnap_script.push_str("while [[ $# -gt 0 ]]; do\n");
        barrnap_script.push_str("  case $1 in\n");
        barrnap_script.push_str("    -*) shift 2;;\n");
        barrnap_script.push_str("    *) if [[ -z $infile ]]; then infile=\"$1\"; fi; shift;;\n");
        barrnap_script.push_str("  esac\n");
        barrnap_script.push_str("done\n");
        barrnap_script.push_str("stem=$(basename \"$infile\"); stem=${stem%.*}\n");
        barrnap_script.push_str("case \"$stem\" in\n");
        for (name, _, _, r5, r16, r23, _) in genomes.iter() {
            barrnap_script.push_str(&format!("  {name})\n", name = name));
            if *r5 > 0 {
                barrnap_script.push_str("    echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t1\t100\t.\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA'\n");
            }
            if *r16 > 0 {
                barrnap_script.push_str("    echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t200\t300\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA'\n");
            }
            if *r23 > 0 {
                barrnap_script.push_str("    echo -e '##gff-version 3\nmock_contig\tbarrnap\trRNA\t400\t500\t.\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA'\n");
            }
            barrnap_script.push_str("    ;;\n");
        }
        barrnap_script.push_str("  *) ;;\nesac\n");
        let barrnap = dir.join("barrnap");
        fs::write(&barrnap, barrnap_script).unwrap();

        // tRNAscan-SE writes to a file specified by -o
        let common_trnas = [
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
            "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Fake",
        ];
        let mut trnascan_script = String::from("#!/bin/bash\n");
        trnascan_script.push_str("out=\"\"\n");
        trnascan_script.push_str("infile=\"\"\n");
        trnascan_script.push_str("while [[ $# -gt 0 ]]; do\n");
        trnascan_script.push_str("  case $1 in\n");
        trnascan_script.push_str("    -o) out=$2; shift 2;;\n");
        trnascan_script.push_str("    -*) shift;;\n");
        trnascan_script.push_str("    *) if [[ -z $infile ]]; then infile=\"$1\"; fi; shift;;\n");
        trnascan_script.push_str("  esac\n");
        trnascan_script.push_str("done\n");

        trnascan_script.push_str("echo -e 'Sequence                      \t\ttRNA \tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      ' > \"$out\"\n");
        trnascan_script.push_str("echo -e 'Name                          \ttRNA #\tBegin\tEnd  \tType\tCodon\tBegin\tEnd\tScore\tNote' >> \"$out\"\n");
        trnascan_script.push_str("echo -e '--------                      \t------\t-----\t------\t----\t-----\t-----\t----\t------\t------' >> \"$out\"\n");
        trnascan_script.push_str("stem=$(basename \"$infile\"); stem=${stem%.*}\n");
        trnascan_script.push_str("case \"$stem\" in\n");
        for (name, _, _, _, _, _, trnas) in genomes.iter() {
            trnascan_script.push_str(&format!("  {name})\n", name = name));
            for trna in common_trnas.iter().take(*trnas) {
                trnascan_script.push_str(&format!(
                    "    echo -e 'mock_contig\t1\t101\t200\t{trna}\tGCC\t0\t0\t20.0\tNote' >> \"$out\"\n"
                ));
            }
            trnascan_script.push_str("    ;;\n");
        }
        trnascan_script.push_str("  *) ;;\nesac\n");
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
    fn test_process_real() {
        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_quality = tmpdir.path().join("quality_report.tsv");

        let checkm2_db_path = std::env::var("CHECKM2DB")
            .expect("CHECKM2DB environment variable must be set to run this test");
        println!("Using CheckM2 database at {}", checkm2_db_path);

        Assert::main_binary()
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--output-cluster-definition",
                "/dev/stdout",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-quality-report",
                output_quality.to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .is("\
            tests/data/abisko4/73.20110800_S2M.16.fna\ttests/data/abisko4/73.20110800_S2M.16.fna\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\ttests/data/abisko4/73.20120800_S1D.21.fna\n\
            tests/data/set1/500kb.fna\ttests/data/set1/500kb.fna\n\
            tests/data/set1/500kb.fna\ttests/data/set1/1mbp.fna\n")
            .unwrap();

        // Verify analyse outputs were created properly
        assert!(output_mimag.exists());
        let content = fs::read_to_string(&output_mimag).unwrap();
        let expected = "\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t6.35\t0.67\t0\t0\t0\t0\tLow quality\n\
            tests/data/set1/500kb.fna\t4.08\t0.02\t0\t0\t0\t0\tLow quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t82.17\t0.00\t1\t1\t1\t19\tMedium quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t84.95\t0.03\t1\t1\t1\t20\tMedium quality\n";
        assert_eq!(content, expected);

        assert!(output_quality.exists());
    }

    #[test]
    fn test_process_mock() {
        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_quality = tmpdir.path().join("quality_report.tsv");

        setup_mock_bin(
            tmpdir.path(),
            &[
                (String::from("73.20120800_S1D.21"), 95.0, 2.0, 1, 1, 1, 20),
                (String::from("73.20110800_S2M.16"), 90.0, 5.0, 1, 1, 1, 20),
                (String::from("1mbp"), 85.0, 3.0, 1, 1, 1, 15),
                (String::from("500kb"), 80.0, 4.0, 0, 1, 0, 10),
            ],
        );
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[
                ("PATH", new_path),
                ("CHECKM2DB", String::from("/tmp/mockdb")),
            ])
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--output-cluster-definition",
                "/dev/stdout",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-quality-report",
                output_quality.to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .is("\
            tests/data/abisko4/73.20120800_S1D.21.fna\ttests/data/abisko4/73.20120800_S1D.21.fna\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\ttests/data/abisko4/73.20110800_S2M.16.fna\n\
            tests/data/set1/1mbp.fna\ttests/data/set1/1mbp.fna\n\
            tests/data/set1/1mbp.fna\ttests/data/set1/500kb.fna\n")
            .unwrap();

        // Verify analyse outputs were created properly
        assert!(output_mimag.exists());
        let content = fs::read_to_string(&output_mimag).unwrap();
        let expected = "\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t85.00\t3.00\t1\t1\t1\t15\tMedium quality\n\
            tests/data/set1/500kb.fna\t80.00\t4.00\t0\t1\t0\t10\tMedium quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t95.00\t2.00\t1\t1\t1\t20\tHigh quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t90.00\t5.00\t1\t1\t1\t20\tMedium quality\n";
        assert_eq!(content, expected);

        assert!(output_quality.exists());
    }

    #[test]
    fn test_process_mock_invert() {
        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_quality = tmpdir.path().join("quality_report.tsv");

        setup_mock_bin(
            tmpdir.path(),
            &[
                (String::from("1mbp"), 80.0, 4.0, 0, 1, 0, 10),
                (String::from("500kb"), 85.0, 3.0, 1, 1, 1, 15),
                (String::from("73.20120800_S1D.21"), 90.0, 5.0, 1, 1, 1, 20),
                (String::from("73.20110800_S2M.16"), 95.0, 2.0, 1, 1, 1, 20),
            ],
        );
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[
                ("PATH", new_path),
                ("CHECKM2DB", String::from("/tmp/mockdb")),
            ])
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--output-cluster-definition",
                "/dev/stdout",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-quality-report",
                output_quality.to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .is("\
            tests/data/abisko4/73.20110800_S2M.16.fna\ttests/data/abisko4/73.20110800_S2M.16.fna\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\ttests/data/abisko4/73.20120800_S1D.21.fna\n\
            tests/data/set1/500kb.fna\ttests/data/set1/500kb.fna\n\
            tests/data/set1/500kb.fna\ttests/data/set1/1mbp.fna\n")
            .unwrap();

        // Verify analyse outputs were created properly
        assert!(output_mimag.exists());
        let content = fs::read_to_string(&output_mimag).unwrap();
        let expected = "\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/set1/1mbp.fna\t80.00\t4.00\t0\t1\t0\t10\tMedium quality\n\
            tests/data/set1/500kb.fna\t85.00\t3.00\t1\t1\t1\t15\tMedium quality\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\t90.00\t5.00\t1\t1\t1\t20\tMedium quality\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\t95.00\t2.00\t1\t1\t1\t20\tHigh quality\n";
        assert_eq!(content, expected);

        assert!(output_quality.exists());
    }

    #[test]
    fn test_process_mock_with_reference_genomes() {
        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_quality = tmpdir.path().join("quality_report.tsv");

        setup_mock_bin(
            tmpdir.path(),
            &[
                (String::from("1mbp"), 80.0, 4.0, 0, 1, 0, 10),
                (String::from("500kb"), 85.0, 3.0, 1, 1, 1, 15),
                (String::from("73.20120800_S1X.13"), 90.0, 5.0, 1, 1, 1, 20),
                (String::from("73.20120600_S2D.19"), 95.0, 2.0, 1, 1, 1, 20),
            ],
        );
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        Assert::main_binary()
            .with_env(&[
                ("PATH", new_path),
                ("CHECKM2DB", String::from("/tmp/mockdb")),
            ])
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/set1/500kb.fna",
                "--reference-genomes",
                "tests/data/set1/1mbp.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "--precluster-method",
                "skani",
                "--cluster-method",
                "skani",
                "--precluster-ani",
                "90",
                "--ani",
                "95",
                "--output-cluster-definition",
                "/dev/stdout",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-quality-report",
                output_quality.to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .is("\
            tests/data/abisko4/73.20120600_S2D.19.fna	tests/data/abisko4/73.20120600_S2D.19.fna\n\
            tests/data/abisko4/73.20120600_S2D.19.fna	tests/data/abisko4/73.20120800_S1X.13.fna\n\
            tests/data/set1/500kb.fna	tests/data/set1/500kb.fna\n\
            tests/data/set1/500kb.fna	tests/data/set1/1mbp.fna\n")
            .unwrap();

        // Verify analyse outputs are just for non-reference genomes
        assert!(output_mimag.exists());
        let content = fs::read_to_string(&output_mimag).unwrap();
        let expected = "\
            genome\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\ttRNAs\tMIMAG_quality\n\
            tests/data/abisko4/73.20120800_S1X.13.fna\t90.00\t5.00\t1\t1\t1\t20\tMedium quality\n\
            tests/data/set1/500kb.fna\t85.00\t3.00\t1\t1\t1\t15\tMedium quality\n";
        assert_eq!(content, expected);

        // Quality report should exist
        assert!(output_quality.exists());
    }
}
