extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;
    use std::env;
    use std::fs;
    use std::path::Path;
    use tempfile::tempdir;

    /// Write a mock isiteuk binary that classifies every genome as Bacteria.
    fn write_mock_isiteuk(dir: &Path) {
        let script = r#"#!/bin/bash
output=""
genome_list=""
while [[ $# -gt 0 ]]; do
  case $1 in
    --output) output=$2; shift 2;;
    --genome-list) genome_list=$2; shift 2;;
    *) shift;;
  esac
done
printf 'genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n' > "$output"
if [[ -n "$genome_list" && -f "$genome_list" ]]; then
  while IFS= read -r genome; do
    [[ -z "$genome" ]] && continue
    printf '%s\td__Bacteria\t14.2\t0.0\n' "$genome" >> "$output"
  done < "$genome_list"
fi
"#;
        let isiteuk = dir.join("isiteuk");
        fs::write(&isiteuk, script).unwrap();
        let _ = std::process::Command::new("chmod")
            .arg("+x")
            .arg(&isiteuk)
            .status();
    }

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

        // isiteuk mock: classify all genomes as Bacteria
        write_mock_isiteuk(dir);

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
                "--domain-choice",
                "bac",
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t6.35\t0.67\t0\t0\t0\t0\t0\t0\t0\tLow quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t4.08\t0.02\t0\t0\t0\t0\t0\t0\t0\tLow quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t82.17\t0.00\t1\t1\t1\t0\t0\t0\t19\tMedium quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t84.95\t0.03\t1\t1\t1\t0\t0\t0\t18\tMedium quality\t\n";
        assert_eq!(content, expected);

        assert!(output_quality.exists());
    }

    /// Tests the full process pipeline (analyse + cluster) with isiteuk domain classification
    /// on a mixed set of bacterial, archaeal, and eukaryotic genomes, running the full quality
    /// pipeline (CheckM2, EukCC, Barrnap, tRNAscan-SE).
    /// Requires CHECKM2DB, ISITEUK_METAPACKAGE_PATH, and EUKCC2_DB to be set.
    #[test]
    #[ignore]
    fn test_process_real_isiteuk_domain_examples() {
        let _checkm2_db_path = std::env::var("CHECKM2DB")
            .expect("CHECKM2DB environment variable must be set to run this test");
        let _isiteuk_metapackage = std::env::var("ISITEUK_METAPACKAGE_PATH")
            .expect("ISITEUK_METAPACKAGE_PATH environment variable must be set to run this test");
        let _eukcc_db = std::env::var("EUKCC2_DB")
            .expect("EUKCC2_DB environment variable must be set to run this test");

        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_clusters = tmpdir.path().join("clusters.tsv");

        Assert::main_binary()
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/domain_examples/GCF_002008365.1_genomic.fna.gz",
                "tests/data/domain_examples/GCA_003139855.1_genomic.fna.gz",
                "tests/data/domain_examples/binchicken_co8412.34_euk.fna.gz",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-cluster-definition",
                output_clusters.to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();

        let content = fs::read_to_string(&output_mimag).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(
            lines[0],
            "genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes",
            "Output header mismatch"
        );
        assert_eq!(lines.len(), 4, "Expected header + 3 genome lines");

        let bac_line = lines
            .iter()
            .find(|l| l.contains("GCF_002008365.1"))
            .expect("Bacteria genome line not found in output");
        assert_eq!(
            bac_line.split('\t').nth(1),
            Some("Bacteria"),
            "GCF_002008365.1 should be classified as Bacteria"
        );

        let arc_line = lines
            .iter()
            .find(|l| l.contains("GCA_003139855.1"))
            .expect("Archaea genome line not found in output");
        assert_eq!(
            arc_line.split('\t').nth(1),
            Some("Archaea"),
            "GCA_003139855.1 should be classified as Archaea"
        );

        let euk_line = lines
            .iter()
            .find(|l| l.contains("binchicken_co8412"))
            .expect("Eukaryota genome line not found in output");
        assert_eq!(
            euk_line.split('\t').nth(1),
            Some("Eukaryota"),
            "binchicken_co8412.34_euk should be classified as Eukaryota"
        );

        // Each genome from a different domain should form its own cluster
        assert!(
            output_clusters.exists(),
            "Cluster definition file should be created"
        );
        let cluster_content = fs::read_to_string(&output_clusters).unwrap();
        assert_eq!(
            cluster_content.lines().count(),
            3,
            "Expected 3 clusters (one per genome)"
        );
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t85.00\t3.00\t1\t1\t1\t0\t0\t0\t15\tMedium quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t80.00\t4.00\t0\t1\t0\t0\t0\t0\t10\tMedium quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t90.00\t5.00\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n";
        assert_eq!(content, expected);

        assert!(output_quality.exists());
    }

    #[test]
    fn test_process_mock_low_memory() {
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
                "--low-memory",
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t85.00\t3.00\t1\t1\t1\t0\t0\t0\t15\tMedium quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t80.00\t4.00\t0\t1\t0\t0\t0\t0\t10\tMedium quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t90.00\t5.00\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n";
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t80.00\t4.00\t0\t1\t0\t0\t0\t0\t10\tMedium quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t85.00\t3.00\t1\t1\t1\t0\t0\t0\t15\tMedium quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t90.00\t5.00\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n";
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1X.13.fna\tBacteria\t90.00\t5.00\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t85.00\t3.00\t1\t1\t1\t0\t0\t0\t15\tMedium quality\t\n";
        assert_eq!(content, expected);

        // Quality report should exist
        assert!(output_quality.exists());
    }

    #[test]
    fn test_process_mock_combined_prok_euk_quality_report() {
        let tmpdir = tempdir().unwrap();

        // Prokaryotic pair: two distinct real genomes already known to cluster together
        // (see test_process_mock_with_reference_genomes) with a large completeness gap so the
        // Parks2020_reduced quality formula's per-genome penalties can't flip the ranking.
        let prok_hi = "tests/data/abisko4/73.20120600_S2D.19.fna";
        let prok_lo = "tests/data/abisko4/73.20120800_S1X.13.fna";

        // Synthetic eukaryotic pair: identical sequence content (guaranteed 100% ANI, so they
        // always cluster together) but distinct filenames/quality, labelled Eukaryota purely
        // via the injected isiteuk output below. The real content doesn't matter here since
        // completeness/contamination and domain both come from injected reports.
        let euk_hi = tmpdir.path().join("euk_hi.fna");
        let euk_lo = tmpdir.path().join("euk_lo.fna");
        fs::copy("tests/data/set1/500kb.fna", &euk_hi).unwrap();
        fs::copy("tests/data/set1/500kb.fna", &euk_lo).unwrap();
        let euk_hi = euk_hi.to_str().unwrap().to_string();
        let euk_lo = euk_lo.to_str().unwrap().to_string();

        // Pre-computed isiteuk domain classification (num_in_target_domain above the default
        // cutoffs of 10 for bac/arc, 20 for euk).
        fs::write(
            tmpdir.path().join("isiteuk_output.tsv"),
            format!(
                "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n\
                 {prok_hi}\td__Bacteria\t30.0\t0.0\n\
                 {prok_lo}\td__Bacteria\t30.0\t0.0\n\
                 {euk_hi}\td__Eukaryota\t25.0\t0.0\n\
                 {euk_lo}\td__Eukaryota\t25.0\t0.0\n"
            ),
        )
        .unwrap();

        // Mocked checkm2 binary standing in for the real tool (prokaryotes only - no entry for
        // eukaryotes, matching what real CheckM2 would produce).
        setup_mock_bin(
            tmpdir.path(),
            &[
                (String::from("73.20120600_S2D.19"), 95.0, 1.0, 0, 0, 0, 0),
                (String::from("73.20120800_S1X.13"), 60.0, 1.0, 0, 0, 0, 0),
            ],
        );
        let path = env::var("PATH").unwrap();
        let new_path = format!("{}:{}", tmpdir.path().display(), path);

        // Pre-computed EukCC quality report (eukaryotes only - no entry for prokaryotes).
        let eukcc_report = tmpdir.path().join("eukcc_quality_report.tsv");
        fs::write(
            &eukcc_report,
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             euk_hi\t90.0\t2.0\tEukaryota\n\
             euk_lo\t40.0\t2.0\tEukaryota\n",
        )
        .unwrap();

        // Minimal (header-only) Barrnap/tRNAscan-SE outputs, reused for all four genomes -
        // rRNA/tRNA content is irrelevant to this test, only completeness/contamination-driven
        // representative selection is being checked.
        let empty_gff = tmpdir.path().join("empty.gff");
        fs::write(&empty_gff, "##gff-version 3\n").unwrap();
        let empty_trna = tmpdir.path().join("empty.trna.out");
        fs::write(
            &empty_trna,
            "Sequence                      \t\ttRNA \tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      \n\
             Name                          \ttRNA #\tBegin\tEnd  \tType\tCodon\tBegin\tEnd\tScore\tNote\n\
             --------                      \t------\t-----\t------\t----\t-----\t-----\t----\t------\t------\n",
        )
        .unwrap();

        // These are matched against the exact --genome-fasta-files paths (no stem fallback),
        // so use the same literal path strings passed to the command below.
        let empty_gff_str = empty_gff.to_str().unwrap();
        let empty_trna_str = empty_trna.to_str().unwrap();
        let barrnap_list = tmpdir.path().join("barrnap_gff_list.tsv");
        fs::write(
            &barrnap_list,
            format!(
                "{prok_hi}\t{empty_gff_str}\n\
                 {prok_lo}\t{empty_gff_str}\n\
                 {euk_hi}\t{empty_gff_str}\n\
                 {euk_lo}\t{empty_gff_str}\n"
            ),
        )
        .unwrap();
        let trnascan_list = tmpdir.path().join("trnascan_out_list.tsv");
        fs::write(
            &trnascan_list,
            format!(
                "{prok_hi}\t{empty_trna_str}\n\
                 {prok_lo}\t{empty_trna_str}\n\
                 {euk_hi}\t{empty_trna_str}\n\
                 {euk_lo}\t{empty_trna_str}\n"
            ),
        )
        .unwrap();

        let output_clusters = tmpdir.path().join("clusters.tsv");
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");
        let output_quality = tmpdir.path().join("quality_report.tsv");

        Assert::main_binary()
            .with_env(&[
                ("PATH", new_path),
                ("CHECKM2DB", String::from("/tmp/mockdb")),
            ])
            .with_args(&[
                "process",
                "--genome-fasta-files",
                // Lower-quality genome listed first in each pair, so a correct result can only
                // come from genuine quality-based ordering, not input-order fallback.
                prok_lo,
                prok_hi,
                &euk_lo,
                &euk_hi,
                "--isiteuk-output",
                tmpdir.path().join("isiteuk_output.tsv").to_str().unwrap(),
                "--eukcc-quality-report",
                eukcc_report.to_str().unwrap(),
                "--barrnap-gff-list",
                barrnap_list.to_str().unwrap(),
                "--trnascan-out-list",
                trnascan_list.to_str().unwrap(),
                "--precluster-method",
                "skani",
                "--cluster-method",
                "skani",
                "--precluster-ani",
                "90",
                "--ani",
                "95",
                "--output-cluster-definition",
                output_clusters.to_str().unwrap(),
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
                "--output-quality-report",
                output_quality.to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();

        // Paths under tmpdir aren't known ahead of time, so parse the cluster definition
        // rather than asserting on an exact hardcoded string.
        let cluster_content = fs::read_to_string(&output_clusters).unwrap();
        let lines: Vec<&str> = cluster_content.lines().collect();
        assert_eq!(
            lines.len(),
            4,
            "Expected 4 rows (2 prok + 2 euk members): {cluster_content}"
        );

        let prok_members: Vec<&str> = lines
            .iter()
            .filter(|l| l.starts_with(&format!("{prok_hi}\t")))
            .copied()
            .collect();
        assert_eq!(
            prok_members.len(),
            2,
            "Higher-quality (CheckM2) prok genome should be representative of both prok genomes: {cluster_content}"
        );
        assert!(prok_members.iter().any(|l| l.ends_with(prok_hi)));
        assert!(prok_members.iter().any(|l| l.ends_with(prok_lo)));

        let euk_members: Vec<&str> = lines
            .iter()
            .filter(|l| l.starts_with(&format!("{euk_hi}\t")))
            .copied()
            .collect();
        assert_eq!(
            euk_members.len(),
            2,
            "Higher-quality (EukCC) euk genome should be representative of both euk genomes: {cluster_content}"
        );
        assert!(euk_members.iter().any(|l| l.ends_with(euk_hi.as_str())));
        assert!(euk_members.iter().any(|l| l.ends_with(euk_lo.as_str())));

        // Sanity-check the merged per-genome quality made it into the MIMAG summary too.
        let mimag_content = fs::read_to_string(&output_mimag).unwrap();
        assert!(mimag_content.contains(&format!("{prok_hi}\tBacteria\t95.00\t1.00")));
        assert!(mimag_content.contains(&format!("{prok_lo}\tBacteria\t60.00\t1.00")));
        assert!(mimag_content.contains(&format!("{euk_hi}\tEukaryota\t90.00\t2.00")));
        assert!(mimag_content.contains(&format!("{euk_lo}\tEukaryota\t40.00\t2.00")));
    }

    /// `--domain-choice all` produces multiple rows per genome (one per domain), which has no
    /// single quality value to cluster on - `process` must reject it rather than guess.
    #[test]
    fn test_process_domain_choice_all_is_rejected() {
        Assert::main_binary()
            .with_args(&[
                "process",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "--domain-choice",
                "all",
                "--output-mimag-summary",
                "/dev/stdout",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .fails()
            .stderr()
            .contains("--domain-choice all is not supported by `process`")
            .unwrap();
    }
}
