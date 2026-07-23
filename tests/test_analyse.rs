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
    --genomes) shift; while [[ $# -gt 0 && "${1:0:1}" != "-" ]]; do shift; done;;
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

    /// Write mock binaries (checkm2, barrnap, tRNAscan-SE, isiteuk) to `dir`.
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
        // CheckM2 mock
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

        // Barrnap mock: output rRNA GFF lines to stdout
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

        // tRNAscan-SE mock
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
                "--domain-choice",
                "bac",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t6.35\t0.67\t0\t0\t0\t0\t0\t0\t0\tLow quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t4.08\t0.02\t0\t0\t0\t0\t0\t0\t0\tLow quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t82.17\t0.00\t1\t1\t1\t0\t0\t0\t19\tMedium quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t84.95\t0.03\t1\t1\t1\t0\t0\t0\t18\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_mock_fake() {
        // Extra non-standard tRNAs should not be counted
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t0\t0\t0\t0\t0\t15\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t0\t1\t0\t0\t0\t20\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t16\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t89.90\t2.00\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t5.10\t1\t1\t1\t0\t0\t0\t20\tMedium quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t49.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tLow quality\t\n")
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
            .with_env(&[(
                "PATH",
                new_path
            ),(
                "CHECKM2DB",
                String::from("/tmp/mockdb")
            )])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t11.00\t1\t1\t1\t0\t0\t0\t20\tLow quality\t\n")
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
                "--domain-choice",
                "bac",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t95.50\t1.20\t1\t1\t1\t0\t0\t0\t19\tHigh quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t68.37\t2.91\t0\t1\t1\t0\t0\t0\t10\tMedium quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.50\t1.20\t0\t0\t1\t0\t0\t0\t1\tMedium quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t95.37\t2.91\t0\t0\t0\t0\t0\t0\t0\tMedium quality\t\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_ambiguous_domain_prefers_eukcc_when_higher_completeness() {
        let tmpdir = tempdir().unwrap();
        let genome = "tests/data/set1/1mbp.fna";

        // No row for `genome` at all - it fails every domain cutoff, so isiteuk classification
        // is empty for it and it falls back to being a candidate for every domain.
        fs::write(
            tmpdir.path().join("isiteuk_output.tsv"),
            "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n",
        )
        .unwrap();

        fs::write(
            tmpdir.path().join("checkm2_quality_report.tsv"),
            "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t\
             Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t\
             Total_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes\n\
             1mbp\t40.0\t1.0\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3\t355151\t\
             0.33\t446\t75\t24150\tNone\n",
        )
        .unwrap();

        fs::write(
            tmpdir.path().join("eukcc_quality_report.tsv"),
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t90.0\t2.0\tEukaryota\n",
        )
        .unwrap();

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                genome,
                "--isiteuk-output",
                tmpdir.path().join("isiteuk_output.tsv").to_str().unwrap(),
                "--checkm2-quality-report",
                tmpdir
                    .path()
                    .join("checkm2_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--eukcc-quality-report",
                tmpdir
                    .path()
                    .join("eukcc_quality_report.tsv")
                    .to_str()
                    .unwrap(),
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tEukaryota\t90.00\t2.00\t1\t1\t1\t0\t0\t0\t19\tHigh quality\tno confident isiteuk domain call; assessed under Bacteria, Archaea, Eukaryota; domain resolved to Eukaryota via higher completeness (CheckM2 40.00% vs EukCC 90.00%)\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_ambiguous_domain_prefers_checkm2_when_higher_completeness() {
        let tmpdir = tempdir().unwrap();
        let genome = "tests/data/set1/1mbp.fna";

        fs::write(
            tmpdir.path().join("isiteuk_output.tsv"),
            "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n",
        )
        .unwrap();

        fs::write(
            tmpdir.path().join("checkm2_quality_report.tsv"),
            "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t\
             Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t\
             Total_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes\n\
             1mbp\t95.0\t1.0\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3\t355151\t\
             0.33\t446\t75\t24150\tNone\n",
        )
        .unwrap();

        fs::write(
            tmpdir.path().join("eukcc_quality_report.tsv"),
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t30.0\t2.0\tEukaryota\n",
        )
        .unwrap();

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                genome,
                "--isiteuk-output",
                tmpdir.path().join("isiteuk_output.tsv").to_str().unwrap(),
                "--checkm2-quality-report",
                tmpdir
                    .path()
                    .join("checkm2_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--eukcc-quality-report",
                tmpdir
                    .path()
                    .join("eukcc_quality_report.tsv")
                    .to_str()
                    .unwrap(),
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria,Archaea\t95.00\t1.00\t1\t1\t1\t0\t0\t0\t19\tHigh quality\tno confident isiteuk domain call; assessed under Bacteria, Archaea, Eukaryota; domain resolved to Bacteria,Archaea via higher completeness (CheckM2 95.00% vs EukCC 30.00%)\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_ambiguous_domain_restricts_rrna_trna_to_resolved_domain() {
        let tmpdir = tempdir().unwrap();
        let bin_dir = tempdir().unwrap();
        let genome = "tests/data/set1/1mbp.fna";

        fs::write(
            tmpdir.path().join("isiteuk_output.tsv"),
            "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n",
        )
        .unwrap();
        fs::write(
            tmpdir.path().join("checkm2_quality_report.tsv"),
            "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t\
             Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t\
             Total_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes\n\
             1mbp\t40.0\t1.0\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3\t355151\t\
             0.33\t446\t75\t24150\tNone\n",
        )
        .unwrap();
        fs::write(
            tmpdir.path().join("eukcc_quality_report.tsv"),
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t90.0\t2.0\tEukaryota\n",
        )
        .unwrap();

        // Barrnap mock: --kingdom fun succeeds with one 18S hit; bac/arc are hard failures.
        let barrnap_script = "#!/bin/bash\n\
             kingdom=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20--kingdom) kingdom=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             case \"$kingdom\" in\n\
             \x20\x20fun) echo -e '##gff-version 3\\nmock\\tbarrnap\\trRNA\\t1\\t100\\t.\\t+\\t.\\tName=18S_rRNA;product=18S ribosomal RNA' ;;\n\
             \x20\x20*) echo \"Unexpected barrnap --kingdom $kingdom\" >&2; exit 1 ;;\n\
             esac\n";
        let barrnap = bin_dir.path().join("barrnap");
        fs::write(&barrnap, barrnap_script).unwrap();

        // tRNAscan-SE mock: -E succeeds with 18 standard tRNAs; -B/-A are hard failures.
        let mut trnascan_script = String::from(
            "#!/bin/bash\n\
             mode=\"\"\n\
             out=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20-B|-A|-E) mode=$1; shift;;\n\
             \x20\x20\x20\x20-o) out=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             if [[ \"$mode\" != \"-E\" ]]; then\n\
             \x20\x20echo \"Unexpected tRNAscan-SE mode $mode\" >&2\n\
             \x20\x20exit 1\n\
             fi\n\
             echo -e 'Sequence                      \\t\\ttRNA \\tBounds\\ttRNA\\tAnti\\tIntron Bounds\\tInf\\t      ' > \"$out\"\n\
             echo -e 'Name                          \\ttRNA #\\tBegin\\tEnd  \\tType\\tCodon\\tBegin\\tEnd\\tScore\\tNote' >> \"$out\"\n\
             echo -e '--------                      \\t------\\t-----\\t------\\t----\\t-----\\t-----\\t----\\t------\\t------' >> \"$out\"\n",
        );
        for trna in [
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
            "Met", "Phe", "Pro", "Ser", "Thr", "Trp",
        ] {
            trnascan_script.push_str(&format!(
                "echo -e 'mock\\t1\\t101\\t200\\t{trna}\\tGCC\\t0\\t0\\t20.0\\tNote' >> \"$out\"\n"
            ));
        }
        let trnascan = bin_dir.path().join("tRNAscan-SE");
        fs::write(&trnascan, trnascan_script).unwrap();

        for script in [&barrnap, &trnascan] {
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(script)
                .status();
        }

        let new_path = format!("{}:{}", bin_dir.path().display(), env::var("PATH").unwrap());

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                genome,
                "--isiteuk-output",
                tmpdir.path().join("isiteuk_output.tsv").to_str().unwrap(),
                "--checkm2-quality-report",
                tmpdir
                    .path()
                    .join("checkm2_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--eukcc-quality-report",
                tmpdir
                    .path()
                    .join("eukcc_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tEukaryota\t90.00\t2.00\t0\t0\t0\t1\t0\t0\t18\tHigh quality\tno confident isiteuk domain call; assessed under Bacteria, Archaea, Eukaryota; domain resolved to Eukaryota via higher completeness (CheckM2 40.00% vs EukCC 90.00%)\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_ambiguous_domain_picks_rrna_and_trna_together() {
        let tmpdir = tempdir().unwrap();
        let bin_dir = tempdir().unwrap();
        let genome = "tests/data/set1/1mbp.fna";

        fs::write(
            tmpdir.path().join("isiteuk_output.tsv"),
            "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n",
        )
        .unwrap();
        fs::write(
            tmpdir.path().join("checkm2_quality_report.tsv"),
            "Name\tCompleteness\tContamination\tCompleteness_Model_Used\tTranslation_Table_Used\t\
             Coding_Density\tContig_N50\tAverage_Gene_Length\tGenome_Size\tGC_Content\t\
             Total_Coding_Sequences\tTotal_Contigs\tMax_Contig_Length\tAdditional_Notes\n\
             1mbp\t90.0\t1.0\tGradient Boost (General Model)\t11\t0.885\t5745\t235.3\t355151\t\
             0.33\t446\t75\t24150\tNone\n",
        )
        .unwrap();
        fs::write(
            tmpdir.path().join("eukcc_quality_report.tsv"),
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t20.0\t2.0\tEukaryota\n",
        )
        .unwrap();

        // Barrnap mock: bac kingdom finds all 3 prokaryotic rRNA types; arc kingdom finds none.
        let barrnap_script = "#!/bin/bash\n\
             kingdom=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20--kingdom) kingdom=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             case \"$kingdom\" in\n\
             \x20\x20bac) echo -e '##gff-version 3\\n\
mock\\tbarrnap\\trRNA\\t1\\t100\\t.\\t+\\t.\\tName=5S_rRNA;product=5S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t200\\t300\\t.\\t+\\t.\\tName=16S_rRNA;product=16S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t400\\t500\\t.\\t+\\t.\\tName=23S_rRNA;product=23S ribosomal RNA' ;;\n\
             \x20\x20arc) echo -e '##gff-version 3' ;;\n\
             \x20\x20*) echo \"Unexpected barrnap --kingdom $kingdom\" >&2; exit 1 ;;\n\
             esac\n";
        let barrnap = bin_dir.path().join("barrnap");
        fs::write(&barrnap, barrnap_script).unwrap();

        // tRNAscan-SE mock: -B (Bacteria) finds 5 standard tRNAs; -A (Archaea) finds 20.
        let trnascan_script = "#!/bin/bash\n\
             mode=\"\"\n\
             out=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20-B|-A) mode=$1; shift;;\n\
             \x20\x20\x20\x20-o) out=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             {\n\
             \x20\x20echo -e 'Sequence                      \\t\\ttRNA \\tBounds\\ttRNA\\tAnti\\tIntron Bounds\\tInf\\t      '\n\
             \x20\x20echo -e 'Name                          \\ttRNA #\\tBegin\\tEnd  \\tType\\tCodon\\tBegin\\tEnd\\tScore\\tNote'\n\
             \x20\x20echo -e '--------                      \\t------\\t-----\\t------\\t----\\t-----\\t-----\\t----\\t------\\t------'\n\
             \x20\x20if [[ \"$mode\" == \"-B\" ]]; then\n\
             \x20\x20\x20\x20for t in Ala Arg Asn Asp Cys; do echo -e \"mock\\t1\\t101\\t200\\t$t\\tGCC\\t0\\t0\\t20.0\\tNote\"; done\n\
             \x20\x20elif [[ \"$mode\" == \"-A\" ]]; then\n\
             \x20\x20\x20\x20for t in Ala Arg Asn Asp Cys Gln Glu Gly His Ile Leu Lys Met Phe Pro Ser Thr Trp Tyr Val; do echo -e \"mock\\t1\\t101\\t200\\t$t\\tGCC\\t0\\t0\\t20.0\\tNote\"; done\n\
             \x20\x20else\n\
             \x20\x20\x20\x20echo \"Unexpected tRNAscan-SE mode $mode\" >&2\n\
             \x20\x20fi\n\
             } > \"$out\"\n";
        let trnascan = bin_dir.path().join("tRNAscan-SE");
        fs::write(&trnascan, trnascan_script).unwrap();

        for script in [&barrnap, &trnascan] {
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(script)
                .status();
        }

        let new_path = format!("{}:{}", bin_dir.path().display(), env::var("PATH").unwrap());

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                genome,
                "--isiteuk-output",
                tmpdir.path().join("isiteuk_output.tsv").to_str().unwrap(),
                "--checkm2-quality-report",
                tmpdir
                    .path()
                    .join("checkm2_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--eukcc-quality-report",
                tmpdir
                    .path()
                    .join("eukcc_quality_report.tsv")
                    .to_str()
                    .unwrap(),
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria,Archaea\t90.00\t1.00\t0\t0\t0\t0\t0\t0\t20\tMedium quality\tno confident isiteuk domain call; assessed under Bacteria, Archaea, Eukaryota; domain resolved to Bacteria,Archaea via higher completeness (CheckM2 90.00% vs EukCC 20.00%)\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_domain_choice_completeness() {
        let tmpdir = tempdir().unwrap();
        fs::write(
            tmpdir.path().join("eukcc_quality_report.tsv"),
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t75.0\t5.0\tEukaryota\n",
        )
        .unwrap();

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/set1/1mbp.fna",
                "--domain-choice",
                "completeness",
                "--checkm2-quality-report",
                "tests/data/analyse_file_inputs/checkm2_quality_report.tsv",
                "--eukcc-quality-report",
                tmpdir
                    .path()
                    .join("eukcc_quality_report.tsv")
                    .to_str()
                    .unwrap(),
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria,Archaea\t95.50\t1.20\t1\t1\t1\t0\t0\t0\t19\tHigh quality\t--domain-choice completeness: assessed under Bacteria, Archaea, Eukaryota; domain resolved to Bacteria,Archaea via higher completeness (CheckM2 95.50% vs EukCC 75.00%)\n")
            .unwrap();
    }

    #[test]
    fn test_analyse_domain_choice_all_produces_one_row_per_domain() {
        let bin_dir = tempdir().unwrap();
        let genome = "tests/data/set1/1mbp.fna";

        // Barrnap mock: full prokaryotic set for bac, partial (missing 5S) for arc, full
        // eukaryotic set for fun.
        let barrnap_script = "#!/bin/bash\n\
             kingdom=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20--kingdom) kingdom=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             case \"$kingdom\" in\n\
             \x20\x20bac) echo -e '##gff-version 3\\n\
mock\\tbarrnap\\trRNA\\t1\\t100\\t.\\t+\\t.\\tName=5S_rRNA;product=5S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t200\\t300\\t.\\t+\\t.\\tName=16S_rRNA;product=16S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t400\\t500\\t.\\t+\\t.\\tName=23S_rRNA;product=23S ribosomal RNA' ;;\n\
             \x20\x20arc) echo -e '##gff-version 3\\n\
mock\\tbarrnap\\trRNA\\t200\\t300\\t.\\t+\\t.\\tName=16S_rRNA;product=16S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t400\\t500\\t.\\t+\\t.\\tName=23S_rRNA;product=23S ribosomal RNA' ;;\n\
             \x20\x20fun) echo -e '##gff-version 3\\n\
mock\\tbarrnap\\trRNA\\t1\\t100\\t.\\t+\\t.\\tName=18S_rRNA;product=18S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t200\\t300\\t.\\t+\\t.\\tName=28S_rRNA;product=28S ribosomal RNA\\n\
mock\\tbarrnap\\trRNA\\t400\\t500\\t.\\t+\\t.\\tName=5.8S_rRNA;product=5.8S ribosomal RNA' ;;\n\
             \x20\x20*) echo \"Unexpected barrnap --kingdom $kingdom\" >&2; exit 1 ;;\n\
             esac\n";
        let barrnap = bin_dir.path().join("barrnap");
        fs::write(&barrnap, barrnap_script).unwrap();

        // tRNAscan-SE mock: 20 standard tRNAs regardless of mode (sufficient for all rows).
        let mut trnascan_script = String::from(
            "#!/bin/bash\n\
             out=\"\"\n\
             while [[ $# -gt 0 ]]; do\n\
             \x20\x20case $1 in\n\
             \x20\x20\x20\x20-o) out=$2; shift 2;;\n\
             \x20\x20\x20\x20*) shift;;\n\
             \x20\x20esac\n\
             done\n\
             echo -e 'Sequence                      \\t\\ttRNA \\tBounds\\ttRNA\\tAnti\\tIntron Bounds\\tInf\\t      ' > \"$out\"\n\
             echo -e 'Name                          \\ttRNA #\\tBegin\\tEnd  \\tType\\tCodon\\tBegin\\tEnd\\tScore\\tNote' >> \"$out\"\n\
             echo -e '--------                      \\t------\\t-----\\t------\\t----\\t-----\\t-----\\t----\\t------\\t------' >> \"$out\"\n",
        );
        for trna in [
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
            "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
        ] {
            trnascan_script.push_str(&format!(
                "echo -e 'mock\\t1\\t101\\t200\\t{trna}\\tGCC\\t0\\t0\\t20.0\\tNote' >> \"$out\"\n"
            ));
        }
        let trnascan = bin_dir.path().join("tRNAscan-SE");
        fs::write(&trnascan, trnascan_script).unwrap();

        for script in [&barrnap, &trnascan] {
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(script)
                .status();
        }

        let new_path = format!("{}:{}", bin_dir.path().display(), env::var("PATH").unwrap());

        let eukcc_report = bin_dir.path().join("eukcc_quality_report.tsv");
        fs::write(
            &eukcc_report,
            "fasta\tcompleteness\tcontamination\tncbi_lng\n\
             1mbp\t75.0\t5.0\tEukaryota\n",
        )
        .unwrap();

        Assert::main_binary()
            .with_env(&[("PATH", new_path)])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                genome,
                "--domain-choice",
                "all",
                "--checkm2-quality-report",
                "tests/data/analyse_file_inputs/checkm2_quality_report.tsv",
                "--eukcc-quality-report",
                eukcc_report.to_str().unwrap(),
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t95.50\t1.20\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n\
            tests/data/set1/1mbp.fna\tArchaea\t95.50\t1.20\t0\t1\t1\t0\t0\t0\t20\tMedium quality\t\n\
            tests/data/set1/1mbp.fna\tEukaryota\t75.00\t5.00\t0\t0\t0\t1\t1\t1\t20\tMedium quality\t\n")
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
                "--domain-choice",
                "bac",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/set1/1mbp.fna\tBacteria\t95.50\t1.20\t1\t1\t1\t0\t0\t0\t19\tHigh quality\t\n\
            tests/data/set1/500kb.fna\tBacteria\t58.37\t12.91\t0\t1\t1\t0\t0\t0\t10\tLow quality\t\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t48.37\t1.20\t0\t0\t1\t0\t0\t0\t1\tLow quality\t\n\
            tests/data/abisko4/73.20110800_S2M.16.fna\tBacteria\t38.37\t2.91\t0\t0\t0\t0\t0\t0\t0\tLow quality\t\n")
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
            .with_env(&[("PATH", new_path), ("CHECKM2DB", String::from("/tmp/mockdb"))])
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
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n")
            .unwrap();

        let quality_report_contents = std::fs::read_to_string(&quality_report_path).unwrap();
        assert!(quality_report_contents.contains("Name\tCompleteness\tContamination"));
        assert!(quality_report_contents.contains("73.20120800_S1D.21\t95\t2\t"));
    }

    /// Tests isiteuk domain classification on a mixed set of bacterial, archaeal, and eukaryotic
    /// genomes, running the full quality pipeline (CheckM2, EukCC, Barrnap, tRNAscan-SE).
    /// Requires CHECKM2DB, ISITEUK_METAPACKAGE_PATH, and EUKCC2_DB to be set.
    #[test]
    #[ignore]
    fn test_analyse_real_isiteuk_domain_examples() {
        let _checkm2_db_path = std::env::var("CHECKM2DB")
            .expect("CHECKM2DB environment variable must be set to run this test");
        let _isiteuk_metapackage = std::env::var("ISITEUK_METAPACKAGE_PATH")
            .expect("ISITEUK_METAPACKAGE_PATH environment variable must be set to run this test");
        let _eukcc_db = std::env::var("EUKCC2_DB")
            .expect("EUKCC2_DB environment variable must be set to run this test");

        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/domain_examples/GCF_002008365.1_genomic.fna.gz",
                "tests/data/domain_examples/GCA_003139855.1_genomic.fna.gz",
                "tests/data/domain_examples/binchicken_co8412.34_euk.fna.gz",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
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
        let valid_tiers = ["High quality", "Medium quality", "Low quality"];
        assert!(
            valid_tiers.contains(&euk_line.split('\t').nth(11).unwrap_or("")),
            "Eukaryota genome MIMAG quality should be a valid tier"
        );
    }

    /// Tests EukCC quality assessment directly using --domain-choice euk.
    /// Requires EUKCC2_DB to be set.
    #[test]
    #[ignore]
    fn test_analyse_real_eukcc_domain_example() {
        let _eukcc_db = std::env::var("EUKCC2_DB")
            .expect("EUKCC2_DB environment variable must be set to run this test");

        let tmpdir = tempdir().unwrap();
        let output_mimag = tmpdir.path().join("mimag_summary.tsv");

        Assert::main_binary()
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/domain_examples/binchicken_co8412.34_euk.fna.gz",
                "--domain-choice",
                "euk",
                "--output-mimag-summary",
                output_mimag.to_str().unwrap(),
            ])
            .succeeds()
            .unwrap();

        let content = fs::read_to_string(&output_mimag).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines.len(), 2, "Expected header + 1 genome line");
        let euk_line = lines[1];
        assert_eq!(
            euk_line.split('\t').nth(1),
            Some("Eukaryota"),
            "binchicken_co8412.34_euk should have domain Eukaryota"
        );
        // Verify euk rRNA columns are present (indices 7=18S, 8=28S, 9=5.8S)
        let fields: Vec<&str> = euk_line.split('\t').collect();
        assert_eq!(fields.len(), 13, "Expected 13 columns in output");
        let valid_tiers = ["High quality", "Medium quality", "Low quality"];
        assert!(
            valid_tiers.contains(&fields[11]),
            "MIMAG quality should be a valid tier, got: {}",
            fields[11]
        );
    }

    /// Verify that all working-dir skips fire correctly.
    ///
    /// Pre-populates a working dir with cached outputs for every tool
    /// (isiteuk, CheckM2, EukCC, Barrnap, tRNAscan-SE), then runs galah
    /// with GALAH_*_CMD=false and failing stubs on PATH so that any
    /// unexpected tool invocation causes the test to fail.
    #[test]
    fn test_analyse_working_dir_skips() {
        let work_dir = tempdir().unwrap();
        let bin_dir = tempdir().unwrap();

        let bac_genome = "tests/data/abisko4/73.20120800_S1D.21.fna";
        let bac_stem = "73.20120800_S1D.21";
        let euk_genome = "tests/data/set1/1mbp.fna";
        let euk_stem = "1mbp";

        // Failing stubs: barrnap and tRNAscan-SE must NOT be invoked.
        for tool in ["barrnap", "tRNAscan-SE"] {
            let p = bin_dir.path().join(tool);
            fs::write(
                &p,
                "#!/bin/bash\necho \"Unexpected invocation: $0\" >&2\nexit 1\n",
            )
            .unwrap();
            let _ = std::process::Command::new("chmod")
                .arg("+x")
                .arg(&p)
                .status();
        }

        // isiteuk_output.tsv: classify bac→Bacteria (score 30), euk→Eukaryota (score 25).
        // Default cutoffs: bac/arc = 10, euk = 20, so both pass.
        fs::write(
            work_dir.path().join("isiteuk_output.tsv"),
            format!(
                "genome\tdomain\tnum_in_target_domain\tnum_not_in_target_domain\n\
                 {bac_genome}\td__Bacteria\t30.0\t5.0\n\
                 {euk_genome}\td__Eukaryota\t25.0\t3.0\n"
            ),
        )
        .unwrap();

        // checkm2/quality_report.tsv: only the prokaryotic genome.
        fs::create_dir(work_dir.path().join("checkm2")).unwrap();
        fs::write(
            work_dir.path().join("checkm2/quality_report.tsv"),
            format!(
                "Name\tCompleteness\tContamination\tCompleteness_Model_Used\t\
                 Translation_Table_Used\tCoding_Density\tContig_N50\tAverage_Gene_Length\t\
                 Genome_Size\tGC_Content\tTotal_Coding_Sequences\tTotal_Contigs\t\
                 Max_Contig_Length\tAdditional_Notes\n\
                 {bac_stem}\t95.0\t2.0\tGradient Boost (General Model)\t\
                 11\t0.885\t5745\t235.3\t355151\t0.33\t446\t75\t24150\tNone\n"
            ),
        )
        .unwrap();

        // EukCC output: eukcc_<stem>/eukcc.tsv with completeness/contamination columns.
        let eukcc_dir = work_dir.path().join(format!("eukcc_{euk_stem}"));
        fs::create_dir(&eukcc_dir).unwrap();
        fs::write(
            eukcc_dir.join("eukcc.tsv"),
            "completeness\tcontamination\n90.0\t2.0\n",
        )
        .unwrap();

        // Barrnap GFF for bac genome (kingdom "bac"): 5S, 16S, 23S.
        fs::write(
            work_dir.path().join(format!("{bac_stem}.bac.gff")),
            "##gff-version 3\n\
             mock\tbarrnap\trRNA\t1\t100\t.\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n\
             mock\tbarrnap\trRNA\t200\t300\t.\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n\
             mock\tbarrnap\trRNA\t400\t500\t.\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA\n",
        )
        .unwrap();

        // Barrnap GFF for euk genome (kingdom "fun"): 18S.
        fs::write(
            work_dir.path().join(format!("{euk_stem}.fun.gff")),
            "##gff-version 3\n\
             mock\tbarrnap\trRNA\t1\t1800\t.\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA\n",
        )
        .unwrap();

        // tRNAscan-SE output for bac genome (mode "B"): 20 standard tRNAs → High quality.
        let trna_types = [
            "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "Leu", "Lys",
            "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
        ];
        let trna_header = "Sequence                      \t\ttRNA \tBounds\ttRNA\tAnti\t\
                           Intron Bounds\tInf\t      \n\
                           Name                          \ttRNA #\tBegin\tEnd  \tType\tCodon\t\
                           Begin\tEnd\tScore\tNote\n\
                           --------                      \t------\t-----\t------\t----\t-----\t\
                           -----\t----\t------\t------\n";
        let mut bac_trna = trna_header.to_string();
        for t in &trna_types {
            bac_trna.push_str(&format!("mock\t1\t101\t200\t{t}\tGCC\t0\t0\t20.0\tNote\n"));
        }
        fs::write(
            work_dir.path().join(format!("{bac_stem}.B.trna.out")),
            bac_trna,
        )
        .unwrap();

        // tRNAscan-SE output for euk genome (mode "E"): 18 standard tRNAs → High quality.
        let mut euk_trna = trna_header.to_string();
        for t in trna_types.iter().take(18) {
            euk_trna.push_str(&format!("mock\t1\t101\t200\t{t}\tGCC\t0\t0\t20.0\tNote\n"));
        }
        fs::write(
            work_dir.path().join(format!("{euk_stem}.E.trna.out")),
            euk_trna,
        )
        .unwrap();

        let new_path = format!("{}:{}", bin_dir.path().display(), env::var("PATH").unwrap());

        Assert::main_binary()
            .with_env(&[
                ("PATH", new_path),
                ("GALAH_CHECKM2_CMD", String::from("false")),
                ("GALAH_ISITEUK_CMD", String::from("false")),
                ("GALAH_EUKCC_CMD", String::from("false")),
            ])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                bac_genome,
                euk_genome,
                "--working-dir",
                work_dir.path().to_str().unwrap(),
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is(
                "genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\t\
                 rRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
                 tests/data/abisko4/73.20120800_S1D.21.fna\t\
                 Bacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n\
                 tests/data/set1/1mbp.fna\t\
                 Eukaryota\t90.00\t2.00\t0\t0\t0\t1\t0\t0\t18\tHigh quality\t\n",
            )
            .unwrap();
    }

    #[test]
    fn test_analyse_domain_choice_bac() {
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
            .with_env(&[("PATH", new_path), ("CHECKM2DB", String::from("/tmp/mockdb"))])
            .with_args(&[
                "analyse",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "--domain-choice",
                "bac",
                "--output-mimag-summary",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
            genome\tdomain\tcompleteness\tcontamination\trRNA_5S\trRNA_16S\trRNA_23S\trRNA_18S\trRNA_28S\trRNA_5.8S\ttRNAs\tMIMAG_quality\tnotes\n\
            tests/data/abisko4/73.20120800_S1D.21.fna\tBacteria\t95.00\t2.00\t1\t1\t1\t0\t0\t0\t20\tHigh quality\t\n")
            .unwrap();
    }
}
