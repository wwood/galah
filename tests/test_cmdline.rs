extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;

    #[test]
    fn test_completeness_4contamination_quality_score() {
        // 73.20120800_S1D.21	p__Euryarchaeota (UID49)	95	228	153	10	218	0	0	0	0	95.21	0.00	0.00
        // 73.20110800_S2M.16	p__Euryarchaeota (UID49)	95	228	153	8	219	1	0	0	0	95.92	0.65	0.00
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "completeness-4contamination",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20120800_S1D.21.fna	tests/data/abisko4/73.20120800_S1D.21.fna\n\
                tests/data/abisko4/73.20120800_S1D.21.fna	tests/data/abisko4/73.20110800_S2M.16.fna\n")
                .unwrap();
    }

    #[test]
    fn test_parks2020_reduced_quality_score() {
        // 73.20120800_S1D.21	p__Euryarchaeota (UID49)	95	228	153	10	218	0	0	0	0	95.21	0.00	0.00
        // 73.20110800_S2M.16	p__Euryarchaeota (UID49)	95	228	153	8	219	1	0	0	0	95.92	0.65	0.00
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "Parks2020_reduced",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20110800_S2M.16.fna\n\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20120800_S1D.21.fna\n")
                .unwrap();
    }

    #[test]
    fn test_output_symlink_directory_dir_exists() {
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "Parks2020_reduced",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-fasta-directory",
                tdp.to_str().unwrap(),
            ])
            .succeeds()
            .stdout()
            .is("")
            .unwrap();
        let out = tdp.join("500kb.fna");
        assert!(out.exists());
        assert!(std::fs::symlink_metadata(out)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(!tdp.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_symlink_directory_dir_doesnt_exist() {
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        let extra = &format!("{}_and", tdp.to_str().unwrap());
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "Parks2020_reduced",
                "--genome-fasta-files",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-fasta-directory",
                extra,
            ])
            .succeeds()
            .stdout()
            .is("")
            .unwrap();
        let extra_path = std::path::Path::new(extra);
        let out = extra_path.join("500kb.fna");
        assert!(extra_path.exists());
        assert!(std::fs::symlink_metadata(out)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(!extra_path.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_symlink_directory_names_clash() {
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        let extra = &format!("{}_and", tdp.to_str().unwrap());
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "Parks2020_reduced",
                "--genome-fasta-files",
                "tests/data/set1_name_clash/500kb.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-fasta-directory",
                extra,
            ])
            .succeeds()
            .stdout()
            .is("")
            .stderr()
            .contains("One or more sequence files have the same file name")
            .unwrap();
        let extra_path = std::path::Path::new(extra);
        let out = extra_path.join("500kb.fna");
        assert!(extra_path.exists());
        assert!(std::fs::symlink_metadata(out)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(extra_path.join("500kb.fna.1.fna").exists());
        assert!(!extra_path.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_representative_list() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/set1_name_clash/500kb.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-list",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("
                    tests/data/set1/500kb.fna\n\
                    tests/data/set1_name_clash/500kb.fna\n")
            .unwrap();
    }

    #[test]
    fn test_output_symlink_directory_names_clash_copy() {
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        let extra = &format!("{}_and", tdp.to_str().unwrap());
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--quality-formula",
                "Parks2020_reduced",
                "--genome-fasta-files",
                "tests/data/set1_name_clash/500kb.fna",
                "tests/data/set1/500kb.fna",
                "tests/data/set1/1mbp.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-fasta-directory-copy",
                extra,
            ])
            .succeeds()
            .stdout()
            .is("")
            .stderr()
            .contains("One or more sequence files have the same file name")
            .unwrap();
        let extra_path = std::path::Path::new(extra);
        let out = extra_path.join("500kb.fna");
        assert!(extra_path.exists());
        assert!(!std::fs::symlink_metadata(out)
            .unwrap()
            .file_type()
            .is_symlink());
        assert!(extra_path.join("500kb.fna.1.fna").exists());
        assert!(!extra_path.join("1mbp.fna").exists());
    }

    #[test]
    fn test_min_aligned_fraction() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/set2/1mbp.fna",
                "tests/data/set2/1mbp.half_aligned.fna",
                "--min-aligned-fraction",
                "0.2",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-list",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("
                    tests/data/set2/1mbp.fna\n")
            .unwrap();

        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/set2/1mbp.fna",
                "tests/data/set2/1mbp.half_aligned.fna",
                "--min-aligned-fraction",
                "0.6",
                "--precluster-method", // only needed temporarily
                "finch",
                "--output-representative-list",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("
                    tests/data/set2/1mbp.fna\n\
                    tests/data/set2/1mbp.half_aligned.fna\n")
            .unwrap();
    }

    #[test]
    fn test_skani_clusterer() {
        // 73.20120800_S1D.21	p__Euryarchaeota (UID49)	95	228	153	10	218	0	0	0	0	95.21	0.00	0.00
        // 73.20110800_S2M.16	p__Euryarchaeota (UID49)	95	228	153	8	219	1	0	0	0	95.92	0.65	0.00
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--cluster-method",
                "skani",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20110800_S2M.16.fna\n\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20120800_S1D.21.fna\n")
                .unwrap();
    }

    #[test]
    fn test_skani_checkm2() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna",
                "tests/data/abisko4/73.20110800_S2M.16.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--cluster-method",
                "skani",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm2-quality-report",
                "tests/data/abisko4/abisko4_quality_report.tsv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20110800_S2M.16.fna\n\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20120800_S1D.21.fna\n")
                .unwrap();
    }

    #[test]
    fn test_skani_skani_clusterer() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
                "--precluster-method",
                "skani",
                "--cluster-method",
                "skani",
                "--precluster-ani",
                "99",
                "--ani",
                "95",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20120800_S1X.13.fna	tests/data/abisko4/73.20120800_S1X.13.fna\n\
                tests/data/abisko4/73.20120800_S1X.13.fna	tests/data/abisko4/73.20110800_S2D.13.fna\n\
                tests/data/abisko4/73.20120800_S1X.13.fna	tests/data/abisko4/73.20120600_S2D.19.fna\n\
                tests/data/abisko4/73.20120800_S1X.13.fna	tests/data/abisko4/73.20120700_S3X.12.fna\n")
                .unwrap();
    }

    #[test]
    fn test_skani_with_low_ani() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "tests/data/abisko4/73.20120700_S3X.12.fna",
                "tests/data/abisko4/73.20110800_S2D.13.fna",
                "--precluster-method",
                "skani",
                "--cluster-method",
                "skani",
                "--precluster-ani",
                "80",
                "--ani",
                "80",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv",
            ])
            .fails()
            .stderr()
            .contains(
                "Error: skani produces inaccurate results with ANI less than 85%. Provided: 80",
            )
            .unwrap();
    }

    #[test]
    fn test_github7() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/antonio_mags/BE_RX_R2_MAG52.fna",
                "tests/data/antonio_mags/BE_RX_R3_MAG189.fna",
                "--precluster-method", // only needed temporarily
                "finch",
                "--precluster-ani",
                "90",
                "--ani",
                "95",
                "--min-aligned-fraction",
                "60",
                "--output-representative-list",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("tests/data/antonio_mags/BE_RX_R2_MAG52.fna\n")
            .unwrap();
    }

    #[test]
    fn test_genome_cluster_with_small_genomes() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1X.13.fna",
                "tests/data/abisko4/73.20120600_S2D.19.fna",
                "--small-genomes",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .contains("tests/data/abisko4/73.20120800_S1X.13.fna")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs.fna",
                "--cluster-contigs",
                "--large-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024\n\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024_2\n\
                73.20110600_S2D.10_contig_50844	73.20110600_S2D.10_contig_50844\n\
                73.20110600_S2D.10_contig_37820	73.20110600_S2D.10_contig_37820\n")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_specific() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs_specific.fna",
                "--cluster-contigs",
                "--small-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024\n\
                73.20110600_S2D.10_contig_13024	100ANI_100AF\n\
                73.20110600_S2D.10_contig_13024	100ANI_100refAF_90queryAF\n\
                73.20110600_S2D.10_contig_13024	100ANI_90refAF_90queryAF\n\
                73.20110600_S2D.10_contig_13024	100ANI_80refAF_80queryAF\n\
                73.20110600_S2D.10_contig_13024	96ANI_80refAF_80queryAF\n\
                94ANI_80refAF_80queryAF	94ANI_80refAF_80queryAF\n\
                73.20110600_S2D.10_contig_50844	73.20110600_S2D.10_contig_50844\n\
                73.20110600_S2D.10_contig_37820	73.20110600_S2D.10_contig_37820\n")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_without_size_flag_should_error() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs.fna",
                "--cluster-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .fails()
            .stderr()
            .contains("Error: When --cluster-contigs is used, either --small-contigs or --large-contigs must be specified.")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_with_both_flags_should_error() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs.fna",
                "--cluster-contigs",
                "--small-contigs",
                "--large-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .fails()
            .stderr()
            .contains("error: the argument '--small-contigs' cannot be used with '--large-contigs'")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_multiple_files() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs.fna",
                "tests/data/contigs/contigs_extra.fna",
                "--cluster-contigs",
                "--small-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024\n\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024_2\n\
                73.20110600_S2D.10_contig_13024	73.20110600_S2D.10_contig_13024_3\n\
                73.20110600_S2D.10_contig_50844	73.20110600_S2D.10_contig_50844\n\
                73.20110600_S2D.10_contig_37820	73.20110600_S2D.10_contig_37820\n")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_rep_bug_large() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs_rep_bug.fna",
                "--cluster-contigs",
                "--large-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
                k141_313035 flag=1 multi=13.9893 len=27966	k141_313035 flag=1 multi=13.9893 len=27966\n\
                k141_313035 flag=1 multi=13.9893 len=27966	k141_401621 flag=1 multi=12.7497 len=42088\n\
                k141_313035 flag=1 multi=13.9893 len=27966	NODE_1070_length_34582_cov_11.872969\n")
            .unwrap();
    }

    #[test]
    fn test_contig_cluster_rep_bug_small() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/contigs/contigs_rep_bug.fna",
                "--cluster-contigs",
                "--small-contigs",
                "--output-cluster-definition",
                "/dev/stdout",
            ])
            .succeeds()
            .stdout()
            .is("\
                k141_313035 flag=1 multi=13.9893 len=27966	k141_313035 flag=1 multi=13.9893 len=27966\n\
                k141_313035 flag=1 multi=13.9893 len=27966	k141_401621 flag=1 multi=12.7497 len=42088\n\
                NODE_1070_length_34582_cov_11.872969	NODE_1070_length_34582_cov_11.872969\n")
            .unwrap();
    }

    #[test]
    fn test_github53() {
        Assert::main_binary()
            .with_args(&[
                "cluster",
                "--genome-fasta-files",
                "tests/data/abisko4/73.20120800_S1D.21.fna.gz",
                "tests/data/abisko4/73.20110800_S2M.16.fna.gz",
                "--output-cluster-definition",
                "/dev/stdout",
                "--checkm2-quality-report",
                "tests/data/abisko4/abisko4_quality_report.tsv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20110800_S2M.16.fna.gz	tests/data/abisko4/73.20110800_S2M.16.fna.gz\n\
                tests/data/abisko4/73.20110800_S2M.16.fna.gz	tests/data/abisko4/73.20120800_S1D.21.fna.gz\n")
                .unwrap();
    }

    // #[test]
    // fn test_fraglen() {
    //     Assert::main_binary()
    //         .with_args(&[
    //             "cluster",
    //             "--genome-fasta-files",
    //             "tests/data/fraglen_test/sequence2.fna",
    //             "tests/data/fraglen_test/sequence1.fna",
    //             "--min-aligned-fraction",
    //             "0.8",
    //             "--precluster-method", // only needed temporarily
    //             "finch",
    //             "--output-representative-list",
    //             "/dev/stdout",
    //         ])
    //         .succeeds()
    //         .stdout()
    //         .is("
    //             tests/data/fraglen_test/sequence2.fna\n\
    //             tests/data/fraglen_test/sequence1.fna\n")
    //         .unwrap();

    //     Assert::main_binary()
    //         .with_args(&[
    //             "cluster",
    //             "--genome-fasta-files",
    //             "tests/data/fraglen_test/sequence2.fna",
    //             "tests/data/fraglen_test/sequence1.fna",
    //             "--min-aligned-fraction",
    //             "0.8",
    //             "--fragment-length",
    //             "1000",
    //             "--precluster-method", // only needed temporarily
    //             "finch",
    //             "--output-representative-list",
    //             "/dev/stdout",
    //         ])
    //         .succeeds()
    //         .stdout()
    //         .is("
    //             tests/data/fraglen_test/sequence2.fna")
    //         .unwrap();
    // }
}
