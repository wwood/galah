extern crate assert_cli;

#[cfg(test)]
mod tests {
    use assert_cli::Assert;

    #[test]
    fn test_completeness_4contamination_quality_score(){
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
    fn test_parks2020_reduced_quality_score(){
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
    fn test_output_symlink_directory_dir_exists(){
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
                tdp.to_str().unwrap()])
                .succeeds()
                .stdout()
                .is("")
                .unwrap();
        let out = tdp.join("500kb.fna");
        assert!(out.exists());
        assert!(std::fs::symlink_metadata(out).unwrap().file_type().is_symlink());
        assert!(!tdp.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_symlink_directory_dir_doesnt_exist(){
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        let extra = &format!("{}_and",tdp.to_str().unwrap());
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
                extra])
                .succeeds()
                .stdout()
                .is("")
                .unwrap();
        let extra_path = std::path::Path::new(extra);
        let out = extra_path.join("500kb.fna");
        assert!(extra_path.exists());
        assert!(std::fs::symlink_metadata(out).unwrap().file_type().is_symlink());
        assert!(!extra_path.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_symlink_directory_names_clash(){
        let td = tempfile::TempDir::new().unwrap();
        let tdp = td.path();
        let extra = &format!("{}_and",tdp.to_str().unwrap());
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
                extra])
                .succeeds()
                .stdout()
                .is("")
                .stderr()
                .contains("One or more sequence files have the same file name")
                .unwrap();
        let extra_path = std::path::Path::new(extra);
        let out = extra_path.join("500kb.fna");
        assert!(extra_path.exists());
        assert!(std::fs::symlink_metadata(out).unwrap().file_type().is_symlink());
        assert!(extra_path.join("500kb.fna.1.fna").exists());
        assert!(!extra_path.join("1mbp.fna").exists());
    }

    #[test]
    fn test_output_representative_list(){
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
                "/dev/stdout"])
                .succeeds()
                .stdout()
                .is("
                    tests/data/set1/500kb.fna\n\
                    tests/data/set1_name_clash/500kb.fna\n")
                .unwrap();
    }
}