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
                "--checkm-tab-table",
                "tests/data/abisko4/abisko4.csv"])
                .succeeds()
                .stdout()
                .is("\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20110800_S2M.16.fna\n\
                tests/data/abisko4/73.20110800_S2M.16.fna	tests/data/abisko4/73.20120800_S1D.21.fna\n")
                .unwrap();
    }
}