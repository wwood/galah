
Runs both `analyse` and `cluster` in sequence, producing a MIMAG quality summary and an
ANI-based cluster definition in a single command. Supports bacteria, archaea, and eukaryotes
via the same domain-classification workflow as `analyse`.

Determines the MIMAG quality score based on completeness, contamination, rRNA, and tRNA presence.
Completeness and contamination are estimated using CheckM2 by default, unless CheckM1/2 quality reports are provided.

Cluster genomes into ANI-based groups for downstream analysis.

```bash
# Process a mixed set of genomes: isiteuk classifies domains automatically (default)
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg \
galah process \
    --genome-fasta-files genome1.fna genome2.fna genome3.fna \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv

# Process bacterial genomes only (skips isiteuk)
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah process \
    --genome-fasta-files genome1.fna genome2.fna \
    --domain-choice bac \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv

# Process a mixed set using a pre-computed isiteuk classification and EukCC quality report
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah process \
    --genome-fasta-files bacteria.fna archaea.fna euk.fna \
    --isiteuk-output isiteuk_results.tsv \
    --eukcc-quality-report eukcc_quality.tsv \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv

# Save the CheckM2 quality report produced during the run for later use
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah process \
    --genome-fasta-files genome1.fna genome2.fna \
    --domain-choice bac \
    --output-quality-report quality_report.tsv \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv

# Save intermediate outputs to a directory; re-running reuses completed steps
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg \
galah process \
    --genome-fasta-list genomes.txt \
    --working-dir work/ \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv
```

See `galah analyse --full-help` for a full description of domain classification, quality tools,
and output format.
