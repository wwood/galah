
Runs both analyse and cluster in sequence.

Determines the MIMAG quality score based on completeness, contamination, rRNA, and tRNA presence.
Completeness and contamination are estimated using CheckM2 by default, unless CheckM1/2 quality reports are provided.

Cluster genomes into ANI-based groups for downstream analysis.

```bash
# Example: process genomes to produce cluster definition and MIMAG summary
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd galah process \
    --genome-fasta-files genome1.fna genome2.fna \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag.tsv
```
