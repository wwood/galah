
Runs both analyse and cluster in sequence.

```bash
# Example: process genomes to produce cluster definition and MIMAG summary
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd galah process \
    --genome-fasta-files genome1.fna genome2.fna \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag.tsv
```
