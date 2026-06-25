
Determines the MIMAG quality score based on completeness, contamination, rRNA, and tRNA presence.
Completeness and contamination are estimated using CheckM2 by default, unless CheckM1/2 quality reports are provided.

```bash
# Example: determine MIMAG quality scores
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd galah analyse --genome-fasta-files genome1.fna genome2.fna --output-mimag-summary mimag.tsv
# Example: determine MIMAG quality scores for a directory of genomes using CheckM2 database specified by argument
galah analyse --genome-fasta-directory input_genomes/ --checkm2-db-path /path/to/checkm2_db.dmnd --output-mimag-summary mimag_summary.tsv
# Example: determine MIMAG quality scores using precomputed CheckM2, Barrnap, and tRNASCAN-SE results
galah analyse --genome-fasta-list genomes.txt --output-mimag-summary mimag_summary.tsv \
    --checkm2-quality-report quality_report.tsv --barrnap-gff-list barrnap_gff_list.tsv --trnascan-out-list trnascan_out_list.tsv
```
