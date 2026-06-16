
Cluster genomes into ANI-based groups for downstream analysis.

```bash
# Example: cluster genomes at 95% ANI and produce cluster definition file
galah cluster --genome-fasta-files genome1.fna genome2.fna --output-cluster-definition clusters.tsv
# Example: cluster contigs and produce cluster definition file (can be used for viruses/plasmids)
galah cluster --cluster-contigs --small-genomes --genome-fasta-files contigs.fna --output-cluster-definition clusters.tsv
# Example: cluster a directory of genomes and create a new directory of symlinked FASTA files of representatives
galah cluster --genome-fasta-directory input_genomes/ --output-representative-fasta-directory output_directory/ 
# Example: cluster genomes specified in genomes.txt at 95% ANI after preclustering at 90% using finch method
galah cluster --ani 95 --precluster-ani 90 --precluster-method finch --genome-fasta-list genomes.txt --output-cluster-definition clusters.tsv
# Example: cluster a set of genomes and then their representatives against a set of reference genomes (reduces memory usage against clustering all together)
galah cluster --genome-fasta-directory input_genomes/ --output-representative-list genome_reps.txt
galah cluster --genome-fasta-list genome_reps.txt --reference-genomes-list reference_genomes.txt --output-cluster-definition clusters.tsv
# Example: cluster a large set of genomes using low-memory mode
galah cluster --low-memory --genome-fasta-directory input_genomes/ --output-representative-fasta-directory output_directory/
```

### Precluster ANI

Similar to dRep, Galah operates in two stages. In the first, a fast
pre-clustering distance ([finch](https://github.com/onecodex/finch-rs)
or [skani](https://github.com/bluenote-1577/skani)) is
calculated between each pair of genomes. Genome pairs are only considered as
potentially in the same cluster with [skani](https://github.com/bluenote-1577/skani) or
[FastANI](https://github.com/ParBLiSS/FastANI) if the prethreshold ANI is
greater than the specified value. By default, the precluster ANI is set at 95%
and the final ANI is set at 99%.
