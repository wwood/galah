---
title: Galah cluster
---
# galah cluster

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

# GENOME INPUT

**-f**, **\--genome-fasta-files** *PATH ..*

  Path(s) to FASTA files of each genome e.g.
    `pathA/genome1.fna pathB/genome2.fa`.

<!-- -->

**-d**, **\--genome-fasta-directory** *PATH*

  Directory containing FASTA files of each genome.

<!-- -->

**-x**, **\--genome-fasta-extension** *EXT*

  File extension of genomes in the directory specified with
    `-d/--genome-fasta-directory`. [default: `fna`]

<!-- -->

**\--genome-fasta-list** *PATH*

  File containing FASTA file paths, one per line.

# FILTERING PARAMETERS

**\--checkm2-quality-report** *PATH*

  CheckM version 2 quality_report.tsv (i.e. the `quality_report.tsv`
    in the output directory output of `checkm2 predict ..`) for defining
    genome quality, which is used both for filtering and to rank genomes
    during clustering.

<!-- -->

**\--checkm-tab-table** *PATH*

  CheckM tab table (i.e. the output of
    `checkm .. --tab_table -f PATH ..`). The information contained is
    used like `--checkm2-quality-report`.

<!-- -->

**\--genome-info** *PATH*

  dRep style genome info table for defining quality. The information
    contained is used like `--checkm2-quality-report`.

<!-- -->

**\--min-completeness** *FLOAT*

  Ignore genomes with less completeness than this percentage.
    [default: not set]

<!-- -->

**\--max-contamination** *FLOAT*

  Ignore genomes with more contamination than this percentage.
    [default: not set]

<!-- -->

**\--run-checkm2**

  Run CheckM2 to generate quality scoring used for clustering.
    Requires \--checkm2-db-path or CHECKM2DB env variable to be set.

<!-- -->

**\--checkm2-db-path** *DB_PATH*

  Path to CheckM2 database (required for running CheckM2) [default:
    from CHECKM2DB environment variable]

# CLUSTERING PARAMETERS

**\--ani** *FLOAT*

  Overall ANI level to dereplicate at with the primary clusterer.
    [default: `95`]

<!-- -->

**\--min-aligned-fraction** *FLOAT*

  Min aligned fraction of two genomes for clustering. [default:
    `15`]

<!-- -->

**\--small-genomes**

  Use small-genomes settings in skani calculation. Recommended for
    sequences \< 20kb.

<!-- -->

**\--fragment-length** *FLOAT*

  Length of fragment used in FastANI calculation (i.e. `--fragLen`).
    [default: `3000`]

<!-- -->

**\--quality-formula** *FORMULA*

  Scoring function for genome quality [default:
    `Parks2020_reduced`]. One of:

| formula | description |
|:---|:---|
| `Parks2020_reduced` | (default) A quality formula described in Parks et. al. 2020 https://doi.org/10.1038/s41587-020-0501-8 (Supplementary Table 19) but only including those scoring criteria that can be calculated from the sequence without homology searching: `completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000` |
| `completeness-4contamination` | `completeness-4*contamination` |
| `completeness-5contamination` | `completeness-5*contamination` |
| `dRep` | `completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log10(N50)` |

**\--precluster-ani** *FLOAT*

  Require at least this precluster-derived ANI for preclustering and
    to avoid primary clustering on distant lineages within preclusters.
    [default: `90`]

<!-- -->

**\--precluster-method** *NAME*

  method of calculating rough ANI for dereplication. \'`finch`\' for
    finch MinHash, \'`skani`\' for Skani. [default: `skani`]

<!-- -->

**\--cluster-method** *NAME*

  method of calculating ANI. \'`fastani`\' for FastANI, \'`skani`\'
    for Skani. [default: `skani`]

<!-- -->

**\--cluster-contigs**

  Cluster contigs within a fasta file instead of genomes. When used,
    either \--small-contigs or \--large-contigs must be specified.

<!-- -->

**\--small-contigs**

  Use small-genomes settings in skani when clustering contigs.
    Recommended for contigs \< 20kb. Mutually exclusive with
    \--large-contigs.

<!-- -->

**\--large-contigs**

  Do not use small-genomes settings in skani when clustering contigs.
    Recommended for contigs \>= 20kb. Mutually exclusive with
    \--small-contigs.

<!-- -->

**\--reference-genomes** *PATH \...*

  Reference genomes to cluster against. These should be pre-clustered
    at the chosen %ANI. If quality is provided for representative
    selection, values for these genomes must also be provided. Genomes
    within the precluster ANI cutoff of each reference will be placed in
    the same precluster. Mutually exclusive with
    \--reference-genomes-list.

<!-- -->

**\--reference-genomes-list** *PATH*

  File containing paths to reference genomes (one per line). These
    should be pre-clustered at the chosen %ANI. If quality is provided
    for representative selection, values for these genomes must also be
    provided. Genomes within the precluster ANI cutoff of each reference
    will be placed in the same precluster. Mutually exclusive with
    \--reference-genomes.

# OUTPUT

**\--output-cluster-definition** *PATH*

  Output a file of representative\<TAB\>member lines.

<!-- -->

**\--output-representative-fasta-directory** *PATH*

  Symlink representative genomes into this directory.

<!-- -->

**\--output-representative-fasta-directory-copy** *PATH*

  Copy representative genomes into this directory.

<!-- -->

**\--output-representative-list** *PATH*

  Print newline separated list of paths to representatives into this
    file.

# GENERAL PARAMETERS

**-t**, **\--threads** *INT*

  Number of threads. [default: `1`]

<!-- -->

**-v**, **\--verbose**

  Print extra debugging information

<!-- -->

**-q**, **\--quiet**

  Unless there is an error, do not print log messages

<!-- -->

**-h**, **\--help**

  Output a short usage message.

<!-- -->

**\--full-help**

  Output a full help message and display in \'man\'.

<!-- -->

**\--full-help-roff**

  Output a full help message in raw ROFF format for conversion to
    other formats.

# EXIT STATUS

**0**

  Successful program execution.

<!-- -->

**1**

  Unsuccessful program execution.

<!-- -->

**101**

  The program panicked.

# AUTHOR

> Ben J. Woodcroft, Centre for Microbiome Research, Queensland University of Technology <benjwoodcroft near gmail.com>
