---
title: Galah process
---
# galah process

Runs both `analyse` and `cluster` in sequence: determines the MIMAG quality score for each genome,
then dereplicates genomes into ANI-based clusters, using those quality scores to pick each cluster's
representative. Supports bacteria, archaea, and eukaryotes in the same run.

By default, galah runs [isiteuk](https://github.com/wwood/isiteuk) to classify each genome into its
biological domain (Bacteria, Archaea, or Eukaryota), then selects the appropriate quality tool and
marker gene search mode for that domain, exactly as `analyse` does:

* **Bacteria / Archaea**: completeness and contamination estimated by CheckM2; rRNA genes found by
  Barrnap in bacterial/archaeal mode; tRNAs found by tRNAscan-SE in bacterial/archaeal mode.
* **Eukaryota**: completeness and contamination estimated by EukCC; rRNA genes found by Barrnap in
  eukaryotic mode (finds 18S, 28S, 5.8S, 5S); tRNAs found by tRNAscan-SE in eukaryotic mode.

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

### Domain choice

The `--domain-choice` flag controls how galah determines the biological domain for each genome:

| Value | Description |
|---|---|
| `isiteuk` (default) | Run isiteuk to classify each genome automatically |
| `bac` | Treat all genomes as Bacteria |
| `arc` | Treat all genomes as Archaea |
| `euk` | Treat all genomes as Eukaryota |
| `completeness` | Run CheckM2 and EukCC on every genome and keep whichever reports the higher completeness as a single row |
| `all` | Run CheckM2 and EukCC on every genome and report all three domains as separate rows |

`--domain-choice all`, which reports every domain as a separate row, is **not** supported by `process`:
since it produces multiple rows per genome, there is no single quality value left to cluster on. Use
`galah analyse --domain-choice all` instead if you need that view.

### MIMAG quality criteria

**Bacteria / Archaea:**

| Tier | Criteria |
|---|---|
| High quality | Completeness ≥ 90%, contamination < 5%, ≥ 1 each of 5S/16S/23S rRNA, ≥ 18 tRNA types |
| Medium quality | Completeness ≥ 50%, contamination < 10% (and not High quality) |
| Low quality | Completeness < 50% or contamination ≥ 10% |

**Eukaryota:**

| Tier | Criteria |
|---|---|
| High quality | Completeness ≥ 90%, contamination < 5%, ≥ 1 eukaryotic rRNA (18S/28S/5.8S/5S), ≥ 18 tRNA types |
| Medium quality | Completeness ≥ 50%, contamination < 10% (and not High quality) |
| Low quality | Completeness < 50% or contamination ≥ 10% |

See `galah analyse --full-help` for a full description of the `--output-mimag-summary` output format,
and `galah cluster --full-help` for clustering-specific parameters.

# GENOME INPUT

**-f**, **\--genome-fasta-files** *PATH ..*

Path(s) to FASTA files of each genome e.g.
  `pathA/genome1.fna pathB/genome2.fa`.

**-d**, **\--genome-fasta-directory** *PATH*

Directory containing FASTA files of each genome.

**-x**, **\--genome-fasta-extension** *EXT*

File extension of genomes in the directory specified with
  `-d/--genome-fasta-directory`. [default: `fna`]

**\--genome-fasta-list** *PATH*

File containing FASTA file paths, one per line.

# DOMAIN PARAMETERS

**\--domain-choice** *CHOICE*

Method for determining genome domain. \'`isiteuk`\' runs isiteuk to
  classify each genome (default). \'`bac`\', \'`arc`\', \'`euk`\' fix
  all genomes to that domain. \'`completeness`\' runs CheckM2 and EukCC
  (and matching rRNA/tRNA) for every genome and reports only the
  higher-completeness domain\'s result, one row per genome. \'`all`\'
  runs CheckM2 and EukCC for every genome and reports all three domains
  as separate rows. [default: `isiteuk`]

**\--isiteuk-output** *FILE*

Pre-computed isiteuk output TSV. If given, skips running isiteuk when
  \--domain-choice isiteuk.

**\--isiteuk-metapackage** *PATH*

Path to isiteuk metapackage. If not given, uses
  ISITEUK_METAPACKAGE_PATH environment variable.

**\--isiteuk-bacteria-cutoff** *FLOAT*

Minimum isiteuk num_in_target_domain for Bacteria domain assignment.
  [default: `10`]

**\--isiteuk-archaea-cutoff** *FLOAT*

Minimum isiteuk num_in_target_domain for Archaea domain assignment.
  [default: `10`]

**\--isiteuk-eukaryota-cutoff** *FLOAT*

Minimum isiteuk num_in_target_domain for Eukaryota domain assignment.
  [default: `20`]

**\--eukcc-db-path** *PATH*

Path to EukCC database (required for eukaryotic quality assessment).
  If not given, uses EUKCC2_DB environment variable.

**\--eukcc-quality-report** *FILE*

Pre-computed merged EukCC TSV with fasta/completeness/contamination
  columns. If given, skips running EukCC for eukaryotic genomes.

# QUALITY PARAMETERS

**\--quality-method** *NAME*

method for finding genome quality. \'`checkm2`\' for CheckM2.
  [default: `checkm2`]

**\--checkm2-db-path** *PATH*

Path to CheckM2 database (required for CheckM2 quality method). If not
  given, will use CHECKM2DB environment variable if set.

**\--checkm2-quality-report** *PATH*

Path to pre-generated CheckM2 quality_report.tsv file. If given, will
  use this file instead of running quality method.

**\--checkm-tab-table** *PATH*

Path to pre-generated CheckM tab table file. If given, will use this
  file instead of running quality method.

# RNA PARAMETERS

**\--rrna-method** *NAME*

method for finding rRNA genes. \'`barrnap`\' for Barrnap. [default:
  `barrnap`]

**\--trna-method** *NAME*

method for finding tRNA genes. \'`trnascan`\' for tRNAscan-SE.
  [default: `trnascan`]

**\--barrnap-gff-list** *PATH*

Path to two-column TSV file mapping genome paths (as given in input)
  to Barrnap GFF paths (no headers). If given, will use these files
  instead of running rRNA method.

**\--trnascan-out-list** *PATH*

Path to two-column TSV file mapping genome paths (as given in input)
  to tRNAscan-SE output paths (no headers). If given, will use these
  files instead of running tRNA method.

# FILTERING PARAMETERS

**\--checkm2-quality-report** *PATH*

CheckM version 2 quality_report.tsv (i.e. the `quality_report.tsv` in
  the output directory output of `checkm2 predict ..`) for defining
  genome quality, which is used both for filtering and to rank genomes
  during clustering.

**\--checkm-tab-table** *PATH*

CheckM tab table (i.e. the output of
  `checkm .. --tab_table -f PATH ..`). The information contained is used
  like `--checkm2-quality-report`.

**\--genome-info** *PATH*

dRep style genome info table for defining quality. The information
  contained is used like `--checkm2-quality-report`.

**\--min-completeness** *FLOAT*

Ignore genomes with less completeness than this percentage. [default:
  not set]

**\--max-contamination** *FLOAT*

Ignore genomes with more contamination than this percentage.
  [default: not set]

**\--run-checkm2**

Run CheckM2 to generate quality scoring used for clustering. Requires
  \--checkm2-db-path or CHECKM2DB env variable to be set.

**\--checkm2-db-path** *DB_PATH*

Path to CheckM2 database (required for running CheckM2) [default:
  from CHECKM2DB environment variable]

# CLUSTERING PARAMETERS

**\--ani** *FLOAT*

Overall ANI level to dereplicate at with the primary clusterer.
  [default: `95`]

**\--min-aligned-fraction** *FLOAT*

Min aligned fraction of two genomes for clustering. [default: `15`]

**\--small-genomes**

Use small-genomes settings in skani calculation. Recommended for
  sequences \< 20kb.

**\--fragment-length** *FLOAT*

Length of fragment used in FastANI calculation (i.e. `--fragLen`).
  [default: `3000`]

**\--quality-formula** *FORMULA*

Scoring function for genome quality [default: `Parks2020_reduced`].
  One of:

  | formula | description |
  |:---|:---|
  | `Parks2020_reduced` | (default) A quality formula described in Parks et. al. 2020 https://doi.org/10.1038/s41587-020-0501-8 (Supplementary Table 19) but only including those scoring criteria that can be calculated from the sequence without homology searching: `completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000` |
  | `completeness-4contamination` | `completeness-4*contamination` |
  | `completeness-5contamination` | `completeness-5*contamination` |
  | `dRep` | `completeness-5*contamination+contamination*(strain_heterogeneity/100)+0.5*log10(N50)` |

**\--precluster-ani** *FLOAT*

Require at least this precluster-derived ANI for preclustering and to
  avoid primary clustering on distant lineages within preclusters.
  [default: `90`]

**\--precluster-method** *NAME*

method of calculating rough ANI for dereplication. \'`finch`\' for
  finch MinHash, \'`skani`\' for Skani. [default: `skani`]

**\--cluster-method** *NAME*

method of calculating ANI. \'`fastani`\' for FastANI, \'`skani`\' for
  Skani. [default: `skani`]

**\--cluster-contigs**

Cluster contigs within a fasta file instead of genomes. When used,
  either \--small-contigs or \--large-contigs must be specified.

**\--small-contigs**

Use small-genomes settings in skani when clustering contigs.
  Recommended for contigs \< 20kb. Mutually exclusive with
  \--large-contigs.

**\--large-contigs**

Do not use small-genomes settings in skani when clustering contigs.
  Recommended for contigs \>= 20kb. Mutually exclusive with
  \--small-contigs.

**\--low-memory**

Reduce memory use by sketching to file and searching it instead.

**\--skip-sanitise-headers**

Skip checking/rewriting FASTA headers that contain tab characters
  before running skani, passing genome paths straight through unchanged.
  Mainly useful for benchmarking against tools which do not perform this
  sanitizing step. If any input genome actually has a tab character in a
  header line, skani\'s TSV output will be silently corrupted, so only
  use this on genome sets already known not to have tabs in their
  headers.

**\--reference-genomes** *PATH \...*

Reference genomes to cluster against. These should be pre-clustered at
  the chosen %ANI - reference-vs-reference comparisons are never made,
  so this is not checked. Input genomes (\--genome-fasta-files etc.) are
  dereplicated amongst themselves first, as if no reference genomes were
  given at all, and only the resulting representative(s) are then
  compared against these reference genomes - so input genomes do not
  need to be pre-dereplicated beforehand. If quality is provided for
  representative selection, values for these genomes must also be
  provided. Genomes within the precluster ANI cutoff of each reference
  will be placed in the same precluster. Mutually exclusive with
  \--reference-genomes-list.

**\--reference-genomes-list** *PATH*

File containing paths to reference genomes (one per line). These
  should be pre-clustered at the chosen %ANI - reference-vs-reference
  comparisons are never made, so this is not checked. Input genomes
  (\--genome-fasta-files etc.) are dereplicated amongst themselves
  first, as if no reference genomes were given at all, and only the
  resulting representative(s) are then compared against these reference
  genomes - so input genomes do not need to be pre-dereplicated
  beforehand. If quality is provided for representative selection,
  values for these genomes must also be provided. Genomes within the
  precluster ANI cutoff of each reference will be placed in the same
  precluster. Mutually exclusive with \--reference-genomes.

**\--skip-input-dereplication**

When used with \--reference-genomes/\--reference-genomes-list, skip
  dereplicating input genomes amongst themselves first: compare every
  input genome directly against the reference set instead, so
  near-duplicate input genomes are not merged before matching (restores
  the behaviour prior to this flag\'s introduction). Has no effect
  unless reference genomes are given.

# OUTPUT

**\--output-mimag-summary** *PATH*

Output a tsv file summarising the MIMAG status for each genome.

**\--output-quality-report** *PATH*

Output a CheckM2-format quality report TSV file.

**\--output-cluster-definition** *PATH*

Output a file of representative\<TAB\>member lines.

**\--output-representative-fasta-directory** *PATH*

Symlink representative genomes into this directory.

**\--output-representative-fasta-directory-copy** *PATH*

Copy representative genomes into this directory.

**\--output-representative-list** *PATH*

Print newline separated list of paths to representatives into this
  file.

# GENERAL PARAMETERS

**-t**, **\--threads** *INT*

Number of threads. [default: `1`]

**-v**, **\--verbose**

Print extra debugging information

**-q**, **\--quiet**

Unless there is an error, do not print log messages

**-h**, **\--help**

Output a short usage message.

**\--full-help**

Output a full help message and display in \'man\'.

**\--full-help-roff**

Output a full help message in raw ROFF format for conversion to other
  formats.

# EXIT STATUS

**0**

Successful program execution.

**1**

Unsuccessful program execution.

**101**

The program panicked.

# AUTHOR

>     Ben J. Woodcroft, Centre for Microbiome Research, Queensland University of Technology <benjwoodcroft near gmail.com>
