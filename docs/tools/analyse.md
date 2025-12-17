---
title: Galah analyse
---
# galah analyse

Determines the MIMAG quality score based on completeness, contamination, rRNA, and tRNA presence.
Completeness and contamination are estimated using CheckM2 by default, unless CheckM1/2 quality reports are provided.

on CheckM2 completeness/contamination, For determining MIMAG quality scores for a set of genomes with CheckM2:

```bash
# Example: determine MIMAG quality scores
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd galah analyse --genome-fasta-files genome1.fna genome2.fna --output-mimag-summary mimag.tsv
# Example: determine MIMAG quality scores for a directory of genomes using CheckM2 database specified by argument
galah analyse --genome-fasta-directory input_genomes/ --checkm2-db-path /path/to/checkm2_db.dmnd --output-mimag-summary mimag_summary.tsv
# Example: determine MIMAG quality scores using precomputed CheckM2, Barrnap, and tRNASCAN-SE results
galah analyse --genome-fasta-list genomes.txt --output-mimag-summary mimag_summary.tsv \
    --checkm2-quality-report quality_report.tsv --barrnap-gff-list barrnap_gff_list.tsv --trnascan-out-list trnascan_out_list.tsv
```

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

# QUALITY PARAMETERS

**\--quality-method** *NAME*

  method for finding genome quality. \'`checkm2`\' for CheckM2.
    [default: `checkm2`]

<!-- -->

**\--checkm2-db-path** *PATH*

  Path to CheckM2 database (required for CheckM2 quality method). If
    not given, will use CHECKM2DB environment variable if set.

<!-- -->

**\--checkm2-quality-report** *PATH*

  Path to pre-generated CheckM2 quality_report.tsv file. If given,
    will use this file instead of running quality method.

<!-- -->

**\--checkm-tab-table** *PATH*

  Path to pre-generated CheckM tab table file. If given, will use this
    file instead of running quality method.

# RNA PARAMETERS

**\--rrna-method** *NAME*

  method for finding rRNA genes. \'`barrnap`\' for Barrnap. [default:
    `barrnap`]

<!-- -->

**\--trna-method** *NAME*

  method for finding tRNA genes. \'`trnascan`\' for tRNAscan-SE.
    [default: `trnascan`]

<!-- -->

**\--barrnap-gff-list** *PATH*

  Path to two-column TSV file mapping genome paths (as given in input)
    to Barrnap GFF paths (no headers). If given, will use these files
    instead of running rRNA method.

<!-- -->

**\--trnascan-out-list** *PATH*

  Path to two-column TSV file mapping genome paths (as given in input)
    to tRNAscan-SE output paths (no headers). If given, will use these
    files instead of running tRNA method.

# OUTPUT

**\--output-mimag-summary** *PATH*

  Output a tsv file summarising the MIMAG status for each genome.

<!-- -->

**\--output-quality-report** *PATH*

  Output a CheckM2-format quality report TSV file.

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
