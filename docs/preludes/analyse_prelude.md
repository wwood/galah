Determines the MIMAG quality score for each genome based on completeness, contamination, rRNA gene
presence, and tRNA diversity. Supports bacteria, archaea, and eukaryotes.

By default, galah runs [isiteuk](https://github.com/wwood/isiteuk) to classify each genome into its
biological domain (Bacteria, Archaea, or Eukaryota), then selects the appropriate quality tool and
marker gene search mode for that domain:

* **Bacteria / Archaea**: completeness and contamination estimated by CheckM2; rRNA genes found by
  Barrnap in bacterial/archaeal mode; tRNAs found by tRNAscan-SE in bacterial/archaeal mode.
* **Eukaryota**: completeness and contamination estimated by EukCC; rRNA genes found by Barrnap in
  eukaryotic mode (finds 18S, 28S, 5.8S, 5S); tRNAs found by tRNAscan-SE in eukaryotic mode.

```bash
# Classify genomes by domain automatically with isiteuk (default), then run the appropriate tools
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg \
galah analyse \
    --genome-fasta-files genome1.fna genome2.fna \
    --output-mimag-summary mimag_summary.tsv

# Treat all genomes as bacteria (skips isiteuk; useful when domain is known)
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah analyse \
    --genome-fasta-files genome1.fna genome2.fna \
    --domain-choice bac \
    --output-mimag-summary mimag_summary.tsv

# Treat all genomes as eukaryotes
EUKCC2_DB=/path/to/eukcc2_db_ver_1.1 \
galah analyse \
    --genome-fasta-files euk1.fna euk2.fna \
    --domain-choice euk \
    --output-mimag-summary mimag_summary.tsv

# Mixed set: isiteuk classifies automatically; provide pre-computed EukCC report to skip EukCC run
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg \
galah analyse \
    --genome-fasta-files bacteria.fna archaea.fna euk.fna \
    --eukcc-quality-report eukcc_quality.tsv \
    --output-mimag-summary mimag_summary.tsv

# Use a pre-computed isiteuk classification result (avoids re-running isiteuk)
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah analyse \
    --genome-fasta-files genome1.fna genome2.fna \
    --isiteuk-output isiteuk_results.tsv \
    --output-mimag-summary mimag_summary.tsv

# Use a CheckM2 database specified by argument instead of environment variable
galah analyse \
    --genome-fasta-directory input_genomes/ \
    --domain-choice bac \
    --checkm2-db-path /path/to/checkm2_db.dmnd \
    --output-mimag-summary mimag_summary.tsv

# Use pre-computed CheckM2, Barrnap, and tRNAscan-SE results
galah analyse \
    --genome-fasta-list genomes.txt \
    --domain-choice bac \
    --checkm2-quality-report quality_report.tsv \
    --barrnap-gff-list barrnap_gff_list.tsv \
    --trnascan-out-list trnascan_out_list.tsv \
    --output-mimag-summary mimag_summary.tsv

# Save intermediate outputs to a directory; re-running reuses completed steps
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg \
galah analyse \
    --genome-fasta-list genomes.txt \
    --working-dir work/ \
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

### Output format

The `--output-mimag-summary` TSV contains one row per genome with these columns:

| Column | Description |
|---|---|
| `genome` | Path to the genome FASTA file |
| `domain` | Assigned domain (Bacteria, Archaea, or Eukaryota) |
| `completeness` | Completeness estimate (%) |
| `contamination` | Contamination estimate (%) |
| `rRNA_5S` | Number of 5S rRNA genes found |
| `rRNA_16S` | Number of 16S rRNA genes found |
| `rRNA_23S` | Number of 23S rRNA genes found |
| `rRNA_18S` | Number of 18S rRNA genes found (eukaryotes) |
| `rRNA_28S` | Number of 28S rRNA genes found (eukaryotes) |
| `rRNA_5.8S` | Number of 5.8S rRNA genes found (eukaryotes) |
| `tRNAs` | Number of unique standard tRNA types found |
| `MIMAG_quality` | MIMAG quality tier: High quality, Medium quality, or Low quality |
| `notes` | Empty unless the domain call was ambiguous; otherwise records which tools/domains were involved and how the ambiguity was resolved |

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
