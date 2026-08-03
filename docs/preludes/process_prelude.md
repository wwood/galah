
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
