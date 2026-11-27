
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
  If EukCC cannot estimate quality (e.g. no placement marker genes are found), completeness
  and contamination are recorded as 0.

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

### Representative selection and MIMAG quality

By default, cluster representatives are chosen using the `Parks2020_reduced` quality formula, which
ranks genomes on completeness, contamination and assembly statistics, not MIMAG tier. A lower tier
genome can therefore be chosen over a higher tier genome in the same cluster, e.g. a Medium quality
genome that is slightly more complete than a High quality genome but is missing a 16S rRNA gene.
When a representative has a lower MIMAG tier than any member of its cluster, galah logs a warning,
and the `notes` column of the `--output-mimag-summary` records it for both the representative and
the higher quality genome(s) it was chosen over.

To rank genomes by MIMAG tier first (High, then Medium, then Low), and only then by
`Parks2020_reduced`, use `--quality-formula Parks2020_reduced_hqpreferred`. This formula is only
available in `process`, since it needs the MIMAG quality from the analyse step, and it cannot be
used with reference genomes, which are not assessed.

Even with this formula, a genome can occasionally end up in a cluster with a lower quality
representative, which is also reported in the warning and `notes` column. A genome becomes a
representative only if its ANI to every higher-ranked representative is below the threshold, but
each remaining genome is then assigned to the representative it has the highest ANI to. Since ANI
is not transitive, a High quality genome X can be within the threshold of High quality
representative A (so X is not a representative), while a Medium quality genome M is outside the
threshold of A (so M becomes a representative) but closer to X than A is, in which case X joins
M's cluster.

```bash
CHECKM2DB=CheckM2_database/uniref100.KO.1.dmnd \
galah process \
    --genome-fasta-files genome1.fna genome2.fna \
    --domain-choice bac \
    --quality-formula Parks2020_reduced_hqpreferred \
    --output-cluster-definition clusters.tsv \
    --output-mimag-summary mimag_summary.tsv
```

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
| High quality | Completeness ≥ 90%, contamination < 5%, all eukaryotic rRNAs (≥ 1 each of 18S, 28S, 5.8S and 5S), ≥ 18 tRNA types |
| Medium quality | Completeness ≥ 50%, contamination < 10% (and not High quality) |
| Low quality | Completeness < 50% or contamination ≥ 10% |

**Note:** MIMAG defines its quality criteria for bacteria and archaea only. Galah extends MIMAG to eukaryotes by requiring eukaryotic rRNA genes (18S, 28S, 5.8S and 5S) in place of the bacterial/archaeal 5S/16S/23S rRNAs for High quality, and otherwise keeping the thresholds the same (completeness, contamination and tRNA count).

See `galah analyse --full-help` for a full description of the `--output-mimag-summary` output format,
and `galah cluster --full-help` for clustering-specific parameters.
