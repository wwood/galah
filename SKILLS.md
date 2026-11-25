---
name: galah
description: Use Galah to dereplicate genomes/MAGs or contigs by ANI (galah cluster), to assess MIMAG quality of bacterial, archaeal and eukaryotic genomes (galah analyse), or both in one run (galah process).
---

# Using Galah

Galah dereplicates genomes by average nucleotide identity (ANI) and assigns MIMAG quality tiers.
Full documentation: https://wwood.github.io/galah/ (source in `docs/`). Use `galah <subcommand> --full-help`
for the complete, version-accurate argument list.

## Choosing a subcommand

| Goal | Subcommand | Needs databases? |
|---|---|---|
| Dereplicate genomes or contigs by ANI | `galah cluster` | No (unless `--run-checkm2`) |
| MIMAG quality tier per genome (completeness, contamination, rRNA, tRNA) | `galah analyse` | Yes |
| Quality assessment, then dereplication using that quality to pick representatives | `galah process` | Yes |
| Check an existing cluster file against a stricter ANI | `galah cluster-validate` | No |

## Specifying genomes

All subcommands accept one of:

- `-f/--genome-fasta-files a.fna b.fna ...`
- `-d/--genome-fasta-directory dir/` with `-x/--genome-fasta-extension` (default `fna`, so
  `.fa`/`.fasta` directories need `-x`)
- `--genome-fasta-list genomes.txt` (one path per line; best for large sets)

Threads default to 1: always set `-t/--threads`.

## galah cluster

```bash
galah cluster -t 16 --genome-fasta-list genomes.txt \
    --checkm2-quality-report quality_report.tsv \
    --output-cluster-definition clusters.tsv \
    --output-representative-list reps.txt
```

- Defaults: `--ani 95`, `--min-aligned-fraction 15`, `--precluster-ani 90`, skani for both
  preclustering and clustering.
- Representative choice uses genome quality from `--checkm2-quality-report`, `--checkm-tab-table`,
  `--genome-info` or `--run-checkm2`, scored with `--quality-formula` (default `Parks2020_reduced`).
  **Without any quality input, the first genome in input order is chosen as representative** (a
  warning is logged). `--min-completeness`/`--max-contamination` filter genomes before clustering.
- Outputs: `--output-cluster-definition` (TSV, `representative<TAB>member`, no header; singletons
  are their own representative), `--output-representative-list`,
  `--output-representative-fasta-directory` (symlinks) or `...-directory-copy` (copies).
- Contigs (e.g. viruses, plasmids) within one FASTA: `--cluster-contigs` plus exactly one of
  `--small-contigs` (< 20 kb) or `--large-contigs`. Quality filtering is skipped in this mode.
  Consider `--both-min-aligned-fraction` so a small sequence does not cluster with a larger one
  that contains it.
- Small genomes (< 20 kb) given as separate files: `--small-genomes`.
- New genomes against an already-dereplicated reference set: `--reference-genomes-list refs.txt`.
  Inputs are dereplicated amongst themselves first, then representatives are compared to the
  references (add `--skip-input-dereplication` to compare every input directly). If quality is
  supplied, it must cover the reference genomes too.
- Large sets / limited RAM: `--low-memory`.

## galah analyse

```bash
export CHECKM2DB=/path/to/CheckM2_database/uniref100.KO.1.dmnd
export ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk_data/isiteuk-backpack-0.0.1
export EUKCC2_DB=/path/to/eukcc2_db_ver_1.1
galah analyse -t 16 --genome-fasta-list genomes.txt \
    --working-dir work/ \
    --output-mimag-summary mimag_summary.tsv
```

- Domain is set by `--domain-choice`: `isiteuk` (default, classifies each genome), `bac`, `arc`,
  `euk`, `completeness` (runs CheckM2 and EukCC, keeps the higher-completeness result) or `all`
  (one row per domain per genome).
- Bacteria/Archaea: CheckM2 + Barrnap + tRNAscan-SE (prokaryotic modes). Eukaryota: EukCC +
  Barrnap + tRNAscan-SE (eukaryotic modes). If EukCC cannot place a genome, completeness and
  contamination are recorded as 0.
- If the domain is known, `--domain-choice bac` (or `arc`/`euk`) skips isiteuk and the unused
  databases.
- Reuse earlier results instead of re-running tools: `--isiteuk-output`, `--checkm2-quality-report`
  or `--checkm-tab-table`, `--eukcc-quality-report`, `--barrnap-gff-list`, `--trnascan-out-list`
  (the last two are headerless two-column TSVs mapping genome path, as given in input, to output
  path).
- `--working-dir DIR` keeps intermediate outputs and reuses completed steps on re-run, which makes
  interrupted runs resumable. Without it a temporary directory is used and deleted on exit.
- `--output-mimag-summary` columns: `genome`, `domain`, `completeness`, `contamination`,
  `rRNA_5S`, `rRNA_16S`, `rRNA_23S`, `rRNA_18S`, `rRNA_28S`, `rRNA_5.8S`, `tRNAs`,
  `MIMAG_quality`, `notes`. `--output-quality-report` writes a CheckM2-format report.
- MIMAG tiers: High = completeness ≥ 90%, contamination < 5%, rRNAs present (5S/16S/23S for
  prokaryotes; any of 18S/28S/5.8S/5S for eukaryotes) and ≥ 18 tRNA types; Medium = completeness
  ≥ 50%, contamination < 10%; otherwise Low.

## galah process

Takes the union of `analyse` and `cluster` arguments, and outputs both `--output-mimag-summary`
and the cluster outputs. Quality from the analyse step is used to pick representatives.
`--domain-choice all` is **not** supported (no single quality value per genome); use
`galah analyse --domain-choice all` for that.

## Environment and databases

- The Bioconda package on Linux brings skani, FastANI, Barrnap, tRNAscan-SE, isiteuk and pixi.
  CheckM2 and EukCC conflict with isiteuk, so Galah installs each into its own pixi environment the
  first time it is needed.
- Tool lookup order: `GALAH_CHECKM2_CMD` / `GALAH_ISITEUK_CMD` / `GALAH_EUKCC_CMD`, then `PATH`,
  then automatic `pixi run` install (Linux only).
- Databases (only for `analyse`/`process`): `CHECKM2DB` or `--checkm2-db-path`;
  `ISITEUK_METAPACKAGE_PATH` or `--isiteuk-metapackage` (download with
  `isiteuk data --output-directory DIR`); `EUKCC2_DB` or `--eukcc-db-path`.
- macOS (osx-arm64): no isiteuk or EukCC, so use `--domain-choice bac`/`arc` for
  `analyse`/`process`. `cluster` works as normal.

## Caveats

- 95% ANI with 15% minimum aligned fraction approximates species-level clusters; change `--ani`
  for strain-level or genus-level grouping and report the thresholds used.
- Cluster membership depends on input set and order (greedy, quality-ranked); representatives can
  change when genomes are added or the quality source changes.
- Record the Galah version (`galah --version`) along with the CheckM2, EukCC and isiteuk database
  versions used.
- FASTA headers containing tabs are rewritten before skani is run; do not use
  `--skip-sanitise-headers` unless headers are known to be clean.

## Developing Galah

See `CLAUDE.md` for test conventions: `pixi run cargo test -- <test_name>` for a single test,
`pixi run cargo test` for the full fast suite, and `bash tests/run_expensive_tests_at_cmr.sh` for
the expensive tests. Update `docs/` alongside any user-facing change.
