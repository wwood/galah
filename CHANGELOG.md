# Changelog

## Unreleased

### Added
- Multi-domain genome quality assessment: `galah analyse` and `galah process` now non-exclusively classify each genome by domain (Bacteria, Archaea, or Eukaryota) before choosing the appropriate quality tool and RNA criteria
- `--domain-choice isiteuk` (default): runs [isiteuk](https://github.com/wwood/isiteuk) to classify genomes by domain; genomes with no confident domain assignment are assessed under all three domains
- EukCC support for eukaryotic genome quality (completeness/contamination); run automatically for genomes classified as Eukaryota
- Eukaryotic rRNA and tRNA criteria: 18S/28S/5.8S/5S rRNA and ≥18 tRNAs required for high-quality eukaryotic MAGs
- `--isiteuk-bacteria-cutoff`, `--isiteuk-archaea-cutoff`, `--isiteuk-eukaryota-cutoff` arguments to tune the minimum isiteuk marker count for domain assignment (defaults: 10, 10, 20)
- `--isiteuk-output` to supply a pre-computed isiteuk TSV and skip running isiteuk
- `--isiteuk-metapackage` / `ISITEUK_METAPACKAGE_PATH` environment variable support
- `--eukcc-db-path` / `EUKCC2_DB` environment variable support
- `--eukcc-quality-report` to supply a pre-computed EukCC TSV and skip running EukCC
- Bundled pixi environments for CheckM2, isiteuk, and EukCC: tools are installed automatically on first use if not found on `PATH`
- `--skip-input-dereplication` (for `cluster`/`process`): when used with `--reference-genomes`/`--reference-genomes-list`, restores the behaviour prior to this release by comparing every input genome directly against the reference set instead of dereplicating input genomes amongst themselves first
- `--skip-sanitise-headers` (for `cluster`/`process`): skip checking/rewriting FASTA headers containing tabs before running skani, passing genome paths straight through unchanged. Mainly intended for benchmarking against tools which do not perform this sanitizing step; if any input genome actually has a tab in a header line, skani's TSV output will be silently corrupted, so this should only be used on genome sets already known not to have tabs in their headers

### Changed
- Domain output column in MIMAG summary now reflects the assigned domain(s); multi-domain genomes show comma-separated values (e.g. `Bacteria,Archaea`)
- Barrnap and tRNAscan-SE are now run in the mode matching each genome's assigned domain
- `--reference-genomes`/`--reference-genomes-list` no longer require input genomes to be pre-dereplicated amongst themselves before being clustered against the reference set: input genomes are now dereplicated amongst themselves first automatically, and only the resulting representative(s) are then matched against the reference genomes (which must still already be dereplicated, since reference-vs-reference comparisons are never made)
- Significantly sped up `cluster`/`process` on large genome sets by parallelizing FASTA header sanitizing (previously ran single-threaded regardless of `--threads`, dominating runtime for tens of thousands of genomes), buffering the sanitized-copy writes, and skipping the sanitizing step entirely for genomes whose headers contain no tab characters (the common case) so they are passed straight to skani without being copied at all

### Fixed
- `.gz` genome inputs built from multiple concatenated gzip streams (e.g. `cat a.fna.gz b.fna.gz > combined.fna.gz`) were silently truncated to just the first stream before being passed to isiteuk, CheckM2, EukCC or tRNAscan-SE, since decompression used `GzDecoder` (single-member) instead of `MultiGzDecoder` (multi-member, per RFC 1952)

## Version 0.5.2

### Changed
- Release procedure

## [0.5.1] - 2026-06-26

### Fixed
- Re-release for packaging

## [0.5.0] - 2026-06-26

### Added
- `galah analyse` subcommand for determining MIMAG quality scores
- `galah process` subcommand combining analyse and cluster functions
- `--reference-genomes` argument for clustering against existing reference genomes, reducing ANI comparisons
- `--low-memory` flag substantially decreasing memory requirements during clustering

### Fixed
- Rare contig-clustering bug when transitive property is not satisfied

### Removed
- Dashing preclusterer

### Changed
- Whitespace is now stripped from genome paths and contig names

## [0.4.2] - 2024-09-03

### Changed
- Documentation updates (thanks [@solc42](https://github.com/solc42))

## [0.4.1] - 2024-09-02

### Changed
- Documentation updates (thanks [@solc42](https://github.com/solc42))
- Added DOI to citation

## [0.4.0] - 2024-01-18

### Added
- `--checkm2-quality-file` option for cluster subcommand

### Changed
- Updated to clap v4 command line parser
- skani is now the default genome comparison tool instead of fastANI

## [0.3.1] - 2021-11-26

### Fixed
- Improved symlinking procedure for Windows (thanks [@apcamargo](https://github.com/apcamargo))
- Updated coverm-rs to allow absolute paths in CheckM file (reported by [@rhysnewell](https://github.com/rhysnewell))

## [0.3.0] - 2020-12-11

### Fixed
- Fixed argument parsing for `--precluster-ani` (thanks [@apcamargo](https://github.com/apcamargo))

### Changed
- Added troubleshooting guidance for dashing installation (thanks Rafael Laso Pérez)
- Updated to bird_tool_utils v0.3.0

## [0.2.0] - 2020-08-26

### Added
- `--quality-formula` parameter with `Parks2020_reduced` as default
- `--output-representative-fasta-directory-copy` and `--output-representative-list` output options

### Fixed
- Aligned fraction calculation now computed from fragment counts rather than relying on FastANI's thresholds

### Changed
- Updated FastANI dependency to v1.31
- Renamed `--prethreshold-ani` argument to `--precluster-ani`
- Overall full help text for `cluster` mode

## [0.1.0] - 2020-02-20

- Initial release
