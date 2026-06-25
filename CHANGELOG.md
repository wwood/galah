# Changelog

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
