---
title: installation
---

Installation
========

There are several ways to install Galah

## Install from Bioconda via Pixi

Create pixi.toml file:

```toml
[workspace]
channels = ["conda-forge", "bioconda"]
name = "galah"
platforms = ["linux-64"]

[dependencies]
galah = "*"
```

Create pixi environment.

```bash
pixi install

# Either run within your current environment
pixi run galah -h
# Or enter the environment
pixi shell
```

One can see [details of the galah recipe](https://bioconda.github.io/recipes/galah/README.html).

Galah can also be used indirectly through
[CoverM](https://github.com/wwood/CoverM) via its `cluster` subcommand, which is also available on bioconda.

## Install from Bioconda via Conda

Install latest release via conda (or mamba).

```bash
conda create -n galah -c bioconda -c conda-forge galah

# Activate the environment
conda activate galah
```

### Pre-compiled binary

Galah can be installed by downloading statically compiled binaries, available on
the [releases page](https://github.com/wwood/galah/releases).

Third party dependencies listed below are required for this method.

### Compiling from source

Galah can also be installed from source, using the cargo build system after
installing [Rust](https://www.rust-lang.org/).

```
cargo install galah
```

Third party dependencies listed below are required for this method.

### Development

To run an unreleased version of Galah, after installing
[Rust](https://www.rust-lang.org/):

```
git clone https://github.com/wwood/galah
cd galah
pixi run cargo run -- cluster ...etc...
```

### Dependencies

Some usages of Galah require third party tools, which must be installed separately:

* skani v0.2.2 https://github.com/bluenote-1577/skani
* FastANI v1.34 https://github.com/ParBLiSS/FastANI
* Barrnap v0.9 https://github.com/tseemann/barrnap
* tRNAscan-SE v2.0.12 https://github.com/UCSC-LoweLab/tRNAscan-SE
* CheckM2 v1.1.0 https://github.com/chklovski/CheckM2
* isiteuk (for automatic domain classification) https://github.com/wwood/isiteuk
* EukCC v2 (for eukaryotic genome quality assessment) https://github.com/EBI-Metagenomics/EukCC

Most of these tools can be installed via pixi, using the `pixi.toml` file within the github repository.

```
pixi install
```

**Note:** `checkm2`, `isiteuk`, and `eukcc` cannot all be installed in the same conda environment
due to incompatible `diamond` and `zlib` dependency constraints between them.

Galah locates each tool using the following priority order:

1. `GALAH_CHECKM2_CMD` / `GALAH_ISITEUK_CMD` / `GALAH_EUKCC_CMD` environment variable (if set)
2. Tool found on `PATH`
3. Automatic installation via `pixi exec` (requires [pixi](https://pixi.sh) to be installed)

With pixi installed, no manual setup is needed — galah will download and cache each tool in an
isolated environment the first time it is needed.

To pin a specific version or override the command explicitly:

```bash
# Development (pixi workspace)
export GALAH_CHECKM2_CMD="pixi run -e checkm2 checkm2"
export GALAH_ISITEUK_CMD="pixi run -e isiteuk isiteuk"
export GALAH_EUKCC_CMD="pixi run -e eukcc eukcc"

# Production (separate conda envs)
conda create -n galah_checkm2 -c bioconda -c conda-forge checkm2
conda create -n galah_isiteuk -c bioconda -c conda-forge isiteuk
conda create -n galah_eukcc -c bioconda -c conda-forge eukcc
export GALAH_CHECKM2_CMD="conda run -n galah_checkm2 checkm2"
export GALAH_ISITEUK_CMD="conda run -n galah_isiteuk isiteuk"
export GALAH_EUKCC_CMD="conda run -n galah_eukcc eukcc"
```

#### CheckM2 database

CheckM2 requires a database to be set using the environment variable `CHECKM2DB` or the
argument `--checkm2-db-path`. See https://github.com/chklovski/CheckM2 for details.

```bash
export CHECKM2DB=/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

#### isiteuk metapackage

isiteuk uses a SingleM metapackage to classify genomes into biological domains. Set the
path using the environment variable `ISITEUK_METAPACKAGE_PATH` or the argument
`--isiteuk-metapackage`.

```bash
export ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg
```

#### EukCC database

EukCC requires a database for eukaryotic quality assessment. Set the path using the
environment variable `EUKCC2_DB` or the argument `--eukcc-db-path`. The database can be
downloaded from https://github.com/EBI-Metagenomics/EukCC.

```bash
export EUKCC2_DB=/path/to/eukcc2_db_ver_1.1
```
