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

These tools can be installed via pixi, using the `pixi.toml` file within the github repository.

```
pixi install
```

Note that CheckM2 requires a database to be set using the environment variable `CHECKM2DB` or
the argument `--checkm2-db`.
See https://github.com/chklovski/CheckM2 for details.
