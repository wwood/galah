# Galah Installation

Containerised, worked examples of installing [Galah](https://github.com/wwood/galah) on Linux
(x86_64) and smoke-testing each install with `galah --version` and `galah cluster` on the
release's own test genomes (`tests/data/set1`).

This folder is self-contained and can be split into its own repository
(`git subtree split --prefix galah-installation`).

| Method | Dockerfile template | What is installed |
|---|---|---|
| Pixi (Bioconda) | `Dockerfile.pixi.in` + `pixi.toml.in` | `galah` plus bundled third-party tools |
| Conda (Bioconda) | `Dockerfile.conda.in` | as above, via Miniconda |
| Cargo (crates.io) | `Dockerfile.cargo.in` | compiled `galah`; skani/FastANI/pixi from conda |

## Usage

Requires Docker.

```bash
bash compile_and_test_install_methods.bash RELEASE_VERSION   # e.g. 0.7.0
```

Each method is tested in two stages:

1. **Build** (`<method>.build.log`): install Galah, then `galah --version` and `galah cluster`.
2. **Quality** (`<method>.quality.log`, via `test_quality.sh`): download the CheckM2, isiteuk and
   EukCC databases, then run `galah analyse` and `galah process` on the three genomes in
   `tests/data/domain_examples` (bacterium, archaeon, eukaryote), checking each is assigned the
   expected domain. Galah's own automatic pixi installation of Barrnap/tRNAscan-SE (cargo method only) and CheckM2/isiteuk/EukCC is exercised
   here.

Databases are downloaded once into `./databases` (mounted at `/db`) and reused by the later
methods; expect >10 GB of disk and a long first run. Files there are written by root in the
container. Use `--skip-quality` to run stage 1 only.

Generated files (`Dockerfile.*`, `pixi.toml`) are git-ignored.
The release tag `vRELEASE_VERSION` must exist on GitHub, and the version must be on Bioconda and
crates.io for the respective methods.

## Scope and caveats

- `test_quality.sh` downloads databases with the latest CheckM2 and isiteuk from Bioconda (via
  `pixi exec`) and EukCC DB v1.1 from the EBI FTP. The assertions only check domain assignment
  and the output shape, not completeness/contamination values.
- The Cargo method installs only skani, FastANI and pixi explicitly (these must be on `PATH`);
  Barrnap, tRNAscan-SE, CheckM2, isiteuk and EukCC rely on Galah's automatic pixi installation.
- Bioconda builds are unpinned beyond Galah's version, so dependency versions drift over time.
  Record `pixi list` / `conda list` output if you need exact versions.

## Errors

Docker daemon not running: `sudo systemctl start docker` (Linux) or `sudo service docker start` (WSL).

Debug a failed build: `docker images`, then `docker run -it --entrypoint /bin/bash <image_id>`.
Clean up: `docker system prune -a`.
