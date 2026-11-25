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
| Cargo (crates.io) | `Dockerfile.cargo.in` | compiled `galah`; skani/FastANI from conda |

## Usage

Requires Docker.

```bash
bash compile_and_test_install_methods.bash RELEASE_VERSION   # e.g. 0.7.0
```

Generated files (`Dockerfile.*`, `pixi.toml`) are git-ignored; build output goes to `<method>.build.log`.
The release tag `vRELEASE_VERSION` must exist on GitHub, and the version must be on Bioconda and
crates.io for the respective methods.

## Scope and caveats

- The smoke test only runs `galah cluster`, which needs skani/FastANI. `galah analyse`/`process`
  additionally need the CheckM2, isiteuk and EukCC reference databases (large downloads), so
  they are not tested here. See the Galah docs on databases.
- The Cargo method installs only skani and FastANI; Barrnap, tRNAscan-SE, CheckM2, isiteuk and
  EukCC would need installing separately for `analyse`/`process`.
- Bioconda builds are unpinned beyond Galah's version, so dependency versions drift over time.
  Record `pixi list` / `conda list` output if you need exact versions.

## Errors

Docker daemon not running: `sudo systemctl start docker` (Linux) or `sudo service docker start` (WSL).

Debug a failed build: `docker images`, then `docker run -it --entrypoint /bin/bash <image_id>`.
Clean up: `docker system prune -a`.
