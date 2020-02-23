- [Galah](#galah)
  - [Installation](#installation)
    - [Install through the bioconda package](#install-through-the-bioconda-package)
    - [Pre-compiled binary](#pre-compiled-binary)
    - [Compiling from source](#compiling-from-source)
    - [Development](#development)
    - [Dependencies](#dependencies)
  - [Usage](#usage)
    - [Precluster ANI](#precluster-ani)
  - [License](#license)

![Galah image][galah]

By David Cook Wildlife Photography - originally posted to Flickr as Galah (Eolophus roseicapillus), CC BY 2.0, https://commons.wikimedia.org/w/index.php?curid=8388233

# Galah

[![Travis](https://img.shields.io/travis/wwood/galah.svg?style=flat-square)](https://travis-ci.org/wwood/galah)

Galah aims to be a more scalable metagenome assembled genome (MAG) dereplication
method. That is, it clusters microbial genomes together based on their average
nucleotide identity (ANI), and chooses a single member of each cluster as the
representative.

Galah uses a greedy clustering approach to speed up genome dereplication,
relative to e.g. [dRep](https://drep.readthedocs.io/), particularly when there
are many closely related genomes (i.e. >95% ANI). Generated cluster
representatives have 2 properties. If the ANI threshold was set to 99%, then:

1. Each representative is <99% ANI to each other representative.
2. All members are >=99% ANI to the representative.

If [CheckM](https://ecogenomics.github.io/CheckM/) genome qualities were
specified, then the clusters have an additional property:

3. Each representative genome has a better quality score than other members of
   the cluster. Each genome is assigned a quality score based on the formula
   "completeness - 4*contamination".

If CheckM qualities are not used, then:

1. Each representative genome was specified to galah before other members of the
   cluster.

The overall greedy clustering approach was largely inspired by the work of
Donovan Parks, as described in [this
publication](https://www.biorxiv.org/content/10.1101/771964v2.abstract). It
operates in 3 steps. In the first step, genomes are assigned as representative
if no genomes of higher quality are >99% ANI. In the second step, each
non-representative genome is assigned to the representative genome it has the
highest ANI with.

## Installation

### Install through the bioconda package

CoverM can be installed through the [bioconda](https://bioconda.github.io/user/install.html) conda channel. After initial setup of conda and the bioconda channel, it can be installed with

```
conda install galah
```
Galah can also be used indirectly through
[CoverM](https://github.com/wwood/CoverM), which is also available on bioconda.

### Pre-compiled binary

Galah can be installed by downloading statically compiled binaries, available on
the [releases page](https://github.com/wwood/Galah/releases).

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
cargo run -- cluster ...etc...
```
Third party dependencies listed below are required for this method.

### Dependencies

Galah relies on these 3rd party tools, which must be installed separately.

* Dashing v0.4.0 https://github.com/dnbaker/dashing
* FastANI v1.3 https://github.com/ParBLiSS/FastANI

## Usage
For clustering a set of genomes at 99% ANI:
```
galah cluster --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna >clusters
```
There are several other options for specifying genomes, ANI cutoffs, etc. See
`galah cluster --help` for more information.

### Precluster ANI
Similar to dRep, galah operates in two stages. In the first, a fast
pre-clustering distance ([dashing](https://github.com/dnbaker/dashing)) is
calculated between each pair of genomes. Genome pairs are only considered as
potentially in the same cluster with
[FastANI](https://github.com/ParBLiSS/FastANI) if the prethreshold ANI is
greater than the specified value. By default, the precluster ANI is set at 95%
and the final ANI is set at 99%.

## License

Galah is made available under GPL3+. See LICENSE.txt for details. Copyright Ben
Woodcroft.

Developed by Ben Woodcroft at the [Australian Centre for Ecogenomics](http://ecogenomic.org).

[galah]: Eolophus_roseicapilla_-Wamboin,_NSW,_Australia_-juvenile-8.smaller.jpg