<!-- NOTE: This intro should manually be kept in sync between the repo README and the docs README -->

[![Current Build](https://github.com/wwood/galah/actions/workflows/test-galah.yml/badge.svg)](https://github.com/wwood/galah/actions)
[![Conda version](https://img.shields.io/conda/v/bioconda/galah)](https://anaconda.org/bioconda/galah)
[![Conda downloads](https://img.shields.io/conda/d/bioconda/galah)](https://anaconda.org/bioconda/galah)
[![Crates.io version](https://img.shields.io/crates/v/galah)](https://crates.io/crates/galah)
[![Crates.io downloads](https://img.shields.io/crates/d/galah)](https://crates.io/crates/galah)

# Galah

[<img src="docs/_include/galah_logo.png" alt="Galah logo" width="600"/>](https://github.com/wwood/galah/blob/main/docs/_include/galah_logo.png)

Galah - Scalable dereplication and MIMAG calculation for metagenome assembled genomes.

Documentation can be found at [https://wwood.github.io/galah/](https://wwood.github.io/galah/).

Galah aims to be a scalable metagenome assembled genome (MAG) dereplication and quality assessment method.
Dereplication clusters genomes together based on their average nucleotide identity (ANI), and chooses a single member of each cluster as the representative.
Quality assessment results in a MIMAG quality score for each genome, based on its completeness, contamination and the presence of rRNA and tRNA genes.

## Quick install

```bash
# Install latest release via conda.
conda create -n galah -c bioconda -c conda-forge galah
```

## Example usage

For clustering and determining MIMAG quality scores:

```bash
galah process --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna \
  --output-cluster-definition clusters.tsv \
  --output-mimag-summary mimag.tsv
```

For clustering a set of genomes at 95% ANI:

```bash
galah cluster --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna \
  --output-cluster-definition clusters.tsv
```

For clustering a set of contigs at 95% ANI:

```bash
galah cluster --cluster-contigs --small-genomes --genome-fasta-files /path/to/contigs.fna \
  --output-cluster-definition clusters.tsv
```

For determining MIMAG quality scores for a set of genomes with CheckM2, Barrnap, and tRNAscan-SE:

```bash
galah analyse --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna \
  --output-mimag-summary mimag.tsv
```

## Help

If you have any questions or need help, please [open an issue](https://github.com/wwood/galah/issues).

## License
Galah is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions from [Samuel Aroney](https://github.com/AroneyS), [Antônio Camargo](https://github.com/apcamargo), and [Rhys Newell](https://github.com/rhysnewell). It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/galah](https://github.com/wwood/galah).

## Citation
<!-- NOTE: Citations should manually be kept in sync between the repo README and the docs README -->

Aroney, S.T.N., Camargo, A.P., Tyson, G.W. and Woodcroft B.J.
Galah: More scalable dereplication for metagenome assembled genomes.
Zenodo (2024). https://doi.org/10.5281/zenodo.13637856
