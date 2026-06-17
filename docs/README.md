Galah
=====

Galah is a scalable dereplication and MIMAG calculation tool for metagenome assembled genomes.

Galah aims to be a more scalable metagenome assembled genome (MAG) dereplication
method. That is, it clusters microbial genomes together based on their average
nucleotide identity (ANI), and chooses a single member of each cluster as the
representative.

Galah also determines MIMAG quality scores for genomes based on their
completeness, contamination and the presence of rRNA and tRNA genes.

## Clustering

Galah uses a greedy clustering approach to speed up genome dereplication,
relative to e.g. [dRep](https://drep.readthedocs.io/), particularly when there
are many closely related genomes (i.e. >95% ANI). Generated cluster
representatives have 2 properties. If the ANI threshold was set to 95%, then:

1. Each representative is <95% ANI to each other representative.
2. All members are >=95% ANI to the representative.

If `--run-checkm2` was specified, or [CheckM2](https://github.com/chklovski/CheckM2) /
[CheckM](https://ecogenomics.github.io/CheckM/) genome qualities were provided,
then the clusters have an additional property:

3. Each representative genome has a better quality score than other members of
   the cluster. Each genome is assigned a quality score based on the formula
   `completeness-5*contamination-5*num_contigs/100-5*num_ambiguous_bases/100000`,
   which is reduced from a quality formula described in
   Parks et. al. 2020 https://doi.org/10.1038/s41587-020-0501-8.
   Other quality score formula are available via `--quality-formula`.

If instead CheckM1/2 qualities are not available, then the following holds instead:

3. Each representative genome was specified to Galah before other members of the
   cluster.

The overall greedy clustering approach was largely inspired by the work of
Donovan Parks, as described in [Parks et. al. 2020](https://doi.org/10.1038/s41587-020-0501-8).
It operates in 3 steps. In the first step, genomes are assigned as representative
if no genomes of higher quality are >95% ANI. In the second step, each
non-representative genome is assigned to the representative genome with which it
has the highest ANI.

## Example usage

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

For determining MIMAG quality scores for a set of genomes with CheckM2:

```bash
galah analyse --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna \
  --output-mimag-summary mimag.tsv
```

For clustering and determining MIMAG quality scores:

```bash
galah process --genome-fasta-files /path/to/genome1.fna /path/to/genome2.fna \
  --output-cluster-definition clusters.tsv --output-mimag-summary mimag.tsv
```

## Help

If you have any questions or need help, please [open an issue](https://github.com/wwood/galah/issues).

## License

Galah is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions from [Samuel Aroney](https://github.com/AroneyS), [Antônio Camargo](https://github.com/apcamargo), and [Rhys Newell](https://github.com/rhysnewell). It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/galah](https://github.com/wwood/galah).

## Citation

<!-- NOTE: Citation should manually be kept in sync between the repo README and the docs README -->

Aroney, S.T.N., Camargo, A.P., Tyson, G.W. and Woodcroft B.J.
Galah: More scalable dereplication for metagenome assembled genomes.
Zenodo (2024). [https://doi.org/10.5281/zenodo.13637856](https://doi.org/10.5281/zenodo.13637856)
