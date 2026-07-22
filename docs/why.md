---
title: why
---

Why
========

# Why use Galah?

Metagenome-assembled genome (MAG) catalogues have grown from hundreds to hundreds of thousands of genomes.
Dereplicating a catalogue of this size (clustering genomes by average nucleotide identity (ANI) and picking one representative per cluster) is normally an all-vs-all comparison, which becomes impractical well before six-figure genome counts are reached.
Galah solves this issue through either a low-memory mode that avoids storing all pairwise comparisons in memory, solving the memory bottleneck, or by clustering against reference genomes instead of all-vs-all, avoiding the quadratic time complexity of pairwise comparisons entirely.

## Well-defined clusters

Generated cluster representatives have well-defined properties. If clustering at 95% ANI:

1. Each representative is \<95% ANI to every other representative.
2. Every member is \>=95% ANI to its representative.
3. If genome quality is available (e.g. from CheckM2), each representative has a better quality score than the other members of its cluster.

See the [Clustering](/#clustering) section of the docs for the full algorithm description.

## Benchmarking

We benchmarked galah against [skDER](https://github.com/raufs/skDER) (greedy clustering with skani), [dRep](https://drep.readthedocs.io/) (using skani), and [pyani](https://github.com/widdowquinn/pyani) (ANIm) at 100, 500, 1,000, 5,000, 10,000 and 50,000 genomes, using 32 CPUs, 250GB memory, and a 48-hour walltime limit.

![Benchmarking of galah against skDER, dRep and pyani across genome catalogue sizes](/benchmarking_summary_figure.png)

**Wall-clock time (a), CPU time (b) and peak memory (c) across catalogue sizes.**

pyani (ANIm) is dramatically slower than the other tools even at small scale, 1,000 genomes alone took over 20 hours wall-clock time, so it did not complete the larger benchmarks.
dRep also becomes comparatively slow at larger scale, and was unable to complete within the resource limits at 50,000 genomes.
galah, galah's low-memory mode, and skDER remain fast across the full range tested.
At 50,000 genomes, galah's low-memory mode uses substantially less peak memory than either galah's default mode or skDER, at some cost to wall-clock/CPU time.

Note that galah's preclustering step (grouping genomes roughly by ANI before fine-grained comparison) was previously important for performance with the older FastANI backend, but is no longer necessary now that skani is the default comparison method — skani itself is fast enough that comparing all genome pairs directly is practical at these scales.

## Dereplicating against an existing catalogue

Growing a genome catalogue over time usually means re-dereplicating everything from scratch against the new batch of genomes.
Galah's `--reference-genomes` option instead lets new genomes be clustered directly against an already-dereplicated reference set, so only input-vs-reference comparisons are needed, rather than all-vs-all.

We benchmarked this by clustering 100 to 50,000 input genomes against all 346,233 genomes in [GlobDB](https://globdb.org/) r232 as reference genomes, comparing galah's low-memory mode (clustering the input genomes together with all reference genomes) against galah's `--reference-genomes` mode (clustering the new genomes against the pre-existing reference set). No other tool benchmarked above was able to complete this comparison under the same resource limits.

![Benchmarking galah low-memory mode against galah --reference-genomes mode with 346,233 GlobDB r232 reference genomes](/benchmarking_summary_figure_globdb.png)

**Wall-clock time (a), CPU time (b) and peak memory (c) when clustering against 346,233 reference genomes.**
`--reference-genomes` is several-fold faster in both wall-clock and CPU time than clustering all genomes together, and uses substantially less peak memory, since it avoids repeating comparisons among the reference genomes themselves.

## Consistency of cluster assignments

Despite the large resource differences, cluster assignments are highly consistent between tools across both benchmarks, with the exception of pyani which uses a different clustering method:

![Normalised cluster counts relative to the mean across tools](/benchmarking_clusters_summary_figure.png)

**Cluster counts, normalised as (n − mean) / mean across tools, for (a) the benchmark above and (b) the reference-genome benchmark below.**
galah, galah (low memory), skDER and dRep all agree to within a fraction of a percent. pyani is the outlier, consistently producing more clusters than the ANI/skani-based tools. galah (low memory) and galah (reference genomes) are also highly consistent, with only a few clusters differing between the two methods.

## Eukaryotic genomes

Metagenomes increasingly yield eukaryotic MAGs (primarily microeukaryotes like protists and fungi) alongside bacterial and archaeal ones, but dereplication and quality-assessment pipelines are typically built around prokaryotes only: completeness/contamination models trained on bacterial/archaeal marker genes, and MIMAG rRNA/tRNA criteria (5S/16S/23S) that don't apply to eukaryotic genomes.
Mixed-domain catalogues have generally had to be split by domain and run through separate tools before being recombined.

Galah instead classifies each genome's domain automatically with [isiteuk](https://github.com/wwood/isiteuk), then routes it to the appropriate tools: CheckM2, Barrnap and tRNAscan-SE (bacterial/archaeal mode) for Bacteria/Archaea, and [EukCC](https://github.com/EBI-Metagenomics/EukCC) with eukaryotic rRNA (18S/28S/5.8S) and tRNA criteria for Eukaryota.
Genomes with no confident domain assignment are assessed under all applicable domains, and the `--domain-choice` flag can override automatic classification when the domain is already known.
This means a single `galah process` run can dereplicate and quality-assess a catalogue spanning all three domains together, using one cluster ranking based on each genome's domain-appropriate quality score.

## One tool for dereplication and genome quality

MIMAG quality (completeness, contamination, rRNA and tRNA presence) is normally assessed with a separate set of tools before or after dereplication.
Galah's `process` subcommand runs both steps together: it determines MIMAG quality and uses it to choose the best representative of each cluster, then reports both the cluster definition and the MIMAG summary from a single invocation.

## Integration with CoverM

Galah's clustering is also available through [CoverM](https://github.com/wwood/CoverM)'s `cluster` subcommand, so read coverage/relative abundance calculation and genome dereplication can be run from the same toolset.

> Aroney, S.T.N., Camargo, A.P., Tyson, G.W. and Woodcroft B.J. *Galah: More scalable dereplication for metagenome assembled genomes.* Zenodo (2024). <https://doi.org/10.5281/zenodo.13637856>
