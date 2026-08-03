---
title: demo
---

Demo
========

This demo will guide you through running `galah process` on a small, mixed-domain set of genomes: from one bacterium, one archaeon, and one eukaryote.
The only prerequisite is that Galah is fully installed and set up as per the [installation](/installation) instructions, including the CheckM2 database, isiteuk metapackage, and EukCC database.

## Download the demo genomes

We will use three real genomes, taken from Galah's own test suite: a bacterial genome ([GCF_002008365.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002008365.1/)), an archaeal genome ([GCA_003139855.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003139855.1/)), and a eukaryotic genome recovered by coassembly ([binchicken_co8412.34](https://github.com/aroneys/binchicken)).

```bash
mkdir galah_demo
cd galah_demo

wget https://raw.githubusercontent.com/wwood/galah/refs/heads/main/tests/data/domain_examples/GCF_002008365.1_genomic.fna.gz
wget https://raw.githubusercontent.com/wwood/galah/refs/heads/main/tests/data/domain_examples/GCA_003139855.1_genomic.fna.gz
wget https://raw.githubusercontent.com/wwood/galah/refs/heads/main/tests/data/domain_examples/binchicken_co8412.34_euk.fna.gz
```

To demonstrate clustering (rather than just quality assessment), we also make 5 duplicate copies of the bacterial genome, to simulate multiple near-identical strains recovered from different samples:

```bash
for i in 1 2 3 4 5; do
  cp GCF_002008365.1_genomic.fna.gz GCF_002008365.1_genomic.dup${i}.fna.gz
done
```

To see how Galah handles a genome that genuinely contains sequence from two domains at once, we make a synthetic "chimera" genome by concatenating the archaeal and eukaryotic genomes together:

```bash
cat GCA_003139855.1_genomic.fna.gz binchicken_co8412.34_euk.fna.gz > chimera_arc_euk.fna.gz
```

We now have 9 genome files in total: 6 (near-)identical copies of the bacterial genome, 1 archaeal genome, 1 eukaryotic genome, and 1 synthetic archaeal/eukaryotic chimera.

## Run `galah process`

`galah process` runs both quality assessment and clustering in one command.
Set the paths to the CheckM2, isiteuk, and EukCC databases (see [installation](/installation)), then run:

```bash
export CHECKM2DB=/path/to/CheckM2_database/uniref100.KO.1.dmnd
export ISITEUK_METAPACKAGE_PATH=/path/to/isiteuk.smpkg
export EUKCC2_DB=/path/to/eukcc2_db_ver_1.1

galah process \
  --genome-fasta-files \
    GCF_002008365.1_genomic.fna.gz GCF_002008365.1_genomic.dup1.fna.gz \
    GCF_002008365.1_genomic.dup2.fna.gz GCF_002008365.1_genomic.dup3.fna.gz \
    GCF_002008365.1_genomic.dup4.fna.gz GCF_002008365.1_genomic.dup5.fna.gz \
    GCA_003139855.1_genomic.fna.gz \
    binchicken_co8412.34_euk.fna.gz \
    chimera_arc_euk.fna.gz \
  --threads 8 \
  --output-mimag-summary mimag_summary.tsv \
  --output-cluster-definition clusters.tsv
```

By default, Galah runs [isiteuk](https://github.com/wwood/isiteuk) to classify each genome by domain (based on marker gene content), then assesses quality with the tool appropriate to that domain: CheckM2, Barrnap and tRNAscan-SE for Bacteria/Archaea, and EukCC (with eukaryotic rRNA/tRNA criteria) for Eukaryota.
No domain needs to be specified manually — each genome is classified and routed automatically.

## Understanding the output

### `--output-mimag-summary`

`mimag_summary.tsv` reports one row per input genome, with its assigned domain, completeness, contamination, rRNA/tRNA gene counts, and overall MIMAG quality category:

```tsv
genome	domain	completeness	contamination	rRNA_5S	rRNA_16S	rRNA_23S	rRNA_18S	rRNA_28S	rRNA_5.8S	tRNAs	MIMAG_quality	notes
GCF_002008365.1_genomic.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCF_002008365.1_genomic.dup1.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCF_002008365.1_genomic.dup2.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCF_002008365.1_genomic.dup3.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCF_002008365.1_genomic.dup4.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCF_002008365.1_genomic.dup5.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality	
GCA_003139855.1_genomic.fna.gz	Archaea	84.95	0.03	1	1	1	0	0	1	20	Medium quality	
binchicken_co8412.34_euk.fna.gz	Eukaryota	99.25	1.88	7	0	0	0	0	0	17	Medium quality	
chimera_arc_euk.fna.gz	Archaea	100.00	34.26	8	1	1	0	0	1	20	Low quality	isiteuk assigned multiple domains: Eukaryota, Archaea; domain resolved to Archaea via higher completeness (CheckM2 100.00% vs EukCC 99.25%)
```

- `domain`: the domain assigned by isiteuk (`Bacteria`, `Archaea`, or `Eukaryota`).
- `completeness`/`contamination`: from CheckM2 (Bacteria/Archaea) or EukCC (Eukaryota).
- `rRNA_*`/`tRNAs`: counts of the rRNA genes relevant to the assigned domain (5S/16S/23S for prokaryotes; 18S/28S/5.8S/5S for eukaryotes) and the total number of distinct standard tRNAs found, from Barrnap and tRNAscan-SE respectively.
- `MIMAG_quality`: the overall [MIMAG](https://doi.org/10.1038/nbt.3893) quality category determined from the above.
- `notes`: empty unless isiteuk's domain call was multi-domain/ambiguous (no cutoff passed, or several passed at once). When ambiguous, both CheckM2 and EukCC are run and whichever reports higher completeness wins; this column then records that the domain was ambiguous and, if resolved, which tool's result (and completeness scores) won.

The bacterial genome reaches High quality (completeness ≥ 90%, contamination < 5%, all of 5S/16S/23S present, ≥ 18 tRNA types).
The archaeal and eukaryotic genomes both land at Medium quality instead: the archaeal genome's completeness (84.95%) is just under the 90% High-quality threshold, and the eukaryotic genome is one tRNA type short of the ≥ 18 required (17 found).

The chimera genome is exactly the case `--domain-choice isiteuk`'s ambiguity handling is for: isiteuk finds marker genes for both Eukaryota and Archaea in the same file (unsurprising, since it's literally the concatenation of an archaeal and a eukaryotic genome), so Galah falls back to running both CheckM2 and EukCC and keeping whichever reports higher completeness. Here CheckM2 calls it 100.00% complete against the Archaea/Bacteria marker set (the archaeal contigs plus the eukaryote contigs are enough to convince CheckM2 that the genome is complete) versus EukCC's 99.25% (which is actually just the completeness of the embedded eukaryotic genome), so the row is reported as Archaea, but with contamination inflated to 34.26% by the eukaryotic contigs CheckM2 doesn't expect to see — dragging its MIMAG quality down to Low despite the high completeness. The `notes` column records exactly what happened so this isn't mistaken for a straightforward domain call.

### `--output-cluster-definition`

`clusters.tsv` reports one `representative<TAB>member` line per input genome, clustered by ANI (95% by default) using the MIMAG quality scores above to pick the best representative of each cluster:

```tsv
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup1.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup3.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup2.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup4.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup5.fna.gz
GCA_003139855.1_genomic.fna.gz	GCA_003139855.1_genomic.fna.gz
GCA_003139855.1_genomic.fna.gz	chimera_arc_euk.fna.gz
binchicken_co8412.34_euk.fna.gz	binchicken_co8412.34_euk.fna.gz
```

The 6 (near-)identical bacterial genomes collapse into a single cluster, and the eukaryotic genome forms its own singleton — but the chimera genome clusters *with* the archaeal genome rather than forming a new cluster: its file contains the archaeal genome's contigs unmodified, so skani finds them close enough (well within 95% ANI over the aligned region) to group together, with `GCA_003139855.1_genomic.fna.gz` chosen as representative over the chimera itself since it has much lower contamination.
So even with the chimera added, `galah process` still reports **3 representative genomes** from the 9 input genomes — one per distinct organism.

## Domain choice options

`--domain-choice` controls how each genome's domain is decided, and defaults to `isiteuk` (used throughout this demo). Two other values are useful when a genome's domain is uncertain or you want to double-check isiteuk's call:

- `--domain-choice completeness` skips isiteuk and instead runs *both* CheckM2 and EukCC on every genome, keeping whichever reports the higher completeness as a single row. If CheckM2 wins, the domain is reported as `Bacteria,Archaea` (CheckM2 doesn't distinguish the two); rRNA and tRNA are then searched under Bacteria and Archaea kingdoms/modes and the *pair* from whichever single one scores higher overall is kept. This is also what happens automatically to any genome isiteuk can't confidently classify. Works with both `analyse` and `process`.
- `--domain-choice all` also runs both CheckM2 and EukCC on every genome, but reports **all three domains as separate rows** instead of picking a winner — useful for inspecting how a genome scores under every domain's criteria. Since this produces multiple rows per genome, there's no single quality value left to cluster on, so it's supported by `analyse` only; `process --domain-choice all` is rejected with an error.

For example, running `analyse --domain-choice all` on just the archaeal genome from this demo:

```bash
galah analyse \
  --genome-fasta-files GCA_003139855.1_genomic.fna.gz \
  --domain-choice all \
  --output-mimag-summary mimag_summary_all.tsv
```

produces three rows for that one genome, one per domain:

```tsv
genome	domain	completeness	contamination	rRNA_5S	rRNA_16S	rRNA_23S	rRNA_18S	rRNA_28S	rRNA_5.8S	tRNAs	MIMAG_quality	notes
GCA_003139855.1_genomic.fna.gz	Archaea	84.95	0.03	1	1	1	0	0	1	20	Medium quality	
GCA_003139855.1_genomic.fna.gz	Bacteria	84.95	0.03	1	1	1	0	0	0	18	Medium quality	
GCA_003139855.1_genomic.fna.gz	Eukaryota	10.60	0.92	1	0	0	1	0	0	19	Low quality	
```

The Bacteria and Archaea rows share the same completeness/contamination (CheckM2 doesn't distinguish between them), but their rRNA/tRNA counts differ since each row's Barrnap/tRNAscan-SE search is restricted to that row's own kingdom/mode.
The Eukaryota row shows EukCC finding almost no real eukaryotic signal, correctly landing at Low quality.
The exact Eukaryota-row numbers can vary a little between EukCC/database versions: since this genome isn't actually eukaryotic, EukCC's phylogenetic placement step has no confident clade to land on (expect something like `Genome belongs to clade: protozoa (Best TaxID: protist_common)` in its log).

## Dereplicating against a reference catalogue (e.g. GlobDB)

Rather than re-dereplicating an entire growing genome catalogue from scratch every time new genomes arrive, `--reference-genomes`/`--reference-genomes-list` let you cluster only the new genomes against an already-dereplicated reference set (see [Why use Galah?](/why#dereplicating-against-an-existing-catalogue) for benchmarking against [GlobDB](https://globdb.org/), a public catalogue of hundreds of thousands of dereplicated prokaryotic genomes).

The new genomes given via `--genome-fasta-files`/`--genome-fasta-directory`/`--genome-fasta-list` do **not** need to be dereplicated amongst themselves first, Galah does that automatically before comparing against the reference set. The reference genomes themselves must already be dereplicated at the chosen ANI, since reference-vs-reference comparisons are never made.

```bash
# reference_genomes.txt: one path per line to an already-dereplicated reference catalogue,
# e.g. downloaded from https://globdb.org/
galah cluster \
  --genome-fasta-directory new_genomes/ \
  --reference-genomes-list reference_genomes.txt \
  --precluster-ani 90 --ani 95 \
  --threads 8 \
  --output-cluster-definition clusters.tsv
```

Any of `new_genomes/`'s genomes within 95% ANI of each other collapse into one representative first; each resulting representative is then compared only against the reference catalogue (not against every other reference), so the cost scales with the number of *new* genomes rather than the size of the reference catalogue. A `new_genomes/` genome that doesn't match anything in the reference catalogue becomes a new representative of its own.
