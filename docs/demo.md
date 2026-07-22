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

We now have 8 genome files in total: 6 (near-)identical copies of the bacterial genome, 1 archaeal genome, and 1 eukaryotic genome.

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
  --threads 8 \
  --output-mimag-summary mimag_summary.tsv \
  --output-cluster-definition clusters.tsv
```

By default, Galah runs [isiteuk](https://github.com/wwood/isiteuk) to classify each genome by domain, then assesses quality with the tool appropriate to that domain: CheckM2, Barrnap and tRNAscan-SE for Bacteria/Archaea, and EukCC (with eukaryotic rRNA/tRNA criteria) for Eukaryota.
No domain needs to be specified manually — each genome is classified and routed automatically.

## Understanding the output

### `--output-mimag-summary`

`mimag_summary.tsv` reports one row per input genome, with its assigned domain, completeness, contamination, rRNA/tRNA gene counts, and overall MIMAG quality category:

```tsv
genome	domain	completeness	contamination	rRNA_5S	rRNA_16S	rRNA_23S	rRNA_18S	rRNA_28S	rRNA_5.8S	tRNAs	MIMAG_quality
GCF_002008365.1_genomic.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCF_002008365.1_genomic.dup1.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCF_002008365.1_genomic.dup2.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCF_002008365.1_genomic.dup3.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCF_002008365.1_genomic.dup4.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCF_002008365.1_genomic.dup5.fna.gz	Bacteria	100.00	0.85	1	1	1	0	0	0	20	High quality
GCA_003139855.1_genomic.fna.gz	Archaea	84.95	0.03	1	1	1	0	0	1	20	Medium quality
binchicken_co8412.34_euk.fna.gz	Eukaryota	99.25	1.88	7	0	0	0	0	0	17	Medium quality
```

- `domain`: the domain assigned by isiteuk (`Bacteria`, `Archaea`, or `Eukaryota`).
- `completeness`/`contamination`: from CheckM2 (Bacteria/Archaea) or EukCC (Eukaryota).
- `rRNA_*`/`tRNAs`: counts of the rRNA genes relevant to the assigned domain (5S/16S/23S for prokaryotes; 18S/28S/5.8S/5S for eukaryotes) and the total number of distinct standard tRNAs found, from Barrnap and tRNAscan-SE respectively.
- `MIMAG_quality`: the overall [MIMAG](https://doi.org/10.1038/nbt.3893) quality category determined from the above.

The bacterial genome reaches High quality (completeness ≥ 90%, contamination < 5%, all of 5S/16S/23S present, ≥ 18 tRNA types).
The archaeal and eukaryotic genomes both land at Medium quality instead: the archaeal genome's completeness (84.95%) is just under the 90% High-quality threshold, and the eukaryotic genome is one tRNA type short of the ≥ 18 required (17 found) despite otherwise-excellent completeness and contamination.

### `--output-cluster-definition`

`clusters.tsv` reports one `representative<TAB>member` line per input genome, clustered by ANI (95% by default) using the MIMAG quality scores above to pick the best representative of each cluster:

```tsv
binchicken_co8412.34_euk.fna.gz	binchicken_co8412.34_euk.fna.gz
GCA_003139855.1_genomic.fna.gz	GCA_003139855.1_genomic.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup1.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup2.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup3.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup4.fna.gz
GCF_002008365.1_genomic.fna.gz	GCF_002008365.1_genomic.dup5.fna.gz
```

The 6 (near-)identical bacterial genomes collapse into a single cluster (all sharing the same representative in the first column), while the archaeal and eukaryotic genomes each form their own single-genome cluster, since neither is within 95% ANI of anything else in the input.
So from 8 input genomes, `galah process` reports **3 representative genomes** in total — one per distinct organism.
