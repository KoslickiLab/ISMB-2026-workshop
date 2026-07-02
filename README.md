# ISMB 2026 Workshop: k-mer Methods and Sketching

This tutorial covers [sourmash](https://sourmash.readthedocs.io/) and [YACHT](https://github.com/KoslickiLab/YACHT) for k-mer-based analysis of genomic and metagenomic data. It is designed for a one-hour session but includes additional material for self-directed exploration.

The two tutorials share a single dataset and tell one story. We start with a metagenomic sample and a handful of candidate reference genomes, and we work our way from raw sequence to a statistically supported answer to the question *which of these genomes are actually in my sample?*

- [Sourmash tutorial](Sourmash.md): what a sketch is, how to build and compare sketches, and why ANI is a better similarity measure than raw containment.
- [YACHT tutorial](YACHT.md): turning those same sketches into a hypothesis test for genome presence/absence, and how that differs from taxonomic profiling.

---

## Prerequisites

conda or miniconda must be installed. See the [miniconda installation guide](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) if you need to set this up.

---

## Setup

### 1. Create the conda environment

```bash
conda create -n ISMBtutorial -c bioconda -c conda-forge yacht=1.3.2
conda activate ISMBtutorial
```

`yacht` depends on `sourmash`, so this single environment covers both tutorials.

### 2. Create a working directory and download the data

Everything below runs from one working directory. Create it, download the demo dataset into it, and make a folder to hold the sketches we will build:

```bash
mkdir ISMBtutorial && cd ISMBtutorial
yacht download demo --outfolder ./demo
mkdir sketches
```

This gives you:

```
ISMBtutorial/
|-- demo/
|   |-- query_data/query_data.fq      # a metagenomic sample (reads)
|   |-- ref_genomes/                   # 15 candidate reference genomes (GTDB)
|   |-- toy_genome_to_taxid.tsv        # maps each genome to a taxonomy ID
|-- sketches/                          # (empty; we will fill this in)
```

The `demo` dataset is small and self-contained. The `query_data.fq` sample was simulated from 5 of the reference genomes at 0.5x coverage, so we have a known ground truth to check our answers against. The 15 reference genomes are randomly selected GTDB representatives.

All commands in both tutorials assume you are in the top-level `ISMBtutorial/` directory.

---

## Tutorials

Work through the tutorials in order:

1. [Sourmash tutorial](Sourmash.md)
2. [YACHT tutorial](YACHT.md)

---

## Optional: extra data for self-study

The classic sourmash tutorials use an *E. coli* dataset (reads, an assembled genome, and a collection of 50 strain signatures) and a pair of small synthetic genomes. These are not needed for this workshop, but if you want to work through the upstream [sourmash tutorial](https://sourmash.readthedocs.io/en/latest/tutorial-basic.html) on your own afterwards, you can grab them with:

```bash
# (optional, self-study only) run from inside ISMBtutorial/
mkdir -p extras && cd extras
curl -L https://osf.io/ruanf/download -o ecoliMG1655.fa.gz
curl -L https://osf.io/q472x/download -o ecoli_ref-5m.fastq.gz
curl -O -L https://github.com/sourmash-bio/sourmash/raw/latest/data/eschericia-sigs.tar.gz
tar xzf eschericia-sigs.tar.gz && rm eschericia-sigs.tar.gz
cd ..
```
