# ISMB 2026 Workshop: k-mer Methods and Sketching

This tutorial covers [sourmash](https://sourmash.readthedocs.io/) and [YACHT](https://github.com/KoslickiLab/YACHT) for k-mer-based analysis of genomic and metagenomic data. It is designed for a one-hour session but includes additional material for self-directed exploration.

- [Sourmash tutorial](Sourmash.md): sketching, comparing, and searching genomic signatures
- [YACHT tutorial](YACHT.md): hypothesis-testing-based taxonomic profiling of metagenomic samples

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

---

### 2. Download tutorial data

Create a directory where you would like this work to take place and navigate into it. For example:

```bash
mkdir ISMBtutorial && cd ISMBtutorial
```

Then create the required subdirectories:

```bash
mkdir -p sourmash/{data,output,scripts} yacht
```

#### Sourmash tutorial data

From inside `ISMBtutorial/`:

```bash
cd sourmash/data

# Sample genomes for comparison exercises
wget https://raw.githubusercontent.com/Penn-State-Microbiome-Center/KickStart-Workshop-2024/main/Day3-Shotgun/Data/sample_001.fna
wget https://raw.githubusercontent.com/Penn-State-Microbiome-Center/KickStart-Workshop-2024/main/Day3-Shotgun/Data/sample_002.fna

# Metagenomic samples for the MetaGenomQuest section
wget -i https://raw.githubusercontent.com/Penn-State-Microbiome-Center/KickStart-Workshop-2022/main/Day5-Shotgun/Data/file_list.txt
ls *.gz | xargs -P6 -I{} gunzip {}

# E. coli reference and reads for the database search section
curl -L https://osf.io/ruanf/download -o ecoliMG1655.fa.gz
curl -L https://osf.io/q472x/download -o ecoli_ref-5m.fastq.gz

# Collection of 50 E. coli signatures for database construction
mkdir ecoli_many_sigs && cd ecoli_many_sigs
curl -O -L https://github.com/sourmash-bio/sourmash/raw/latest/data/eschericia-sigs.tar.gz
tar xzf eschericia-sigs.tar.gz && rm eschericia-sigs.tar.gz
cd ../..
```

#### YACHT tutorial data

From inside `ISMBtutorial/`:

```bash
cd yacht

# Query metagenome
wget https://github.com/KoslickiLab/YACHT/raw/refs/heads/main/demo/query_data/query_data.fq

# Sketch the sample
yacht sketch sample --infile query_data.fq --kmer 31 --scaled 1000 --outfile query_data.sig.zip

# Demo dataset for the YACHT tutorial walkthrough
yacht download demo --outfolder ./demo
cd ..
```

---

### 3. Download the GTDB reference database

From `ISMBtutorial/yacht/`, choose one option.

#### Option A: Download pre-sketched reference and train YACHT yourself

This is somewhat computationally intensive (just takes some time to iterate over all ~85K genomes; adjust threads as your system allows):

```bash
cd yacht

# Download the pre-sketched GTDB representative genomes (k=31)
yacht download default_ref_db --database gtdb --db_version rs214 --gtdb_type reps --k 31 --outfolder ./

# Train YACHT (adjust --num_threads to match your system)
yacht train --ref_file gtdb-rs214-reps.k31.zip --ksize 31 --num_threads 4 --ani_thresh 0.95 --prefix 'gtdb_ani_thresh_0.95' --outdir ./

# Run YACHT
yacht run --json gtdb_ani_thresh_0.95_config.json --sample_file query_data.sig.zip --num_threads 4 --min_coverage_list 0.01 --out test.xlsx
```

#### Option B: Download a pre-trained database

This is computationally lightweight, but Zenodo is painfully slow, so will likely take the same time as Option A, just with less resources:

```bash
cd yacht

# Download and unpack the pre-trained database
yacht download pretrained_ref_db --database gtdb --db_version rs214 --k 31 --ani_thresh 0.95 --outfolder ./
unzip gtdb-rs214-reps.k31_0.95_pretrained.zip

# The config JSON contains placeholder paths that must be replaced with the actual path on your system
PRETRAINED_DIR=$(realpath gtdb-rs214-reps.k31_0.95_pretrained)
sed -i "s|/Enter/Path/Here/|${PRETRAINED_DIR}/|g" \
    gtdb-rs214-reps.k31_0.95_pretrained/gtdb-rs214-reps.k31_0.95_config.json

# Run YACHT
yacht run \
    --json gtdb-rs214-reps.k31_0.95_pretrained/gtdb-rs214-reps.k31_0.95_config.json \
    --sample_file query_data.sig.zip \
    --num_threads 4 \
    --min_coverage_list 0.01 \
    --out test.xlsx
```

---

## Tutorials

Work through the tutorials in order:

1. [Sourmash tutorial](Sourmash.md)
2. [YACHT tutorial](YACHT.md)
