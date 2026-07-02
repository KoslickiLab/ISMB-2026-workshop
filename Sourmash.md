# Sourmash Tutorial

## Setup

Activate the workshop environment and navigate to your working directory:

```bash
conda activate ISMBtutorial
cd ISMBtutorial
```

Everything below runs from the top-level `ISMBtutorial/` directory. If you have not downloaded the data yet, see the [README](README.md).

Our running example is a single dataset: a metagenomic sample (`demo/query_data/query_data.fq`) together with 15 candidate reference genomes (`demo/ref_genomes/`). The question we are building toward is: *which of those 15 genomes are actually present in the sample?* In this first half we build the sketching machinery with sourmash; in the [YACHT tutorial](YACHT.md) we turn it into a statistical test.

---

## What is a sketch?

A genome or a metagenome is a very large set of k-mers (all substrings of length `k`). Comparing these sets directly is exact but expensive. A **sketch** is a small, representative subsample of that set that preserves the quantities we care about (similarity, containment) while being orders of magnitude smaller.

Sourmash uses **FracMinHash**: hash every k-mer, and keep only those hashes that fall below a fraction `1/scaled` of the hash space. With `scaled=1000` we keep, on average, 1 out of every 1000 distinct k-mers. Because the same k-mer always hashes to the same value, two samples that share a k-mer will both keep it (or both drop it), so overlap between sketches estimates overlap between the full sets.

![FracMinHash sketch](Data/Picture6.png)

Let's sketch a single reference genome and look at what we get:

```bash
sourmash sketch dna demo/ref_genomes/GCF_018918045.1_genomic.fna.gz \
    -p k=31,scaled=1000 -o sketches/one_genome.sig.zip
```

| Parameter | Description |
|-----------|-------------|
| `sourmash sketch dna` | Build a sketch from DNA sequence (use `protein` for amino-acid sequence). Input can be FASTA or FASTQ, gzipped or not. |
| `k=31` | k-mer size. Larger k is more specific; smaller k is more sensitive. |
| `scaled=1000` | Compression factor. Keep ~1/1000 of the distinct k-mers. Larger `scaled` = smaller sketch. |
| `-o` | Output signature file. |

Now inspect it:

```bash
sourmash sig describe sketches/one_genome.sig.zip
```

```
k=31 molecule=DNA num=0 scaled=1000 seed=42 track_abundance=0
size: 2452
```

That genome is about 2.45 Mbp, so its sketch of **2452 hashes** is roughly a 1000-fold compression, yet it is enough to accurately estimate how this genome relates to any other sketch built the same way.

**When is a sketch the right tool?** When you need fast set-level comparisons at scale: is this genome contained in that metagenome, how similar are these two samples, which of a million references does my sample match. Sketching is *not* the right tool when you need base-level detail, for example calling SNPs, recovering exact coordinates, or assembling. For those you keep the reads. The rule of thumb: sketch when you are asking "how much overlap," keep the reads when you are asking "what exactly is at this position."

**Two sketches can only be compared if they were built with the same `k` and are compatible in `scaled`.** Keep those consistent across a project.

---

## Sketch the whole reference set and the sample

Sketch all 15 reference genomes into one file, and sketch the metagenomic sample:

```bash
# 15 reference genomes -> one signature file
sourmash sketch dna demo/ref_genomes/*.fna.gz \
    -p k=31,scaled=1000 --name-from-first -o sketches/refs.sig.zip

# the metagenomic sample (reads work exactly the same as genomes)
sourmash sketch dna demo/query_data/query_data.fq \
    -p k=31,scaled=1000 --name query_metagenome -o sketches/query.sig.zip
```

Note that sketching reads is no different from sketching a genome: FracMinHash just hashes k-mers, whatever their source. The `--name-from-first` flag labels each reference signature by its first sequence header (a readable organism name), and `--name` gives the sample a label of our choosing.

---

## Comparing sketches: containment and the ANI lesson

`sourmash compare` estimates the pairwise similarity between sketches. Two useful measures:

**Jaccard** ( `--jaccard` ) is the size of the intersection over the size of the union. It is symmetric, but it degrades when the two sets are very different in size (a small genome inside a large metagenome will have a tiny Jaccard even if it is fully present).

**Containment** ( `--containment` ) is the size of the intersection over the size of *one* of the sets: what fraction of A is found in B. It handles size asymmetry gracefully, which is exactly what we need for metagenomics.

![Jaccard and containment](Data/Picture8.png)

### A vignette: why ANI beats raw containment

Two of our reference genomes are isolates of the same species, *Cloacibacterium caeni*. They are similar but not identical, which makes them a good illustration. Sketch each at three k-mer sizes and measure containment at each:

```bash
sourmash sketch dna demo/ref_genomes/GCF_907163105.1_genomic.fna.gz \
    -p k=21,k=31,k=51,scaled=1000 -o sketches/cc_105.sig.zip
sourmash sketch dna demo/ref_genomes/GCF_907163125.1_genomic.fna.gz \
    -p k=21,k=31,k=51,scaled=1000 -o sketches/cc_125.sig.zip

sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --containment -k 21
sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --containment -k 31
sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --containment -k 51
```

The containment of one isolate in the other drops sharply as k grows:

```
k=21:  37.3%
k=31:  27.9%
k=51:  16.7%
```

Nothing about the two genomes changed; only k did. The reason is that containment depends on k: it falls off roughly as `containment ~ ANI^k`, because a single mismatch destroys every k-mer that spans it, and longer k-mers are more likely to span a mismatch. So a raw containment value is not comparable across runs that used different k, and on its own it is hard to interpret biologically.

The fix is to invert that relationship. **Average Nucleotide Identity (ANI)** estimates the per-base identity between the two genomes as `ANI ~ containment^(1/k)`, which is largely independent of k. Add `--estimate-ani`:

```bash
sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --estimate-ani -k 21
sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --estimate-ani -k 31
sourmash compare sketches/cc_105.sig.zip sketches/cc_125.sig.zip --estimate-ani -k 51
```

```
k=21:  95.6%
k=31:  96.1%
k=51:  96.6%
```

The containment more than halved from k=21 to k=51, but the ANI stays put at about **96%** at every k. This is why we report and reason about ANI rather than raw containment: it is stable across k and it means something concrete, namely two genomes that are about 96% identical, i.e. two members of the same species. (Sourmash can also report confidence intervals on ANI with `--estimate-ani-ci`; the theory behind these estimates is in Hera, Pierce-Ward, and Koslicki 2023 [[3]](#3).)

Hold on to that 96% number. When we get to YACHT, its training step will decide whether two references are "the same organism" using an ANI threshold, and these two *Cloacibacterium* isolates will be exactly the kind of pair that gets collapsed.

---

## Searching the sample against the references

Now the payoff question: what is in the metagenome? The right tool is `sourmash gather`, which greedily explains the sample as a combination of reference genomes (a minimum-metagenome-cover):

```bash
sourmash gather sketches/query.sig.zip sketches/refs.sig.zip
```

```
overlap     p_query p_match
---------   ------- -------
0.8 Mbp       30.5%   34.4%    NZ_JAHLQE010000140.1 Eubacterium sp. MSJ-13 ...
0.6 Mbp       21.2%   18.4%    NZ_JAHLPV010000001.1 Acetivibrio sp. MSJd-27 ...
0.5 Mbp       19.2%   20.5%    NZ_JAHLPU010000070.1 Pseudoflavonifractor sp. MSJ-30 ...
314.0 kbp     12.0%    9.0%    NZ_JAHLQA010000007.1 Blautia sp. MSJ-19 ...
156.0 kbp      5.9%    5.4%    NZ_JAHLPZ010000001.1 Eubacterium sp. MSJ-21 ...

found 5 matches total;
the recovered matches hit 88.7% of the query k-mers (unweighted).
```

Gather recovers exactly 5 of the 15 references, and these are precisely the 5 genomes the sample was simulated from. `p_query` is the fraction of the sample explained by each match; `p_match` is the fraction of that reference seen in the sample.

`sourmash search` answers a slightly different question, ranking references by how much of the sample each one contains individually:

```bash
sourmash search sketches/query.sig.zip sketches/refs.sig.zip --containment
```

```
4 matches above threshold 0.080; showing first 3:
similarity   match
----------   -----
 30.5%       NZ_JAHLQE010000140.1 Eubacterium sp. MSJ-13 ...
 21.2%       NZ_JAHLPV010000001.1 Acetivibrio sp. MSJd-27 ...
 19.2%       NZ_JAHLPU010000070.1 Pseudoflavonifractor sp. MSJ-30 ...
```

Notice that `search` applies a default containment threshold of 0.08, which is enough to drop the lowest-abundance organism (*Eubacterium* sp. MSJ-21, at 5.9%) entirely, and it shows only the top few hits unless you ask for more. A hard threshold like this is convenient, but it is a blunt instrument: it treats 8% as the line between present and absent regardless of genome size, coverage, or how much confidence you actually have.

This is a good place to stop and notice what we do *not* yet have. Our sample was simulated at only 0.5x coverage, so every present genome is seen only partially, and the containment values are correspondingly low (6% to 31%). Gather and search will happily rank matches, but they lean on hard thresholds (a 50 kbp overlap cutoff, a 0.08 containment cutoff), and they give us no principled statement about how confident we should be that a genome is truly present versus a stray k-mer coincidence. Answering *that* question, with control over the false positive and false negative rates, is exactly what YACHT is for.

---

## Optional: comparing many sketches at once

`sourmash compare` with a single argument builds the full all-vs-all matrix, and `sourmash plot` renders it as a dendrogram and heatmap:

```bash
sourmash compare sketches/refs.sig.zip --estimate-ani -k 31 -o sketches/refs_cmp
sourmash plot --pdf --labels sketches/refs_cmp
```

For our 15 references the matrix is mostly near-zero (they are distinct GTDB species), with a single hot off-diagonal pair: the two *Cloacibacterium caeni* isolates from the ANI vignette, at about 96% ANI. Everything else is a different species.

---

## Further exploration

- The canonical [sourmash tutorial](https://sourmash.readthedocs.io/en/latest/tutorial-basic.html) works through reads-vs-assembly containment and searching a 50-strain *E. coli* database (see the optional data in the [README](README.md)).
- `sourmash sketch` can track k-mer abundance ( `-p ...,abund` ), which lets `gather` report relative abundances. We use this in the YACHT tutorial to decorate presence/absence calls with abundances.
- FracMinHash sketches are a general-purpose tool: the same machinery underlies ANI estimation, cosine-similarity and other metric estimation, dN/dS estimation, and functional profiling. See the workshop slides for references.

---

## References

<a id="1">[1]</a> Irber, L., et al. (2024). sourmash v4: A multitool to quickly search, compare, and analyze genomic and metagenomic data sets. Journal of Open Source Software, 9(98), 6830. https://doi.org/10.21105/joss.06830

<a id="2">[2]</a> Koslicki, D., & Zabeti, H. (2019). Improving MinHash via the containment index with applications to metagenomic analysis. Applied Mathematics and Computation, 354, 206-215.

<a id="3">[3]</a> Hera, M. R., Pierce-Ward, N. T., & Koslicki, D. (2023). Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash. Genome Research, 33(7), 1061-1068.

---

# Please proceed to the [YACHT tutorial](YACHT.md)
