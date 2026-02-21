# dstats_poptests_report.sh

A bash + R reporting tool that summarises ANGSD block-jackknifed D-statistics into population-level topology tests, distinct-population (“real”) triplet summaries, and an optional cladogram-style topology summary based on aggregated “who tends to be sister?” votes across population triplets.

This tool is designed for datasets with few individuals per population, where classic phylogenetic branch lengths and traditional bootstraps can be misleading. It focuses on robust summaries and topology tendencies rather than divergence time inference.

---

## Overview

Given:

1. An ANGSD block-jackknifed D-statistics file (with columns including H1, H2, H3, SE, Z and jackEst/Dstat/D)
2. A sample-to-population mapping file

The script produces:

- Population-level topology tests (Right / Wrong / Within patterns)
- Summary tables with D and Z statistics
- Diagnostic plots (Z distributions and heatmaps)
- Optional distinct-population triple summaries
- Optional cladogram based on aggregated triplet sister votes

---

## Conceptual summary

### Topology test (Right / Wrong / Within)

For each ordered population pair:

- One population is treated as the “double” population (H1/H2).
- Another is treated as the “single” population (H3).

The script extracts rows from the ANGSD file that match:

- Right: H1,H2 in DoubPop and H3 in SingPop
- Wrong: configurations implying the alternative topology (with sign correction applied)
- Within: all three individuals from the same population

It summarises Z and D distributions per population pair and test type.

### Real (distinct-population) triplets

If enabled, the script extracts only comparisons where H1pop, H2pop and H3pop are all distinct.

- H1pop/H2pop are canonicalised (alphabetically sorted).
- D and Z are sign-flipped when H1/H2 are swapped.
- Each population triple is summarised consistently.

### Cladogram: “who tends to be sister?”

For each population triplet (A,B,C):

- Evaluate the three possible resolutions:
  - (A,B)|C
  - (A,C)|B
  - (B,C)|A
- Select the resolution with the smallest |mean_Z| (i.e. most tree-like for that trio).
- Record the winning sister pair for that trio.

These per-trio “sister” votes are aggregated across all triplets to build a cladogram-style topology.

Important:
- This is not a phylogenetic distance tree.
- Branch lengths are not meaningful.
- It is a topology summary based on triplet tendencies.

---

## Interpretation cautions

1. Z values can be inflated if SE is very small (often due to missing data or low informative site counts).
2. |Z| > threshold should be treated as a heuristic diagnostic.
3. The optional resampling is triplet resampling, not classical sequence bootstrapping.
4. Branch lengths in the tree output must not be interpreted biologically.

---

## Requirements

Command line:
- bash
- awk
- sort / coreutils

R:
- dplyr
- tidyr
- ggplot2
- scales
- viridisLite
- ape (if tree mode enabled)

The script attempts to install missing packages from CRAN if needed.

---

## Inputs

### 1. ANGSD D-statistics file (-i)

Must include header columns:
- H1
- H2
- H3
- SE
- Z
- jackEst or Dstat or D

### 2. Sample-to-population map (-m)

Whitespace-delimited 2-column file:

SampleID  Population

Example:

Ind01  Saitama  
Ind02  Saitama  
Ind03  Hokuriku  
Ind04  Shikoku  

Lines beginning with # are ignored.

---

## Usage

dstats_poptests_report.sh -i dstats_jackknifed.txt -m sample_to_pop.txt -p out_prefix [options]

Required:
- -i  ANGSD D-stat file
- -m  Sample-to-population map
- -p  Output prefix

Optional:
- -z  Z threshold (default 3)
- -R  Also output distinct-population triplet summaries
- -T  Build cladogram from triplet votes (requires -R)
- -B  Number of triplet resampling replicates (requires -T)

---

## Example workflows

Basic topology summaries:

./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p results/topology \
  -z 3

With real triplet summaries:

./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p results/real \
  -R \
  -z 3

With cladogram:

./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p results/tree \
  -R -T \
  -z 3

With triplet resampling:

./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p results/tree_boot \
  -R -T \
  -B 200 \
  -z 3

---

## Output files

Topology outputs:
- <prefix>.topology.all.tsv
- <prefix>.topology.summary.tsv
- <prefix>.topology.wilcox.tsv
- <prefix>.topology.Z_box_ALLpairs.pdf
- <prefix>.topology.heatmap_meanZ.pdf
- <prefix>.topology.heatmap_propSig.pdf

Real outputs (with -R):
- <prefix>.real.tsv
- <prefix>.real.summary.tsv
- <prefix>.real.heatmap_meanZ.pdf

Tree outputs (with -T):
- <prefix>.quartet_votes.tsv
- <prefix>.quartet_support.tsv
- <prefix>.quartet.nwk
- <prefix>.quartet.pdf

Resampling outputs (with -B):
- <prefix>.quartet.bootstrap_support.tsv
- <prefix>.quartet.bootstrap_trees.nwk

---

## Recommended reporting wording

Instead of:
“Bootstrap support”

Use:
“Triplet-vote resampling support”

Instead of:
“Phylogenetic tree”

Use:
“Cladogram summarising aggregated triplet sister tendencies”

---


---

## Issues

When reporting issues, include:
- Command used
- First ~20 lines of input files (anonymised if necessary)
- Full error message
- Relevant summary TSV outputs
