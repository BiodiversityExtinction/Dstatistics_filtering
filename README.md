Dstats Population Structure Toolkit

This repository provides a command-line script for filtering, summarising, and visualising D-statistics (ABBA–BABA) output from ANGSD.

The tool is designed to work directly with the block-jackknifed output produced by angsd -doAbbababa (i.e. tables containing H1, H2, H3, D-statistic estimates, standard errors, and Z-scores). It implements a topology-based framework to distinguish signal caused by incorrect population topology from signal potentially caused by introgression, and to summarise D-statistic results across multiple individuals.

────────────────────
Input files
────────────────────

1) ANGSD D-statistics output  
   A block-jackknifed D-statistics file produced using ANGSD (e.g. via jackKnife.R).  
   The file must contain at least the following columns:
   - H1, H2, H3 (sample IDs)
   - D or Dstat (or jackEst)
   - SE
   - Z

2) Sample-to-population mapping file  
   A plain-text file with two columns:
   - sample ID
   - population name  
   (whitespace- or tab-delimited)

   This mapping is used to assign samples to populations and removes any dependence on sample naming conventions.

────────────────────
Command-line parameters
────────────────────

Required:
- -i  Input ANGSD block-jackknifed D-statistics file
- -m  Sample-to-population mapping file
- -p  Output prefix

Optional:
- -z  Z-score threshold used to define “significant” results (default: 3)
- -R  Enable extraction and summarisation of standard (“real”) D-statistic comparisons
      where H1, H2, and H3 all belong to different populations

All analyses are run from a single command. Required R packages are automatically checked
and installed if missing.

────────────────────
Topology-based analyses
────────────────────

Using the population assignments, the script automatically tests all population
combinations and classifies each D-statistic into one of three categories:

- Correct topologies:
  H1 and H2 belong to the same population, H3 belongs to a different population.

- Incorrect topologies:
  H2 and H3 belong to the same population, H1 belongs to a different population.
  D and Z values are sign-flipped where necessary to enforce a consistent orientation.

- Within-population comparisons:
  H1, H2, and H3 all belong to the same population.

These categories are used to summarise and contrast D-statistic behaviour across
replicate individual-level tests.

────────────────────
Optional “real” D-statistic comparisons
────────────────────

When the -R flag is used, the script additionally extracts all D-statistics where
H1, H2, and H3 belong to three different populations.

For these comparisons:
- H1 and H2 are canonicalised (e.g. alphabetically) so that equivalent tests
  such as (A,B,C) and (B,A,C) are collapsed into a single orientation.
- D and Z values are sign-flipped when H1/H2 are swapped.
- Results are summarised across individuals for each population triple.

────────────────────
Outputs
────────────────────

Tabular outputs:
- <prefix>.topology.all.tsv
  All topology-classified D-statistic results.

- <prefix>.topology.summary.tsv
  Summary statistics per population pair and topology class, including:
  number of tests, mean and median D, mean and median Z, and proportion of tests
  exceeding the Z threshold.

- <prefix>.topology.wilcox.tsv
  Wilcoxon tests comparing Z distributions between correct and incorrect topologies
  for each population pair (with BH correction).

- <prefix>.real.tsv (optional)
  Canonicalised “real” D-statistic comparisons where all populations are distinct.

- <prefix>.real.summary.tsv (optional)
  Summary statistics for real D-statistic population triples.

────────────────────
Plots
────────────────────

All plots are written as PDF files suitable for inspection or inclusion as supplementary figures.

Topology-based plots:
- <prefix>.topology.Z_box_ALLpairs.pdf  
  Boxplots of Z-score distributions for all population pairs (correct, incorrect,
  and within-population). The y-axis is capped at the 99th percentile to prevent
  extreme values from obscuring structure.

- <prefix>.topology.heatmap_meanZ.pdf  
  A multi-page PDF containing mean Z heatmaps:
    page 1: correct topologies
    page 2: incorrect topologies
    page 3: within-population comparisons  
  Each page uses its own symmetric colour scale to preserve contrast even when
  some comparisons have very large values.

- <prefix>.topology.heatmap_propSig.pdf  
  Heatmaps showing the proportion of tests with |Z| greater than the chosen threshold,
  faceted by topology class.

Real D-statistic plots (optional):
- <prefix>.real.heatmap_meanZ.pdf  
  A multi-page PDF with one page per H3 population (top H3 populations by number of tests),
  showing mean Z values for canonicalised H1/H2 comparisons. Each page has its own
  colour scale.

────────────────────
Intended use
────────────────────

This tool is intended for exploratory population structure analyses, summarising
multiple individual-level D-statistic results, identifying consistent signals of
gene flow, and diagnosing topology-driven artefacts in D-statistics derived from ANGSD.
