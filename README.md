dstats_poptests_report.sh

A small, self-contained pipeline to summarise ANGSD block-jackknifed D-statistics at the population level, produce diagnostic plots, and (optionally) build a cladogram-style topology summary based on “who tends to be sister?” votes across population triplets.

This tool is designed for situations where you have few individuals per population and you want robust, interpretable population-level summaries rather than over-interpreting per-comparison p-values, classic “bootstrap” support, or branch lengths.

What this tool does

Given:

An ANGSD block-jackknifed D-statistics file (with columns like H1 H2 H3 ... SE Z jackEst or Dstat or D)

A sample → population mapping file (2 columns)

It will:

A) Topology-style population tests (all ordered population pairs)

For every ordered pair of populations:

Treat one population as the “double” population (H1/H2)

Treat another as the “single” population (H3)

Then, from the ANGSD D-stat file it extracts any rows that match either:

Right pattern: H1,H2 in DoubPop and H3 in SingPop

Wrong patterns: H1 in SingPop and H2,H3 in DoubPop (including swapped variants)

Where needed, the script sign-flips D and Z so that “wrong” patterns are comparable to “right” patterns.

It writes all extracted rows to a tidy TSV and summarises them to:

per-pair distributions of Z and D

proportion of “significant” tests using a threshold (default |Z| > 3)

a Wilcoxon comparison of Right vs Wrong Z distributions (BH-adjusted), as a diagnostic (not a definitive hypothesis test)

It also produces:

a boxplot PDF of Z distributions for all population pairs

heatmaps of mean Z and proportion significant by pair and test type

B) “Real” distinct-population triple comparisons (optional)

Optionally, the tool can also extract only “real” comparisons where:

H1pop, H2pop, H3pop are all distinct

H1pop/H2pop are canonicalised (sorted) so the same triple is represented consistently

D and Z are sign-flipped if H1/H2 are swapped, so summaries are consistent

It then summarises these “real” triples and produces a multi-page heatmap PDF (one page per H3pop for the most common H3 values).

C) Optional cladogram-style topology summary from triplet votes (optional)

If enabled, the script builds a cladogram-style topology summary that answers:

Across many population triplets (A,B,C), which pair most often behaves like the “sister” pair?

Important: this is not a phylogenetic tree with meaningful branch lengths. It is a summary of local triplet resolutions (“within triplets”) aggregated across all triplets.

Internal logic:

For each trio of populations (A,B,C), there are three possible “sister” resolutions:

(A,B)|C

(A,C)|B

(B,C)|A

Each resolution corresponds to a set of D-stat comparisons summarised in the real-summary table.

The tool assigns a score to each resolution (often |mean_Z|) and picks the smallest score as that trio’s “winner”.

The winner becomes a vote for the sister pair in that trio.

Votes are aggregated to construct a global cladogram-like summary.

Because Z can be distorted by missingness / tiny SE, a minimum-sites filter is strongly recommended if your script version supports it.

Key interpretation notes (please read)
1) Z “significance” here is a heuristic, not a gold standard

ANGSD’s block jackknife Z values are often treated like a test statistic, but in practice:

Z can become huge if SE becomes tiny (often due to missingness / uneven informative site counts)

Many comparisons are not independent

Many “replicates” are effectively the same genomic signal re-counted across different triplets

So treat |Z| > threshold as a screening / diagnostic rule, not “proof”.

2) If the tool reports “bootstrap/resampling” support, it is triplet resampling, not classic phylogenetic bootstrap

If the script includes resampling for the tree, it is:

resampling triplets / votes with replacement to check stability of the vote aggregation,

not a standard phylogenetic bootstrap of alignment sites.

If you describe this in a paper/report, call it “triplet-vote resampling support” to avoid confusion.

3) Branch lengths are not meaningful

Any tree plot produced here should be treated as a cladogram:

topology reflects aggregated triplet vote tendencies

branch lengths are arbitrary / algorithm-dependent

do not interpret divergence times, drift lengths, etc.

4) Any derived “distance” is not a genetic distance

If you derive a distance from vote proportions or Z magnitudes, it is:

not additive

not a measure of genetic divergence

sensitive to missingness and uneven site counts

Prefer the vote table and sister-support summaries over interpreting distances biologically.

Inputs
1) ANGSD block-jackknifed D-stat file (-i)

This must be the text output from ANGSD’s D-statistics with a header row including at least:

H1, H2, H3

SE

Z

and one of: jackEst, Dstat, or D

If your file uses different column names, rename columns upstream or adjust the script.

2) Sample-to-population map (-m)

A 2-column whitespace-delimited file:

SampleID  Population

Lines beginning with # are ignored.

Example:

# Sample to population map
Ind01  Saitama
Ind02  Saitama
Ind03  Hokuriku
Ind04  Shikoku
Installation / Requirements

This is a bash script that uses standard Unix tools plus R.

Required command-line tools

bash (uses set -euo pipefail)

awk

sort, coreutils

Required R + packages

The script uses R packages such as:

dplyr

tidyr

ggplot2

scales

viridisLite

ape (only if tree mode is enabled)

The script attempts to install missing packages automatically from CRAN (https://cloud.r-project.org). On clusters without internet access, install packages in advance (module/conda/renv/shared library).

Usage
./dstats_poptests_report.sh -i dstats_jackknifed.txt -m sample_to_pop.txt -p results/outprefix [options]
Required arguments

-i ANGSD D-stat file (block-jackknifed)

-m sample → population map

-p output prefix

Options (typical)

-z <float> Z threshold for “significance” summaries (default: 3)

-R also output “real” distinct-population triple comparisons

-T build triplet-vote cladogram summary (requires -R)

(optional) a min-sites / min-informative filter if supported by your script version

Recommended workflow
1) Run topology summaries and plots
./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p out/topo \
  -z 3
2) Also extract “real” comparisons (distinct populations)
./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p out/real \
  -R \
  -z 3
3) Build the triplet-vote cladogram
./dstats_poptests_report.sh \
  -i dstats_jackknifed.txt \
  -m sample_to_pop.txt \
  -p out/tree \
  -R -T \
  -z 3

If your results show extreme Z values (e.g. 50–200+) and you have a min-sites/min-informative-sites option, use it.

Outputs

All outputs are prefixed by -p <prefix>.

Topology test outputs (always)

<prefix>.topology.all.tsv
All extracted rows across all population pairs and types (Right/Wrong/Within).

<prefix>.topology.summary.tsv
Summary per (DoubPop, SingPop, Type) including:

n

prop_sig (using |Z| > zthr)

mean_Z, median_Z

mean_D, median_D

<prefix>.topology.wilcox.tsv
Wilcoxon test of Z ~ Type comparing Right vs Wrong per pair (BH adjusted). Diagnostic only.

<prefix>.topology.Z_box_ALLpairs.pdf
Boxplot of Z distributions by pair and Type.

<prefix>.topology.heatmap_meanZ.pdf
Multi-page PDF: heatmaps of mean Z for Right / Wrong / Within (separate page per Type).

<prefix>.topology.heatmap_propSig.pdf
Faceted heatmap of proportion significant by Type (fixed 0–1 scale).

“Real” outputs (only with -R)

<prefix>.real.tsv
Distinct-population triples, canonicalised H1/H2 with sign flips applied.

<prefix>.real.summary.tsv
Summary by (H1pop, H2pop, H3pop) including:

n

prop_sig_positive / prop_sig_negative / prop_sig_total (using zthr)

mean/median Z and D

<prefix>.real.heatmap_meanZ.pdf
Multi-page heatmap PDF: one page per H3pop (top H3 by number of comparisons).

Tree / cladogram outputs (only with -T, requires -R)

Names vary slightly by script version, but typically:

<prefix>.quartet_votes.tsv
One row per trio (A,B,C) showing which sister pair won and the score margin.

<prefix>.quartet_support.tsv (or similar)
Pairwise support: how often each pair wins as sisters (normalised by opportunities).

<prefix>.quartetNJ.nwk
Newick topology (treat as cladogram; branch lengths not meaningful).

<prefix>.quartetNJ.pdf
Quick tree plot for visualisation.

If resampling support is implemented:

<prefix>.quartetNJ.resample_support.tsv
Triplet-vote resampling support for internal nodes.

<prefix>.quartetNJ.resample_trees.nwk
Resampled trees.

How to interpret the “who tends to be sister?” cladogram

For each trio (A,B,C), the tool picks a “winner” among the three possible sister resolutions. That is a within-triplet decision.

The global tree is then a summary of those within-triplet winners across many triplets.

So yes: it is fundamentally within triplets, aggregated across the dataset.

Common issues & troubleshooting
“could not find function %>%”

This means the R packages providing %>% (usually dplyr/magrittr) were not loaded or not installed.

Fix:

ensure library(dplyr) runs (and succeeds)

install packages ahead of time if your environment blocks CRAN installs

Wilcoxon warning: “cannot compute exact p-value with ties”

This is common and usually harmless. It means many identical Z values exist and R uses an approximation. It does not necessarily indicate a biological problem with that population.

Extremely high |Z| (e.g. 50–200+)

Often indicates tiny SE due to:

small number of informative sites

uneven missingness

a small subset of blocks dominating

Mitigations:

apply a min-sites / min-(ABBA+BABA) filter if available

treat Z as diagnostic not truth

inspect site-count distributions per population/triple

Tree looks very different from PCA

Not necessarily an error:

PCA summarises global covariance structure

D-stat triplet votes summarise specific asymmetry signals (admixture-like patterns)

missingness can distort which triplets are “decidable”

Checks:

use min-sites filtering

restrict to well-covered populations

inspect the votes table to see which triplets dominate

Best practices

Treat Z thresholds as screening, not proof.

Always sanity-check site/informative counts (ABBA+BABA, or equivalent).

Prefer vote/support summaries over branch lengths.

When you have few individuals per population, avoid overstating resampling support.

Citation / attribution

If you use this tool in a paper or report, cite:

ANGSD (for D-statistics computation)

the ABBA-BABA / Patterson’s D-statistics methodology

This script is a wrapper to organise and visualise ANGSD outputs; it does not implement D-statistics from scratch.
