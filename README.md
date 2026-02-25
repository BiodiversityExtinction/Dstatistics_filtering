# DstatWorkbench.R

`DstatWorkbench.R` is an R CLI for filtering and summarizing ANGSD block-jackknifed D-statistics, with optional plotting, topology consistency checks, ranking, and PCA/ordination outputs.

This repository is documented for the R CLI only. Bash variants are not part of the supported interface.

## What it does

From an ANGSD D-stat table and a sample-to-population map, the script can:

- Extract canonical distinct-population triples (`--triples`)
- Summarize triples by population combination
- Build trio-vote support matrices (both Z-based and D-based)
- Run population ordination (`--pcaPop`) from support matrices
- Run individual PCA (`--pcaInd`) from H3-focused affinity profiles
- Plot topology diagnostics (`--topology`)
- Rank likely sister populations for target H3 populations/samples (`--rankH3Pop`, `--rankH3Ind`)
- Filter triple comparisons by a supplied rooted tree (`--plotRealByTree`)
- Plot D-score distributions for selected H3 individuals (`--plotH3Ind`)

## Requirements

R packages:

- `dplyr`
- `tidyr`
- `ggplot2`
- `ape` (only required for `--plotRealByTree`)

The script does not auto-install packages; it fails with a clear message if missing.

## Inputs

### 1. D-stat file (`--input`)

ANGSD block-jackknifed table with header containing:

- `H1`, `H2`, `H3`
- `SE`, `Z`
- one of `jackEst` / `Dstat` / `D`

For site filtering, the script uses:

- `nABBA + nBABA` if both columns exist
- otherwise `nSites` if present
- otherwise `0`

### 2. Population map (`--Popfile`)

Two-column whitespace-delimited file:

- column 1: `SampleID`
- column 2: `PopulationID`

## Usage

```bash
Rscript DstatWorkbench.R \
  --input dstats.jack.txt \
  --Popfile sample_to_pop.txt \
  --out run_prefix \
  [options]
```

## Parameters

### Required

- `--input FILE`: D-stat input table
- `--Popfile FILE`: sample-to-population map
- `--out PREFIX`: output prefix

### Core filters

- `--minsites INT`: minimum informative sites to retain a row (default `0`)
- `-z FLOAT`: significance threshold for `|Z|` labels (default `3`)

### Main analysis switches

- `--triples`: build distinct-population triple tables (`.real.tsv`, `.real.summary.tsv`)
- `--topology`: produce topology pairplot TSV/PDF
- `--votes`: build trio-vote tables and support matrices (requires `--triples`)
- `--pcaPop[=POP_IDS]`: population ordination from support matrix (requires `--triples`)
  - optional comma-separated population subset
- `--pcaInd[=SAMPLE_IDS]`: individual PCA from H3-only affinity (requires `--triples`)
  - optional comma-separated sample subset
- `--rankH3Pop POP`: rank closest populations where `H3pop == POP`
- `--rankH3Ind SAMPLE`: rank closest populations where `H3 == SAMPLE`
- `--plotH3Ind LIST`: comma-separated H3 sample IDs for D-score plotting
- `--plotRealByTree FILE.nwk`: keep only triples consistent with rooted tree topology (requires `--triples`)
  - tree must contain tip `Outgroup`

### PCA controls

- `--pcaK INT`: number of PCs written to TSV (default `6`, minimum `2`)
- `--pcaPlotPC A,B`: PCs plotted in PDF (default `1,2`)

## Z-based vs D-based outputs

For `--votes`, `--pcaPop`, and `--pcaInd`, the script writes both:

- Z-based outputs (default filenames)
- D-based sensitivity outputs with `.d.` in filename

This lets you compare robustness when low site counts may inflate Z.

## Output files

Only files relevant to selected options are written.

### `--triples`

- `<prefix>.real.tsv`
  - Row-level distinct-population comparisons (`H1pop/H2pop/H3pop` all different), with canonicalized `H1pop/H2pop` and sign-adjusted `D`,`Z`.
- `<prefix>.real.summary.tsv`
  - Grouped summaries per population triple: counts, site depth summaries, proportions above `|Z|` threshold, and central tendency for `Z` and `D`.

### `--plotRealByTree`

- `<prefix>.real.byTree.tsv`
  - Subset of triple rows consistent with the supplied rooted tree topology, including significance label (`sig`/`non`).
- `<prefix>.real.byTree.pdf`
  - Jitter plot of `D` for tree-consistent topologies (`x=D`, `y=topology`, point shape indicates significance).

### `--topology`

- `<prefix>.topology.mapped_rows.tsv`
  - All mapped/filter-passing input rows with sample/pop labels and numeric fields (`D`,`SE`,`Z`,`nSites`).
- `<prefix>.topology.pairplot.tsv`
  - Plot-ready table of pairwise topology contrasts (`(Pop1,Pop1),Pop2` vs `(Pop1,Pop2),Pop1`) with `D`.
- `<prefix>.topology.pairplot.pdf`
  - Jitter plot comparing topology contrast `D` distributions across population pairs.

### `--votes` (also produced when `--pcaPop` or `--pcaInd` is requested)

Z-based:

- `<prefix>.votes.tsv`
  - One row per population trio (`A,B,C`) with winning sister pair and trio-resolution scores (from `abs(mean_Z)`).
- `<prefix>.support_matrix.tsv`
  - Symmetric pairwise sister-support matrix (0-1) derived from Z-based trio votes.

D-based:

- `<prefix>.votes.d.tsv`
  - Same as `.votes.tsv`, but trio-resolution scores use `abs(mean_D)`.
- `<prefix>.support_matrix.d.tsv`
  - Same as `.support_matrix.tsv`, but derived from D-based votes.

### `--pcaPop`

Z-based:

- `<prefix>.pca.pop.tsv`
  - Population coordinates from MDS on `1 - support_matrix` (Z-based support).
- `<prefix>.pca.pop.pdf`
  - Population ordination plot from Z-based support matrix.

D-based:

- `<prefix>.pca.pop.d.tsv`
  - Population coordinates from MDS on `1 - support_matrix.d` (D-based support).
- `<prefix>.pca.pop.d.pdf`
  - Population ordination plot from D-based support matrix.

### `--pcaInd`

Z-based:

- `<prefix>.pca.ind.tsv`
  - Individual PC coordinates from H3-only signed affinity profiles built with `Z` (`H1=-Z`, `H2=+Z`).
- `<prefix>.pca.ind.pdf`
  - Individual PCA scatter plot from Z-based affinities.

D-based:

- `<prefix>.pca.ind.d.tsv`
  - Individual PC coordinates from H3-only signed affinity profiles built with `D` (`H1=-D`, `H2=+D`).
- `<prefix>.pca.ind.d.pdf`
  - Individual PCA scatter plot from D-based affinities.

### `--rankH3Pop`

- `<prefix>.rankH3Pop.tsv`
  - Ranked candidate sister populations for target `H3pop`, with weighted wins (`abs(Z)`) and win counts.

### `--rankH3Ind`

- `<prefix>.rankH3Ind.tsv`
  - Ranked candidate sister populations for target H3 sample, with weighted wins (`abs(Z)`) and win counts.

### `--plotH3Ind`

- `<prefix>.plotH3Ind.mapped_rows.tsv`
  - Full mapped/filter-passing rows used as the source table for H3 plotting.
- `<prefix>.plotH3Ind.tsv`
  - Filtered plot table for selected H3 individuals: `D`,`Z`, comparison label (`H1pop_H2pop`), significance label, and metadata.
- `<prefix>.plotH3Ind.pdf`
  - Jitter plot of `D` for selected H3 individuals (`x=D`, `y=H3`, color=`H1pop_H2pop`, shape by significance).

## Example commands

All main analyses:

```bash
Rscript DstatWorkbench.R \
  --input Japan_mod_anc_all_Dstats_specout.jack.txt \
  --Popfile Bamlist_Dstats_wanc.pops \
  --out test \
  --triples --topology --votes \
  --pcaPop --pcaInd --pcaK 3 --pcaPlotPC 1,2
```

Tree-filtered triple plot + selected H3 plotting:

```bash
Rscript DstatWorkbench.R \
  --input Japan_mod_anc_all_Dstats_specout.jack.txt \
  --Popfile Bamlist_Dstats_wanc.pops \
  --out test \
  --triples \
  --plotRealByTree tree.nwk \
  --plotH3Ind Kumamoto_his,Saitama-217_his,Saitama-218_his,Saitama-220_his
```

## Interpretation notes

- `|Z| > threshold` is a diagnostic heuristic based on block-jackknife Z.
- Very low informative-site counts can inflate Z; use `--minsites` and compare `.d.` outputs.
- If Z- and D-based PCA/ordination disagree strongly, treat low-coverage samples with caution.
