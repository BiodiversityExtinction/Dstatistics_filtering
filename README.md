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
  - distinct-population triple rows with canonicalized `H1pop/H2pop`
- `<prefix>.real.summary.tsv`
  - grouped summary per `H1pop,H2pop,H3pop`

### `--plotRealByTree`

- `<prefix>.real.byTree.tsv`
- `<prefix>.real.byTree.pdf`

### `--topology`

- `<prefix>.topology.mapped_rows.tsv`
- `<prefix>.topology.pairplot.tsv`
- `<prefix>.topology.pairplot.pdf`

### `--votes` (also produced when `--pcaPop` or `--pcaInd` is requested)

Z-based:

- `<prefix>.votes.tsv`
- `<prefix>.support_matrix.tsv`

D-based:

- `<prefix>.votes.d.tsv`
- `<prefix>.support_matrix.d.tsv`

### `--pcaPop`

Z-based:

- `<prefix>.pca.pop.tsv`
- `<prefix>.pca.pop.pdf`

D-based:

- `<prefix>.pca.pop.d.tsv`
- `<prefix>.pca.pop.d.pdf`

### `--pcaInd`

Z-based:

- `<prefix>.pca.ind.tsv`
- `<prefix>.pca.ind.pdf`

D-based:

- `<prefix>.pca.ind.d.tsv`
- `<prefix>.pca.ind.d.pdf`

### `--rankH3Pop`

- `<prefix>.rankH3Pop.tsv`

### `--rankH3Ind`

- `<prefix>.rankH3Ind.tsv`

### `--plotH3Ind`

- `<prefix>.plotH3Ind.mapped_rows.tsv`
- `<prefix>.plotH3Ind.tsv`
- `<prefix>.plotH3Ind.pdf`

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
