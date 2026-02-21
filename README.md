
Choose `-s` based on your dataset. If you see extreme Z values (e.g. 50–200+) and warnings that suggest few informative sites, increase `-s`.

---

## Outputs

All outputs are prefixed by `-p <prefix>`.

### Topology test outputs (always produced)
- `<prefix>.topology.all.tsv`  
  All extracted rows across all population pairs and types (Right/Wrong/Within).

- `<prefix>.topology.summary.tsv`  
  Summary per `(DoubPop, SingPop, Type)` including:
  - `n`
  - `prop_sig` (using `|Z| > zthr`)
  - `mean_Z`, `median_Z`
  - `mean_D`, `median_D`
  - (and any added diagnostics such as site counts if enabled)

- `<prefix>.topology.wilcox.tsv`  
  Wilcoxon test of `Z ~ Type` comparing Right vs Wrong per pair (BH adjusted).  
  This is a diagnostic to see if Right and Wrong behave differently, not a definitive test.

- `<prefix>.topology.Z_box_ALLpairs.pdf`  
  Boxplot of Z distributions by pair and Type.

- `<prefix>.topology.heatmap_meanZ.pdf`  
  Multi-page PDF: heatmaps of mean Z for Right / Wrong / Within (separate page per Type).

- `<prefix>.topology.heatmap_propSig.pdf`  
  Faceted heatmap of proportion significant by Type (fixed 0–1 scale).

### “Real” outputs (only with `-R`)
- `<prefix>.real.tsv`  
  Distinct-population triples, canonicalised H1/H2 with sign flips applied.

- `<prefix>.real.summary.tsv`  
  Summary by `(H1pop, H2pop, H3pop)` including:
  - `n`
  - prop_sig_positive / negative / total (using zthr)
  - mean/median Z and D
  - (and optionally site diagnostics)

- `<prefix>.real.heatmap_meanZ.pdf`  
  Multi-page heatmap PDF: one page per H3pop (top H3 by number of comparisons).

### Tree / cladogram outputs (only with `-T`, and requires `-R`)
- `<prefix>.quartet_votes.tsv`  
  One row per trio (A,B,C) showing which sister pair won and the score margin.

- `<prefix>.quartet_support.tsv` (or similar)  
  Pairwise support table: how often each pair wins as sisters (normalised by opportunities).

- `<prefix>.quartetNJ.nwk`  
  Newick topology. Treat as a cladogram (branch lengths not meaningful).

- `<prefix>.quartetNJ.pdf`  
  Tree plot for quick visualisation.

If resampling support is implemented:
- `<prefix>.quartetNJ.resample_support.tsv`  
  “Triplet vote resampling support” for internal nodes.
- `<prefix>.quartetNJ.resample_trees.nwk`  
  Set of resampled trees.

---

## Understanding the “tree” output

### What “who tends to be sister?” means
For each trio of populations (A,B,C):

- compute a score for each resolution (A,B)|C, (A,C)|B, (B,C)|A
- pick the winner (smallest score)
- that is a **within-trio** decision

Then across all trios, aggregate winners to see:

- which pairs most often win as sisters
- and build a global cladogram-like summary

### What it does NOT mean
- It is not a rooted phylogeny
- It is not a coalescent species tree
- It does not produce meaningful branch lengths
- It can be biased by missingness / uneven site counts

---

## Common issues & troubleshooting

### “could not find function %>%”
This means the R packages providing `%>%` (usually **dplyr** or magrittr) were not loaded.

- The script should `library(dplyr)` early.  
- If the cluster blocks package installs, pre-install dependencies:
  - `install.packages("dplyr")` (and the others) in an interactive session
  - or use conda / renv / a shared R library

### “cannot compute exact p-value with ties” (Wilcoxon warning)
This is normal when there are ties in Z values. R falls back to an approximation.

This warning does **not** mean your population is “weird”. It usually means:
- many identical Z values (often from repeated or discrete outcomes), or
- small sample size per group

If Wilcoxon is not central to your analysis, you can ignore these warnings.

### Extremely high |Z| (e.g., 50–200+)
Often indicates tiny SE due to:
- small number of informative sites
- uneven missingness
- a small subset of blocks dominating

Fixes:
- use a minimum-sites filter (`-s`)
- inspect distributions of site counts / ABBA+BABA
- remove or merge problematic populations with too little data
- interpret Z as a diagnostic, not as a truth signal

### Tree looks very different from PCA
This can happen because:
- PCA captures overall covariance structure across many loci
- D-stat triplet votes capture specific **asymmetry signals** (often admixture-like)
- missingness can distort which triplets are “decidable”

Suggested checks:
- verify populations have comparable site counts across comparisons
- increase `-s`
- compare only a subset of well-covered populations
- inspect `<prefix>.quartet_votes.tsv` to see which trios drive the topology

### “replacement has length zero” / Newick manipulation errors
This typically indicates a bug in tree-string editing or that the tree-building algorithm received an empty or malformed set of votes.

Fixes:
- confirm `-R` outputs have enough rows
- confirm at least 4 populations with sufficient data
- check that after filtering there are still complete trios (all three resolutions present)
- raise/lower `-s` depending on whether filtering is too strict

---

## Best practices

- Treat Z-based “significance” as **screening**, not confirmation.
- Always look at data quantity diagnostics (ABBA/BABA or site counts).
- Prefer **vote tables and support summaries** over branch lengths.
- When you have few individuals per population, avoid overstating bootstrap-like supports.
- Keep outputs and script version together for reproducibility.

---

## Citation / attribution

If you use this tool in a paper or report, cite:
- ANGSD (for D-statistics computation)
- the relevant D-statistics methodology (ABBA-BABA / Patterson’s D)

This script itself is a convenience wrapper to organise and visualise ANGSD outputs; it does not implement D-statistics from scratch.

---

## License

Add your license here (e.g., MIT, GPL-3.0, etc.).

---

## Contact / Issues

Open a GitHub issue with:
- your command line used
- the first 30 lines of the input files (with sensitive IDs removed)
- the exact error output
- the generated `<prefix>.*.tsv` summaries if available
