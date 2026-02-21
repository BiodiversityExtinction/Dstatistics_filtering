#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  dstats_poptests_report.sh -i dstats_jackknifed.txt -m sample_to_pop.txt -p out_prefix [options]

Required:
  -i  ANGSD block-jackknifed D-stats file (must have header with H1 H2 H3, SE, Z and Dstat/jackEst/D)
  -m  Sample-to-population map (2 columns: SampleID  Population)
  -p  Output prefix

Optional:
  -z  Z threshold used for "significance" summaries (default: 3)
  -R  Also output "real" D-stat comparisons (distinct populations, canonicalised H1/H2)
  -T  Build a quartet sister-vote tree from real comparisons (requires -R)
      Tree is a cladogram based on "who tends to be sister?" votes (not distances; branch lengths removed).
  -B  Bootstrap replicates for quartet sister-vote tree uncertainty (requires -T). Default: 0
      This is trio-resampling (resampling trios with replacement), not a classical sequence bootstrap.
      Example: -B 200
  -s  Minimum median nSites per trio-resolution to include a trio in voting (default: 0 = no filter)
      nSites is computed as nABBA + nBABA when those columns are present (otherwise NA).
EOF
  exit 1
}

# ------------------------
# Parse arguments
# ------------------------
infile=""
mapfile=""
prefix=""
zthr="3"
do_real="0"
do_tree="0"
boot_reps="0"
min_sites="0"

while getopts "i:m:p:z:s:RTB:" opt; do
  case "$opt" in
    i) infile="$OPTARG" ;;
    m) mapfile="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    z) zthr="$OPTARG" ;;
    s) min_sites="$OPTARG" ;;
    R) do_real="1" ;;
    T) do_tree="1" ;;
    B) boot_reps="$OPTARG" ;;
    *) usage ;;
  esac
done

[[ -z "${infile:-}" || -z "${mapfile:-}" || -z "${prefix:-}" ]] && usage
[[ ! -f "$infile" ]] && { echo "ERROR: infile not found: $infile" >&2; exit 2; }
[[ ! -f "$mapfile" ]] && { echo "ERROR: mapfile not found: $mapfile" >&2; exit 2; }

if [[ "$do_tree" == "1" && "$do_real" != "1" ]]; then
  echo "ERROR: -T requires -R (tree is built from real comparisons)." >&2
  exit 2
fi

if ! [[ "$boot_reps" =~ ^[0-9]+$ ]]; then
  echo "ERROR: -B must be a non-negative integer." >&2
  exit 2
fi
if [[ "$boot_reps" -gt 0 && "$do_tree" != "1" ]]; then
  echo "ERROR: -B requires -T." >&2
  exit 2
fi

# allow min_sites numeric (integer-ish)
if ! [[ "$min_sites" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "ERROR: -s must be numeric (>= 0)." >&2
  exit 2
fi

# ------------------------
# Output files
# ------------------------
allout="${prefix}.topology.all.tsv"
sumout="${prefix}.topology.summary.tsv"
wilout="${prefix}.topology.wilcox.tsv"

plot_zpairs="${prefix}.topology.Z_box_ALLpairs.pdf"
plot_heat_meanZ="${prefix}.topology.heatmap_meanZ.pdf"   # multi-page: Right/Wrong/Within (per-page scaling)
plot_heat_prop="${prefix}.topology.heatmap_propSig.pdf"  # faceted; fixed 0..1 scale

realout="${prefix}.real.tsv"
realsum="${prefix}.real.summary.tsv"
realplot="${prefix}.real.heatmap_meanZ.pdf"              # multi-page: one page per H3 (per-page scaling)

tree_prefix="${prefix}.quartetSisterVote"
tree_votes="${tree_prefix}.quartet_votes.tsv"
tree_newick="${tree_prefix}.nwk"
tree_pdf="${tree_prefix}.pdf"
tree_boot_support="${tree_prefix}.bootstrap_support.tsv" # node support per internal node (if -B > 0)
tree_boot_trees="${tree_prefix}.bootstrap_trees.nwk"     # all bootstrap trees (if -B > 0)

# ------------------------
# Extract populations
# ------------------------
map_pops=$(awk 'NF>=2 && $0 !~ /^[[:space:]]*#/ {print $2}' "$mapfile" | sort -u)
[[ -z "${map_pops:-}" ]] && { echo "ERROR: No populations found in map file." >&2; exit 3; }

# ------------------------
# Write headers
# ------------------------
printf "Type\tTopology\tDoubPop\tSingPop\tH1\tH2\tH3\tPop1\tPop2\tPop3\tD\tSE\tZ\tnSites\n" > "$allout"
if [[ "$do_real" == "1" ]]; then
  printf "H1pop\tH2pop\tH3pop\tH1\tH2\tH3\tD\tSE\tZ\tnSites\n" > "$realout"
fi

# ------------------------
# AWK core: mapping-based, header-aware
# - mode=pair:   topology Right/Wrong for (doub,sing)
# - mode=within: within-pop
# - mode=real:   distinct-pop triples; canonicalise H1pop/H2pop and flip D,Z if swapped
# Also captures nSites if nABBA & nBABA present (nSites = nABBA + nBABA), else NA.
# ------------------------
run_awk='
BEGIN { OFS="\t" }

FNR == NR {
  if ($0 ~ /^[[:space:]]*$/ || $0 ~ /^[[:space:]]*#/) next
  if (NF < 2) next
  samp2pop[$1] = $2
  next
}

FNR == 1 {
  for (i=1; i<=NF; i++) idx[$i]=i
  h1=idx["H1"]; h2=idx["H2"]; h3=idx["H3"]
  se=idx["SE"]; z=idx["Z"]
  d=(idx["jackEst"]?idx["jackEst"]:(idx["Dstat"]?idx["Dstat"]:idx["D"]))

  nabba = idx["nABBA"]
  nbaba = idx["nBABA"]

  if (!h1||!h2||!h3||!se||!z||!d) {
    print "ERROR: Missing required columns in header (need H1 H2 H3 SE Z and Dstat/jackEst/D)." > "/dev/stderr"; exit 10
  }
  next
}

function pop(s){return (s in samp2pop ? samp2pop[s] : "NA")}

{
  H1=$(h1); H2=$(h2); H3=$(h3)
  P1=pop(H1); P2=pop(H2); P3=pop(H3)
  if (P1=="NA"||P2=="NA"||P3=="NA") {
    print "ERROR: Missing population mapping:",H1,H2,H3 > "/dev/stderr"; exit 11
  }
  Dv=$(d)+0; SEv=$(se)+0; Zv=$(z)+0

  nSites="NA"
  if (nabba && nbaba) {
    # treat as numeric if possible; otherwise NA
    a=$(nabba); b=$(nbaba)
    if (a ~ /^-?[0-9.]+$/ && b ~ /^-?[0-9.]+$/) nSites = (a+0) + (b+0)
  }

  if (mode=="pair") {
    if (P1==doub && P2==doub && P3==sing)
      print "Right",doub"_"doub"_"sing,doub,sing,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv,nSites
    else if (P1==sing && P2==doub && P3==doub)
      print "Wrong",sing"_"doub"_"doub,doub,sing,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv,nSites
    else if (P2==sing && P1==doub && P3==doub)
      print "Wrong",sing"_"doub"_"doub,doub,sing,H2,H1,H3,P2,P1,P3,-Dv,SEv,-Zv,nSites
  }
  else if (mode=="within") {
    if (P1==doub && P2==doub && P3==doub)
      print "Within",doub"_"doub"_"doub,doub,doub,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv,nSites
  }
  else if (mode=="real") {
    if (P1!=P2 && P1!=P3 && P2!=P3) {
      # canonicalise H1pop/H2pop; flip D and Z if swapped
      if (P1 < P2)
        print P1,P2,P3,H1,H2,H3,Dv,SEv,Zv,nSites
      else
        print P2,P1,P3,H2,H1,H3,-Dv,SEv,-Zv,nSites
    }
  }
}
'

# ------------------------
# Run topology tests (all ordered pairs) + within
# ------------------------
while read -r doub; do
  while read -r sing; do
    [[ "$doub" == "$sing" ]] && continue
    awk -v doub="$doub" -v sing="$sing" -v mode="pair" "$run_awk" \
      "$mapfile" "$infile" >> "$allout"
  done <<< "$map_pops"
done <<< "$map_pops"

while read -r pop; do
  awk -v doub="$pop" -v mode="within" "$run_awk" \
    "$mapfile" "$infile" >> "$allout"
done <<< "$map_pops"

# ------------------------
# Real comparisons (optional)
# ------------------------
if [[ "$do_real" == "1" ]]; then
  awk -v mode="real" "$run_awk" "$mapfile" "$infile" >> "$realout"
fi

# ------------------------
# R: topology summaries + plots
# ------------------------
Rscript --vanilla - "$allout" "$sumout" "$wilout" \
  "$plot_zpairs" "$plot_heat_meanZ" "$plot_heat_prop" \
  "$zthr" <<'RSCRIPT'

required_pkgs <- c("dplyr","tidyr","ggplot2","scales","viridisLite")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing missing R packages: ", paste(missing, collapse=", "))
  install.packages(missing, repos="https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
sumout <- args[2]
wilout <- args[3]
plot_zpairs <- args[4]
plot_heat_meanZ <- args[5]
plot_heat_prop <- args[6]
zthr <- as.numeric(args[7])

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)

df <- df |>
  dplyr::mutate(
    Pair = paste0(DoubPop, "_vs_", SingPop),
    Sig = abs(Z) > zthr
  )

# Summary (include nSites, SE, |Z| diagnostics)
sumtab <- df |>
  dplyr::group_by(DoubPop, SingPop, Pair, Type) |>
  dplyr::summarise(
    n = dplyr::n(),
    prop_sig = mean(Sig, na.rm=TRUE),
    mean_Z = mean(Z, na.rm=TRUE),
    median_Z = median(Z, na.rm=TRUE),
    mean_D = mean(D, na.rm=TRUE),
    median_D = median(D, na.rm=TRUE),
    med_SE = median(SE, na.rm=TRUE),
    min_SE = suppressWarnings(min(SE, na.rm=TRUE)),
    med_absZ = median(abs(Z), na.rm=TRUE),
    max_absZ = suppressWarnings(max(abs(Z), na.rm=TRUE)),
    med_nSites = if (all(is.na(nSites))) NA_real_ else median(nSites, na.rm=TRUE),
    min_nSites = if (all(is.na(nSites))) NA_real_ else suppressWarnings(min(nSites, na.rm=TRUE)),
    .groups = "drop"
  )

write.table(sumtab, sumout, sep="\t", row.names=FALSE, quote=FALSE)

# Wilcoxon Right vs Wrong per Pair (BH correction)
# Use exact=FALSE to avoid "cannot compute exact p-value with ties" warnings.
wil <- df |>
  dplyr::filter(Type %in% c("Right","Wrong")) |>
  dplyr::group_by(Pair, DoubPop, SingPop) |>
  dplyr::summarise(
    n_right = sum(Type=="Right"),
    n_wrong = sum(Type=="Wrong"),
    p_value = if (n_right >= 3 && n_wrong >= 3) {
      suppressWarnings(stats::wilcox.test(Z ~ Type, exact=FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) |>
  dplyr::mutate(p_adj_BH = p.adjust(p_value, "BH"))

write.table(wil, wilout, sep="\t", row.names=FALSE, quote=FALSE)

# Z plot for ALL pairs, includes Within. Cap y-axis to 99th percentile.
zcap <- stats::quantile(abs(df$Z), probs=0.99, na.rm=TRUE)
zcap <- max(zcap, zthr*2)
npairs <- length(unique(df$Pair))
plot_width <- max(11, min(90, 0.30 * npairs))

pZ <- ggplot2::ggplot(df, ggplot2::aes(x=Pair, y=Z, fill=Type)) +
  ggplot2::geom_boxplot(outlier.shape=NA, alpha=0.85) +
  ggplot2::geom_hline(yintercept=c(-zthr, zthr), linetype="dashed") +
  ggplot2::coord_cartesian(ylim=c(-zcap, zcap)) +
  ggplot2::theme_bw(base_size=11) +
  ggplot2::theme(
    legend.position="top",
    axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1),
    panel.grid.major=ggplot2::element_blank(),
    panel.grid.minor=ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title="Topology test: Z distributions for ALL population pairs",
    subtitle=paste0("Pair = DoubPop_vs_SingPop (DoubPop=H1/H2; SingPop=H3). Y capped at 99th percentile: ±", signif(zcap,4)),
    x="", y="Z"
  )
ggplot2::ggsave(plot_zpairs, pZ, width=plot_width, height=8, limitsize=FALSE)

# ---- Mean Z heatmaps: ONE PDF, one page per Type, per-page scaling ----
heat <- sumtab |>
  dplyr::mutate(Type = factor(Type, levels=c("Right","Wrong","Within")))

grDevices::pdf(plot_heat_meanZ, width=11, height=7)
for (tt in levels(heat$Type)) {
  sub <- heat |> dplyr::filter(Type == tt)
  if (nrow(sub) == 0) next

  L <- max(abs(sub$mean_Z), na.rm=TRUE)
  L <- max(L, 1e-6)

  pMean <- ggplot2::ggplot(sub, ggplot2::aes(x=SingPop, y=DoubPop, fill=mean_Z)) +
    ggplot2::geom_tile(color="white", linewidth=0.2) +
    ggplot2::scale_fill_gradient2(
      low="#3B4CC0", mid="white", high="#B40426", midpoint=0,
      limits=c(-L, L), oob=scales::squish
    ) +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(
      panel.grid=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=45, hjust=1)
    ) +
    ggplot2::labs(
      title=paste("Topology test mean Z:", tt),
      subtitle=paste0("Colour scale symmetric per page; limits ±", signif(L,3)),
      x="SingPop (H3)", y="DoubPop (H1/H2)", fill="Mean Z"
    )
  print(pMean)
}
grDevices::dev.off()

# propSig heatmap: fixed 0..1 scale is comparable across Types
pProp <- ggplot2::ggplot(heat, ggplot2::aes(x=SingPop, y=DoubPop, fill=prop_sig)) +
  ggplot2::geom_tile(color="white", linewidth=0.2) +
  ggplot2::scale_fill_viridis_c(option="C", limits=c(0,1)) +
  ggplot2::facet_wrap(~Type) +
  ggplot2::theme_bw(base_size=11) +
  ggplot2::theme(
    panel.grid=ggplot2::element_blank(),
    axis.text.x=ggplot2::element_text(angle=45, hjust=1)
  ) +
  ggplot2::labs(
    title=paste0("Topology test: proportion with |Z| > ", zthr),
    x="SingPop (H3)", y="DoubPop (H1/H2)", fill="Prop. sig"
  )
ggplot2::ggsave(plot_heat_prop, pProp, width=11, height=7)

RSCRIPT

# ------------------------
# R: real comparisons summary + plots (optional)
# ------------------------
if [[ "$do_real" == "1" ]]; then
Rscript --vanilla - "$realout" "$realsum" "$realplot" "$zthr" <<'RSCRIPT2'

required_pkgs <- c("dplyr","ggplot2","scales","viridisLite")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing missing R packages: ", paste(missing, collapse=", "))
  install.packages(missing, repos="https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outsum <- args[2]
outplot <- args[3]
zthr <- as.numeric(args[4])

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)

sumtab <- df |>
  dplyr::group_by(H1pop, H2pop, H3pop) |>
  dplyr::summarise(
    n = dplyr::n(),
    prop_sig_positive = mean(Z >  zthr, na.rm=TRUE),
    prop_sig_negative = mean(Z < -zthr, na.rm=TRUE),
    prop_sig_total    = mean(abs(Z) > zthr, na.rm=TRUE),
    mean_Z = mean(Z, na.rm=TRUE),
    median_Z = median(Z, na.rm=TRUE),
    mean_D = mean(D, na.rm=TRUE),
    median_D = median(D, na.rm=TRUE),
    med_SE = median(SE, na.rm=TRUE),
    min_SE = suppressWarnings(min(SE, na.rm=TRUE)),
    med_absZ = median(abs(Z), na.rm=TRUE),
    max_absZ = suppressWarnings(max(abs(Z), na.rm=TRUE)),
    med_nSites = if (all(is.na(nSites))) NA_real_ else median(nSites, na.rm=TRUE),
    min_nSites = if (all(is.na(nSites))) NA_real_ else suppressWarnings(min(nSites, na.rm=TRUE)),
    .groups="drop"
  )

write.table(sumtab, outsum, sep="\t", row.names=FALSE, quote=FALSE)

# Plot: one page per H3, per-page scaling (use top H3 by total n)
topH3 <- sumtab |>
  dplyr::count(H3pop, wt=n, name="total_n") |>
  dplyr::arrange(dplyr::desc(total_n)) |>
  dplyr::slice_head(n=12) |>
  dplyr::pull(H3pop)

plotdat <- sumtab |>
  dplyr::filter(H3pop %in% topH3)

grDevices::pdf(outplot, width=11, height=7)
for (h3 in unique(plotdat$H3pop)) {
  sub <- plotdat |> dplyr::filter(H3pop == h3)
  if (nrow(sub) == 0) next

  L <- max(abs(sub$mean_Z), na.rm=TRUE)
  L <- max(L, 1e-6)

  p <- ggplot2::ggplot(sub, ggplot2::aes(x=H2pop, y=H1pop, fill=mean_Z)) +
    ggplot2::geom_tile(color="white", linewidth=0.2) +
    ggplot2::scale_fill_gradient2(
      low="#3B4CC0", mid="white", high="#B40426", midpoint=0,
      limits=c(-L, L), oob=scales::squish
    ) +
    ggplot2::theme_bw(base_size=11) +
    ggplot2::theme(
      panel.grid=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=45, hjust=1)
    ) +
    ggplot2::labs(
      title="Real D-stat mean Z (canonical H1/H2)",
      subtitle=paste0("H3 = ", h3, " | Colour scale symmetric per page; limits ±", signif(L,3)),
      x="H2 population (canonical)", y="H1 population (canonical)", fill="Mean Z"
    )
  print(p)
}
grDevices::dev.off()

RSCRIPT2
fi

# ------------------------
# R: quartet sister-vote tree (NOT distance-based)
# - Uses real.summary.tsv (mean_Z + med_nSites) to score trio resolutions and determine sister pairs
# - For each trio (A,B,C), choose the resolution with smallest |mean_Z| as the "best sister" (A-B, A-C, or B-C)
# - Tree is built by greedy agglomeration on sister-vote *counts* between clusters (pure voting, no distances)
# - Bootstrap: resample trios with replacement, rebuild tree; node support via prop.clades (trees rooted at an arbitrary tip for consistency)
# ------------------------
if [[ "$do_tree" == "1" ]]; then
Rscript --vanilla - "$realsum" "$tree_votes" "$tree_newick" "$tree_pdf" \
  "$tree_boot_support" "$tree_boot_trees" "$boot_reps" "$min_sites" <<'RSCRIPT3'

required_pkgs <- c("dplyr","tidyr","ape")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
  message("Installing missing R packages: ", paste(missing, collapse=", "))
  install.packages(missing, repos="https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only=TRUE))

args <- commandArgs(trailingOnly=TRUE)
realsum <- args[1]
out_votes <- args[2]
out_nwk <- args[3]
out_pdf <- args[4]
out_boot_support <- args[5]
out_boot_trees <- args[6]
B <- as.integer(args[7])
min_sites <- as.numeric(args[8])

tab <- read.delim(realsum, sep="\t", header=TRUE, stringsAsFactors=FALSE)

need_cols <- c("H1pop","H2pop","H3pop","mean_Z")
if (!all(need_cols %in% names(tab))) {
  stop("real.summary.tsv must contain columns: H1pop, H2pop, H3pop, mean_Z")
}

has_sites <- "med_nSites" %in% names(tab)

# score: smaller |mean_Z| => more tree-like for that trio resolution
tab <- tab |> dplyr::mutate(score = abs(mean_Z))

# key by canonical (H1pop < H2pop) already; still trust input
key <- paste(tab$H1pop, tab$H2pop, tab$H3pop, sep="|")
score_map <- setNames(tab$score, key)
sites_map <- if (has_sites) setNames(tab$med_nSites, key) else NULL

get_key <- function(A, B, C) {
  x1 <- ifelse(A < B, A, B)
  x2 <- ifelse(A < B, B, A)
  paste(x1, x2, C, sep="|")
}

get_score <- function(A, B, C) {
  k <- get_key(A, B, C)
  v <- score_map[[k]]
  if (is.null(v)) return(NA_real_)
  as.numeric(v)
}

get_sites <- function(A, B, C) {
  if (!has_sites) return(NA_real_)
  k <- get_key(A, B, C)
  v <- sites_map[[k]]
  if (is.null(v)) return(NA_real_)
  as.numeric(v)
}

pops <- sort(unique(c(tab$H1pop, tab$H2pop, tab$H3pop)))
n <- length(pops)
if (n < 4) stop("Need >= 4 populations to build a quartet-based tree.")

# Enumerate all complete trios, pick best sister per trio
trios <- list()
votes <- list()
ti <- 1L
vi <- 1L

for (i in 1:(n-2)) for (j in (i+1):(n-1)) for (k in (j+1):n) {
  A <- pops[i]; Bp <- pops[j]; C <- pops[k]

  sAB_C <- get_score(A,  Bp, C)   # ((A,B),C)
  sAC_B <- get_score(A,  C,  Bp)  # ((A,C),B)
  sBC_A <- get_score(Bp, C,  A)   # ((B,C),A)

  if (any(is.na(c(sAB_C, sAC_B, sBC_A)))) next

  # Optional sites filter: require each of the 3 resolutions has med_nSites >= min_sites
  if (min_sites > 0 && has_sites) {
    nsAB_C <- get_sites(A,  Bp, C)
    nsAC_B <- get_sites(A,  C,  Bp)
    nsBC_A <- get_sites(Bp, C,  A)

    if (any(is.na(c(nsAB_C, nsAC_B, nsBC_A)))) next
    if (any(c(nsAB_C, nsAC_B, nsBC_A) < min_sites)) next
  }

  scores <- c(sAB_C, sAC_B, sBC_A)
  w <- which.min(scores)

  ord <- sort(scores)
  winner <- ord[1]
  runnerup <- ord[2]
  margin <- runnerup - winner

  if (w == 1) sister <- c(A, Bp)
  if (w == 2) sister <- c(A, C)
  if (w == 3) sister <- c(Bp, C)

  trios[[ti]] <- c(A, Bp, C); ti <- ti + 1L

  votes[[vi]] <- data.frame(
    A = A, B = Bp, C = C,
    best_sister1 = sister[1],
    best_sister2 = sister[2],
    score_AB_C = sAB_C,
    score_AC_B = sAC_B,
    score_BC_A = sBC_A,
    winner = winner,
    runnerup = runnerup,
    margin = margin,
    stringsAsFactors = FALSE
  )
  vi <- vi + 1L
}

if (length(trios) < 1) stop("No complete trios found (need all three resolutions per trio; and pass min_sites filter if set).")

votes_df <- dplyr::bind_rows(votes)
write.table(votes_df, out_votes, sep="\t", quote=FALSE, row.names=FALSE)

# Build sister vote count matrix from a vector of trio indices
build_pair_votes <- function(idx) {
  wins <- matrix(0, n, n, dimnames=list(pops, pops))
  for (t in idx) {
    A <- trios[[t]][1]; Bp <- trios[[t]][2]; C <- trios[[t]][3]

    sAB_C <- get_score(A,  Bp, C)
    sAC_B <- get_score(A,  C,  Bp)
    sBC_A <- get_score(Bp, C,  A)
    scores <- c(sAB_C, sAC_B, sBC_A)
    w <- which.min(scores)
    if (w == 1) sister <- c(A, Bp)
    if (w == 2) sister <- c(A, C)
    if (w == 3) sister <- c(Bp, C)

    wins[sister[1], sister[2]] <- wins[sister[1], sister[2]] + 1
    wins[sister[2], sister[1]] <- wins[sister[2], sister[1]] + 1
  }
  wins
}

# Greedy agglomeration into a cladogram using ONLY sister-vote counts.
# Each step merges the two clusters with maximum total wins between them.
greedy_tree_from_votes <- function(wins) {
  # clusters: named list. Each cluster has members + newick label
  clusters <- lapply(pops, function(x) list(members=x, label=x))
  names(clusters) <- pops

  cluster_score <- function(ci, cj) {
    mi <- clusters[[ci]]$members
    mj <- clusters[[cj]]$members
    # total wins between members, symmetric, so sum upper only (or full /2)
    sum(wins[mi, mj, drop=FALSE], na.rm=TRUE)
  }

  while (length(clusters) > 1) {
    cn <- names(clusters)
    best <- NULL
    best_score <- -Inf

    for (i in 1:(length(cn)-1)) for (j in (i+1):length(cn)) {
      sc <- cluster_score(cn[i], cn[j])
      if (sc > best_score) {
        best_score <- sc
        best <- c(cn[i], cn[j])
      }
    }

    # if everything is 0, still merge deterministically by alphabetical order
    if (is.infinite(best_score) || is.na(best_score)) best <- cn[1:2]

    a <- best[1]; b <- best[2]
    new_members <- c(clusters[[a]]$members, clusters[[b]]$members)
    new_label <- paste0("(", clusters[[a]]$label, ",", clusters[[b]]$label, ")")

    # new cluster name should not collide
    new_name <- paste(sort(new_members), collapse="|")

    # delete old clusters safely
    clusters[[a]] <- NULL
    clusters[[b]] <- NULL

    clusters[[new_name]] <- list(members=new_members, label=new_label)
  }

  # return final Newick string (rooted); we will unroot for plotting and strip branch lengths
  paste0(clusters[[1]]$label, ";")
}

all_idx <- seq_along(trios)
wins_all <- build_pair_votes(all_idx)
nwk <- greedy_tree_from_votes(wins_all)

# write Newick
ape::write.tree(ape::read.tree(text=nwk), file=out_nwk)

# Bootstrap: trio-resampling (not classical sequence bootstrap)
boot_vals <- NULL
boot_trees <- NULL

# Rooting for consistent clade comparison in prop.clades
root_tip <- pops[1]

if (B > 0) {
  set.seed(1)
  boot_trees <- vector("list", B)

  for (b in 1:B) {
    idx_b <- sample(all_idx, length(all_idx), replace=TRUE)
    wins_b <- build_pair_votes(idx_b)
    nwk_b <- greedy_tree_from_votes(wins_b)
    tb <- ape::read.tree(text=nwk_b)

    # root consistently for prop.clades
    tb <- ape::root(tb, outgroup=root_tip, resolve.root=TRUE)

    # remove lengths (cladogram)
    tb$edge.length <- NULL

    boot_trees[[b]] <- tb
  }

  # Write bootstrap trees
  ape::write.tree(boot_trees, file=out_boot_trees)
}

# Final tree object
final_tree <- ape::read.tree(text=nwk)
final_tree <- ape::root(final_tree, outgroup=root_tip, resolve.root=TRUE)
final_tree$edge.length <- NULL

# Bootstrap node supports on rooted trees for prop.clades
if (B > 0) {
  counts <- ape::prop.clades(final_tree, boot_trees)
  boot_vals <- round(100 * counts / B, 1)

  outtab <- data.frame(
    node = which(!is.na(counts)),
    bootstrap = boot_vals,
    stringsAsFactors = FALSE
  )
  write.table(outtab, out_boot_support, sep="\t", quote=FALSE, row.names=FALSE)
}

# Plot unrooted cladogram (display only; Newick is rooted by construction)
grDevices::pdf(out_pdf, width=10, height=8)
plot(ape::unroot(final_tree), type="unrooted",
     main="Quartet sister-vote cladogram from D-statistics (winner = smallest |mean Z| per trio)")
mtext("Tree is built from 'who tends to be sister?' votes (NOT distances; branch lengths removed).", side=3, line=0.2, cex=0.85)
mtext("Bootstrap (if enabled) = trio-resampling with replacement, not classical sequence bootstrapping.", side=3, line=1.2, cex=0.75)

if (!is.null(boot_vals)) {
  # nodelabels expects node numbers in the *rooted* tree; place labels on rooted, but plot is unrooted.
  # So: add labels before unrooting is hard to align; simplest is to plot rooted with "fan" style.
  # Instead: show rooted cladogram with bootstrap labels on a second page for clarity.
  plot(final_tree, main="Same tree (rooted for node labels only)")
  nodelabels(text=boot_vals, frame="none", adj=c(1.1, -0.2), cex=0.8)
  legend("topleft", legend=paste0("Bootstrap: ", B, " trio-resampling replicates"), bty="n")
}
grDevices::dev.off()

RSCRIPT3
fi

echo "Done."
echo "Wrote:"
echo "  $allout"
echo "  $sumout"
echo "  $wilout"
echo "  $plot_zpairs"
echo "  $plot_heat_meanZ   (multi-page: Right/Wrong/Within)"
echo "  $plot_heat_prop"
if [[ "$do_real" == "1" ]]; then
  echo "  $realout"
  echo "  $realsum  (includes med/min nSites if available; plus SE and |Z| diagnostics)"
  echo "  $realplot (multi-page: top H3 by n)"
fi
if [[ "$do_tree" == "1" ]]; then
  echo "  $tree_votes"
  echo "  $tree_newick"
  echo "  $tree_pdf"
  if [[ "$boot_reps" -gt 0 ]]; then
    echo "  $tree_boot_support"
    echo "  $tree_boot_trees"
  fi
fi
