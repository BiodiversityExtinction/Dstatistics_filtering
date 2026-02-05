#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage:
  dstats_poptests_report.sh -i dstats_jackknifed.txt -m sample_to_pop.txt -p out_prefix [options]

Required:
  -i  ANGSD jackknife output (header required; e.g. H1 H2 H3 ... Dstat/jackEst SE Z)
  -m  Sample-to-population mapping file (2 columns: SampleID  Population)
  -p  Output prefix

Options:
  -z  Z threshold for significance (default: 3)
  -R  Also output "real" D-stat comparisons (distinct populations, canonicalised H1/H2)
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

while getopts "i:m:p:z:R" opt; do
  case "$opt" in
    i) infile="$OPTARG" ;;
    m) mapfile="$OPTARG" ;;
    p) prefix="$OPTARG" ;;
    z) zthr="$OPTARG" ;;
    R) do_real="1" ;;
    *) usage ;;
  esac
done

[[ -z "${infile:-}" || -z "${mapfile:-}" || -z "${prefix:-}" ]] && usage
[[ ! -f "$infile" ]] && { echo "ERROR: infile not found: $infile" >&2; exit 2; }
[[ ! -f "$mapfile" ]] && { echo "ERROR: mapfile not found: $mapfile" >&2; exit 2; }

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

# ------------------------
# Extract populations
# ------------------------
map_pops=$(awk 'NF>=2 && $0 !~ /^[[:space:]]*#/ {print $2}' "$mapfile" | sort -u)
[[ -z "${map_pops:-}" ]] && { echo "ERROR: No populations found in map file." >&2; exit 3; }

# ------------------------
# Write headers
# ------------------------
printf "Type\tTopology\tDoubPop\tSingPop\tH1\tH2\tH3\tPop1\tPop2\tPop3\tD\tSE\tZ\n" > "$allout"
if [[ "$do_real" == "1" ]]; then
  printf "H1pop\tH2pop\tH3pop\tH1\tH2\tH3\tD\tSE\tZ\n" > "$realout"
fi

# ------------------------
# AWK core: mapping-based, header-aware
# - mode=pair:   topology Right/Wrong for (doub,sing)
# - mode=within: within-pop
# - mode=real:   distinct-pop triples; canonicalise H1pop/H2pop and flip D,Z if swapped
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
  if (!h1||!h2||!h3||!se||!z||!d) {
    print "ERROR: Missing required columns in header" > "/dev/stderr"; exit 10
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

  if (mode=="pair") {
    if (P1==doub && P2==doub && P3==sing)
      print "Right",doub"_"doub"_"sing,doub,sing,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv
    else if (P1==sing && P2==doub && P3==doub)
      print "Wrong",sing"_"doub"_"doub,doub,sing,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv
    else if (P2==sing && P1==doub && P3==doub)
      print "Wrong",sing"_"doub"_"doub,doub,sing,H2,H1,H3,P2,P1,P3,-Dv,SEv,-Zv
  }
  else if (mode=="within") {
    if (P1==doub && P2==doub && P3==doub)
      print "Within",doub"_"doub"_"doub,doub,doub,H1,H2,H3,P1,P2,P3,Dv,SEv,Zv
  }
  else if (mode=="real") {
    if (P1!=P2 && P1!=P3 && P2!=P3) {
      # canonicalise H1pop/H2pop; flip D and Z if swapped
      if (P1 < P2)
        print P1,P2,P3,H1,H2,H3,Dv,SEv,Zv
      else
        print P2,P1,P3,H2,H1,H3,-Dv,SEv,-Zv
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
# - auto-installs missing packages
# - Z plot uses ALL pairs and includes Within
# - meanZ heatmap is ONE multi-page PDF (Right/Wrong/Within), per-page scaling
# - propSig heatmap fixed 0..1 scale, faceted
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

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>%
  mutate(
    Pair = paste0(DoubPop, "_vs_", SingPop),
    Sig = abs(Z) > zthr
  )

# Summary
sumtab <- df %>%
  group_by(DoubPop, SingPop, Pair, Type) %>%
  summarise(
    n = n(),
    prop_sig = mean(Sig),
    mean_Z = mean(Z),
    median_Z = median(Z),
    mean_D = mean(D),
    median_D = median(D),
    .groups = "drop"
  )
write.table(sumtab, sumout, sep="\t", row.names=FALSE, quote=FALSE)

# Wilcoxon Right vs Wrong per Pair (BH correction)
wil <- df %>%
  filter(Type %in% c("Right","Wrong")) %>%
  group_by(Pair, DoubPop, SingPop) %>%
  summarise(
    n_right = sum(Type=="Right"),
    n_wrong = sum(Type=="Wrong"),
    p_value = if (n_right >= 3 && n_wrong >= 3) wilcox.test(Z ~ Type)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, "BH"))
write.table(wil, wilout, sep="\t", row.names=FALSE, quote=FALSE)

# Z plot for ALL pairs, includes Within. Cap y-axis to 99th percentile.
zcap <- quantile(abs(df$Z), probs=0.99, na.rm=TRUE)
zcap <- max(zcap, zthr*2)
npairs <- length(unique(df$Pair))
plot_width <- max(11, min(90, 0.30 * npairs)) # scale width with number of pairs

pZ <- ggplot(df, aes(x=Pair, y=Z, fill=Type)) +
  geom_boxplot(outlier.shape=NA, alpha=0.85) +
  geom_hline(yintercept=c(-zthr, zthr), linetype="dashed") +
  coord_cartesian(ylim=c(-zcap, zcap)) +
  theme_bw(base_size=11) +
  theme(
    legend.position="top",
    axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank()
  ) +
  labs(
    title="Topology test: Z distributions for ALL population pairs",
    subtitle=paste0("Pair = DoubPop_vs_SingPop (DoubPop=H1/H2; SingPop=H3). Y capped at 99th percentile: ±", signif(zcap,4)),
    x="", y="Z"
  )
ggsave(plot_zpairs, pZ, width=plot_width, height=8, limitsize=FALSE)

# ---- Mean Z heatmaps: ONE PDF, one page per Type, per-page scaling ----
heat <- sumtab %>% mutate(Type = factor(Type, levels=c("Right","Wrong","Within")))

pdf(plot_heat_meanZ, width=11, height=7)
for (tt in levels(heat$Type)) {

  sub <- heat %>% filter(Type == tt)
  if (nrow(sub) == 0) next

  L <- max(abs(sub$mean_Z), na.rm=TRUE)
  L <- max(L, 1e-6)

  pMean <- ggplot(sub, aes(x=SingPop, y=DoubPop, fill=mean_Z)) +
    geom_tile(color="white", linewidth=0.2) +
    scale_fill_gradient2(
      low="#3B4CC0", mid="white", high="#B40426", midpoint=0,
      limits=c(-L, L), oob=scales::squish
    ) +
    theme_bw(base_size=11) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=45, hjust=1)
    ) +
    labs(
      title=paste("Topology test mean Z:", tt),
      subtitle=paste0("Colour scale symmetric per page; limits ±", signif(L,3)),
      x="SingPop (H3)", y="DoubPop (H1/H2)", fill="Mean Z"
    )

  print(pMean)
}
dev.off()

# propSig heatmap: fixed 0..1 scale is comparable across Types
pProp <- ggplot(heat, aes(x=SingPop, y=DoubPop, fill=prop_sig)) +
  geom_tile(color="white", linewidth=0.2) +
  scale_fill_viridis_c(option="C", limits=c(0,1)) +
  facet_wrap(~Type) +
  theme_bw(base_size=11) +
  theme(
    panel.grid=element_blank(),
    axis.text.x=element_text(angle=45, hjust=1)
  ) +
  labs(
    title=paste0("Topology test: proportion with |Z| > ", zthr),
    x="SingPop (H3)", y="DoubPop (H1/H2)", fill="Prop. sig"
  )
ggsave(plot_heat_prop, pProp, width=11, height=7)

RSCRIPT

# ------------------------
# R: real comparisons summary + plots (optional)
# - Canonicalised already in AWK (H1pop<=H2pop)
# - Writes ONE multi-page PDF: one page per H3 (top 12 by n), per-page scaling
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

sumtab <- df %>%
  group_by(H1pop, H2pop, H3pop) %>%
  summarise(
    n = n(),
    prop_sig = mean(abs(Z) > zthr),
    mean_Z = mean(Z),
    median_Z = median(Z),
    mean_D = mean(D),
    median_D = median(D),
    .groups="drop"
  )
write.table(sumtab, outsum, sep="\t", row.names=FALSE, quote=FALSE)

topH3 <- sumtab %>%
  count(H3pop, wt=n, name="total_n") %>%
  arrange(desc(total_n)) %>%
  slice_head(n=12) %>%
  pull(H3pop)

plotdat <- sumtab %>% filter(H3pop %in% topH3)

pdf(outplot, width=11, height=7)
for (h3 in unique(plotdat$H3pop)) {

  sub <- plotdat %>% filter(H3pop == h3)
  if (nrow(sub) == 0) next

  L <- max(abs(sub$mean_Z), na.rm=TRUE)
  L <- max(L, 1e-6)

  p <- ggplot(sub, aes(x=H2pop, y=H1pop, fill=mean_Z)) +
    geom_tile(color="white", linewidth=0.2) +
    scale_fill_gradient2(
      low="#3B4CC0", mid="white", high="#B40426", midpoint=0,
      limits=c(-L, L), oob=scales::squish
    ) +
    theme_bw(base_size=11) +
    theme(
      panel.grid=element_blank(),
      axis.text.x=element_text(angle=45, hjust=1)
    ) +
    labs(
      title="Real D-stat mean Z (canonical H1/H2)",
      subtitle=paste0("H3 = ", h3, " | Colour scale symmetric per page; limits ±", signif(L,3)),
      x="H2 population (canonical)", y="H1 population (canonical)", fill="Mean Z"
    )

  print(p)
}
dev.off()

RSCRIPT2
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
  echo "  $realsum"
  echo "  $realplot  (multi-page: top H3 by n)"
fi
