#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
${script_name}: D-stat filtering, plotting, topology tests, ranking, and PCA.

Usage:
  ${script_name} -i dstats.jack.txt -m sample_to_pop.txt -p out_prefix [options]

Required:
  -i FILE                 ANGSD block-jackknifed D-stats table.
                          Header must include: H1 H2 H3 SE Z and one of jackEst/Dstat/D.
  -m FILE                 Sample-to-pop map (2 columns): SampleID Population
  -p PREFIX               Output prefix

Core filters:
  -s INT                  Minimum informative sites (default: 0)
                          Uses nABBA+nBABA when available, else nSites, else 0.
  -z FLOAT                |Z| threshold for "significant" points (default: 3)

Main output modes (all off by default):
  -R, --triples           Write distinct-population triple comparisons:
                          <prefix>.real.tsv, <prefix>.real.summary.tsv (legacy filenames)
  --topology              Write topology pairplot table + PDF
  --votes                 Write trio votes + support matrix (requires -R/--triples)
  --pcaPop[=POP_IDS]      Population-level ordination from support matrix (requires -R/--triples)
                          Optional value filters which populations to include by Pop ID.
                          Example: --pcaPop=BlackBear,BrownBear,PolarBear
  --pcaInd[=SAMPLE_IDS]   Individual-level PCA from H3-only affinities (requires -R/--triples)
                          Optional value filters which H3 individuals to include by sample ID.
                          Example: --pcaInd=Kumamoto_his,Saitama-218_his
  --rankH3Pop POP         Rank closest populations for rows where H3 is population POP
  --rankH3Ind SAMPLE      Rank closest populations for rows where H3 is SAMPLE
  --plotH3Ind LIST        Plot D by selected H3 individuals (comma-separated IDs)
  --plotRealByTree FILE   Filter distinct-population triple comparisons by rooted tree topology (requires -R/--triples)
                          Tree must include a tip named exactly: Outgroup

PCA controls:
  --pcaK INT              Number of PCs to write (default: 6, minimum: 2)
  --pcaPlotPC A,B         PCs to plot (default: 1,2; example: --pcaPlotPC 2,3)

Examples:
  ${script_name} -i in.jack.txt -m sample_to_pop.txt -p run1 --triples --votes --pcaPop
  ${script_name} -i in.jack.txt -m sample_to_pop.txt -p run2 --triples --pcaInd --pcaPlotPC 1,2

Notes:
  * Significance uses ANGSD block-jackknife Z (not IID bootstrap).
  * Low effective site counts can inflate Z. Use -s to reduce this.
EOF
  exit 1
}

# ------------------------
# Defaults
# ------------------------
infile=""
mapfile=""
prefix=""
script_name="$(basename "$0")"
minsites="0"
zthr="3"
do_real="0"
do_topology="0"
do_votes="0"
rankH3Pop=""
rankH3Ind=""
do_pca_pop="0"
pca_pop_ids_csv=""
pca_ind_mode="0"
pca_ind_samples_csv=""
plotH3Ind_csv=""
plotRealByTree_nwk=""

pcaK="6"
pcaPlotPC="1,2"

# ------------------------
# Arg parse
# ------------------------
if [[ $# -eq 0 ]]; then usage; fi

require_arg() {
  local opt="$1"
  local val="${2-}"
  if [[ -z "$val" ]]; then
    echo "ERROR: $opt requires a value." >&2
    usage
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) require_arg "$1" "${2-}"; infile="$2"; shift 2 ;;
    -m) require_arg "$1" "${2-}"; mapfile="$2"; shift 2 ;;
    -p) require_arg "$1" "${2-}"; prefix="$2"; shift 2 ;;
    -s) require_arg "$1" "${2-}"; minsites="$2"; shift 2 ;;
    -z) require_arg "$1" "${2-}"; zthr="$2"; shift 2 ;;
    -R|--triples) do_real="1"; shift 1 ;;
    --topology) do_topology="1"; shift 1 ;;
    --votes) do_votes="1"; shift 1 ;;
    --rankH3Pop) require_arg "$1" "${2-}"; rankH3Pop="$2"; shift 2 ;;
    --rankH3Ind) require_arg "$1" "${2-}"; rankH3Ind="$2"; shift 2 ;;
    --pcaPop)
      do_pca_pop="1"
      if [[ "${2:-}" =~ ^[^-].* ]] && [[ "${2:-}" != "" ]]; then
        pca_pop_ids_csv="$2"
        shift 2
      else
        shift 1
      fi
      ;;
    --pcaPop=*)
      do_pca_pop="1"
      pca_pop_ids_csv="${1#*=}"
      shift 1
      ;;
    --pcaK) require_arg "$1" "${2-}"; pcaK="$2"; shift 2 ;;
    --pcaPlotPC) require_arg "$1" "${2-}"; pcaPlotPC="$2"; shift 2 ;;
    --pcaInd)
      pca_ind_mode="1"
      if [[ "${2:-}" =~ ^[^-].* ]] && [[ "${2:-}" != "" ]]; then
        pca_ind_samples_csv="$2"
        shift 2
      else
        shift 1
      fi
      ;;
    --pcaInd=*)
      pca_ind_mode="1"
      pca_ind_samples_csv="${1#*=}"
      shift 1
      ;;
    --plotH3Ind)
      require_arg "$1" "${2-}"
      plotH3Ind_csv="$2"
      shift 2
      ;;
    --plotH3Ind=*)
      plotH3Ind_csv="${1#*=}"
      shift 1
      ;;
    --plotRealByTree)
      require_arg "$1" "${2-}"
      plotRealByTree_nwk="$2"
      shift 2
      ;;
    --plotRealByTree=*)
      plotRealByTree_nwk="${1#*=}"
      shift 1
      ;;
    -h|--help) usage ;;
    *)
      echo "ERROR: Unknown option: $1" >&2
      usage
      ;;
  esac
done

# ------------------------
# Validate
# ------------------------
[[ -z "${infile:-}" || -z "${mapfile:-}" || -z "${prefix:-}" ]] && usage
[[ ! -f "$infile" ]] && { echo "ERROR: infile not found: $infile" >&2; exit 2; }
[[ ! -f "$mapfile" ]] && { echo "ERROR: mapfile not found: $mapfile" >&2; exit 2; }

for cmd in awk Rscript; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "ERROR: Required command not found: $cmd" >&2
    exit 2
  fi
done

if ! [[ "$minsites" =~ ^[0-9]+$ ]]; then
  echo "ERROR: -s must be a non-negative integer." >&2
  exit 2
fi
if ! [[ "$zthr" =~ ^[-+]?[0-9]*\.?[0-9]+$ ]]; then
  echo "ERROR: -z must be numeric." >&2
  exit 2
fi

if ! [[ "$pcaK" =~ ^[0-9]+$ ]] || [[ "$pcaK" -lt 2 ]]; then
  echo "ERROR: --pcaK must be an integer >= 2." >&2
  exit 2
fi

if ! [[ "$pcaPlotPC" =~ ^[0-9]+,[0-9]+$ ]]; then
  echo "ERROR: --pcaPlotPC must be like A,B (e.g. 1,2 or 2,3)." >&2
  exit 2
fi

if [[ "$do_votes" == "1" || "$do_pca_pop" == "1" || "$pca_ind_mode" == "1" || -n "${plotRealByTree_nwk:-}" ]]; then
  if [[ "$do_real" != "1" ]]; then
    echo "ERROR: --votes/--pcaPop/--pcaInd/--plotRealByTree require -R or --triples." >&2
    exit 2
  fi
fi

if [[ -n "${plotRealByTree_nwk:-}" && ! -f "$plotRealByTree_nwk" ]]; then
  echo "ERROR: --plotRealByTree file not found: $plotRealByTree_nwk" >&2
  exit 2
fi

# ------------------------
# Outputs
# ------------------------
topo_raw="${prefix}.topology.mapped_rows.tsv"
topo_pairplot_tsv="${prefix}.topology.pairplot.tsv"
topo_pairplot_pdf="${prefix}.topology.pairplot.pdf"

realout="${prefix}.real.tsv"
realsum="${prefix}.real.summary.tsv"

votes_out="${prefix}.votes.tsv"
support_out="${prefix}.support_matrix.tsv"

pca_pop_pdf="${prefix}.pca.pop.pdf"
pca_pop_tsv="${prefix}.pca.pop.tsv"

pca_ind_pdf="${prefix}.pca.ind.pdf"
pca_ind_tsv="${prefix}.pca.ind.tsv"

rankH3Pop_out="${prefix}.rankH3Pop.tsv"
rankH3Ind_out="${prefix}.rankH3Ind.tsv"

plotH3Ind_tsv="${prefix}.plotH3Ind.tsv"
plotH3Ind_pdf="${prefix}.plotH3Ind.pdf"

realByTree_tsv="${prefix}.real.byTree.tsv"
realByTree_pdf="${prefix}.real.byTree.pdf"

# ------------------------
# AWK core
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
  nabba=idx["nABBA"]; nbaba=idx["nBABA"]; nsites=idx["nSites"]
  if (!h1||!h2||!h3||!se||!z||!d) {
    print "ERROR: Missing required columns in header (need H1 H2 H3 SE Z and jackEst/Dstat/D)" > "/dev/stderr"; exit 10
  }
  next
}

function pop(s){return (s in samp2pop ? samp2pop[s] : "NA")}
function get_nsites(){
  if (nabba && nbaba) return ($(nabba)+0)+($(nbaba)+0)
  if (nsites) return $(nsites)+0
  return 0
}

{
  H1=$(h1); H2=$(h2); H3=$(h3)
  P1=pop(H1); P2=pop(H2); P3=pop(H3)
  if (P1=="NA"||P2=="NA"||P3=="NA") next
  Dv=$(d)+0; SEv=$(se)+0; Zv=$(z)+0
  Ns=get_nsites()
  if (Ns < minsites) next

  if (mode=="real") {
    if (P1!=P2 && P1!=P3 && P2!=P3) {
      if (P1 < P2)
        print P1,P2,P3,H1,H2,H3,Dv,SEv,Zv,Ns
      else
        print P2,P1,P3,H2,H1,H3,-Dv,SEv,-Zv,Ns
    }
  }
  else if (mode=="topology_raw") {
    print H1,H2,H3,P1,P2,P3,Dv,SEv,Zv,Ns
  }
  else if (mode=="plotH3Ind_raw") {
    print H1,H2,H3,P1,P2,P3,Dv,SEv,Zv,Ns
  }
}
'

# ------------------------
# Distinct-population triples (-R/--triples)
# ------------------------
if [[ "$do_real" == "1" ]]; then
  printf "H1pop\tH2pop\tH3pop\tH1\tH2\tH3\tD\tSE\tZ\tnSites\n" > "$realout"
  awk -v mode="real" -v minsites="$minsites" "$run_awk" "$mapfile" "$infile" >> "$realout"

  Rscript --vanilla - "$realout" "$realsum" "$zthr" <<'RS_REALSUM'
required_pkgs <- c("dplyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outsum <- args[2]
zthr <- as.numeric(args[3])

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)

sumtab <- df %>%
  group_by(H1pop, H2pop, H3pop) %>%
  summarise(
    n = n(),
    med_nSites = suppressWarnings(median(nSites, na.rm=TRUE)),
    min_nSites = suppressWarnings(min(nSites, na.rm=TRUE)),
    prop_sig_positive = mean(Z >  zthr),
    prop_sig_negative = mean(Z < -zthr),
    prop_sig_total    = mean(abs(Z) > zthr),
    mean_Z = mean(Z),
    median_Z = median(Z),
    mean_D = mean(D),
    median_D = median(D),
    .groups="drop"
  )

write.table(sumtab, outsum, sep="\t", row.names=FALSE, quote=FALSE)
RS_REALSUM
fi

# ------------------------
# REAL-by-tree
# ------------------------
if [[ -n "${plotRealByTree_nwk:-}" ]]; then
  Rscript --vanilla - "$realout" "$plotRealByTree_nwk" "$realByTree_tsv" "$realByTree_pdf" "$zthr" <<'RS_REALBYTREE'
required_pkgs <- c("ape","dplyr","ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
realfile <- args[1]
nwkfile  <- args[2]
outtsv   <- args[3]
outpdf   <- args[4]
zthr <- as.numeric(args[5])

real <- read.delim(realfile, sep="\t", header=TRUE, stringsAsFactors=FALSE)

tr <- ape::read.tree(nwkfile)
if (is.null(tr)) stop("Could not read tree.")
if (!ape::is.binary(tr)) stop("Tree must be strictly binary (fully resolved).")
if (!("Outgroup" %in% tr$tip.label)) stop("Tree must contain a tip named exactly: Outgroup")

tr <- ape::root(tr, outgroup="Outgroup", resolve.root=TRUE)

ntip <- length(tr$tip.label)
internal_nodes <- (ntip+1):(ntip+tr$Nnode)

get_tips <- function(tree, node, ntip) {
  if (node <= ntip) return(tree$tip.label[node])
  ape::extract.clade(tree, node)$tip.label
}

tops <- list()
for (node in internal_nodes) {
  kids <- tr$edge[tr$edge[,1]==node, 2]
  if (length(kids) != 2) next

  tips1 <- get_tips(tr, kids[1], ntip)
  tips2 <- get_tips(tr, kids[2], ntip)

  if (length(tips1) == 2) {
    sis <- sort(tips1)
    others <- setdiff(tr$tip.label, tips1)
    others <- setdiff(others, "Outgroup")
    for (h3 in others) tops[[length(tops)+1]] <- data.frame(sis1=sis[1], sis2=sis[2], H3=h3)
  }
  if (length(tips2) == 2) {
    sis <- sort(tips2)
    others <- setdiff(tr$tip.label, tips2)
    others <- setdiff(others, "Outgroup")
    for (h3 in others) tops[[length(tops)+1]] <- data.frame(sis1=sis[1], sis2=sis[2], H3=h3)
  }
}

if (length(tops) == 0) stop("Tree produced no (sister pair)|H3 topologies with a 2-tip sister clade.")

topdf <- dplyr::bind_rows(tops) %>% distinct()

real2 <- real %>% mutate(H1c=pmin(H1pop,H2pop), H2c=pmax(H1pop,H2pop))

m <- real2 %>%
  inner_join(topdf, by=c("H1c"="sis1","H2c"="sis2","H3pop"="H3")) %>%
  mutate(
    Topology = paste0("(", H1c, ",", H2c, ")|", H3pop),
    sig = ifelse(abs(Z) > zthr, "sig", "non")
  )

if (nrow(m) == 0) stop("No distinct-population triple comparisons matched the supplied tree topology. Check labels.")

out <- m %>% select(Topology, H1pop, H2pop, H3pop, H1, H2, H3, D, SE, Z, nSites, sig)
write.table(out, outtsv, sep="\t", quote=FALSE, row.names=FALSE)

lev <- out %>% distinct(Topology) %>% arrange(Topology) %>% pull(Topology)
out$Topology <- factor(out$Topology, levels=rev(lev))

pdf(outpdf, width=10, height=max(6, 0.18*length(lev) + 2))
p <- ggplot(out, aes(x=D, y=Topology)) +
  geom_jitter(aes(shape=sig), size=2.2, height=0.15, width=0,
              stroke=0.9, color="black", fill=NA) +
  scale_shape_manual(values=c("sig"=19, "non"=21)) +
  theme_bw(base_size=11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(linetype="dashed")) +
  labs(
    title="Real D-stat comparisons consistent with supplied tree",
    subtitle=paste0("Filled = |Z|>", zthr, " ; Open = |Z|<=", zthr, "  |  Label: (sister1,sister2)|H3"),
    x="D-score", y=""
  )
print(p)
dev.off()
RS_REALBYTREE
fi

# ------------------------
# TOPOLOGY plot
# ------------------------
if [[ "$do_topology" == "1" ]]; then
  printf "H1\tH2\tH3\tH1pop\tH2pop\tH3pop\tD\tSE\tZ\tnSites\n" > "$topo_raw"
  awk -v mode="topology_raw" -v minsites="$minsites" "$run_awk" "$mapfile" "$infile" >> "$topo_raw"

  Rscript --vanilla - "$topo_raw" "$topo_pairplot_tsv" "$topo_pairplot_pdf" <<'RS_TOPOPAIR'
required_pkgs <- c("dplyr","ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outtsv <- args[2]
outpdf <- args[3]

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)

blue <- df %>%
  filter(H1pop==H2pop, H3pop!=H1pop) %>%
  transmute(Pop1=H1pop, Pop2=H3pop, D=D, Source="(Pop1,Pop1),Pop2")

red <- df %>%
  filter(H1pop!=H2pop, H3pop==H1pop) %>%
  transmute(Pop1=H1pop, Pop2=H2pop, D=D, Source="(Pop1,Pop2),Pop1")

plotdat <- bind_rows(blue, red) %>% mutate(Pair=paste0(Pop1," - ",Pop2))
write.table(plotdat, outtsv, sep="\t", quote=FALSE, row.names=FALSE)

pair_levels <- plotdat %>% distinct(Pair,Pop1,Pop2) %>% arrange(Pop1,Pop2) %>% pull(Pair)
plotdat$Pair <- factor(plotdat$Pair, levels=rev(pair_levels))

pdf(outpdf, width=10, height=max(6, 0.18*length(pair_levels) + 2))
p <- ggplot(plotdat, aes(x=D, y=Pair)) +
  geom_jitter(aes(color=Source), height=0.15, width=0, size=2.2, shape=21, fill=NA, stroke=0.9) +
  theme_bw(base_size=11) +
  theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
        panel.grid.major.x=element_line(linetype="dashed")) +
  labs(
    title="Topology comparison per population pair",
    subtitle="Blue open circles: (Pop1,Pop1),Pop2   |   Red open circles: (Pop1,Pop2),Pop1",
    x="D-score", y=""
  ) +
  scale_color_manual(values=c("(Pop1,Pop1),Pop2"="blue", "(Pop1,Pop2),Pop1"="red"))
print(p)
dev.off()
RS_TOPOPAIR
fi

# ------------------------
# Votes + support matrix (needed for PCA)
# ------------------------
if [[ "$do_votes" == "1" || "$do_pca_pop" == "1" || "$pca_ind_mode" == "1" ]]; then
  Rscript --vanilla - "$realsum" "$votes_out" "$support_out" <<'RS_VOTES'
required_pkgs <- c("dplyr","tidyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))

args <- commandArgs(trailingOnly=TRUE)
realsum <- args[1]
out_votes <- args[2]
out_support <- args[3]

tab <- read.delim(realsum, sep="\t", header=TRUE, stringsAsFactors=FALSE)
stopifnot(all(c("H1pop","H2pop","H3pop","mean_Z") %in% names(tab)))

tab <- tab %>% mutate(score = abs(mean_Z))

pops <- sort(unique(c(tab$H1pop, tab$H2pop, tab$H3pop)))
n <- length(pops)
if (n < 4) stop("Need >= 4 populations for trio voting / support matrix.")

key <- paste(tab$H1pop, tab$H2pop, tab$H3pop, sep="|")
score_map <- setNames(tab$score, key)

get_score <- function(A, B, C) {
  x1 <- ifelse(A < B, A, B)
  x2 <- ifelse(A < B, B, A)
  k <- paste(x1, x2, C, sep="|")
  v <- score_map[[k]]
  if (is.null(v) || length(v)==0) return(NA_real_)
  as.numeric(v)
}

trios <- list(); votes <- list()
ti <- 1; vi <- 1

for (i in 1:(n-2)) for (j in (i+1):(n-1)) for (k in (j+1):n) {
  A <- pops[i]; Bp <- pops[j]; C <- pops[k]
  sAB_C <- get_score(A,  Bp, C)
  sAC_B <- get_score(A,  C,  Bp)
  sBC_A <- get_score(Bp, C,  A)
  if (any(is.na(c(sAB_C, sAC_B, sBC_A)))) next

  scores <- c(sAB_C, sAC_B, sBC_A)
  w <- which.min(scores)

  ord <- sort(scores)
  winner <- ord[1]; runnerup <- ord[2]
  margin <- runnerup - winner

  if (w == 1) sister <- c(A, Bp)
  if (w == 2) sister <- c(A, C)
  if (w == 3) sister <- c(Bp, C)

  trios[[ti]] <- c(A, Bp, C); ti <- ti + 1

  votes[[vi]] <- data.frame(
    A=A, B=Bp, C=C,
    best_sister1=sister[1], best_sister2=sister[2],
    score_AB_C=sAB_C, score_AC_B=sAC_B, score_BC_A=sBC_A,
    winner=winner, runnerup=runnerup, margin=margin,
    stringsAsFactors=FALSE
  ); vi <- vi + 1
}

if (length(trios) < 1) stop("No complete trios found (need all 3 resolutions per trio).")

votes_df <- bind_rows(votes)
write.table(votes_df, out_votes, sep="\t", quote=FALSE, row.names=FALSE)

wins <- matrix(0, n, n, dimnames=list(pops, pops))
den  <- matrix(0, n, n, dimnames=list(pops, pops))

for (t in seq_along(trios)) {
  A <- trios[[t]][1]; Bp <- trios[[t]][2]; C <- trios[[t]][3]
  scores <- c(get_score(A,Bp,C), get_score(A,C,Bp), get_score(Bp,C,A))

  for (pp in list(c(A,Bp), c(A,C), c(Bp,C))) {
    den[pp[1], pp[2]] <- den[pp[1], pp[2]] + 1
    den[pp[2], pp[1]] <- den[pp[2], pp[1]] + 1
  }

  w <- which.min(scores)
  if (w == 1) sister <- c(A, Bp)
  if (w == 2) sister <- c(A, C)
  if (w == 3) sister <- c(Bp, C)

  wins[sister[1], sister[2]] <- wins[sister[1], sister[2]] + 1
  wins[sister[2], sister[1]] <- wins[sister[2], sister[1]] + 1
}

support <- wins / den
support[is.na(support)] <- 0
diag(support) <- 1

write.table(support, out_support, sep="\t", quote=FALSE, col.names=NA)
RS_VOTES
fi

# ------------------------
# pcaPop (MDS)
# ------------------------
if [[ "$do_pca_pop" == "1" ]]; then
  Rscript --vanilla - "$support_out" "$pca_pop_pdf" "$pca_pop_tsv" "$pca_pop_ids_csv" <<'RS_PCAPOP'
required_pkgs <- c("ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
inmat <- args[1]
outpdf <- args[2]
outtsv <- args[3]
pops_csv <- args[4]

mat <- read.delim(inmat, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
S <- as.matrix(mat)
if (!is.null(pops_csv) && nchar(pops_csv) > 0) {
  keep_pops <- trimws(unlist(strsplit(pops_csv, ",", fixed=TRUE)))
  unknown_pops <- setdiff(keep_pops, rownames(S))
  if (length(unknown_pops) > 0) {
    stop(paste0("Unknown population ID(s) in --pcaPop: ", paste(unknown_pops, collapse=", ")))
  }
  S <- S[keep_pops, keep_pops, drop=FALSE]
}
if (nrow(S) < 3) stop("pcaPop requires at least 3 populations after filtering.")
D <- 1 - S
fit <- cmdscale(as.dist(D), k=2, eig=TRUE)
coords <- data.frame(Pop=rownames(S), Axis1=fit$points[,1], Axis2=fit$points[,2], stringsAsFactors=FALSE)
write.table(coords, outtsv, sep="\t", quote=FALSE, row.names=FALSE)

pdf(outpdf, width=7, height=6)
p <- ggplot(coords, aes(x=Axis1, y=Axis2, label=Pop)) +
  geom_point(size=2.5) +
  geom_text(vjust=-0.6, size=3) +
  theme_bw(base_size=11) +
  labs(title="Population ordination from trio-vote support (classical MDS)", x="Axis 1", y="Axis 2")
print(p)
dev.off()
RS_PCAPOP
fi

# ------------------------
# pcaInd (fixed aes_string + constant columns)
# ------------------------
if [[ "$pca_ind_mode" == "1" ]]; then
  Rscript --vanilla - "$infile" "$mapfile" "$pca_ind_pdf" "$pca_ind_tsv" "$minsites" "$pca_ind_samples_csv" "$pcaK" "$pcaPlotPC" <<'RS_PCAIND'
required_pkgs <- c("dplyr","ggplot2","tidyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
dstats <- args[1]
mapfile <- args[2]
outpdf <- args[3]
outtsv <- args[4]
minsites <- as.numeric(args[5])
samples_csv <- args[6]
K <- as.integer(args[7])
plotPC <- args[8]

pcs <- as.integer(strsplit(plotPC, ",", fixed=TRUE)[[1]])
pcA <- pcs[1]; pcB <- pcs[2]

map <- read.table(mapfile, header=FALSE, stringsAsFactors=FALSE)
colnames(map) <- c("Sample","Pop")

df <- read.delim(dstats, sep="", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
need <- c("H1","H2","H3","Z")
miss <- setdiff(need, names(df))
if (length(miss)>0) stop(paste("Missing required columns:", paste(miss, collapse=", ")))

if (all(c("nABBA","nBABA") %in% names(df))) {
  df$nSites <- df$nABBA + df$nBABA
} else if ("nSites" %in% names(df)) {
  df$nSites <- df$nSites
} else {
  df$nSites <- 0
}
df <- df %>% filter(nSites >= minsites)

pop_of <- function(x) map$Pop[match(x, map$Sample)]
df$H1pop <- pop_of(df$H1); df$H2pop <- pop_of(df$H2); df$H3pop <- pop_of(df$H3)
df <- df %>% filter(!is.na(H1pop), !is.na(H2pop), !is.na(H3pop))

if (!is.null(samples_csv) && nchar(samples_csv)>0) {
  keep_samps <- trimws(unlist(strsplit(samples_csv, ",", fixed=TRUE)))
  unknown_samps <- setdiff(keep_samps, unique(map$Sample))
  if (length(unknown_samps) > 0) {
    stop(paste0("Unknown sample ID(s) in --pcaInd: ", paste(unknown_samps, collapse=", ")))
  }
  if (length(keep_samps) < 1) {
    stop("No samples selected for --pcaInd filter.")
  }
  pops <- sort(unique(map$Pop[match(keep_samps, map$Sample)]))
} else {
  keep_samps <- map$Sample
  pops <- sort(unique(map$Pop))
}

sub <- df %>% filter(H3 %in% keep_samps)
if (nrow(sub) < 1) stop("No D-stat rows left for selected H3 samples after filtering.")

pair <- sub %>%
  filter(H1pop != H2pop)

# Build per-H3 population affinities symmetrically:
# for each row, H2 gets +Z support, H1 gets -Z support.
aff_long <- bind_rows(
  pair %>% transmute(H3, popX=H1pop, score=-Z),
  pair %>% transmute(H3, popX=H2pop, score= Z)
) %>%
  filter(popX %in% pops)

aff <- aff_long %>%
  group_by(H3, popX) %>%
  summarise(aff = mean(score), .groups="drop") %>%
  tidyr::pivot_wider(names_from=popX, values_from=aff, values_fill=0)

for (p in pops) if (!(p %in% names(aff))) aff[[p]] <- 0

mat <- as.matrix(aff[, pops, drop=FALSE])
rownames(mat) <- aff$H3

sds <- apply(mat, 2, sd)
keep_cols <- which(is.finite(sds) & sds > 0)
if (length(keep_cols) < 2) stop("Too few non-constant population columns for PCA (need >=2).")
mat2 <- mat[, keep_cols, drop=FALSE]

mat_s <- scale(mat2)
pc <- prcomp(mat_s, center=TRUE, scale.=FALSE)

K_use <- min(K, ncol(pc$x))
if (pcA > K_use || pcB > K_use) stop(paste0("Requested PCs ", pcA, ",", pcB, " but only ", K_use, " available."))

coords <- data.frame(Sample=rownames(mat_s), stringsAsFactors=FALSE)
for (k in 1:K_use) coords[[paste0("PC",k)]] <- pc$x[,k]
coords$Pop <- map$Pop[match(coords$Sample, map$Sample)]

write.table(coords, outtsv, sep="\t", quote=FALSE, row.names=FALSE)

pdf(outpdf, width=7, height=6)
p <- ggplot(coords, aes(x=.data[[paste0("PC",pcA)]], y=.data[[paste0("PC",pcB)]], color=Pop)) +
  geom_point(size=2.4) +
  theme_bw(base_size=11) +
  labs(title="Individuals PCA from H3-only affinity (mean Z vs populations)",
       x=paste0("PC",pcA), y=paste0("PC",pcB))
print(p)
dev.off()
RS_PCAIND
fi

# ------------------------
# rankH3*
# ------------------------
if [[ -n "${rankH3Pop:-}" || -n "${rankH3Ind:-}" ]]; then
  Rscript --vanilla - "$infile" "$mapfile" "$minsites" "$rankH3Pop" "$rankH3Ind" "$rankH3Pop_out" "$rankH3Ind_out" <<'RS_RANK'
required_pkgs <- c("dplyr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))

args <- commandArgs(trailingOnly=TRUE)
dstats <- args[1]
mapfile <- args[2]
minsites <- as.numeric(args[3])
targetPop <- args[4]
targetInd <- args[5]
outPop <- args[6]
outInd <- args[7]

map <- read.table(mapfile, header=FALSE, stringsAsFactors=FALSE)
colnames(map) <- c("Sample","Pop")

df <- read.delim(dstats, sep="", header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
need <- c("H1","H2","H3","Z")
miss <- setdiff(need, names(df))
if (length(miss)>0) stop(paste("Missing required columns:", paste(miss, collapse=", ")))

if (all(c("nABBA","nBABA") %in% names(df))) {
  df$nSites <- df$nABBA + df$nBABA
} else if ("nSites" %in% names(df)) {
  df$nSites <- df$nSites
} else {
  df$nSites <- 0
}

df <- df %>% filter(nSites >= minsites)

pop_of <- function(x) map$Pop[match(x, map$Sample)]
df$H1pop <- pop_of(df$H1); df$H2pop <- pop_of(df$H2); df$H3pop <- pop_of(df$H3)
df <- df %>% filter(!is.na(H1pop), !is.na(H2pop), !is.na(H3pop))

rank_from_H3 <- function(subdf) {
  subdf %>%
    filter(H1pop != H2pop) %>%
    transmute(winner = ifelse(Z > 0, H2pop, H1pop), w = abs(Z)) %>%
    group_by(winner) %>%
    summarise(wins_wt=sum(w), n_wins=n(), .groups="drop") %>%
    rename(Candidate=winner) %>%
    arrange(desc(wins_wt), desc(n_wins))
}

if (nchar(targetPop)>0) {
  out <- rank_from_H3(df %>% filter(H3pop == targetPop))
  write.table(out, outPop, sep="\t", quote=FALSE, row.names=FALSE)
}
if (nchar(targetInd)>0) {
  out <- rank_from_H3(df %>% filter(H3 == targetInd))
  write.table(out, outInd, sep="\t", quote=FALSE, row.names=FALSE)
}
RS_RANK
fi

# ------------------------
# plotH3Ind (FIXED layout)
# ------------------------
if [[ -n "${plotH3Ind_csv:-}" ]]; then
  plot_raw="${prefix}.plotH3Ind.mapped_rows.tsv"
  printf "H1\tH2\tH3\tH1pop\tH2pop\tH3pop\tD\tSE\tZ\tnSites\n" > "$plot_raw"
  awk -v mode="plotH3Ind_raw" -v minsites="$minsites" "$run_awk" "$mapfile" "$infile" >> "$plot_raw"

  Rscript --vanilla - "$plot_raw" "$plotH3Ind_tsv" "$plotH3Ind_pdf" "$plotH3Ind_csv" "$zthr" <<'RS_PLOTH3'
required_pkgs <- c("dplyr","ggplot2")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing) > 0) stop(paste0("Missing R package(s): ", paste(missing, collapse=", "), ". Install them before running this script."))
invisible(lapply(required_pkgs, library, character.only=TRUE))
try(suppressWarnings(Sys.setlocale("LC_CTYPE","C")), silent=TRUE)

args <- commandArgs(trailingOnly=TRUE)
infile <- args[1]
outtsv <- args[2]
outpdf <- args[3]
h3csv <- args[4]
zthr <- as.numeric(args[5])

df <- read.delim(infile, sep="\t", header=TRUE, stringsAsFactors=FALSE)
h3list <- trimws(unlist(strsplit(h3csv, ",", fixed=TRUE)))

sub <- df %>%
  filter(H3 %in% h3list) %>%
  filter(H1pop != H2pop) %>%
  filter(!(H1 %in% h3list), !(H2 %in% h3list)) %>%
  mutate(
    H1H2 = paste0(H1pop, "_", H2pop),
    sig = ifelse(abs(Z) > zthr, "sig", "non")
  )

outtab <- sub %>% select(H3, H3pop, H1H2, D, Z, H1, H2, H1pop, H2pop, nSites, sig)
write.table(outtab, outtsv, sep="\t", quote=FALSE, row.names=FALSE)

# y = H3 individual; colour = H1H2 comparison; shape = sig/non (filled/open)
outtab$H3 <- factor(outtab$H3, levels=rev(unique(outtab$H3)))

pdf(outpdf, width=11, height=max(6, 0.25*length(levels(outtab$H3)) + 2))
p <- ggplot(outtab, aes(x=D, y=H3)) +
  geom_jitter(aes(color=H1H2, shape=sig), size=2.4, height=0.15, width=0) +
  scale_shape_manual(values=c("sig"=19, "non"=21)) +
  theme_bw(base_size=11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype="dashed"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title="D-score for selected H3 individuals",
    subtitle=paste0("H3 in: ", paste(h3list, collapse=", "), " | filled=|Z|>", zthr, ", open=|Z|<=", zthr),
    x="D-score", y="H3 individual", color="H1pop_H2pop"
  )
print(p)
dev.off()
RS_PLOTH3
fi

echo "Done."
echo "Wrote (depending on options):"
[[ "$do_real" == "1" ]] && echo "  $realout" && echo "  $realsum"
[[ -n "${plotRealByTree_nwk:-}" ]] && echo "  $realByTree_tsv" && echo "  $realByTree_pdf"
[[ "$do_topology" == "1" ]] && echo "  $topo_raw" && echo "  $topo_pairplot_tsv" && echo "  $topo_pairplot_pdf"
[[ "$do_votes" == "1" || "$do_pca_pop" == "1" || "$pca_ind_mode" == "1" ]] && echo "  $votes_out" && echo "  $support_out"
[[ "$do_pca_pop" == "1" ]] && echo "  $pca_pop_pdf" && echo "  $pca_pop_tsv"
[[ "$pca_ind_mode" == "1" ]] && echo "  $pca_ind_pdf" && echo "  $pca_ind_tsv"
[[ -n "${rankH3Pop:-}" ]] && echo "  $rankH3Pop_out"
[[ -n "${rankH3Ind:-}" ]] && echo "  $rankH3Ind_out"
[[ -n "${plotH3Ind_csv:-}" ]] && echo "  $plotH3Ind_pdf" && echo "  $plotH3Ind_tsv"
