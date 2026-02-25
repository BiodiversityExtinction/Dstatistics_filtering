#!/usr/bin/env Rscript

fail <- function(msg, code = 2L) {
  cat(sprintf("ERROR: %s\n", msg), file = stderr())
  quit(save = "no", status = code)
}

usage <- function(script_name, code = 1L) {
  txt <- sprintf(paste(
    "%s: D-stat filtering, plotting, topology tests, ranking, and PCA.",
    "",
    "Usage:",
    "  %s --input dstats.jack.txt --Popfile sample_to_pop.txt --out out_prefix [options]",
    "",
    "Required:",
    "  --input FILE            ANGSD block-jackknifed D-stats table.",
    "                          Header must include: H1 H2 H3 SE Z and one of jackEst/Dstat/D.",
    "  --Popfile FILE          Sample-to-pop map (2 columns): SampleID Population",
    "  --out PREFIX            Output prefix",
    "",
    "Core filters:",
    "  --minsites INT          Minimum informative sites (default: 0)",
    "                          Uses nABBA+nBABA when available, else nSites, else 0.",
    "  -z FLOAT                |Z| threshold for significant points (default: 3)",
    "",
    "Main output modes (all off by default):",
    "  --triples               Write distinct-population triple comparisons:",
    "                          <prefix>.real.tsv, <prefix>.real.summary.tsv",
    "  --topology              Write topology pairplot table + PDF",
    "  --votes                 Write trio votes + support matrix (requires --triples; writes Z- and D-based outputs)",
    "  --pcaPop[=POP_IDS]      Population-level ordination from support matrix (requires --triples; writes Z- and D-based outputs)",
    "                          Optional value filters which populations to include by Pop ID.",
    "  --pcaInd[=SAMPLE_IDS]   Individual-level PCA from H3-only affinities (requires --triples; writes Z- and D-based outputs)",
    "                          Optional value filters which H3 individuals to include by sample ID.",
    "  --rankH3Pop POP         Rank closest populations for rows where H3 is population POP",
    "  --rankH3Ind SAMPLE      Rank closest populations for rows where H3 is SAMPLE",
    "  --plotH3Ind LIST        Plot D by selected H3 individuals (comma-separated IDs)",
    "  --plotRealByTree FILE   Filter distinct-population triple comparisons by rooted tree topology (requires --triples)",
    "                          Tree must include a tip named exactly: Outgroup",
    "",
    "PCA controls:",
    "  --pcaK INT              Number of PCs to write (default: 6, minimum: 2)",
    "  --pcaPlotPC A,B         PCs to plot (default: 1,2; example: --pcaPlotPC 2,3)",
    "",
    "Examples:",
    "  %s --input in.jack.txt --Popfile sample_to_pop.txt --out run1 --triples --votes --pcaPop",
    "  %s --input in.jack.txt --Popfile sample_to_pop.txt --out run2 --triples --pcaInd --pcaPlotPC 1,2",
    "",
    "Notes:",
    "  * Significance uses ANGSD block-jackknife Z (not IID bootstrap).",
    "  * Low effective site counts can inflate Z. Use --minsites to reduce this.",
    sep = "\n"
  ), script_name, script_name, script_name, script_name)
  cat(txt, "\n")
  quit(save = "no", status = code)
}

require_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    fail(sprintf("Missing R package(s): %s. Install them before running this script.", paste(missing, collapse = ", ")))
  }
}

parse_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(strsplit(x, ",", fixed = TRUE)[[1]])
}

cmd_all <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd_all[grep("^--file=", cmd_all)]
script_name <- if (length(file_arg) > 0) basename(sub("^--file=", "", file_arg[1])) else "DstatWorkbench.R"
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) usage(script_name, code = 1L)

opt <- list(
  infile = "",
  mapfile = "",
  prefix = "",
  minsites = 0L,
  zthr = 3,
  do_real = FALSE,
  do_topology = FALSE,
  do_votes = FALSE,
  rankH3Pop = "",
  rankH3Ind = "",
  do_pca_pop = FALSE,
  pca_pop_ids_csv = "",
  pca_ind_mode = FALSE,
  pca_ind_samples_csv = "",
  plotH3Ind_csv = "",
  plotRealByTree_nwk = "",
  pcaK = 6L,
  pcaPlotPC = "1,2"
)

need_value <- function(i, name) {
  if (i >= length(args) || startsWith(args[i + 1], "-")) fail(sprintf("%s requires a value.", name))
}

i <- 1L
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("-h", "--help")) usage(script_name, code = 0L)

  if (a == "--input") {
    need_value(i, a)
    opt$infile <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--Popfile") {
    need_value(i, a)
    opt$mapfile <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--out") {
    need_value(i, a)
    opt$prefix <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--minsites") {
    need_value(i, a)
    opt$minsites <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "-z") {
    need_value(i, a)
    opt$zthr <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--topology") { opt$do_topology <- TRUE; i <- i + 1L; next }
  if (a == "--triples") { opt$do_real <- TRUE; i <- i + 1L; next }
  if (a == "--votes") { opt$do_votes <- TRUE; i <- i + 1L; next }

  if (a == "--rankH3Pop") {
    need_value(i, a)
    opt$rankH3Pop <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--rankH3Ind") {
    need_value(i, a)
    opt$rankH3Ind <- args[i + 1]
    i <- i + 2L
    next
  }

  if (a == "--pcaPop") {
    opt$do_pca_pop <- TRUE
    if (i < length(args) && !startsWith(args[i + 1], "-")) {
      opt$pca_pop_ids_csv <- args[i + 1]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
    next
  }
  if (startsWith(a, "--pcaPop=")) {
    opt$do_pca_pop <- TRUE
    opt$pca_pop_ids_csv <- sub("^--pcaPop=", "", a)
    i <- i + 1L
    next
  }

  if (a == "--pcaInd") {
    opt$pca_ind_mode <- TRUE
    if (i < length(args) && !startsWith(args[i + 1], "-")) {
      opt$pca_ind_samples_csv <- args[i + 1]
      i <- i + 2L
    } else {
      i <- i + 1L
    }
    next
  }
  if (startsWith(a, "--pcaInd=")) {
    opt$pca_ind_mode <- TRUE
    opt$pca_ind_samples_csv <- sub("^--pcaInd=", "", a)
    i <- i + 1L
    next
  }

  if (a == "--pcaK") {
    need_value(i, a)
    opt$pcaK <- args[i + 1]
    i <- i + 2L
    next
  }
  if (a == "--pcaPlotPC") {
    need_value(i, a)
    opt$pcaPlotPC <- args[i + 1]
    i <- i + 2L
    next
  }

  if (a == "--plotH3Ind") {
    need_value(i, a)
    opt$plotH3Ind_csv <- args[i + 1]
    i <- i + 2L
    next
  }
  if (startsWith(a, "--plotH3Ind=")) {
    opt$plotH3Ind_csv <- sub("^--plotH3Ind=", "", a)
    i <- i + 1L
    next
  }

  if (a == "--plotRealByTree") {
    need_value(i, a)
    opt$plotRealByTree_nwk <- args[i + 1]
    i <- i + 2L
    next
  }
  if (startsWith(a, "--plotRealByTree=")) {
    opt$plotRealByTree_nwk <- sub("^--plotRealByTree=", "", a)
    i <- i + 1L
    next
  }

  fail(sprintf("Unknown option: %s", a))
}

if (!nzchar(opt$infile) || !nzchar(opt$mapfile) || !nzchar(opt$prefix)) usage(script_name, code = 1L)
if (!file.exists(opt$infile)) fail(sprintf("infile not found: %s", opt$infile))
if (!file.exists(opt$mapfile)) fail(sprintf("mapfile not found: %s", opt$mapfile))
if (nzchar(opt$plotRealByTree_nwk) && !file.exists(opt$plotRealByTree_nwk)) {
  fail(sprintf("--plotRealByTree file not found: %s", opt$plotRealByTree_nwk))
}

if (!grepl("^[0-9]+$", as.character(opt$minsites))) fail("--minsites must be a non-negative integer.")
opt$minsites <- as.integer(opt$minsites)
if (!grepl("^[-+]?[0-9]*\\.?[0-9]+$", as.character(opt$zthr))) fail("-z must be numeric.")
opt$zthr <- as.numeric(opt$zthr)
if (!grepl("^[0-9]+$", as.character(opt$pcaK)) || as.integer(opt$pcaK) < 2L) fail("--pcaK must be an integer >= 2.")
opt$pcaK <- as.integer(opt$pcaK)
if (!grepl("^[0-9]+,[0-9]+$", opt$pcaPlotPC)) fail("--pcaPlotPC must be like A,B (e.g. 1,2 or 2,3).")

if ((opt$do_votes || opt$do_pca_pop || opt$pca_ind_mode || nzchar(opt$plotRealByTree_nwk)) && !opt$do_real) {
  fail("--votes/--pcaPop/--pcaInd/--plotRealByTree require --triples.")
}

out <- list(
  topo_raw = sprintf("%s.topology.mapped_rows.tsv", opt$prefix),
  topo_pairplot_tsv = sprintf("%s.topology.pairplot.tsv", opt$prefix),
  topo_pairplot_pdf = sprintf("%s.topology.pairplot.pdf", opt$prefix),
  realout = sprintf("%s.real.tsv", opt$prefix),
  realsum = sprintf("%s.real.summary.tsv", opt$prefix),
  votes_out = sprintf("%s.votes.tsv", opt$prefix),
  votes_out_d = sprintf("%s.votes.d.tsv", opt$prefix),
  support_out = sprintf("%s.support_matrix.tsv", opt$prefix),
  support_out_d = sprintf("%s.support_matrix.d.tsv", opt$prefix),
  pca_pop_pdf = sprintf("%s.pca.pop.pdf", opt$prefix),
  pca_pop_tsv = sprintf("%s.pca.pop.tsv", opt$prefix),
  pca_pop_pdf_d = sprintf("%s.pca.pop.d.pdf", opt$prefix),
  pca_pop_tsv_d = sprintf("%s.pca.pop.d.tsv", opt$prefix),
  pca_ind_pdf = sprintf("%s.pca.ind.pdf", opt$prefix),
  pca_ind_tsv = sprintf("%s.pca.ind.tsv", opt$prefix),
  pca_ind_pdf_d = sprintf("%s.pca.ind.d.pdf", opt$prefix),
  pca_ind_tsv_d = sprintf("%s.pca.ind.d.tsv", opt$prefix),
  rankH3Pop_out = sprintf("%s.rankH3Pop.tsv", opt$prefix),
  rankH3Ind_out = sprintf("%s.rankH3Ind.tsv", opt$prefix),
  plotH3Ind_tsv = sprintf("%s.plotH3Ind.tsv", opt$prefix),
  plotH3Ind_pdf = sprintf("%s.plotH3Ind.pdf", opt$prefix),
  realByTree_tsv = sprintf("%s.real.byTree.tsv", opt$prefix),
  realByTree_pdf = sprintf("%s.real.byTree.pdf", opt$prefix),
  plot_raw = sprintf("%s.plotH3Ind.mapped_rows.tsv", opt$prefix)
)

read_map <- function(mapfile) {
  map <- read.table(mapfile, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(map) < 2) fail("map file must have at least 2 columns: SampleID Population")
  map <- map[, 1:2, drop = FALSE]
  names(map) <- c("Sample", "Pop")
  map
}

read_mapped <- function(dstats_file, map, minsites) {
  df <- read.delim(dstats_file, sep = "", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  dcol <- if ("jackEst" %in% names(df)) "jackEst" else if ("Dstat" %in% names(df)) "Dstat" else if ("D" %in% names(df)) "D" else ""
  need <- c("H1", "H2", "H3", "SE", "Z")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0 || !nzchar(dcol)) {
    fail("Missing required columns in header (need H1 H2 H3 SE Z and jackEst/Dstat/D)")
  }

  if (all(c("nABBA", "nBABA") %in% names(df))) {
    df$nSites <- suppressWarnings(as.numeric(df$nABBA) + as.numeric(df$nBABA))
  } else if ("nSites" %in% names(df)) {
    df$nSites <- suppressWarnings(as.numeric(df$nSites))
  } else {
    df$nSites <- 0
  }

  df$D <- suppressWarnings(as.numeric(df[[dcol]]))
  df$SE <- suppressWarnings(as.numeric(df$SE))
  df$Z <- suppressWarnings(as.numeric(df$Z))
  df$nSites[!is.finite(df$nSites)] <- 0

  pop_of <- function(x) map$Pop[match(x, map$Sample)]
  df$H1pop <- pop_of(df$H1)
  df$H2pop <- pop_of(df$H2)
  df$H3pop <- pop_of(df$H3)

  keep <- !is.na(df$H1pop) & !is.na(df$H2pop) & !is.na(df$H3pop) & (df$nSites >= minsites)
  df <- df[keep, c("H1", "H2", "H3", "H1pop", "H2pop", "H3pop", "D", "SE", "Z", "nSites"), drop = FALSE]
  rownames(df) <- NULL
  df
}

map <- read_map(opt$mapfile)
mapped <- read_mapped(opt$infile, map, opt$minsites)

real_df <- NULL
realsum_df <- NULL
support <- NULL
support_d <- NULL

require_pkgs(c("dplyr"))

if (opt$do_real) {
  real_keep <- mapped$H1pop != mapped$H2pop & mapped$H1pop != mapped$H3pop & mapped$H2pop != mapped$H3pop
  real <- mapped[real_keep, , drop = FALSE]

  swap <- real$H1pop > real$H2pop
  H1pop2 <- ifelse(swap, real$H2pop, real$H1pop)
  H2pop2 <- ifelse(swap, real$H1pop, real$H2pop)
  H1_2 <- ifelse(swap, real$H2, real$H1)
  H2_2 <- ifelse(swap, real$H1, real$H2)
  D2 <- ifelse(swap, -real$D, real$D)
  Z2 <- ifelse(swap, -real$Z, real$Z)

  real_df <- data.frame(
    H1pop = H1pop2, H2pop = H2pop2, H3pop = real$H3pop,
    H1 = H1_2, H2 = H2_2, H3 = real$H3,
    D = D2, SE = real$SE, Z = Z2, nSites = real$nSites,
    stringsAsFactors = FALSE
  )
  write.table(real_df, out$realout, sep = "\t", row.names = FALSE, quote = FALSE)

  realsum_df <- dplyr::group_by(real_df, H1pop, H2pop, H3pop)
  realsum_df <- dplyr::summarise(
    realsum_df,
    n = dplyr::n(),
    med_nSites = suppressWarnings(stats::median(nSites, na.rm = TRUE)),
    min_nSites = suppressWarnings(min(nSites, na.rm = TRUE)),
    prop_sig_positive = mean(Z > opt$zthr),
    prop_sig_negative = mean(Z < -opt$zthr),
    prop_sig_total = mean(abs(Z) > opt$zthr),
    mean_Z = mean(Z),
    median_Z = stats::median(Z),
    mean_D = mean(D),
    median_D = stats::median(D),
    .groups = "drop"
  )
  write.table(realsum_df, out$realsum, sep = "\t", row.names = FALSE, quote = FALSE)
}

if (nzchar(opt$plotRealByTree_nwk)) {
  require_pkgs(c("ape", "dplyr", "ggplot2"))
  tr <- ape::read.tree(opt$plotRealByTree_nwk)
  if (is.null(tr)) fail("Could not read tree.")
  if (!ape::is.binary(tr)) fail("Tree must be strictly binary (fully resolved).")
  if (!("Outgroup" %in% tr$tip.label)) fail("Tree must contain a tip named exactly: Outgroup")

  tr <- ape::root(tr, outgroup = "Outgroup", resolve.root = TRUE)
  ntip <- length(tr$tip.label)
  internal_nodes <- (ntip + 1):(ntip + tr$Nnode)

  get_tips <- function(tree, node, ntip) {
    if (node <= ntip) return(tree$tip.label[node])
    ape::extract.clade(tree, node)$tip.label
  }

  tops <- list()
  idx <- 1L
  for (node in internal_nodes) {
    kids <- tr$edge[tr$edge[, 1] == node, 2]
    if (length(kids) != 2) next

    tips1 <- get_tips(tr, kids[1], ntip)
    tips2 <- get_tips(tr, kids[2], ntip)

    if (length(tips1) == 2) {
      sis <- sort(tips1)
      others <- setdiff(setdiff(tr$tip.label, tips1), "Outgroup")
      for (h3 in others) {
        tops[[idx]] <- data.frame(sis1 = sis[1], sis2 = sis[2], H3 = h3, stringsAsFactors = FALSE)
        idx <- idx + 1L
      }
    }
    if (length(tips2) == 2) {
      sis <- sort(tips2)
      others <- setdiff(setdiff(tr$tip.label, tips2), "Outgroup")
      for (h3 in others) {
        tops[[idx]] <- data.frame(sis1 = sis[1], sis2 = sis[2], H3 = h3, stringsAsFactors = FALSE)
        idx <- idx + 1L
      }
    }
  }
  if (length(tops) == 0) fail("Tree produced no (sister pair)|H3 topologies with a 2-tip sister clade.")

  topdf <- dplyr::distinct(dplyr::bind_rows(tops))
  real2 <- dplyr::mutate(real_df, H1c = pmin(H1pop, H2pop), H2c = pmax(H1pop, H2pop))

  m <- dplyr::inner_join(real2, topdf, by = c("H1c" = "sis1", "H2c" = "sis2", "H3pop" = "H3"))
  m <- dplyr::mutate(m, Topology = paste0("(", H1c, ",", H2c, ")|", H3pop), sig = ifelse(abs(Z) > opt$zthr, "sig", "non"))

  if (nrow(m) == 0) fail("No distinct-population triple comparisons matched the supplied tree topology. Check labels.")

  outtab <- dplyr::select(m, Topology, H1pop, H2pop, H3pop, H1, H2, H3, D, SE, Z, nSites, sig)
  write.table(outtab, out$realByTree_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  lev <- dplyr::arrange(dplyr::distinct(outtab, Topology), Topology)$Topology
  outtab$Topology <- factor(outtab$Topology, levels = rev(lev))

  grDevices::pdf(out$realByTree_pdf, width = 10, height = max(6, 0.18 * length(lev) + 2))
  p <- ggplot2::ggplot(outtab, ggplot2::aes(x = D, y = Topology)) +
    ggplot2::geom_jitter(ggplot2::aes(shape = sig), size = 2.2, height = 0.15, width = 0, stroke = 0.9, color = "black", fill = NA) +
    ggplot2::scale_shape_manual(values = c(sig = 19, non = 21)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linetype = "dashed")
    ) +
    ggplot2::labs(
      title = "Real D-stat comparisons consistent with supplied tree",
      subtitle = paste0("Filled = |Z|>", opt$zthr, " ; Open = |Z|<=", opt$zthr, "  |  Label: (sister1,sister2)|H3"),
      x = "D-score", y = ""
    )
  print(p)
  grDevices::dev.off()
}

if (opt$do_topology) {
  require_pkgs(c("dplyr", "ggplot2"))

  write.table(mapped, out$topo_raw, sep = "\t", row.names = FALSE, quote = FALSE)

  blue <- dplyr::transmute(
    dplyr::filter(mapped, H1pop == H2pop, H3pop != H1pop),
    Pop1 = H1pop, Pop2 = H3pop, D = D, Source = "(Pop1,Pop1),Pop2"
  )
  red <- dplyr::transmute(
    dplyr::filter(mapped, H1pop != H2pop, H3pop == H1pop),
    Pop1 = H1pop, Pop2 = H2pop, D = D, Source = "(Pop1,Pop2),Pop1"
  )

  plotdat <- dplyr::bind_rows(blue, red)
  plotdat$Pair <- paste0(plotdat$Pop1, " - ", plotdat$Pop2)
  write.table(plotdat, out$topo_pairplot_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  pair_levels <- dplyr::arrange(dplyr::distinct(plotdat, Pair, Pop1, Pop2), Pop1, Pop2)$Pair
  plotdat$Pair <- factor(plotdat$Pair, levels = rev(pair_levels))

  grDevices::pdf(out$topo_pairplot_pdf, width = 10, height = max(6, 0.18 * length(pair_levels) + 2))
  p <- ggplot2::ggplot(plotdat, ggplot2::aes(x = D, y = Pair)) +
    ggplot2::geom_jitter(ggplot2::aes(color = Source), height = 0.15, width = 0, size = 2.2, shape = 21, fill = NA, stroke = 0.9) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linetype = "dashed")
    ) +
    ggplot2::labs(
      title = "Topology comparison per population pair",
      subtitle = "Blue open circles: (Pop1,Pop1),Pop2   |   Red open circles: (Pop1,Pop2),Pop1",
      x = "D-score", y = ""
    ) +
    ggplot2::scale_color_manual(values = c("(Pop1,Pop1),Pop2" = "blue", "(Pop1,Pop2),Pop1" = "red"))
  print(p)
  grDevices::dev.off()
}

if (opt$do_votes || opt$do_pca_pop || opt$pca_ind_mode) {
  require_pkgs(c("dplyr", "tidyr"))

  tab <- realsum_df
  if (!all(c("H1pop", "H2pop", "H3pop", "mean_Z", "mean_D") %in% names(tab))) fail("real summary missing required columns.")

  pops <- sort(unique(c(tab$H1pop, tab$H2pop, tab$H3pop)))
  n <- length(pops)
  if (n < 4) fail("Need >= 4 populations for trio voting / support matrix.")

  build_votes_support <- function(score_vec, out_votes_file, out_support_file) {
    key <- paste(tab$H1pop, tab$H2pop, tab$H3pop, sep = "|")
    score_map <- stats::setNames(score_vec, key)

    get_score <- function(A, B, C) {
      x1 <- if (A < B) A else B
      x2 <- if (A < B) B else A
      k <- paste(x1, x2, C, sep = "|")
      v <- score_map[[k]]
      if (is.null(v) || length(v) == 0) return(NA_real_)
      as.numeric(v)
    }

    trios <- list(); votes <- list(); ti <- 1L; vi <- 1L
    for (ii in 1:(n - 2)) {
      for (jj in (ii + 1):(n - 1)) {
        for (kk in (jj + 1):n) {
          A <- pops[ii]; Bp <- pops[jj]; C <- pops[kk]
          sAB_C <- get_score(A, Bp, C)
          sAC_B <- get_score(A, C, Bp)
          sBC_A <- get_score(Bp, C, A)
          if (any(is.na(c(sAB_C, sAC_B, sBC_A)))) next

          scores <- c(sAB_C, sAC_B, sBC_A)
          w <- which.min(scores)
          ord <- sort(scores)
          winner <- ord[1]
          runnerup <- ord[2]
          margin <- runnerup - winner

          sister <- if (w == 1) c(A, Bp) else if (w == 2) c(A, C) else c(Bp, C)

          trios[[ti]] <- c(A, Bp, C); ti <- ti + 1L
          votes[[vi]] <- data.frame(
            A = A, B = Bp, C = C,
            best_sister1 = sister[1], best_sister2 = sister[2],
            score_AB_C = sAB_C, score_AC_B = sAC_B, score_BC_A = sBC_A,
            winner = winner, runnerup = runnerup, margin = margin,
            stringsAsFactors = FALSE
          )
          vi <- vi + 1L
        }
      }
    }

    if (length(trios) < 1) fail("No complete trios found (need all 3 resolutions per trio).")

    votes_df <- dplyr::bind_rows(votes)
    write.table(votes_df, out_votes_file, sep = "\t", row.names = FALSE, quote = FALSE)

    wins <- matrix(0, n, n, dimnames = list(pops, pops))
    den <- matrix(0, n, n, dimnames = list(pops, pops))

    for (t in seq_along(trios)) {
      A <- trios[[t]][1]; Bp <- trios[[t]][2]; C <- trios[[t]][3]
      scores <- c(get_score(A, Bp, C), get_score(A, C, Bp), get_score(Bp, C, A))

      for (pp in list(c(A, Bp), c(A, C), c(Bp, C))) {
        den[pp[1], pp[2]] <- den[pp[1], pp[2]] + 1
        den[pp[2], pp[1]] <- den[pp[2], pp[1]] + 1
      }

      w <- which.min(scores)
      sister <- if (w == 1) c(A, Bp) else if (w == 2) c(A, C) else c(Bp, C)
      wins[sister[1], sister[2]] <- wins[sister[1], sister[2]] + 1
      wins[sister[2], sister[1]] <- wins[sister[2], sister[1]] + 1
    }

    support_local <- wins / den
    support_local[is.na(support_local)] <- 0
    diag(support_local) <- 1
    write.table(support_local, out_support_file, sep = "\t", quote = FALSE, col.names = NA)
    support_local
  }

  support <- build_votes_support(abs(tab$mean_Z), out$votes_out, out$support_out)
  support_d <- build_votes_support(abs(tab$mean_D), out$votes_out_d, out$support_out_d)
}

if (opt$do_pca_pop) {
  require_pkgs(c("ggplot2"))
  run_pca_pop <- function(S, out_tsv, out_pdf, metric_label) {
    keep_pops <- parse_csv(opt$pca_pop_ids_csv)
    if (length(keep_pops) > 0) {
      unknown <- setdiff(keep_pops, rownames(S))
      if (length(unknown) > 0) fail(sprintf("Unknown population ID(s) in --pcaPop: %s", paste(unknown, collapse = ", ")))
      S <- S[keep_pops, keep_pops, drop = FALSE]
    }
    if (nrow(S) < 3) fail("pcaPop requires at least 3 populations after filtering.")

    Dm <- 1 - S
    fit <- stats::cmdscale(stats::as.dist(Dm), k = 2, eig = TRUE)
    coords <- data.frame(Pop = rownames(S), Axis1 = fit$points[, 1], Axis2 = fit$points[, 2], stringsAsFactors = FALSE)
    write.table(coords, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

    grDevices::pdf(out_pdf, width = 7, height = 6)
    p <- ggplot2::ggplot(coords, ggplot2::aes(x = Axis1, y = Axis2, label = Pop)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_text(vjust = -0.6, size = 3) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(title = sprintf("Population ordination from trio-vote support (%s-based)", metric_label), x = "Axis 1", y = "Axis 2")
    print(p)
    grDevices::dev.off()
  }

  run_pca_pop(support, out$pca_pop_tsv, out$pca_pop_pdf, "Z")
  run_pca_pop(support_d, out$pca_pop_tsv_d, out$pca_pop_pdf_d, "D")
}

if (opt$pca_ind_mode) {
  require_pkgs(c("dplyr", "tidyr", "ggplot2"))

  pcs <- as.integer(strsplit(opt$pcaPlotPC, ",", fixed = TRUE)[[1]])
  pcA <- pcs[1]
  pcB <- pcs[2]

  keep_samps <- parse_csv(opt$pca_ind_samples_csv)
  if (length(keep_samps) > 0) {
    unknown <- setdiff(keep_samps, unique(map$Sample))
    if (length(unknown) > 0) fail(sprintf("Unknown sample ID(s) in --pcaInd: %s", paste(unknown, collapse = ", ")))
    pops <- sort(unique(map$Pop[match(keep_samps, map$Sample)]))
  } else {
    keep_samps <- map$Sample
    pops <- sort(unique(map$Pop))
  }

  sub <- dplyr::filter(mapped, H3 %in% keep_samps)
  if (nrow(sub) < 1) fail("No D-stat rows left for selected H3 samples after filtering.")

  run_pca_ind <- function(metric_col, out_tsv, out_pdf, metric_label) {
    pair <- dplyr::filter(sub, H1pop != H2pop)
    score_base <- if (metric_col == "Z") pair$Z else pair$D
    aff_long <- dplyr::bind_rows(
      data.frame(H3 = pair$H3, popX = pair$H1pop, score = -score_base, stringsAsFactors = FALSE),
      data.frame(H3 = pair$H3, popX = pair$H2pop, score =  score_base, stringsAsFactors = FALSE)
    )
    aff_long <- dplyr::filter(aff_long, popX %in% pops)

    aff <- dplyr::summarise(dplyr::group_by(aff_long, H3, popX), aff = mean(score), .groups = "drop")
    aff <- tidyr::pivot_wider(aff, names_from = popX, values_from = aff, values_fill = 0)

    for (p in pops) if (!(p %in% names(aff))) aff[[p]] <- 0
    mat <- as.matrix(aff[, pops, drop = FALSE])
    rownames(mat) <- aff$H3

    sds <- apply(mat, 2, stats::sd)
    keep_cols <- which(is.finite(sds) & sds > 0)
    if (length(keep_cols) < 2) fail("Too few non-constant population columns for PCA (need >=2).")

    mat2 <- mat[, keep_cols, drop = FALSE]
    mat_s <- scale(mat2)
    pc <- stats::prcomp(mat_s, center = TRUE, scale. = FALSE)

    K_use <- min(opt$pcaK, ncol(pc$x))
    if (pcA > K_use || pcB > K_use) {
      fail(sprintf("Requested PCs %d,%d but only %d available.", pcA, pcB, K_use))
    }

    coords <- data.frame(Sample = rownames(mat_s), stringsAsFactors = FALSE)
    for (k in seq_len(K_use)) coords[[paste0("PC", k)]] <- pc$x[, k]
    coords$Pop <- map$Pop[match(coords$Sample, map$Sample)]
    write.table(coords, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

    grDevices::pdf(out_pdf, width = 7, height = 6)
    p <- ggplot2::ggplot(coords, ggplot2::aes(x = .data[[paste0("PC", pcA)]], y = .data[[paste0("PC", pcB)]], color = Pop)) +
      ggplot2::geom_point(size = 2.4) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::labs(
        title = sprintf("Individuals PCA from H3-only affinity (mean signed %s vs populations)", metric_label),
        x = paste0("PC", pcA), y = paste0("PC", pcB)
      )
    print(p)
    grDevices::dev.off()
  }

  run_pca_ind("Z", out$pca_ind_tsv, out$pca_ind_pdf, "Z")
  run_pca_ind("D", out$pca_ind_tsv_d, out$pca_ind_pdf_d, "D")
}

if (nzchar(opt$rankH3Pop) || nzchar(opt$rankH3Ind)) {
  require_pkgs(c("dplyr"))

  rank_from_H3 <- function(subdf) {
    out <- dplyr::transmute(dplyr::filter(subdf, H1pop != H2pop), winner = ifelse(Z > 0, H2pop, H1pop), w = abs(Z))
    out <- dplyr::summarise(dplyr::group_by(out, winner), wins_wt = sum(w), n_wins = dplyr::n(), .groups = "drop")
    out <- dplyr::rename(out, Candidate = winner)
    dplyr::arrange(out, dplyr::desc(wins_wt), dplyr::desc(n_wins))
  }

  if (nzchar(opt$rankH3Pop)) {
    outtab <- rank_from_H3(dplyr::filter(mapped, H3pop == opt$rankH3Pop))
    write.table(outtab, out$rankH3Pop_out, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  if (nzchar(opt$rankH3Ind)) {
    outtab <- rank_from_H3(dplyr::filter(mapped, H3 == opt$rankH3Ind))
    write.table(outtab, out$rankH3Ind_out, sep = "\t", row.names = FALSE, quote = FALSE)
  }
}

if (nzchar(opt$plotH3Ind_csv)) {
  require_pkgs(c("dplyr", "ggplot2"))
  h3list <- parse_csv(opt$plotH3Ind_csv)

  write.table(mapped, out$plot_raw, sep = "\t", row.names = FALSE, quote = FALSE)

  sub <- mapped
  sub <- dplyr::filter(sub, H3 %in% h3list)
  sub <- dplyr::filter(sub, H1pop != H2pop)
  sub <- dplyr::filter(sub, !(H1 %in% h3list), !(H2 %in% h3list))
  sub <- dplyr::mutate(sub, H1H2 = paste0(H1pop, "_", H2pop), sig = ifelse(abs(Z) > opt$zthr, "sig", "non"))

  outtab <- dplyr::select(sub, H3, H3pop, H1H2, D, Z, H1, H2, H1pop, H2pop, nSites, sig)
  write.table(outtab, out$plotH3Ind_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  outtab$H3 <- factor(outtab$H3, levels = rev(unique(outtab$H3)))
  grDevices::pdf(out$plotH3Ind_pdf, width = 11, height = max(6, 0.25 * length(levels(outtab$H3)) + 2))
  p <- ggplot2::ggplot(outtab, ggplot2::aes(x = D, y = H3)) +
    ggplot2::geom_jitter(ggplot2::aes(color = H1H2, shape = sig), size = 2.4, height = 0.15, width = 0) +
    ggplot2::scale_shape_manual(values = c(sig = 19, non = 21)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(linetype = "dashed"),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = "D-score for selected H3 individuals",
      subtitle = paste0("H3 in: ", paste(h3list, collapse = ", "), " | filled=|Z|>", opt$zthr, ", open=|Z|<=", opt$zthr),
      x = "D-score", y = "H3 individual", color = "H1pop_H2pop"
    )
  print(p)
  grDevices::dev.off()
}

cat("Done.\n")
cat("Wrote (depending on options):\n")
if (opt$do_real) {
  cat(" ", out$realout, "\n", sep = "")
  cat(" ", out$realsum, "\n", sep = "")
}
if (nzchar(opt$plotRealByTree_nwk)) {
  cat(" ", out$realByTree_tsv, "\n", sep = "")
  cat(" ", out$realByTree_pdf, "\n", sep = "")
}
if (opt$do_topology) {
  cat(" ", out$topo_raw, "\n", sep = "")
  cat(" ", out$topo_pairplot_tsv, "\n", sep = "")
  cat(" ", out$topo_pairplot_pdf, "\n", sep = "")
}
if (opt$do_votes || opt$do_pca_pop || opt$pca_ind_mode) {
  cat(" ", out$votes_out, "\n", sep = "")
  cat(" ", out$votes_out_d, "\n", sep = "")
  cat(" ", out$support_out, "\n", sep = "")
  cat(" ", out$support_out_d, "\n", sep = "")
}
if (opt$do_pca_pop) {
  cat(" ", out$pca_pop_pdf, "\n", sep = "")
  cat(" ", out$pca_pop_tsv, "\n", sep = "")
  cat(" ", out$pca_pop_pdf_d, "\n", sep = "")
  cat(" ", out$pca_pop_tsv_d, "\n", sep = "")
}
if (opt$pca_ind_mode) {
  cat(" ", out$pca_ind_pdf, "\n", sep = "")
  cat(" ", out$pca_ind_tsv, "\n", sep = "")
  cat(" ", out$pca_ind_pdf_d, "\n", sep = "")
  cat(" ", out$pca_ind_tsv_d, "\n", sep = "")
}
if (nzchar(opt$rankH3Pop)) cat(" ", out$rankH3Pop_out, "\n", sep = "")
if (nzchar(opt$rankH3Ind)) cat(" ", out$rankH3Ind_out, "\n", sep = "")
if (nzchar(opt$plotH3Ind_csv)) {
  cat(" ", out$plotH3Ind_pdf, "\n", sep = "")
  cat(" ", out$plotH3Ind_tsv, "\n", sep = "")
}
