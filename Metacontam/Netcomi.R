# Load necessary library
library(NetCoMi)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
sparseMethod_threshold <- as.numeric(args[3])

# Argument 4 = BH-FDR alpha. Absent or <= 0 means no FDR, which is the original
# v0.0.1 behaviour. The CLI passes 0.1 by default. The adaptive prevalence filter
# (utils.find_best_threshold) calls this script with three arguments only, so no
# FDR is applied there -- that keeps the taxon set passing the filter independent
# of the edge rule, so the effect of r and the effect of the FDR do not mix.
fdr_alpha <- if (length(args) >= 4 && nzchar(args[4])) as.numeric(args[4]) else NA_real_

# Argument 5 = zero handling before the CLR. Absent keeps the published setting.
#   multRepl          NetCoMi's default. With no dl given it injects dl = 1e-3 and
#                     zCompositions::multRepl fills zeros with 0.65*dl.
#   multRepl:<dl>     set the detection limit directly (proportion of the composition).
#   bayesMult         Bayesian-multiplicative; uses no detection limit.
#   pseudo:<count>    add a pseudocount at the count level (e.g. 0.5).
zero_spec <- if (length(args) >= 5 && nzchar(args[5])) args[5] else "multRepl"
zs        <- strsplit(zero_spec, ":", fixed = TRUE)[[1]]
zero_method <- zs[1]
zero_value  <- if (length(zs) >= 2) as.numeric(zs[2]) else NA_real_

zero_par <- if (zero_method == "multRepl") {
  if (is.na(zero_value)) list(z.delete = FALSE)                       # published
  else list(z.delete = FALSE, z.warning = 1, dl = zero_value)
} else if (zero_method == "bayesMult") {
  list(z.delete = FALSE, z.warning = 1)
} else if (zero_method == "pseudo") {
  list(pseudocount = if (is.na(zero_value)) 0.5 else zero_value)
} else {
  stop(sprintf("unknown zero method: %s", zero_method))
}
cat(sprintf("[zero] method=%s%s\n", zero_method,
            if (is.na(zero_value)) "" else sprintf(" value=%g", zero_value)))

# Load and preprocess the data
df <- t(read.csv(input_file, sep="\t", row.names = 1))

# Construct the network
net_pears <- netConstruct(df,
                          measure = "pearson",
                          normMethod = "clr",
                          zeroMethod = zero_method,
                          zeroPar = zero_par,
                          sparsMethod = "threshold",
                          dataType = "counts",
                          thresh = sparseMethod_threshold,
                          verbose = 3)

edge_list <- net_pears$edgelist1

# -- BH-FDR ------------------------------------------------------------------
# A fixed correlation threshold is more or less stringent depending on the sample
# size, so Benjamini-Hochberg is applied on top of it, to the edges that already
# passed the r threshold.
#
#   * P values come from a two-sided t test on |r| (df = n-2). NetCoMi thresholds
#     on |r| as well, so negative edges reach the graph and every edge is tested
#     regardless of sign. Testing only the positive edges would halve i in the BH
#     bound alpha*i/m and make it far stricter.
#   * m is the number of all taxon pairs, not the number of edges that passed.
#     The edges kept by the r threshold carry the smallest P values, so sorting
#     only those still gives the correct BH result.
#   * The original row order is restored afterwards. community_detection takes
#     np.unique over v1/v2, so order does not affect it, but keeping the original
#     output format makes runs easier to diff.
if (!is.na(fdr_alpha) && fdr_alpha > 0) {
  n_samples <- nrow(df)
  n_taxa    <- ncol(df)
  m         <- n_taxa * (n_taxa - 1) / 2

  r  <- pmin(abs(edge_list$asso), 0.999999)
  tt <- r * sqrt((n_samples - 2) / (1 - r^2))
  pv <- 2 * pt(tt, df = n_samples - 2, lower.tail = FALSE)

  ord  <- order(pv)
  pass <- pv[ord] <= fdr_alpha * seq_along(ord) / m
  k    <- if (any(pass)) max(which(pass)) else 0L

  keep <- rep(FALSE, length(pv))
  if (k > 0) keep[ord[seq_len(k)]] <- TRUE

  cat(sprintf(
    "[FDR] alpha=%.3g  samples=%d  taxa=%d  tests=%.0f  edges %d -> %d (removed %d)\n",
    fdr_alpha, n_samples, n_taxa, m, length(pv), k, length(pv) - k))

  edge_list <- edge_list[keep, , drop = FALSE]
}

# Save the edge information to the output file
write.table(edge_list, file = output_file, sep="\t")
