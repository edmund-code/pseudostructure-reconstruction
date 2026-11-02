# Fig. 3e Slingshot baseline, called by fig3_baselines.py (protocol: docs/paper/fig3-validation-experiments.md, 3e).
# Usage: Rscript --vanilla fig3_slingshot.R <folder> <approx_points>
# Reads <folder>/embedding.f64 (units x dims, row-major little-endian float64), dims.txt (n_units, n_dims),
# units.txt, clusters.txt and start.txt, written by fig3_baselines.py.
# Writes <folder>/slingshot.csv (unit, averaged pseudotime, per-lineage pseudotime and curve weight),
# lineages.txt (the cluster path of each lineage) and session.txt.
# Needs an R with slingshot and DelayedMatrixStats, which getLineages loads (conda-forge + bioconda: r-base 4.4,
# bioconductor-slingshot, bioconductor-delayedmatrixstats). Never install into kidney-pseudospace.
suppressPackageStartupMessages(library(slingshot))
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
approx_points <- as.integer(args[2])
dims <- as.integer(readLines(file.path(folder, "dims.txt")))
units <- readLines(file.path(folder, "units.txt"))
clusters <- readLines(file.path(folder, "clusters.txt"))
start <- readLines(file.path(folder, "start.txt"))
con <- file(file.path(folder, "embedding.f64"), "rb")
X <- matrix(readBin(con, "double", n = dims[1] * dims[2], size = 8, endian = "little"), nrow = dims[1], byrow = TRUE)
close(con)
stopifnot(nrow(X) == length(units), length(clusters) == length(units), start %in% clusters, all(is.finite(X)))
rownames(X) <- units

elapsed <- system.time(
  sds <- slingshot(X, clusterLabels = clusters, start.clus = start, approx_points = approx_points)
)[["elapsed"]]
pt <- slingPseudotime(sds, na = FALSE)
w <- slingCurveWeights(sds, as.probs = TRUE)
avg <- slingAvgPseudotime(sds)  # rowSums(pseudotime * curve weights as probabilities)
out <- data.frame(unit = units, slingshot = avg, check.names = FALSE)
for (j in seq_len(ncol(pt))) {
  out[[paste0("pt_", colnames(pt)[j])]] <- pt[, j]
  out[[paste0("w_", colnames(w)[j])]] <- w[, j]
}
write.csv(out, file.path(folder, "slingshot.csv"), row.names = FALSE)
lin <- slingLineages(sds)
writeLines(vapply(names(lin), function(n) paste0(n, ": ", paste(lin[[n]], collapse = " > ")), ""),
           file.path(folder, "lineages.txt"))
writeLines(c(sprintf("slingshot_seconds %.1f", elapsed), sprintf("approx_points %d", approx_points),
             paste("slingshot", as.character(packageVersion("slingshot"))),
             paste("princurve", as.character(packageVersion("princurve"))),
             paste("TrajectoryUtils", as.character(packageVersion("TrajectoryUtils"))),
             R.version.string), file.path(folder, "session.txt"))
