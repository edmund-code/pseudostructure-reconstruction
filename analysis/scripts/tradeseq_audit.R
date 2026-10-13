# Notebook 87, section A2: tradeSeq run on the units, genes and partitions of notebook 71.
# Usage: Rscript tradeseq_audit.R <folder> <workers>
# Reads <folder>/counts.mtx (genes x units), genes.txt, units.csv, subset.txt and grid.txt, written by notebook 87.
# Writes condition_<arm>.csv, association_mouse.csv, gam_subset.csv, curves_subset.csv, offsets.csv and session.txt.
# Run from the tradeseq-audit environment (analysis/scripts/tradeseq_audit_environment.yml) with one BLAS thread per
# worker (OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1).
suppressPackageStartupMessages({
  library(tradeSeq)
  library(mgcv)
  library(Matrix)
  library(BiocParallel)
  library(SingleCellExperiment)
})
args <- commandArgs(trailingOnly = TRUE)
folder <- args[1]
workers <- as.integer(args[2])
set.seed(87)
counts <- as(Matrix::readMM(file.path(folder, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(folder, "genes.txt"))
units <- read.csv(file.path(folder, "units.csv"), stringsAsFactors = FALSE)
subset <- readLines(file.path(folder, "subset.txt"))
grid <- as.numeric(readLines(file.path(folder, "grid.txt")))
rownames(counts) <- genes
colnames(counts) <- units$unit
counts <- as.matrix(counts)
stopifnot(ncol(counts) == nrow(units), all(subset %in% genes))

within_contrasts <- function(specimen, group) {
  # Within-group sum-to-zero specimen contrasts (notebook 71's coding), one column per group with two specimens.
  out <- NULL
  for (g in sort(unique(group))) {
    names <- sort(unique(specimen[group == g]))
    for (s in names[-length(names)]) {
      out <- cbind(out, (specimen == s) - (specimen == names[length(names)]))
    }
  }
  out
}

offset_all <- tradeSeq:::.get_offset(NULL, counts)            # tradeSeq's default: log(library x TMM factor)
write.csv(data.frame(unit = units$unit, offset_tmm = offset_all, library = colSums(counts)),
          file.path(folder, "offsets.csv"), row.names = FALSE)

CHUNK <- 1500

fit_chunked <- function(counts_in, test, ...) {
  # fitGAM on gene chunks (bounded memory per forked worker); a chunk that loses a worker is retried once with half
  # the workers. Returns the rbind of test(sce) over chunks, with tradeSeq's convergence flags.
  ids <- split(seq_len(nrow(counts_in)), ceiling(seq_len(nrow(counts_in)) / CHUNK))
  out <- lapply(ids, function(id) {
    run <- function(w) {
      sce <- fitGAM(counts_in, genes = id, ..., verbose = FALSE, parallel = TRUE, BPPARAM = MulticoreParam(workers = w))
      res <- test(sce)
      res$gene <- rownames(res)
      res$converged <- rowData(sce)$tradeSeq$converged
      res
    }
    tryCatch(run(workers), error = function(e) {
      message("retrying a chunk with ", workers %/% 2, " workers: ", conditionMessage(e))
      run(workers %/% 2)
    })
  })
  do.call(rbind, out)
}

fit_and_test <- function(arm, condition, use_u) {
  cond <- factor(condition)
  U <- NULL
  if (use_u) U <- cbind(1, within_contrasts(units$specimen, condition))
  t0 <- Sys.time()
  res <- fit_chunked(counts, function(sce) conditionTest(sce, global = TRUE, l2fc = 0),
                     pseudotime = matrix(units$position), cellWeights = matrix(1, ncol(counts), 1),
                     conditions = cond, U = U, nknots = 6, offset = offset_all)
  write.csv(res, file.path(folder, paste0("condition_", arm, ".csv")), row.names = FALSE)
  message(arm, ": ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

arms <- list(c("species_default", "species", "FALSE"), c("relabel1_default", "relabel1", "FALSE"),
             c("relabel2_default", "relabel2", "FALSE"), c("species_U", "species", "TRUE"),
             c("relabel1_U", "relabel1", "TRUE"), c("relabel2_U", "relabel2", "TRUE"))
for (a in arms) fit_and_test(a[1], units[[a[2]]], as.logical(a[3]))

# Within mouse: no conditions, tradeSeq's default offset on the mouse units.
mouse <- units$species == "mouse"
cm <- counts[, mouse]
t0 <- Sys.time()
assoc <- fit_chunked(cm, function(sce) associationTest(sce, global = TRUE),
                     pseudotime = matrix(units$position[mouse]), cellWeights = matrix(1, sum(mouse), 1), nknots = 6)
write.csv(assoc, file.path(folder, "association_mouse.csv"), row.names = FALSE)
message("mouse association: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")

# Gene subset: the species fit as a list of gam objects, for the penalty (edf), theta, AIC and the curves.
gams <- fitGAM(counts[subset, , drop = FALSE], pseudotime = matrix(units$position),
               cellWeights = matrix(1, ncol(counts), 1), conditions = factor(units$species), nknots = 6,
               offset = offset_all, verbose = FALSE, sce = FALSE, parallel = TRUE,
               BPPARAM = MulticoreParam(workers = workers))
levels_cond <- levels(factor(units$species))
rows <- list()
curves <- list()
for (g in names(gams)) {
  m <- gams[[g]]
  if (inherits(m, "try-error")) next
  row <- data.frame(gene = g, theta = m$family$getTheta(TRUE), aic = m$aic, sp = m$sp[1], converged = m$converged)
  b <- coef(m)
  for (k in seq_along(m$smooth)) {
    sm <- m$smooth[[k]]
    idx <- sm$first.para:sm$last.para
    lev <- levels_cond[k]
    row[[paste0("edf_", lev)]] <- sum(m$edf[idx])
    nd <- data.frame(t1 = grid)
    nd[[sm$by]] <- 1
    curves[[paste(g, lev)]] <- data.frame(gene = g, condition = lev, point = seq_along(grid),
                                          value = as.numeric(b[1] + PredictMat(sm, nd) %*% b[idx]))
  }
  rows[[g]] <- row
}
write.csv(do.call(rbind, rows), file.path(folder, "gam_subset.csv"), row.names = FALSE)
write.csv(do.call(rbind, curves), file.path(folder, "curves_subset.csv"), row.names = FALSE)
writeLines(c(capture.output(sessionInfo()), paste("tradeSeq", packageVersion("tradeSeq")),
             paste("mgcv", packageVersion("mgcv")), paste("edgeR", packageVersion("edgeR"))),
           file.path(folder, "session.txt"))
