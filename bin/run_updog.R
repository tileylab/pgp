#!/usr/bin/env Rscript
# Genotype one ploidy group with updog::multidog and write a tidy table of
# per-SNP per-individual reference-allele dosage calls.
#
# Usage: run_updog.R <ref.txt> <size.txt> <ploidy> <model> <prefix> [ncores] [extra]
#   ref.txt / size.txt : TSV, rows = SNPs (row names), cols = individuals (header)
#   extra              : optional R named-argument string forwarded to multidog,
#                        e.g. "seq = 0.01, update_bias = FALSE" (from ext.args)
#   output             : <prefix>.updog.tsv  (columns: snp, ind, geno)
# updog's 'geno' is the reference-allele dosage (0..ploidy).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: run_updog.R <ref.txt> <size.txt> <ploidy> <model> <prefix> [ncores] [extra]")
}
ref_file  <- args[1]
size_file <- args[2]
ploidy    <- as.integer(args[3])
model     <- args[4]
prefix    <- args[5]
ncores    <- if (length(args) >= 6) as.integer(args[6]) else 1L

# Extra multidog options from the module's ext.args (trusted pipeline config).
extra <- list()
if (length(args) >= 7 && nzchar(trimws(args[7]))) {
  extra <- eval(parse(text = paste0("list(", args[7], ")")))
}

suppressMessages(library(updog))

read_mat <- function(path) {
  as.matrix(read.table(path, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE))
}

refmat  <- read_mat(ref_file)
sizemat <- read_mat(size_file)

mout <- do.call(multidog, c(
  list(refmat = refmat, sizemat = sizemat, ploidy = ploidy, model = model, nc = ncores),
  extra
))

ind <- mout$inddf[, c("snp", "ind", "geno")]
write.table(
  ind,
  file      = paste0(prefix, ".updog.tsv"),
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
