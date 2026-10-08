# Export raw RNA counts + metadata from a VIPerturb-seq Seurat object, losslessly.
# Usage: Rscript export_seurat.R <in.RDS> <out_dir>
# Writes: i.int32, p.int32, x.int32 (dgCMatrix slots, genes x cells, CSC), genes.txt,
# barcodes.txt, meta.csv, object.json (structure + validation statistics).
suppressPackageStartupMessages({library(SeuratObject); library(Matrix); library(jsonlite); library(digest)})
args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
obj <- readRDS(inp)
assays <- Assays(obj)
layers <- lapply(assays, function(a) tryCatch(Layers(obj[[a]]), error = function(e) character(0)))
names(layers) <- assays
other_dims <- lapply(setdiff(assays, "RNA"), function(a) dim(obj[[a]]))
cls <- class(obj)[1]; dflt <- DefaultAssay(obj)
md <- obj[[]]
m <- LayerData(obj, assay = "RNA", layer = "counts")
rm(obj); invisible(gc())
if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
stopifnot(all(m@x == round(m@x)), all(m@x >= 0))
con <- file(file.path(out, "i.int32"), "wb"); writeBin(as.integer(m@i), con, size = 4); close(con)
con <- file(file.path(out, "p.int32"), "wb"); writeBin(as.integer(m@p), con, size = 4); close(con)
con <- file(file.path(out, "x.int32"), "wb"); writeBin(as.integer(m@x), con, size = 4); close(con)
writeLines(rownames(m), file.path(out, "genes.txt"))
writeLines(colnames(m), file.path(out, "barcodes.txt"))
write.csv(md, file.path(out, "meta.csv"), row.names = TRUE)
set.seed(20261006)
rc <- sample(ncol(m), min(2000, ncol(m))); rg <- sample(nrow(m), min(500, nrow(m)))
re <- cbind(sample(nrow(m), 5000, replace = TRUE), sample(ncol(m), 5000, replace = TRUE))
info <- list(
  input = basename(inp), md5_input = digest(file = inp, algo = "md5"),
  r_version = R.version.string, SeuratObject = as.character(packageVersion("SeuratObject")),
  Matrix = as.character(packageVersion("Matrix")), class = cls,
  assays = assays, layers = layers, default_assay = dflt,
  dims_genes_cells = dim(m), nnz = length(m@x), total_counts = sum(m@x),
  meta_columns = colnames(md), n_meta_rows = nrow(md),
  check_cols = rc, check_col_sums = as.numeric(colSums(m[, rc])),
  check_rows = rg, check_row_sums = as.numeric(rowSums(m[rg, ])),
  check_entries = re, check_entry_values = as.numeric(m[re]),
  other_assay_dims = other_dims
)
write(toJSON(info, auto_unbox = TRUE, digits = NA), file.path(out, "object.json"))
cat("exported", inp, "\n")
