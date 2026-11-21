# Transport to the unmodified official CytoTRACE 0.3.1 source.
args <- commandArgs(trailingOnly=TRUE)
folder <- args[1]
.libPaths(args[3], include.site=FALSE)
required <- c("Matrix", "HiClimR", "ccaPP", "nnls")
for (pkg in required) loadNamespace(pkg)
source(args[2])
mat <- as.matrix(Matrix::readMM(file.path(folder, "counts.mtx")))
rownames(mat) <- readLines(file.path(folder, "genes.txt"))
colnames(mat) <- readLines(file.path(folder, "cells.txt"))
stopifnot(all(is.finite(mat)), all(mat >= 0), all(mat == floor(mat)),
          all(colSums(mat) == scan(file.path(folder, "cell_sums.txt"), quiet=TRUE)))
# The official function interprets max(mat)<50 as log2 input. Full-library CPM
# preserves detection counts and avoids that heuristic for low integer counts.
input_unit <- "raw_integer_counts"
if (max(mat) < 50) {
  mat <- sweep(mat, 2, colSums(mat), "/") * 1e6
  input_unit <- "CPM_from_full_raw_library_to_avoid_official_log2_input_heuristic"
}
set.seed(as.integer(args[4]))
res <- CytoTRACE(mat, batch=NULL, enableFast=as.logical(args[5]),
                 ncores=as.integer(args[6]), subsamplesize=as.integer(args[7]))
stopifnot(identical(names(res$CytoTRACE), colnames(mat)))
score <- as.numeric(res$CytoTRACE)
qc <- is.finite(score) & score >= 0 & score <= 1
write.table(data.frame(cell=colnames(mat), cytotrace_score=score,
                       cytotrace_qc=qc,
                       cytotrace_detected_genes=as.numeric(res$Counts)),
            file.path(folder, "scores.tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="NA")
writeLines(c(paste("input_unit", input_unit), capture.output(sessionInfo()),
             paste(required, vapply(required, function(p) as.character(packageVersion(p)), ""))),
           file.path(folder, "runtime.txt"))
