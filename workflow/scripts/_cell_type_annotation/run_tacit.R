# TACIT cell type annotation of one feature table (quantify_cells output, any segmentation method).
#
# The signature markers are normalized with arcsinh(x / cofactor) followed by a z-score per marker
# (shifted to min 0, see below), then TACIT (TACIT_function.R, next to this script) labels every cell. Writes one row per cell:
# cell_id, cell_type ("Others" when unassigned), the first two of TACIT's p UMAP dimensions, the cell
# types whose CTR passes its threshold and whether TACIT treated the cell as mixed. Also writes the
# per-cell CTR scores and the inferred thresholds (all in the normalized space below), and for the mixed
# cells the cell type fractions among their k nearest clean neighbours (TACIT's kNN label is the top one).
#
# Signature CSV: first column cell_type, then one column per marker (weights in [0, 1], empty = 0).
# Negative weights mark markers a cell type should NOT express (e.g. CD3E for Fibroblast): they lower
# its CTR score in cells positive for that marker. Not part of TACIT's design; they only count as
# "not a marker" in TACIT's subset rules (weight > 0). TACIT_function.R fits each CTR threshold on the
# positive-weight score and applies it to the penalized score (written to --out_ctr), so negative
# weights can only remove positive calls; see LOCAL CHANGE there.
# The min-0 shift below means a cell at background level of an excluded marker gets ~no penalty.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

options_list <- list(
  make_option("--feature_table", help = "feature_table.csv of quantify_cells"),
  make_option("--signature", help = "cell type x marker signature CSV"),
  make_option("--cofactor", type = "double", help = "arcsinh cofactor"),
  make_option("--r", type = "double", help = "TACIT clustering resolution"),
  make_option("--p", type = "integer", help = "TACIT number of UMAP dimensions"),
  make_option("--threads", type = "integer", default = 1, help = "cores TACIT may use"),
  make_option("--out_csv", help = "output CSV: cell_id, cell_type, UMAP1, UMAP2, positive_cell_types, mixed"),
  make_option("--out_ctr", help = "output CSV: cell_id, then the CTR score of every cell type"),
  make_option("--out_thresholds", help = "output CSV: kind (cell_type = CTR threshold, marker = subset marker threshold), name, threshold"),
  make_option("--out_knn", help = "output CSV, one row per mixed cell: cell_id, knn_k, fraction of the k nearest clean neighbours of every cell type")
)
args <- parse_args(OptionParser(option_list = options_list))

# TACIT runs availableCores() - 1 workers/threads everywhere: report --threads + 1 cores so that it
# uses exactly --threads (the main process only waits while the workers run)
n_cores <- max(args$threads, 1) + 1
options(
  parallelly.availableCores.methods = "custom",
  parallelly.availableCores.custom = function() n_cores,
  # ScaleData runs under plan(multisession) and ships the whole matrix to each worker
  future.globals.maxSize = Inf
)
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "TACIT_function.R"))

# read inputs, keeping marker and cell type names as they are (no make.names)
features <- fread(args$feature_table, check.names = FALSE, data.table = FALSE)
signature <- fread(args$signature, check.names = FALSE, data.table = FALSE)

# check the signature: cell_type column, at least one cell type, markers present in the table
if (colnames(signature)[1] != "cell_type") stop("the first signature column must be named cell_type")
if (nrow(signature) == 0) stop(paste("no cell types in", args$signature, "- fill in the signature first"))
markers <- colnames(signature)[-1]
missing <- setdiff(markers, colnames(features))
if (length(missing) > 0) stop(paste("signature markers missing from the feature table:", paste(missing, collapse = ", ")))
if (anyDuplicated(signature$cell_type)) stop("duplicated cell types in the signature")
cat(sprintf("%d cells, %d cell types, %d markers\n", nrow(features), nrow(signature), length(markers)))

# arcsinh + z-score of every signature marker
expr <- asinh(as.matrix(features[, markers, drop = FALSE]) / args$cofactor)
expr <- scale(expr)
# a constant marker has sd 0 and becomes NaN: set it to 0
expr[is.nan(expr)] <- 0
# shift each marker to min 0: TACIT's microclusters use Seurat's CLR, which breaks on negative
# values (synthetic test: 70% vs 90% accuracy); its scores, thresholds and kNN are shift-invariant
expr <- sweep(expr, 2, apply(expr, 2, min))

# TACIT looks columns up by name (get(), formulas, data.frame()), so markers and cell types get
# syntactic internal names (m1.., ct1..), mapped back after the run
marker_ids <- paste0("m", seq_along(markers))
ct_ids <- paste0("ct", seq_len(nrow(signature)))
colnames(expr) <- marker_ids
expr <- as.data.frame(expr)
sig <- data.frame(cell_type = ct_ids, signature[, markers, drop = FALSE], check.names = FALSE)
colnames(sig)[-1] <- marker_ids
# empty weights mean "not a marker of this cell type"
sig[is.na(sig)] <- 0

# run TACIT, all arguments named (the README's positional call does not match the signature)
res <- TACIT(data_expression = expr, r = args$r, p = args$p, Signature = sig)
# one label per input cell, otherwise data.frame() below would silently recycle cell_id
if (nrow(res) != nrow(features)) stop(sprintf("TACIT returned %d rows for %d cells", nrow(res), nrow(features)))

# labels back to the signature's cell type names; "Others" stays
labels <- res$TACIT
labels <- ifelse(labels %in% ct_ids, signature$cell_type[match(labels, ct_ids)], labels)

# cell types whose CTR passes its threshold in each cell (after TACIT's subset rules), "::"-joined
positive <- as.matrix(attr(res, "positive")) == 1
colnames(positive) <- signature$cell_type[match(colnames(positive), ct_ids)]
positive_cell_types <- apply(positive, 1, function(x) paste(colnames(positive)[x], collapse = "::"))
# mixed: cell of a mixed microcluster positive for several cell types, labelled by TACIT's kNN
mixed <- seq_len(nrow(features)) %in% attr(res, "mixed_ids")

# one row per cell, in feature table order (TACIT keeps the input order)
out <- data.frame(cell_id = features$cell_id, cell_type = labels, UMAP1 = res$UMAP1, UMAP2 = res$UMAP2,
                  positive_cell_types = positive_cell_types, mixed = mixed)
fwrite(out, args$out_csv)
cat("written:", args$out_csv, "\n")
print(sort(table(out$cell_type), decreasing = TRUE))
cat(sprintf("%d mixed cells (%.1f%%)\n", sum(mixed), 100 * mean(mixed)))

# CTR scores (normalized markers %*% signature weights), one column per cell type
ctr <- as.data.frame(attr(res, "ctr"))
colnames(ctr) <- signature$cell_type[match(colnames(ctr), ct_ids)]
fwrite(data.frame(cell_id = features$cell_id, ctr, check.names = FALSE), args$out_ctr)
cat("written:", args$out_ctr, "\n")

# CTR threshold per cell type, plus the marker thresholds TACIT fits for subset markers (if any)
ctr_thr <- attr(res, "ctr_thresholds")
thresholds <- data.frame(kind = "cell_type", name = signature$cell_type[match(ctr_thr$Name, ct_ids)],
                         threshold = ctr_thr$value)
marker_thr <- attr(res, "marker_thresholds")
if (!is.null(marker_thr)) {
  thresholds <- rbind(thresholds, data.frame(kind = "marker", name = markers[match(marker_thr$Name, marker_ids)],
                                             threshold = marker_thr$value))
}
fwrite(thresholds, args$out_thresholds)
cat("written:", args$out_thresholds, "\n")

# kNN of the mixed cells: fraction of the k nearest clean neighbours per cell type (NA = not in the cell's
# mix of positive cell types, all NA = not resolved by kNN), cell types back to the signature names
knn <- attr(res, "knn_fractions")
frac_cols <- setdiff(colnames(knn), c("ID", "knn_k"))
colnames(knn)[match(frac_cols, colnames(knn))] <- signature$cell_type[match(frac_cols, ct_ids)]
fwrite(data.frame(cell_id = features$cell_id[knn$ID], knn[, -1], check.names = FALSE), args$out_knn)
cat("written:", args$out_knn, "\n")
