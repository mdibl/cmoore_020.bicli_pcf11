# ============================================================
# DEXSeq-based Alternative Polyadenylation (APA) Analysis
# Generic script — configure the USER CONFIGURATION section only.
#
# Usage:
#   Set BASE_DIR, input file paths, condition labels, and
#   GROUPING_VAR (NULL for a simple two-condition design; set to
#   a colData column name to run one analysis per level of that
#   variable, e.g., GROUPING_VAR <- "timepoint").
#
# Input:
#   - PAS count matrix (rows = PAS IDs, columns = samples)
#   - PAS annotation file (PolyA_DB v3.2 or v4.1 auto-detected)
#   - Design file (tab-delimited; rows = samples, columns = metadata)
#   - (optional) sample name list
# ============================================================

# -- Libraries ---------------------------------------------------------------
library(DEXSeq)
library(SummarizedExperiment)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(scales)
library(limma)
library(dplyr)
library(tidyr)
library(tibble)
library(matrixStats)
library(stringr)

# ============================================================
#  USER CONFIGURATION — EDIT THESE
# ============================================================

BASE_DIR <- "/compbio/analysis/ClaireMoore/cmoore_020.bicli_pcf11/"   # project root

# Input files (relative to BASE_DIR)
COUNTS_CSV  <- "data/cluster.all.reads.csv"
SAMPLES_TXT <- "data/sample_name.txt"       # set to NULL to use all count_mat columns
ANNO_FILE   <- "data/PolyA_DB_v4.1/hg38.PAS.main.tsv"
DESIGN_FILE <- "data/design.txt"

# Output base directory
OUT_BASE <- file.path(BASE_DIR, "output/DEXseq/")

# Experimental design
CTRL_LABEL  <- "Control"      # reference condition (fold-change denominator)
TRTMT_LABEL <- "Experimental"    # treatment/experimental condition

# Grouping variable:
#   NULL  — run a single analysis on all samples (simple two-condition design)
#   "timepoint" — run one analysis per level of colData$timepoint
GROUPING_VAR    <- NULL
GROUPING_LEVELS <- NULL       # NULL = auto-detect from colData (alphabetical order)

# PAS type filter — keep rows whose PAS_type matches this regex
PAS_TYPE_REGEX <- "^3'UTR"

# Analysis thresholds
MIN_TOTAL_READS   <- 10    # min "this"-read row-sum per group to keep a PAS
MIN_PER_CONDITION <- 1     # min samples with count > 0 per condition; 0 = disable
PADJ_CUT          <- 0.05  # FDR threshold for significance calls
LFC_CUT           <- 0.5   # |log2FC| threshold for "big" hits in gene summary

# PCA settings
NTOP_PCA      <- 2000        # top-N variable PAS for PCA (by row variance)
PCA_LABEL_COL <- "replicate" # colData column to use as point labels (NULL = none)

# Genes to plot PSI bar charts; set to character(0) to skip
GENES_OF_INTEREST <- c("PCF11", "TAB2", "ICAM1")

# RUVseq batch correction
#   USE_RUV <- TRUE  when replicates cluster by batch in PCA rather than condition.
#   RUVs estimates W factors from within-condition replicate variation; W_1:exon is
#   then added to the DEXSeq design to correct for batch effects on relative PAS usage.
#   RUV_K: number of factors to estimate (1 is almost always sufficient; raise to 2
#   only if the first factor does not account for the outlier clustering).
USE_RUV <- TRUE
RUV_K   <- 1

# ============================================================
#  DERIVED PATHS — do not edit below this line
# ============================================================

path_counts  <- file.path(BASE_DIR, COUNTS_CSV)
path_anno    <- file.path(BASE_DIR, ANNO_FILE)
path_design  <- file.path(BASE_DIR, DESIGN_FILE)
path_samples <- if (!is.null(SAMPLES_TXT)) file.path(BASE_DIR, SAMPLES_TXT) else NULL

dir_results     <- OUT_BASE
dir_plots       <- file.path(OUT_BASE, "plots")
dir_usage_plots <- file.path(OUT_BASE, "plots/usage")
for (d in c(dir_results, dir_plots, dir_usage_plots)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

CONDITION_LEVELS <- c(CTRL_LABEL, TRTMT_LABEL)
FC_COL           <- paste0("log2fold_", TRTMT_LABEL, "_", CTRL_LABEL)
FC_EXPORT        <- paste0("log2fold_", TRTMT_LABEL, "_v_", CTRL_LABEL)
USAGE_CTRL_COL   <- paste0("meanUsage_", CTRL_LABEL)
USAGE_TRTMT_COL  <- paste0("meanUsage_", TRTMT_LABEL)
SE_CTRL_COL      <- paste0("seUsage_", CTRL_LABEL)
SE_TRTMT_COL     <- paste0("seUsage_", TRTMT_LABEL)

# ============================================================
#  BUILD COUNT MATRIX
# ============================================================

message("Reading count matrix...")
counts_df <- read.csv(path_counts, check.names = FALSE)
stopifnot("hit_PAS_ID" %in% names(counts_df))

if (!is.null(path_samples)) {
  samples <- scan(path_samples, what = character(), quiet = TRUE)
} else {
  samples <- setdiff(names(counts_df), "hit_PAS_ID")
}
stopifnot(all(samples %in% names(counts_df)))

count_mat_raw <- as.matrix(counts_df[, samples])            # unfiltered — used for library sizes
rownames(count_mat_raw) <- counts_df$hit_PAS_ID
storage.mode(count_mat_raw) <- "integer"

count_mat <- count_mat_raw   # will be row-filtered after annotation alignment below

# ============================================================
#  BUILD PAS ANNOTATION
# ============================================================

message("Reading PAS annotation...")
pas_anno <- read.delim(path_anno, stringsAsFactors = FALSE, check.names = FALSE)

# Auto-detect PolyA_DB format:
#   v3.2: explicit Chromosome / Position / Strand columns
#   v4.1: all three encoded in PAS_ID as "chr:strand:pos"
has_coord_cols <- all(c("Chromosome", "Position", "Strand") %in% names(pas_anno))

if (!has_coord_cols) {
  if ("PAS_ID" %in% names(pas_anno)) {
    message("PolyA_DB v4.1 format detected: parsing coordinates from PAS_ID.")
    parsed              <- strsplit(as.character(pas_anno$PAS_ID), ":", fixed = TRUE)
    pas_anno$Chromosome <- vapply(parsed, `[`, character(1), 1)
    pas_anno$Strand     <- vapply(parsed, `[`, character(1), 2)
    pas_anno$Position   <- as.integer(vapply(parsed, `[`, character(1), 3))
  } else {
    message("No coordinate columns found — GRanges will not be built; APA direction will be NA.")
  }
}

pas_anno <- pas_anno[grepl(PAS_TYPE_REGEX, pas_anno$PAS_type), ]

blank_to_na <- function(x) {
  x <- stringr::str_squish(x)
  x[tolower(x) %in% c("", "na")] <- NA_character_
  x
}

id_cols <- intersect(c("Gene Symbol", "GeneSymbol", "Ensembl ID", "RefSeq Gene ID", "FANTOM ID"),
                     names(pas_anno))
if (length(id_cols) == 0) stop("No recognized gene-ID columns in annotation file.")
for (col in id_cols) pas_anno[[col]] <- blank_to_na(pas_anno[[col]])
pas_anno$groupID <- do.call(dplyr::coalesce, as.list(pas_anno[, id_cols, drop = FALSE]))

# Rescue PAS IDs that are in the count matrix but have no gene name
unmatched_ids <- setdiff(counts_df$hit_PAS_ID, pas_anno$PAS_ID[!is.na(pas_anno$groupID)])
rescue <- pas_anno$PAS_ID %in% unmatched_ids & is.na(pas_anno$groupID)
pas_anno$groupID[rescue] <- pas_anno$PAS_ID[rescue]
pas_anno <- pas_anno[!is.na(pas_anno$groupID), ]

if (anyDuplicated(pas_anno[c("PAS_ID", "groupID")])) {
  pas_anno <- pas_anno[!duplicated(pas_anno[c("PAS_ID", "groupID")]), ]
}

# ============================================================
#  ALIGN MATRIX & ANNOTATION; SET groupID / featureID
# ============================================================

common_pas <- intersect(rownames(count_mat), pas_anno$PAS_ID)
count_mat  <- count_mat[common_pas, , drop = FALSE]
pas_anno   <- pas_anno[match(common_pas, pas_anno$PAS_ID), ]

groupID   <- pas_anno$groupID
# PolyA_DB v4.1 PAS_IDs contain ':' (chr:strand:pos). DEXSeq strips ':' from
# featureIDs, which would desync GRanges names from internal rownames and crash.
# Sanitize here so DEXSeq sees IDs it doesn't need to modify.
featureID           <- gsub("[: ]", "_", pas_anno$PAS_ID)
rownames(count_mat) <- featureID          # count_mat rows are in pas_anno order here
pas_anno$featureID_safe <- featureID      # retained for annotation merge below

# ============================================================
#  BUILD GRanges (enables APA direction labelling)
# ============================================================

if (all(c("Chromosome", "Position", "Strand") %in% names(pas_anno))) {
  gr <- GRanges(
    seqnames = pas_anno$Chromosome,
    ranges   = IRanges(start = as.integer(pas_anno$Position), width = 1L),
    strand   = pas_anno$Strand
  )
  # DEXSeq composite key — allows gr[rrn] to subset correctly after per-group filtering
  names(gr) <- paste(groupID, featureID, sep = ":")
} else {
  gr <- NULL
}

# ============================================================
#  LOAD DESIGN / colData
# ============================================================

colData <- read.delim(path_design, row.names = 1, check.names = FALSE)
if (!"condition" %in% names(colData)) stop("Design file must have a 'condition' column.")
colData$condition <- factor(colData$condition, levels = CONDITION_LEVELS)
if (!is.null(GROUPING_VAR) && GROUPING_VAR %in% names(colData)) {
  colData[[GROUPING_VAR]] <- factor(colData[[GROUPING_VAR]])
}
stopifnot(all(colnames(count_mat) == rownames(colData)))

# ============================================================
#  RUVseq BATCH CORRECTION (optional — set USE_RUV <- TRUE)
# ============================================================

if (USE_RUV) {
  if (!requireNamespace("RUVSeq", quietly = TRUE))
    stop("RUVSeq is required when USE_RUV = TRUE.  Install with: BiocManager::install('RUVSeq')")
  library(RUVSeq)

  message("Estimating unwanted variation with RUVs (k = ", RUV_K, ")...")

  # scIdx: each row is a set of within-condition replicates — any variation among
  # them is treated as unwanted.  Pad uneven groups with -1.
  ctrl_idx <- which(colData$condition == CTRL_LABEL)
  trt_idx  <- which(colData$condition == TRTMT_LABEL)
  nc       <- max(length(ctrl_idx), length(trt_idx))
  scIdx    <- matrix(-1L, nrow = 2L, ncol = nc)
  scIdx[1L, seq_along(ctrl_idx)] <- ctrl_idx
  scIdx[2L, seq_along(trt_idx)]  <- trt_idx

  ruv_set <- RUVSeq::newSeqExpressionSet(
    counts    = count_mat,
    phenoData = data.frame(condition = colData$condition, row.names = rownames(colData))
  )
  ruv_fit <- RUVSeq::RUVs(ruv_set,
                           cIdx  = seq_len(nrow(count_mat)),
                           k     = RUV_K,
                           scIdx = scIdx)

  # Attach W factors to colData so run_dexseq_group picks them up automatically
  w_cols <- paste0("W_", seq_len(RUV_K))
  for (wc in w_cols) {
    colData[[wc]] <- pData(ruv_fit)[rownames(colData), wc]
  }
  message("  W factors attached to colData: ", paste(w_cols, collapse = ", "))
  message("  W_1 per sample:")
  for (s in rownames(colData))
    message(sprintf("    %-25s  %+.4f  [%s]", s, colData[s, "W_1"], colData[s, "condition"]))
}

# ============================================================
#  GROUPING SETUP
# ============================================================

if (is.null(GROUPING_VAR)) {
  group_list   <- list("all" = rownames(colData))
  group_suffix <- c("all" = "")
} else {
  if (is.null(GROUPING_LEVELS)) {
    GROUPING_LEVELS <- sort(unique(as.character(colData[[GROUPING_VAR]])))
  }
  group_list <- setNames(
    lapply(GROUPING_LEVELS, function(lv) {
      rownames(colData)[as.character(colData[[GROUPING_VAR]]) == lv]
    }),
    GROUPING_LEVELS
  )
  group_suffix <- setNames(
    vapply(GROUPING_LEVELS, function(lv) paste0(".", sub("day$", "d", lv)), character(1)),
    GROUPING_LEVELS
  )
}

# ============================================================
#  FUNCTION: global PCA on raw VST PAS counts (pre-DEXSeq)
# ============================================================

make_pca_global <- function(count_mat, col_data,
                             color_col    = "condition",
                             shape_col    = GROUPING_VAR,
                             label_col    = PCA_LABEL_COL,
                             ntop         = NTOP_PCA,
                             batch_covars = NULL,
                             outfile      = NULL) {
  dds     <- DESeq2::DESeqDataSetFromMatrix(countData = count_mat,
                                             colData  = col_data,
                                             design   = ~ 1)
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))

  # Optionally remove batch covariates from VST matrix for visualization only.
  # This does not affect the statistical model — it is purely for the PCA plot.
  if (!is.null(batch_covars) && length(batch_covars) > 0) {
    stopifnot(all(batch_covars %in% colnames(col_data)))
    vst_mat <- limma::removeBatchEffect(
      vst_mat,
      covariates = as.matrix(col_data[, batch_covars, drop = FALSE])
    )
  }

  if (!is.null(ntop) && ntop > 0 && nrow(vst_mat) > ntop) {
    rv  <- matrixStats::rowVars(vst_mat)
    vst_mat <- vst_mat[order(rv, decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
  }

  p   <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pct <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1)

  pcs <- as.data.frame(p$x[, 1:2]) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(tibble::rownames_to_column(as.data.frame(col_data), "sample_id"),
                     by = "sample_id")

  aes_base <- if (!is.null(color_col) && color_col %in% names(pcs))
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]])
  else
    ggplot2::aes(x = PC1, y = PC2)

  pca_title <- if (!is.null(batch_covars))
    "Global PCA — VST 3'UTR PAS counts (RUV-corrected, visual only)"
  else
    "Global PCA — VST 3'UTR PAS counts"

  g <- ggplot2::ggplot(pcs, aes_base) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                  y = paste0("PC2 (", pct[2], "%)"),
                  title = pca_title) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(shape_col) && shape_col %in% names(pcs)) {
    g <- g + ggplot2::aes(shape = .data[[shape_col]])
  }
  if (!is.null(label_col) && label_col %in% names(pcs)) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_col]]), size = 3)
  }

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 7, height = 5, dpi = 300)
  invisible(list(plot = g, pcs = pcs, prcomp = p))
}

# ============================================================
#  FUNCTION: library size bar chart
# ============================================================

make_library_sizes <- function(count_mat_raw, col_data, outfile = NULL) {
  df <- data.frame(
    sample    = colnames(count_mat_raw),
    lib_size  = colSums(count_mat_raw),
    condition = as.character(col_data[colnames(count_mat_raw), "condition"]),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$condition, df$sample), ]
  df$sample <- factor(df$sample, levels = df$sample)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = lib_size / 1e6, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Sample", y = "Total PAS reads (M)",
                  title = "Library sizes (total PAS read counts)") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(outfile)) {
    w <- max(4, 0.5 * nrow(df) + 2)
    ggplot2::ggsave(outfile, g, width = w, height = 4, dpi = 300)
  }
  invisible(list(plot = g, data = df))
}

# ============================================================
#  FUNCTION: per-group PCA using DEXSeqDataSet (post-normalization)
# ============================================================

make_pca_dxd <- function(dxd, grp_label = "",
                          color_col = "condition",
                          label_col = PCA_LABEL_COL,
                          ntop = NTOP_PCA,
                          outfile = NULL) {
  cd_full <- as.data.frame(SummarizedExperiment::colData(dxd))
  this_idx <- which(cd_full$exon == "this")
  cd  <- cd_full[this_idx, , drop = FALSE]
  mat <- SummarizedExperiment::assay(dxd)[, this_idx, drop = FALSE]

  dds     <- DESeq2::DESeqDataSetFromMatrix(countData = mat, colData = cd, design = ~ 1)
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))

  if (!is.null(ntop) && ntop > 0 && nrow(vst_mat) > ntop) {
    rv  <- matrixStats::rowVars(vst_mat)
    vst_mat <- vst_mat[order(rv, decreasing = TRUE)[seq_len(ntop)], , drop = FALSE]
  }

  p   <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
  pct <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1)

  pcs <- as.data.frame(p$x[, 1:2]) |>
    tibble::rownames_to_column("sample_id") |>
    dplyr::left_join(tibble::rownames_to_column(cd, "sample_id"), by = "sample_id")

  aes_base <- if (!is.null(color_col) && color_col %in% names(pcs))
    ggplot2::aes(x = PC1, y = PC2, color = .data[[color_col]])
  else
    ggplot2::aes(x = PC1, y = PC2)

  ttl <- paste0("PCA — normalized PAS counts",
                if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")

  g <- ggplot2::ggplot(pcs, aes_base) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                  y = paste0("PC2 (", pct[2], "%)"),
                  title = ttl) +
    ggplot2::theme_minimal(base_size = 12)

  if (!is.null(label_col) && label_col %in% names(pcs)) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_col]]), size = 3)
  }

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 7, height = 5, dpi = 300)
  invisible(list(plot = g, pcs = pcs))
}

# ============================================================
#  FUNCTION: size factor bar chart
# ============================================================

make_size_factors <- function(dxd, grp_label = "", outfile = NULL) {
  cd_full  <- as.data.frame(SummarizedExperiment::colData(dxd))
  this_idx <- which(cd_full$exon == "this")
  cd       <- cd_full[this_idx, , drop = FALSE]
  sf       <- sizeFactors(dxd)[this_idx]

  df <- data.frame(
    sample      = as.character(cd$sample),
    size_factor = as.numeric(sf),
    condition   = as.character(cd$condition),
    stringsAsFactors = FALSE
  )
  df <- df[order(df$condition, df$sample), ]
  df$sample <- factor(df$sample, levels = df$sample)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = size_factor, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    ggplot2::labs(
      x     = "Sample",
      y     = "Size factor",
      title = paste0("DEXSeq size factors",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(outfile)) {
    w <- max(4, 0.5 * nrow(df) + 2)
    ggplot2::ggsave(outfile, g, width = w, height = 4, dpi = 300)
  }
  invisible(list(plot = g, data = df))
}

# ============================================================
#  FUNCTION: dispersion plot (base R via DEXSeq)
# ============================================================

make_dispersion_plot <- function(dxd, grp_label = "", outfile = NULL) {
  ttl <- paste0("Dispersion estimates",
                if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
  if (!is.null(outfile)) png(outfile, width = 800, height = 600, res = 100)
  DESeq2::plotDispEsts(dxd, main = ttl)
  if (!is.null(outfile)) { dev.off(); invisible(NULL) }
}

# ============================================================
#  FUNCTION: p-value histogram
# ============================================================

make_pvalue_hist <- function(res_df, grp_label = "", outfile = NULL) {
  pv <- res_df$pvalue[!is.na(res_df$pvalue)]
  df <- data.frame(pvalue = pv)
  n_bins <- 20
  expected_height <- length(pv) / n_bins

  g <- ggplot2::ggplot(df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = n_bins, fill = "steelblue", color = "white", boundary = 0) +
    ggplot2::geom_hline(yintercept = expected_height, linetype = "dashed",
                        color = "firebrick", linewidth = 0.8) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    ggplot2::labs(
      x     = "p-value",
      y     = "Count",
      title = paste0("p-value distribution",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else ""),
      caption = "Dashed red line = uniform expectation"
    ) +
    ggplot2::theme_minimal(base_size = 11)

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 5, height = 4, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: MA plot
# ============================================================

make_ma_plot <- function(res_df, grp_label = "",
                          padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                          fc_col = FC_COL, outfile = NULL) {
  if (!all(c("exonBaseMean", fc_col, "padj") %in% names(res_df))) {
    warning("MA plot skipped — required columns missing.")
    return(invisible(NULL))
  }

  df <- res_df %>%
    dplyr::filter(!is.na(exonBaseMean), !is.na(.data[[fc_col]]), exonBaseMean > 0) %>%
    dplyr::mutate(
      logMean = log10(exonBaseMean),
      sig     = !is.na(padj) & padj <= padj_cut
    )

  n_up   <- sum(df$sig & df[[fc_col]] > 0, na.rm = TRUE)
  n_down <- sum(df$sig & df[[fc_col]] < 0, na.rm = TRUE)

  g <- ggplot2::ggplot(df, ggplot2::aes(x = logMean, y = .data[[fc_col]], color = sig)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::scale_color_manual(values = c("FALSE" = "gray60", "TRUE" = "firebrick"),
                                labels = c("Not sig", paste0("padj<", padj_cut)),
                                name = NULL) +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid",  color = "black",   linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::annotate("text", x = -Inf, y = Inf,  hjust = -0.1, vjust = 1.5,
                      label = paste0("Up: ", n_up), color = "firebrick", size = 3.5) +
    ggplot2::annotate("text", x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.5,
                      label = paste0("Down: ", n_down), color = "steelblue", size = 3.5) +
    ggplot2::labs(
      x     = expression(log[10](mean~PAS~count)),
      y     = expression(log[2]~fold~change),
      title = paste0("MA plot",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 6, height = 5, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: Volcano plot
# ============================================================

make_volcano_plot <- function(res_df, grp_label = "",
                               padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                               fc_col = FC_COL, n_label = 10, outfile = NULL) {
  if (!all(c(fc_col, "padj", "gene") %in% names(res_df))) {
    warning("Volcano plot skipped — required columns missing.")
    return(invisible(NULL))
  }

  df <- res_df %>%
    dplyr::filter(!is.na(.data[[fc_col]]), !is.na(padj)) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      category = dplyr::case_when(
        padj <= padj_cut & abs(.data[[fc_col]]) >= lfc_cut & .data[[fc_col]] > 0 ~ "Up (sig + |L2FC|)",
        padj <= padj_cut & abs(.data[[fc_col]]) >= lfc_cut & .data[[fc_col]] < 0 ~ "Down (sig + |L2FC|)",
        padj <= padj_cut ~ "Sig (small FC)",
        TRUE ~ "Not sig"
      )
    )

  top_labs <- df %>%
    dplyr::filter(padj <= padj_cut) %>%
    dplyr::slice_min(padj, n = n_label, with_ties = FALSE)

  cols <- c("Up (sig + |L2FC|)"   = "firebrick",
             "Down (sig + |L2FC|)" = "steelblue",
             "Sig (small FC)"      = "goldenrod3",
             "Not sig"             = "gray70")

  g <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[fc_col]],
                                         y = neg_log10_padj,
                                         color = category)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.6) +
    ggplot2::scale_color_manual(values = cols, name = NULL) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut),
                        linetype = "dashed", color = "gray40", linewidth = 0.4) +
    ggrepel::geom_text_repel(data = top_labs,
                              ggplot2::aes(label = gene),
                              size = 2.5, max.overlaps = 20) +
    ggplot2::labs(
      x     = expression(log[2]~fold~change),
      y     = expression(-log[10](padj)),
      title = paste0("Volcano plot",
                     if (nchar(grp_label) > 0) paste0(" (", grp_label, ")") else "")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(outfile)) ggplot2::ggsave(outfile, g, width = 6, height = 6, dpi = 300)
  invisible(g)
}

# ============================================================
#  FUNCTION: run DEXSeq for one group of samples
# ============================================================

run_dexseq_group <- function(grp_label, sub_samples, count_mat, colData,
                              featureID, groupID,
                              gr = NULL, pas_anno = NULL,
                              min_total = 10, min_per_condition = 0,
                              use_ruv = FALSE, ruv_k = 1) {

  message(sprintf("[%s] Subsetting samples...", grp_label))
  sub_cd  <- droplevels(colData[sub_samples, , drop = FALSE])
  sub_cnt <- count_mat[, sub_samples, drop = FALSE]
  stopifnot(all(colnames(sub_cnt) == rownames(sub_cd)))

  sub_cd$sample    <- factor(rownames(sub_cd))
  sub_cd$condition <- droplevels(factor(sub_cd$condition, levels = CONDITION_LEVELS))
  char_cols <- vapply(sub_cd, is.character, logical(1))
  sub_cd[char_cols] <- lapply(sub_cd[char_cols], factor)

  if (length(levels(sub_cd$condition)) < 2 || any(table(sub_cd$condition) == 0)) {
    stop(sprintf("Group '%s' is missing one or both conditions.", grp_label))
  }

  # W_k:exon corrects for batch effects on *relative* PAS usage (PSI); it is not
  # collinear with the `sample` term, which only absorbs gene-level expression totals.
  dex_design <- if (use_ruv && ruv_k >= 1) {
    w_terms <- paste(paste0("W_", seq_len(ruv_k)), "exon", sep = ":")
    as.formula(paste("~ sample + exon +", paste(c(w_terms, "condition:exon"), collapse = " + ")))
  } else {
    ~ sample + exon + condition:exon
  }
  message(sprintf("[%s] Creating DEXSeqDataSet (design: %s)...",
                  grp_label, deparse(dex_design)))
  dxd <- DEXSeqDataSet(
    countData     = sub_cnt,
    sampleData    = sub_cd,
    design        = dex_design,
    featureID     = featureID,
    groupID       = groupID,
    featureRanges = gr
  )

  # Filter on "this" columns only — the doubled matrix inflates apparent counts
  message(sprintf("[%s] Filtering low-count PAS...", grp_label))
  this_idx <- which(colData(dxd)$exon == "this")
  mat_this <- counts(dxd)[, this_idx, drop = FALSE]
  dxd <- dxd[rowSums(mat_this) >= min_total, , drop = FALSE]

  if (min_per_condition > 0) {
    this_idx  <- which(colData(dxd)$exon == "this")
    mat_this  <- counts(dxd)[, this_idx, drop = FALSE]
    cond_this <- colData(dxd)$condition[this_idx]
    keep2 <- rowSums(mat_this[, cond_this == CTRL_LABEL,  drop = FALSE] > 0) >= min_per_condition &
              rowSums(mat_this[, cond_this == TRTMT_LABEL, drop = FALSE] > 0) >= min_per_condition
    dxd <- dxd[keep2, , drop = FALSE]
  }

  # Drop genes with < 2 PAS remaining (no within-gene contrast)
  tab <- table(rowData(dxd)$groupID)
  dxd <- dxd[rowData(dxd)$groupID %in% names(tab[tab >= 2]), ]

  message(sprintf("[%s] Normalizing and testing (%d PAS, %d genes)...",
                  grp_label, nrow(dxd), length(unique(rowData(dxd)$groupID))))
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)
  dxr <- DEXSeqResults(dxd)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition",
                                  denominator = CTRL_LABEL)

  message(sprintf("[%s] Building results table...", grp_label))
  res     <- S4Vectors::as.data.frame(dxr)
  rd      <- S4Vectors::mcols(dxd)
  fc_cols <- grep("^log2fold_", colnames(rd), value = TRUE)
  if (length(fc_cols) == 0) stop("No fold-change columns found in rowData(dxd).")

  fc_df <- base::as.data.frame(
    lapply(S4Vectors::as.list(rd[, fc_cols, drop = FALSE]), as.numeric),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  key       <- rownames(dxd)
  res$key   <- rownames(res)
  fc_df$key <- key
  res_tp    <- merge(res, fc_df, by = "key", all.x = TRUE, sort = FALSE)
  rownames(res_tp) <- res_tp$key
  res_tp$key       <- NULL

  res_tp$gene    <- rowData(dxd)$groupID[match(rownames(res_tp), rownames(dxd))]
  res_tp$feature <- rowData(dxd)$featureID[match(rownames(res_tp), rownames(dxd))]

  if (!is.null(pas_anno)) {
    # Join on featureID_safe (sanitized); bring original PAS_ID through for readability
    ann_cols <- intersect(c("featureID_safe", "PAS_ID", "Intron_exon_location", "PAS_type"),
                          names(pas_anno))
    if (length(ann_cols) > 1) {
      ann <- pas_anno[, ann_cols, drop = FALSE]
      names(ann)[names(ann) == "featureID_safe"] <- "featureID"
      ann <- ann[!duplicated(ann$featureID), , drop = FALSE]
      ann$featureID    <- as.character(ann$featureID)
      res_tp$featureID <- as.character(res_tp$featureID)
      res_tp <- merge(res_tp, ann, by = "featureID", all.x = TRUE, sort = FALSE)
    }
  }

  # APA direction: distal PAS = furthest 3' per gene
  #   + strand: highest coordinate; - strand: lowest coordinate
  # Lengthened = increased distal usage; Shortened = increased proximal usage
  message(sprintf("[%s] Assigning APA direction...", grp_label))
  res_tp$APA_direction <- NA_character_

  if (!is.null(gr)) {
    rrn <- rownames(dxd)
    if (!is.null(names(gr)) && all(rrn %in% names(gr))) {
      pos      <- as.integer(GenomicRanges::start(gr[rrn]))
      strand   <- as.character(GenomicRanges::strand(gr[rrn]))
      pos_meta <- data.frame(row = rrn, pos = pos, strand = strand,
                             stringsAsFactors = FALSE)
      res_tp$row <- rownames(res_tp)
      res_tp     <- dplyr::left_join(res_tp, pos_meta, by = "row")

      res_tp <- res_tp |>
        dplyr::group_by(gene) |>
        dplyr::mutate(
          g_max     = suppressWarnings(max(pos, na.rm = TRUE)),
          g_min     = suppressWarnings(min(pos, na.rm = TRUE)),
          g_max     = ifelse(is.finite(g_max), g_max, NA_real_),
          g_min     = ifelse(is.finite(g_min), g_min, NA_real_),
          is_distal = dplyr::case_when(
            strand == "+" & !is.na(pos) & !is.na(g_max) & pos == g_max ~ TRUE,
            strand == "-" & !is.na(pos) & !is.na(g_min) & pos == g_min ~ TRUE,
            TRUE ~ FALSE
          )
        ) |>
        dplyr::ungroup()

      fc_col_dir <- grep("^log2fold_", names(res_tp), value = TRUE)[1]
      if (!is.na(fc_col_dir)) {
        res_tp <- dplyr::mutate(res_tp,
          APA_direction = dplyr::case_when(
            is_distal  & .data[[fc_col_dir]] > 0  ~ "Lengthened",
            is_distal  & .data[[fc_col_dir]] < 0  ~ "Shortened",
            !is_distal & .data[[fc_col_dir]] > 0  ~ "Shortened",
            !is_distal & .data[[fc_col_dir]] < 0  ~ "Lengthened",
            TRUE ~ "Ambiguous"
          )
        )
        bad <- is.na(res_tp[[fc_col_dir]]) | is.na(res_tp$is_distal)
        res_tp$APA_direction[bad] <- "Ambiguous"
      }
      res_tp <- res_tp[, setdiff(names(res_tp), c("row", "g_max", "g_min")), drop = FALSE]
    } else {
      warning(sprintf("[%s] GRanges names don't cover all dxd features — APA direction skipped.", grp_label))
    }
  }

  fc_out_cols <- grep("^log2fold_", colnames(res_tp), value = TRUE)
  count_cols  <- grep("^countData\\.", colnames(res_tp), value = TRUE)
  ctrl_cols   <- grep(paste0("^countData\\.", CTRL_LABEL),  count_cols, value = TRUE)
  treat_cols  <- grep(paste0("^countData\\.", TRTMT_LABEL), count_cols, value = TRUE)

  desired_order <- c(
    "groupID", "featureID", "PAS_ID", "gene", "feature",
    "Intron_exon_location", "PAS_type", "APA_direction",
    "exonBaseMean", "dispersion", "stat", "pvalue", "padj",
    fc_out_cols, ctrl_cols, treat_cols,
    "genomicData.seqnames", "genomicData.start", "genomicData.end",
    "genomicData.width", "genomicData.strand"
  )
  desired_order <- desired_order[desired_order %in% colnames(res_tp)]
  res_tp <- res_tp[, desired_order, drop = FALSE]

  list(dxd = dxd, dxr = dxr, results = res_tp)
}

# ============================================================
#  FUNCTION: collapse per-PAS results to gene-level summary
# ============================================================

collapse_gene <- function(res_df, gene_q = NULL,
                           padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
                           fc_col = FC_COL) {
  safe_min  <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x)  else NA_real_ }
  safe_max  <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x)  else NA_real_ }
  safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  has_dir   <- "APA_direction" %in% names(res_df)

  res_df %>%
    dplyr::mutate(
      padj_clean = ifelse(is.infinite(padj), NA_real_, padj),
      lfc_clean  = ifelse(is.infinite(.data[[fc_col]]), NA_real_, .data[[fc_col]])
    ) %>%
    dplyr::mutate(
      sig = !is.na(padj_clean) & padj_clean <= padj_cut,
      big = sig & abs(lfc_clean) >= lfc_cut
    ) %>%
    dplyr::group_by(groupID) %>%
    dplyr::summarise(
      n_PAS        = dplyr::n(),
      n_sig        = sum(sig, na.rm = TRUE),
      n_big        = sum(big, na.rm = TRUE),
      min_padj     = safe_min(padj_clean),
      max_absL2FC  = safe_max(abs(lfc_clean[sig])),
      mean_absL2FC = safe_mean(abs(lfc_clean[sig])),
      dir_consensus = if (has_dir) {
        d <- APA_direction[sig & !is.na(APA_direction)]
        if (!length(d))                  NA_character_
        else if (all(d == "Lengthened")) "Lengthened_only"
        else if (all(d == "Shortened"))  "Shortened_only"
        else                             "Mixed"
      } else NA_character_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(perGeneQ = if (!is.null(gene_q)) unname(gene_q[groupID]) else NA_real_) %>%
    dplyr::arrange(is.na(perGeneQ), perGeneQ, min_padj, dplyr::desc(max_absL2FC))
}

# ============================================================
#  FUNCTION: per-condition mean PSI + SEM
# ============================================================

usage_by_condition <- function(res_df, sub_cd) {
  cnt_cols <- grep("^countData\\.", names(res_df), value = TRUE)
  counts   <- as.matrix(res_df[, cnt_cols, drop = FALSE])
  colnames(counts) <- sub("^countData\\.", "", colnames(counts))
  stopifnot(all(colnames(counts) %in% rownames(sub_cd)))

  gene_totals <- rowsum(counts, res_df$groupID)
  totals_row  <- gene_totals[res_df$groupID, colnames(counts), drop = FALSE]
  usage       <- counts / pmax(totals_row, 1)

  cond     <- sub_cd[colnames(counts), "condition", drop = TRUE]
  ctrl_idx <- which(cond == CTRL_LABEL)
  trt_idx  <- which(cond == TRTMT_LABEL)

  se_fn <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) > 1) sd(x) / sqrt(length(x)) else NA_real_
  }

  out <- res_df
  out[["meanUsage_All"]] <- rowMeans(usage, na.rm = TRUE)
  out[[USAGE_CTRL_COL]]  <- if (length(ctrl_idx)) rowMeans(usage[, ctrl_idx, drop = FALSE]) else NA_real_
  out[[USAGE_TRTMT_COL]] <- if (length(trt_idx))  rowMeans(usage[, trt_idx,  drop = FALSE]) else NA_real_
  out[[SE_CTRL_COL]]     <- if (length(ctrl_idx) > 1) apply(usage[, ctrl_idx, drop = FALSE], 1, se_fn) else rep(NA_real_, nrow(usage))
  out[[SE_TRTMT_COL]]    <- if (length(trt_idx)  > 1) apply(usage[, trt_idx,  drop = FALSE], 1, se_fn) else rep(NA_real_, nrow(usage))
  out
}

# ============================================================
#  FUNCTION: PSI bar chart with ΔPSI labels, SEM error bars,
#            significance markers
# ============================================================

plot_gene_usage <- function(res_u, gene_symbol, title_suffix = NULL,
                             ctrl_col  = USAGE_CTRL_COL,
                             trt_col   = USAGE_TRTMT_COL,
                             se_ctrl   = SE_CTRL_COL,
                             se_trt    = SE_TRTMT_COL,
                             fc_col    = FC_COL) {
  stopifnot(all(c("gene", ctrl_col, trt_col) %in% names(res_u)))

  df <- res_u %>%
    dplyr::filter(gene == gene_symbol) %>%
    dplyr::mutate(
      PAS_raw = dplyr::coalesce(.data$feature, .data$featureID),
      pos     = dplyr::coalesce(.data$genomicData.start, NA_integer_),
      strand  = as.character(dplyr::coalesce(.data$genomicData.strand, NA))
    ) %>%
    dplyr::filter(!(is.na(.data[[ctrl_col]]) & is.na(.data[[trt_col]])))

  if (nrow(df) == 0) stop(sprintf("No rows for gene '%s'.", gene_symbol))

  if (!all(is.na(df$pos)) && !all(is.na(df$strand))) {
    df <- df %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(is_distal = dplyr::case_when(
        strand == "+" & pos == max(pos, na.rm = TRUE) ~ TRUE,
        strand == "-" & pos == min(pos, na.rm = TRUE) ~ TRUE,
        TRUE ~ FALSE
      )) %>%
      dplyr::ungroup()
  } else {
    df$is_distal <- NA
  }

  if (!all(is.na(df$pos))) {
    df <- if (all(na.omit(df$strand) == "-")) dplyr::arrange(df, dplyr::desc(pos)) else dplyr::arrange(df, pos)
  } else {
    df <- dplyr::arrange(df, PAS_raw)
  }

  # Build per-PAS labels with ΔPSI and significance marker
  sig_mark <- dplyr::case_when(
    !is.na(df$padj) & df$padj < 0.001 ~ "***",
    !is.na(df$padj) & df$padj < 0.01  ~ "**",
    !is.na(df$padj) & df$padj < 0.05  ~ "*",
    TRUE ~ ""
  )
  base_lab  <- ifelse(!is.na(df$is_distal) & df$is_distal,
                      paste0(df$PAS_raw, " (distal)"), df$PAS_raw)
  delta_psi <- df[[trt_col]] - df[[ctrl_col]]
  PAS_label <- sprintf("%s\nΔPSI=%+.2f%s", base_lab, delta_psi, sig_mark)
  if (any(duplicated(PAS_label))) PAS_label <- make.unique(PAS_label, sep = "_")

  df$PAS_label <- PAS_label

  # Reshape to long, carrying SE through pivot
  long <- df %>%
    dplyr::transmute(
      PAS_label,
      ctrl_u  = .data[[ctrl_col]],
      trt_u   = .data[[trt_col]],
      ctrl_se = if (se_ctrl  %in% names(df)) .data[[se_ctrl]]  else NA_real_,
      trt_se  = if (se_trt   %in% names(df)) .data[[se_trt]]   else NA_real_
    ) %>%
    dplyr::rename(!!CTRL_LABEL := ctrl_u, !!TRTMT_LABEL := trt_u) %>%
    tidyr::pivot_longer(
      cols      = dplyr::all_of(c(CTRL_LABEL, TRTMT_LABEL)),
      names_to  = "condition",
      values_to = "mean_usage"
    ) %>%
    dplyr::mutate(
      se_val    = ifelse(condition == CTRL_LABEL, ctrl_se, trt_se),
      PAS_label = factor(PAS_label, levels = unique(PAS_label)),
      mean_usage = pmax(pmin(mean_usage, 1), 0)
    )

  ttl <- paste0("Relative PAS usage (PSI): ", gene_symbol,
                if (!is.null(title_suffix)) paste0(" — ", title_suffix) else "")

  ggplot2::ggplot(long, ggplot2::aes(x = PAS_label, y = mean_usage, fill = condition)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(mean_usage - se_val, 0),
                   ymax = pmin(mean_usage + se_val, 1)),
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.25, na.rm = TRUE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.1))
    ) +
    ggplot2::labs(
      x     = "PAS site (5' → 3')",
      y     = "Mean fractional usage (PSI)",
      title = ttl,
      fill  = "Condition",
      caption = "ΔPSI = Treatment − Control; error bars = ±1 SEM; * padj<0.05, ** <0.01, *** <0.001"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# ============================================================
#  PRE-ANALYSIS QC
# ============================================================

message("--- Pre-analysis QC ---")

make_library_sizes(
  count_mat_raw, colData,
  outfile = file.path(dir_plots, "library_sizes.png")
)

if (USE_RUV) {
  # W factor bar chart — shows how much each sample is pulled by the batch factor
  w_cols <- paste0("W_", seq_len(RUV_K))
  w_df   <- tibble::rownames_to_column(
               as.data.frame(colData[, c("condition", w_cols), drop = FALSE]), "sample")
  w_long <- tidyr::pivot_longer(w_df, cols = dplyr::all_of(w_cols),
                                 names_to = "factor", values_to = "value")
  gw <- ggplot2::ggplot(w_long, ggplot2::aes(x = sample, y = value, fill = condition)) +
    ggplot2::geom_col() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::facet_wrap(~ factor, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = "Sample", y = "W (RUV factor)",
                  title = "RUVs estimated batch factors") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(dir_plots, "RUV_W_factors.png"), gw,
                  width = max(4, 0.5 * nrow(colData) + 2),
                  height = 3 * RUV_K + 1, dpi = 300)

  # Uncorrected PCA — shows the raw batch structure you are correcting for
  make_pca_global(count_mat, colData,
                  outfile = file.path(dir_plots, "PCA_global_uncorrected.png"))
  # RUV-corrected PCA — W factors regressed from VST matrix for visualization only
  make_pca_global(count_mat, colData,
                  batch_covars = w_cols,
                  outfile = file.path(dir_plots, "PCA_global_RUV_corrected.png"))
} else {
  make_pca_global(count_mat, colData,
                  outfile = file.path(dir_plots, "PCA_global.png"))
}

# ============================================================
#  MAIN LOOP
# ============================================================

message("--- Running DEXSeq analyses ---")

runs <- lapply(names(group_list), function(grp) {
  run_dexseq_group(
    grp_label         = grp,
    sub_samples       = group_list[[grp]],
    count_mat         = count_mat,
    colData           = colData,
    featureID         = featureID,
    groupID           = groupID,
    gr                = gr,
    pas_anno          = pas_anno,
    min_total         = MIN_TOTAL_READS,
    min_per_condition = MIN_PER_CONDITION,
    use_ruv           = USE_RUV,
    ruv_k             = RUV_K
  )
})
names(runs) <- names(group_list)

# ============================================================
#  POST-ANALYSIS QC (per group)
# ============================================================

message("--- Post-analysis QC plots ---")

for (grp in names(runs)) {
  sfx <- group_suffix[grp]
  dxd <- runs[[grp]]$dxd
  res <- runs[[grp]]$results

  make_pca_dxd(
    dxd, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("PCA_normalized%s.png", sfx))
  )
  make_size_factors(
    dxd, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("size_factors%s.png", sfx))
  )
  make_dispersion_plot(
    dxd, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("dispersion%s.png", sfx))
  )
  make_pvalue_hist(
    res, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("pvalue_hist%s.png", sfx))
  )
  make_ma_plot(
    res, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("MA%s.png", sfx))
  )
  make_volcano_plot(
    res, grp_label = grp,
    outfile = file.path(dir_plots, sprintf("volcano%s.png", sfx))
  )
}

# ============================================================
#  EXPORT PER-PAS RESULTS
# ============================================================

message("--- Exporting results ---")

for (grp in names(runs)) {
  sfx <- group_suffix[grp]
  out <- runs[[grp]]$results
  names(out) <- sub(FC_COL, FC_EXPORT, names(out), fixed = TRUE)
  write.csv(out,
            file.path(dir_results, sprintf("pas_results%s.csv", sfx)),
            row.names = FALSE, quote = FALSE)
}

# ============================================================
#  GENE-LEVEL SUMMARIES
# ============================================================

for (grp in names(runs)) {
  sfx    <- group_suffix[grp]
  gene_q <- perGeneQValue(runs[[grp]]$dxr)
  sum_df <- collapse_gene(runs[[grp]]$results, gene_q)
  write.csv(sum_df,
            file.path(dir_results, sprintf("gene_summary%s.csv", sfx)),
            row.names = FALSE, quote = FALSE)
}

# ============================================================
#  USAGE TABLES (mean PSI + SEM per condition)
# ============================================================

usage_results <- lapply(names(group_list), function(grp) {
  cd_sub <- colData[group_list[[grp]], , drop = FALSE]
  usage_by_condition(runs[[grp]]$results, cd_sub)
})
names(usage_results) <- names(group_list)

for (grp in names(usage_results)) {
  sfx <- group_suffix[grp]
  out <- usage_results[[grp]]
  names(out) <- sub(FC_COL, FC_EXPORT, names(out), fixed = TRUE)
  write.csv(out,
            file.path(dir_results, sprintf("pas_usage%s.csv", sfx)),
            row.names = FALSE, quote = FALSE)
}

# ============================================================
#  RELATIVE USAGE BAR CHARTS
# ============================================================

if (length(GENES_OF_INTEREST) > 0) {
  message("--- Plotting PSI bar charts ---")

  for (g in GENES_OF_INTEREST) {
    for (grp in names(usage_results)) {
      sfx     <- group_suffix[grp]
      ttl_sfx <- if (nchar(sfx) > 0) sub("^\\.", "", sfx) else NULL
      outfile <- file.path(dir_usage_plots, sprintf("%s%s.png", g, sfx))
      tryCatch({
        png(outfile, width = 1000, height = 800, res = 100)
        print(plot_gene_usage(usage_results[[grp]], g, title_suffix = ttl_sfx))
        dev.off()
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        warning(sprintf("Could not plot '%s' (%s): %s", g, grp, conditionMessage(e)))
      })
    }
  }
}

message("Done. Outputs written to: ", OUT_BASE)
