#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({

  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("R package 'DESeq2' is required.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }

  library(DESeq2)
  library(ggplot2)
})

cat("\n")
cat("============================================================\n")
cat("05 2026 TIMEPOINT-vs-BASELINE DIFFERENTIAL EXPRESSION\n")
cat("DESeq2 analysis of the 2026 primary cohort\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed analysis parameters
# ============================================================

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

EXPECTED_SAMPLES <- 77
EXPECTED_GENES <- 18364

group_levels <- c(
  "Baseline",
  "R0h",
  "R1h",
  "R2h",
  "R6h",
  "R10h",
  "R14h",
  "R18h",
  "R36h",
  "R60h",
  "R72h"
)

injury_groups <- group_levels[
  group_levels != "Baseline"
]

group_labels <- c(
  Baseline = "Baseline",
  R0h = "0 h",
  R1h = "1 h",
  R2h = "2 h",
  R6h = "6 h",
  R10h = "10 h",
  R14h = "14 h",
  R18h = "18 h",
  R36h = "36 h",
  R60h = "60 h",
  R72h = "72 h"
)


# ============================================================
# 1. Input/output
# ============================================================

count_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_filtered_counts.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_metadata_aligned.tsv"
)

outdir <- "05_2026_timepoint_vs_baseline_DESeq2_result"

result_dir <- file.path(
  outdir,
  "DE_full_results"
)

significant_dir <- file.path(
  outdir,
  "DE_significant_genes"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  result_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  significant_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Plot helper
# ============================================================

save_plot_both <- function(
    plot_object,
    filename_stem,
    width,
    height) {

  ggsave(
    file.path(
      outdir,
      paste0(filename_stem, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in"
  )

  tiff_file <- file.path(
    outdir,
    paste0(filename_stem, "_600dpi.tiff")
  )

  if (requireNamespace("ragg", quietly = TRUE)) {

    ragg::agg_tiff(
      filename = tiff_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

    print(plot_object)
    dev.off()

  } else {

    grDevices::tiff(
      filename = tiff_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      type = "cairo",
      compression = "lzw"
    )

    print(plot_object)
    dev.off()
  }
}


# ============================================================
# 3. Input validation
# ============================================================

required_files <- c(
  count_file,
  meta_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing input files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("Input files found: PASS\n\n")


# ============================================================
# 4. Read filtered counts
# ============================================================

cat("Reading filtered raw counts...\n")

count_df <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(count_df)[1] <- "GeneID"

if (anyDuplicated(count_df$GeneID)) {
  stop("Duplicated GeneIDs detected.")
}

count_mat <- as.matrix(
  count_df[, -1, drop = FALSE]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- count_df$GeneID

if (anyNA(count_mat) ||
    any(!is.finite(count_mat)) ||
    any(count_mat < 0)) {

  stop("Invalid values detected in count matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Noninteger values detected in count matrix.")
}

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"


# ============================================================
# 5. Read metadata
# ============================================================

cat("Reading metadata...\n")

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

required_meta <- c(
  "sampleID",
  "group",
  "time",
  "condition"
)

missing_meta <- setdiff(
  required_meta,
  colnames(meta)
)

if (length(missing_meta) > 0) {

  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop("Duplicated sampleIDs detected.")
}


# ============================================================
# 6. Align metadata to count matrix
# ============================================================

sample_ids <- colnames(count_mat)

if (!setequal(
  sample_ids,
  meta$sampleID
)) {

  stop(
    "Count and metadata sample sets differ."
  )
}

meta <- meta[
  match(sample_ids, meta$sampleID),
  ,
  drop = FALSE
]

if (!identical(
  sample_ids,
  meta$sampleID
)) {

  stop(
    "Metadata alignment failed."
  )
}

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)

if (anyNA(meta$group)) {
  stop("Unexpected group detected.")
}


# ============================================================
# 7. Dimension validation
# ============================================================

cat("\nInput dimensions:\n")

cat(
  "Counts   :",
  nrow(count_mat),
  "genes x",
  ncol(count_mat),
  "samples\n"
)

cat(
  "Metadata :",
  nrow(meta),
  "samples\n\n"
)

if (nrow(count_mat) != EXPECTED_GENES) {
  stop(
    "Expected ",
    EXPECTED_GENES,
    " genes but observed ",
    nrow(count_mat)
  )
}

if (ncol(count_mat) != EXPECTED_SAMPLES) {
  stop(
    "Expected ",
    EXPECTED_SAMPLES,
    " samples but observed ",
    ncol(count_mat)
  )
}

cat("Dimension validation: PASS\n")
cat("All 77 samples retained: PASS\n\n")


# ============================================================
# 8. Group sample numbers
# ============================================================

group_summary <- data.frame(
  Group = group_levels,
  Stage = unname(
    group_labels[group_levels]
  ),
  N = as.integer(
    table(
      factor(
        meta$group,
        levels = group_levels
      )
    )
  ),
  stringsAsFactors = FALSE
)

write.table(
  group_summary,
  file = file.path(
    outdir,
    "05_2026_group_sample_numbers.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Sample distribution:\n")
print(
  group_summary,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 9. Create DESeq2 dataset
# ============================================================

cat("Creating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ group
)

# Explicitly confirm Baseline reference
dds$group <- relevel(
  dds$group,
  ref = "Baseline"
)


# ============================================================
# 10. Fit DESeq2 model
# ============================================================

cat("Running DESeq2 model fitting...\n")

dds <- DESeq(
  dds,
  test = "Wald",
  minReplicatesForReplace = Inf,
  quiet = FALSE
)

cat("\nDESeq2 fitting completed.\n\n")

cat("Available coefficient names:\n")
print(resultsNames(dds))
cat("\n")


# ============================================================
# 11. Size-factor output
# ============================================================

size_factor_table <- data.frame(
  sampleID = colnames(dds),
  Group = as.character(colData(dds)$group),
  Size_factor = sizeFactors(dds),
  stringsAsFactors = FALSE
)

write.table(
  size_factor_table,
  file = file.path(
    outdir,
    "05_2026_DESeq2_size_factors.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 12. Contrast analysis
# ============================================================

cat("Extracting timepoint-vs-Baseline contrasts...\n\n")

all_full_results <- list()
summary_list <- list()

for (grp in injury_groups) {

  stage_label <- unname(
    group_labels[grp]
  )

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Contrast: ",
    grp,
    " vs Baseline\n",
    sep = ""
  )

  # Independent filtering is disabled because genes were
  # already filtered using the predefined count-based rule.
  res <- results(
    dds,
    contrast = c(
      "group",
      grp,
      "Baseline"
    ),
    alpha = PADJ_THRESHOLD,
    independentFiltering = FALSE
  )

  res_df <- as.data.frame(res)

  res_df$GeneID <- rownames(res_df)

  res_df <- res_df[
    ,
    c(
      "GeneID",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    )
  ]

  # Classification
  res_df$DE_status <- "Not_significant"

  res_df$DE_status[
    !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD &
    res_df$log2FoldChange >= LFC_THRESHOLD
  ] <- "Up"

  res_df$DE_status[
    !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD &
    res_df$log2FoldChange <= -LFC_THRESHOLD
  ] <- "Down"

  # Add contrast information
  res_df$Contrast <- paste0(
    grp,
    "_vs_Baseline"
  )

  res_df$Stage <- stage_label

  # Sort full results by padj, then pvalue
  ord <- order(
    is.na(res_df$padj),
    res_df$padj,
    res_df$pvalue
  )

  res_df <- res_df[
    ord,
    ,
    drop = FALSE
  ]

  # Save full result
  full_file <- file.path(
    result_dir,
    paste0(
      "05_",
      grp,
      "_vs_Baseline_full.tsv"
    )
  )

  write.table(
    res_df,
    file = full_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # Significant DE genes
  sig_df <- res_df[
    res_df$DE_status %in% c(
      "Up",
      "Down"
    ),
    ,
    drop = FALSE
  ]

  sig_file <- file.path(
    significant_dir,
    paste0(
      "05_",
      grp,
      "_vs_Baseline_DEG.tsv"
    )
  )

  write.table(
    sig_df,
    file = sig_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  n_up <- sum(
    res_df$DE_status == "Up",
    na.rm = TRUE
  )

  n_down <- sum(
    res_df$DE_status == "Down",
    na.rm = TRUE
  )

  n_sig_padj <- sum(
    !is.na(res_df$padj) &
    res_df$padj < PADJ_THRESHOLD
  )

  n_na_padj <- sum(
    is.na(res_df$padj)
  )

  summary_list[[grp]] <- data.frame(
    Group = grp,
    Stage = stage_label,
    N_group = sum(meta$group == grp),
    N_baseline = sum(meta$group == "Baseline"),
    Tested_genes = nrow(res_df),
    Genes_with_PADJ_NA = n_na_padj,
    PADJ_LT_0.05 = n_sig_padj,
    Up = n_up,
    Down = n_down,
    Total_DEG = n_up + n_down,
    PADJ_threshold = PADJ_THRESHOLD,
    Absolute_LFC_threshold = LFC_THRESHOLD,
    stringsAsFactors = FALSE
  )

  all_full_results[[grp]] <- res_df

  cat(
    "padj < 0.05          : ",
    n_sig_padj,
    "\n",
    sep = ""
  )

  cat(
    "Up                  : ",
    n_up,
    "\n",
    sep = ""
  )

  cat(
    "Down                : ",
    n_down,
    "\n",
    sep = ""
  )

  cat(
    "Total DEG           : ",
    n_up + n_down,
    "\n",
    sep = ""
  )
}

cat(
  "------------------------------------------------------------\n\n"
)


# ============================================================
# 13. Combined DEG summary
# ============================================================

de_summary <- do.call(
  rbind,
  summary_list
)

rownames(de_summary) <- NULL

write.table(
  de_summary,
  file = file.path(
    outdir,
    "05_2026_DEG_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 14. Combined full-statistics long table
# ============================================================

all_full_long <- do.call(
  rbind,
  all_full_results
)

rownames(all_full_long) <- NULL

write.table(
  all_full_long,
  file = file.path(
    outdir,
    "05_2026_all_timepoints_full_statistics_long.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 15. Create genome-wide log2FC matrix
# This will later be used for 2026-vs-2021 concordance
# ============================================================

gene_order <- rownames(dds)

logfc_matrix <- data.frame(
  GeneID = gene_order,
  stringsAsFactors = FALSE
)

padj_matrix <- data.frame(
  GeneID = gene_order,
  stringsAsFactors = FALSE
)

for (grp in injury_groups) {

  x <- all_full_results[[grp]]

  idx <- match(
    gene_order,
    x$GeneID
  )

  if (anyNA(idx)) {
    stop(
      "Gene matching failed while constructing log2FC matrix."
    )
  }

  logfc_matrix[[grp]] <-
    x$log2FoldChange[idx]

  padj_matrix[[grp]] <-
    x$padj[idx]
}

write.table(
  logfc_matrix,
  file = file.path(
    outdir,
    "05_2026_genomewide_log2FC_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  padj_matrix,
  file = file.path(
    outdir,
    "05_2026_genomewide_PADJ_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 16. DEG-number plot
# ============================================================

plot_df <- rbind(
  data.frame(
    Stage = de_summary$Stage,
    Direction = "Up",
    Number = de_summary$Up,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Stage = de_summary$Stage,
    Direction = "Down",
    Number = -de_summary$Down,
    stringsAsFactors = FALSE
  )
)

plot_df$Stage <- factor(
  plot_df$Stage,
  levels = unname(
    group_labels[injury_groups]
  )
)

plot_df$Direction <- factor(
  plot_df$Direction,
  levels = c(
    "Up",
    "Down"
  )
)

direction_palette <- c(
  Up = "#D95F02",
  Down = "#1B9E77"
)

p_deg <- ggplot(
  plot_df,
  aes(
    x = Stage,
    y = Number,
    fill = Direction
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  scale_fill_manual(
    values = direction_palette
  ) +
  scale_y_continuous(
    labels = function(x) abs(x)
  ) +
  labs(
    x = "Sampling stage",
    y = "Number of differentially expressed genes",
    fill = "Direction"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

save_plot_both(
  p_deg,
  "05_2026_DEG_number_over_time",
  width = 9,
  height = 6
)


# ============================================================
# 17. Total DEG plot
# ============================================================

total_plot <- de_summary

total_plot$Stage <- factor(
  total_plot$Stage,
  levels = unname(
    group_labels[injury_groups]
  )
)

p_total <- ggplot(
  total_plot,
  aes(
    x = Stage,
    y = Total_DEG
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_point(
    size = 2
  ) +
  labs(
    x = "Sampling stage",
    y = "Total number of differentially expressed genes"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

save_plot_both(
  p_total,
  "05_2026_total_DEG_number_over_time",
  width = 9,
  height = 5.8
)


# ============================================================
# 18. Save DESeq2 object
# ============================================================

saveRDS(
  dds,
  file = file.path(
    outdir,
    "05_2026_DESeq2_fitted_object.rds"
  )
)


# ============================================================
# 19. Session information
# ============================================================

sink(
  file.path(
    outdir,
    "05_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 20. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("05 2026 DIFFERENTIAL EXPRESSION COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analyzed              : ",
  ncol(count_mat),
  "\n",
  sep = ""
)

cat(
  "Genes analyzed                : ",
  nrow(count_mat),
  "\n",
  sep = ""
)

cat(
  "Baseline samples              : ",
  sum(meta$group == "Baseline"),
  "\n",
  sep = ""
)

cat(
  "Injury samples                : ",
  sum(meta$group != "Baseline"),
  "\n",
  sep = ""
)

cat(
  "Contrasts                     : ",
  length(injury_groups),
  "\n",
  sep = ""
)

cat(
  "Adjusted-P threshold          : ",
  PADJ_THRESHOLD,
  "\n",
  sep = ""
)

cat(
  "Absolute log2FC threshold     : ",
  LFC_THRESHOLD,
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

print(
  de_summary[
    ,
    c(
      "Stage",
      "PADJ_LT_0.05",
      "Up",
      "Down",
      "Total_DEG"
    )
  ],
  row.names = FALSE
)

cat("------------------------------------------------------------\n")
cat("All 77 samples were retained.\n")
cat("No CPM or logCPM was used.\n")
cat("Raw counts were analyzed with DESeq2.\n")
cat("Automatic Cook's outlier count replacement was disabled.\n")
cat("Independent filtering in results() was disabled because\n")
cat("genes had already passed the predefined count filter.\n")
cat("No pathway enrichment analysis was performed.\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 05 is complete. ",
  "The next step is Step 06: ",
  "continuous post-injury temporal modeling using ",
  "DESeq2 natural-spline likelihood-ratio testing.\n"
)
