#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({

  pkgs <- c(
    "DESeq2",
    "ggplot2",
    "pheatmap",
    "vegan"
  )

  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required R package is missing: ", p)
    }
  }

  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(vegan)
})

cat("\n")
cat("============================================================\n")
cat("08 2021 HISTORICAL-REFERENCE COHORT REANALYSIS\n")
cat("Internal structure and Baseline-relative temporal responses\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_SAMPLES <- 60
EXPECTED_GENES <- 17328

PADJ_THRESHOLD <- 0.05
LFC_THRESHOLD <- 1

N_PERMUTATIONS <- 9999
RANDOM_SEED <- 20260910

set.seed(RANDOM_SEED)

group_levels <- c(
  "Baseline",
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)

injury_groups <- group_levels[
  group_levels != "Baseline"
]

group_labels <- c(
  Baseline = "Baseline",
  R4h = "4 h",
  R8h = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)


# ============================================================
# 1. Input/output
# ============================================================

count_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2021_historical_reference/",
  "02_2021_filtered_counts.tsv"
)

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2021_historical_reference/",
  "02_2021_VST_matrix.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2021_historical_reference/",
  "02_2021_metadata_aligned.tsv"
)

outdir <- "08_2021_historical_reference_reanalysis_result"

de_full_dir <- file.path(
  outdir,
  "DE_full_results"
)

de_sig_dir <- file.path(
  outdir,
  "DE_significant_genes"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  de_full_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  de_sig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Plot helper
# ============================================================

save_plot_both <- function(
    p,
    stem,
    width,
    height) {

  ggsave(
    filename = file.path(
      outdir,
      paste0(stem, ".pdf")
    ),
    plot = p,
    width = width,
    height = height,
    units = "in"
  )

  tifile <- file.path(
    outdir,
    paste0(stem, "_600dpi.tiff")
  )

  if (requireNamespace("ragg", quietly = TRUE)) {

    ragg::agg_tiff(
      filename = tifile,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

  } else {

    grDevices::tiff(
      filename = tifile,
      width = width,
      height = height,
      units = "in",
      res = 600,
      type = "cairo",
      compression = "lzw"
    )
  }

  print(p)
  dev.off()
}


# ============================================================
# 3. Check input files
# ============================================================

required_files <- c(
  count_file,
  vst_file,
  meta_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("Missing input files:\n")

  for (f in missing_files) {
    cat("  ", f, "\n")
  }

  quit(status = 1)
}

cat("All input files found: PASS\n\n")


# ============================================================
# 4. Read counts
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
  stop("Duplicated GeneIDs in count matrix.")
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

if (any(
  abs(count_mat - round(count_mat)) > 1e-8
)) {
  stop("Noninteger values detected.")
}

count_mat <- round(count_mat)
storage.mode(count_mat) <- "integer"


# ============================================================
# 5. Read VST
# ============================================================

cat("Reading VST matrix...\n")

vst_df <- read.delim(
  vst_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

colnames(vst_df)[1] <- "GeneID"

if (anyDuplicated(vst_df$GeneID)) {
  stop("Duplicated GeneIDs in VST matrix.")
}

vst_mat <- as.matrix(
  vst_df[, -1, drop = FALSE]
)

storage.mode(vst_mat) <- "numeric"
rownames(vst_mat) <- vst_df$GeneID

if (anyNA(vst_mat) ||
    any(!is.finite(vst_mat))) {
  stop("Invalid values detected in VST matrix.")
}


# ============================================================
# 6. Read metadata
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
# 7. Explicit alignment
# ============================================================

if (!setequal(
  colnames(count_mat),
  colnames(vst_mat)
)) {
  stop("Count/VST sample sets differ.")
}

if (!setequal(
  colnames(count_mat),
  meta$sampleID
)) {
  stop("Count/metadata sample sets differ.")
}

vst_mat <- vst_mat[
  match(
    rownames(count_mat),
    rownames(vst_mat)
  ),
  ,
  drop = FALSE
]

if (anyNA(match(
  rownames(count_mat),
  rownames(vst_mat)
))) {
  stop("Count/VST GeneID alignment failed.")
}

vst_mat <- vst_mat[
  rownames(count_mat),
  colnames(count_mat),
  drop = FALSE
]

meta <- meta[
  match(
    colnames(count_mat),
    meta$sampleID
  ),
  ,
  drop = FALSE
]

if (!identical(
  colnames(count_mat),
  colnames(vst_mat)
)) {
  stop("Count/VST sample-order alignment failed.")
}

if (!identical(
  colnames(count_mat),
  meta$sampleID
)) {
  stop("Count/metadata alignment failed.")
}

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)

if (anyNA(meta$group)) {
  stop("Unexpected 2021 group detected.")
}


# ============================================================
# 8. Dimension validation
# ============================================================

cat("\nInput dimensions:\n")

cat(
  "Counts   : ",
  nrow(count_mat),
  " genes x ",
  ncol(count_mat),
  " samples\n",
  sep = ""
)

cat(
  "VST      : ",
  nrow(vst_mat),
  " genes x ",
  ncol(vst_mat),
  " samples\n",
  sep = ""
)

cat(
  "Metadata : ",
  nrow(meta),
  " samples\n\n",
  sep = ""
)

if (nrow(count_mat) != EXPECTED_GENES ||
    nrow(vst_mat) != EXPECTED_GENES) {
  stop("Unexpected gene number.")
}

if (ncol(count_mat) != EXPECTED_SAMPLES ||
    ncol(vst_mat) != EXPECTED_SAMPLES ||
    nrow(meta) != EXPECTED_SAMPLES) {
  stop("Unexpected sample number.")
}

cat("Dimension validation: PASS\n\n")


# ============================================================
# 9. Sample distribution
# ============================================================

sample_summary <- data.frame(
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
  sample_summary,
  file = file.path(
    outdir,
    "08_2021_group_sample_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Sample distribution:\n")
print(sample_summary, row.names = FALSE)
cat("\n")


# ============================================================
# 10. PCA
# ============================================================

cat("Running 2021 PCA...\n")

pca <- prcomp(
  t(vst_mat),
  center = TRUE,
  scale. = FALSE
)

variance <- pca$sdev^2

variance_percent <- 100 *
  variance /
  sum(variance)

pca_variance <- data.frame(
  PC = paste0(
    "PC",
    seq_along(variance_percent)
  ),
  Variance_percent = variance_percent,
  Cumulative_percent = cumsum(
    variance_percent
  ),
  stringsAsFactors = FALSE
)

write.table(
  pca_variance,
  file = file.path(
    outdir,
    "08_2021_PCA_variance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

pca_df <- data.frame(
  sampleID = rownames(pca$x),
  pca$x,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pca_df$Group <- as.character(
  meta[
    pca_df$sampleID,
    "group"
  ]
)

pca_df$Stage <- factor(
  unname(
    group_labels[
      pca_df$Group
    ]
  ),
  levels = unname(
    group_labels[
      group_levels
    ]
  )
)

write.table(
  pca_df,
  file = file.path(
    outdir,
    "08_2021_PCA_coordinates.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

stage_palette <- setNames(
  grDevices::hcl.colors(
    length(group_levels),
    palette = "Dynamic"
  ),
  unname(group_labels[group_levels])
)

p_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = Stage
  )
) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_colour_manual(
    values = stage_palette
  ) +
  labs(
    x = paste0(
      "PC1 (",
      sprintf("%.2f", variance_percent[1]),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      sprintf("%.2f", variance_percent[2]),
      "%)"
    ),
    colour = "Sampling stage"
  ) +
  theme_classic(
    base_size = 12
  )

save_plot_both(
  p_pca,
  "08_2021_PCA_PC1_PC2",
  8.5,
  6.5
)

cat(
  "PCA variance: PC1 = ",
  sprintf("%.2f", variance_percent[1]),
  "%; PC2 = ",
  sprintf("%.2f", variance_percent[2]),
  "%\n\n",
  sep = ""
)


# ============================================================
# 11. Spearman sample correlation
# ============================================================

cat("Calculating 2021 sample Spearman correlations...\n")

cor_mat <- cor(
  vst_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

diag(cor_mat) <- 1

if (anyNA(cor_mat)) {
  stop("NA detected in correlation matrix.")
}

write.table(
  data.frame(
    sampleID = rownames(cor_mat),
    cor_mat,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "08_2021_sample_Spearman_correlation_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 12. Correlation-based clustering
# ============================================================

cor_distance_mat <- 1 - cor_mat

cor_distance_mat[
  cor_distance_mat < 0
] <- 0

diag(cor_distance_mat) <- 0

if (nrow(cor_distance_mat) !=
    ncol(cor_distance_mat)) {
  stop("Correlation distance matrix is not square.")
}

cor_distance <- as.dist(
  cor_distance_mat
)

sample_hclust <- hclust(
  cor_distance,
  method = "average"
)

annotation <- data.frame(
  Stage = factor(
    unname(
      group_labels[
        as.character(meta$group)
      ]
    ),
    levels = unname(
      group_labels[group_levels]
    )
  ),
  row.names = meta$sampleID
)

pdf(
  file.path(
    outdir,
    "08_2021_sample_Spearman_correlation_heatmap.pdf"
  ),
  width = 12,
  height = 11
)

pheatmap(
  cor_mat,
  cluster_rows = sample_hclust,
  cluster_cols = sample_hclust,
  annotation_row = annotation,
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 5,
  fontsize_col = 5,
  border_color = NA,
  main = "2021 historical-reference cohort"
)

dev.off()

cat("Spearman correlation and clustering: PASS\n\n")


# ============================================================
# 13. Euclidean VST distance
# ============================================================

cat("Calculating full-expression Euclidean distances...\n")

sample_distance <- dist(
  t(vst_mat),
  method = "euclidean"
)


# ============================================================
# 14. PERMANOVA
# ============================================================

cat(
  "Running PERMANOVA with ",
  N_PERMUTATIONS,
  " permutations...\n",
  sep = ""
)

set.seed(RANDOM_SEED)

permanova <- vegan::adonis2(
  sample_distance ~ group,
  data = meta,
  permutations = N_PERMUTATIONS
)

write.table(
  data.frame(
    Term = rownames(permanova),
    as.data.frame(permanova),
    row.names = NULL,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "08_2021_PERMANOVA_stage.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nPERMANOVA:\n")
print(permanova)
cat("\n")


# ============================================================
# 15. PERMDISP
# ============================================================

cat(
  "Running PERMDISP with ",
  N_PERMUTATIONS,
  " permutations...\n",
  sep = ""
)

bd <- vegan::betadisper(
  sample_distance,
  group = meta$group,
  type = "median",
  bias.adjust = TRUE
)

set.seed(RANDOM_SEED)

bd_perm <- vegan::permutest(
  bd,
  permutations = N_PERMUTATIONS
)

write.table(
  data.frame(
    Term = rownames(bd_perm$tab),
    as.data.frame(bd_perm$tab),
    row.names = NULL,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "08_2021_PERMDISP_permutation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nPERMDISP:\n")
print(bd_perm)
cat("\n")


# ============================================================
# 16. DESeq2 model
# ============================================================

cat("Creating 2021 DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ group
)

dds$group <- relevel(
  dds$group,
  ref = "Baseline"
)

cat("Running DESeq2 Wald model...\n")

dds <- DESeq(
  dds,
  test = "Wald",
  minReplicatesForReplace = Inf,
  quiet = FALSE
)

cat("\nDESeq2 fitting completed.\n\n")


# ============================================================
# 17. Extract Baseline-relative contrasts
# ============================================================

de_summary_list <- list()
full_results_list <- list()

for (grp in injury_groups) {

  cat(
    "------------------------------------------------------------\n"
  )

  cat(
    "Contrast: ",
    grp,
    " vs Baseline\n",
    sep = ""
  )

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

  x <- as.data.frame(res)

  x$GeneID <- rownames(x)

  x <- x[
    ,
    c(
      "GeneID",
      "baseMean",
      "log2FoldChange",
      "lfcSE",
      "stat",
      "pvalue",
      "padj"
    ),
    drop = FALSE
  ]

  x$DE_status <- "Not_significant"

  x$DE_status[
    !is.na(x$padj) &
    x$padj < PADJ_THRESHOLD &
    x$log2FoldChange >= LFC_THRESHOLD
  ] <- "Up"

  x$DE_status[
    !is.na(x$padj) &
    x$padj < PADJ_THRESHOLD &
    x$log2FoldChange <= -LFC_THRESHOLD
  ] <- "Down"

  x$Contrast <- paste0(
    grp,
    "_vs_Baseline"
  )

  x$Stage <- unname(
    group_labels[grp]
  )

  x <- x[
    order(
      is.na(x$padj),
      x$padj,
      x$pvalue
    ),
    ,
    drop = FALSE
  ]

  write.table(
    x,
    file = file.path(
      de_full_dir,
      paste0(
        "08_",
        grp,
        "_vs_Baseline_full.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  sig <- x[
    x$DE_status %in% c(
      "Up",
      "Down"
    ),
    ,
    drop = FALSE
  ]

  write.table(
    sig,
    file = file.path(
      de_sig_dir,
      paste0(
        "08_",
        grp,
        "_vs_Baseline_DEG.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  n_up <- sum(
    x$DE_status == "Up"
  )

  n_down <- sum(
    x$DE_status == "Down"
  )

  n_padj <- sum(
    !is.na(x$padj) &
    x$padj < PADJ_THRESHOLD
  )

  de_summary_list[[grp]] <- data.frame(
    Group = grp,
    Stage = unname(
      group_labels[grp]
    ),
    N_group = sum(
      meta$group == grp
    ),
    N_baseline = sum(
      meta$group == "Baseline"
    ),
    Tested_genes = nrow(x),
    Genes_with_PADJ_NA = sum(
      is.na(x$padj)
    ),
    PADJ_LT_0.05 = n_padj,
    Up = n_up,
    Down = n_down,
    Total_DEG = n_up + n_down,
    stringsAsFactors = FALSE
  )

  full_results_list[[grp]] <- x

  cat(
    "padj < 0.05 : ",
    n_padj,
    "\n",
    sep = ""
  )

  cat(
    "Up          : ",
    n_up,
    "\n",
    sep = ""
  )

  cat(
    "Down        : ",
    n_down,
    "\n",
    sep = ""
  )

  cat(
    "Total DEG   : ",
    n_up + n_down,
    "\n",
    sep = ""
  )
}

cat(
  "------------------------------------------------------------\n\n"
)


# ============================================================
# 18. DEG summary
# ============================================================

de_summary <- do.call(
  rbind,
  de_summary_list
)

rownames(de_summary) <- NULL

write.table(
  de_summary,
  file = file.path(
    outdir,
    "08_2021_DEG_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Genome-wide log2FC/PADJ matrices
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

  x <- full_results_list[[grp]]

  idx <- match(
    gene_order,
    x$GeneID
  )

  if (anyNA(idx)) {
    stop(
      "Gene matching failed for ",
      grp
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
    "08_2021_genomewide_log2FC_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  padj_matrix,
  file = file.path(
    outdir,
    "08_2021_genomewide_PADJ_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. DEG overview plot
# ============================================================

plot_df <- rbind(
  data.frame(
    Stage = de_summary$Stage,
    Direction = "Up",
    Number = de_summary$Up
  ),
  data.frame(
    Stage = de_summary$Stage,
    Direction = "Down",
    Number = -de_summary$Down
  )
)

plot_df$Stage <- factor(
  plot_df$Stage,
  levels = unname(
    group_labels[injury_groups]
  )
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
  )

save_plot_both(
  p_deg,
  "08_2021_DEG_number_over_time",
  8.5,
  5.8
)


# ============================================================
# 21. Extract global statistics
# ============================================================

if (!"Model" %in% rownames(permanova)) {
  stop("Unexpected PERMANOVA output.")
}

permanova_r2 <- as.numeric(
  permanova[
    "Model",
    "R2"
  ]
)

permanova_f <- as.numeric(
  permanova[
    "Model",
    "F"
  ]
)

permanova_p <- as.numeric(
  permanova[
    "Model",
    "Pr(>F)"
  ]
)

if (!"Groups" %in% rownames(bd_perm$tab)) {
  stop("Unexpected PERMDISP output.")
}

permdisp_p <- as.numeric(
  bd_perm$tab[
    "Groups",
    "Pr(>F)"
  ]
)


# ============================================================
# 22. Summary
# ============================================================

summary_table <- data.frame(
  Metric = c(
    "Samples",
    "Genes",
    "Baseline_samples",
    "Injury_samples",
    "PC1_variance_percent",
    "PC2_variance_percent",
    "PERMANOVA_R2",
    "PERMANOVA_F",
    "PERMANOVA_P",
    "PERMDISP_P",
    "DESeq2_contrasts"
  ),
  Value = c(
    ncol(count_mat),
    nrow(count_mat),
    sum(meta$group == "Baseline"),
    sum(meta$group != "Baseline"),
    variance_percent[1],
    variance_percent[2],
    permanova_r2,
    permanova_f,
    permanova_p,
    permdisp_p,
    length(injury_groups)
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    outdir,
    "08_2021_historical_reference_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Save objects/session info
# ============================================================

saveRDS(
  dds,
  file = file.path(
    outdir,
    "08_2021_DESeq2_fitted_object.rds"
  )
)

saveRDS(
  pca,
  file = file.path(
    outdir,
    "08_2021_PCA_object.rds"
  )
)

sink(
  file.path(
    outdir,
    "08_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("08 2021 HISTORICAL-REFERENCE REANALYSIS COMPLETED\n")
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
  "PC1 variance                  : ",
  sprintf("%.2f%%", variance_percent[1]),
  "\n",
  sep = ""
)

cat(
  "PC2 variance                  : ",
  sprintf("%.2f%%", variance_percent[2]),
  "\n",
  sep = ""
)

cat("------------------------------------------------------------\n")

cat(
  "PERMANOVA stage R2            : ",
  sprintf("%.6f", permanova_r2),
  "\n",
  sep = ""
)

cat(
  "PERMANOVA stage F             : ",
  sprintf("%.6f", permanova_f),
  "\n",
  sep = ""
)

cat(
  "PERMANOVA stage P             : ",
  format(
    permanova_p,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "PERMDISP stage P              : ",
  format(
    permdisp_p,
    scientific = TRUE
  ),
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
cat("2021 is treated as a previously published historical reference.\n")
cat("No pooling or batch correction with the 2026 cohort was performed.\n")
cat("Raw counts were analyzed with DESeq2.\n")
cat("Automatic Cook's count replacement was disabled.\n")
cat("Independent filtering in results() was disabled.\n")
cat("Genome-wide Baseline-relative log2FC values were retained\n")
cat("for subsequent cross-cohort temporal concordance analysis.\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 08 is complete. ",
  "The next step is Step 09: ",
  "2026-versus-2021 cross-cohort temporal-response concordance.\n"
)
