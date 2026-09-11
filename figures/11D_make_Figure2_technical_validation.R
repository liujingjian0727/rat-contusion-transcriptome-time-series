#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

cat("\n")
cat("============================================================\n")
cat("11D MAKE FIGURE 2\n")
cat("Internal technical validation of the 2026 transcriptome resource\n")
cat("Scientific Data final main figure\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed settings
# ============================================================

groups <- c(
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

stage_labels <- c(
  Baseline = "Baseline",
  R0h  = "0 h",
  R1h  = "1 h",
  R2h  = "2 h",
  R6h  = "6 h",
  R10h = "10 h",
  R14h = "14 h",
  R18h = "18 h",
  R36h = "36 h",
  R60h = "60 h",
  R72h = "72 h"
)

stage_colors <- c(
  Baseline = "#7A7A7A",
  R0h  = "#4C78A8",
  R1h  = "#5A88B5",
  R2h  = "#6998C1",
  R6h  = "#72A5B5",
  R10h = "#79B39C",
  R14h = "#8DBB78",
  R18h = "#B4B85D",
  R36h = "#D1A354",
  R60h = "#D68155",
  R72h = "#C85A5A"
)


# ============================================================
# 1. Input files
# ============================================================

raw_count_file <-
  "rattus_counts_2026_primary_77samples.txt"

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_VST_matrix.tsv"
)

meta_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2026_primary/",
  "02_2026_metadata_aligned.tsv"
)

review_file <- paste0(
  "03_2026_internal_expression_structure_result/",
  "03_2026_expression_outlier_diagnostics.tsv"
)

required_files <- c(
  raw_count_file,
  vst_file,
  meta_file,
  review_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing required file(s):\n")

  for (f in missing_files) {
    cat("  - ", f, "\n", sep = "")
  }

  quit(status = 1)
}

cat("All required input files found: PASS\n\n")


# ============================================================
# 2. Output directories
# ============================================================

outdir <- "final_submission_figures"

panel_dir <- file.path(
  outdir,
  "panels"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  panel_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Generic matrix reader
# ============================================================

read_gene_matrix <- function(file) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  colnames(x)[1] <- "GeneID"

  if (anyDuplicated(x$GeneID)) {
    stop(
      "Duplicated GeneID in ",
      file
    )
  }

  x
}


# ============================================================
# 4. Read raw counts
# ============================================================

cat("Reading raw count matrix...\n")

raw_df <- read_gene_matrix(
  raw_count_file
)

raw_counts <- as.matrix(
  raw_df[, -1, drop = FALSE]
)

storage.mode(raw_counts) <- "numeric"

rownames(raw_counts) <- raw_df$GeneID

if (anyNA(raw_counts) ||
    any(!is.finite(raw_counts)) ||
    any(raw_counts < 0)) {

  stop(
    "Invalid values detected in raw count matrix."
  )
}

if (ncol(raw_counts) != 77) {

  stop(
    "Expected 77 samples in raw counts; observed ",
    ncol(raw_counts)
  )
}

cat(
  "Raw counts : ",
  nrow(raw_counts),
  " genes x ",
  ncol(raw_counts),
  " samples\n",
  sep = ""
)


# ============================================================
# 5. Read VST
# ============================================================

cat("Reading VST matrix...\n")

vst_df <- read_gene_matrix(
  vst_file
)

vst_mat <- as.matrix(
  vst_df[, -1, drop = FALSE]
)

storage.mode(vst_mat) <- "numeric"

rownames(vst_mat) <- vst_df$GeneID

if (nrow(vst_mat) != 18364) {

  warning(
    "Expected 18,364 VST genes; observed ",
    nrow(vst_mat)
  )
}

if (ncol(vst_mat) != 77) {

  stop(
    "Expected 77 VST samples; observed ",
    ncol(vst_mat)
  )
}

cat(
  "VST matrix : ",
  nrow(vst_mat),
  " genes x ",
  ncol(vst_mat),
  " samples\n",
  sep = ""
)


# ============================================================
# 6. Read metadata
# ============================================================

meta <- read.delim(
  meta_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!all(
  c(
    "sampleID",
    "group"
  ) %in%
    colnames(meta)
)) {

  stop(
    "Metadata must contain sampleID and group."
  )
}

if (nrow(meta) != 77) {

  stop(
    "Expected 77 metadata rows; observed ",
    nrow(meta)
  )
}


# ============================================================
# 7. Validate sample sets
# ============================================================

if (!setequal(
  colnames(raw_counts),
  meta$sampleID
)) {

  stop(
    "Raw-count and metadata sample sets differ."
  )
}

if (!setequal(
  colnames(vst_mat),
  meta$sampleID
)) {

  stop(
    "VST and metadata sample sets differ."
  )
}

meta <- meta[
  match(
    colnames(vst_mat),
    meta$sampleID
  ),
  ,
  drop = FALSE
]

rownames(meta) <- meta$sampleID

raw_counts <- raw_counts[
  ,
  meta$sampleID,
  drop = FALSE
]

vst_mat <- vst_mat[
  ,
  meta$sampleID,
  drop = FALSE
]

meta$group <- factor(
  meta$group,
  levels = groups
)

cat("Sample alignment: PASS\n\n")


# ============================================================
# 8. Read Step 03 review diagnostics
# ============================================================

review <- read.delim(
  review_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!all(
  c(
    "sampleID",
    "Review_level"
  ) %in%
    colnames(review)
)) {

  stop(
    "Review table must contain sampleID and Review_level."
  )
}

if (!setequal(
  review$sampleID,
  meta$sampleID
)) {

  stop(
    "Review-diagnostic and metadata sample sets differ."
  )
}

review <- review[
  match(
    meta$sampleID,
    review$sampleID
  ),
  ,
  drop = FALSE
]


# ============================================================
# 9. Normalize Review_level labels
# ============================================================

normalize_review <- function(x) {

  out <- rep(
    NA_character_,
    length(x)
  )

  out[
    grepl(
      "multi",
      x,
      ignore.case = TRUE
    )
  ] <- "Multi-metric review"

  out[
    grepl(
      "single",
      x,
      ignore.case = TRUE
    )
  ] <- "Single-metric review"

  remaining <- is.na(out)

  out[remaining] <- "No flag"

  out
}

review$Review_status <-
  normalize_review(
    review$Review_level
  )

review$Review_status <- factor(
  review$Review_status,
  levels = c(
    "No flag",
    "Single-metric review",
    "Multi-metric review"
  )
)

review_counts <- table(
  review$Review_status
)

cat("Expression-level review status:\n")
print(review_counts)
cat("\n")


# ============================================================
# 10. Join metadata/review
# ============================================================

meta$Review_status <-
  review$Review_status[
    match(
      meta$sampleID,
      review$sampleID
    )
  ]


# ============================================================
# 11. Publication theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 10
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    axis.title = element_text(
      size = 10
    ),
    axis.text = element_text(
      size = 8.5
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# 12. Saving functions
# ============================================================

open_tiff <- function(
    file,
    width,
    height) {

  if (requireNamespace(
    "ragg",
    quietly = TRUE
  )) {

    ragg::agg_tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

  } else {

    grDevices::tiff(
      filename = file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      type = "cairo",
      compression = "lzw"
    )
  }
}


save_panel <- function(
    p,
    stem,
    width,
    height) {

  grDevices::cairo_pdf(
    filename = file.path(
      panel_dir,
      paste0(
        stem,
        ".pdf"
      )
    ),
    width = width,
    height = height
  )

  print(p)
  dev.off()


  open_tiff(
    file.path(
      panel_dir,
      paste0(
        stem,
        "_600dpi.tiff"
      )
    ),
    width,
    height
  )

  print(p)
  dev.off()
}


save_final_figure <- function(
    pA,
    pB,
    pC,
    pD,
    width = 13,
    height = 10) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 2,
          ncol = 2,
          widths = unit(
            c(
              1,
              1
            ),
            "null"
          ),
          heights = unit(
            c(
              0.92,
              1.08
            ),
            "null"
          )
        )
      )
    )

    print(
      pA,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    print(
      pB,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )

    print(
      pC,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
      )
    )

    print(
      pD,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 2
      )
    )
  }


  pdf_file <- file.path(
    outdir,
    paste0(
      "Figure_2_internal_technical_",
      "validation.pdf"
    )
  )

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width,
    height = height
  )

  draw()
  dev.off()


  tif_file <- file.path(
    outdir,
    paste0(
      "Figure_2_internal_technical_",
      "validation_600dpi.tiff"
    )
  )

  open_tiff(
    tif_file,
    width,
    height
  )

  draw()
  dev.off()
}


# ============================================================
# PANEL A
# Count-level library depth and gene detection
# ============================================================

cat("Constructing Figure 2A...\n")

library_size <- colSums(
  raw_counts
)

detected_genes <- colSums(
  raw_counts > 0
)

genes_count10 <- colSums(
  raw_counts >= 10
)

qc <- data.frame(
  sampleID = colnames(raw_counts),
  Library_size = library_size,
  Library_size_million =
    library_size / 1e6,
  Detected_genes =
    detected_genes,
  Genes_count_ge10 =
    genes_count10,
  stringsAsFactors = FALSE
)

qc$group <- meta[
  qc$sampleID,
  "group"
]

qc$Stage <- unname(
  stage_labels[
    as.character(
      qc$group
    )
  ]
)

qc$Review_status <- meta[
  qc$sampleID,
  "Review_status"
]


p2a <- ggplot(
  qc,
  aes(
    x = Library_size_million,
    y = Detected_genes,
    colour = group
  )
) +

  geom_point(
    size = 2.7,
    alpha = 0.85
  ) +

  scale_colour_manual(
    values = stage_colors,
    breaks = groups,
    labels = unname(
      stage_labels[groups]
    ),
    name = "Stage"
  ) +

  labs(
    tag = "A",
    title = "Count-level library depth and gene detection",
    x = "Gene-assigned library size (million counts)",
    y = "Detected genes (count > 0)"
  ) +

  theme_pub +

  theme(
    legend.position = "right",
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 7.5
    ),
    legend.key.height = unit(
      0.34,
      "cm"
    )
  )


# ============================================================
# PANEL B
# 77-sample PCA
# ============================================================

cat("Constructing Figure 2B...\n")

pca <- prcomp(
  t(vst_mat),
  center = TRUE,
  scale. = FALSE
)

pca_var <- 100 *
  pca$sdev^2 /
  sum(
    pca$sdev^2
  )

pca_df <- data.frame(
  sampleID = rownames(
    pca$x
  ),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$group <- meta[
  pca_df$sampleID,
  "group"
]

pca_df$Stage <- unname(
  stage_labels[
    as.character(
      pca_df$group
    )
  ]
)


p2b <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    colour = group
  )
) +

  geom_point(
    size = 2.8,
    alpha = 0.88
  ) +

  scale_colour_manual(
    values = stage_colors,
    breaks = groups,
    labels = unname(
      stage_labels[groups]
    ),
    name = "Stage"
  ) +

  labs(
    tag = "B",
    title = "Global VST expression structure",
    x = paste0(
      "PC1 (",
      sprintf(
        "%.2f",
        pca_var[1]
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      sprintf(
        "%.2f",
        pca_var[2]
      ),
      "%)"
    )
  ) +

  theme_pub +

  theme(
    legend.position = "right",
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 7.5
    ),
    legend.key.height = unit(
      0.34,
      "cm"
    )
  )


# ============================================================
# 13. Validate PCA variance against frozen Step 03
# ============================================================

cat(
  "PCA variance: PC1 = ",
  sprintf(
    "%.2f",
    pca_var[1]
  ),
  "%; PC2 = ",
  sprintf(
    "%.2f",
    pca_var[2]
  ),
  "%\n",
  sep = ""
)


# ============================================================
# PANEL C
# 77 x 77 Spearman correlation heatmap
# ============================================================

cat("Constructing Figure 2C...\n")

spearman_mat <- cor(
  vst_mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)


# ------------------------------------------------------------
# Order samples by sampling stage
# ------------------------------------------------------------

sample_order <- unlist(
  lapply(
    groups,
    function(g) {

      x <- meta$sampleID[
        meta$group == g
      ]

      sort(x)
    }
  ),
  use.names = FALSE
)

spearman_ord <- spearman_mat[
  sample_order,
  sample_order,
  drop = FALSE
]


# ------------------------------------------------------------
# Stage boundaries and midpoint positions
# ------------------------------------------------------------

group_sizes <- sapply(
  groups,
  function(g) {
    sum(
      meta$group == g
    )
  }
)

group_ends <- cumsum(
  group_sizes
)

group_starts <- c(
  1,
  head(
    group_ends,
    -1
  ) + 1
)

group_midpoints <- (
  group_starts +
  group_ends
) / 2

boundary_positions <- head(
  group_ends,
  -1
) + 0.5


# ------------------------------------------------------------
# Long-format correlation matrix
# ------------------------------------------------------------

cor_long <- as.data.frame(
  as.table(
    spearman_ord
  ),
  stringsAsFactors = FALSE
)

colnames(
  cor_long
) <- c(
  "Sample_y",
  "Sample_x",
  "Spearman"
)

cor_long$x <- match(
  cor_long$Sample_x,
  sample_order
)

cor_long$y <- match(
  cor_long$Sample_y,
  sample_order
)

min_cor <- min(
  spearman_ord,
  na.rm = TRUE
)

legend_min <- max(
  0,
  floor(
    min_cor * 10
  ) / 10
)

cat(
  "Minimum sample Spearman correlation: ",
  sprintf(
    "%.4f",
    min_cor
  ),
  "\n",
  sep = ""
)


p2c <- ggplot(
  cor_long,
  aes(
    x = x,
    y = y,
    fill = Spearman
  )
) +

  geom_tile() +

  geom_vline(
    xintercept =
      boundary_positions,
    colour = "white",
    linewidth = 0.32
  ) +

  geom_hline(
    yintercept =
      boundary_positions,
    colour = "white",
    linewidth = 0.32
  ) +

  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#08306B",
    limits = c(
      legend_min,
      1
    ),
    oob = scales::squish,
    name = "Spearman\nρ"
  ) +

  scale_x_continuous(
    breaks =
      group_midpoints,
    labels =
      unname(
        stage_labels[
          groups
        ]
      ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    breaks =
      group_midpoints,
    labels =
      unname(
        stage_labels[
          groups
        ]
      ),
    expand = c(
      0,
      0
    )
  ) +

  coord_fixed() +

  labs(
    tag = "C",
    title = "Sample-to-sample VST Spearman correlation",
    x = "Sampling stage",
    y = "Sampling stage"
  ) +

  theme_minimal(
    base_size = 10
  ) +

  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    panel.grid = element_blank(),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 7
    ),
    axis.title = element_text(
      size = 9.5
    ),
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# PANEL D
# Sample-level expression review diagnostics
# ============================================================

cat("Constructing Figure 2D...\n")


# ------------------------------------------------------------
# Order samples exactly like correlation heatmap
# ------------------------------------------------------------

review_plot <- data.frame(
  sampleID = sample_order,
  stringsAsFactors = FALSE
)

review_plot$group <- meta[
  review_plot$sampleID,
  "group"
]

review_plot$Review_status <- meta[
  review_plot$sampleID,
  "Review_status"
]

review_plot$x <- seq_len(
  nrow(
    review_plot
  )
)

review_plot$Review_numeric <- c(
  "No flag" = 0,
  "Single-metric review" = 1,
  "Multi-metric review" = 2
)[
  as.character(
    review_plot$Review_status
  )
]


review_colors <- c(
  "No flag" = "#BDBDBD",
  "Single-metric review" = "#E69F00",
  "Multi-metric review" = "#D55E00"
)


p2d <- ggplot(
  review_plot,
  aes(
    x = x,
    y = Review_numeric,
    colour = Review_status
  )
) +

  geom_segment(
    aes(
      xend = x,
      y = 0,
      yend = Review_numeric
    ),
    colour = "grey82",
    linewidth = 0.45
  ) +

  geom_point(
    size = 2.6
  ) +

  geom_vline(
    xintercept =
      boundary_positions,
    colour = "grey80",
    linewidth = 0.35,
    linetype = "dashed"
  ) +

  scale_colour_manual(
    values = review_colors,
    drop = FALSE,
    name = "Review status"
  ) +

  scale_x_continuous(
    breaks =
      group_midpoints,
    labels =
      unname(
        stage_labels[
          groups
        ]
      ),
    expand = expansion(
      mult = c(
        0.01,
        0.01
      )
    )
  ) +

  scale_y_continuous(
    breaks = c(
      0,
      1,
      2
    ),
    labels = c(
      "No flag",
      "Single-metric",
      "Multi-metric"
    ),
    limits = c(
      -0.15,
      2.42
    ),
    expand = c(
      0,
      0
    )
  ) +

  annotate(
    "label",
    x = 58,
    y = 2.31,
    label = paste0(
      "No flag: ",
      unname(
        review_counts[
          "No flag"
        ]
      ),
      "   •   Single: ",
      unname(
        review_counts[
          "Single-metric review"
        ]
      ),
      "   •   Multi: ",
      unname(
        review_counts[
          "Multi-metric review"
        ]
      )
    ),
    hjust = 0.5,
    size = 2.75,
    fill = "white",
    colour = "grey30",
    linewidth = 0.25
  ) +

  annotate(
    "text",
    x = 38.5,
    y = -0.11,
    label =
      "All 77 samples retained in the primary analysis",
    size = 3.0,
    fontface = "italic",
    colour = "grey35"
  ) +

  labs(
    tag = "D",
    title = "Quantitative sample-level expression diagnostics",
    x = "Sampling stage",
    y = "Review level"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    axis.text.y = element_text(
      size = 8
    ),
    legend.position = "top",
    legend.title = element_text(
      size = 8.5
    ),
    legend.text = element_text(
      size = 8
    )
  )


# ============================================================
# 14. Save individual panels
# ============================================================

save_panel(
  p2a,
  "Figure_2A_count_level_QC",
  width = 6.5,
  height = 5
)

save_panel(
  p2b,
  "Figure_2B_VST_PCA",
  width = 6.5,
  height = 5
)

save_panel(
  p2c,
  "Figure_2C_sample_Spearman_heatmap",
  width = 6.5,
  height = 6
)

save_panel(
  p2d,
  "Figure_2D_sample_review_diagnostics",
  width = 6.5,
  height = 5.5
)


# ============================================================
# 15. Save final Figure 2
# ============================================================

save_final_figure(
  p2a,
  p2b,
  p2c,
  p2d,
  width = 13,
  height = 10
)


# ============================================================
# 16. Save source tables
# ============================================================

write.table(
  qc,
  file = file.path(
    panel_dir,
    "Figure_2A_count_level_QC_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  pca_df,
  file = file.path(
    panel_dir,
    "Figure_2B_PCA_coordinates.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    sampleID = rownames(
      spearman_ord
    ),
    spearman_ord,
    check.names = FALSE
  ),
  file = file.path(
    panel_dir,
    "Figure_2C_Spearman_matrix.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  review_plot,
  file = file.path(
    panel_dir,
    "Figure_2D_sample_review_plot_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 17. Save plot objects
# ============================================================

saveRDS(
  list(
    Figure2A = p2a,
    Figure2B = p2b,
    Figure2C = p2c,
    Figure2D = p2d
  ),
  file = file.path(
    outdir,
    "Figure_2_all_plot_objects.rds"
  )
)


# ============================================================
# 18. Summary file
# ============================================================

summary_table <- data.frame(
  Metric = c(
    "Samples",
    "Raw_count_genes",
    "VST_genes",
    "PCA_PC1_percent",
    "PCA_PC2_percent",
    "Minimum_sample_Spearman",
    "No_flag_samples",
    "Single_metric_review_samples",
    "Multi_metric_review_samples",
    "Samples_retained_primary_analysis"
  ),

  Value = c(
    ncol(
      vst_mat
    ),
    nrow(
      raw_counts
    ),
    nrow(
      vst_mat
    ),
    pca_var[1],
    pca_var[2],
    min_cor,
    unname(
      review_counts[
        "No flag"
      ]
    ),
    unname(
      review_counts[
        "Single-metric review"
      ]
    ),
    unname(
      review_counts[
        "Multi-metric review"
      ]
    ),
    77
  ),

  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    panel_dir,
    "Figure_2_technical_validation_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Session info
# ============================================================

sink(
  file.path(
    panel_dir,
    "Figure_2_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 20. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 2 COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples                     : ",
  ncol(vst_mat),
  "\n",
  sep = ""
)

cat(
  "Raw count genes             : ",
  nrow(raw_counts),
  "\n",
  sep = ""
)

cat(
  "VST genes                   : ",
  nrow(vst_mat),
  "\n",
  sep = ""
)

cat(
  "PC1 variance                : ",
  sprintf(
    "%.2f%%",
    pca_var[1]
  ),
  "\n",
  sep = ""
)

cat(
  "PC2 variance                : ",
  sprintf(
    "%.2f%%",
    pca_var[2]
  ),
  "\n",
  sep = ""
)

cat(
  "Minimum sample Spearman     : ",
  sprintf(
    "%.4f",
    min_cor
  ),
  "\n",
  sep = ""
)

cat(
  "No flag                     : ",
  unname(
    review_counts[
      "No flag"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Single-metric review        : ",
  unname(
    review_counts[
      "Single-metric review"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Multi-metric review         : ",
  unname(
    review_counts[
      "Multi-metric review"
    ]
  ),
  "\n",
  sep = ""
)

cat(
  "Samples retained            : 77\n"
)

cat("\nOutput:\n")

cat(
  file.path(
    outdir,
    "Figure_2_internal_technical_validation.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    outdir,
    "Figure_2_internal_technical_validation_600dpi.tiff"
  ),
  "\n",
  sep = ""
)

cat("\nPASS\n")
