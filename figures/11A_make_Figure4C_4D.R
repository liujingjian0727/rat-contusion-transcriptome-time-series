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
cat("11A MAKE FIGURE 4C AND FIGURE 4D\n")
cat("Scientific Data final figure preparation\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed settings
# ============================================================

injury_groups <- c(
  "R0h", "R1h", "R2h", "R6h", "R10h",
  "R14h", "R18h", "R36h", "R60h", "R72h"
)

stage_labels <- c(
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

EXPECTED_GENES <- 476
Z_LIMIT <- 2.5


# ============================================================
# 1. Input files
# ============================================================

expr_file <- paste0(
  "07_2026_TPM_tau_SPM_peak_stage_result/",
  "07_2026_group_median_log2TPMplus1_all11stages.tsv"
)

annot_file <- paste0(
  "07_2026_TPM_tau_SPM_peak_stage_result/",
  "07_2026_high_confidence_dynamic_",
  "injury_time_specific_tau085_spm050.tsv"
)

count_file <- paste0(
  "07_2026_TPM_tau_SPM_peak_stage_result/",
  "07_2026_high_confidence_dynamic_",
  "injury_time_specific_gene_numbers.tsv"
)

outdir <- "final_submission_figures"
panel_dir <- file.path(
  outdir,
  "panels"
)

dir.create(
  panel_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

required_files <- c(
  expr_file,
  annot_file,
  count_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat("Input files found: PASS\n")


# ============================================================
# 2. Helper for publication TIFF/PDF
# ============================================================

save_single_plot <- function(
    p,
    stem,
    width,
    height) {

  pdf_file <- file.path(
    panel_dir,
    paste0(stem, ".pdf")
  )

  grDevices::cairo_pdf(
    filename = pdf_file,
    width = width,
    height = height
  )

  print(p)
  dev.off()


  tif_file <- file.path(
    panel_dir,
    paste0(stem, "_600dpi.tiff")
  )

  if (requireNamespace(
    "ragg",
    quietly = TRUE
  )) {

    ragg::agg_tiff(
      filename = tif_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

  } else {

    grDevices::tiff(
      filename = tif_file,
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


save_4CD <- function(
    pC,
    pD,
    stem,
    width = 12,
    height = 7) {

  draw_layout <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 1,
          ncol = 4,
          widths = unit(
            c(1, 1, 1, 1.15),
            "null"
          )
        )
      )
    )

    print(
      pC,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1:3
      )
    )

    print(
      pD,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 4
      )
    )
  }

  pdf_file <- file.path(
    outdir,
    paste0(stem, ".pdf")
  )

  grDevices::cairo_pdf(
    pdf_file,
    width = width,
    height = height
  )

  draw_layout()
  dev.off()


  tif_file <- file.path(
    outdir,
    paste0(stem, "_600dpi.tiff")
  )

  if (requireNamespace(
    "ragg",
    quietly = TRUE
  )) {

    ragg::agg_tiff(
      filename = tif_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      compression = "lzw"
    )

  } else {

    grDevices::tiff(
      filename = tif_file,
      width = width,
      height = height,
      units = "in",
      res = 600,
      type = "cairo",
      compression = "lzw"
    )
  }

  draw_layout()
  dev.off()
}


# ============================================================
# 3. Read expression matrix
# ============================================================

cat("Reading median log2(TPM+1)...\n")

expr <- read.delim(
  expr_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

colnames(expr)[1] <- "GeneID"

if (anyDuplicated(expr$GeneID)) {
  stop("Duplicated GeneID in expression matrix.")
}

missing_stage_cols <- setdiff(
  injury_groups,
  colnames(expr)
)

if (length(missing_stage_cols) > 0) {
  stop(
    "Missing stage columns: ",
    paste(
      missing_stage_cols,
      collapse = ", "
    )
  )
}


# ============================================================
# 4. Read 476-gene annotation
# ============================================================

cat("Reading high-confidence temporal-specific genes...\n")

annot <- read.delim(
  annot_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!"GeneID" %in% colnames(annot)) {
  stop("GeneID missing from annotation file.")
}

if (nrow(annot) != EXPECTED_GENES) {
  warning(
    "Expected ",
    EXPECTED_GENES,
    " genes but observed ",
    nrow(annot)
  )
}

peak_col <- if (
  "Peak_stage_injury10" %in%
  colnames(annot)
) {
  "Peak_stage_injury10"
} else {
  stop(
    "Peak_stage_injury10 column not found."
  )
}

spm_col <- if (
  "Peak_SPM_injury10" %in%
  colnames(annot)
) {
  "Peak_SPM_injury10"
} else {
  NULL
}


# ============================================================
# 5. Align 476 genes
# ============================================================

idx <- match(
  annot$GeneID,
  expr$GeneID
)

if (anyNA(idx)) {
  stop(
    "Some 476 genes were not found ",
    "in expression matrix."
  )
}

mat <- as.matrix(
  expr[
    idx,
    injury_groups,
    drop = FALSE
  ]
)

storage.mode(mat) <- "numeric"

rownames(mat) <- annot$GeneID


# ============================================================
# 6. Row Z-score
# ============================================================

row_mean <- rowMeans(mat)

row_sd <- apply(
  mat,
  1,
  sd
)

z_mat <- sweep(
  mat,
  1,
  row_mean,
  "-"
)

z_mat <- sweep(
  z_mat,
  1,
  row_sd,
  "/"
)

z_mat[
  !is.finite(z_mat)
] <- 0

z_mat <- pmax(
  pmin(
    z_mat,
    Z_LIMIT
  ),
  -Z_LIMIT
)


# ============================================================
# 7. Sort genes by peak stage
# ============================================================

annot$Peak_factor <- factor(
  annot[[peak_col]],
  levels = injury_groups
)

if (!is.null(spm_col)) {

  order_index <- order(
    annot$Peak_factor,
    -annot[[spm_col]],
    annot$GeneID
  )

} else {

  order_index <- order(
    annot$Peak_factor,
    annot$GeneID
  )
}

annot_ord <- annot[
  order_index,
  ,
  drop = FALSE
]

z_mat <- z_mat[
  annot_ord$GeneID,
  ,
  drop = FALSE
]

ordered_genes <- rownames(
  z_mat
)


# ============================================================
# 8. Convert heatmap to long table
# ============================================================

heat_df <- as.data.frame(
  as.table(z_mat),
  stringsAsFactors = FALSE
)

colnames(heat_df) <- c(
  "GeneID",
  "Group",
  "Z"
)

heat_df$Gene_index <- match(
  heat_df$GeneID,
  ordered_genes
)

heat_df$Group <- factor(
  heat_df$Group,
  levels = injury_groups,
  labels = unname(
    stage_labels[injury_groups]
  )
)


# ============================================================
# 9. Peak-stage block boundaries
# ============================================================

peak_counts <- table(
  factor(
    annot_ord[[peak_col]],
    levels = injury_groups
  )
)

boundaries <- cumsum(
  as.numeric(
    peak_counts
  )
)

boundaries <- boundaries[
  -length(boundaries)
] + 0.5


# ============================================================
# 10. Figure 4C
# ============================================================

p4c <- ggplot(
  heat_df,
  aes(
    x = Group,
    y = Gene_index,
    fill = Z
  )
) +
  geom_tile() +
  geom_hline(
    yintercept = boundaries,
    linewidth = 0.18,
    colour = "grey75"
  ) +
  scale_y_reverse(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -Z_LIMIT,
      Z_LIMIT
    ),
    oob = scales::squish,
    name = "Row Z-score"
  ) +
  labs(
    tag = "C",
    x = "Post-injury sampling stage",
    y = "Genes ordered by peak stage",
    title = paste0(
      nrow(annot_ord),
      " high-confidence injury-stage-specific genes"
    )
  ) +
  theme_classic(
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
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# 11. Figure 4D
# ============================================================

counts <- read.delim(
  count_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_cols <- c(
  "Peak_stage",
  "Gene_number"
)

if (!all(
  required_cols %in%
  colnames(counts)
)) {
  stop(
    "Expected columns missing from ",
    count_file
  )
}

counts$Peak_stage <- factor(
  counts$Peak_stage,
  levels = injury_groups
)

counts$Stage_label <- unname(
  stage_labels[
    as.character(
      counts$Peak_stage
    )
  ]
)

counts$Stage_label <- factor(
  counts$Stage_label,
  levels = unname(
    stage_labels[injury_groups]
  )
)

p4d <- ggplot(
  counts,
  aes(
    x = Stage_label,
    y = Gene_number
  )
) +
  geom_col(
    width = 0.72,
    fill = "#4C78A8"
  ) +
  geom_text(
    aes(
      label = Gene_number
    ),
    vjust = -0.35,
    size = 3
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  labs(
    tag = "D",
    x = "Peak stage",
    y = "Number of genes",
    title = "Peak-stage distribution"
  ) +
  theme_classic(
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
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# 12. Save individual panels
# ============================================================

save_single_plot(
  p4c,
  "Figure_4C_temporal_heatmap",
  width = 8.5,
  height = 7
)

save_single_plot(
  p4d,
  "Figure_4D_peak_stage_distribution",
  width = 5,
  height = 6.5
)


# ============================================================
# 13. Save plot objects for final Figure 4
# ============================================================

saveRDS(
  p4c,
  file.path(
    panel_dir,
    "Figure_4C_plot.rds"
  )
)

saveRDS(
  p4d,
  file.path(
    panel_dir,
    "Figure_4D_plot.rds"
  )
)


# ============================================================
# 14. Save combined 4C + 4D preview
# ============================================================

save_4CD(
  p4c,
  p4d,
  "Figure_4CD_temporal_specificity_preview",
  width = 12,
  height = 7
)


# ============================================================
# 15. Save ordering information
# ============================================================

ordering_table <- data.frame(
  Gene_order = seq_along(
    ordered_genes
  ),
  GeneID = ordered_genes,
  Peak_stage = annot_ord[[peak_col]],
  stringsAsFactors = FALSE
)

if (!is.null(spm_col)) {
  ordering_table$Peak_SPM <-
    annot_ord[[spm_col]]
}

write.table(
  ordering_table,
  file = file.path(
    panel_dir,
    "Figure_4C_gene_order.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


cat("\n")
cat("============================================================\n")
cat("FIGURE 4C/4D COMPLETED\n")
cat("============================================================\n")
cat("Genes in heatmap : ", nrow(annot_ord), "\n", sep = "")
cat("Expected genes   : ", EXPECTED_GENES, "\n", sep = "")
cat("Peak-stage counts:\n")
print(
  counts[
    ,
    c(
      "Stage_label",
      "Gene_number"
    )
  ],
  row.names = FALSE
)
cat("\nOutput directory:\n")
cat(outdir, "\n")
cat("\nPASS\n")
