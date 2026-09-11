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
cat("14A SUPPLEMENTARY FIGURE S4\n")
cat("Sensitivity of transcriptomic structure to fastp QC variation\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Input / output paths
# ============================================================

input_dir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13D_fastp_QC_adjusted_expression_space_sensitivity"
)

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

pca_file <- file.path(
  input_dir,
  "13D_PCA_coordinates_original_and_adjusted.tsv"
)

variance_file <- file.path(
  input_dir,
  "13D_PCA_variance_summary.tsv"
)

permanova_file <- file.path(
  input_dir,
  "13D_PERMANOVA_PERMDISP_summary.tsv"
)

sample_distance_file <- file.path(
  input_dir,
  "13D_sample_distance_concordance.tsv"
)

centroid_distance_file <- file.path(
  input_dir,
  "13D_stage_centroid_distance_concordance.tsv"
)

adjacent_rms_file <- file.path(
  input_dir,
  "13D_adjacent_stage_RMS_concordance.tsv"
)

qc_pc_file <- file.path(
  input_dir,
  "13D_QC_vs_PC1_to_PC5_after_adjustment.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S4_fastp_QC_sensitivity"
)

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
# 1. Fixed stage information
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

version_levels <- c(
  "Original",
  "Stage_preserving",
  "QC_only"
)

version_display <- c(
  Original = "Original",
  Stage_preserving = "Stage-preserving",
  QC_only = "QC-only"
)


# ============================================================
# 2. Validate all inputs
# ============================================================

required_files <- c(
  metadata_file,
  pca_file,
  variance_file,
  permanova_file,
  sample_distance_file,
  centroid_distance_file,
  adjacent_rms_file,
  qc_pc_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Missing required file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat("All required 13D outputs found: PASS\n")


# ============================================================
# 3. Read canonical metadata
# ============================================================

meta <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!all(
  c("sampleID", "group") %in%
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

if (anyDuplicated(meta$sampleID)) {
  stop(
    "Duplicated sampleID in metadata."
  )
}

meta$group <- factor(
  meta$group,
  levels = groups
)

final_samples <- meta$sampleID


# ============================================================
# 4. Read PCA coordinates
# ============================================================

pca <- read.delim(
  pca_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_pca_cols <- c(
  "sampleID",
  "PC1",
  "PC2",
  "Version"
)

if (!all(
  required_pca_cols %in%
    colnames(pca)
)) {
  stop(
    "PCA coordinate table lacks required columns."
  )
}

pca$Version <- factor(
  pca$Version,
  levels = version_levels
)

if (anyNA(pca$Version)) {
  stop(
    "Unexpected Version value in PCA coordinate table."
  )
}


# ============================================================
# 5. Verify 3 x 77 sample structure
# ============================================================

for (v in version_levels) {

  tmp_samples <- pca$sampleID[
    pca$Version == v
  ]

  if (length(tmp_samples) != 77) {
    stop(
      v,
      ": expected 77 samples; observed ",
      length(tmp_samples)
    )
  }

  if (!setequal(
    tmp_samples,
    final_samples
  )) {
    stop(
      v,
      ": sample set differs from final metadata."
    )
  }
}

if (nrow(pca) != 231) {
  stop(
    "Expected 231 PCA coordinate rows; observed ",
    nrow(pca)
  )
}

cat(
  "PCA structure verified : 3 expression spaces x 77 samples\n"
)


# ============================================================
# 6. Attach canonical group annotation
# ============================================================

pca$group <- meta$group[
  match(
    pca$sampleID,
    meta$sampleID
  )
]

pca$Stage <- unname(
  stage_labels[
    as.character(
      pca$group
    )
  ]
)

if (anyNA(pca$group)) {
  stop(
    "Failed to attach group annotation."
  )
}


# ============================================================
# 7. Read PCA variance
# ============================================================

variance <- read.delim(
  variance_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_variance_cols <- c(
  "Version",
  "PC1_percent",
  "PC2_percent"
)

if (!all(
  required_variance_cols %in%
    colnames(variance)
)) {
  stop(
    "Variance summary lacks required columns."
  )
}

variance$Version <- factor(
  variance$Version,
  levels = version_levels
)

if (!setequal(
  as.character(variance$Version),
  version_levels
)) {
  stop(
    "Unexpected versions in PCA variance summary."
  )
}


# ============================================================
# 8. Read PERMANOVA / PERMDISP
# ============================================================

perm <- read.delim(
  permanova_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_perm_cols <- c(
  "Version",
  "PERMANOVA_R2",
  "PERMANOVA_F",
  "PERMANOVA_P",
  "PERMDISP_F",
  "PERMDISP_P"
)

if (!all(
  required_perm_cols %in%
    colnames(perm)
)) {
  stop(
    "PERMANOVA/PERMDISP table lacks required columns."
  )
}

perm$Version <- factor(
  perm$Version,
  levels = version_levels
)

if (nrow(perm) != 3) {
  stop(
    "Expected three PERMANOVA rows."
  )
}


# ============================================================
# 9. Read supporting sensitivity summaries
# ============================================================

sample_distance <- read.delim(
  sample_distance_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

centroid_distance <- read.delim(
  centroid_distance_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

adjacent_rms <- read.delim(
  adjacent_rms_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

qc_pc <- read.delim(
  qc_pc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 10. Reorder rows consistently
# ============================================================

variance <- variance[
  match(
    version_levels,
    as.character(
      variance$Version
    )
  ),
  ,
  drop = FALSE
]

perm <- perm[
  match(
    version_levels,
    as.character(
      perm$Version
    )
  ),
  ,
  drop = FALSE
]


# ============================================================
# 11. Sanity check against frozen 13D values
# ============================================================

expected_original_r2 <- 0.41549
expected_qc_only_r2 <- 0.27311

original_r2 <- perm$PERMANOVA_R2[
  perm$Version == "Original"
]

qc_only_r2 <- perm$PERMANOVA_R2[
  perm$Version == "QC_only"
]

if (abs(
  original_r2 -
    expected_original_r2
) > 0.001) {
  warning(
    "Original PERMANOVA R2 differs from expected 13D value."
  )
}

if (abs(
  qc_only_r2 -
    expected_qc_only_r2
) > 0.001) {
  warning(
    "QC-only PERMANOVA R2 differs from expected 13D value."
  )
}


# ============================================================
# 12. Calculate PCA-stage centroids
# ============================================================

make_centroids <- function(
    version_name) {

  tmp <- pca[
    pca$Version == version_name,
    ,
    drop = FALSE
  ]

  out_list <- vector(
    "list",
    length(groups)
  )

  for (i in seq_along(groups)) {

    g <- groups[i]

    x <- tmp[
      tmp$group == g,
      ,
      drop = FALSE
    ]

    out_list[i] <- list(
      data.frame(
        group = g,
        PC1 = mean(x$PC1),
        PC2 = mean(x$PC2),
        stringsAsFactors = FALSE
      )
    )
  }

  result <- do.call(
    rbind,
    out_list
  )

  result$group <- factor(
    result$group,
    levels = groups
  )

  result$Stage <- unname(
    stage_labels[
      as.character(
        result$group
      )
    ]
  )

  result
}


cent_original <- make_centroids(
  "Original"
)

cent_stage_preserving <- make_centroids(
  "Stage_preserving"
)

cent_qc_only <- make_centroids(
  "QC_only"
)


# ============================================================
# 13. Publication theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 10.5
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),

    plot.title = element_text(
      face = "bold",
      size = 11
    ),

    plot.subtitle = element_text(
      size = 8.3,
      colour = "grey35"
    ),

    axis.title = element_text(
      size = 9.7
    ),

    axis.text = element_text(
      size = 8.3
    ),

    legend.title = element_text(
      size = 8.4
    ),

    legend.text = element_text(
      size = 7.5
    ),

    legend.key.height = unit(
      0.34,
      "cm"
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 14. PCA plotting function
# ============================================================

make_pca_panel <- function(
    version_name,
    centroid_df,
    tag,
    title,
    subtitle,
    show_legend = FALSE) {

  tmp <- pca[
    pca$Version == version_name,
    ,
    drop = FALSE
  ]

  var_row <- variance[
    variance$Version == version_name,
    ,
    drop = FALSE
  ]

  if (nrow(var_row) != 1) {
    stop(
      "Could not identify variance row for ",
      version_name
    )
  }


  out <- ggplot(
    tmp,
    aes(
      x = PC1,
      y = PC2,
      colour = group
    )
  ) +

    geom_path(
      data = centroid_df,
      aes(
        x = PC1,
        y = PC2,
        group = 1
      ),
      inherit.aes = FALSE,
      colour = "grey35",
      linewidth = 0.72,
      arrow = arrow(
        length = unit(
          0.11,
          "cm"
        ),
        type = "closed"
      )
    ) +

    geom_point(
      size = 2.15,
      alpha = 0.72
    ) +

    geom_point(
      data = centroid_df,
      aes(
        x = PC1,
        y = PC2,
        fill = group
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.7,
      stroke = 0.65,
      colour = "black"
    ) +

    scale_colour_manual(
      values = stage_colors,
      breaks = groups,
      labels = unname(
        stage_labels[
          groups
        ]
      ),
      name = "Sampling stage"
    ) +

    scale_fill_manual(
      values = stage_colors,
      guide = "none"
    ) +

    labs(
      tag = tag,
      title = title,
      subtitle = subtitle,

      x = paste0(
        "PC1 (",
        sprintf(
          "%.2f",
          var_row$PC1_percent
        ),
        "%)"
      ),

      y = paste0(
        "PC2 (",
        sprintf(
          "%.2f",
          var_row$PC2_percent
        ),
        "%)"
      )
    ) +

    theme_pub


  if (show_legend) {

    out <- out +
      theme(
        legend.position = "right"
      )

  } else {

    out <- out +
      theme(
        legend.position = "none"
      )
  }

  out
}


# ============================================================
# 15. Panels A-C
# ============================================================

pA <- make_pca_panel(
  version_name =
    "Original",

  centroid_df =
    cent_original,

  tag = "A",

  title =
    "Original VST expression space",

  subtitle =
    "No fastp-derived QC adjustment",

  show_legend = FALSE
)


pB <- make_pca_panel(
  version_name =
    "Stage_preserving",

  centroid_df =
    cent_stage_preserving,

  tag = "B",

  title =
    "Stage-preserving fastp-QC adjustment",

  subtitle =
    "Q30-, GC- and duplication-associated effects removed conditional on stage",

  show_legend = TRUE
)


pC <- make_pca_panel(
  version_name =
    "QC_only",

  centroid_df =
    cent_qc_only,

  tag = "C",

  title =
    "QC-only adjustment",

  subtitle =
    "All linear signal associated with Q30, GC and duplication removed",

  show_legend = FALSE
)


# ============================================================
# 16. Prepare Panel D
# ============================================================

perm$Display_version <- factor(
  unname(
    version_display[
      as.character(
        perm$Version
      )
    ]
  ),
  levels = unname(
    version_display[
      version_levels
    ]
  )
)


format_p <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf(
    "=%.3f",
    p
  )
}


perm$Top_annotation <- paste0(
  "R² = ",
  sprintf(
    "%.3f",
    perm$PERMANOVA_R2
  ),
  "\nPERMANOVA P ",
  vapply(
    perm$PERMANOVA_P,
    format_p,
    character(1)
  )
)


axis_labels <- setNames(
  paste0(
    as.character(
      perm$Display_version
    ),
    "\nPERMDISP P = ",
    sprintf(
      "%.3f",
      perm$PERMDISP_P
    )
  ),
  as.character(
    perm$Display_version
  )
)


# ============================================================
# 17. Panel D
# ============================================================

pD <- ggplot(
  perm,
  aes(
    x = Display_version,
    y = PERMANOVA_R2
  )
) +

  geom_col(
    width = 0.62,
    fill = "grey67",
    colour = "grey20",
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = Top_annotation
    ),
    vjust = -0.35,
    size = 3.0,
    lineheight = 0.95,
    colour = "grey15"
  ) +

  scale_x_discrete(
    labels = axis_labels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      max(
        perm$PERMANOVA_R2
      ) * 1.22
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "D",

    title =
      "Sampling-stage structure after fastp-QC adjustment",

    subtitle =
      "PERMANOVA based on full VST Euclidean expression space",

    x = NULL,

    y =
      "PERMANOVA R² for sampling stage"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 7.8,
      face = "bold",
      lineheight = 1.05
    ),

    plot.margin = margin(
      8,
      10,
      14,
      8
    )
  )


# ============================================================
# 18. TIFF helper
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


# ============================================================
# 19. Save individual panels
# ============================================================

save_panel <- function(
    plot_object,
    stem,
    width = 6.3,
    height = 5.1) {

  pdf_name <- file.path(
    panel_dir,
    paste0(
      stem,
      ".pdf"
    )
  )

  grDevices::cairo_pdf(
    filename = pdf_name,
    width = width,
    height = height
  )

  print(
    plot_object
  )

  dev.off()


  tif_name <- file.path(
    panel_dir,
    paste0(
      stem,
      "_600dpi.tiff"
    )
  )

  open_tiff(
    tif_name,
    width,
    height
  )

  print(
    plot_object
  )

  dev.off()
}


save_panel(
  pA,
  "Supplementary_Figure_S4A_original_PCA"
)

save_panel(
  pB,
  "Supplementary_Figure_S4B_stage_preserving_PCA"
)

save_panel(
  pC,
  "Supplementary_Figure_S4C_QC_only_PCA"
)

save_panel(
  pD,
  "Supplementary_Figure_S4D_PERMANOVA_sensitivity"
)


# ============================================================
# 20. Draw combined figure
# ============================================================

draw_combined <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2,
        widths = unit(
          c(1, 1),
          "null"
        ),
        heights = unit(
          c(1, 1),
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


# ============================================================
# 21. Save combined publication figure
# ============================================================

pdf_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S4_",
    "fastp_QC_adjusted_expression_space_sensitivity.pdf"
  )
)

grDevices::cairo_pdf(
  filename = pdf_file,
  width = 13,
  height = 10
)

draw_combined()

dev.off()


tiff_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S4_",
    "fastp_QC_adjusted_expression_space_",
    "sensitivity_600dpi.tiff"
  )
)

open_tiff(
  tiff_file,
  width = 13,
  height = 10
)

draw_combined()

dev.off()


# ============================================================
# 22. Publication summary
# ============================================================

publication_summary <- data.frame(
  Expression_space = c(
    "Original",
    "Stage-preserving fastp-QC adjustment",
    "QC-only adjustment"
  ),

  PC1_percent =
    variance$PC1_percent,

  PC2_percent =
    variance$PC2_percent,

  PERMANOVA_R2 =
    perm$PERMANOVA_R2,

  PERMANOVA_F =
    perm$PERMANOVA_F,

  PERMANOVA_P =
    perm$PERMANOVA_P,

  PERMDISP_F =
    perm$PERMDISP_F,

  PERMDISP_P =
    perm$PERMDISP_P,

  stringsAsFactors = FALSE
)


write.table(
  publication_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_publication_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Source-data summary
# ============================================================

write.table(
  pca,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_PCA_source_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sample_distance,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_sample_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  centroid_distance,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_stage_centroid_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  adjacent_rms,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_adjacent_stage_RMS_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  qc_pc,
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_QC_PC_associations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 24. Final caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S4 | Sensitivity of global transcriptomic ",
  "structure to fastp-derived sequencing-quality variation. ",
  "a, Principal component analysis of the original variance-stabilized ",
  "expression matrix for the 77-sample 2026 primary cohort. ",
  "b, PCA after removal of expression components associated with ",
  "Q30 base-call quality, GC content and estimated read duplication ",
  "conditional on sampling stage. This stage-preserving adjustment ",
  "retains stage-associated expression effects while removing residual ",
  "QC-associated variation. ",
  "c, PCA following a more conservative QC-only adjustment in which ",
  "all linear expression components associated with Q30, GC content ",
  "and duplication were removed without explicitly preserving sampling ",
  "stage. Because these QC metrics were themselves associated with ",
  "sampling stage, the QC-only adjustment represents a sensitivity bound ",
  "and may also remove genuine stage-related biological variation. ",
  "d, PERMANOVA effect sizes for sampling stage in the original and ",
  "adjusted full-expression spaces. Sampling-stage structure remained ",
  "significant in all three expression spaces (all PERMANOVA P < 0.001), ",
  "with R-squared values of 0.415, 0.539 and 0.273 for the original, ",
  "stage-preserving and QC-only expression spaces, respectively. ",
  "Multivariate dispersion remained nonsignificant in all analyses ",
  "(PERMDISP P = 0.722, 0.819 and 0.521, respectively). ",
  "PCA was recalculated independently for each expression matrix; ",
  "therefore, the orientation and sign of individual principal-component ",
  "axes are not directly comparable among panels."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S4_caption.txt"
  )
)


# ============================================================
# 25. Technical interpretation note
# ============================================================

interpretation_text <- c(
  "Technical interpretation for manuscript preparation:",
  "",
  "1. fastp-derived Q30, GC content and duplication characteristics",
  "   contribute measurably to global transcriptomic structure.",
  "",
  "2. The stage-preserving adjustment deliberately retains sampling-stage",
  "   effects and therefore should not be interpreted as an independent",
  "   estimate of biological effect size.",
  "",
  "3. The QC-only adjustment is the more conservative stress test.",
  "   Sampling stage remained significant after this adjustment",
  "   (PERMANOVA R2 = 0.273, P < 0.001).",
  "",
  "4. Because QC characteristics covary with sampling stage, QC-only",
  "   adjustment may remove both technical and genuine biological signal.",
  "",
  "5. Therefore, adjusted matrices are sensitivity-analysis products and",
  "   should not replace the original VST matrix in the primary analysis.",
  "",
  "6. Fine-scale stage geometry and adjacent-stage distance rankings should",
  "   remain descriptive rather than mechanistic interpretations."
)

writeLines(
  interpretation_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S4_interpretation_note.txt"
  )
)


# ============================================================
# 26. Save plot objects / session
# ============================================================

saveRDS(
  list(
    S4A = pA,
    S4B = pB,
    S4C = pC,
    S4D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S4_plot_objects.rds"
  )
)

sink(
  file.path(
    outdir,
    "Supplementary_Figure_S4_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 27. Final checks
# ============================================================

if ("R2h_1" %in%
    pca$sampleID) {
  stop(
    "R2h_1 unexpectedly present."
  )
}

if (!("R2h_8" %in%
      pca$sampleID)) {
  stop(
    "R2h_8 unexpectedly absent."
  )
}

if (nrow(pca) != 231) {
  stop(
    "Expected 231 PCA rows."
  )
}


# ============================================================
# 28. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S4 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Final cohort samples per expression space : 77\n"
)

cat(
  "Expression spaces                       : 3\n"
)

cat(
  "R2h_1 excluded                          : TRUE\n"
)

cat(
  "R2h_8 retained                          : TRUE\n\n"
)

cat("Publication summary:\n")

print(
  publication_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n")
cat("Sample-distance concordance:\n")

print(
  sample_distance,
  row.names = FALSE,
  digits = 5
)

cat("\n")
cat("Stage-centroid concordance:\n")

print(
  centroid_distance,
  row.names = FALSE,
  digits = 5
)

cat("\n")
cat("Adjacent-stage RMS concordance:\n")

print(
  adjacent_rms,
  row.names = FALSE,
  digits = 5
)

cat("\nPrimary outputs:\n")

cat(
  pdf_file,
  "\n",
  sep = ""
)

cat(
  tiff_file,
  "\n",
  sep = ""
)

cat("\nPASS\n")
