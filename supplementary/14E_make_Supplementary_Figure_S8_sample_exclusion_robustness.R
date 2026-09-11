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
cat("14E SUPPLEMENTARY FIGURE S8\n")
cat("Sample-exclusion and leave-one-out robustness\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Output directories
# ============================================================

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S8_sample_exclusion_robustness"
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
# 1. Frozen Step10:
# Multi-sample exclusion sensitivity
#
# Primary_all77:
#   all final 77 samples
#
# Exclude_multi_metric:
#   RB_2, R0h_6, R6h_8, R10h_1, R18h_5
#
# Exclude_R2h4:
#   sequencing-review sample R2h_4
#
# Exclude_union:
#   union of the above six reviewed samples
# ============================================================

multi_global <- data.frame(
  Scenario = c(
    "Primary_all77",
    "Exclude_multi_metric",
    "Exclude_R2h4",
    "Exclude_union"
  ),

  Display = c(
    "Primary",
    "Exclude 5\nexpression-review",
    "Exclude\nR2h_4",
    "Exclude\nall 6"
  ),

  N_samples = c(
    77,
    72,
    76,
    71
  ),

  PERMANOVA_R2 = c(
    0.415491,
    0.466386,
    0.416445,
    0.468294
  ),

  PERMANOVA_P = c(
    1e-04,
    1e-04,
    1e-04,
    1e-04
  ),

  PERMDISP_P = c(
    0.7220,
    0.0575,
    0.7371,
    0.0646
  ),

  stringsAsFactors = FALSE
)


multi_cross <- data.frame(
  Scenario = c(
    "Primary_all77",
    "Exclude_multi_metric",
    "Exclude_R2h4",
    "Exclude_union"
  ),

  Display = c(
    "Primary",
    "Exclude 5\nexpression-review",
    "Exclude\nR2h_4",
    "Exclude\nall 6"
  ),

  Median_cross_cohort_Spearman = c(
    0.423107,
    0.469166,
    0.421734,
    0.469171
  ),

  Correlation_with_primary_70cell_matrix = c(
    1.000000,
    0.961368,
    0.999965,
    0.961648
  ),

  Mean_absolute_matrix_difference = c(
    0.000000,
    0.054870,
    0.000433,
    0.054600
  ),

  Maximum_absolute_matrix_difference = c(
    0.000000,
    0.140510,
    0.006060,
    0.140540
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 2. Frozen Step10B:
# Single-sample leave-one-out sensitivity
# ============================================================

loo <- data.frame(
  Removed_sample = c(
    "RB_2",
    "R0h_6",
    "R6h_8",
    "R10h_1",
    "R18h_5",
    "R2h_4"
  ),

  Review_type = c(
    "Expression review",
    "Expression review",
    "Expression review",
    "Expression review",
    "Expression review",
    "Sequencing review"
  ),

  N_samples = rep(
    76,
    6
  ),

  PERMANOVA_R2 = c(
    0.432401,
    0.426524,
    0.420773,
    0.417850,
    0.425804,
    0.416445
  ),

  PERMANOVA_P = rep(
    1e-04,
    6
  ),

  PERMDISP_P = c(
    0.3466,
    0.5132,
    0.7808,
    0.6362,
    0.5821,
    0.7371
  ),

  Minimum_genomewide_LFC_Spearman = c(
    0.89220,
    0.95450,
    0.93773,
    0.96378,
    0.90111,
    0.98733
  ),

  Cross_cohort_matrix_correlation = c(
    0.978515,
    0.999020,
    0.995241,
    0.999230,
    0.987333,
    0.999965
  ),

  Cross_cohort_matrix_MAE = c(
    0.060680,
    0.006310,
    0.005260,
    0.002350,
    0.005340,
    0.000433
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 3. Frozen primary references
# ============================================================

primary_permanova_r2 <- 0.415491
primary_permdisp_p <- 0.7220
primary_cross_median <- 0.423107


# ============================================================
# 4. Internal validation
# ============================================================

if (nrow(
  multi_global
) != 4) {
  stop(
    "Expected four Step10 scenarios."
  )
}

if (nrow(
  loo
) != 6) {
  stop(
    "Expected six Step10B leave-one-out scenarios."
  )
}

if (!all(
  loo$N_samples == 76
)) {
  stop(
    "All leave-one-out scenarios must contain 76 samples."
  )
}

if (abs(
  multi_global$PERMANOVA_R2[
    multi_global$Scenario ==
      "Primary_all77"
  ] -
    primary_permanova_r2
) > 1e-10) {
  stop(
    "Primary PERMANOVA reference mismatch."
  )
}

if (abs(
  multi_cross$
    Correlation_with_primary_70cell_matrix[
      multi_cross$Scenario ==
        "Primary_all77"
    ] -
    1
) > 1e-10) {
  stop(
    "Primary cross-cohort matrix reference mismatch."
  )
}

cat("Frozen Step10 / Step10B summaries loaded: PASS\n\n")


# ============================================================
# 5. Factors / plotting labels
# ============================================================

multi_levels <- c(
  "Primary_all77",
  "Exclude_multi_metric",
  "Exclude_R2h4",
  "Exclude_union"
)

multi_global$Scenario <- factor(
  multi_global$Scenario,
  levels = multi_levels
)

multi_cross$Scenario <- factor(
  multi_cross$Scenario,
  levels = multi_levels
)


loo_levels <- c(
  "RB_2",
  "R0h_6",
  "R6h_8",
  "R10h_1",
  "R18h_5",
  "R2h_4"
)

loo$Removed_sample <- factor(
  loo$Removed_sample,
  levels = loo_levels
)


# ============================================================
# 6. Publication theme
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
      size = 9.6
    ),

    axis.text = element_text(
      size = 8
    ),

    legend.title = element_text(
      size = 8.5
    ),

    legend.text = element_text(
      size = 8
    ),

    plot.margin = margin(
      8,
      8,
      10,
      8
    )
  )


# ============================================================
# 7. Panel A:
# Global temporal structure after multi-sample exclusion
# ============================================================

multi_global$Top_label <- paste0(
  "R² = ",
  sprintf(
    "%.3f",
    multi_global$PERMANOVA_R2
  ),
  "\nP < 0.001"
)


multi_axis_labels <- setNames(
  paste0(
    multi_global$Display,
    "\nN = ",
    multi_global$N_samples,
    "\nPERMDISP P = ",
    sprintf(
      "%.3f",
      multi_global$PERMDISP_P
    )
  ),
  as.character(
    multi_global$Scenario
  )
)


pA <- ggplot(
  multi_global,
  aes(
    x = Scenario,
    y = PERMANOVA_R2
  )
) +

  geom_hline(
    yintercept =
      primary_permanova_r2,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey55"
  ) +

  geom_col(
    width = 0.64,
    fill = "grey68",
    colour = "grey20",
    linewidth = 0.45
  ) +

  geom_text(
    aes(
      label = Top_label
    ),
    vjust = -0.35,
    size = 2.9,
    lineheight = 0.95
  ) +

  scale_x_discrete(
    labels = multi_axis_labels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      0.55
    ),
    breaks = seq(
      0,
      0.5,
      0.1
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "A",

    title =
      "Global temporal structure under multi-sample exclusion",

    subtitle =
      "Sampling-stage PERMANOVA in the full VST expression space",

    x = NULL,

    y =
      "PERMANOVA R² for sampling stage"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 7.3,
      face = "bold",
      lineheight = 1.03
    )
  )


# ============================================================
# 8. Panel B:
# Cross-cohort response matrix robustness
# after multi-sample exclusion
#
# Use categorical point display rather than XY scatter because
# Exclude_multi_metric and Exclude_union have nearly identical
# numerical values and overlap in the original XY layout.
# ============================================================

multi_cross_plot <- multi_cross[
  multi_cross$Scenario !=
    "Primary_all77",
  ,
  drop = FALSE
]

multi_cross_plot$Display_plot <- factor(
  multi_cross_plot$Display,
  levels = multi_cross_plot$Display
)

multi_cross_plot$Point_label <- paste0(
  "r = ",
  sprintf(
    "%.4f",
    multi_cross_plot$
      Correlation_with_primary_70cell_matrix
  ),
  "\nMAE = ",
  sprintf(
    "%.4f",
    multi_cross_plot$
      Mean_absolute_matrix_difference
  )
)


pB <- ggplot(
  multi_cross_plot,
  aes(
    x = Display_plot,
    y =
      Correlation_with_primary_70cell_matrix
  )
) +

  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey50"
  ) +

  geom_segment(
    aes(
      xend = Display_plot,
      y = 0.95,
      yend =
        Correlation_with_primary_70cell_matrix
    ),
    linewidth = 0.75,
    colour = "grey65"
  ) +

  geom_point(
    size = 4.2,
    shape = 21,
    fill = "white",
    colour = "grey20",
    stroke = 0.85
  ) +

  geom_text(
    aes(
      label = Point_label
    ),
    vjust = -0.75,
    size = 2.75,
    lineheight = 0.92,
    colour = "grey20"
  ) +

  scale_y_continuous(
    limits = c(
      0.95,
      1.006
    ),
    breaks = c(
      0.95,
      0.96,
      0.97,
      0.98,
      0.99,
      1.00
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "B",

    title =
      "Cross-cohort response structure after multi-sample exclusion",

    subtitle =
      "Similarity of each 10 x 7 Spearman concordance matrix to the primary analysis",

    x = NULL,

    y =
      "Correlation with primary concordance matrix"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 7.4,
      face = "bold",
      lineheight = 1.02
    )
  )


# 9. Panel C:
# Single-sample leave-one-out temporal structure
# ============================================================

loo$Top_label <- paste0(
  "R²=",
  sprintf(
    "%.3f",
    loo$PERMANOVA_R2
  )
)


loo_axis_labels <- setNames(
  paste0(
    as.character(
      loo$Removed_sample
    ),
    "\nPERMDISP P=",
    sprintf(
      "%.3f",
      loo$PERMDISP_P
    )
  ),
  as.character(
    loo$Removed_sample
  )
)


pC <- ggplot(
  loo,
  aes(
    x = Removed_sample,
    y = PERMANOVA_R2,
    fill = Review_type
  )
) +

  geom_hline(
    yintercept =
      primary_permanova_r2,
    linetype = "dashed",
    linewidth = 0.60,
    colour = "grey45"
  ) +

  geom_col(
    width = 0.63,
    colour = "grey20",
    linewidth = 0.4
  ) +

  geom_text(
    aes(
      label = Top_label
    ),
    vjust = -0.45,
    size = 2.75
  ) +

  scale_fill_manual(
    values = c(
      "Expression review" = "grey62",
      "Sequencing review" = "grey82"
    ),
    name = "Review category"
  ) +

  scale_x_discrete(
    labels = loo_axis_labels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      0.48
    ),
    breaks = seq(
      0,
      0.4,
      0.1
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "C",

    title =
      "Sampling-stage structure in single-sample leave-one-out analyses",

    subtitle = paste0(
      "Dashed line: primary 77-sample PERMANOVA R² = ",
      sprintf(
        "%.3f",
        primary_permanova_r2
      ),
      "; all leave-one-out PERMANOVA P < 0.001"
    ),

    x =
      "Removed reviewed sample",

    y =
      "PERMANOVA R² for sampling stage"
  ) +

  theme_pub +

  theme(
    legend.position = "top",

    axis.text.x = element_text(
      size = 7.1,
      lineheight = 1.02
    )
  )


# ============================================================
# 10. Panel D:
# Single-sample leave-one-out cross-cohort robustness
#
# Point/lollipop display is used instead of geom_col().
# Correlations are close to 1, so a zero-baseline bar chart is
# neither visually efficient nor compatible with a zoomed y-axis.
# ============================================================

loo$Cross_label <- paste0(
  sprintf(
    "%.4f",
    loo$Cross_cohort_matrix_correlation
  ),
  "\nMAE ",
  sprintf(
    "%.4f",
    loo$Cross_cohort_matrix_MAE
  )
)


pD <- ggplot(
  loo,
  aes(
    x = Removed_sample,
    y =
      Cross_cohort_matrix_correlation
  )
) +

  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey45"
  ) +

  geom_segment(
    aes(
      xend = Removed_sample,
      y = 0.965,
      yend =
        Cross_cohort_matrix_correlation
    ),
    linewidth = 0.75,
    colour = "grey68"
  ) +

  geom_point(
    aes(
      fill = Review_type
    ),
    shape = 21,
    size = 4.1,
    stroke = 0.75,
    colour = "grey20"
  ) +

  geom_text(
    aes(
      label = Cross_label
    ),
    vjust = -0.80,
    size = 2.55,
    lineheight = 0.90,
    colour = "grey20"
  ) +

  scale_fill_manual(
    values = c(
      "Expression review" = "grey62",
      "Sequencing review" = "grey85"
    ),
    guide = "none"
  ) +

  scale_y_continuous(
    limits = c(
      0.965,
      1.005
    ),
    breaks = c(
      0.97,
      0.98,
      0.99,
      1.00
    ),
    expand = c(
      0,
      0
    )
  ) +

  labs(
    tag = "D",

    title =
      "Cross-cohort concordance after single-sample removal",

    subtitle =
      "Correlation of each leave-one-out 10 x 7 response-concordance matrix with the primary matrix",

    x =
      "Removed reviewed sample",

    y =
      "Correlation with primary concordance matrix"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 7.6
    )
  )


# 11. TIFF helper
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
# 12. Save individual panels
# ============================================================

save_panel <- function(
    plot_object,
    stem,
    width = 6.4,
    height = 5.1) {

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

  print(
    plot_object
  )

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

  print(
    plot_object
  )

  dev.off()
}


save_panel(
  pA,
  "Supplementary_Figure_S8A_multi_sample_global_structure"
)

save_panel(
  pB,
  "Supplementary_Figure_S8B_multi_sample_cross_cohort_robustness"
)

save_panel(
  pC,
  "Supplementary_Figure_S8C_single_sample_PERMANOVA"
)

save_panel(
  pD,
  "Supplementary_Figure_S8D_single_sample_cross_cohort_robustness"
)


# ============================================================
# 13. Combined S8
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
# 14. Save final publication figure
# ============================================================

pdf_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S8_",
    "sample_exclusion_and_leave_one_out_robustness.pdf"
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
    "Supplementary_Figure_S8_",
    "sample_exclusion_and_leave_one_out_",
    "robustness_600dpi.tiff"
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
# 15. Save source-data tables
# ============================================================

write.table(
  multi_global,
  file = file.path(
    outdir,
    "Supplementary_Figure_S8A_multi_sample_global_structure.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  multi_cross,
  file = file.path(
    outdir,
    "Supplementary_Figure_S8B_multi_sample_cross_cohort_robustness.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  loo,
  file = file.path(
    outdir,
    "Supplementary_Figure_S8C_D_single_sample_leave_one_out.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 16. Publication summary table
# ============================================================

publication_summary <- data.frame(
  Analysis = c(
    "Primary",
    "Exclude five expression-review samples",
    "Exclude R2h_4 only",
    "Exclude union of six reviewed samples",
    "Single-sample leave-one-out minimum matrix correlation",
    "Single-sample leave-one-out maximum matrix correlation"
  ),

  Value = c(
    primary_permanova_r2,

    multi_global$PERMANOVA_R2[
      multi_global$Scenario ==
        "Exclude_multi_metric"
    ],

    multi_global$PERMANOVA_R2[
      multi_global$Scenario ==
        "Exclude_R2h4"
    ],

    multi_global$PERMANOVA_R2[
      multi_global$Scenario ==
        "Exclude_union"
    ],

    min(
      loo$Cross_cohort_matrix_correlation
    ),

    max(
      loo$Cross_cohort_matrix_correlation
    )
  ),

  stringsAsFactors = FALSE
)


write.table(
  publication_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S8_publication_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 17. Caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S8 | Robustness of temporal and cross-cohort ",
  "transcriptomic structure to exclusion of reviewed samples. ",
  "a, Sampling-stage PERMANOVA effect sizes in the complete 77-sample ",
  "2026 primary cohort and after exclusion of the five samples flagged ",
  "for multi-metric expression-level review, exclusion of the ",
  "sequencing-review sample R2h_4, or exclusion of the union of all six ",
  "reviewed samples. PERMDISP P values are shown beneath each scenario. ",
  "b, Stability of the complete 10-by-7 cross-cohort Spearman ",
  "response-concordance matrix after multi-sample exclusion, quantified ",
  "by the correlation with the primary 70-cell matrix and its mean ",
  "absolute difference. ",
  "c, Sampling-stage PERMANOVA effect sizes after leave-one-out removal ",
  "of each reviewed sample individually. The dashed horizontal line ",
  "indicates the primary 77-sample result. ",
  "d, Correlation between each leave-one-out cross-cohort response-",
  "concordance matrix and the primary 70-cell matrix; labels additionally ",
  "report the mean absolute matrix difference. ",
  "All leave-one-out analyses retained significant sampling-stage ",
  "structure, and no single reviewed sample dominated the cross-cohort ",
  "response-concordance pattern. Multi-sample exclusion changed some ",
  "quantitative estimates and exact differential-expression counts and ",
  "is therefore interpreted as a sensitivity analysis rather than ",
  "evidence that sample removal had no effect."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S8_caption.txt"
  )
)


# ============================================================
# 18. Interpretation note
# ============================================================

interpretation_text <- c(
  "Supplementary Figure S8 interpretation notes",
  "",
  "1. Review flags were diagnostic indicators and were not automatic",
  "   exclusion criteria.",
  "",
  "2. The primary analysis retains all final 77 samples.",
  "",
  "3. Excluding the five multi-metric expression-review samples changes",
  "   some exact DEG numbers and quantitative estimates; therefore the",
  "   correct interpretation is robustness of global structure, not",
  "   invariance of all results.",
  "",
  "4. R2h_4 removal has negligible influence on the global and",
  "   cross-cohort response structure.",
  "",
  "5. RB_2 produces the largest single-sample leave-one-out quantitative",
  "   effect because Baseline serves as the common reference for all",
  "   post-injury contrasts.",
  "",
  "6. Nevertheless, the cross-cohort 70-cell matrix remains highly",
  "   concordant with the primary analysis after RB_2 removal.",
  "",
  "7. No single reviewed sample provides evidence of dominating the",
  "   major temporal or cross-cohort transcriptomic organization.",
  "",
  "8. These analyses support retention of all 77 final samples in the",
  "   primary dataset and primary analyses."
)

writeLines(
  interpretation_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S8_interpretation_note.txt"
  )
)


# ============================================================
# 19. Plot objects / session info
# ============================================================

saveRDS(
  list(
    S8A = pA,
    S8B = pB,
    S8C = pC,
    S8D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S8_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "Supplementary_Figure_S8_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 20. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S8 COMPLETED\n")
cat("============================================================\n\n")


cat("Multi-sample exclusion global structure:\n")

print(
  multi_global,
  row.names = FALSE,
  digits = 6
)


cat("\nMulti-sample exclusion cross-cohort robustness:\n")

print(
  multi_cross,
  row.names = FALSE,
  digits = 6
)


cat("\nSingle-sample leave-one-out summary:\n")

print(
  loo,
  row.names = FALSE,
  digits = 6
)


cat("\nKey ranges:\n")

cat(
  "LOO PERMANOVA R2 range                 : ",
  sprintf(
    "%.6f - %.6f",
    min(
      loo$PERMANOVA_R2
    ),
    max(
      loo$PERMANOVA_R2
    )
  ),
  "\n",
  sep = ""
)

cat(
  "LOO PERMDISP P range                  : ",
  sprintf(
    "%.4f - %.4f",
    min(
      loo$PERMDISP_P
    ),
    max(
      loo$PERMDISP_P
    )
  ),
  "\n",
  sep = ""
)

cat(
  "LOO cross-cohort matrix correlation   : ",
  sprintf(
    "%.6f - %.6f",
    min(
      loo$Cross_cohort_matrix_correlation
    ),
    max(
      loo$Cross_cohort_matrix_correlation
    )
  ),
  "\n",
  sep = ""
)

cat(
  "LOO cross-cohort matrix MAE           : ",
  sprintf(
    "%.6f - %.6f",
    min(
      loo$Cross_cohort_matrix_MAE
    ),
    max(
      loo$Cross_cohort_matrix_MAE
    )
  ),
  "\n",
  sep = ""
)


cat("\nPrimary outputs:\n")

cat(
  pdf_file,
  "\n",
  sep = ""
)

cat(
  tiff_file,
  "\n\n",
  sep = ""
)

cat("PASS\n")
