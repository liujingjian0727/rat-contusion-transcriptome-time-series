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
cat("12D SUPPLEMENTARY FIGURE S2\n")
cat("Gene-body coverage sensitivity analysis\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Paths
# ============================================================

input_dir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage/",
  "12C_coverage_adjusted_PCA_sensitivity"
)

pca_file <- file.path(
  input_dir,
  "12C_PCA_coordinates_original_and_adjusted.tsv"
)

variance_file <- file.path(
  input_dir,
  "12C_PCA_variance_summary.tsv"
)

permanova_file <- file.path(
  input_dir,
  "12C_PERMANOVA_PERMDISP_summary.tsv"
)

distance_file <- file.path(
  input_dir,
  "12C_sample_distance_concordance.tsv"
)

centroid_file <- file.path(
  input_dir,
  "12C_stage_centroid_distance_concordance.tsv"
)

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S2_coverage_sensitivity"
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
  "Conditional",
  "Aggressive"
)

version_labels <- c(
  Original = "Original",
  Conditional = "Stage-preserving",
  Aggressive = "Coverage-only"
)


# ============================================================
# 2. Input validation
# ============================================================

required_files <- c(
  pca_file,
  variance_file,
  permanova_file,
  distance_file,
  centroid_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "Missing required input file(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

cat("All required 12C inputs found: PASS\n")


# ============================================================
# 3. Read metadata
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
    "Expected 77 metadata samples; observed ",
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
# 4. Read 12C PCA coordinates
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
    "Unexpected Version value in PCA table."
  )
}


# ============================================================
# 5. Verify each PCA version contains exactly final 77
# ============================================================

for (v in version_levels) {

  samples_v <- pca$sampleID[
    pca$Version == v
  ]

  if (length(samples_v) != 77) {
    stop(
      v,
      ": expected 77 samples; observed ",
      length(samples_v)
    )
  }

  if (!setequal(
    samples_v,
    final_samples
  )) {
    stop(
      v,
      ": sample set differs from final metadata."
    )
  }
}

cat(
  "PCA versions verified: 3 x 77 samples\n"
)


# ============================================================
# 6. Reattach canonical stage labels
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
    "Stage annotation failed."
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

if (!all(
  c(
    "Version",
    "PC1_percent",
    "PC2_percent"
  ) %in%
    colnames(variance)
)) {
  stop(
    "PCA variance table lacks required columns."
  )
}

variance$Version <- factor(
  variance$Version,
  levels = version_levels
)


# ============================================================
# 8. Read PERMANOVA / PERMDISP results
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
    "PERMANOVA summary lacks required columns."
  )
}

perm$Version <- factor(
  perm$Version,
  levels = version_levels
)

perm$Display_version <- factor(
  version_labels[
    as.character(
      perm$Version
    )
  ],
  levels = unname(
    version_labels[
      version_levels
    ]
  )
)


# ============================================================
# 9. Read concordance summaries
# ============================================================

distance_concordance <- read.delim(
  distance_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

centroid_concordance <- read.delim(
  centroid_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)


# ============================================================
# 10. Calculate PCA centroids
# ============================================================

make_centroids <- function(version_name) {

  tmp <- pca[
    pca$Version == version_name,
    ,
    drop = FALSE
  ]

  result <- do.call(
    rbind,
    lapply(
      groups,
      function(g) {

        x <- tmp[
          tmp$group == g,
          ,
          drop = FALSE
        ]

        data.frame(
          group = g,
          PC1 = mean(
            x$PC1
          ),
          PC2 = mean(
            x$PC2
          ),
          stringsAsFactors = FALSE
        )
      }
    )
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

centroids_original <-
  make_centroids(
    "Original"
  )

centroids_conditional <-
  make_centroids(
    "Conditional"
  )

centroids_aggressive <-
  make_centroids(
    "Aggressive"
  )


# ============================================================
# 11. Publication theme
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
      size = 8.5,
      colour = "grey35"
    ),

    axis.title = element_text(
      size = 10
    ),

    axis.text = element_text(
      size = 8.5
    ),

    legend.title = element_text(
      size = 8.5
    ),

    legend.text = element_text(
      size = 7.7
    ),

    legend.key.height = unit(
      0.35,
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
# 12. PCA panel helper
# ============================================================

make_pca_panel <- function(
    version_name,
    centroid_df,
    tag,
    title,
    subtitle) {

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

  ggplot(
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
      linewidth = 0.70,
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

    theme_pub +

    theme(
      legend.position = "right"
    )
}


# ============================================================
# 13. Panels A-C
# ============================================================

pA <- make_pca_panel(
  version_name = "Original",
  centroid_df = centroids_original,
  tag = "A",
  title = "Original VST expression space",
  subtitle =
    "No gene-body coverage adjustment"
)


pB <- make_pca_panel(
  version_name = "Conditional",
  centroid_df = centroids_conditional,
  tag = "B",
  title = "Stage-preserving adjustment",
  subtitle =
    "Coverage-associated effects removed conditional on sampling stage"
)


pC <- make_pca_panel(
  version_name = "Aggressive",
  centroid_df = centroids_aggressive,
  tag = "C",
  title = "Coverage-only adjustment",
  subtitle =
    "All linear signal associated with coverage metrics removed"
)


# ============================================================
# 14. Prepare Panel D annotations
# ============================================================

format_p <- function(p) {

  if (is.na(p)) {
    return("NA")
  }

  if (p < 0.001) {
    return("<0.001")
  }

  sprintf(
    "%.3f",
    p
  )
}


# Top-of-bar annotation:
# PERMANOVA effect size + significance
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


# Put PERMDISP information into x-axis labels
x_axis_labels <- setNames(
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
# 15. Panel D
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
    labels = x_axis_labels
  ) +

  scale_y_continuous(
    limits = c(
      0,
      0.50
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
    tag = "D",

    title =
      "Sampling-stage structure after coverage adjustment",

    subtitle =
      "PERMANOVA based on full VST Euclidean expression space",

    x = NULL,

    y =
      "PERMANOVA R² for sampling stage"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      size = 8.0,
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
# 16. Save individual panels
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
    plot_object,
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
  "Supplementary_Figure_S2A_original_PCA",
  6.3,
  5.1
)

save_panel(
  pB,
  "Supplementary_Figure_S2B_stage_preserving_PCA",
  6.3,
  5.1
)

save_panel(
  pC,
  "Supplementary_Figure_S2C_coverage_only_PCA",
  6.3,
  5.1
)

save_panel(
  pD,
  "Supplementary_Figure_S2D_PERMANOVA_sensitivity",
  6.3,
  5.1
)


# ============================================================
# 17. Combined Supplementary Figure S2
# ============================================================

draw_combined <- function() {

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
            1,
            1
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
  "Supplementary_Figure_S2_coverage_adjusted_sensitivity.pdf"
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
    "Supplementary_Figure_S2_",
    "coverage_adjusted_sensitivity_600dpi.tiff"
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
# 18. Publication summary table
# ============================================================

publication_summary <- data.frame(
  Expression_space = c(
    "Original",
    "Stage-preserving coverage adjustment",
    "Coverage-only adjustment"
  ),

  PC1_percent = variance$PC1_percent[
    match(
      version_levels,
      as.character(
        variance$Version
      )
    )
  ],

  PC2_percent = variance$PC2_percent[
    match(
      version_levels,
      as.character(
        variance$Version
      )
    )
  ],

  PERMANOVA_R2 = perm$PERMANOVA_R2[
    match(
      version_levels,
      as.character(
        perm$Version
      )
    )
  ],

  PERMANOVA_P = perm$PERMANOVA_P[
    match(
      version_levels,
      as.character(
        perm$Version
      )
    )
  ],

  PERMDISP_P = perm$PERMDISP_P[
    match(
      version_levels,
      as.character(
        perm$Version
      )
    )
  ],

  stringsAsFactors = FALSE
)


write.table(
  publication_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S2_publication_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Save concordance summary for supplementary source data
# ============================================================

write.table(
  distance_concordance,
  file = file.path(
    outdir,
    "Supplementary_Figure_S2_sample_distance_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  centroid_concordance,
  file = file.path(
    outdir,
    "Supplementary_Figure_S2_stage_centroid_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Create caption text file
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S2 | Sensitivity of global ",
  "transcriptomic structure to gene-body coverage variation. ",
  "a, Principal component analysis of the original ",
  "variance-stabilized expression matrix for the 77-sample ",
  "2026 primary cohort. ",
  "b, PCA after removal of coverage-associated expression ",
  "components conditional on sampling stage, thereby preserving ",
  "the stage-associated component of the expression structure. ",
  "c, PCA following a more conservative coverage-only adjustment ",
  "in which all linear expression components associated with the ",
  "two gene-body coverage metrics were removed without explicitly ",
  "preserving sampling-stage effects. Because gene-body coverage ",
  "characteristics were themselves associated with sampling stage, ",
  "this adjustment represents a sensitivity bound and may also ",
  "remove genuine stage-related biological variation. ",
  "d, PERMANOVA effect sizes for sampling stage in the original ",
  "and adjusted expression spaces. Stage-associated transcriptomic ",
  "structure remained significant after both adjustments, although ",
  "the PERMANOVA R-squared decreased under the coverage-only ",
  "adjustment."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S2_caption.txt"
  )
)


# ============================================================
# 21. Save R objects / session
# ============================================================

saveRDS(
  list(
    S2A = pA,
    S2B = pB,
    S2C = pC,
    S2D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S2_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "Supplementary_Figure_S2_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 22. Final validation
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

if (nrow(pca) !=
    77 * 3) {
  stop(
    "Expected 231 PCA coordinate rows; observed ",
    nrow(pca)
  )
}


# ============================================================
# 23. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S2 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Final cohort samples per PCA space : 77\n"
)

cat(
  "Expression spaces                 : 3\n"
)

cat(
  "R2h_1 excluded                    : TRUE\n"
)

cat(
  "R2h_8 retained                    : TRUE\n\n"
)

cat("Publication summary:\n")

print(
  publication_summary,
  row.names = FALSE,
  digits = 5
)

cat("\n")
cat("Sample-distance sensitivity:\n")

print(
  distance_concordance,
  row.names = FALSE,
  digits = 5
)

cat("\n")
cat("Stage-centroid sensitivity:\n")

print(
  centroid_concordance,
  row.names = FALSE,
  digits = 5
)

cat("\nOutput files:\n")

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
