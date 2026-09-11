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
cat("11B FINAL FIGURES 3, 4 AND 5\n")
cat("Scientific Data final main figures\n")
cat("============================================================\n\n")


# ============================================================
# 0. Global stage definitions
# ============================================================

groups11 <- c(
  "Baseline",
  "R0h", "R1h", "R2h", "R6h", "R10h",
  "R14h", "R18h", "R36h", "R60h", "R72h"
)

labels11 <- c(
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

groups26 <- groups11[
  groups11 != "Baseline"
]

labels26 <- labels11[
  groups26
]

groups21 <- c(
  "R4h", "R8h", "R12h",
  "R16h", "R20h", "R24h", "R48h"
)

labels21 <- c(
  R4h  = "4 h",
  R8h  = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)

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
# 1. Helper functions
# ============================================================

read_tsv <- function(file) {

  if (!file.exists(file)) {
    stop("Missing file: ", file)
  }

  read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )
}


pick_col <- function(
    df,
    candidates,
    regex = NULL) {

  hit <- candidates[
    candidates %in%
    colnames(df)
  ]

  if (length(hit) > 0) {
    return(hit[1])
  }

  if (!is.null(regex)) {

    hit <- grep(
      regex,
      colnames(df),
      value = TRUE,
      ignore.case = TRUE
    )

    if (length(hit) > 0) {
      return(hit[1])
    }
  }

  stop(
    "Cannot identify required column.\nColumns are:\n",
    paste(
      colnames(df),
      collapse = ", "
    )
  )
}


get_metric <- function(
    df,
    regex) {

  if (!all(
    c(
      "Metric",
      "Value"
    ) %in%
    colnames(df)
  )) {
    stop(
      "Summary table must contain Metric and Value."
    )
  }

  idx <- grep(
    regex,
    df$Metric,
    ignore.case = TRUE
  )

  if (length(idx) == 0) {
    stop(
      "Cannot find metric matching: ",
      regex,
      "\nAvailable metrics:\n",
      paste(
        df$Metric,
        collapse = "\n"
      )
    )
  }

  as.numeric(
    df$Value[
      idx[1]
    ]
  )
}


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
# 2. Figure saving helpers
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


save_2x2 <- function(
    plots,
    stem,
    width = 12,
    height = 9) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          2,
          2,
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
      plots[[1]],
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    print(
      plots[[2]],
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )

    print(
      plots[[3]],
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
      )
    )

    print(
      plots[[4]],
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 2
      )
    )
  }

  grDevices::cairo_pdf(
    file.path(
      outdir,
      paste0(stem, ".pdf")
    ),
    width = width,
    height = height
  )

  draw()
  dev.off()

  open_tiff(
    file.path(
      outdir,
      paste0(
        stem,
        "_600dpi.tiff"
      )
    ),
    width,
    height
  )

  draw()
  dev.off()
}


save_figure4 <- function(
    pA,
    pB,
    pC,
    pD,
    stem,
    width = 13,
    height = 9) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          2,
          4,
          widths = unit(
            c(1, 1, 1, 1),
            "null"
          ),
          heights = unit(
            c(0.85, 1.35),
            "null"
          )
        )
      )
    )

    print(
      pA,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1:2
      )
    )

    print(
      pB,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 3:4
      )
    )

    print(
      pC,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1:3
      )
    )

    print(
      pD,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 4
      )
    )
  }

  grDevices::cairo_pdf(
    file.path(
      outdir,
      paste0(stem, ".pdf")
    ),
    width = width,
    height = height
  )

  draw()
  dev.off()

  open_tiff(
    file.path(
      outdir,
      paste0(
        stem,
        "_600dpi.tiff"
      )
    ),
    width,
    height
  )

  draw()
  dev.off()
}


save_figure5 <- function(
    pA,
    pB,
    pC,
    pD,
    stem,
    width = 13,
    height = 10) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          3,
          3,
          widths = unit(
            c(1.15, 1.15, 1),
            "null"
          ),
          heights = unit(
            c(1, 1, 1),
            "null"
          )
        )
      )
    )

    print(
      pA,
      vp = viewport(
        layout.pos.row = 1:3,
        layout.pos.col = 1:2
      )
    )

    print(
      pB,
      vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 3
      )
    )

    print(
      pC,
      vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 3
      )
    )

    print(
      pD,
      vp = viewport(
        layout.pos.row = 3,
        layout.pos.col = 3
      )
    )
  }

  grDevices::cairo_pdf(
    file.path(
      outdir,
      paste0(stem, ".pdf")
    ),
    width = width,
    height = height
  )

  draw()
  dev.off()

  open_tiff(
    file.path(
      outdir,
      paste0(
        stem,
        "_600dpi.tiff"
      )
    ),
    width,
    height
  )

  draw()
  dev.off()
}


# ============================================================
# FIGURE 3
# ============================================================

cat("Constructing Figure 3...\n")


# ------------------------------------------------------------
# Figure 3A: centroid PCA trajectory
# ------------------------------------------------------------

centroid_file <- paste0(
  "04_2026_temporal_structure_result/",
  "04_2026_PCA_group_centroid_coordinates.tsv"
)

cent <- read_tsv(
  centroid_file
)

group_col <- pick_col(
  cent,
  c(
    "Group",
    "group"
  )
)

pc1_col <- pick_col(
  cent,
  c("PC1")
)

pc2_col <- pick_col(
  cent,
  c("PC2")
)

cent$Group_fixed <- as.character(
  cent[[group_col]]
)

cent$Stage_label <- unname(
  labels11[
    cent$Group_fixed
  ]
)

cent$Order <- match(
  cent$Group_fixed,
  groups11
)

cent <- cent[
  order(cent$Order),
  ,
  drop = FALSE
]

p3a <- ggplot(
  cent,
  aes(
    x = .data[[pc1_col]],
    y = .data[[pc2_col]]
  )
) +
  geom_path(
    aes(group = 1),
    linewidth = 0.7,
    colour = "grey45"
  ) +
  geom_point(
    size = 3,
    colour = "#3366A6"
  ) +
  geom_text(
    aes(
      label = Stage_label
    ),
    vjust = -0.8,
    size = 3
  ) +
  labs(
    tag = "A",
    x = "PC1 (28.28%)",
    y = "PC2 (12.85%)",
    title = "Temporal centroid trajectory"
  ) +
  theme_pub


# ------------------------------------------------------------
# Figure 3B: PERMANOVA/PERMDISP
# ------------------------------------------------------------

summary4_file <- paste0(
  "04_2026_temporal_structure_result/",
  "04_2026_temporal_structure_summary.tsv"
)

sum4 <- read_tsv(
  summary4_file
)

r2 <- get_metric(
  sum4,
  "PERMANOVA.*R2|PERMANOVA_stage_R2"
)

fstat <- get_metric(
  sum4,
  "PERMANOVA.*F"
)

perm_p <- get_metric(
  sum4,
  "PERMANOVA.*P"
)

disp_p <- get_metric(
  sum4,
  "^PERMDISP_stage_P$|^PERMDISP_P$"
)

stat_text <- paste0(
  "PERMANOVA\n",
  "R² = ",
  sprintf("%.3f", r2),
  "\nF = ",
  sprintf("%.3f", fstat),
  "\nP < 0.001\n\n",
  "PERMDISP\n",
  "P = ",
  sprintf("%.3f", disp_p),
  "\n\n",
  "Stage-dependent separation\n",
  "without significant dispersion\n",
  "heterogeneity"
)

p3b <- ggplot() +
  annotate(
    "rect",
    xmin = 0,
    xmax = 1,
    ymin = 0,
    ymax = 1,
    fill = "grey97",
    colour = "grey70"
  ) +
  annotate(
    "text",
    x = 0.08,
    y = 0.90,
    label = stat_text,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1)
  ) +
  labs(
    tag = "B",
    title = "Global temporal structure"
  ) +
  theme_void() +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 16
    ),
    plot.title = element_text(
      face = "bold",
      size = 11
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ------------------------------------------------------------
# Figure 3C: adjacent-stage RMS distance
# ------------------------------------------------------------

distance_file <- paste0(
  "04_2026_temporal_structure_result/",
  "04_2026_adjacent_stage_centroid_distances.tsv"
)

adj <- read_tsv(
  distance_file
)

adj$Transition <- gsub(
  " -> ",
  " → ",
  adj$Transition,
  fixed = TRUE
)

adj$Transition <- factor(
  adj$Transition,
  levels = adj$Transition
)

p3c <- ggplot(
  adj,
  aes(
    x = Transition,
    y = Full_VST_RMS_distance
  )
) +
  geom_col(
    width = 0.70,
    fill = "#4C78A8"
  ) +
  geom_point(
    size = 1.8
  ) +
  labs(
    tag = "C",
    x = NULL,
    y = "Full-VST RMS centroid distance",
    title = "Adjacent-stage transcriptomic displacement"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# ------------------------------------------------------------
# Figure 3D: distance of samples to group centroid
# ------------------------------------------------------------

sample_dist_file <- paste0(
  "04_2026_temporal_structure_result/",
  "04_2026_sample_distance_to_fullVST_centroid.tsv"
)

sdist <- read_tsv(
  sample_dist_file
)

stage_col <- pick_col(
  sdist,
  c(
    "Stage",
    "stage",
    "Group",
    "group"
  )
)

distance_col <- pick_col(
  sdist,
  c(
    "Distance_to_group_centroid",
    "Distance_to_centroid",
    "Full_VST_distance_to_group_centroid",
    "Euclidean_distance_to_group_centroid"
  ),
  regex = "centroid.*distance|distance.*centroid"
)

stage_raw <- as.character(
  sdist[[stage_col]]
)

sdist$Stage_fixed <- ifelse(
  stage_raw %in% groups11,
  unname(labels11[stage_raw]),
  stage_raw
)

sdist$Stage_fixed <- factor(
  sdist$Stage_fixed,
  levels = unname(
    labels11[groups11]
  )
)

p3d <- ggplot(
  sdist,
  aes(
    x = Stage_fixed,
    y = .data[[distance_col]]
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.62,
    fill = "grey90"
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.3,
    alpha = 0.70
  ) +
  labs(
    tag = "D",
    x = "Sampling stage",
    y = "Distance to stage centroid",
    title = "Within-stage expression dispersion"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


save_2x2(
  list(
    p3a,
    p3b,
    p3c,
    p3d
  ),
  "Figure_3_global_temporal_transcriptomic_structure",
  width = 12,
  height = 9
)


# ============================================================
# FIGURE 4
# ============================================================

cat("Constructing Figure 4...\n")


# ------------------------------------------------------------
# Figure 4A: DEG numbers
# ------------------------------------------------------------

deg_file <- paste0(
  "05_2026_timepoint_vs_baseline_DESeq2_result/",
  "05_2026_DEG_summary.tsv"
)

deg <- read_tsv(
  deg_file
)

deg_long <- rbind(
  data.frame(
    Stage = deg$Stage,
    Direction = "Up",
    Number = deg$Up
  ),
  data.frame(
    Stage = deg$Stage,
    Direction = "Down",
    Number = -deg$Down
  )
)

deg_long$Stage <- factor(
  deg_long$Stage,
  levels = unname(
    labels26
  )
)

p4a <- ggplot(
  deg_long,
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
    linewidth = 0.35
  ) +
  scale_fill_manual(
    values = c(
      Up = "#D55E00",
      Down = "#0072B2"
    )
  ) +
  scale_y_continuous(
    labels = function(x) abs(x)
  ) +
  labs(
    tag = "A",
    x = "Sampling stage",
    y = "Number of DEGs",
    fill = NULL,
    title = "Timepoint-versus-Baseline differential expression"
  ) +
  theme_pub +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "top"
  )


# ------------------------------------------------------------
# Figure 4B: dynamic-gene hierarchy
# ------------------------------------------------------------

sum6_file <- paste0(
  "06_2026_continuous_time_spline_LRT_result/",
  "06_2026_continuous_time_model_summary.tsv"
)

sum7_file <- paste0(
  "07_2026_TPM_tau_SPM_peak_stage_result/",
  "07_2026_Tau_SPM_summary.tsv"
)

sum6 <- read_tsv(
  sum6_file
)

sum7 <- read_tsv(
  sum7_file
)

n_filtered <- get_metric(
  sum7,
  "^Filtered_genes$"
)

n_eligible <- get_metric(
  sum6,
  "Genes_retained_for_time_model"
)

n_lrt <- get_metric(
  sum6,
  "Dynamic_genes_FDR_LT_0.05"
)

n_high <- get_metric(
  sum7,
  "^High_confidence_dynamic$"
)

n_specific <- get_metric(
  sum7,
  "High_confidence_and_injury_specific"
)

hier <- data.frame(
  Category = c(
    "Expressed genes",
    "Injury-time eligible",
    "LRT dynamic",
    "Amplitude-filtered dynamic",
    "Injury-stage specific"
  ),
  Genes = c(
    n_filtered,
    n_eligible,
    n_lrt,
    n_high,
    n_specific
  ),
  stringsAsFactors = FALSE
)

hier$Percent <- 100 *
  hier$Genes /
  n_filtered

hier$Label <- paste0(
  format(
    hier$Genes,
    big.mark = ",",
    scientific = FALSE
  ),
  "\n(",
  sprintf(
    "%.1f%%",
    hier$Percent
  ),
  ")"
)

hier$Category <- factor(
  hier$Category,
  levels = rev(
    hier$Category
  )
)

p4b <- ggplot(
  hier,
  aes(
    x = Category,
    y = Genes
  )
) +
  geom_col(
    width = 0.66,
    fill = "#6A51A3"
  ) +
  geom_text(
    aes(
      label = Label
    ),
    hjust = -0.08,
    size = 3
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.24
      )
    ),
    labels = scales::comma
  ) +
  labs(
    tag = "B",
    x = NULL,
    y = "Number of genes",
    title = "Temporal dynamic-gene hierarchy"
  ) +
  theme_pub


# ------------------------------------------------------------
# Figure 4C/D: load publication objects from Step 11A
# ------------------------------------------------------------

p4c_file <- file.path(
  panel_dir,
  "Figure_4C_plot.rds"
)

p4d_file <- file.path(
  panel_dir,
  "Figure_4D_plot.rds"
)

if (!file.exists(p4c_file) ||
    !file.exists(p4d_file)) {

  stop(
    "Figure 4C/4D RDS files missing. ",
    "Run 11A_make_Figure4C_4D.R first."
  )
}

p4c <- readRDS(
  p4c_file
)

p4d <- readRDS(
  p4d_file
)


save_figure4(
  p4a,
  p4b,
  p4c,
  p4d,
  "Figure_4_temporal_expression_dynamics_resource",
  width = 13,
  height = 9
)


# ============================================================
# FIGURE 5
# ============================================================

cat("Constructing Figure 5...\n")


# ------------------------------------------------------------
# Figure 5A: genome-wide Spearman heatmap
# ------------------------------------------------------------

spearman_file <- paste0(
  "09_cross_cohort_temporal_response_concordance_result/",
  "09_cross_cohort_Spearman_matrix.tsv"
)

sp <- read_tsv(
  spearman_file
)

rownames(sp) <- sp[[1]]

sp_mat <- as.matrix(
  sp[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(sp_mat) <- "numeric"

sp_mat <- sp_mat[
  groups26,
  groups21,
  drop = FALSE
]

sp_long <- as.data.frame(
  as.table(sp_mat),
  stringsAsFactors = FALSE
)

colnames(sp_long) <- c(
  "Group2026",
  "Group2021",
  "Spearman"
)

sp_long$Stage2026 <- factor(
  unname(
    labels26[
      sp_long$Group2026
    ]
  ),
  levels = rev(
    unname(labels26)
  )
)

sp_long$Stage2021 <- factor(
  unname(
    labels21[
      sp_long$Group2021
    ]
  ),
  levels = unname(
    labels21
  )
)

p5a <- ggplot(
  sp_long,
  aes(
    x = Stage2021,
    y = Stage2026,
    fill = Spearman
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman
      )
    ),
    size = 3
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -0.2,
      0.8
    ),
    oob = scales::squish,
    name = "Spearman\nρ"
  ) +
  labs(
    tag = "A",
    x = "2021 historical-reference cohort",
    y = "2026 primary cohort",
    title = "Genome-wide Baseline-relative response concordance"
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
      hjust = 1
    ),
    axis.title = element_text(
      size = 10
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ------------------------------------------------------------
# Figure 5B/C: nearest-time pairs
# ------------------------------------------------------------

nearest_file <- paste0(
  "10_2026_sample_exclusion_sensitivity_result/",
  "10_nearest_time_cross_cohort_pairs.tsv"
)

near <- read_tsv(
  nearest_file
)

near$Pair <- paste0(
  near$Cohort2026_stage,
  " ↔ ",
  near$Cohort2021_stage
)

near$Pair <- factor(
  near$Pair,
  levels = rev(
    unique(
      near$Pair
    )
  )
)

p5b <- ggplot(
  near,
  aes(
    x = Spearman_rho,
    y = Pair
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Spearman_rho,
      yend = Pair
    ),
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.5,
    colour = "#3366A6"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman_rho
      )
    ),
    hjust = -0.35,
    size = 2.6
  ) +
  scale_x_continuous(
    limits = c(
      0,
      0.72
    ),
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    tag = "B",
    x = "Genome-wide Spearman ρ",
    y = NULL,
    title = "Chronologically nearest-stage concordance"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(
      size = 7
    )
  )


p5c <- ggplot(
  near,
  aes(
    x = Shared_DEG_direction_concordance,
    y = Pair
  )
) +
  geom_segment(
    aes(
      x = 0.80,
      xend =
        Shared_DEG_direction_concordance,
      yend = Pair
    ),
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    size = 2.5,
    colour = "#009E73"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.1f%%",
        100 *
          Shared_DEG_direction_concordance
      )
    ),
    hjust = -0.25,
    size = 2.5
  ) +
  scale_x_continuous(
    limits = c(
      0.80,
      1.015
    ),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    tag = "C",
    x = "Direction concordance",
    y = NULL,
    title = "Shared-DEG directional concordance"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(
      size = 7
    )
  )


# ------------------------------------------------------------
# Figure 5D: leave-one-out robustness
# ------------------------------------------------------------

loo_file <- paste0(
  "10B_single_sample_leave_one_out_result/",
  "10B_single_sample_influence_summary.tsv"
)

loo <- read_tsv(
  loo_file
)

loo$Removed_sample <- factor(
  loo$Removed_sample,
  levels = rev(
    loo$Removed_sample
  )
)

p5d <- ggplot(
  loo,
  aes(
    x =
      Cross_cohort_70cell_matrix_Spearman,
    y = Removed_sample
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.5,
    colour = "grey45"
  ) +
  geom_segment(
    aes(
      x = 0.95,
      xend =
        Cross_cohort_70cell_matrix_Spearman,
      yend = Removed_sample
    ),
    linewidth = 0.7,
    colour = "grey70"
  ) +
  geom_point(
    size = 2.8,
    colour = "#CC79A7"
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        Cross_cohort_70cell_matrix_Spearman
      )
    ),
    hjust = -0.28,
    size = 2.6
  ) +
  scale_x_continuous(
    limits = c(
      0.95,
      1.008
    )
  ) +
  labs(
    tag = "D",
    x = "Correlation with primary\n70-cell concordance matrix",
    y = "Sample removed",
    title = "Single-sample leave-one-out robustness"
  ) +
  theme_pub +
  theme(
    axis.text.y = element_text(
      size = 8
    )
  )


save_figure5(
  p5a,
  p5b,
  p5c,
  p5d,
  "Figure_5_cross_cohort_reproducibility_and_robustness",
  width = 13,
  height = 10
)


# ============================================================
# Save all component RDS objects
# ============================================================

plots <- list(
  Figure3A = p3a,
  Figure3B = p3b,
  Figure3C = p3c,
  Figure3D = p3d,
  Figure4A = p4a,
  Figure4B = p4b,
  Figure4C = p4c,
  Figure4D = p4d,
  Figure5A = p5a,
  Figure5B = p5b,
  Figure5C = p5c,
  Figure5D = p5d
)

saveRDS(
  plots,
  file.path(
    outdir,
    "Figures_3_4_5_all_plot_objects.rds"
  )
)


# ============================================================
# Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FINAL FIGURES 3, 4 AND 5 COMPLETED\n")
cat("============================================================\n")

cat(
  "Figure 3: global temporal transcriptomic structure\n"
)

cat(
  "Figure 4: temporal expression dynamics resource\n"
)

cat(
  "Figure 5: cross-cohort reproducibility and robustness\n"
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
