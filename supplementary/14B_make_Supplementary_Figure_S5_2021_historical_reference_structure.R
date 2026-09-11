#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(vegan)
})

cat("\n")
cat("============================================================\n")
cat("14B SUPPLEMENTARY FIGURE S5\n")
cat("2021 historical-reference cohort expression structure\n")
cat("Independently reprocessed previously published cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed paths
# ============================================================

vst_file <- paste0(
  "02_count_filtering_and_VST_result/",
  "2021_historical_reference/",
  "02_2021_VST_matrix.tsv"
)

metadata_candidates <- c(
  paste0(
    "02_count_filtering_and_VST_result/",
    "2021_historical_reference/",
    "02_2021_metadata_aligned.tsv"
  ),
  "rattus_meta_2021.tsv"
)

metadata_existing <- metadata_candidates[
  file.exists(metadata_candidates)
]

if (length(metadata_existing) == 0) {
  stop(
    "No 2021 metadata file found."
  )
}

metadata_file <- metadata_existing[1]


outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S5_2021_historical_reference_structure"
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

RANDOM_SEED <- 20260910


# ============================================================
# 1. Frozen 2021 stage definitions
# ============================================================

groups <- c(
  "Baseline",
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)

stage_labels <- c(
  Baseline = "Baseline",
  R4h  = "4 h",
  R8h  = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)

expected_stage_n <- c(
  Baseline = 6,
  R4h = 7,
  R8h = 8,
  R12h = 8,
  R16h = 8,
  R20h = 7,
  R24h = 8,
  R48h = 8
)

stage_colors <- c(
  Baseline = "#7A7A7A",
  R4h  = "#4C78A8",
  R8h  = "#5F8DB8",
  R12h = "#72A1B3",
  R16h = "#7FB28F",
  R20h = "#A3B76B",
  R24h = "#D2A455",
  R48h = "#C85A5A"
)


# ============================================================
# 2. Validate inputs
# ============================================================

if (!file.exists(vst_file)) {
  stop(
    "VST file not found:\n",
    vst_file
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found:\n",
    metadata_file
  )
}

cat(
  "VST file      : ",
  vst_file,
  "\n",
  sep = ""
)

cat(
  "Metadata file : ",
  metadata_file,
  "\n\n",
  sep = ""
)


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


# ------------------------------------------------------------
# Detect sample-ID column
# ------------------------------------------------------------

sample_col_candidates <- c(
  "sampleID",
  "SampleID",
  "sample",
  "Sample"
)

sample_col <- sample_col_candidates[
  sample_col_candidates %in%
    colnames(meta)
]

if (length(sample_col) == 0) {
  stop(
    "Could not identify sample-ID column in metadata."
  )
}

sample_col <- sample_col[1]

if (sample_col != "sampleID") {
  colnames(meta)[
    colnames(meta) == sample_col
  ] <- "sampleID"
}


# ------------------------------------------------------------
# Detect / derive group column
# ------------------------------------------------------------

if (!"group" %in%
    colnames(meta)) {

  if ("Group" %in%
      colnames(meta)) {

    colnames(meta)[
      colnames(meta) == "Group"
    ] <- "group"

  } else {

    stop(
      "Metadata lacks a recognizable group column."
    )
  }
}


if (nrow(meta) != 60) {
  stop(
    "Expected 60 historical-reference samples; observed ",
    nrow(meta)
  )
}

if (anyDuplicated(
  meta$sampleID
)) {
  stop(
    "Duplicated sampleID in 2021 metadata."
  )
}


unknown_groups <- setdiff(
  unique(
    meta$group
  ),
  groups
)

if (length(unknown_groups) > 0) {
  stop(
    "Unexpected 2021 group(s): ",
    paste(
      unknown_groups,
      collapse = ", "
    )
  )
}


meta$group <- factor(
  meta$group,
  levels = groups
)


# ============================================================
# 4. Verify frozen sample numbers
# ============================================================

observed_stage_n <- table(
  meta$group
)

cat("2021 stage sample numbers:\n")

print(
  observed_stage_n
)

cat("\n")


for (g in groups) {

  observed <- as.integer(
    observed_stage_n[g]
  )

  expected <- as.integer(
    expected_stage_n[g]
  )

  if (observed != expected) {
    stop(
      "Unexpected sample number for ",
      g,
      ": expected ",
      expected,
      ", observed ",
      observed
    )
  }
}

cat(
  "Stage sample-number verification: PASS\n\n"
)


# ============================================================
# 5. Read frozen 2021 VST matrix
# ============================================================

vst_df <- read.delim(
  vst_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

colnames(vst_df)[1] <- "GeneID"

if (anyDuplicated(
  vst_df$GeneID
)) {
  stop(
    "Duplicated GeneID in VST matrix."
  )
}


vst_mat <- as.matrix(
  vst_df[
    ,
    -1,
    drop = FALSE
  ]
)

storage.mode(
  vst_mat
) <- "numeric"

rownames(vst_mat) <-
  vst_df$GeneID


if (nrow(vst_mat) != 17328) {

  warning(
    "Expected 17,328 filtered genes; observed ",
    nrow(vst_mat)
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


# Align metadata to VST columns
meta <- meta[
  match(
    colnames(vst_mat),
    meta$sampleID
  ),
  ,
  drop = FALSE
]

if (!identical(
  meta$sampleID,
  colnames(vst_mat)
)) {
  stop(
    "Metadata/VST alignment failed."
  )
}


if (anyNA(vst_mat) ||
    any(!is.finite(vst_mat))) {
  stop(
    "Invalid numeric values in VST matrix."
  )
}


cat(
  "VST matrix : ",
  nrow(vst_mat),
  " genes x ",
  ncol(vst_mat),
  " samples\n\n",
  sep = ""
)


# ============================================================
# 6. PCA
# ============================================================

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


pca_coords <- data.frame(
  sampleID =
    rownames(
      pca$x
    ),

  PC1 =
    pca$x[, 1],

  PC2 =
    pca$x[, 2],

  PC3 =
    pca$x[, 3],

  group =
    meta$group,

  stringsAsFactors = FALSE
)


cat(
  "PCA variance:\n",
  "PC1 = ",
  sprintf(
    "%.2f%%",
    pca_var[1]
  ),
  "\n",
  "PC2 = ",
  sprintf(
    "%.2f%%",
    pca_var[2]
  ),
  "\n\n",
  sep = ""
)


# Frozen Step08 check
if (abs(
  pca_var[1] -
    29.52
) > 0.10) {

  warning(
    "PC1 differs from frozen Step08 value 29.52%."
  )
}

if (abs(
  pca_var[2] -
    23.97
) > 0.10) {

  warning(
    "PC2 differs from frozen Step08 value 23.97%."
  )
}


# ============================================================
# 7. PCA stage centroids
# ============================================================

centroid_list <- vector(
  "list",
  length(groups)
)

for (i in seq_along(groups)) {

  g <- groups[i]

  tmp <- pca_coords[
    pca_coords$group == g,
    ,
    drop = FALSE
  ]

  centroid_list[i] <- list(
    data.frame(
      group = g,
      PC1 = mean(
        tmp$PC1
      ),
      PC2 = mean(
        tmp$PC2
      ),
      stringsAsFactors = FALSE
    )
  )
}

pca_centroids <- do.call(
  rbind,
  centroid_list
)

pca_centroids$group <- factor(
  pca_centroids$group,
  levels = groups
)


# ============================================================
# 8. Sample Spearman correlation
# ============================================================

cor_mat <- cor(
  vst_mat,
  method = "spearman"
)

if (nrow(cor_mat) != 60 ||
    ncol(cor_mat) != 60) {
  stop(
    "Unexpected sample-correlation matrix dimension."
  )
}


cor_values <- cor_mat[
  upper.tri(
    cor_mat
  )
]

cor_summary <- data.frame(
  N_samples = 60,
  N_pairs =
    length(cor_values),

  Minimum_Spearman =
    min(cor_values),

  Median_Spearman =
    median(cor_values),

  Mean_Spearman =
    mean(cor_values),

  Maximum_Spearman =
    max(cor_values),

  stringsAsFactors = FALSE
)


# ============================================================
# 9. Order samples by stage
# ============================================================

sample_order <- unlist(
  lapply(
    groups,
    function(g) {

      meta$sampleID[
        meta$group == g
      ]
    }
  ),
  use.names = FALSE
)


if (length(sample_order) != 60) {
  stop(
    "Unexpected ordered sample count."
  )
}


cor_ordered <- cor_mat[
  sample_order,
  sample_order,
  drop = FALSE
]


# ============================================================
# 10. Long-format correlation table
# ============================================================

cor_long <- as.data.frame(
  as.table(
    cor_ordered
  ),
  stringsAsFactors = FALSE
)

colnames(cor_long) <- c(
  "Sample_Y",
  "Sample_X",
  "Spearman"
)

cor_long$X <- match(
  cor_long$Sample_X,
  sample_order
)

cor_long$Y <- match(
  cor_long$Sample_Y,
  sample_order
)


# ============================================================
# 11. Stage boundaries / centers for heatmap
# ============================================================

stage_sizes <- as.integer(
  table(
    factor(
      meta$group,
      levels = groups
    )
  )
)

stage_end <- cumsum(
  stage_sizes
)

stage_start <- c(
  1,
  head(
    stage_end,
    -1
  ) + 1
)

stage_center <- (
  stage_start +
    stage_end
) / 2

stage_boundaries <- stage_end[
  -length(
    stage_end
  )
] + 0.5


# ============================================================
# 12. Full-VST Euclidean distance
# ============================================================

sample_distance <- dist(
  t(vst_mat),
  method = "euclidean"
)


# ============================================================
# 13. PERMANOVA
# ============================================================

set.seed(
  RANDOM_SEED
)

perm <- vegan::adonis2(
  sample_distance ~ group,
  data = meta,
  permutations = 9999,
  by = "margin"
)


permanova_summary <- data.frame(
  R2 =
    perm[
      "group",
      "R2"
    ],

  F =
    perm[
      "group",
      "F"
    ],

  P =
    perm[
      "group",
      "Pr(>F)"
    ],

  stringsAsFactors = FALSE
)


# ============================================================
# 14. PERMDISP
# ============================================================

bd <- vegan::betadisper(
  sample_distance,
  group = meta$group,
  type = "median",
  bias.adjust = TRUE
)

set.seed(
  RANDOM_SEED
)

bd_perm <- vegan::permutest(
  bd,
  permutations = 9999
)


permdisp_summary <- data.frame(
  F =
    bd_perm$tab[
      "Groups",
      "F"
    ],

  P =
    bd_perm$tab[
      "Groups",
      "Pr(>F)"
    ],

  stringsAsFactors = FALSE
)


cat(
  "PERMANOVA:\n",
  "R2 = ",
  sprintf(
    "%.6f",
    permanova_summary$R2
  ),
  "\n",
  "F  = ",
  sprintf(
    "%.6f",
    permanova_summary$F
  ),
  "\n",
  "P  = ",
  format(
    permanova_summary$P,
    scientific = TRUE
  ),
  "\n\n",
  sep = ""
)

cat(
  "PERMDISP:\n",
  "F = ",
  sprintf(
    "%.6f",
    permdisp_summary$F
  ),
  "\n",
  "P = ",
  sprintf(
    "%.4f",
    permdisp_summary$P
  ),
  "\n\n",
  sep = ""
)


# Frozen Step08 check
if (abs(
  permanova_summary$R2 -
    0.519648
) > 0.005) {

  warning(
    "PERMANOVA R2 differs from frozen Step08 value."
  )
}

if (abs(
  permdisp_summary$P -
    0.2373
) > 0.05) {

  warning(
    "PERMDISP P differs noticeably from frozen Step08 value."
  )
}


# ============================================================
# 15. Full-expression-space stage centroids
# ============================================================

full_centroids <- matrix(
  NA_real_,
  nrow = length(groups),
  ncol = nrow(vst_mat),
  dimnames = list(
    groups,
    rownames(vst_mat)
  )
)


for (g in groups) {

  samples_g <- meta$sampleID[
    meta$group == g
  ]

  full_centroids[
    g,
  ] <- rowMeans(
    vst_mat[
      ,
      samples_g,
      drop = FALSE
    ]
  )
}


# ============================================================
# 16. Adjacent-stage centroid RMS distances
# ============================================================

adjacent_list <- vector(
  "list",
  length(groups) - 1
)


for (i in seq_len(
  length(groups) - 1
)) {

  g1 <- groups[i]
  g2 <- groups[i + 1]

  delta <- full_centroids[
    g2,
  ] -
    full_centroids[
      g1,
    ]

  rms <- sqrt(
    mean(
      delta^2
    )
  )

  adjacent_list[i] <- list(
    data.frame(
      From = g1,
      To = g2,

      Transition = paste0(
        unname(
          stage_labels[g1]
        ),
        " \u2192 ",
        unname(
          stage_labels[g2]
        )
      ),

      Order = i,

      RMS_distance =
        rms,

      stringsAsFactors = FALSE
    )
  )
}

adjacent_rms <- do.call(
  rbind,
  adjacent_list
)

adjacent_rms$Transition <- factor(
  adjacent_rms$Transition,
  levels = adjacent_rms$Transition
)


# ============================================================
# 17. Sample distance to own stage centroid
# ============================================================

sample_centroid_list <- vector(
  "list",
  nrow(meta)
)


for (i in seq_len(
  nrow(meta)
)) {

  sample_id <- meta$sampleID[i]

  g <- as.character(
    meta$group[i]
  )

  sample_vector <- vst_mat[
    ,
    sample_id
  ]

  centroid_vector <- full_centroids[
    g,
  ]

  euclidean_distance <- sqrt(
    sum(
      (
        sample_vector -
          centroid_vector
      )^2
    )
  )

  sample_centroid_list[i] <- list(
    data.frame(
      sampleID =
        sample_id,

      group =
        g,

      Distance_to_stage_centroid =
        euclidean_distance,

      stringsAsFactors = FALSE
    )
  )
}

sample_centroid_distance <- do.call(
  rbind,
  sample_centroid_list
)

sample_centroid_distance$group <- factor(
  sample_centroid_distance$group,
  levels = groups
)


# ============================================================
# 18. Publication theme
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
      size = 8
    ),

    legend.title = element_text(
      size = 8.5
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
# 19. Panel A: PCA + centroid trajectory
# ============================================================

pA <- ggplot(
  pca_coords,
  aes(
    x = PC1,
    y = PC2,
    colour = group
  )
) +

  geom_path(
    data = pca_centroids,
    aes(
      x = PC1,
      y = PC2,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "grey35",
    linewidth = 0.75,
    arrow = arrow(
      length = unit(
        0.11,
        "cm"
      ),
      type = "closed"
    )
  ) +

  geom_point(
    size = 2.3,
    alpha = 0.78
  ) +

  geom_point(
    data = pca_centroids,
    aes(
      x = PC1,
      y = PC2,
      fill = group
    ),
    inherit.aes = FALSE,
    shape = 21,
    colour = "black",
    stroke = 0.7,
    size = 4.0
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
    tag = "A",

    title =
      "Global expression structure of the 2021 historical-reference cohort",

    subtitle =
      "Variance-stabilized expression profiles; arrows connect stage centroids",

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
    legend.position = "right"
  )


# ============================================================
# 20. Panel B: sample-correlation heatmap
# ============================================================

pB <- ggplot(
  cor_long,
  aes(
    x = X,
    y = Y,
    fill = Spearman
  )
) +

  geom_tile() +

  geom_vline(
    xintercept =
      stage_boundaries,
    colour = "white",
    linewidth = 0.45
  ) +

  geom_hline(
    yintercept =
      stage_boundaries,
    colour = "white",
    linewidth = 0.45
  ) +

  scale_fill_gradientn(
    colours = c(
      "#2C4C7C",
      "#5E81AC",
      "#D8E2ED",
      "#F4E7CF",
      "#C77C4A"
    ),
    limits = c(
      min(
        cor_long$Spearman
      ),
      1
    ),
    name = "Spearman\ncorrelation"
  ) +

  scale_x_continuous(
    breaks =
      stage_center,

    labels = unname(
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
      stage_center,

    labels = unname(
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
    tag = "B",

    title =
      "Pairwise sample-expression correlation",

    subtitle =
      "Spearman correlations across 17,328 filtered genes",

    x = NULL,
    y = NULL
  ) +

  theme_minimal(
    base_size = 9.5
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

    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.5
    ),

    axis.text.y = element_text(
      size = 7.5
    ),

    legend.title = element_text(
      size = 8
    ),

    legend.text = element_text(
      size = 7.5
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 21. Panel C: adjacent-stage centroid RMS
# ============================================================

pC <- ggplot(
  adjacent_rms,
  aes(
    x = Transition,
    y = RMS_distance,
    group = 1
  )
) +

  geom_line(
    colour = "grey35",
    linewidth = 0.8
  ) +

  geom_point(
    size = 2.8,
    shape = 21,
    fill = "white",
    stroke = 0.7
  ) +

  labs(
    tag = "C",

    title =
      "Adjacent-stage transcriptomic transition magnitude",

    subtitle =
      "RMS distance between consecutive full-VST stage centroids",

    x = NULL,

    y =
      "Centroid RMS distance"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.5
    )
  )


# ============================================================
# 22. Panel D: sample distance to stage centroid
# ============================================================

pD <- ggplot(
  sample_centroid_distance,
  aes(
    x = group,
    y = Distance_to_stage_centroid
  )
) +

  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = "white",
    colour = "grey35",
    linewidth = 0.45
  ) +

  geom_jitter(
    aes(
      colour = group
    ),
    width = 0.15,
    height = 0,
    size = 2.1,
    alpha = 0.82
  ) +

  scale_colour_manual(
    values = stage_colors,
    guide = "none"
  ) +

  scale_x_discrete(
    labels = unname(
      stage_labels[
        groups
      ]
    )
  ) +

  labs(
    tag = "D",

    title =
      "Within-stage expression dispersion",

    subtitle = paste0(
      "Stage PERMANOVA R² = ",
      sprintf(
        "%.3f",
        permanova_summary$R2
      ),
      ", P < 0.001; PERMDISP P = ",
      sprintf(
        "%.3f",
        permdisp_summary$P
      )
    ),

    x = NULL,

    y =
      "Euclidean distance to stage centroid"
  ) +

  theme_pub +

  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.5
    )
  )


# ============================================================
# 23. TIFF helper
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
# 24. Save individual panels
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
  "Supplementary_Figure_S5A_2021_PCA"
)

save_panel(
  pB,
  "Supplementary_Figure_S5B_2021_sample_correlation"
)

save_panel(
  pC,
  "Supplementary_Figure_S5C_2021_adjacent_stage_RMS"
)

save_panel(
  pD,
  "Supplementary_Figure_S5D_2021_within_stage_dispersion"
)


# ============================================================
# 25. Draw combined S5
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


# ============================================================
# 26. Save final publication figure
# ============================================================

pdf_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S5_",
    "2021_historical_reference_expression_structure.pdf"
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
    "Supplementary_Figure_S5_",
    "2021_historical_reference_expression_structure_600dpi.tiff"
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
# 27. Save source-data tables
# ============================================================

write.table(
  pca_coords,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5A_PCA_source_data.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  cor_ordered,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5B_sample_Spearman_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  adjacent_rms,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5C_adjacent_stage_RMS.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  sample_centroid_distance,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5D_sample_centroid_distance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  cor_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5_correlation_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  permanova_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5_PERMANOVA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  permdisp_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S5_PERMDISP.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 28. Publication caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S5 | Internal transcriptomic structure of the ",
  "independently reprocessed 2021 historical-reference cohort. ",
  "a, Principal component analysis of variance-stabilized expression ",
  "profiles across 60 samples and 17,328 genes retained after the ",
  "predefined count-based filtering procedure. Points represent ",
  "individual samples and arrows connect sampling-stage centroids in ",
  "chronological order. ",
  "b, Pairwise Spearman correlation matrix of variance-stabilized ",
  "expression profiles across all 60 samples, ordered by sampling stage. ",
  "White boundaries delineate sampling-stage groups. ",
  "c, Adjacent-stage transcriptomic transition magnitudes quantified as ",
  "root-mean-square distances between consecutive stage centroids in the ",
  "full variance-stabilized expression space. These distances are treated ",
  "as descriptive measures rather than biological velocities. ",
  "d, Euclidean distance of individual samples from their corresponding ",
  "sampling-stage centroid. Global transcriptomic profiles differed ",
  "significantly among sampling stages (PERMANOVA, R-squared = ",
  sprintf(
    "%.3f",
    permanova_summary$R2
  ),
  ", P < 0.001), whereas multivariate dispersion did not differ ",
  "significantly among stages (PERMDISP, P = ",
  sprintf(
    "%.3f",
    permdisp_summary$P
  ),
  "). The 2021 dataset was analysed independently as a previously ",
  "published historical-reference cohort and was not pooled or ",
  "batch-corrected with the 2026 primary cohort."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S5_caption.txt"
  )
)


# ============================================================
# 29. Save plot objects / session
# ============================================================

saveRDS(
  list(
    S5A = pA,
    S5B = pB,
    S5C = pC,
    S5D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S5_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "Supplementary_Figure_S5_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 30. Console summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S5 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Historical-reference samples : ",
  ncol(vst_mat),
  "\n",
  sep = ""
)

cat(
  "Filtered genes               : ",
  nrow(vst_mat),
  "\n",
  sep = ""
)

cat(
  "PC1 variance                 : ",
  sprintf(
    "%.2f%%",
    pca_var[1]
  ),
  "\n",
  sep = ""
)

cat(
  "PC2 variance                 : ",
  sprintf(
    "%.2f%%",
    pca_var[2]
  ),
  "\n",
  sep = ""
)

cat(
  "Minimum sample Spearman      : ",
  sprintf(
    "%.4f",
    cor_summary$Minimum_Spearman
  ),
  "\n",
  sep = ""
)

cat(
  "Median sample Spearman       : ",
  sprintf(
    "%.4f",
    cor_summary$Median_Spearman
  ),
  "\n",
  sep = ""
)

cat(
  "PERMANOVA R2                 : ",
  sprintf(
    "%.6f",
    permanova_summary$R2
  ),
  "\n",
  sep = ""
)

cat(
  "PERMANOVA P                  : ",
  format(
    permanova_summary$P,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "PERMDISP P                   : ",
  sprintf(
    "%.4f",
    permdisp_summary$P
  ),
  "\n\n",
  sep = ""
)

cat("Adjacent-stage centroid RMS:\n")

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
  "\n\n",
  sep = ""
)

cat("PASS\n")
