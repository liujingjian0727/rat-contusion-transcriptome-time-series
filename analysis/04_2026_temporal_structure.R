#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  bitmapType = "cairo"
)

suppressPackageStartupMessages({

  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("R package 'vegan' is required.")
  }

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("R package 'ggplot2' is required.")
  }

  library(vegan)
  library(ggplot2)
})

cat("\n")
cat("============================================================\n")
cat("04 2026 TEMPORAL TRANSCRIPTOMIC STRUCTURE\n")
cat("PERMANOVA, PERMDISP, centroid trajectory,\n")
cat("and adjacent-stage transcriptomic distances\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")


# ============================================================
# 0. Fixed parameters
# ============================================================

EXPECTED_SAMPLES <- 77
EXPECTED_GENES <- 18364

N_PERMUTATIONS <- 9999
RANDOM_SEED <- 20260910

set.seed(RANDOM_SEED)

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

group_time <- c(
  Baseline = NA,
  R0h = 0,
  R1h = 1,
  R2h = 2,
  R6h = 6,
  R10h = 10,
  R14h = 14,
  R18h = 18,
  R36h = 36,
  R60h = 60,
  R72h = 72
)


# ============================================================
# 1. Input/output files
# ============================================================

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
  "03_2026_samples_requiring_expression_review.tsv"
)

outdir <- "04_2026_temporal_structure_result"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. Plot saving helper
# ============================================================

save_plot_both <- function(
    plot_object,
    filename_stem,
    width,
    height) {

  ggsave(
    filename = file.path(
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
# 3. Check files
# ============================================================

required_files <- c(
  vst_file,
  meta_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("Input files found: PASS\n\n")


# ============================================================
# 4. Read VST matrix
# ============================================================

cat("Reading 2026 VST matrix...\n")

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
  stop("Duplicated GeneIDs detected.")
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
# 5. Read metadata
# ============================================================

cat("Reading 2026 metadata...\n")

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
# 6. Explicit sample alignment
# ============================================================

sample_ids <- colnames(vst_mat)

if (!setequal(sample_ids, meta$sampleID)) {
  stop("VST and metadata sample sets differ.")
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
  stop("Metadata alignment failed.")
}

rownames(meta) <- meta$sampleID

meta$group <- factor(
  meta$group,
  levels = group_levels
)

if (anyNA(meta$group)) {
  stop("Unexpected group detected.")
}

meta$stage <- factor(
  group_labels[as.character(meta$group)],
  levels = unname(group_labels)
)


# ============================================================
# 7. Dimension validation
# ============================================================

cat("\nInput dimensions:\n")
cat(
  "VST matrix :",
  nrow(vst_mat),
  "genes x",
  ncol(vst_mat),
  "samples\n"
)
cat(
  "Metadata   :",
  nrow(meta),
  "samples\n\n"
)

if (ncol(vst_mat) != EXPECTED_SAMPLES) {
  stop(
    "Expected ",
    EXPECTED_SAMPLES,
    " samples, observed ",
    ncol(vst_mat)
  )
}

if (nrow(vst_mat) != EXPECTED_GENES) {
  stop(
    "Expected ",
    EXPECTED_GENES,
    " genes, observed ",
    nrow(vst_mat)
  )
}

cat("Dimension validation: PASS\n")
cat("No sample exclusion applied: PASS\n\n")


# ============================================================
# 8. Report Step 03 review samples without excluding them
# ============================================================

if (file.exists(review_file)) {

  review <- read.delim(
    review_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  cat(
    "Expression-level review samples from Step 03 :",
    nrow(review),
    "\n"
  )

  cat(
    "These samples remain INCLUDED in Step 04.\n\n"
  )
}


# ============================================================
# 9. PCA
# ============================================================

cat("Running PCA for centroid visualization...\n")

pca <- prcomp(
  t(vst_mat),
  center = TRUE,
  scale. = FALSE
)

variance <- pca$sdev^2
variance_percent <- 100 * variance / sum(variance)

pca_scores <- data.frame(
  sampleID = rownames(pca$x),
  pca$x,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pca_scores$group <- as.character(
  meta[pca_scores$sampleID, "group"]
)

pca_scores$stage <- factor(
  group_labels[pca_scores$group],
  levels = unname(group_labels)
)


# ============================================================
# 10. Full VST Euclidean distance
# ============================================================

cat("Calculating full-expression Euclidean distance matrix...\n")

sample_distance <- dist(
  t(vst_mat),
  method = "euclidean"
)

if (length(sample_distance) !=
    EXPECTED_SAMPLES * (EXPECTED_SAMPLES - 1) / 2) {

  stop("Unexpected sample-distance length.")
}

cat(
  "Distance-vector length:",
  length(sample_distance),
  ": PASS\n\n"
)


# ============================================================
# 11. PERMANOVA
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

permanova_table <- data.frame(
  Term = rownames(permanova),
  as.data.frame(permanova),
  row.names = NULL,
  check.names = FALSE
)

write.table(
  permanova_table,
  file = file.path(
    outdir,
    "04_2026_PERMANOVA_stage.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nPERMANOVA result:\n")
print(permanova)
cat("\n")


# ============================================================
# 12. PERMDISP
# ============================================================

cat(
  "Running PERMDISP with ",
  N_PERMUTATIONS,
  " permutations...\n",
  sep = ""
)

# Spatial median is robust to individual extreme samples.
# bias.adjust=TRUE is useful because group sizes are not identical.
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

cat("\nPERMDISP permutation result:\n")
print(bd_perm)
cat("\n")

bd_anova <- anova(bd)

permdisp_anova_table <- data.frame(
  Term = rownames(bd_anova),
  as.data.frame(bd_anova),
  row.names = NULL,
  check.names = FALSE
)

write.table(
  permdisp_anova_table,
  file = file.path(
    outdir,
    "04_2026_PERMDISP_ANOVA.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

bd_perm_table <- data.frame(
  Term = rownames(bd_perm$tab),
  as.data.frame(bd_perm$tab),
  row.names = NULL,
  check.names = FALSE
)

write.table(
  bd_perm_table,
  file = file.path(
    outdir,
    "04_2026_PERMDISP_permutation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 13. Save sample-level dispersion distances
# ============================================================

dispersion_samples <- data.frame(
  sampleID = names(bd$distances),
  Distance_to_group_spatial_median =
    as.numeric(bd$distances),
  stringsAsFactors = FALSE
)

dispersion_samples$group <- as.character(
  meta[
    dispersion_samples$sampleID,
    "group"
  ]
)

dispersion_samples$stage <- group_labels[
  dispersion_samples$group
]

write.table(
  dispersion_samples,
  file = file.path(
    outdir,
    "04_2026_sample_dispersion_distances.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 14. Group dispersion summary
# ============================================================

dispersion_group_summary <- do.call(
  rbind,
  lapply(
    group_levels,
    function(g) {

      x <- dispersion_samples[
        dispersion_samples$group == g,
        "Distance_to_group_spatial_median"
      ]

      data.frame(
        Group = g,
        Stage = unname(group_labels[g]),
        N = length(x),
        Mean_distance = mean(x),
        Median_distance = median(x),
        SD_distance = sd(x),
        Minimum_distance = min(x),
        Maximum_distance = max(x),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  dispersion_group_summary,
  file = file.path(
    outdir,
    "04_2026_group_dispersion_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 15. Dispersion plot
# ============================================================

dispersion_samples$stage <- factor(
  dispersion_samples$stage,
  levels = unname(group_labels)
)

p_disp <- ggplot(
  dispersion_samples,
  aes(
    x = stage,
    y = Distance_to_group_spatial_median
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.8,
    alpha = 0.75
  ) +
  labs(
    x = "Sampling stage",
    y = "Distance to group spatial median"
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
  p_disp,
  "04_2026_PERMDISP_group_dispersion",
  width = 9,
  height = 5.8
)


# ============================================================
# 16. Calculate full-expression group centroids
# ============================================================

cat("Calculating group centroids in full VST expression space...\n")

centroid_mat <- sapply(
  group_levels,
  function(g) {

    group_samples <- meta$sampleID[
      meta$group == g
    ]

    rowMeans(
      vst_mat[
        ,
        group_samples,
        drop = FALSE
      ]
    )
  }
)

colnames(centroid_mat) <- group_levels
rownames(centroid_mat) <- rownames(vst_mat)

centroid_output <- data.frame(
  GeneID = rownames(centroid_mat),
  centroid_mat,
  check.names = FALSE
)

write.table(
  centroid_output,
  file = file.path(
    outdir,
    "04_2026_fullVST_group_centroids.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 17. PCA-space group centroids
# ============================================================

pca_centroids <- do.call(
  rbind,
  lapply(
    group_levels,
    function(g) {

      x <- pca_scores[
        pca_scores$group == g,
        ,
        drop = FALSE
      ]

      data.frame(
        Group = g,
        Stage = unname(group_labels[g]),
        Time_h = group_time[g],
        N = nrow(x),
        PC1 = mean(x$PC1),
        PC2 = mean(x$PC2),
        stringsAsFactors = FALSE
      )
    }
  )
)

pca_centroids$Stage <- factor(
  pca_centroids$Stage,
  levels = unname(group_labels)
)

write.table(
  pca_centroids,
  file = file.path(
    outdir,
    "04_2026_PCA_group_centroid_coordinates.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 18. PCA temporal centroid trajectory plot
# ============================================================

stage_palette <- setNames(
  grDevices::hcl.colors(
    length(group_levels),
    palette = "Dynamic"
  ),
  unname(group_labels)
)

p_trajectory <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2
  )
) +
  geom_point(
    aes(colour = stage),
    size = 2,
    alpha = 0.40
  ) +
  geom_path(
    data = pca_centroids,
    aes(
      x = PC1,
      y = PC2,
      group = 1
    ),
    linewidth = 0.9,
    arrow = grid::arrow(
      length = grid::unit(
        0.15,
        "inches"
      ),
      type = "closed"
    )
  ) +
  geom_point(
    data = pca_centroids,
    aes(
      x = PC1,
      y = PC2,
      colour = Stage
    ),
    size = 4
  ) +
  geom_text(
    data = pca_centroids,
    aes(
      x = PC1,
      y = PC2,
      label = Stage
    ),
    size = 3.2,
    vjust = -1,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = stage_palette,
    drop = FALSE
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
  ) +
  theme(
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    )
  )

save_plot_both(
  p_trajectory,
  "04_2026_temporal_centroid_trajectory",
  width = 9,
  height = 7
)


# ============================================================
# 19. Adjacent-stage centroid distances
# ============================================================

cat("Calculating adjacent-stage centroid distances...\n")

transition_from <- group_levels[
  seq_len(length(group_levels) - 1)
]

transition_to <- group_levels[
  2:length(group_levels)
]

adjacent_distance <- do.call(
  rbind,
  lapply(
    seq_along(transition_from),
    function(i) {

      g1 <- transition_from[i]
      g2 <- transition_to[i]

      delta <- centroid_mat[, g2] -
        centroid_mat[, g1]

      euclidean_full <- sqrt(
        sum(delta^2)
      )

      rms_full <- sqrt(
        mean(delta^2)
      )

      pca2_delta <- as.numeric(
        pca_centroids[
          pca_centroids$Group == g2,
          c("PC1", "PC2")
        ]
      ) -
        as.numeric(
          pca_centroids[
            pca_centroids$Group == g1,
            c("PC1", "PC2")
          ]
        )

      pca2_distance <- sqrt(
        sum(pca2_delta^2)
      )

      if (g1 == "Baseline") {

        elapsed_hours <- NA_real_
        distance_rate <- NA_real_

      } else {

        elapsed_hours <-
          group_time[g2] -
          group_time[g1]

        distance_rate <-
          rms_full / elapsed_hours
      }

      data.frame(
        From_group = g1,
        To_group = g2,
        From_stage = group_labels[g1],
        To_stage = group_labels[g2],
        Transition = paste0(
          group_labels[g1],
          " -> ",
          group_labels[g2]
        ),
        Elapsed_hours = elapsed_hours,
        Full_VST_Euclidean_distance =
          euclidean_full,
        Full_VST_RMS_distance =
          rms_full,
        PCA2_centroid_distance =
          pca2_distance,
        RMS_distance_per_hour =
          distance_rate,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  adjacent_distance,
  file = file.path(
    outdir,
    "04_2026_adjacent_stage_centroid_distances.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Adjacent-stage distance plot
# ============================================================

adjacent_distance$Transition <- factor(
  adjacent_distance$Transition,
  levels = adjacent_distance$Transition
)

p_distance <- ggplot(
  adjacent_distance,
  aes(
    x = Transition,
    y = Full_VST_RMS_distance
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_point(
    size = 2.2
  ) +
  labs(
    x = "Adjacent sampling-stage transition",
    y = "Centroid RMS distance in VST expression space"
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
  p_distance,
  "04_2026_adjacent_stage_centroid_distance",
  width = 10,
  height = 6
)


# ============================================================
# 21. Distance-per-hour plot for injury stages only
# ============================================================

distance_rate_data <- adjacent_distance[
  !is.na(adjacent_distance$RMS_distance_per_hour),
  ,
  drop = FALSE
]

distance_rate_data$Transition <- factor(
  distance_rate_data$Transition,
  levels = distance_rate_data$Transition
)

p_rate <- ggplot(
  distance_rate_data,
  aes(
    x = Transition,
    y = RMS_distance_per_hour
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_point(
    size = 2.2
  ) +
  labs(
    x = "Injury-stage transition",
    y = "Centroid RMS distance per hour"
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
  p_rate,
  "04_2026_adjacent_stage_distance_per_hour",
  width = 9.5,
  height = 5.8
)


# ============================================================
# 22. Group-level distance to own centroid
# ============================================================

cat("Calculating sample distances to full-expression group centroids...\n")

sample_centroid_distance <- do.call(
  rbind,
  lapply(
    group_levels,
    function(g) {

      samples <- meta$sampleID[
        meta$group == g
      ]

      centroid <- centroid_mat[, g]

      dist_values <- sapply(
        samples,
        function(s) {

          sqrt(
            mean(
              (
                vst_mat[, s] -
                centroid
              )^2
            )
          )
        }
      )

      data.frame(
        sampleID = samples,
        Group = g,
        Stage = unname(group_labels[g]),
        RMS_distance_to_group_centroid =
          as.numeric(dist_values),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.table(
  sample_centroid_distance,
  file = file.path(
    outdir,
    "04_2026_sample_distance_to_fullVST_centroid.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Analysis summary
# ============================================================

# adonis2() with a single predictor reports the tested term
# as "Model", followed by "Residual" and "Total".
if (!"Model" %in% rownames(permanova)) {
  stop(
    "Unexpected PERMANOVA output row names: ",
    paste(rownames(permanova), collapse = ", ")
  )
}

permanova_r2 <- as.numeric(
  permanova["Model", "R2"]
)

permanova_f <- as.numeric(
  permanova["Model", "F"]
)

permanova_p <- as.numeric(
  permanova["Model", "Pr(>F)"]
)

if (!"Groups" %in% rownames(bd_perm$tab)) {
  stop(
    "Unexpected PERMDISP output row names: ",
    paste(rownames(bd_perm$tab), collapse = ", ")
  )
}

permdisp_p <- as.numeric(
  bd_perm$tab["Groups", "Pr(>F)"]
)

if (length(permanova_r2) != 1 ||
    length(permanova_f) != 1 ||
    length(permanova_p) != 1 ||
    length(permdisp_p) != 1 ||
    any(!is.finite(c(
      permanova_r2,
      permanova_f,
      permanova_p,
      permdisp_p
    )))) {

  stop(
    "Failed to extract PERMANOVA/PERMDISP summary statistics."
  )
}

summary_table <- data.frame(
  Metric = c(
    "Samples",
    "Genes",
    "PC1_variance_percent",
    "PC2_variance_percent",
    "PERMANOVA_permutations",
    "PERMANOVA_stage_R2",
    "PERMANOVA_stage_F",
    "PERMANOVA_stage_P",
    "PERMDISP_permutations",
    "PERMDISP_stage_P"
  ),
  Value = c(
    ncol(vst_mat),
    nrow(vst_mat),
    variance_percent[1],
    variance_percent[2],
    N_PERMUTATIONS,
    permanova_r2,
    permanova_f,
    permanova_p,
    N_PERMUTATIONS,
    permdisp_p
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_table,
  file = file.path(
    outdir,
    "04_2026_temporal_structure_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 24. Save objects
# ============================================================

saveRDS(
  permanova,
  file = file.path(
    outdir,
    "04_2026_PERMANOVA_object.rds"
  )
)

saveRDS(
  bd,
  file = file.path(
    outdir,
    "04_2026_PERMDISP_object.rds"
  )
)

saveRDS(
  pca,
  file = file.path(
    outdir,
    "04_2026_PCA_object.rds"
  )
)


# ============================================================
# 25. Session info
# ============================================================

sink(
  file.path(
    outdir,
    "04_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 26. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("04 2026 TEMPORAL STRUCTURE COMPLETED\n")
cat("============================================================\n")

cat(
  "Samples analyzed              : ",
  ncol(vst_mat),
  "\n",
  sep = ""
)

cat(
  "Genes analyzed                : ",
  nrow(vst_mat),
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
cat("All 77 samples were retained.\n")
cat("No sample exclusion was performed.\n")
cat("PERMANOVA used Euclidean distances of VST expression.\n")
cat("PERMDISP used distance to group spatial median.\n")
cat("Centroids were arithmetic means in VST expression space.\n")
cat("============================================================\n\n")

cat("Adjacent-stage centroid distances:\n")

print(
  adjacent_distance[
    ,
    c(
      "Transition",
      "Elapsed_hours",
      "Full_VST_RMS_distance",
      "RMS_distance_per_hour"
    )
  ],
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 04 is complete. ",
  "The next step is Step 05: ",
  "2026 timepoint-versus-Baseline differential-expression overview.\n"
)
