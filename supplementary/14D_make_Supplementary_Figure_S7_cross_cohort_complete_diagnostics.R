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
cat("14D SUPPLEMENTARY FIGURE S7\n")
cat("Complete cross-cohort response-concordance diagnostics\n")
cat("2026 primary cohort versus 2021 historical-reference cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Output paths
# ============================================================

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S7_cross_cohort_complete_diagnostics"
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
# 1. Frozen stage definitions
# ============================================================

stages_2026 <- c(
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

times_2026 <- c(
  0,
  1,
  2,
  6,
  10,
  14,
  18,
  36,
  60,
  72
)

names(times_2026) <- stages_2026


labels_2026 <- c(
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


stages_2021 <- c(
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)

times_2021 <- c(
  4,
  8,
  12,
  16,
  20,
  24,
  48
)

names(times_2021) <- stages_2021


labels_2021 <- c(
  R4h  = "4 h",
  R8h  = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)


# ============================================================
# 2. File-location helper
# ============================================================

locate_unique_file <- function(
    preferred_basenames,
    regex_pattern = NULL) {

  all_files <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE
  )

  basenames <- basename(
    all_files
  )


  for (target in preferred_basenames) {

    direct_hit <- which(
      basenames == target
    )

    if (length(direct_hit) == 1) {

      return(
        normalizePath(
          all_files[direct_hit],
          mustWork = TRUE
        )
      )
    }

    if (length(direct_hit) > 1) {

      cat(
        "\nMultiple files found for basename:\n",
        target,
        "\n",
        sep = ""
      )

      cat(
        paste(
          all_files[direct_hit],
          collapse = "\n"
        ),
        "\n"
      )

      stop(
        "Ambiguous input file."
      )
    }
  }


  if (!is.null(
    regex_pattern
  )) {

    regex_hit <- grep(
      regex_pattern,
      basenames,
      ignore.case = TRUE
    )

    if (length(regex_hit) == 1) {

      return(
        normalizePath(
          all_files[regex_hit],
          mustWork = TRUE
        )
      )
    }

    if (length(regex_hit) > 1) {

      cat(
        "\nMultiple regex-matched files found:\n"
      )

      cat(
        paste(
          all_files[regex_hit],
          collapse = "\n"
        ),
        "\n"
      )

      stop(
        "Ambiguous regex input."
      )
    }
  }


  stop(
    "Required input file not found."
  )
}


# ============================================================
# 3. Locate frozen genome-wide matrices
# ============================================================

lfc_2026_file <- locate_unique_file(
  preferred_basenames = c(
    "05_2026_genomewide_log2FC_matrix.tsv"
  ),
  regex_pattern =
    "^05_2026.*genomewide.*log2fc.*matrix.*\\.tsv$"
)


padj_2026_file <- locate_unique_file(
  preferred_basenames = c(
    "05_2026_genomewide_PADJ_matrix.tsv",
    "05_2026_genomewide_padj_matrix.tsv"
  ),
  regex_pattern =
    "^05_2026.*genomewide.*padj.*matrix.*\\.tsv$"
)


lfc_2021_file <- locate_unique_file(
  preferred_basenames = c(
    "08_2021_genomewide_log2FC_matrix.tsv"
  ),
  regex_pattern =
    "^08_2021.*genomewide.*log2fc.*matrix.*\\.tsv$"
)


padj_2021_file <- locate_unique_file(
  preferred_basenames = c(
    "08_2021_genomewide_PADJ_matrix.tsv",
    "08_2021_genomewide_padj_matrix.tsv"
  ),
  regex_pattern =
    "^08_2021.*genomewide.*padj.*matrix.*\\.tsv$"
)


cat(
  "2026 log2FC : ",
  lfc_2026_file,
  "\n",
  sep = ""
)

cat(
  "2026 PADJ   : ",
  padj_2026_file,
  "\n",
  sep = ""
)

cat(
  "2021 log2FC : ",
  lfc_2021_file,
  "\n",
  sep = ""
)

cat(
  "2021 PADJ   : ",
  padj_2021_file,
  "\n\n",
  sep = ""
)


# ============================================================
# 4. Matrix-reading helper
# ============================================================

read_gene_matrix <- function(
    file,
    expected_columns,
    matrix_name) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )


  if (ncol(x) < 2) {
    stop(
      matrix_name,
      " has fewer than two columns."
    )
  }


  colnames(x)[1] <- "GeneID"


  if (anyDuplicated(
    x$GeneID
  )) {
    stop(
      "Duplicated GeneID in ",
      matrix_name
    )
  }


  missing_columns <- setdiff(
    expected_columns,
    colnames(x)
  )


  if (length(
    missing_columns
  ) > 0) {

    cat(
      "\nAvailable columns in ",
      matrix_name,
      ":\n",
      sep = ""
    )

    cat(
      paste(
        colnames(x),
        collapse = "\n"
      ),
      "\n"
    )

    stop(
      "Missing expected column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  }


  mat <- as.matrix(
    x[
      ,
      expected_columns,
      drop = FALSE
    ]
  )

  storage.mode(
    mat
  ) <- "numeric"

  rownames(mat) <-
    x$GeneID


  list(
    dataframe = x,
    matrix = mat
  )
}


# ============================================================
# 5. Read frozen matrices
# ============================================================

lfc_2026_obj <- read_gene_matrix(
  lfc_2026_file,
  stages_2026,
  "2026 log2FC matrix"
)

padj_2026_obj <- read_gene_matrix(
  padj_2026_file,
  stages_2026,
  "2026 PADJ matrix"
)

lfc_2021_obj <- read_gene_matrix(
  lfc_2021_file,
  stages_2021,
  "2021 log2FC matrix"
)

padj_2021_obj <- read_gene_matrix(
  padj_2021_file,
  stages_2021,
  "2021 PADJ matrix"
)


# ============================================================
# 6. Within-cohort gene-set consistency
# ============================================================

if (!setequal(
  rownames(
    lfc_2026_obj$matrix
  ),
  rownames(
    padj_2026_obj$matrix
  )
)) {
  stop(
    "2026 log2FC and PADJ gene sets differ."
  )
}


if (!setequal(
  rownames(
    lfc_2021_obj$matrix
  ),
  rownames(
    padj_2021_obj$matrix
  )
)) {
  stop(
    "2021 log2FC and PADJ gene sets differ."
  )
}


# Align PADJ matrices to LFC gene order
padj_2026_mat <- padj_2026_obj$matrix[
  match(
    rownames(
      lfc_2026_obj$matrix
    ),
    rownames(
      padj_2026_obj$matrix
    )
  ),
  ,
  drop = FALSE
]


padj_2021_mat <- padj_2021_obj$matrix[
  match(
    rownames(
      lfc_2021_obj$matrix
    ),
    rownames(
      padj_2021_obj$matrix
    )
  ),
  ,
  drop = FALSE
]


lfc_2026_mat <- lfc_2026_obj$matrix
lfc_2021_mat <- lfc_2021_obj$matrix


# ============================================================
# 7. Common independently filtered gene universe
# ============================================================

common_genes <- intersect(
  rownames(
    lfc_2026_mat
  ),
  rownames(
    lfc_2021_mat
  )
)


cat(
  "2026 genes         : ",
  nrow(lfc_2026_mat),
  "\n",
  sep = ""
)

cat(
  "2021 genes         : ",
  nrow(lfc_2021_mat),
  "\n",
  sep = ""
)

cat(
  "Common genes       : ",
  length(common_genes),
  "\n\n",
  sep = ""
)


if (length(
  common_genes
) != 16918) {

  warning(
    "Expected 16,918 common genes from frozen Step09; observed ",
    length(common_genes)
  )
}


# Align all matrices to exactly the same common-gene order
lfc_2026_common <- lfc_2026_mat[
  common_genes,
  ,
  drop = FALSE
]

padj_2026_common <- padj_2026_mat[
  common_genes,
  ,
  drop = FALSE
]

lfc_2021_common <- lfc_2021_mat[
  common_genes,
  ,
  drop = FALSE
]

padj_2021_common <- padj_2021_mat[
  common_genes,
  ,
  drop = FALSE
]


# ============================================================
# 8. Pairwise 10 x 7 cross-cohort diagnostics
# ============================================================

n_2026 <- length(
  stages_2026
)

n_2021 <- length(
  stages_2021
)


spearman_mat <- matrix(
  NA_real_,
  nrow = n_2026,
  ncol = n_2021,
  dimnames = list(
    stages_2026,
    stages_2021
  )
)


pearson_mat <- matrix(
  NA_real_,
  nrow = n_2026,
  ncol = n_2021,
  dimnames = list(
    stages_2026,
    stages_2021
  )
)


shared_deg_mat <- matrix(
  0L,
  nrow = n_2026,
  ncol = n_2021,
  dimnames = list(
    stages_2026,
    stages_2021
  )
)


direction_mat <- matrix(
  NA_real_,
  nrow = n_2026,
  ncol = n_2021,
  dimnames = list(
    stages_2026,
    stages_2021
  )
)


pair_list <- list()
counter <- 1


for (s26 in stages_2026) {

  for (s21 in stages_2021) {

    x <- lfc_2026_common[
      ,
      s26
    ]

    y <- lfc_2021_common[
      ,
      s21
    ]


    keep <- is.finite(x) &
      is.finite(y)


    rho <- cor(
      x[keep],
      y[keep],
      method = "spearman"
    )


    r <- cor(
      x[keep],
      y[keep],
      method = "pearson"
    )


    deg26 <-
      !is.na(
        padj_2026_common[
          ,
          s26
        ]
      ) &
      padj_2026_common[
        ,
        s26
      ] < 0.05 &
      abs(
        x
      ) >= 1


    deg21 <-
      !is.na(
        padj_2021_common[
          ,
          s21
        ]
      ) &
      padj_2021_common[
        ,
        s21
      ] < 0.05 &
      abs(
        y
      ) >= 1


    shared <- deg26 &
      deg21


    n_shared <- sum(
      shared
    )


    if (n_shared > 0) {

      same_direction <- sum(
        sign(
          x[
            shared
          ]
        ) ==
          sign(
            y[
              shared
            ]
          )
      )


      direction_concordance <-
        same_direction /
        n_shared

    } else {

      same_direction <- 0L
      direction_concordance <- NA_real_
    }


    spearman_mat[
      s26,
      s21
    ] <- rho


    pearson_mat[
      s26,
      s21
    ] <- r


    shared_deg_mat[
      s26,
      s21
    ] <- n_shared


    direction_mat[
      s26,
      s21
    ] <- direction_concordance


    time26 <- unname(
      times_2026[
        s26
      ]
    )

    time21 <- unname(
      times_2021[
        s21
      ]
    )


    pair_list[counter] <- list(
      data.frame(
        Stage_2026 = s26,
        Time_2026_h = time26,

        Stage_2021 = s21,
        Time_2021_h = time21,

        Absolute_time_gap_h =
          abs(
            time26 -
              time21
          ),

        Spearman_rho =
          rho,

        Pearson_r =
          r,

        Shared_DEG =
          n_shared,

        Same_direction_shared_DEG =
          same_direction,

        Direction_concordance =
          direction_concordance,

        N_common_genes =
          sum(keep),

        stringsAsFactors = FALSE
      )
    )


    counter <- counter + 1
  }
}


pair_table <- do.call(
  rbind,
  pair_list
)

rownames(
  pair_table
) <- NULL


# ============================================================
# 9. Frozen Step09 validation
# ============================================================

spearman_values <-
  pair_table$Spearman_rho

pearson_values <-
  pair_table$Pearson_r


step09_summary <- data.frame(
  Metric = c(
    "Mean_Spearman",
    "Median_Spearman",
    "Minimum_Spearman",
    "Maximum_Spearman",
    "Mean_Pearson",
    "Median_Pearson"
  ),

  Observed = c(
    mean(
      spearman_values
    ),
    median(
      spearman_values
    ),
    min(
      spearman_values
    ),
    max(
      spearman_values
    ),
    mean(
      pearson_values
    ),
    median(
      pearson_values
    )
  ),

  Frozen_expected = c(
    0.37198,
    0.42311,
    -0.11664,
    0.72230,
    0.42502,
    0.48611
  ),

  stringsAsFactors = FALSE
)


step09_summary$Difference <-
  step09_summary$Observed -
  step09_summary$Frozen_expected


step09_summary$PASS <-
  abs(
    step09_summary$Difference
  ) < 0.002


cat(
  "Frozen Step09 concordance validation:\n"
)

print(
  step09_summary,
  row.names = FALSE,
  digits = 6
)

cat("\n")


if (!all(
  step09_summary$PASS
)) {

  warning(
    "One or more overall concordance metrics differ ",
    "from frozen Step09 values."
  )
}


# ============================================================
# 10. Validate max / min Spearman pair identities
# ============================================================

max_idx <- which.max(
  pair_table$Spearman_rho
)

min_idx <- which.min(
  pair_table$Spearman_rho
)


max_pair <- pair_table[
  max_idx,
  ,
  drop = FALSE
]

min_pair <- pair_table[
  min_idx,
  ,
  drop = FALSE
]


cat(
  "Maximum Spearman pair:\n"
)

print(
  max_pair,
  row.names = FALSE,
  digits = 5
)

cat("\n")


cat(
  "Minimum Spearman pair:\n"
)

print(
  min_pair,
  row.names = FALSE,
  digits = 5
)

cat("\n")


if (!(
  max_pair$Stage_2026 == "R72h" &&
    max_pair$Stage_2021 == "R48h"
)) {

  warning(
    "Maximum Spearman pair differs from frozen Step09."
  )
}


if (!(
  min_pair$Stage_2026 == "R1h" &&
    min_pair$Stage_2021 == "R48h"
)) {

  warning(
    "Minimum Spearman pair differs from frozen Step09."
  )
}


# ============================================================
# 11. Alternative-metric concordance summary
# ============================================================

metric_cor_pearson <- cor(
  pair_table$Spearman_rho,
  pair_table$Pearson_r,
  method = "pearson"
)

metric_cor_spearman <- cor(
  pair_table$Spearman_rho,
  pair_table$Pearson_r,
  method = "spearman"
)


metric_concordance_summary <- data.frame(
  N_stage_pairs =
    nrow(pair_table),

  Pearson_correlation_between_metrics =
    metric_cor_pearson,

  Spearman_correlation_between_metrics =
    metric_cor_spearman,

  Mean_absolute_metric_difference =
    mean(
      abs(
        pair_table$Pearson_r -
          pair_table$Spearman_rho
      )
    ),

  Maximum_absolute_metric_difference =
    max(
      abs(
        pair_table$Pearson_r -
          pair_table$Spearman_rho
      )
    ),

  stringsAsFactors = FALSE
)


# ============================================================
# 12. Prepare long-format heatmap data
# ============================================================

matrix_to_long <- function(
    mat,
    value_name) {

  x <- as.data.frame(
    as.table(
      mat
    ),
    stringsAsFactors = FALSE
  )

  colnames(x) <- c(
    "Stage_2026",
    "Stage_2021",
    value_name
  )


  x$Stage_2026 <- factor(
    x$Stage_2026,
    levels = rev(
      stages_2026
    )
  )

  x$Stage_2021 <- factor(
    x$Stage_2021,
    levels = stages_2021
  )


  x
}


pearson_long <- matrix_to_long(
  pearson_mat,
  "Pearson_r"
)

shared_long <- matrix_to_long(
  shared_deg_mat,
  "Shared_DEG"
)

direction_long <- matrix_to_long(
  direction_mat,
  "Direction_concordance"
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
      size = 8.2,
      colour = "grey35"
    ),

    axis.title = element_text(
      size = 9.5
    ),

    axis.text = element_text(
      size = 7.8
    ),

    legend.title = element_text(
      size = 8.3
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


heatmap_theme <- theme_minimal(
  base_size = 10
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
      size = 8.2,
      colour = "grey35"
    ),

    panel.grid = element_blank(),

    axis.title = element_text(
      size = 9
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.8
    ),

    axis.text.y = element_text(
      size = 7.8
    ),

    legend.title = element_text(
      size = 8.2
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
# 14. Panel A
# Complete Pearson response-concordance heatmap
# ============================================================

pA <- ggplot(
  pearson_long,
  aes(
    x = Stage_2021,
    y = Stage_2026,
    fill = Pearson_r
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.65
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Pearson_r
      )
    ),
    size = 2.45
  ) +

  scale_fill_gradient2(
    low = "#3B6FB6",
    mid = "white",
    high = "#C9504D",
    midpoint = 0,
    limits = c(
      -1,
      1
    ),
    name = "Pearson\nr"
  ) +

  scale_x_discrete(
    labels = unname(
      labels_2021[
        stages_2021
      ]
    )
  ) +

  scale_y_discrete(
    labels = function(x) {
      unname(
        labels_2026[
          as.character(x)
        ]
      )
    }
  ) +

  labs(
    tag = "A",

    title =
      "Pearson cross-cohort response concordance",

    subtitle =
      "Baseline-relative genome-wide log2FC vectors across 16,918 common genes",

    x =
      "2021 historical-reference stage",

    y =
      "2026 primary-cohort stage"
  ) +

  heatmap_theme


# ============================================================
# 15. Panel B
# Complete shared-DEG-number heatmap
# ============================================================

max_shared <- max(
  shared_long$Shared_DEG,
  na.rm = TRUE
)


pB <- ggplot(
  shared_long,
  aes(
    x = Stage_2021,
    y = Stage_2026,
    fill = Shared_DEG
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.65
  ) +

  geom_text(
    aes(
      label = Shared_DEG
    ),
    size = 2.35
  ) +

  scale_fill_gradient(
    low = "#F5F5F5",
    high = "#4C78A8",
    limits = c(
      0,
      max_shared
    ),
    name = "Shared\nDEGs"
  ) +

  scale_x_discrete(
    labels = unname(
      labels_2021[
        stages_2021
      ]
    )
  ) +

  scale_y_discrete(
    labels = function(x) {
      unname(
        labels_2026[
          as.character(x)
        ]
      )
    }
  ) +

  labs(
    tag = "B",

    title =
      "Shared differentially expressed genes",

    subtitle =
      "DEG defined independently in each cohort: adjusted P < 0.05 and |log2FC| \u2265 1",

    x =
      "2021 historical-reference stage",

    y =
      "2026 primary-cohort stage"
  ) +

  heatmap_theme


# ============================================================
# 16. Panel C
# Shared-DEG direction-concordance heatmap
# ============================================================

direction_min <- min(
  direction_long$Direction_concordance,
  na.rm = TRUE
)

direction_floor <- max(
  0,
  floor(
    direction_min * 20
  ) / 20
)


pC <- ggplot(
  direction_long,
  aes(
    x = Stage_2021,
    y = Stage_2026,
    fill = Direction_concordance
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.65
  ) +

  geom_text(
    aes(
      label = ifelse(
        is.na(
          Direction_concordance
        ),
        "NA",
        sprintf(
          "%.1f%%",
          100 *
            Direction_concordance
        )
      )
    ),
    size = 2.25
  ) +

  scale_fill_gradient(
    low = "#F3F3F3",
    high = "#4E8C62",
    limits = c(
      direction_floor,
      1
    ),
    name = "Direction\nconcordance"
  ) +

  scale_x_discrete(
    labels = unname(
      labels_2021[
        stages_2021
      ]
    )
  ) +

  scale_y_discrete(
    labels = function(x) {
      unname(
        labels_2026[
          as.character(x)
        ]
      )
    }
  ) +

  labs(
    tag = "C",

    title =
      "Direction concordance among shared DEGs",

    subtitle =
      "Fraction of shared DEGs with matching log2FC direction",

    x =
      "2021 historical-reference stage",

    y =
      "2026 primary-cohort stage"
  ) +

  heatmap_theme


# ============================================================
# 17. Panel D
# Spearman versus Pearson across all 70 pairs
# ============================================================

pair_table$Pair_label <- paste0(
  unname(
    labels_2026[
      pair_table$Stage_2026
    ]
  ),
  " / ",
  unname(
    labels_2021[
      pair_table$Stage_2021
    ]
  )
)


metric_label <- paste0(
  "Across 70 stage pairs: Pearson r = ",
  sprintf(
    "%.3f",
    metric_cor_pearson
  ),
  "; Spearman rho = ",
  sprintf(
    "%.3f",
    metric_cor_spearman
  )
)


pD <- ggplot(
  pair_table,
  aes(
    x = Spearman_rho,
    y = Pearson_r,
    colour = Absolute_time_gap_h
  )
) +

  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey55"
  ) +

  geom_point(
    size = 2.6,
    alpha = 0.85
  ) +

  scale_colour_gradient(
    low = "#5B8DB8",
    high = "#C85A5A",
    name =
      "Absolute time\ngap (h)"
  ) +

  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = metric_label,
    hjust = -0.03,
    vjust = 1.08,
    size = 2.7,
    fill = "white",
    colour = "grey20",
    linewidth = 0.25
  ) +

  labs(
    tag = "D",

    title =
      "Concordance between correlation metrics",

    subtitle =
      "Each point represents one of the 10 \u00d7 7 cross-cohort stage comparisons",

    x =
      "Spearman response concordance",

    y =
      "Pearson response concordance"
  ) +

  theme_pub +

  theme(
    legend.position = "right"
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
  "Supplementary_Figure_S7A_Pearson_concordance"
)

save_panel(
  pB,
  "Supplementary_Figure_S7B_shared_DEG_number"
)

save_panel(
  pC,
  "Supplementary_Figure_S7C_direction_concordance"
)

save_panel(
  pD,
  "Supplementary_Figure_S7D_Spearman_vs_Pearson"
)


# ============================================================
# 20. Combined S7
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
# 21. Save final publication figure
# ============================================================

pdf_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S7_",
    "cross_cohort_complete_diagnostics.pdf"
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
    "Supplementary_Figure_S7_",
    "cross_cohort_complete_diagnostics_600dpi.tiff"
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
# 22. Save all source data
# ============================================================

write.table(
  pair_table,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7_all_70_stage_pair_diagnostics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  pearson_mat,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7A_Pearson_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  shared_deg_mat,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7B_shared_DEG_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  direction_mat,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7C_direction_concordance_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  spearman_mat,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7_primary_Spearman_matrix_for_reference.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  step09_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7_Step09_validation.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  metric_concordance_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S7_correlation_metric_concordance.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


writeLines(
  common_genes,
  con = file.path(
    outdir,
    "Supplementary_Figure_S7_common_16918_geneIDs.txt"
  )
)


# ============================================================
# 23. Caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S7 | Complete cross-cohort diagnostics of ",
  "Baseline-relative transcriptomic responses in the 2026 primary ",
  "cohort and the independently reprocessed 2021 historical-reference ",
  "cohort. ",
  "a, Pearson correlations between genome-wide Baseline-relative log2 ",
  "fold-change vectors for all 70 combinations of the ten 2026 ",
  "post-injury stages and seven 2021 post-injury stages, calculated ",
  "across 16,918 genes retained after independent count filtering in ",
  "both cohorts. ",
  "b, Numbers of shared differentially expressed genes for each ",
  "cross-cohort stage pair. Differentially expressed genes were defined ",
  "independently within each cohort using adjusted P < 0.05 and an ",
  "absolute log2 fold change of at least 1. ",
  "c, Direction concordance among shared differentially expressed genes, ",
  "defined as the proportion showing the same sign of log2 fold change ",
  "in both cohorts. ",
  "d, Comparison of Spearman and Pearson genome-wide response ",
  "concordance across all 70 stage pairs; colours denote the absolute ",
  "difference in sampling time between the paired stages. ",
  "The two cohorts were analysed independently relative to their own ",
  "uninjured Baseline groups and were not pooled or batch-corrected. ",
  "Pairwise stage comparisons are descriptive and are not statistically ",
  "independent; unconstrained high-concordance pairs should therefore ",
  "not be interpreted as evidence that the corresponding chronological ",
  "stages are biologically equivalent."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S7_caption.txt"
  )
)


# ============================================================
# 24. Interpretation note
# ============================================================

interpretation_text <- c(
  "Supplementary Figure S7 interpretation notes",
  "",
  "1. Spearman correlation remains the primary cross-cohort response",
  "   concordance metric used in the main analysis.",
  "",
  "2. Pearson correlation is presented here as an alternative metric",
  "   to demonstrate that the overall cross-cohort pattern is not",
  "   dependent on rank-based correlation alone.",
  "",
  "3. Shared-DEG counts and direction concordance are auxiliary",
  "   threshold-dependent diagnostics.",
  "",
  "4. Each cohort was analysed independently relative to its own",
  "   Baseline group.",
  "",
  "5. No pooling, ComBat correction, or cross-cohort batch correction",
  "   was performed.",
  "",
  "6. The 70 stage-pair cells are dependent because individual",
  "   stage-response vectors participate in multiple comparisons.",
  "",
  "7. High concordance between chronologically distant stages does not",
  "   establish temporal equivalence.",
  "",
  "8. Chronologically constrained comparisons in the main figure should",
  "   remain the preferred interpretation of cross-cohort temporal",
  "   correspondence."
)

writeLines(
  interpretation_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S7_interpretation_note.txt"
  )
)


# ============================================================
# 25. Plot objects / session info
# ============================================================

saveRDS(
  list(
    S7A = pA,
    S7B = pB,
    S7C = pC,
    S7D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S7_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "Supplementary_Figure_S7_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 26. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S7 COMPLETED\n")
cat("============================================================\n\n")


cat(
  "Common independently filtered genes : ",
  length(common_genes),
  "\n",
  sep = ""
)

cat(
  "Cross-cohort stage pairs            : ",
  nrow(pair_table),
  "\n\n",
  sep = ""
)


cat("Overall concordance summary:\n")

print(
  step09_summary,
  row.names = FALSE,
  digits = 6
)


cat("\nMaximum Spearman pair:\n")

print(
  max_pair,
  row.names = FALSE,
  digits = 5
)


cat("\nMinimum Spearman pair:\n")

print(
  min_pair,
  row.names = FALSE,
  digits = 5
)


cat("\nSpearman-versus-Pearson metric concordance:\n")

print(
  metric_concordance_summary,
  row.names = FALSE,
  digits = 5
)


cat("\nShared-DEG range across 70 pairs:\n")

cat(
  "Minimum : ",
  min(
    pair_table$Shared_DEG
  ),
  "\n",
  sep = ""
)

cat(
  "Median  : ",
  median(
    pair_table$Shared_DEG
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum : ",
  max(
    pair_table$Shared_DEG
  ),
  "\n",
  sep = ""
)


cat("\nDirection-concordance range across 70 pairs:\n")

cat(
  "Minimum : ",
  sprintf(
    "%.4f",
    min(
      pair_table$Direction_concordance,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Median  : ",
  sprintf(
    "%.4f",
    median(
      pair_table$Direction_concordance,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum : ",
  sprintf(
    "%.4f",
    max(
      pair_table$Direction_concordance,
      na.rm = TRUE
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
