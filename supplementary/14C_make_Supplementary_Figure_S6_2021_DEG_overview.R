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
cat("14C SUPPLEMENTARY FIGURE S6\n")
cat("2021 historical-reference cohort DEG overview\n")
cat("Timepoint-versus-Baseline transcriptomic responses\n")
cat("============================================================\n\n")


# ============================================================
# 0. Output directory
# ============================================================

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S6_2021_DEG_overview"
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
# 1. Frozen 2021 contrast definitions
# ============================================================

stage_order <- c(
  "R4h",
  "R8h",
  "R12h",
  "R16h",
  "R20h",
  "R24h",
  "R48h"
)

stage_labels <- c(
  R4h  = "4 h",
  R8h  = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)

stage_colors <- c(
  R4h  = "#4C78A8",
  R8h  = "#5F8DB8",
  R12h = "#72A1B3",
  R16h = "#7FB28F",
  R20h = "#A3B76B",
  R24h = "#D2A455",
  R48h = "#C85A5A"
)


# ============================================================
# 2. Frozen Step08 DEG summary
#
# Threshold:
# padj < 0.05 and |log2FC| >= 1
# ============================================================

deg_summary <- data.frame(
  Stage = stage_order,

  Up = c(
    2363,
    2516,
    1729,
    1549,
    978,
    1781,
    2173
  ),

  Down = c(
    1675,
    2314,
    963,
    1368,
    538,
    1847,
    1403
  ),

  stringsAsFactors = FALSE
)

deg_summary$Total <-
  deg_summary$Up +
  deg_summary$Down

expected_total <- c(
  4038,
  4830,
  2692,
  2917,
  1516,
  3628,
  3576
)

if (!identical(
  as.integer(deg_summary$Total),
  as.integer(expected_total)
)) {
  stop(
    "Internal frozen DEG-summary check failed."
  )
}


# ============================================================
# 3. Helper: locate existing Step08 files
# ============================================================

locate_file <- function(
    preferred_names,
    regex_pattern = NULL,
    required = TRUE) {

  for (candidate in preferred_names) {

    if (file.exists(candidate)) {
      return(
        normalizePath(
          candidate,
          mustWork = TRUE
        )
      )
    }
  }


  all_files <- list.files(
    ".",
    recursive = TRUE,
    full.names = TRUE
  )


  base_names <- basename(
    all_files
  )


  for (preferred in preferred_names) {

    hit <- which(
      base_names ==
        basename(preferred)
    )

    if (length(hit) == 1) {

      return(
        normalizePath(
          all_files[hit],
          mustWork = TRUE
        )
      )
    }

    if (length(hit) > 1) {

      stop(
        "Multiple files found with basename: ",
        basename(preferred),
        "\n",
        paste(
          all_files[hit],
          collapse = "\n"
        )
      )
    }
  }


  if (!is.null(
    regex_pattern
  )) {

    hit <- grep(
      regex_pattern,
      base_names,
      ignore.case = TRUE
    )

    if (length(hit) == 1) {

      return(
        normalizePath(
          all_files[hit],
          mustWork = TRUE
        )
      )
    }

    if (length(hit) > 1) {

      cat(
        "Multiple regex candidate files found:\n"
      )

      cat(
        paste(
          all_files[hit],
          collapse = "\n"
        ),
        "\n"
      )

      stop(
        "Please retain a unique Step08 source file."
      )
    }
  }


  if (required) {

    stop(
      "Required Step08 file not found."
    )

  } else {

    return(
      NA_character_
    )
  }
}


# ============================================================
# 4. Locate genome-wide log2FC matrix
# ============================================================

lfc_file <- locate_file(
  preferred_names = c(
    "08_2021_genomewide_log2FC_matrix.tsv",
    paste0(
      "08_2021_historical_reference_reanalysis_result/",
      "08_2021_genomewide_log2FC_matrix.tsv"
    )
  ),
  regex_pattern =
    "^08_2021.*genomewide.*log2fc.*matrix.*\\.tsv$",
  required = TRUE
)


# ============================================================
# 5. Locate PADJ matrix if available
# ============================================================

padj_file <- locate_file(
  preferred_names = c(
    "08_2021_genomewide_padj_matrix.tsv",
    "08_2021_genomewide_PADJ_matrix.tsv",
    paste0(
      "08_2021_historical_reference_reanalysis_result/",
      "08_2021_genomewide_padj_matrix.tsv"
    ),
    paste0(
      "08_2021_historical_reference_reanalysis_result/",
      "08_2021_genomewide_PADJ_matrix.tsv"
    )
  ),
  regex_pattern =
    "^08_2021.*genomewide.*padj.*matrix.*\\.tsv$",
  required = FALSE
)


cat(
  "Genome-wide log2FC matrix:\n",
  lfc_file,
  "\n\n",
  sep = ""
)

if (!is.na(padj_file)) {

  cat(
    "Genome-wide PADJ matrix:\n",
    padj_file,
    "\n\n",
    sep = ""
  )

} else {

  cat(
    "Genome-wide PADJ matrix not found.\n",
    "Panel A will use the frozen Step08 DEG summary.\n\n",
    sep = ""
  )
}


# ============================================================
# 6. Read log2FC matrix
# ============================================================

lfc_df <- read.delim(
  lfc_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (ncol(lfc_df) < 8) {
  stop(
    "log2FC matrix has too few columns."
  )
}

colnames(lfc_df)[1] <- "GeneID"


if (anyDuplicated(
  lfc_df$GeneID
)) {
  stop(
    "Duplicated GeneID in log2FC matrix."
  )
}


# ============================================================
# 7. Detect contrast columns
# ============================================================

missing_stage_columns <- setdiff(
  stage_order,
  colnames(lfc_df)
)

if (length(missing_stage_columns) > 0) {

  cat(
    "Available log2FC columns:\n"
  )

  cat(
    paste(
      colnames(lfc_df),
      collapse = "\n"
    ),
    "\n"
  )

  stop(
    "Missing expected log2FC contrast column(s): ",
    paste(
      missing_stage_columns,
      collapse = ", "
    )
  )
}


lfc_mat <- as.matrix(
  lfc_df[
    ,
    stage_order,
    drop = FALSE
  ]
)

storage.mode(
  lfc_mat
) <- "numeric"

rownames(lfc_mat) <-
  lfc_df$GeneID


if (nrow(lfc_mat) != 17328) {

  warning(
    "Expected 17,328 genes in 2021 log2FC matrix; observed ",
    nrow(lfc_mat)
  )
}


cat(
  "log2FC matrix : ",
  nrow(lfc_mat),
  " genes x ",
  ncol(lfc_mat),
  " contrasts\n",
  sep = ""
)


cat(
  "NA log2FC values : ",
  sum(
    is.na(lfc_mat)
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 8. Optional DEG-count validation from PADJ matrix
# ============================================================

deg_validation <- NULL

if (!is.na(
  padj_file
)) {

  padj_df <- read.delim(
    padj_file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  colnames(padj_df)[1] <- "GeneID"


  if (anyDuplicated(
    padj_df$GeneID
  )) {
    stop(
      "Duplicated GeneID in PADJ matrix."
    )
  }


  missing_padj_cols <- setdiff(
    stage_order,
    colnames(padj_df)
  )

  if (length(
    missing_padj_cols
  ) > 0) {
    stop(
      "PADJ matrix lacks contrast column(s): ",
      paste(
        missing_padj_cols,
        collapse = ", "
      )
    )
  }


  if (!setequal(
    padj_df$GeneID,
    lfc_df$GeneID
  )) {
    stop(
      "Gene sets differ between log2FC and PADJ matrices."
    )
  }


  padj_df <- padj_df[
    match(
      lfc_df$GeneID,
      padj_df$GeneID
    ),
    ,
    drop = FALSE
  ]


  padj_mat <- as.matrix(
    padj_df[
      ,
      stage_order,
      drop = FALSE
    ]
  )

  storage.mode(
    padj_mat
  ) <- "numeric"


  validation_list <- vector(
    "list",
    length(stage_order)
  )


  for (i in seq_along(
    stage_order
  )) {

    stage <- stage_order[i]

    lfc <- lfc_mat[
      ,
      stage
    ]

    padj <- padj_mat[
      ,
      stage
    ]


    up <- sum(
      !is.na(padj) &
        padj < 0.05 &
        lfc >= 1
    )

    down <- sum(
      !is.na(padj) &
        padj < 0.05 &
        lfc <= -1
    )


    validation_list[i] <- list(
      data.frame(
        Stage = stage,
        Recalculated_Up = up,
        Recalculated_Down = down,
        Recalculated_Total =
          up + down,

        Frozen_Up =
          deg_summary$Up[
            deg_summary$Stage == stage
          ],

        Frozen_Down =
          deg_summary$Down[
            deg_summary$Stage == stage
          ],

        Frozen_Total =
          deg_summary$Total[
            deg_summary$Stage == stage
          ],

        stringsAsFactors = FALSE
      )
    )
  }


  deg_validation <- do.call(
    rbind,
    validation_list
  )


  deg_validation$PASS <- with(
    deg_validation,

    Recalculated_Up ==
      Frozen_Up &

      Recalculated_Down ==
      Frozen_Down &

      Recalculated_Total ==
      Frozen_Total
  )


  cat(
    "\nDEG validation against frozen Step08 summary:\n"
  )

  print(
    deg_validation,
    row.names = FALSE
  )

  cat("\n")


  if (!all(
    deg_validation$PASS
  )) {

    stop(
      "Recalculated DEG counts do not match frozen Step08 results."
    )
  }


  cat(
    "Frozen DEG count validation: PASS\n\n"
  )
}


# ============================================================
# 9. Panel A source data
# ============================================================

deg_long <- rbind(
  data.frame(
    Stage =
      deg_summary$Stage,
    Direction =
      "Up",
    Count =
      deg_summary$Up,
    Signed_count =
      deg_summary$Up,
    stringsAsFactors = FALSE
  ),

  data.frame(
    Stage =
      deg_summary$Stage,
    Direction =
      "Down",
    Count =
      deg_summary$Down,
    Signed_count =
      -deg_summary$Down,
    stringsAsFactors = FALSE
  )
)

deg_long$Stage <- factor(
  deg_long$Stage,
  levels = stage_order
)

deg_long$Direction <- factor(
  deg_long$Direction,
  levels = c(
    "Up",
    "Down"
  )
)


# ============================================================
# 10. Genome-wide absolute response magnitude
# ============================================================

response_summary_list <- vector(
  "list",
  length(stage_order)
)


for (i in seq_along(
  stage_order
)) {

  stage <- stage_order[i]

  x <- abs(
    lfc_mat[
      ,
      stage
    ]
  )

  x <- x[
    is.finite(x)
  ]


  response_summary_list[i] <- list(
    data.frame(
      Stage = stage,
      N_genes =
        length(x),

      Q25_abs_log2FC =
        unname(
          quantile(
            x,
            0.25
          )
        ),

      Median_abs_log2FC =
        median(x),

      Q75_abs_log2FC =
        unname(
          quantile(
            x,
            0.75
          )
        ),

      P90_abs_log2FC =
        unname(
          quantile(
            x,
            0.90
          )
        ),

      RMS_log2FC =
        sqrt(
          mean(
            x^2
          )
        ),

      stringsAsFactors = FALSE
    )
  )
}


response_summary <- do.call(
  rbind,
  response_summary_list
)

response_summary$Stage <- factor(
  response_summary$Stage,
  levels = stage_order
)


# ============================================================
# 11. Genome-wide response-vector correlations
# ============================================================

response_cor <- cor(
  lfc_mat,
  use = "pairwise.complete.obs",
  method = "spearman"
)


response_cor_long <- as.data.frame(
  as.table(
    response_cor
  ),
  stringsAsFactors = FALSE
)

colnames(
  response_cor_long
) <- c(
  "Stage_Y",
  "Stage_X",
  "Spearman_rho"
)

response_cor_long$Stage_X <- factor(
  response_cor_long$Stage_X,
  levels = stage_order
)

response_cor_long$Stage_Y <- factor(
  response_cor_long$Stage_Y,
  levels = rev(
    stage_order
  )
)


# ============================================================
# 12. Adjacent-stage response-vector correlations
# ============================================================

adjacent_cor_list <- vector(
  "list",
  length(stage_order) - 1
)


for (i in seq_len(
  length(stage_order) - 1
)) {

  s1 <- stage_order[i]
  s2 <- stage_order[i + 1]

  x <- lfc_mat[
    ,
    s1
  ]

  y <- lfc_mat[
    ,
    s2
  ]

  keep <- is.finite(x) &
    is.finite(y)


  rho <- cor(
    x[keep],
    y[keep],
    method = "spearman"
  )


  adjacent_cor_list[i] <- list(
    data.frame(
      Stage_1 = s1,
      Stage_2 = s2,

      Comparison = paste0(
        unname(
          stage_labels[s1]
        ),
        " \u2194 ",
        unname(
          stage_labels[s2]
        )
      ),

      Order = i,

      N_genes =
        sum(keep),

      Spearman_rho =
        rho,

      stringsAsFactors = FALSE
    )
  )
}


adjacent_cor <- do.call(
  rbind,
  adjacent_cor_list
)

adjacent_cor$Comparison <- factor(
  adjacent_cor$Comparison,
  levels = adjacent_cor$Comparison
)


# ============================================================
# 13. Global response-correlation summary
# ============================================================

pairwise_cor_values <- response_cor[
  upper.tri(
    response_cor
  )
]

response_cor_summary <- data.frame(
  N_stage_pairs =
    length(
      pairwise_cor_values
    ),

  Minimum_Spearman =
    min(
      pairwise_cor_values
    ),

  Median_Spearman =
    median(
      pairwise_cor_values
    ),

  Mean_Spearman =
    mean(
      pairwise_cor_values
    ),

  Maximum_Spearman =
    max(
      pairwise_cor_values
    ),

  stringsAsFactors = FALSE
)


# ============================================================
# 14. Publication theme
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
      size = 8
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 15. Panel A
# Timepoint-vs-Baseline DEG counts
# ============================================================

pA <- ggplot(
  deg_long,
  aes(
    x = Stage,
    y = Signed_count,
    fill = Direction
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    colour = "grey35"
  ) +

  geom_col(
    width = 0.66,
    position = "identity",
    colour = "grey25",
    linewidth = 0.25
  ) +

  geom_text(
    aes(
      label = Count
    ),
    vjust = ifelse(
      deg_long$Signed_count >= 0,
      -0.35,
      1.25
    ),
    size = 2.8
  ) +

  scale_fill_manual(
    values = c(
      Up = "#C85A5A",
      Down = "#4C78A8"
    ),
    name = "Direction"
  ) +

  scale_x_discrete(
    labels = unname(
      stage_labels[
        stage_order
      ]
    )
  ) +

  scale_y_continuous(
    labels = function(x) {
      abs(x)
    },
    expand = expansion(
      mult = c(
        0.10,
        0.10
      )
    )
  ) +

  labs(
    tag = "A",

    title =
      "Differential expression relative to Baseline",

    subtitle =
      "DESeq2; adjusted P < 0.05 and |log2 fold change| \u2265 1",

    x = NULL,

    y =
      "Number of differentially expressed genes"
  ) +

  theme_pub +

  theme(
    legend.position = "top"
  )


# ============================================================
# 16. Panel B
# Genome-wide response magnitude
# ============================================================

pB <- ggplot(
  response_summary,
  aes(
    x = Stage,
    y = Median_abs_log2FC
  )
) +

  geom_errorbar(
    aes(
      ymin = Q25_abs_log2FC,
      ymax = Q75_abs_log2FC
    ),
    width = 0.12,
    linewidth = 0.75,
    colour = "grey35"
  ) +

  geom_line(
    aes(
      group = 1
    ),
    colour = "grey45",
    linewidth = 0.75
  ) +

  geom_point(
    aes(
      fill = Stage
    ),
    shape = 21,
    size = 3.4,
    stroke = 0.6,
    colour = "black"
  ) +

  geom_point(
    aes(
      y = P90_abs_log2FC
    ),
    shape = 17,
    size = 2.6,
    colour = "grey25"
  ) +

  scale_fill_manual(
    values = stage_colors,
    guide = "none"
  ) +

  scale_x_discrete(
    labels = unname(
      stage_labels[
        stage_order
      ]
    )
  ) +

  labs(
    tag = "B",

    title =
      "Genome-wide baseline-relative response magnitude",

    subtitle =
      "Circles: median |log2FC| with IQR; triangles: 90th percentile",

    x = NULL,

    y =
      "Absolute log2 fold change"
  ) +

  theme_pub


# ============================================================
# 17. Panel C
# Genome-wide response-vector correlation
# ============================================================

pC <- ggplot(
  response_cor_long,
  aes(
    x = Stage_X,
    y = Stage_Y,
    fill = Spearman_rho
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman_rho
      )
    ),
    size = 2.8
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
    name = "Spearman\nrho"
  ) +

  scale_x_discrete(
    labels = unname(
      stage_labels[
        stage_order
      ]
    )
  ) +

  scale_y_discrete(
    labels = function(x) {
      unname(
        stage_labels[
          as.character(x)
        ]
      )
    }
  ) +

  labs(
    tag = "C",

    title =
      "Similarity of genome-wide injury-response vectors",

    subtitle =
      "Spearman correlation of Baseline-relative log2FC across 17,328 genes",

    x = "Sampling stage",
    y = "Sampling stage"
  ) +

  theme_minimal(
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
      size = 8.3,
      colour = "grey35"
    ),

    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7.8
    ),

    axis.text.y = element_text(
      size = 7.8
    ),

    axis.title = element_text(
      size = 9.5
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
# 18. Panel D
# Adjacent response-vector correlation
# ============================================================

pD <- ggplot(
  adjacent_cor,
  aes(
    x = Comparison,
    y = Spearman_rho,
    group = 1
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +

  geom_line(
    colour = "grey35",
    linewidth = 0.8
  ) +

  geom_point(
    size = 3.0,
    shape = 21,
    fill = "white",
    stroke = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf(
        "%.2f",
        Spearman_rho
      )
    ),
    vjust = -0.75,
    size = 2.8
  ) +

  scale_y_continuous(
    limits = c(
      min(
        0,
        min(
          adjacent_cor$Spearman_rho
        ) - 0.12
      ),
      min(
        1,
        max(
          adjacent_cor$Spearman_rho
        ) + 0.12
      )
    )
  ) +

  labs(
    tag = "D",

    title =
      "Similarity between consecutive injury-response vectors",

    subtitle =
      "Genome-wide Spearman correlation between adjacent sampled stages",

    x = NULL,

    y =
      "Spearman correlation"
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
# 19. TIFF helper
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
# 20. Save individual panels
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
  "Supplementary_Figure_S6A_2021_DEG_counts"
)

save_panel(
  pB,
  "Supplementary_Figure_S6B_2021_response_magnitude"
)

save_panel(
  pC,
  "Supplementary_Figure_S6C_2021_response_correlation"
)

save_panel(
  pD,
  "Supplementary_Figure_S6D_2021_adjacent_response_correlation"
)


# ============================================================
# 21. Combined S6
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
# 22. Save final publication figure
# ============================================================

pdf_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Figure_S6_",
    "2021_historical_reference_DEG_overview.pdf"
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
    "Supplementary_Figure_S6_",
    "2021_historical_reference_DEG_overview_600dpi.tiff"
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
# 23. Source data
# ============================================================

write.table(
  deg_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S6A_DEG_counts.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  response_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S6B_response_magnitude.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  response_cor,
  file = file.path(
    outdir,
    "Supplementary_Figure_S6C_response_Spearman_matrix.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


write.table(
  adjacent_cor,
  file = file.path(
    outdir,
    "Supplementary_Figure_S6D_adjacent_response_correlations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  response_cor_summary,
  file = file.path(
    outdir,
    "Supplementary_Figure_S6_response_correlation_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


if (!is.null(
  deg_validation
)) {

  write.table(
    deg_validation,
    file = file.path(
      outdir,
      "Supplementary_Figure_S6_DEG_count_validation.tsv"
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}


# ============================================================
# 24. Caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S6 | Baseline-relative differential-expression ",
  "responses in the independently reprocessed 2021 historical-reference ",
  "cohort. ",
  "a, Numbers of up- and downregulated genes at each post-injury stage ",
  "relative to the corresponding uninjured Baseline group. Differentially ",
  "expressed genes were defined using DESeq2 with adjusted P < 0.05 and ",
  "an absolute log2 fold change of at least 1. ",
  "b, Genome-wide magnitude of Baseline-relative expression responses ",
  "across 17,328 genes. Circles indicate the median absolute log2 fold ",
  "change, error bars indicate the interquartile range, and triangles ",
  "indicate the 90th percentile. ",
  "c, Pairwise Spearman correlations between genome-wide Baseline-relative ",
  "log2 fold-change vectors for the seven post-injury stages. ",
  "d, Spearman correlations between genome-wide response vectors for ",
  "consecutive sampled stages. These analyses describe the internal ",
  "temporal organization of the historical-reference cohort and are not ",
  "interpreted as evidence of temporal equivalence with the 2026 primary ",
  "cohort. The 2021 dataset was analysed independently and was not pooled ",
  "or batch-corrected with the 2026 cohort."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S6_caption.txt"
  )
)


# ============================================================
# 25. Plot objects / session info
# ============================================================

saveRDS(
  list(
    S6A = pA,
    S6B = pB,
    S6C = pC,
    S6D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S6_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "Supplementary_Figure_S6_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 26. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S6 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Genes in genome-wide response matrix : ",
  nrow(lfc_mat),
  "\n",
  sep = ""
)

cat(
  "Post-injury contrasts                : ",
  ncol(lfc_mat),
  "\n\n",
  sep = ""
)


cat("Frozen DEG summary:\n")

print(
  deg_summary,
  row.names = FALSE
)


cat("\nGenome-wide response magnitude:\n")

print(
  response_summary,
  row.names = FALSE,
  digits = 5
)


cat("\nGenome-wide response correlation summary:\n")

print(
  response_cor_summary,
  row.names = FALSE,
  digits = 5
)


cat("\nAdjacent-stage response correlations:\n")

print(
  adjacent_cor,
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
