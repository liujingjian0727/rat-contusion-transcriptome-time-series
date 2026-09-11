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
cat("12A SUPPLEMENTARY FIGURE S1\n")
cat("Gene-body coverage technical validation\n")
cat("Final 77-sample 2026 primary cohort\n")
cat("============================================================\n\n")


# ============================================================
# 0. Portable data-root configuration
# ============================================================

DATA_ROOT <- Sys.getenv("SCIDATA_DATA_ROOT")

if (DATA_ROOT == "") {
  stop(
    paste0(
      "Environment variable SCIDATA_DATA_ROOT is not set.\n",
      "Set it to the local bulk_RNA_2026 directory, for example:\n",
      "export SCIDATA_DATA_ROOT=/path/to/bulk_RNA_2026"
    )
  )
}

DATA_ROOT <- normalizePath(
  DATA_ROOT,
  mustWork = TRUE
)

gene_body_dir <- file.path(
  DATA_ROOT,
  "GeneBodyCoverage"
)

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S1_gene_body_coverage"
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
# 1. Stage definitions
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
# 2. Input validation
# ============================================================

if (!dir.exists(gene_body_dir)) {
  stop(
    "GeneBodyCoverage directory not found:\n",
    gene_body_dir
  )
}

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found:\n",
    metadata_file
  )
}

cat(
  "Gene-body directory : ",
  gene_body_dir,
  "\n",
  sep = ""
)

cat(
  "Metadata whitelist  : ",
  metadata_file,
  "\n\n",
  sep = ""
)


# ============================================================
# 3. Read final 77-sample metadata
# ============================================================

meta <- read.delim(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

required_meta_cols <- c(
  "sampleID",
  "group"
)

if (!all(
  required_meta_cols %in%
    colnames(meta)
)) {
  stop(
    "Metadata must contain sampleID and group."
  )
}

if (anyDuplicated(meta$sampleID)) {
  stop(
    "Duplicated sampleID in metadata."
  )
}

if (nrow(meta) != 77) {
  stop(
    "Expected exactly 77 final samples; observed ",
    nrow(meta)
  )
}

unknown_groups <- setdiff(
  unique(meta$group),
  groups
)

if (length(unknown_groups) > 0) {
  stop(
    "Unexpected group(s) in metadata: ",
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

final_samples <- meta$sampleID

cat(
  "Final metadata samples : ",
  length(final_samples),
  "\n",
  sep = ""
)

cat(
  "Final R2h samples      : ",
  paste(
    meta$sampleID[
      meta$group == "R2h"
    ],
    collapse = ", "
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 4. Discover all RSeQC GeneBodyCoverage files
# ============================================================

coverage_files <- list.files(
  gene_body_dir,
  pattern = "\\.geneBodyCoverage\\.txt$",
  full.names = TRUE,
  recursive = FALSE
)

if (length(coverage_files) == 0) {
  stop(
    "No *.geneBodyCoverage.txt files found."
  )
}

file_samples <- sub(
  "\\.geneBodyCoverage\\.txt$",
  "",
  basename(coverage_files)
)

if (anyDuplicated(file_samples)) {
  dup <- unique(
    file_samples[
      duplicated(file_samples)
    ]
  )

  stop(
    "Duplicated GeneBodyCoverage sample files: ",
    paste(
      dup,
      collapse = ", "
    )
  )
}

cat(
  "Coverage files discovered : ",
  length(coverage_files),
  "\n",
  sep = ""
)


# ============================================================
# 5. Compare coverage files against final 77-sample whitelist
# ============================================================

missing_final_samples <- setdiff(
  final_samples,
  file_samples
)

extra_coverage_samples <- setdiff(
  file_samples,
  final_samples
)

cat(
  "Final samples missing coverage files : ",
  length(missing_final_samples),
  "\n",
  sep = ""
)

if (length(missing_final_samples) > 0) {

  cat(
    "  ",
    paste(
      missing_final_samples,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Extra coverage files outside final cohort : ",
  length(extra_coverage_samples),
  "\n",
  sep = ""
)

if (length(extra_coverage_samples) > 0) {

  cat(
    "  ",
    paste(
      extra_coverage_samples,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat("\n")

if (length(missing_final_samples) > 0) {
  stop(
    "One or more final 77 samples lack ",
    "gene-body coverage results."
  )
}


# ============================================================
# 6. Build sample audit table
# ============================================================

audit_samples <- union(
  final_samples,
  file_samples
)

audit <- data.frame(
  sampleID = audit_samples,
  In_final_77_metadata =
    audit_samples %in%
      final_samples,
  GeneBodyCoverage_file_present =
    audit_samples %in%
      file_samples,
  stringsAsFactors = FALSE
)

audit$Status <- ifelse(
  audit$In_final_77_metadata &
    audit$GeneBodyCoverage_file_present,
  "Included_final77",
  ifelse(
    !audit$In_final_77_metadata &
      audit$GeneBodyCoverage_file_present,
    "Excluded_not_in_final_metadata",
    "Missing_coverage_file"
  )
)

write.table(
  audit,
  file = file.path(
    outdir,
    "S1_sample_whitelist_audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 7. RSeQC GeneBodyCoverage parser
# Standard format:
# Percentile  1 2 ... 100
# sample      coverage1 ... coverage100
# ============================================================

read_gene_body_file <- function(
    file,
    expected_sample) {

  # First try standard tab-delimited RSeQC format
  x <- try(
    read.delim(
      file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    ),
    silent = TRUE
  )

  # Fallback to generic whitespace-separated format
  if (inherits(x, "try-error") ||
      ncol(x) < 2 ||
      nrow(x) < 1) {

    x <- read.table(
      file,
      header = TRUE,
      sep = "",
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    )
  }


  if (nrow(x) < 1) {
    stop(
      "No coverage data found in: ",
      file
    )
  }


  positions <- suppressWarnings(
    as.numeric(
      colnames(x)[-1]
    )
  )

  values <- suppressWarnings(
    as.numeric(
      x[
        1,
        -1,
        drop = TRUE
      ]
    )
  )

  embedded_sample <- as.character(
    x[
      1,
      1
    ]
  )


  if (length(positions) !=
      length(values)) {

    stop(
      "Position/value length mismatch in ",
      file
    )
  }


  if (length(values) != 100) {

    stop(
      "Expected 100 gene-body percentile values in ",
      file,
      "; observed ",
      length(values)
    )
  }


  if (anyNA(positions) ||
      anyNA(values) ||
      any(!is.finite(values))) {

    stop(
      "Invalid numeric values in ",
      file
    )
  }


  if (!identical(
    as.integer(positions),
    1:100
  )) {

    stop(
      "Expected percentiles 1:100 in ",
      file
    )
  }


  # Filename is treated as canonical sample ID.
  # Embedded RSeQC label is retained only for audit.
  data.frame(
    sampleID = expected_sample,
    Embedded_RSeQC_sample = embedded_sample,
    Percentile = positions,
    Raw_coverage = values,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 8. Select only final 77 files in metadata order
# ============================================================

selected_idx <- match(
  final_samples,
  file_samples
)

if (anyNA(selected_idx)) {
  stop(
    "Internal file matching error."
  )
}

selected_files <- coverage_files[
  selected_idx
]

selected_file_samples <- file_samples[
  selected_idx
]

stopifnot(
  identical(
    selected_file_samples,
    final_samples
  )
)


# ============================================================
# 9. Parse all 77 profiles
# ============================================================

cat(
  "Reading final 77 GeneBodyCoverage profiles...\n"
)

profile_list <- vector(
  "list",
  length(final_samples)
)

for (i in seq_along(final_samples)) {

  profile_list[[i]] <-
    read_gene_body_file(
      file = selected_files[i],
      expected_sample =
        final_samples[i]
    )
}

profile_raw <- do.call(
  rbind,
  profile_list
)

rownames(profile_raw) <- NULL


# ============================================================
# 10. Verify all final samples
# ============================================================

parsed_samples <- unique(
  profile_raw$sampleID
)

if (length(parsed_samples) != 77) {
  stop(
    "Parsed sample number is not 77."
  )
}

if (!setequal(
  parsed_samples,
  final_samples
)) {
  stop(
    "Parsed coverage sample set differs ",
    "from final metadata."
  )
}

profile_raw$group <- meta$group[
  match(
    profile_raw$sampleID,
    meta$sampleID
  )
]

profile_raw$Stage <- unname(
  stage_labels[
    as.character(
      profile_raw$group
    )
  ]
)


# ============================================================
# 11. RSeQC-style per-sample min-max normalization
#
# normalized = (x - min(x)) / (max(x) - min(x))
#
# This reproduces the transformation used in the
# RSeQC-generated geneBodyCoverage curves.
# ============================================================

normalize_profile <- function(x) {

  x_min <- min(
    x,
    na.rm = TRUE
  )

  x_max <- max(
    x,
    na.rm = TRUE
  )

  denominator <- x_max - x_min

  if (!is.finite(denominator) ||
      denominator <= 0) {
    stop(
      "A sample has zero or invalid ",
      "gene-body coverage range."
    )
  }

  (
    x - x_min
  ) / denominator
}


profile_raw$Relative_coverage <- NA_real_

for (s in final_samples) {

  idx <- which(
    profile_raw$sampleID == s
  )

  profile_raw$Relative_coverage[
    idx
  ] <- normalize_profile(
    profile_raw$Raw_coverage[
      idx
    ]
  )
}


if (anyNA(
  profile_raw$Relative_coverage
)) {
  stop(
    "NA generated during normalization."
  )
}


# ============================================================
# 12. Create normalized 77 x 100 matrix
# ============================================================

normalized_matrix <- matrix(
  NA_real_,
  nrow = length(final_samples),
  ncol = 100,
  dimnames = list(
    final_samples,
    paste0(
      "P",
      1:100
    )
  )
)

raw_matrix <- normalized_matrix

for (s in final_samples) {

  x <- profile_raw[
    profile_raw$sampleID == s,
    ,
    drop = FALSE
  ]

  x <- x[
    order(
      x$Percentile
    ),
    ,
    drop = FALSE
  ]

  normalized_matrix[
    s,
  ] <- x$Relative_coverage

  raw_matrix[
    s,
  ] <- x$Raw_coverage
}


# ============================================================
# 13. Overall median and IQR across 77 samples
# ============================================================

overall_summary <- data.frame(
  Percentile = 1:100,

  Q25 = apply(
    normalized_matrix,
    2,
    quantile,
    probs = 0.25,
    na.rm = TRUE
  ),

  Median = apply(
    normalized_matrix,
    2,
    median,
    na.rm = TRUE
  ),

  Q75 = apply(
    normalized_matrix,
    2,
    quantile,
    probs = 0.75,
    na.rm = TRUE
  )
)


# ============================================================
# 14. Stage median profiles
# ============================================================

stage_median_list <- list()

stage_counter <- 1

for (g in groups) {

  samples_g <- meta$sampleID[
    meta$group == g
  ]

  mat_g <- normalized_matrix[
    samples_g,
    ,
    drop = FALSE
  ]

  stage_median_list[[stage_counter]] <- data.frame(
    Group = g,
    Stage = unname(
      stage_labels[g]
    ),
    N = length(samples_g),
    Percentile = 1:100,
    Median_relative_coverage =
      apply(
        mat_g,
        2,
        median,
        na.rm = TRUE
      ),
    stringsAsFactors = FALSE
  )

  stage_counter <-
    stage_counter + 1
}

stage_median <- do.call(
  rbind,
  stage_median_list
)

rownames(stage_median) <- NULL

stage_median$Group <- factor(
  stage_median$Group,
  levels = groups
)


# ============================================================
# 15. Optional descriptive per-sample 3'/5' coverage metric
#
# Calculated from RAW aggregate coverage:
# mean positions 81-100 / mean positions 1-20
#
# This is descriptive only; no pass/fail cutoff is applied.
# ============================================================

coverage_metrics <- data.frame(
  sampleID = final_samples,
  Group = as.character(
    meta$group[
      match(
        final_samples,
        meta$sampleID
      )
    ]
  ),
  Mean_5prime_raw_coverage = NA_real_,
  Mean_3prime_raw_coverage = NA_real_,
  Ratio_3prime_to_5prime = NA_real_,
  Spearman_to_overall_median =
    NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_along(final_samples)) {

  s <- final_samples[i]

  raw_x <- raw_matrix[
    s,
  ]

  norm_x <- normalized_matrix[
    s,
  ]

  mean5 <- mean(
    raw_x[
      1:20
    ],
    na.rm = TRUE
  )

  mean3 <- mean(
    raw_x[
      81:100
    ],
    na.rm = TRUE
  )

  coverage_metrics$
    Mean_5prime_raw_coverage[i] <-
    mean5

  coverage_metrics$
    Mean_3prime_raw_coverage[i] <-
    mean3

  coverage_metrics$
    Ratio_3prime_to_5prime[i] <-
    ifelse(
      mean5 > 0,
      mean3 / mean5,
      NA_real_
    )

  coverage_metrics$
    Spearman_to_overall_median[i] <-
    suppressWarnings(
      cor(
        norm_x,
        overall_summary$Median,
        method = "spearman"
      )
    )
}


# ============================================================
# 16. Stage-level descriptive metric summary
# ============================================================

stage_metric_summary <- do.call(
  rbind,
  lapply(
    groups,
    function(g) {

      x <- coverage_metrics[
        coverage_metrics$Group == g,
        ,
        drop = FALSE
      ]

      data.frame(
        Group = g,
        Stage = unname(
          stage_labels[g]
        ),
        N = nrow(x),

        Median_3prime_to_5prime_ratio =
          median(
            x$Ratio_3prime_to_5prime,
            na.rm = TRUE
          ),

        Minimum_profile_Spearman_to_overall_median =
          min(
            x$Spearman_to_overall_median,
            na.rm = TRUE
          ),

        Median_profile_Spearman_to_overall_median =
          median(
            x$Spearman_to_overall_median,
            na.rm = TRUE
          ),

        stringsAsFactors = FALSE
      )
    }
  )
)


# ============================================================
# 17. Figure theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 11
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    plot.title = element_text(
      face = "bold",
      size = 11.5
    ),
    axis.title = element_text(
      size = 10.5
    ),
    axis.text = element_text(
      size = 9
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# 18. Panel A
# All 77 final samples + overall median/IQR
# ============================================================

pS1a <- ggplot() +

  geom_ribbon(
    data = overall_summary,
    aes(
      x = Percentile,
      ymin = Q25,
      ymax = Q75
    ),
    fill = "grey75",
    alpha = 0.45
  ) +

  geom_line(
    data = profile_raw,
    aes(
      x = Percentile,
      y = Relative_coverage,
      group = sampleID
    ),
    colour = "grey55",
    alpha = 0.22,
    linewidth = 0.38
  ) +

  geom_line(
    data = overall_summary,
    aes(
      x = Percentile,
      y = Median
    ),
    colour = "black",
    linewidth = 1.1
  ) +

  scale_x_continuous(
    breaks = c(
      1,
      25,
      50,
      75,
      100
    ),
    limits = c(
      1,
      100
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.2
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  labs(
    tag = "A",
    title =
      "Gene-body coverage across all 77 primary-cohort samples",
    x =
      "Gene body percentile (5' to 3')",
    y =
      "Relative coverage"
  ) +

  annotate(
    "text",
    x = 4,
    y = 0.96,
    label =
      "Thin lines: individual samples",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    colour = "grey35"
  ) +

  annotate(
    "text",
    x = 4,
    y = 0.89,
    label =
      "Black line: cohort median; shaded area: IQR",
    hjust = 0,
    vjust = 1,
    size = 3.1,
    colour = "grey35"
  ) +

  theme_pub


# ============================================================
# 19. Panel B
# Stage-level median profiles
# ============================================================

pS1b <- ggplot(
  stage_median,
  aes(
    x = Percentile,
    y = Median_relative_coverage,
    colour = Group,
    group = Group
  )
) +

  geom_line(
    linewidth = 0.95
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

  scale_x_continuous(
    breaks = c(
      1,
      25,
      50,
      75,
      100
    ),
    limits = c(
      1,
      100
    ),
    expand = c(
      0,
      0
    )
  ) +

  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.2
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +

  labs(
    tag = "B",
    title =
      "Stage-level median gene-body coverage profiles",
    x =
      "Gene body percentile (5' to 3')",
    y =
      "Median relative coverage"
  ) +

  theme_pub +

  theme(
    legend.position = "right",
    legend.title = element_text(
      size = 9
    ),
    legend.text = element_text(
      size = 8
    ),
    legend.key.height = unit(
      0.38,
      "cm"
    )
  )


# ============================================================
# 20. Save helpers
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


save_combined <- function(
    pA,
    pB,
    width = 12,
    height = 5.7) {

  draw <- function() {

    grid.newpage()

    pushViewport(
      viewport(
        layout = grid.layout(
          nrow = 1,
          ncol = 2,
          widths = unit(
            c(
              1,
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
  }


  pdf_file <- file.path(
    outdir,
    "Supplementary_Figure_S1_gene_body_coverage.pdf"
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
      "Supplementary_Figure_S1_",
      "gene_body_coverage_600dpi.tiff"
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
# 21. Save figures
# ============================================================

save_panel(
  pS1a,
  "Supplementary_Figure_S1A_all_77_samples",
  width = 6.3,
  height = 5.2
)

save_panel(
  pS1b,
  "Supplementary_Figure_S1B_stage_median_profiles",
  width = 6.8,
  height = 5.2
)

save_combined(
  pS1a,
  pS1b,
  width = 12,
  height = 5.7
)


# ============================================================
# 22. Save source tables
# ============================================================

write.table(
  profile_raw,
  file = file.path(
    outdir,
    "S1_gene_body_coverage_long_final77.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  data.frame(
    sampleID = rownames(
      normalized_matrix
    ),
    normalized_matrix,
    check.names = FALSE
  ),
  file = file.path(
    outdir,
    "S1_gene_body_coverage_normalized_matrix_final77.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  overall_summary,
  file = file.path(
    outdir,
    "S1_overall_gene_body_coverage_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_median,
  file = file.path(
    outdir,
    "S1_stage_median_gene_body_coverage.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  coverage_metrics,
  file = file.path(
    outdir,
    "S1_per_sample_gene_body_coverage_metrics.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_metric_summary,
  file = file.path(
    outdir,
    "S1_stage_gene_body_coverage_metric_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Save R objects and session info
# ============================================================

saveRDS(
  list(
    FigureS1A = pS1a,
    FigureS1B = pS1b
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S1_plot_objects.rds"
  )
)


sink(
  file.path(
    outdir,
    "S1_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 24. Final validation
# ============================================================

n_included <- sum(
  audit$Status ==
    "Included_final77"
)

if (n_included != 77) {
  stop(
    "Final included coverage sample count ",
    "is not 77."
  )
}

if ("R2h_1" %in%
    final_samples) {
  stop(
    "Unexpected: R2h_1 is present in ",
    "the final 77-sample metadata."
  )
}

if (!"R2h_8" %in%
    final_samples) {
  stop(
    "Unexpected: R2h_8 is absent from ",
    "the final 77-sample metadata."
  )
}


# ============================================================
# 25. Final console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S1 COMPLETED\n")
cat("============================================================\n")

cat(
  "GeneBodyCoverage files discovered       : ",
  length(coverage_files),
  "\n",
  sep = ""
)

cat(
  "Final metadata samples                  : ",
  length(final_samples),
  "\n",
  sep = ""
)

cat(
  "Final samples with coverage profiles    : ",
  n_included,
  "\n",
  sep = ""
)

cat(
  "Extra coverage files excluded           : ",
  length(extra_coverage_samples),
  "\n",
  sep = ""
)

if (length(extra_coverage_samples) > 0) {

  cat(
    "Excluded extra sample(s)               : ",
    paste(
      extra_coverage_samples,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "R2h_1 in final cohort                   : ",
  "NO\n"
)

cat(
  "R2h_8 in final cohort                   : ",
  "YES\n"
)

cat(
  "Gene-body percentile positions          : 100\n"
)

cat(
  "Normalization                           : ",
  "RSeQC-style per-sample min-max (0-1)\n"
)

cat(
  "Median 3'/5' raw-coverage ratio         : ",
  sprintf(
    "%.4f",
    median(
      coverage_metrics$
        Ratio_3prime_to_5prime,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Minimum sample-profile Spearman vs median: ",
  sprintf(
    "%.4f",
    min(
      coverage_metrics$
        Spearman_to_overall_median,
      na.rm = TRUE
    )
  ),
  "\n",
  sep = ""
)

cat("\nStage summary:\n")

print(
  stage_metric_summary,
  row.names = FALSE
)

cat("\nOutput directory:\n")
cat(outdir, "\n\n")

cat(
  "Primary output:\n",
  file.path(
    outdir,
    "Supplementary_Figure_S1_gene_body_coverage.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    outdir,
    paste0(
      "Supplementary_Figure_S1_",
      "gene_body_coverage_600dpi.tiff"
    )
  ),
  "\n\n",
  sep = ""
)

cat("PASS\n")
