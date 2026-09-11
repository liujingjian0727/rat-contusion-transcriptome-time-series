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
cat("13B FASTP QC INTEGRATION AND STAGE DIAGNOSTIC\n")
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

metadata_file <-
  "rattus_meta_2026_77samples.tsv"

fastp_file <- file.path(
  DATA_ROOT,
  "fastp",
  "fastp_QC_Cohort_2026.tsv"
)

existing_s2_file <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "Supplementary_Table_S2_2026_primary_77sample_",
  "sequencing_mapping_assignment_QC.tsv"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC/",
  "13B_fastp_QC"
)

dir.create(
  outdir,
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
# 2. Validate inputs
# ============================================================

required_files <- c(
  metadata_file,
  fastp_file,
  existing_s2_file
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

cat("Required files found: PASS\n\n")


# ============================================================
# 3. Metadata
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
    "Expected 77 final samples; observed ",
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

cat(
  "Final metadata samples : ",
  length(final_samples),
  "\n",
  sep = ""
)

cat(
  "R2h_1 present          : ",
  "R2h_1" %in% final_samples,
  "\n",
  sep = ""
)

cat(
  "R2h_8 present          : ",
  "R2h_8" %in% final_samples,
  "\n\n",
  sep = ""
)


# ============================================================
# 4. Read existing sequencing / mapping / assignment S2
# ============================================================

s2 <- read.delim(
  existing_s2_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

if (!"sampleID" %in%
    colnames(s2)) {
  stop(
    "Existing S2 table lacks sampleID."
  )
}

if (nrow(s2) != 77) {
  stop(
    "Expected 77 rows in existing S2; observed ",
    nrow(s2)
  )
}

if (!setequal(
  s2$sampleID,
  final_samples
)) {
  stop(
    "Existing S2 sample set differs from metadata."
  )
}

s2 <- s2[
  match(
    final_samples,
    s2$sampleID
  ),
  ,
  drop = FALSE
]

cat(
  "Existing sequencing/mapping/assignment S2 : PASS\n\n"
)


# ============================================================
# 5. Read consolidated fastp QC
# ============================================================

fastp <- read.delim(
  fastp_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  stringsAsFactors = FALSE
)

cat(
  "fastp QC dimensions : ",
  nrow(fastp),
  " rows x ",
  ncol(fastp),
  " columns\n",
  sep = ""
)

required_fastp_cols <- c(
  "SampleID",
  "Raw_reads",
  "Raw_read_pairs",
  "Clean_reads",
  "Clean_read_pairs",
  "Read_retention_percent",
  "Raw_bases",
  "Clean_bases",
  "Base_retention_percent",
  "Q20_clean_percent",
  "Q30_clean_percent",
  "GC_clean_percent",
  "Mean_length_R1_after_bp",
  "Mean_length_R2_after_bp",
  "Duplication_rate_percent",
  "Adapter_trimmed_reads_percent",
  "Passed_filter_reads",
  "Failed_low_quality_reads",
  "Failed_too_many_N_reads",
  "Failed_too_short_reads"
)

missing_fastp_cols <- setdiff(
  required_fastp_cols,
  colnames(fastp)
)

if (length(missing_fastp_cols) > 0) {
  stop(
    "Missing fastp column(s): ",
    paste(
      missing_fastp_cols,
      collapse = ", "
    )
  )
}

if (anyDuplicated(
  fastp$SampleID
)) {
  stop(
    "Duplicated SampleID in fastp table."
  )
}


# ============================================================
# 6. Final-77 whitelist audit
# ============================================================

fastp_samples <- fastp$SampleID

fastp_extra <- setdiff(
  fastp_samples,
  final_samples
)

fastp_missing <- setdiff(
  final_samples,
  fastp_samples
)

cat(
  "fastp samples discovered : ",
  length(fastp_samples),
  "\n",
  sep = ""
)

cat(
  "Extra fastp samples      : ",
  length(fastp_extra),
  "\n",
  sep = ""
)

if (length(fastp_extra) > 0) {
  cat(
    "  ",
    paste(
      fastp_extra,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Final samples missing    : ",
  length(fastp_missing),
  "\n\n",
  sep = ""
)

if (length(fastp_missing) > 0) {
  stop(
    "One or more final samples missing from fastp QC."
  )
}


fastp77 <- fastp[
  match(
    final_samples,
    fastp$SampleID
  ),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    fastp77$SampleID,
    final_samples
  )
)


# ============================================================
# 7. Merge fastp metrics into Supplementary Table S2
# ============================================================

fastp_for_merge <- fastp77[
  ,
  required_fastp_cols,
  drop = FALSE
]

colnames(fastp_for_merge)[
  colnames(fastp_for_merge) == "SampleID"
] <- "sampleID"


# Avoid duplicate column names already present in S2
duplicate_metric_names <- intersect(
  setdiff(
    colnames(fastp_for_merge),
    "sampleID"
  ),
  colnames(s2)
)

if (length(duplicate_metric_names) > 0) {

  cat(
    "Existing duplicate metric names in S2:\n  ",
    paste(
      duplicate_metric_names,
      collapse = "\n  "
    ),
    "\n",
    sep = ""
  )

  s2 <- s2[
    ,
    !colnames(s2) %in%
      duplicate_metric_names,
    drop = FALSE
  ]
}


full_qc <- merge(
  s2,
  fastp_for_merge,
  by = "sampleID",
  all.x = TRUE,
  sort = FALSE
)

full_qc <- full_qc[
  match(
    final_samples,
    full_qc$sampleID
  ),
  ,
  drop = FALSE
]

full_qc$group <- meta$group[
  match(
    full_qc$sampleID,
    meta$sampleID
  )
]

full_qc$Stage <- unname(
  stage_labels[
    as.character(
      full_qc$group
    )
  ]
)

full_qc$group <- factor(
  full_qc$group,
  levels = groups
)


if (nrow(full_qc) != 77) {
  stop(
    "Merged table does not contain 77 samples."
  )
}

cat(
  "Integrated full sequencing QC table : PASS\n\n"
)


# ============================================================
# 8. FASTQ stats vs fastp raw-read consistency
# ============================================================

raw_consistency <- data.frame(
  sampleID =
    full_qc$sampleID,

  FASTQ_stats_raw_reads =
    full_qc$Raw_reads,

  fastp_raw_reads =
    full_qc$Raw_reads,
  stringsAsFactors = FALSE
)


# Existing 13A S2 Raw_reads may have been overwritten above.
# Recover 13A raw-read count from original table.

original_s2 <- read.delim(
  existing_s2_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

original_s2 <- original_s2[
  match(
    final_samples,
    original_s2$sampleID
  ),
  ,
  drop = FALSE
]

if ("Raw_reads" %in%
    colnames(original_s2)) {

  raw_consistency$
    FASTQ_stats_raw_reads <-
    original_s2$Raw_reads

  raw_consistency$
    fastp_raw_reads <-
    fastp77$Raw_reads

  raw_consistency$
    Difference <-
    raw_consistency$
      fastp_raw_reads -
    raw_consistency$
      FASTQ_stats_raw_reads

  raw_consistency$
    Exact_match <-
    raw_consistency$Difference == 0

  exact_matches <- sum(
    raw_consistency$Exact_match
  )

  max_abs_difference <- max(
    abs(
      raw_consistency$Difference
    )
  )

} else {

  exact_matches <- NA_integer_
  max_abs_difference <- NA_real_
}


# ============================================================
# 9. Cohort-level summary
# ============================================================

metric_map <- data.frame(
  Variable = c(
    "Q20_clean_percent",
    "Q30_clean_percent",
    "GC_clean_percent",
    "Duplication_rate_percent",
    "Read_retention_percent",
    "Base_retention_percent",
    "Adapter_trimmed_reads_percent"
  ),

  Label = c(
    "Q20 clean bases (%)",
    "Q30 clean bases (%)",
    "GC clean bases (%)",
    "Duplication rate (%)",
    "Read retention (%)",
    "Base retention (%)",
    "Adapter-trimmed reads (%)"
  ),

  stringsAsFactors = FALSE
)


cohort_summary_list <- vector(
  "list",
  nrow(metric_map)
)

for (i in seq_len(
  nrow(metric_map)
)) {

  variable <- metric_map$Variable[i]

  x <- full_qc[
    ,
    variable
  ]

  cohort_summary_list[[i]] <- data.frame(
    Metric =
      metric_map$Label[i],

    N = sum(
      is.finite(x)
    ),

    Minimum = min(
      x,
      na.rm = TRUE
    ),

    Median = median(
      x,
      na.rm = TRUE
    ),

    Mean = mean(
      x,
      na.rm = TRUE
    ),

    Maximum = max(
      x,
      na.rm = TRUE
    ),

    stringsAsFactors = FALSE
  )
}

cohort_summary <- do.call(
  rbind,
  cohort_summary_list
)

rownames(cohort_summary) <- NULL


# ============================================================
# 10. Stage-level summaries
# ============================================================

stage_summary_list <- list()
counter <- 1

for (g in groups) {

  tmp <- full_qc[
    full_qc$group == g,
    ,
    drop = FALSE
  ]

  for (i in seq_len(
    nrow(metric_map)
  )) {

    variable <- metric_map$Variable[i]

    x <- tmp[
      ,
      variable
    ]

    stage_summary_list[[counter]] <- data.frame(
      Group = g,
      Stage = unname(
        stage_labels[g]
      ),
      N = nrow(tmp),
      Metric =
        metric_map$Label[i],
      Variable =
        variable,
      Minimum = min(
        x,
        na.rm = TRUE
      ),
      Median = median(
        x,
        na.rm = TRUE
      ),
      Mean = mean(
        x,
        na.rm = TRUE
      ),
      Maximum = max(
        x,
        na.rm = TRUE
      ),
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

stage_summary <- do.call(
  rbind,
  stage_summary_list
)

rownames(stage_summary) <- NULL


# ============================================================
# 11. Kruskal-Wallis stage tests
# ============================================================

stage_tests_list <- vector(
  "list",
  nrow(metric_map)
)

for (i in seq_len(
  nrow(metric_map)
)) {

  variable <- metric_map$Variable[i]

  test_data <- data.frame(
    y = full_qc[
      ,
      variable
    ],
    group = full_qc$group
  )

  test_data <- test_data[
    is.finite(test_data$y) &
      !is.na(test_data$group),
    ,
    drop = FALSE
  ]

  kw <- kruskal.test(
    y ~ group,
    data = test_data
  )

  H <- unname(
    kw$statistic
  )

  n <- nrow(test_data)

  k <- length(
    unique(
      test_data$group
    )
  )

  epsilon_squared <- max(
    0,
    (
      H -
      k +
      1
    ) /
      (
        n -
        k
      )
  )

  stage_tests_list[[i]] <- data.frame(
    Metric =
      metric_map$Label[i],
    Variable =
      variable,
    N = n,
    Kruskal_Wallis_H = H,
    df =
      unname(
        kw$parameter
      ),
    P_value =
      kw$p.value,
    Epsilon_squared =
      epsilon_squared,
    stringsAsFactors = FALSE
  )
}

stage_tests <- do.call(
  rbind,
  stage_tests_list
)

stage_tests$BH_FDR <- p.adjust(
  stage_tests$P_value,
  method = "BH"
)


# ============================================================
# 12. fastp metrics vs downstream metrics
# ============================================================

fastp_predictors <- c(
  "Q30_clean_percent",
  "GC_clean_percent",
  "Duplication_rate_percent",
  "Read_retention_percent",
  "Base_retention_percent"
)

downstream_metrics <- c(
  "HISAT2_overall_alignment_rate_percent",
  "featureCounts_assignment_rate_percent"
)

cor_list <- list()
counter <- 1

for (x_name in fastp_predictors) {

  for (y_name in downstream_metrics) {

    x <- full_qc[
      ,
      x_name
    ]

    y <- full_qc[
      ,
      y_name
    ]

    keep <- is.finite(x) &
      is.finite(y)

    ct <- suppressWarnings(
      cor.test(
        x[keep],
        y[keep],
        method = "spearman",
        exact = FALSE
      )
    )

    cor_list[[counter]] <- data.frame(
      Fastp_metric = x_name,
      Downstream_metric = y_name,
      N = sum(keep),
      Spearman_rho =
        unname(
          ct$estimate
        ),
      P_value =
        ct$p.value,
      stringsAsFactors = FALSE
    )

    counter <- counter + 1
  }
}

correlation_summary <- do.call(
  rbind,
  cor_list
)

correlation_summary$BH_FDR <- p.adjust(
  correlation_summary$P_value,
  method = "BH"
)


# ============================================================
# 13. Extreme-value descriptive tables
#
# Not exclusion rules.
# ============================================================

lowest_q30 <- full_qc[
  order(
    full_qc$Q30_clean_percent
  ),
  c(
    "sampleID",
    "group",
    "Q30_clean_percent",
    "GC_clean_percent",
    "Duplication_rate_percent",
    "Read_retention_percent",
    "HISAT2_overall_alignment_rate_percent",
    "featureCounts_assignment_rate_percent"
  ),
  drop = FALSE
][1:10, ]


highest_duplication <- full_qc[
  order(
    full_qc$Duplication_rate_percent,
    decreasing = TRUE
  ),
  c(
    "sampleID",
    "group",
    "Q30_clean_percent",
    "GC_clean_percent",
    "Duplication_rate_percent",
    "Read_retention_percent",
    "HISAT2_overall_alignment_rate_percent",
    "featureCounts_assignment_rate_percent"
  ),
  drop = FALSE
][1:10, ]


lowest_retention <- full_qc[
  order(
    full_qc$Read_retention_percent
  ),
  c(
    "sampleID",
    "group",
    "Q30_clean_percent",
    "GC_clean_percent",
    "Duplication_rate_percent",
    "Read_retention_percent",
    "HISAT2_overall_alignment_rate_percent",
    "featureCounts_assignment_rate_percent"
  ),
  drop = FALSE
][1:10, ]


# ============================================================
# 14. Plot theme
# ============================================================

theme_pub <- theme_classic(
  base_size = 10
) +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),

    plot.title = element_text(
      face = "bold",
      size = 10.8
    ),

    plot.subtitle = element_text(
      size = 8.3,
      colour = "grey35"
    ),

    axis.title = element_text(
      size = 9.5
    ),

    axis.text = element_text(
      size = 7.8
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )


# ============================================================
# 15. Generic stage-panel function
# ============================================================

make_stage_panel <- function(
    variable,
    y_label,
    title,
    subtitle,
    tag) {

  plot_data <- data.frame(
    group = full_qc$group,
    y = full_qc[
      ,
      variable
    ],
    stringsAsFactors = FALSE
  )

  plot_data$group <- factor(
    plot_data$group,
    levels = groups
  )

  cohort_median <- median(
    plot_data$y,
    na.rm = TRUE
  )

  ggplot(
    plot_data,
    aes(
      x = group,
      y = y
    )
  ) +

    geom_hline(
      yintercept =
        cohort_median,
      linetype = "dashed",
      linewidth = 0.55,
      colour = "grey45"
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
      width = 0.16,
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
      tag = tag,
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = y_label
    ) +

    theme_pub +

    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      )
    )
}


# ============================================================
# 16. Four-panel diagnostic
# ============================================================

pA <- make_stage_panel(
  variable =
    "Q30_clean_percent",
  y_label =
    "Q30 clean bases (%)",
  title =
    "Base-call quality",
  subtitle =
    "Q30 proportion after fastp processing",
  tag = "A"
)


pB <- make_stage_panel(
  variable =
    "GC_clean_percent",
  y_label =
    "GC content (%)",
  title =
    "GC content",
  subtitle =
    "GC proportion after fastp processing",
  tag = "B"
)


pC <- make_stage_panel(
  variable =
    "Duplication_rate_percent",
  y_label =
    "Duplication rate (%)",
  title =
    "Estimated read duplication",
  subtitle =
    "Duplication rate reported by fastp",
  tag = "C"
)


pD <- make_stage_panel(
  variable =
    "Read_retention_percent",
  y_label =
    "Retained reads (%)",
  title =
    "Read retention",
  subtitle =
    "Clean reads / raw reads",
  tag = "D"
)


# ============================================================
# 17. Save figure
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


draw_combined <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 2
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
  "13B_fastp_stage_diagnostic.pdf"
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
  "13B_fastp_stage_diagnostic_600dpi.tiff"
)

open_tiff(
  tiff_file,
  13,
  10
)

draw_combined()
dev.off()


# ============================================================
# 18. Save final Supplementary Table S2
# ============================================================

final_s2_file <- file.path(
  outdir,
  paste0(
    "Supplementary_Table_S2_",
    "2026_primary_77sample_",
    "full_sequencing_QC.tsv"
  )
)

write.table(
  full_qc,
  file = final_s2_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 19. Save diagnostics
# ============================================================

write.table(
  cohort_summary,
  file = file.path(
    outdir,
    "13B_fastp_cohort_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  stage_summary,
  file = file.path(
    outdir,
    "13B_fastp_stage_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  stage_tests,
  file = file.path(
    outdir,
    "13B_fastp_stage_Kruskal_Wallis_tests.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  correlation_summary,
  file = file.path(
    outdir,
    "13B_fastp_vs_mapping_assignment_correlations.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lowest_q30,
  file = file.path(
    outdir,
    "13B_lowest_Q30_samples.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  highest_duplication,
  file = file.path(
    outdir,
    "13B_highest_duplication_samples.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  lowest_retention,
  file = file.path(
    outdir,
    "13B_lowest_retention_samples.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  raw_consistency,
  file = file.path(
    outdir,
    "13B_FASTQ_stats_vs_fastp_raw_read_consistency.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 20. Sample audit
# ============================================================

audit_samples <- union(
  final_samples,
  fastp_samples
)

audit <- data.frame(
  sampleID = audit_samples,
  In_final_77_metadata =
    audit_samples %in%
      final_samples,
  In_fastp_QC =
    audit_samples %in%
      fastp_samples,
  stringsAsFactors = FALSE
)

audit$Status <- ifelse(
  audit$In_final_77_metadata &
    audit$In_fastp_QC,
  "Included_final77",
  ifelse(
    !audit$In_final_77_metadata &
      audit$In_fastp_QC,
    "Excluded_not_in_final_metadata",
    "Missing_fastp_QC"
  )
)

write.table(
  audit,
  file = file.path(
    outdir,
    "13B_fastp_sample_whitelist_audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 21. Session info
# ============================================================

saveRDS(
  list(
    A_Q30 = pA,
    B_GC = pB,
    C_duplication = pC,
    D_retention = pD
  ),
  file = file.path(
    outdir,
    "13B_fastp_plot_objects.rds"
  )
)

sink(
  file.path(
    outdir,
    "13B_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 22. Final validation
# ============================================================

if ("R2h_1" %in%
    full_qc$sampleID) {
  stop(
    "R2h_1 unexpectedly present."
  )
}

if (!("R2h_8" %in%
      full_qc$sampleID)) {
  stop(
    "R2h_8 unexpectedly absent."
  )
}

if (nrow(full_qc) != 77) {
  stop(
    "Final S2 table is not 77 samples."
  )
}


# ============================================================
# 23. Console report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FASTP COHORT SUMMARY\n")
cat("============================================================\n\n")

print(
  cohort_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("STAGE-LEVEL KRUSKAL-WALLIS TESTS\n")
cat("============================================================\n\n")

print(
  stage_tests,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("FASTP VS MAPPING / ASSIGNMENT\n")
cat("============================================================\n\n")

print(
  correlation_summary,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("LOWEST Q30 SAMPLES\n")
cat("============================================================\n\n")

print(
  lowest_q30,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("HIGHEST DUPLICATION SAMPLES\n")
cat("============================================================\n\n")

print(
  highest_duplication,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("LOWEST READ-RETENTION SAMPLES\n")
cat("============================================================\n\n")

print(
  lowest_retention,
  row.names = FALSE,
  digits = 5
)


cat("\n")
cat("============================================================\n")
cat("WHITELIST / CONSISTENCY AUDIT\n")
cat("============================================================\n\n")

cat(
  "fastp table samples               : ",
  length(fastp_samples),
  "\n",
  sep = ""
)

cat(
  "Extra fastp sample(s)             : ",
  ifelse(
    length(fastp_extra) == 0,
    "None",
    paste(
      fastp_extra,
      collapse = ", "
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Final samples missing             : ",
  length(fastp_missing),
  "\n",
  sep = ""
)

cat(
  "R2h_1 excluded                    : ",
  !("R2h_1" %in%
      full_qc$sampleID),
  "\n",
  sep = ""
)

cat(
  "R2h_8 retained                    : ",
  "R2h_8" %in%
    full_qc$sampleID,
  "\n",
  sep = ""
)

cat(
  "Exact raw-read-count matches      : ",
  exact_matches,
  "/77\n",
  sep = ""
)

cat(
  "Maximum raw-read count difference : ",
  max_abs_difference,
  "\n",
  sep = ""
)


cat("\n")
cat("============================================================\n")
cat("OUTPUTS\n")
cat("============================================================\n\n")

cat(
  "Final Supplementary Table S2:\n",
  final_s2_file,
  "\n\n",
  sep = ""
)

cat(
  "Diagnostic figure:\n",
  pdf_file,
  "\n",
  tiff_file,
  "\n\n",
  sep = ""
)

cat("PASS\n")
