#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "R package 'DESeq2' is not installed.\n",
      "Please install DESeq2 before running this script."
    )
  }
  library(DESeq2)
})

cat("\n")
cat("============================================================\n")
cat("02 COUNT FILTERING AND DESeq2 VST\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")

# ============================================================
# 0. Fixed analysis parameters
# ============================================================

MIN_COUNT <- 10
MIN_SAMPLES <- 6

cat("Fixed filtering criterion:\n")
cat("  Raw count >=", MIN_COUNT,
    "in >=", MIN_SAMPLES, "samples\n\n")

# ============================================================
# 1. Input files
# ============================================================

cohorts <- list(

  cohort2026 = list(
    name = "2026_primary",
    short = "2026",
    counts = "rattus_counts_2026_primary_77samples.txt",
    tpm = "rattus_77sample.TMM.TPM.tsv",
    meta = "rattus_meta_2026_77samples.tsv",
    expected_n = 77,
    group_levels = c(
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
  ),

  cohort2021 = list(
    name = "2021_historical_reference",
    short = "2021",
    counts = "rattus_2021.counts.tsv",
    tpm = "rattus_2021.TMM.TPM.tsv",
    meta = "rattus_meta_2021.tsv",
    expected_n = 60,
    group_levels = c(
      "Baseline",
      "R4h",
      "R8h",
      "R12h",
      "R16h",
      "R20h",
      "R24h",
      "R48h"
    )
  )
)

outdir <- "02_count_filtering_and_VST_result"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# 2. Helper functions
# ============================================================

read_expression_matrix <- function(file) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  if (ncol(x) < 2) {
    stop("Invalid expression matrix: ", file)
  }

  colnames(x)[1] <- "GeneID"

  if (anyDuplicated(x$GeneID)) {
    stop("Duplicated GeneIDs detected in: ", file)
  }

  if (any(is.na(x$GeneID) | trimws(x$GeneID) == "")) {
    stop("Missing GeneIDs detected in: ", file)
  }

  return(x)
}


read_metadata <- function(file) {

  meta <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  required <- c(
    "sampleID",
    "group",
    "time",
    "condition"
  )

  missing_columns <- setdiff(required, colnames(meta))

  if (length(missing_columns) > 0) {
    stop(
      "Missing metadata columns in ",
      file,
      ": ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (anyDuplicated(meta$sampleID)) {
    stop("Duplicated sampleIDs detected in: ", file)
  }

  return(meta)
}


write_matrix_with_geneid <- function(mat, file) {

  out <- data.frame(
    GeneID = rownames(mat),
    mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.table(
    out,
    file = file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


# ============================================================
# 3. Main processing function
# ============================================================

process_cohort <- function(info) {

  cohort_name <- info$name
  short_name <- info$short

  cohort_outdir <- file.path(
    outdir,
    cohort_name
  )

  dir.create(
    cohort_outdir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  cat("------------------------------------------------------------\n")
  cat("Processing:", cohort_name, "\n")
  cat("------------------------------------------------------------\n\n")

  # ----------------------------------------------------------
  # 3.1 Read data
  # ----------------------------------------------------------

  cat("Reading raw counts...\n")
  counts_df <- read_expression_matrix(info$counts)

  cat("Reading TMM.TPM...\n")
  tpm_df <- read_expression_matrix(info$tpm)

  cat("Reading metadata...\n")
  meta <- read_metadata(info$meta)

  count_samples <- colnames(counts_df)[-1]
  tpm_samples <- colnames(tpm_df)[-1]

  # ----------------------------------------------------------
  # 3.2 Confirm sample consistency
  # ----------------------------------------------------------

  if (!setequal(count_samples, tpm_samples)) {
    stop(
      cohort_name,
      ": counts and TPM sample sets are different."
    )
  }

  if (!setequal(count_samples, meta$sampleID)) {
    stop(
      cohort_name,
      ": counts and metadata sample sets are different."
    )
  }

  if (length(count_samples) != info$expected_n) {
    stop(
      cohort_name,
      ": unexpected sample number."
    )
  }

  # ----------------------------------------------------------
  # 3.3 Reorder metadata to exactly match counts
  # ----------------------------------------------------------

  meta <- meta[
    match(count_samples, meta$sampleID),
    ,
    drop = FALSE
  ]

  if (!identical(meta$sampleID, count_samples)) {
    stop(
      cohort_name,
      ": failed to reorder metadata."
    )
  }

  # ----------------------------------------------------------
  # 3.4 Set biological group order
  # ----------------------------------------------------------

  unexpected_groups <- setdiff(
    unique(meta$group),
    info$group_levels
  )

  if (length(unexpected_groups) > 0) {
    stop(
      cohort_name,
      ": unexpected groups detected: ",
      paste(unexpected_groups, collapse = ", ")
    )
  }

  missing_groups <- setdiff(
    info$group_levels,
    unique(meta$group)
  )

  if (length(missing_groups) > 0) {
    stop(
      cohort_name,
      ": expected groups missing: ",
      paste(missing_groups, collapse = ", ")
    )
  }

  meta$group <- factor(
    meta$group,
    levels = info$group_levels
  )

  meta$condition <- factor(
    meta$condition,
    levels = c("Control", "Injury")
  )

  rownames(meta) <- meta$sampleID

  # ----------------------------------------------------------
  # 3.5 Convert counts to numeric matrix
  # ----------------------------------------------------------

  count_mat <- as.matrix(
    counts_df[, -1, drop = FALSE]
  )

  suppressWarnings(
    storage.mode(count_mat) <- "numeric"
  )

  rownames(count_mat) <- counts_df$GeneID

  if (anyNA(count_mat)) {
    stop(cohort_name, ": NA values detected in count matrix.")
  }

  if (any(!is.finite(count_mat))) {
    stop(cohort_name, ": non-finite counts detected.")
  }

  if (any(count_mat < 0)) {
    stop(cohort_name, ": negative counts detected.")
  }

  if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
    stop(cohort_name, ": noninteger counts detected.")
  }

  count_mat <- round(count_mat)
  storage.mode(count_mat) <- "integer"

  # ----------------------------------------------------------
  # 3.6 Explicitly align TPM to count GeneID order
  # ----------------------------------------------------------

  if (!setequal(
    rownames(count_mat),
    tpm_df$GeneID
  )) {
    stop(
      cohort_name,
      ": counts and TPM GeneID sets differ."
    )
  }

  tpm_match <- match(
    rownames(count_mat),
    tpm_df$GeneID
  )

  if (anyNA(tpm_match)) {
    stop(
      cohort_name,
      ": GeneID matching between counts and TPM failed."
    )
  }

  tpm_aligned_df <- tpm_df[
    tpm_match,
    ,
    drop = FALSE
  ]

  if (!identical(
    rownames(count_mat),
    tpm_aligned_df$GeneID
  )) {
    stop(
      cohort_name,
      ": TPM GeneID alignment failed."
    )
  }

  # Reorder TPM sample columns to count sample order
  tpm_aligned_df <- tpm_aligned_df[
    ,
    c("GeneID", count_samples),
    drop = FALSE
  ]

  cat("\nInput dimensions:\n")
  cat(
    "Counts :",
    nrow(count_mat),
    "genes x",
    ncol(count_mat),
    "samples\n"
  )

  cat(
    "TPM    :",
    nrow(tpm_aligned_df),
    "genes x",
    ncol(tpm_aligned_df) - 1,
    "samples\n"
  )

  cat(
    "Meta   :",
    nrow(meta),
    "samples\n\n"
  )

  cat("Counts/TPM GeneIDs explicitly aligned by GeneID: PASS\n")
  cat("Counts/TPM sample order explicitly aligned: PASS\n\n")

  # ----------------------------------------------------------
  # 3.7 Gene-level filtering statistics
  # ----------------------------------------------------------

  cat("Calculating gene-level filtering metrics...\n")

  samples_ge_min_count <- rowSums(
    count_mat >= MIN_COUNT
  )

  total_count <- rowSums(count_mat)

  mean_count <- rowMeans(count_mat)

  median_count <- apply(
    count_mat,
    1,
    median
  )

  max_count <- apply(
    count_mat,
    1,
    max
  )

  keep <- samples_ge_min_count >= MIN_SAMPLES

  filtering_metrics <- data.frame(
    GeneID = rownames(count_mat),
    Samples_with_count_GE_10 = samples_ge_min_count,
    Total_raw_count = total_count,
    Mean_raw_count = mean_count,
    Median_raw_count = median_count,
    Maximum_raw_count = max_count,
    Keep = keep,
    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # 3.8 Filtering summary
  # ----------------------------------------------------------

  genes_before <- nrow(count_mat)
  genes_after <- sum(keep)
  genes_removed <- genes_before - genes_after

  retained_percent <- 100 * genes_after / genes_before

  zero_genes <- sum(total_count == 0)

  cat("\nGene filtering:\n")
  cat(
    "Genes before filtering :",
    genes_before, "\n"
  )
  cat(
    "Genes retained         :",
    genes_after, "\n"
  )
  cat(
    "Genes removed          :",
    genes_removed, "\n"
  )
  cat(
    "Genes retained (%)     :",
    sprintf("%.2f", retained_percent), "\n"
  )
  cat(
    "All-zero genes         :",
    zero_genes, "\n"
  )
  cat(
    "Filtering criterion    : count >=",
    MIN_COUNT,
    "in >=",
    MIN_SAMPLES,
    "samples\n\n"
  )

  if (genes_after < 1000) {
    stop(
      cohort_name,
      ": fewer than 1,000 genes remained after filtering."
    )
  }

  # ----------------------------------------------------------
  # 3.9 Apply gene filter
  # ----------------------------------------------------------

  filtered_counts <- count_mat[
    keep,
    ,
    drop = FALSE
  ]

  filtered_geneids <- rownames(filtered_counts)

  filtered_tpm <- tpm_aligned_df[
    match(filtered_geneids, tpm_aligned_df$GeneID),
    ,
    drop = FALSE
  ]

  if (!identical(
    filtered_geneids,
    filtered_tpm$GeneID
  )) {
    stop(
      cohort_name,
      ": filtered TPM alignment failed."
    )
  }

  # ----------------------------------------------------------
  # 3.10 Create DESeq2 object
  # ----------------------------------------------------------

  cat("Creating DESeq2 object...\n")

  dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = meta,
    design = ~ group
  )

  # ----------------------------------------------------------
  # 3.11 Estimate DESeq2 size factors
  # ----------------------------------------------------------

  cat("Estimating DESeq2 size factors...\n")

  dds <- estimateSizeFactors(dds)

  size_factor_table <- data.frame(
    sampleID = colnames(dds),
    group = as.character(colData(dds)$group),
    condition = as.character(colData(dds)$condition),
    size_factor = sizeFactors(dds),
    stringsAsFactors = FALSE
  )

  if (anyNA(size_factor_table$size_factor) ||
      any(!is.finite(size_factor_table$size_factor)) ||
      any(size_factor_table$size_factor <= 0)) {

    stop(
      cohort_name,
      ": invalid DESeq2 size factors detected."
    )
  }

  # ----------------------------------------------------------
  # 3.12 Variance stabilizing transformation
  # ----------------------------------------------------------

  cat("Running DESeq2 variance-stabilizing transformation...\n")

  vsd <- varianceStabilizingTransformation(
    dds,
    blind = TRUE
  )

  vst_mat <- assay(vsd)

  if (anyNA(vst_mat) ||
      any(!is.finite(vst_mat))) {

    stop(
      cohort_name,
      ": invalid values detected in VST matrix."
    )
  }

  cat(
    "VST matrix :",
    nrow(vst_mat),
    "genes x",
    ncol(vst_mat),
    "samples\n\n"
  )

  # ----------------------------------------------------------
  # 3.13 Save filtering metrics
  # ----------------------------------------------------------

  write.table(
    filtering_metrics,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_gene_filtering_metrics.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # 3.14 Save retained gene list
  # ----------------------------------------------------------

  retained_gene_table <- data.frame(
    GeneID = filtered_geneids,
    stringsAsFactors = FALSE
  )

  write.table(
    retained_gene_table,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_filtered_gene_list.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # 3.15 Save filtered raw count matrix
  # ----------------------------------------------------------

  write_matrix_with_geneid(
    filtered_counts,
    file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_filtered_counts.tsv"
      )
    )
  )

  # ----------------------------------------------------------
  # 3.16 Save filtered TPM matrix
  # ----------------------------------------------------------

  write.table(
    filtered_tpm,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_filtered_TMM_TPM.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # 3.17 Save VST matrix
  # ----------------------------------------------------------

  write_matrix_with_geneid(
    vst_mat,
    file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_VST_matrix.tsv"
      )
    )
  )

  # ----------------------------------------------------------
  # 3.18 Save size factors
  # ----------------------------------------------------------

  write.table(
    size_factor_table,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_DESeq2_size_factors.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # 3.19 Save reordered metadata
  # ----------------------------------------------------------

  meta_output <- as.data.frame(meta)

  meta_output$group <- as.character(meta_output$group)
  meta_output$condition <- as.character(meta_output$condition)

  write.table(
    meta_output,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_metadata_aligned.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # 3.20 Save DESeq2/VST objects for exact reproducibility
  # ----------------------------------------------------------

  saveRDS(
    dds,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_DESeq2_sizefactor_object.rds"
      )
    )
  )

  saveRDS(
    vsd,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_VST_object.rds"
      )
    )
  )

  # ----------------------------------------------------------
  # 3.21 Cohort summary
  # ----------------------------------------------------------

  summary_table <- data.frame(
    Cohort = cohort_name,
    Samples = ncol(count_mat),
    Genes_before_filtering = genes_before,
    Genes_after_filtering = genes_after,
    Genes_removed = genes_removed,
    Percent_genes_retained = round(
      retained_percent,
      3
    ),
    All_zero_genes = zero_genes,
    Minimum_count_threshold = MIN_COUNT,
    Minimum_sample_threshold = MIN_SAMPLES,
    VST_genes = nrow(vst_mat),
    VST_samples = ncol(vst_mat),
    Status = "PASS",
    stringsAsFactors = FALSE
  )

  write.table(
    summary_table,
    file = file.path(
      cohort_outdir,
      paste0(
        "02_",
        short_name,
        "_processing_summary.tsv"
      )
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  cat("Final cohort status: PASS\n")
  cat("------------------------------------------------------------\n\n")

  return(
    list(
      summary = summary_table,
      genes = filtered_geneids,
      samples = colnames(filtered_counts)
    )
  )
}

# ============================================================
# 4. Check all required files
# ============================================================

required_files <- c(
  cohorts$cohort2026$counts,
  cohorts$cohort2026$tpm,
  cohorts$cohort2026$meta,
  cohorts$cohort2021$counts,
  cohorts$cohort2021$tpm,
  cohorts$cohort2021$meta
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  cat("ERROR: Missing input files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("All six required input files were found.\n\n")

# ============================================================
# 5. Process 2026 primary cohort
# ============================================================

result2026 <- process_cohort(
  cohorts$cohort2026
)

# ============================================================
# 6. Process 2021 historical reference cohort
# ============================================================

result2021 <- process_cohort(
  cohorts$cohort2021
)

# ============================================================
# 7. Cross-cohort filtered GeneID comparison
# ============================================================

cat("============================================================\n")
cat("CROSS-COHORT FILTERED GENE COMPARISON\n")
cat("============================================================\n\n")

common_genes <- intersect(
  result2026$genes,
  result2021$genes
)

only_2026 <- setdiff(
  result2026$genes,
  result2021$genes
)

only_2021 <- setdiff(
  result2021$genes,
  result2026$genes
)

cat(
  "2026 retained genes      :",
  length(result2026$genes),
  "\n"
)

cat(
  "2021 retained genes      :",
  length(result2021$genes),
  "\n"
)

cat(
  "Common retained genes    :",
  length(common_genes),
  "\n"
)

cat(
  "2026-only retained genes :",
  length(only_2026),
  "\n"
)

cat(
  "2021-only retained genes :",
  length(only_2021),
  "\n\n"
)

# Preserve 2026 order for common genes
common_genes_ordered <- result2026$genes[
  result2026$genes %in% common_genes
]

write.table(
  data.frame(
    GeneID = common_genes_ordered,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "02_cross_cohort_common_filtered_genes.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = only_2026,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "02_genes_retained_only_in_2026.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  data.frame(
    GeneID = only_2021,
    stringsAsFactors = FALSE
  ),
  file = file.path(
    outdir,
    "02_genes_retained_only_in_2021.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# 8. Combined filtering summary
# ============================================================

combined_summary <- rbind(
  result2026$summary,
  result2021$summary
)

write.table(
  combined_summary,
  file = file.path(
    outdir,
    "02_count_filtering_and_VST_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cross_summary <- data.frame(
  Metric = c(
    "2026_retained_genes",
    "2021_retained_genes",
    "Common_retained_genes",
    "2026_only_retained_genes",
    "2021_only_retained_genes"
  ),
  Value = c(
    length(result2026$genes),
    length(result2021$genes),
    length(common_genes),
    length(only_2026),
    length(only_2021)
  ),
  stringsAsFactors = FALSE
)

write.table(
  cross_summary,
  file = file.path(
    outdir,
    "02_cross_cohort_filtered_gene_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# ============================================================
# 9. Record R session information
# ============================================================

sink(
  file.path(
    outdir,
    "02_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()

# ============================================================
# 10. Final report
# ============================================================

cat("============================================================\n")
cat("02 COUNT FILTERING AND VST COMPLETED\n")
cat("============================================================\n")

cat(
  "2026 primary cohort       : PASS\n"
)

cat(
  "2021 historical reference: PASS\n"
)

cat(
  "Filtering criterion       : count >=",
  MIN_COUNT,
  "in >=",
  MIN_SAMPLES,
  "samples\n"
)

cat(
  "2026 retained genes       :",
  length(result2026$genes),
  "\n"
)

cat(
  "2021 retained genes       :",
  length(result2021$genes),
  "\n"
)

cat(
  "Common retained genes     :",
  length(common_genes),
  "\n"
)

cat("------------------------------------------------------------\n")
cat("No CPM or logCPM matrices were generated.\n")
cat("No differential-expression analysis was performed.\n")
cat("No PCA or clustering analysis was performed.\n")
cat("Original input files were not modified.\n")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("PASS\n")
cat(
  "Step 02 is complete. The project can proceed to ",
  "Step 03: 2026 internal expression structure ",
  "(PCA, sample correlation, clustering and outlier diagnostics).\n"
)
