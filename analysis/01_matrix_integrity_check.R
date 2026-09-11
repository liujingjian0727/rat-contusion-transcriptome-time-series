#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

cat("\n")
cat("============================================================\n")
cat("01 MATRIX INTEGRITY CHECK\n")
cat("Scientific Data rat skeletal muscle transcriptome project\n")
cat("============================================================\n\n")

# ============================================================
# 0. Input files
# ============================================================

files <- list(
  cohort2026 = list(
    counts = "rattus_counts_2026_primary_77samples.txt",
    tpm    = "rattus_77sample.TMM.TPM.tsv",
    meta   = "rattus_meta_2026_77samples.tsv",
    expected_n = 77
  ),
  cohort2021 = list(
    counts = "rattus_2021.counts.tsv",
    tpm    = "rattus_2021.TMM.TPM.tsv",
    meta   = "rattus_meta_2021.tsv",
    expected_n = 60
  )
)

outdir <- "01_matrix_integrity_check_result"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. Helper functions
# ============================================================

read_matrix_file <- function(file) {

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
    stop("Matrix has fewer than 2 columns: ", file)
  }

  colnames(x)[1] <- "GeneID"

  return(x)
}


read_meta_file <- function(file) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE
  )

  required <- c("sampleID", "group", "time", "condition")

  missing_cols <- setdiff(required, colnames(x))

  if (length(missing_cols) > 0) {
    stop(
      "Metadata file ",
      file,
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  return(x)
}


check_numeric_matrix <- function(x, matrix_name, integer_required = FALSE) {

  sample_cols <- colnames(x)[-1]

  invalid_cells <- 0L
  na_cells <- 0L
  inf_cells <- 0L
  negative_cells <- 0L
  noninteger_cells <- 0L

  for (nm in sample_cols) {

    original <- x[[nm]]

    suppressWarnings(numeric_values <- as.numeric(original))

    invalid <- is.na(numeric_values) & !is.na(original)

    # Empty strings should also count as invalid/missing
    if (is.character(original)) {
      invalid <- invalid | trimws(original) == ""
    }

    invalid_cells <- invalid_cells + sum(invalid, na.rm = TRUE)
    na_cells <- na_cells + sum(is.na(numeric_values))
    inf_cells <- inf_cells + sum(is.infinite(numeric_values), na.rm = TRUE)

    negative_cells <- negative_cells +
      sum(numeric_values < 0, na.rm = TRUE)

    if (integer_required) {
      noninteger_cells <- noninteger_cells +
        sum(
          is.finite(numeric_values) &
          abs(numeric_values - round(numeric_values)) > 1e-8,
          na.rm = TRUE
        )
    }
  }

  data.frame(
    Matrix = matrix_name,
    Invalid_nonnumeric_cells = invalid_cells,
    NA_cells = na_cells,
    Infinite_cells = inf_cells,
    Negative_cells = negative_cells,
    Noninteger_cells = if (integer_required) noninteger_cells else NA_integer_,
    stringsAsFactors = FALSE
  )
}


safe_text <- function(x) {
  if (length(x) == 0) "None" else paste(x, collapse = ",")
}


bool_text <- function(x) {
  ifelse(isTRUE(x), "TRUE", "FALSE")
}


# ============================================================
# 2. Main cohort checking function
# ============================================================

check_cohort <- function(cohort_name,
                         counts_file,
                         tpm_file,
                         meta_file,
                         expected_n) {

  cat("------------------------------------------------------------\n")
  cat("Processing:", cohort_name, "\n")
  cat("------------------------------------------------------------\n")

  # ----------------------------------------------------------
  # Read data
  # ----------------------------------------------------------

  cat("Reading counts ...\n")
  counts <- read_matrix_file(counts_file)

  cat("Reading TMM.TPM ...\n")
  tpm <- read_matrix_file(tpm_file)

  cat("Reading metadata ...\n")
  meta <- read_meta_file(meta_file)

  count_samples <- colnames(counts)[-1]
  tpm_samples <- colnames(tpm)[-1]
  meta_samples <- meta$sampleID

  count_genes <- counts$GeneID
  tpm_genes <- tpm$GeneID

  cat("\nDimensions:\n")
  cat(
    "Counts :",
    nrow(counts), "genes x",
    length(count_samples), "samples\n"
  )
  cat(
    "TPM    :",
    nrow(tpm), "genes x",
    length(tpm_samples), "samples\n"
  )
  cat(
    "Meta   :",
    nrow(meta), "samples\n\n"
  )

  # ----------------------------------------------------------
  # Sample checks
  # ----------------------------------------------------------

  duplicate_count_samples <- unique(
    count_samples[duplicated(count_samples)]
  )

  duplicate_tpm_samples <- unique(
    tpm_samples[duplicated(tpm_samples)]
  )

  duplicate_meta_samples <- unique(
    meta_samples[duplicated(meta_samples)]
  )

  counts_vs_tpm_set <- setequal(count_samples, tpm_samples)
  counts_vs_meta_set <- setequal(count_samples, meta_samples)
  tpm_vs_meta_set <- setequal(tpm_samples, meta_samples)

  counts_vs_tpm_order <- identical(count_samples, tpm_samples)
  counts_vs_meta_order <- identical(count_samples, meta_samples)

  missing_in_tpm <- setdiff(count_samples, tpm_samples)
  extra_in_tpm <- setdiff(tpm_samples, count_samples)

  missing_in_meta <- setdiff(count_samples, meta_samples)
  extra_in_meta <- setdiff(meta_samples, count_samples)

  # ----------------------------------------------------------
  # Gene checks
  # ----------------------------------------------------------

  duplicate_count_genes <- unique(
    count_genes[duplicated(count_genes)]
  )

  duplicate_tpm_genes <- unique(
    tpm_genes[duplicated(tpm_genes)]
  )

  counts_vs_tpm_gene_set <- setequal(count_genes, tpm_genes)
  counts_vs_tpm_gene_order <- identical(count_genes, tpm_genes)

  missing_genes_in_tpm <- setdiff(count_genes, tpm_genes)
  extra_genes_in_tpm <- setdiff(tpm_genes, count_genes)

  # ----------------------------------------------------------
  # Missing GeneID checks
  # ----------------------------------------------------------

  count_missing_geneid <- sum(
    is.na(count_genes) | trimws(count_genes) == ""
  )

  tpm_missing_geneid <- sum(
    is.na(tpm_genes) | trimws(tpm_genes) == ""
  )

  # ----------------------------------------------------------
  # Numeric checks
  # ----------------------------------------------------------

  cat("Checking numerical integrity of counts ...\n")

  counts_numeric <- check_numeric_matrix(
    counts,
    paste0(cohort_name, "_counts"),
    integer_required = TRUE
  )

  cat("Checking numerical integrity of TMM.TPM ...\n")

  tpm_numeric <- check_numeric_matrix(
    tpm,
    paste0(cohort_name, "_TMM_TPM"),
    integer_required = FALSE
  )

  # ----------------------------------------------------------
  # Metadata checks
  # ----------------------------------------------------------

  missing_group <- sum(
    is.na(meta$group) | trimws(meta$group) == ""
  )

  missing_condition <- sum(
    is.na(meta$condition) | trimws(meta$condition) == ""
  )

  # time is allowed to be NA for Baseline
  injury_missing_time <- sum(
    meta$condition == "Injury" &
    (
      is.na(meta$time) |
      trimws(as.character(meta$time)) == ""
    ),
    na.rm = TRUE
  )

  control_nonmissing_time <- sum(
    meta$condition == "Control" &
    !is.na(meta$time) &
    trimws(as.character(meta$time)) != "",
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Expected sample number
  # ----------------------------------------------------------

  expected_counts_n <- length(count_samples) == expected_n
  expected_tpm_n <- length(tpm_samples) == expected_n
  expected_meta_n <- nrow(meta) == expected_n

  # ----------------------------------------------------------
  # Sample matching table
  # ----------------------------------------------------------

  all_samples <- unique(
    c(count_samples, tpm_samples, meta_samples)
  )

  sample_match <- data.frame(
    SampleID = all_samples,
    In_counts = all_samples %in% count_samples,
    In_TMM_TPM = all_samples %in% tpm_samples,
    In_metadata = all_samples %in% meta_samples,
    stringsAsFactors = FALSE
  )

  write.table(
    sample_match,
    file = file.path(
      outdir,
      paste0(cohort_name, "_sample_matching.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # Gene matching table
  # Only write mismatches to avoid huge redundant files
  # ----------------------------------------------------------

  gene_mismatch <- data.frame(
    GeneID = unique(c(
      missing_genes_in_tpm,
      extra_genes_in_tpm
    )),
    stringsAsFactors = FALSE
  )

  if (nrow(gene_mismatch) > 0) {

    gene_mismatch$In_counts <-
      gene_mismatch$GeneID %in% count_genes

    gene_mismatch$In_TMM_TPM <-
      gene_mismatch$GeneID %in% tpm_genes

  } else {

    gene_mismatch <- data.frame(
      GeneID = character(0),
      In_counts = logical(0),
      In_TMM_TPM = logical(0)
    )
  }

  write.table(
    gene_mismatch,
    file = file.path(
      outdir,
      paste0(cohort_name, "_gene_mismatches.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # Metadata group summary
  # ----------------------------------------------------------

  group_summary <- as.data.frame(
    table(
      Group = meta$group,
      Condition = meta$condition,
      useNA = "ifany"
    )
  )

  group_summary <- group_summary[
    group_summary$Freq > 0,
    ,
    drop = FALSE
  ]

  write.table(
    group_summary,
    file = file.path(
      outdir,
      paste0(cohort_name, "_group_summary.tsv")
    ),
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  # ----------------------------------------------------------
  # Determine overall PASS / FAIL
  # ----------------------------------------------------------

  overall_pass <- all(
    expected_counts_n,
    expected_tpm_n,
    expected_meta_n,

    length(duplicate_count_samples) == 0,
    length(duplicate_tpm_samples) == 0,
    length(duplicate_meta_samples) == 0,

    length(duplicate_count_genes) == 0,
    length(duplicate_tpm_genes) == 0,

    counts_vs_tpm_set,
    counts_vs_meta_set,
    tpm_vs_meta_set,

    counts_vs_tpm_gene_set,

    count_missing_geneid == 0,
    tpm_missing_geneid == 0,

    counts_numeric$Invalid_nonnumeric_cells == 0,
    counts_numeric$NA_cells == 0,
    counts_numeric$Infinite_cells == 0,
    counts_numeric$Negative_cells == 0,
    counts_numeric$Noninteger_cells == 0,

    tpm_numeric$Invalid_nonnumeric_cells == 0,
    tpm_numeric$NA_cells == 0,
    tpm_numeric$Infinite_cells == 0,
    tpm_numeric$Negative_cells == 0,

    missing_group == 0,
    missing_condition == 0,
    injury_missing_time == 0
  )

  # ----------------------------------------------------------
  # Summary table
  # ----------------------------------------------------------

  summary <- data.frame(
    Cohort = cohort_name,

    Expected_samples = expected_n,
    Count_samples = length(count_samples),
    TPM_samples = length(tpm_samples),
    Metadata_samples = nrow(meta),

    Count_genes = nrow(counts),
    TPM_genes = nrow(tpm),

    Duplicate_count_samples = length(duplicate_count_samples),
    Duplicate_TPM_samples = length(duplicate_tpm_samples),
    Duplicate_metadata_samples = length(duplicate_meta_samples),

    Duplicate_count_genes = length(duplicate_count_genes),
    Duplicate_TPM_genes = length(duplicate_tpm_genes),

    Count_TPM_sample_sets_identical = counts_vs_tpm_set,
    Count_metadata_sample_sets_identical = counts_vs_meta_set,
    TPM_metadata_sample_sets_identical = tpm_vs_meta_set,

    Count_TPM_sample_order_identical = counts_vs_tpm_order,
    Count_metadata_sample_order_identical = counts_vs_meta_order,

    Count_TPM_gene_sets_identical = counts_vs_tpm_gene_set,
    Count_TPM_gene_order_identical = counts_vs_tpm_gene_order,

    Missing_GeneID_counts = count_missing_geneid,
    Missing_GeneID_TPM = tpm_missing_geneid,

    Count_invalid_cells =
      counts_numeric$Invalid_nonnumeric_cells,

    Count_NA_cells =
      counts_numeric$NA_cells,

    Count_infinite_cells =
      counts_numeric$Infinite_cells,

    Count_negative_cells =
      counts_numeric$Negative_cells,

    Count_noninteger_cells =
      counts_numeric$Noninteger_cells,

    TPM_invalid_cells =
      tpm_numeric$Invalid_nonnumeric_cells,

    TPM_NA_cells =
      tpm_numeric$NA_cells,

    TPM_infinite_cells =
      tpm_numeric$Infinite_cells,

    TPM_negative_cells =
      tpm_numeric$Negative_cells,

    Metadata_missing_group = missing_group,
    Metadata_missing_condition = missing_condition,
    Injury_samples_missing_time = injury_missing_time,
    Control_samples_with_numeric_time = control_nonmissing_time,

    Overall_status = ifelse(overall_pass, "PASS", "FAIL"),

    stringsAsFactors = FALSE
  )

  # ----------------------------------------------------------
  # Print detailed console report
  # ----------------------------------------------------------

  cat("\nSample integrity:\n")
  cat(
    "Expected sample number       :",
    expected_n, "\n"
  )
  cat(
    "Counts sample number         :",
    length(count_samples), "\n"
  )
  cat(
    "TMM.TPM sample number        :",
    length(tpm_samples), "\n"
  )
  cat(
    "Metadata sample number       :",
    nrow(meta), "\n"
  )

  cat("\nSample matching:\n")
  cat(
    "Counts vs TPM same set       :",
    bool_text(counts_vs_tpm_set), "\n"
  )
  cat(
    "Counts vs metadata same set  :",
    bool_text(counts_vs_meta_set), "\n"
  )
  cat(
    "TPM vs metadata same set     :",
    bool_text(tpm_vs_meta_set), "\n"
  )
  cat(
    "Counts vs TPM same order     :",
    bool_text(counts_vs_tpm_order), "\n"
  )
  cat(
    "Counts vs metadata same order:",
    bool_text(counts_vs_meta_order), "\n"
  )

  cat("\nGene integrity:\n")
  cat(
    "Counts genes                 :",
    length(count_genes), "\n"
  )
  cat(
    "TMM.TPM genes                :",
    length(tpm_genes), "\n"
  )
  cat(
    "Counts vs TPM same gene set  :",
    bool_text(counts_vs_tpm_gene_set), "\n"
  )
  cat(
    "Counts vs TPM same gene order:",
    bool_text(counts_vs_tpm_gene_order), "\n"
  )

  cat("\nDuplicates:\n")
  cat(
    "Duplicate count samples      :",
    length(duplicate_count_samples), "\n"
  )
  cat(
    "Duplicate TPM samples        :",
    length(duplicate_tpm_samples), "\n"
  )
  cat(
    "Duplicate metadata samples   :",
    length(duplicate_meta_samples), "\n"
  )
  cat(
    "Duplicate count GeneIDs      :",
    length(duplicate_count_genes), "\n"
  )
  cat(
    "Duplicate TPM GeneIDs        :",
    length(duplicate_tpm_genes), "\n"
  )

  cat("\nCounts numerical integrity:\n")
  cat(
    "Invalid nonnumeric cells     :",
    counts_numeric$Invalid_nonnumeric_cells, "\n"
  )
  cat(
    "NA cells                     :",
    counts_numeric$NA_cells, "\n"
  )
  cat(
    "Infinite cells               :",
    counts_numeric$Infinite_cells, "\n"
  )
  cat(
    "Negative cells               :",
    counts_numeric$Negative_cells, "\n"
  )
  cat(
    "Noninteger cells             :",
    counts_numeric$Noninteger_cells, "\n"
  )

  cat("\nTMM.TPM numerical integrity:\n")
  cat(
    "Invalid nonnumeric cells     :",
    tpm_numeric$Invalid_nonnumeric_cells, "\n"
  )
  cat(
    "NA cells                     :",
    tpm_numeric$NA_cells, "\n"
  )
  cat(
    "Infinite cells               :",
    tpm_numeric$Infinite_cells, "\n"
  )
  cat(
    "Negative cells               :",
    tpm_numeric$Negative_cells, "\n"
  )

  cat("\nMetadata integrity:\n")
  cat(
    "Missing group                :",
    missing_group, "\n"
  )
  cat(
    "Missing condition            :",
    missing_condition, "\n"
  )
  cat(
    "Injury samples missing time  :",
    injury_missing_time, "\n"
  )
  cat(
    "Control samples with time    :",
    control_nonmissing_time, "\n"
  )

  if (length(missing_in_tpm) > 0) {
    cat(
      "\nSamples in counts but not TPM:",
      safe_text(missing_in_tpm), "\n"
    )
  }

  if (length(extra_in_tpm) > 0) {
    cat(
      "Samples in TPM but not counts:",
      safe_text(extra_in_tpm), "\n"
    )
  }

  if (length(missing_in_meta) > 0) {
    cat(
      "Samples in counts but not metadata:",
      safe_text(missing_in_meta), "\n"
    )
  }

  if (length(extra_in_meta) > 0) {
    cat(
      "Samples in metadata but not counts:",
      safe_text(extra_in_meta), "\n"
    )
  }

  cat("\nGroup summary:\n")
  print(group_summary, row.names = FALSE)

  cat("\n")
  cat("FINAL STATUS:", ifelse(overall_pass, "PASS", "FAIL"), "\n")
  cat("------------------------------------------------------------\n\n")

  return(
    list(
      summary = summary,
      numeric = rbind(counts_numeric, tpm_numeric),
      pass = overall_pass
    )
  )
}


# ============================================================
# 3. Confirm all files exist
# ============================================================

all_input_files <- c(
  files$cohort2026$counts,
  files$cohort2026$tpm,
  files$cohort2026$meta,
  files$cohort2021$counts,
  files$cohort2021$tpm,
  files$cohort2021$meta
)

missing_files <- all_input_files[!file.exists(all_input_files)]

if (length(missing_files) > 0) {

  cat("ERROR: Missing input files:\n")

  for (f in missing_files) {
    cat("  -", f, "\n")
  }

  quit(status = 1)
}

cat("All six input files were found.\n\n")


# ============================================================
# 4. Run 2026
# ============================================================

res2026 <- check_cohort(
  cohort_name = "2026_primary",
  counts_file = files$cohort2026$counts,
  tpm_file = files$cohort2026$tpm,
  meta_file = files$cohort2026$meta,
  expected_n = files$cohort2026$expected_n
)


# ============================================================
# 5. Run 2021
# ============================================================

res2021 <- check_cohort(
  cohort_name = "2021_historical_reference",
  counts_file = files$cohort2021$counts,
  tpm_file = files$cohort2021$tpm,
  meta_file = files$cohort2021$meta,
  expected_n = files$cohort2021$expected_n
)


# ============================================================
# 6. Combined summary
# ============================================================

combined_summary <- rbind(
  res2026$summary,
  res2021$summary
)

write.table(
  combined_summary,
  file = file.path(
    outdir,
    "01_matrix_integrity_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

numeric_summary <- rbind(
  res2026$numeric,
  res2021$numeric
)

write.table(
  numeric_summary,
  file = file.path(
    outdir,
    "01_numeric_integrity_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 7. Session information
# ============================================================

sink(
  file.path(
    outdir,
    "01_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 8. Final project-level report
# ============================================================

project_pass <- res2026$pass && res2021$pass

cat("\n")
cat("============================================================\n")
cat("PROJECT-LEVEL MATRIX INTEGRITY CHECK\n")
cat("============================================================\n")
cat("2026 primary cohort             : ",
    ifelse(res2026$pass, "PASS", "FAIL"), "\n", sep = "")
cat("2021 historical reference      : ",
    ifelse(res2021$pass, "PASS", "FAIL"), "\n", sep = "")
cat("------------------------------------------------------------\n")
cat("OVERALL                         : ",
    ifelse(project_pass, "PASS", "FAIL"), "\n", sep = "")
cat("============================================================\n\n")

cat("Output directory:\n")
cat(outdir, "\n\n")

cat("Generated files:\n")
cat(
  paste0(
    "  - ",
    list.files(outdir),
    collapse = "\n"
  ),
  "\n\n"
)

if (project_pass) {

  cat("PASS\n")
  cat(
    "Both cohorts passed matrix, metadata, sample-ID, ",
    "GeneID and numerical-integrity checks.\n",
    sep = ""
  )
  cat(
    "The project can proceed to Step 02: ",
    "count-based gene filtering and DESeq2 VST.\n"
  )

  quit(status = 0)

} else {

  cat("FAIL\n")
  cat(
    "At least one integrity issue was detected. ",
    "Review 01_matrix_integrity_summary.tsv before Step 02.\n"
  )

  quit(status = 1)
}
