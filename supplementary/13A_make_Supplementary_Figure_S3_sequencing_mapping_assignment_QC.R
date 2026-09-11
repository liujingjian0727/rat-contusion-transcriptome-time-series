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
cat("13A SUPPLEMENTARY FIGURE S3\n")
cat("Raw-read sequencing, mapping and assignment QC\n")
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

raw_stats_file <- file.path(
  DATA_ROOT,
  "CleanData",
  "fastq_stats.txt"
)

hisat2_dir <- file.path(
  DATA_ROOT,
  "hisat2"
)

featurecounts_dir <- file.path(
  DATA_ROOT,
  "run_featurecounts"
)

outdir <- paste0(
  "supplementary_figures/",
  "Supplementary_Figure_S3_sequencing_mapping_assignment_QC"
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
# 1. Stage information
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
# 2. Validate major inputs
# ============================================================

if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found: ",
    metadata_file
  )
}

if (!file.exists(raw_stats_file)) {
  stop(
    "Raw FASTQ statistics file not found: ",
    raw_stats_file
  )
}

if (!dir.exists(hisat2_dir)) {
  stop(
    "HISAT2 directory not found: ",
    hisat2_dir
  )
}

if (!dir.exists(featurecounts_dir)) {
  stop(
    "featureCounts directory not found: ",
    featurecounts_dir
  )
}

cat("Major input paths found: PASS\n\n")


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

unknown_groups <- setdiff(
  unique(meta$group),
  groups
)

if (length(unknown_groups) > 0) {
  stop(
    "Unexpected group(s): ",
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
  "R2h samples            : ",
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
# 4. Helper: map filename/string to sample ID
# ============================================================

find_sample_id <- function(
    text_string,
    sample_ids) {

  x <- basename(
    text_string
  )

  hits <- sample_ids[
    vapply(
      sample_ids,
      function(s) {
        grepl(
          s,
          x,
          fixed = TRUE
        )
      },
      logical(1)
    )
  ]

  if (length(hits) == 0) {
    return(NA_character_)
  }

  # Prefer longest exact sample string if multiple matches
  hits <- hits[
    order(
      nchar(hits),
      decreasing = TRUE
    )
  ]

  hits[1]
}


# ============================================================
# 5. RAW FASTQ statistics
# ============================================================

cat("Reading raw FASTQ statistics...\n")

raw_lines <- readLines(
  raw_stats_file,
  warn = FALSE
)

header_idx <- grep(
  "^file[[:space:]]+format[[:space:]]+type[[:space:]]+num_seqs",
  raw_lines
)

if (length(header_idx) != 1) {
  stop(
    "Could not uniquely identify fastq_stats header."
  )
}

raw_tmp <- tempfile(
  fileext = ".txt"
)

writeLines(
  raw_lines[
    header_idx:length(raw_lines)
  ],
  raw_tmp
)

raw_fastq <- read.table(
  raw_tmp,
  header = TRUE,
  sep = "",
  check.names = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

unlink(raw_tmp)

required_raw_cols <- c(
  "file",
  "num_seqs",
  "sum_len",
  "min_len",
  "avg_len",
  "max_len"
)

if (!all(
  required_raw_cols %in%
    colnames(raw_fastq)
)) {
  stop(
    "FASTQ statistics table lacks required columns."
  )
}

clean_number <- function(x) {
  as.numeric(
    gsub(
      ",",
      "",
      as.character(x),
      fixed = TRUE
    )
  )
}

raw_fastq$num_seqs <-
  clean_number(
    raw_fastq$num_seqs
  )

raw_fastq$sum_len <-
  clean_number(
    raw_fastq$sum_len
  )

raw_fastq$min_len <-
  clean_number(
    raw_fastq$min_len
  )

raw_fastq$avg_len <-
  clean_number(
    raw_fastq$avg_len
  )

raw_fastq$max_len <-
  clean_number(
    raw_fastq$max_len
  )

if (anyNA(
  raw_fastq[
    ,
    c(
      "num_seqs",
      "sum_len",
      "avg_len"
    )
  ]
)) {
  stop(
    "Failed to parse numeric FASTQ statistics."
  )
}


# Sample ID is the directory component before FASTQ filename
raw_fastq$sampleID <- sub(
  "/.*$",
  "",
  raw_fastq$file
)

raw_fastq$mate <- ifelse(
  grepl(
    "_1\\.fq",
    raw_fastq$file
  ),
  "R1",
  ifelse(
    grepl(
      "_2\\.fq",
      raw_fastq$file
    ),
    "R2",
    NA_character_
  )
)


raw_sample_ids <- unique(
  raw_fastq$sampleID
)

raw_extra <- setdiff(
  raw_sample_ids,
  final_samples
)

raw_missing <- setdiff(
  final_samples,
  raw_sample_ids
)

cat(
  "Raw FASTQ sample IDs discovered : ",
  length(raw_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Raw FASTQ extra sample(s)       : ",
  length(raw_extra),
  "\n",
  sep = ""
)

if (length(raw_extra) > 0) {
  cat(
    "  ",
    paste(
      raw_extra,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Final samples missing FASTQ stats: ",
  length(raw_missing),
  "\n\n",
  sep = ""
)

if (length(raw_missing) > 0) {
  stop(
    "Final sample(s) missing from raw FASTQ statistics."
  )
}


raw_summary_list <- vector(
  "list",
  length(final_samples)
)

for (i in seq_along(final_samples)) {

  s <- final_samples[i]

  x <- raw_fastq[
    raw_fastq$sampleID == s,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 2) {
    stop(
      "Expected exactly two FASTQ mates for ",
      s,
      "; observed ",
      nrow(x)
    )
  }

  if (!setequal(
    x$mate,
    c(
      "R1",
      "R2"
    )
  )) {
    stop(
      "Could not identify R1/R2 mates for ",
      s
    )
  }

  r1 <- x[
    x$mate == "R1",
    ,
    drop = FALSE
  ]

  r2 <- x[
    x$mate == "R2",
    ,
    drop = FALSE
  ]

  if (r1$num_seqs !=
      r2$num_seqs) {

    warning(
      s,
      ": R1/R2 read numbers are not identical."
    )
  }

  raw_summary_list[[i]] <- data.frame(
    sampleID = s,

    Raw_R1_reads =
      r1$num_seqs,

    Raw_R2_reads =
      r2$num_seqs,

    Raw_read_pairs =
      min(
        r1$num_seqs,
        r2$num_seqs
      ),

    Raw_reads =
      r1$num_seqs +
      r2$num_seqs,

    Raw_bases =
      r1$sum_len +
      r2$sum_len,

    Raw_bases_Gb =
      (
        r1$sum_len +
        r2$sum_len
      ) / 1e9,

    Read_length_bp =
      weighted.mean(
        c(
          r1$avg_len,
          r2$avg_len
        ),
        w = c(
          r1$num_seqs,
          r2$num_seqs
        )
      ),

    stringsAsFactors = FALSE
  )
}

raw_summary <- do.call(
  rbind,
  raw_summary_list
)

rownames(raw_summary) <- NULL


# ============================================================
# 6. HISAT2 mapping statistics
# ============================================================

cat("Reading HISAT2 alignment logs...\n")

hisat_logs <- list.files(
  hisat2_dir,
  pattern = "\\.log$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(hisat_logs) == 0) {
  stop(
    "No HISAT2 .log files discovered."
  )
}

all_candidate_samples <- union(
  final_samples,
  raw_extra
)

parse_hisat2_log <- function(file) {

  lines <- readLines(
    file,
    warn = FALSE
  )

  align_line <- grep(
    "overall alignment rate",
    lines,
    value = TRUE
  )

  if (length(align_line) == 0) {
    return(NULL)
  }

  align_rate <- suppressWarnings(
    as.numeric(
      sub(
        ".*?([0-9]+\\.?[0-9]*)% overall alignment rate.*",
        "\\1",
        align_line[
          length(align_line)
        ]
      )
    )
  )

  if (!is.finite(align_rate)) {
    return(NULL)
  }

  first_reads_line <- grep(
    "^[[:space:]]*[0-9,]+ reads; of these:",
    lines,
    value = TRUE
  )

  input_fragments <- NA_real_

  if (length(first_reads_line) > 0) {

    input_fragments <- suppressWarnings(
      as.numeric(
        gsub(
          ",",
          "",
          sub(
            "^[[:space:]]*([0-9,]+) reads;.*",
            "\\1",
            first_reads_line[1]
          ),
          fixed = TRUE
        )
      )
    )
  }

  sample_id <- find_sample_id(
    file,
    all_candidate_samples
  )

  if (is.na(sample_id)) {
    return(NULL)
  }

  data.frame(
    sampleID = sample_id,
    HISAT2_input_records =
      input_fragments,
    HISAT2_overall_alignment_rate_percent =
      align_rate,
    HISAT2_log =
      file,
    stringsAsFactors = FALSE
  )
}


hisat_parsed <- lapply(
  hisat_logs,
  parse_hisat2_log
)

hisat_parsed <- Filter(
  Negate(is.null),
  hisat_parsed
)

if (length(hisat_parsed) == 0) {
  stop(
    "No parsable HISAT2 alignment logs found."
  )
}

hisat_table <- do.call(
  rbind,
  hisat_parsed
)

rownames(hisat_table) <- NULL


# Handle duplicate logs per sample
duplicate_hisat <- unique(
  hisat_table$sampleID[
    duplicated(
      hisat_table$sampleID
    )
  ]
)

if (length(duplicate_hisat) > 0) {

  cat(
    "Samples represented by multiple parsable HISAT2 logs:\n"
  )

  cat(
    paste(
      duplicate_hisat,
      collapse = ", "
    ),
    "\n"
  )

  cat(
    "Using the last parsable log for each duplicated sample.\n"
  )
}

hisat_table <- hisat_table[
  !duplicated(
    hisat_table$sampleID,
    fromLast = TRUE
  ),
  ,
  drop = FALSE
]


hisat_samples <- unique(
  hisat_table$sampleID
)

hisat_extra <- setdiff(
  hisat_samples,
  final_samples
)

hisat_missing <- setdiff(
  final_samples,
  hisat_samples
)

cat(
  "Parsable HISAT2 samples           : ",
  length(hisat_samples),
  "\n",
  sep = ""
)

cat(
  "HISAT2 extra sample(s)            : ",
  length(hisat_extra),
  "\n",
  sep = ""
)

if (length(hisat_extra) > 0) {
  cat(
    "  ",
    paste(
      hisat_extra,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Final samples missing HISAT2 rate : ",
  length(hisat_missing),
  "\n\n",
  sep = ""
)

if (length(hisat_missing) > 0) {

  cat(
    "Missing:\n",
    paste(
      hisat_missing,
      collapse = ", "
    ),
    "\n"
  )

  stop(
    "One or more final samples lack parsable HISAT2 mapping rates."
  )
}


hisat_final <- hisat_table[
  match(
    final_samples,
    hisat_table$sampleID
  ),
  ,
  drop = FALSE
]


# ============================================================
# 7. featureCounts assignment statistics
#
# Preferred source:
#   featureCounts_QC_Cohort_2026_cleaned.tsv
#
# Fallback source:
#   per-sample *_featurecounts.log
#
# No featureCounts re-analysis is performed.
# ============================================================

cat("Reading featureCounts assignment statistics...\n")

fc_clean_file <- file.path(
  featurecounts_dir,
  "featureCounts_QC_Cohort_2026_cleaned.tsv"
)


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

canonical_name <- function(x) {

  x <- tolower(x)

  x <- gsub(
    "[^a-z0-9]+",
    "_",
    x
  )

  x <- gsub(
    "^_+|_+$",
    "",
    x
  )

  x
}


numeric_from_text <- function(x) {

  x <- as.character(x)

  x <- gsub(
    ",",
    "",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "%",
    "",
    x,
    fixed = TRUE
  )

  suppressWarnings(
    as.numeric(x)
  )
}


# ------------------------------------------------------------
# Try consolidated cleaned QC table first
# ------------------------------------------------------------

fc_selected <- NULL

if (file.exists(fc_clean_file)) {

  cat(
    "Found consolidated featureCounts QC table:\n",
    fc_clean_file,
    "\n",
    sep = ""
  )

  fc_clean <- try(
    read.delim(
      fc_clean_file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE
    ),
    silent = TRUE
  )


  if (!inherits(
    fc_clean,
    "try-error"
  ) &&
      nrow(fc_clean) > 0 &&
      ncol(fc_clean) >= 2) {

    canonical_cols <- canonical_name(
      colnames(fc_clean)
    )

    cat(
      "featureCounts QC table dimensions : ",
      nrow(fc_clean),
      " rows x ",
      ncol(fc_clean),
      " columns\n",
      sep = ""
    )

    cat(
      "Detected columns:\n  ",
      paste(
        colnames(fc_clean),
        collapse = "\n  "
      ),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # Detect sample column
    # --------------------------------------------------------

    sample_candidates <- which(
      canonical_cols %in%
        c(
          "sample",
          "sampleid",
          "sample_id",
          "sample_name",
          "samplename"
        )
    )


    # If no obvious sample column name exists,
    # identify the column containing the largest number
    # of recognizable sample IDs.
    if (length(sample_candidates) == 0) {

      recognition_scores <- vapply(
        seq_len(
          ncol(fc_clean)
        ),
        function(j) {

          values <- as.character(
            fc_clean[
              ,
              j
            ]
          )

          sum(
            vapply(
              values,
              function(z) {

                !is.na(
                  find_sample_id(
                    z,
                    all_candidate_samples
                  )
                )
              },
              logical(1)
            )
          )
        },
        numeric(1)
      )

      if (max(
        recognition_scores
      ) >= 50) {

        sample_candidates <-
          which.max(
            recognition_scores
          )
      }
    }


    # --------------------------------------------------------
    # Detect assignment-rate column
    # --------------------------------------------------------

    rate_candidates <- grep(
      paste0(
        "assign.*",
        "(rate|percent|percentage|pct)",
        "|",
        "(rate|percent|percentage|pct).*assign"
      ),
      canonical_cols,
      perl = TRUE
    )


    # --------------------------------------------------------
    # Detect assigned-count column
    # --------------------------------------------------------

    assigned_candidates <- grep(
      paste0(
        "^assigned$",
        "|^assigned_",
        "|successfully_assigned"
      ),
      canonical_cols,
      perl = TRUE
    )


    # --------------------------------------------------------
    # Detect total-count column
    # --------------------------------------------------------

    total_candidates <- grep(
      paste0(
        "^total$",
        "|total_.*(read|alignment|fragment|record)",
        "|(read|alignment|fragment|record).*total"
      ),
      canonical_cols,
      perl = TRUE
    )


    if (length(sample_candidates) >= 1) {

      sample_col <- sample_candidates[1]

      sample_values <- as.character(
        fc_clean[
          ,
          sample_col
        ]
      )

      sample_ids_detected <- vapply(
        sample_values,
        function(z) {

          find_sample_id(
            z,
            all_candidate_samples
          )
        },
        character(1)
      )


      # ------------------------------------------------------
      # Assignment rate
      # ------------------------------------------------------

      assignment_rate <- rep(
        NA_real_,
        nrow(fc_clean)
      )

      if (length(rate_candidates) >= 1) {

        rate_col <- rate_candidates[1]

        assignment_rate <-
          numeric_from_text(
            fc_clean[
              ,
              rate_col
            ]
          )

        # Convert fractions to percentages if required.
        finite_rate <- assignment_rate[
          is.finite(
            assignment_rate
          )
        ]

        if (length(finite_rate) > 0 &&
            max(
              finite_rate
            ) <= 1.2) {

          assignment_rate <-
            100 *
            assignment_rate
        }
      }


      # ------------------------------------------------------
      # Assigned count
      # ------------------------------------------------------

      assigned_count <- rep(
        NA_real_,
        nrow(fc_clean)
      )

      if (length(
        assigned_candidates
      ) >= 1) {

        assigned_count <-
          numeric_from_text(
            fc_clean[
              ,
              assigned_candidates[1]
            ]
          )
      }


      # ------------------------------------------------------
      # Total count
      # ------------------------------------------------------

      total_count <- rep(
        NA_real_,
        nrow(fc_clean)
      )

      if (length(
        total_candidates
      ) >= 1) {

        total_count <-
          numeric_from_text(
            fc_clean[
              ,
              total_candidates[1]
            ]
          )
      }


      # If no explicit rate column exists,
      # calculate rate only when assigned and total counts exist.
      if (all(
        is.na(
          assignment_rate
        )
      ) &&
          any(
            is.finite(
              assigned_count
            )
          ) &&
          any(
            is.finite(
              total_count
            )
          )) {

        assignment_rate <-
          100 *
          assigned_count /
          total_count
      }


      fc_candidate <- data.frame(
        sampleID =
          sample_ids_detected,

        featureCounts_Assigned =
          assigned_count,

        featureCounts_Total_records =
          total_count,

        featureCounts_assignment_rate_percent =
          assignment_rate,

        Unassigned_Unmapped =
          NA_real_,

        Unassigned_NoFeatures =
          NA_real_,

        Unassigned_Ambiguity =
          NA_real_,

        Unassigned_MultiMapping =
          NA_real_,

        featureCounts_summary_file =
          fc_clean_file,

        stringsAsFactors = FALSE
      )


      fc_candidate <- fc_candidate[
        !is.na(
          fc_candidate$sampleID
        ) &
          is.finite(
            fc_candidate$
              featureCounts_assignment_rate_percent
          ),
        ,
        drop = FALSE
      ]


      # Remove duplicate representations, if any.
      fc_candidate <- fc_candidate[
        !duplicated(
          fc_candidate$sampleID
        ),
        ,
        drop = FALSE
      ]


      final_recognized <- intersect(
        final_samples,
        fc_candidate$sampleID
      )


      if (length(
        final_recognized
      ) == 77) {

        fc_selected <- fc_candidate

        cat(
          "Consolidated featureCounts QC table: PASS\n"
        )

      } else {

        cat(
          "Consolidated table did not provide all 77 ",
          "final assignment rates.\n"
        )

        cat(
          "Recognized final samples : ",
          length(
            final_recognized
          ),
          "/77\n",
          sep = ""
        )

        cat(
          "Falling back to per-sample ",
          "*_featurecounts.log files.\n"
        )
      }
    }
  }
}


# ============================================================
# 8. Fallback: parse per-sample featureCounts logs
# ============================================================

if (is.null(
  fc_selected
)) {

  cat(
    "Parsing per-sample featureCounts logs...\n"
  )

  fc_logs <- list.files(
    featurecounts_dir,
    pattern =
      "_featurecounts\\.log$",
    full.names = TRUE,
    recursive = FALSE
  )


  if (length(fc_logs) == 0) {
    stop(
      "No *_featurecounts.log files found."
    )
  }


  parse_featurecounts_log <- function(file) {

    sample_id <- sub(
      "_featurecounts\\.log$",
      "",
      basename(file)
    )

    if (!sample_id %in%
        all_candidate_samples) {

      return(NULL)
    }


    lines <- readLines(
      file,
      warn = FALSE
    )


    # Standard featureCounts output:
    #
    # Successfully assigned alignments : 12345678 (85.3%)
    #
    # Depending on version / paired-end mode,
    # wording may refer to alignments, reads or fragments.

    assigned_line <- grep(
      "Successfully assigned",
      lines,
      value = TRUE,
      ignore.case = TRUE
    )


    if (length(
      assigned_line
    ) == 0) {

      return(NULL)
    }


    assigned_line <-
      assigned_line[
        length(
          assigned_line
        )
      ]


    # --------------------------------------------------------
    # Parse assignment percentage
    # --------------------------------------------------------

    pct_match <- regexec(
      "\\(([0-9]+\\.?[0-9]*)%\\)",
      assigned_line,
      perl = TRUE
    )

    pct_parts <- regmatches(
      assigned_line,
      pct_match
    )[[1]]


    assignment_rate <- NA_real_

    if (length(
      pct_parts
    ) >= 2) {

      assignment_rate <-
        suppressWarnings(
          as.numeric(
            pct_parts[2]
          )
        )
    }


    # --------------------------------------------------------
    # Parse successfully assigned count
    # --------------------------------------------------------

    count_match <- regexec(
      paste0(
        "Successfully assigned[^:]*",
        ":[[:space:]]*([0-9,]+)"
      ),
      assigned_line,
      perl = TRUE,
      ignore.case = TRUE
    )

    count_parts <- regmatches(
      assigned_line,
      count_match
    )[[1]]


    assigned_count <- NA_real_

    if (length(
      count_parts
    ) >= 2) {

      assigned_count <-
        suppressWarnings(
          as.numeric(
            gsub(
              ",",
              "",
              count_parts[2],
              fixed = TRUE
            )
          )
        )
    }


    # --------------------------------------------------------
    # Try to parse total input count
    # --------------------------------------------------------

    total_lines <- grep(
      paste0(
        "Total ",
        "(alignments|reads|fragments)"
      ),
      lines,
      value = TRUE,
      ignore.case = TRUE
    )


    total_count <- NA_real_

    if (length(
      total_lines
    ) > 0) {

      total_line <-
        total_lines[
          length(
            total_lines
          )
        ]

      total_match <- regexec(
        ":[[:space:]]*([0-9,]+)",
        total_line,
        perl = TRUE
      )

      total_parts <- regmatches(
        total_line,
        total_match
      )[[1]]


      if (length(
        total_parts
      ) >= 2) {

        total_count <-
          suppressWarnings(
            as.numeric(
              gsub(
                ",",
                "",
                total_parts[2],
                fixed = TRUE
              )
            )
          )
      }
    }


    # Percentage is essential.
    if (!is.finite(
      assignment_rate
    )) {

      return(NULL)
    }


    data.frame(
      sampleID = sample_id,

      featureCounts_Assigned =
        assigned_count,

      featureCounts_Total_records =
        total_count,

      featureCounts_assignment_rate_percent =
        assignment_rate,

      Unassigned_Unmapped =
        NA_real_,

      Unassigned_NoFeatures =
        NA_real_,

      Unassigned_Ambiguity =
        NA_real_,

      Unassigned_MultiMapping =
        NA_real_,

      featureCounts_summary_file =
        file,

      stringsAsFactors = FALSE
    )
  }


  fc_list <- lapply(
    fc_logs,
    parse_featurecounts_log
  )

  fc_list <- Filter(
    Negate(is.null),
    fc_list
  )


  if (length(
    fc_list
  ) == 0) {

    stop(
      "No parsable featureCounts assignment logs found."
    )
  }


  fc_selected <- do.call(
    rbind,
    fc_list
  )

  rownames(fc_selected) <- NULL


  if (anyDuplicated(
    fc_selected$sampleID
  )) {

    stop(
      "Duplicated sample IDs among parsed featureCounts logs."
    )
  }
}


# ============================================================
# 8B. Final featureCounts whitelist audit
# ============================================================

fc_samples <- unique(
  fc_selected$sampleID
)

fc_extra <- setdiff(
  fc_samples,
  final_samples
)

fc_missing <- setdiff(
  final_samples,
  fc_samples
)

cat(
  "Parsable featureCounts samples           : ",
  length(fc_samples),
  "\n",
  sep = ""
)

cat(
  "featureCounts extra sample(s)            : ",
  length(fc_extra),
  "\n",
  sep = ""
)

if (length(
  fc_extra
) > 0) {

  cat(
    "  ",
    paste(
      fc_extra,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
}

cat(
  "Final samples missing assignment metrics : ",
  length(fc_missing),
  "\n",
  sep = ""
)


if (length(
  fc_missing
) > 0) {

  cat(
    "Missing:\n",
    paste(
      fc_missing,
      collapse = ", "
    ),
    "\n"
  )

  stop(
    "One or more final samples lack ",
    "featureCounts assignment metrics."
  )
}


fc_final <- fc_selected[
  match(
    final_samples,
    fc_selected$sampleID
  ),
  ,
  drop = FALSE
]


if (anyNA(
  fc_final$
    featureCounts_assignment_rate_percent
)) {

  stop(
    "NA assignment rate after final featureCounts alignment."
  )
}


# ============================================================
# 9. Merge final 77-sample QC table
# ============================================================

qc <- data.frame(
  sampleID = final_samples,
  group = meta$group,
  Stage = unname(
    stage_labels[
      as.character(
        meta$group
      )
    ]
  ),
  stringsAsFactors = FALSE
)


qc <- merge(
  qc,
  raw_summary,
  by = "sampleID",
  all.x = TRUE,
  sort = FALSE
)

qc <- merge(
  qc,
  hisat_final[
    ,
    c(
      "sampleID",
      "HISAT2_input_records",
      "HISAT2_overall_alignment_rate_percent",
      "HISAT2_log"
    )
  ],
  by = "sampleID",
  all.x = TRUE,
  sort = FALSE
)

qc <- merge(
  qc,
  fc_final[
    ,
    c(
      "sampleID",
      "featureCounts_Assigned",
      "featureCounts_Total_records",
      "featureCounts_assignment_rate_percent",
      "Unassigned_Unmapped",
      "Unassigned_NoFeatures",
      "Unassigned_Ambiguity",
      "Unassigned_MultiMapping",
      "featureCounts_summary_file"
    )
  ],
  by = "sampleID",
  all.x = TRUE,
  sort = FALSE
)


# Restore canonical metadata order
qc <- qc[
  match(
    final_samples,
    qc$sampleID
  ),
  ,
  drop = FALSE
]

qc$group <- factor(
  qc$group,
  levels = groups
)


essential_cols <- c(
  "Raw_read_pairs",
  "Raw_reads",
  "Raw_bases_Gb",
  "Read_length_bp",
  "HISAT2_overall_alignment_rate_percent",
  "featureCounts_assignment_rate_percent"
)

if (anyNA(
  qc[
    ,
    essential_cols
  ]
)) {
  stop(
    "NA detected in final essential QC metrics."
  )
}


# ============================================================
# 10. Final cohort checks
# ============================================================

if (nrow(qc) != 77) {
  stop(
    "Final QC table does not contain 77 samples."
  )
}

if ("R2h_1" %in%
    qc$sampleID) {
  stop(
    "R2h_1 unexpectedly present."
  )
}

if (!("R2h_8" %in%
      qc$sampleID)) {
  stop(
    "R2h_8 unexpectedly absent."
  )
}

cat(
  "Merged final 77-sample QC table: PASS\n\n"
)


# ============================================================
# 11. Additional derived columns
# ============================================================

qc$Raw_read_pairs_million <-
  qc$Raw_read_pairs / 1e6

qc$Raw_reads_million <-
  qc$Raw_reads / 1e6


# ============================================================
# 12. Stage-level summaries
# ============================================================

stage_summary_list <- vector(
  "list",
  length(groups)
)

for (i in seq_along(groups)) {

  g <- groups[i]

  x <- qc[
    qc$group == g,
    ,
    drop = FALSE
  ]

  stage_summary_list[[i]] <- data.frame(
    Group = g,
    Stage = unname(
      stage_labels[g]
    ),
    N = nrow(x),

    Median_raw_read_pairs_million =
      median(
        x$Raw_read_pairs_million
      ),

    Minimum_raw_read_pairs_million =
      min(
        x$Raw_read_pairs_million
      ),

    Maximum_raw_read_pairs_million =
      max(
        x$Raw_read_pairs_million
      ),

    Median_alignment_rate_percent =
      median(
        x$HISAT2_overall_alignment_rate_percent
      ),

    Minimum_alignment_rate_percent =
      min(
        x$HISAT2_overall_alignment_rate_percent
      ),

    Maximum_alignment_rate_percent =
      max(
        x$HISAT2_overall_alignment_rate_percent
      ),

    Median_assignment_rate_percent =
      median(
        x$featureCounts_assignment_rate_percent
      ),

    Minimum_assignment_rate_percent =
      min(
        x$featureCounts_assignment_rate_percent
      ),

    Maximum_assignment_rate_percent =
      max(
        x$featureCounts_assignment_rate_percent
      ),

    stringsAsFactors = FALSE
  )
}

stage_summary <- do.call(
  rbind,
  stage_summary_list
)

rownames(stage_summary) <- NULL


# ============================================================
# 13. Cohort summary
# ============================================================

cohort_summary <- data.frame(
  Metric = c(
    "Raw read pairs per sample (million)",
    "Raw reads per sample (million)",
    "Raw bases per sample (Gb)",
    "Read length (bp)",
    "HISAT2 overall alignment rate (%)",
    "featureCounts assignment rate (%)"
  ),

  Minimum = c(
    min(
      qc$Raw_read_pairs_million
    ),
    min(
      qc$Raw_reads_million
    ),
    min(
      qc$Raw_bases_Gb
    ),
    min(
      qc$Read_length_bp
    ),
    min(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    min(
      qc$featureCounts_assignment_rate_percent
    )
  ),

  Median = c(
    median(
      qc$Raw_read_pairs_million
    ),
    median(
      qc$Raw_reads_million
    ),
    median(
      qc$Raw_bases_Gb
    ),
    median(
      qc$Read_length_bp
    ),
    median(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    median(
      qc$featureCounts_assignment_rate_percent
    )
  ),

  Mean = c(
    mean(
      qc$Raw_read_pairs_million
    ),
    mean(
      qc$Raw_reads_million
    ),
    mean(
      qc$Raw_bases_Gb
    ),
    mean(
      qc$Read_length_bp
    ),
    mean(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    mean(
      qc$featureCounts_assignment_rate_percent
    )
  ),

  Maximum = c(
    max(
      qc$Raw_read_pairs_million
    ),
    max(
      qc$Raw_reads_million
    ),
    max(
      qc$Raw_bases_Gb
    ),
    max(
      qc$Read_length_bp
    ),
    max(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    max(
      qc$featureCounts_assignment_rate_percent
    )
  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 14. Mapping-assignment relationship
# ============================================================

pearson_test <- cor.test(
  qc$HISAT2_overall_alignment_rate_percent,
  qc$featureCounts_assignment_rate_percent,
  method = "pearson"
)

spearman_test <- suppressWarnings(
  cor.test(
    qc$HISAT2_overall_alignment_rate_percent,
    qc$featureCounts_assignment_rate_percent,
    method = "spearman",
    exact = FALSE
  )
)

mapping_assignment_relationship <- data.frame(
  N = nrow(qc),

  Pearson_r =
    unname(
      pearson_test$estimate
    ),

  Pearson_P =
    pearson_test$p.value,

  Spearman_rho =
    unname(
      spearman_test$estimate
    ),

  Spearman_P =
    spearman_test$p.value,

  stringsAsFactors = FALSE
)


# ============================================================
# 15. Publication theme
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
      size = 8.0
    ),

    legend.title = element_text(
      size = 8.5
    ),

    legend.text = element_text(
      size = 7.5
    ),

    legend.key.height = unit(
      0.33,
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
# 16. Stage QC panel helper
# ============================================================

make_stage_qc_panel <- function(
    y_column,
    y_label,
    title,
    subtitle,
    tag) {

  plot_data <- data.frame(
    sampleID = qc$sampleID,
    group = qc$group,
    y = qc[
      ,
      y_column
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
        vjust = 1,
        size = 7.5
      )
    )
}


# ============================================================
# 17. Panels A-C
# ============================================================

pA <- make_stage_qc_panel(
  y_column =
    "Raw_read_pairs_million",

  y_label =
    "Raw read pairs (million)",

  title =
    "Raw sequencing depth",

  subtitle =
    "Paired-end FASTQ files supplied by the sequencing provider",

  tag = "A"
)


pB <- make_stage_qc_panel(
  y_column =
    "HISAT2_overall_alignment_rate_percent",

  y_label =
    "Overall alignment rate (%)",

  title =
    "HISAT2 mapping efficiency",

  subtitle =
    "Overall alignment rate reported by HISAT2",

  tag = "B"
)


pC <- make_stage_qc_panel(
  y_column =
    "featureCounts_assignment_rate_percent",

  y_label =
    "Assignment rate (%)",

  title =
    "Gene-level read assignment",

  subtitle =
    "Assignment rate reported by featureCounts",

  tag = "C"
)


# ============================================================
# 18. Panel D
# ============================================================

format_p <- function(p) {

  if (p < 0.001) {
    return("<0.001")
  }

  paste0(
    "= ",
    sprintf(
      "%.3f",
      p
    )
  )
}


relationship_label <- paste0(
  "Spearman \u03c1 = ",
  sprintf(
    "%.3f",
    mapping_assignment_relationship$
      Spearman_rho
  ),
  ", P ",
  format_p(
    mapping_assignment_relationship$
      Spearman_P
  )
)


pD <- ggplot(
  qc,
  aes(
    x = HISAT2_overall_alignment_rate_percent,
    y = featureCounts_assignment_rate_percent,
    colour = group
  )
) +

  geom_point(
    size = 2.5,
    alpha = 0.85
  ) +

  geom_smooth(
    aes(
      group = 1
    ),
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = "black",
    fill = "grey75",
    linewidth = 0.65,
    alpha = 0.22
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

  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label =
      relationship_label,
    hjust = -0.05,
    vjust = 1.10,
    size = 2.8,
    fill = "white",
    colour = "grey15",
    linewidth = 0.25
  ) +

  labs(
    tag = "D",
    title =
      "Relationship between mapping and gene assignment",
    subtitle =
      "Each point represents one final primary-cohort sample",
    x =
      "HISAT2 overall alignment rate (%)",
    y =
      "featureCounts assignment rate (%)"
  ) +

  theme_pub +

  theme(
    legend.position = "right"
  )


# ============================================================
# 19. Saving helpers
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
    width = 6.3,
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

  print(plot_object)
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

  print(plot_object)
  dev.off()
}


save_panel(
  pA,
  "Supplementary_Figure_S3A_raw_read_pairs"
)

save_panel(
  pB,
  "Supplementary_Figure_S3B_HISAT2_alignment"
)

save_panel(
  pC,
  "Supplementary_Figure_S3C_featureCounts_assignment"
)

save_panel(
  pD,
  "Supplementary_Figure_S3D_mapping_vs_assignment"
)


# ============================================================
# 20. Combined Supplementary Figure S3
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
  paste0(
    "Supplementary_Figure_S3_",
    "sequencing_mapping_assignment_QC.pdf"
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
    "Supplementary_Figure_S3_",
    "sequencing_mapping_assignment_QC_600dpi.tiff"
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
# 21. Sample audit table
# ============================================================

all_audit_samples <- unique(
  c(
    final_samples,
    raw_sample_ids,
    hisat_samples,
    fc_samples
  )
)

audit <- data.frame(
  sampleID =
    all_audit_samples,

  In_final_77_metadata =
    all_audit_samples %in%
      final_samples,

  Raw_FASTQ_stats_present =
    all_audit_samples %in%
      raw_sample_ids,

  HISAT2_mapping_metric_present =
    all_audit_samples %in%
      hisat_samples,

  featureCounts_assignment_metric_present =
    all_audit_samples %in%
      fc_samples,

  stringsAsFactors = FALSE
)

audit$Final_status <- ifelse(
  audit$In_final_77_metadata,
  "Included_final77",
  "Excluded_not_in_final_metadata"
)


# ============================================================
# 22. Save data tables
# ============================================================

write.table(
  qc,
  file = file.path(
    outdir,
    paste0(
      "Supplementary_Table_S2_",
      "2026_primary_77sample_",
      "sequencing_mapping_assignment_QC.tsv"
    )
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  stage_summary,
  file = file.path(
    outdir,
    "S3_stage_level_QC_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  cohort_summary,
  file = file.path(
    outdir,
    "S3_cohort_level_QC_summary.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  mapping_assignment_relationship,
  file = file.path(
    outdir,
    "S3_mapping_assignment_relationship.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


write.table(
  audit,
  file = file.path(
    outdir,
    "S3_sample_whitelist_audit.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 23. Caption
# ============================================================

caption_text <- paste0(
  "Supplementary Figure S3 | Raw-read sequencing depth, ",
  "alignment and gene-assignment quality control for the ",
  "2026 primary cohort. ",
  "a, Number of paired-end raw read pairs per sample. ",
  "The FASTQ files were supplied by the sequencing provider ",
  "as the original sequencing reads, irrespective of the ",
  "directory name used for data delivery. ",
  "b, HISAT2 overall alignment rates across the 77 final samples. ",
  "c, Gene-level assignment rates reported by featureCounts ",
  "for each sample. ",
  "d, Relationship between HISAT2 overall alignment rate and ",
  "featureCounts assignment rate. ",
  "Boxplots summarize within-stage distributions and points ",
  "represent individual samples. Dashed horizontal lines in ",
  "panels a-c indicate cohort medians. Only samples included ",
  "in the final 77-sample primary cohort are shown."
)

writeLines(
  caption_text,
  con = file.path(
    outdir,
    "Supplementary_Figure_S3_caption.txt"
  )
)


# ============================================================
# 24. Save plot objects / session info
# ============================================================

saveRDS(
  list(
    S3A = pA,
    S3B = pB,
    S3C = pC,
    S3D = pD
  ),
  file = file.path(
    outdir,
    "Supplementary_Figure_S3_plot_objects.rds"
  )
)

sink(
  file.path(
    outdir,
    "Supplementary_Figure_S3_R_sessionInfo.txt"
  )
)

sessionInfo()

sink()


# ============================================================
# 25. Final report
# ============================================================

cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE S3 COMPLETED\n")
cat("============================================================\n\n")

cat(
  "Final samples                         : 77\n"
)

cat(
  "R2h_1 excluded                       : ",
  !("R2h_1" %in%
      qc$sampleID),
  "\n",
  sep = ""
)

cat(
  "R2h_8 retained                       : ",
  "R2h_8" %in%
    qc$sampleID,
  "\n\n",
  sep = ""
)


cat("Source audit:\n")

cat(
  "Raw FASTQ sample IDs                 : ",
  length(raw_sample_ids),
  "\n",
  sep = ""
)

cat(
  "Raw FASTQ extra                     : ",
  ifelse(
    length(raw_extra) == 0,
    "None",
    paste(
      raw_extra,
      collapse = ", "
    )
  ),
  "\n",
  sep = ""
)

cat(
  "HISAT2 sample IDs                    : ",
  length(hisat_samples),
  "\n",
  sep = ""
)

cat(
  "HISAT2 extra                        : ",
  ifelse(
    length(hisat_extra) == 0,
    "None",
    paste(
      hisat_extra,
      collapse = ", "
    )
  ),
  "\n",
  sep = ""
)

cat(
  "featureCounts sample IDs             : ",
  length(fc_samples),
  "\n",
  sep = ""
)

cat(
  "featureCounts extra                 : ",
  ifelse(
    length(fc_extra) == 0,
    "None",
    paste(
      fc_extra,
      collapse = ", "
    )
  ),
  "\n\n",
  sep = ""
)


cat(
  "Raw read pairs/sample (million)\n"
)

cat(
  "  min / median / mean / max : ",
  sprintf(
    "%.3f / %.3f / %.3f / %.3f\n",
    min(
      qc$Raw_read_pairs_million
    ),
    median(
      qc$Raw_read_pairs_million
    ),
    mean(
      qc$Raw_read_pairs_million
    ),
    max(
      qc$Raw_read_pairs_million
    )
  ),
  sep = ""
)


cat(
  "Raw bases/sample (Gb)\n"
)

cat(
  "  min / median / mean / max : ",
  sprintf(
    "%.3f / %.3f / %.3f / %.3f\n",
    min(
      qc$Raw_bases_Gb
    ),
    median(
      qc$Raw_bases_Gb
    ),
    mean(
      qc$Raw_bases_Gb
    ),
    max(
      qc$Raw_bases_Gb
    )
  ),
  sep = ""
)


cat(
  "HISAT2 overall alignment rate (%)\n"
)

cat(
  "  min / median / mean / max : ",
  sprintf(
    "%.3f / %.3f / %.3f / %.3f\n",
    min(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    median(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    mean(
      qc$HISAT2_overall_alignment_rate_percent
    ),
    max(
      qc$HISAT2_overall_alignment_rate_percent
    )
  ),
  sep = ""
)


cat(
  "featureCounts assignment rate (%)\n"
)

cat(
  "  min / median / mean / max : ",
  sprintf(
    "%.3f / %.3f / %.3f / %.3f\n",
    min(
      qc$featureCounts_assignment_rate_percent
    ),
    median(
      qc$featureCounts_assignment_rate_percent
    ),
    mean(
      qc$featureCounts_assignment_rate_percent
    ),
    max(
      qc$featureCounts_assignment_rate_percent
    )
  ),
  sep = ""
)


cat("\nMapping-assignment relationship:\n")

print(
  mapping_assignment_relationship,
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
  "\n",
  sep = ""
)

cat("\nPASS\n")
