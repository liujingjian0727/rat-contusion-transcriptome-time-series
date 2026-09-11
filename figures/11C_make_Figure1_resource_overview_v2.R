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
cat("11C v2 MAKE FIGURE 1\n")
cat("Experimental design and data resource overview\n")
cat("Scientific Data final main figure\n")
cat("============================================================\n\n")


# ============================================================
# 0. INPUT / OUTPUT
# ============================================================

meta26_file <- "rattus_meta_2026_77samples.tsv"
meta21_file <- "rattus_meta_2021.tsv"

for (f in c(meta26_file, meta21_file)) {
  if (!file.exists(f)) {
    stop("Missing input file: ", f)
  }
}

outdir <- "final_submission_figures"
panel_dir <- file.path(outdir, "panels")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 1. FIXED STAGE DEFINITIONS
# ============================================================

groups26 <- c(
  "Baseline",
  "R0h", "R1h", "R2h", "R6h",
  "R10h", "R14h", "R18h",
  "R36h", "R60h", "R72h"
)

labels26 <- c(
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

groups21 <- c(
  "Baseline",
  "R4h", "R8h", "R12h",
  "R16h", "R20h", "R24h", "R48h"
)

labels21 <- c(
  Baseline = "Baseline",
  R4h = "4 h",
  R8h = "8 h",
  R12h = "12 h",
  R16h = "16 h",
  R20h = "20 h",
  R24h = "24 h",
  R48h = "48 h"
)


# ============================================================
# 2. READ METADATA
# ============================================================

read_meta <- function(file) {

  x <- read.delim(
    file,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    quote = "",
    stringsAsFactors = FALSE
  )

  if (!all(c("sampleID", "group") %in% colnames(x))) {
    stop(
      "Metadata must contain sampleID and group: ",
      file
    )
  }

  x
}

meta26 <- read_meta(meta26_file)
meta21 <- read_meta(meta21_file)


make_summary <- function(meta, groups, labels) {

  x <- table(
    factor(
      meta$group,
      levels = groups
    )
  )

  data.frame(
    Group = names(x),
    Stage = unname(labels[names(x)]),
    N = as.integer(x),
    stringsAsFactors = FALSE
  )
}

sum26 <- make_summary(
  meta26,
  groups26,
  labels26
)

sum21 <- make_summary(
  meta21,
  groups21,
  labels21
)

if (sum(sum26$N) != 77) {
  stop(
    "Expected 77 samples in 2026; observed ",
    sum(sum26$N)
  )
}

if (sum(sum21$N) != 60) {
  stop(
    "Expected 60 samples in 2021; observed ",
    sum(sum21$N)
  )
}

n26_baseline <- sum26$N[
  sum26$Group == "Baseline"
]

n21_baseline <- sum21$N[
  sum21$Group == "Baseline"
]

n26_injury <- sum(sum26$N) - n26_baseline
n21_injury <- sum(sum21$N) - n21_baseline

cat("Metadata validation: PASS\n")
cat("2026 samples: 77\n")
cat("2021 samples: 60\n\n")


# ============================================================
# 3. GLOBAL STYLE
# ============================================================

blue <- "#4C78A8"
orange <- "#D18455"
green <- "#72A66A"
purple <- "#8C78A8"

light_blue <- "#E7F0F8"
light_orange <- "#F8EBDD"
light_green <- "#E8F2E5"
light_purple <- "#EEEAF4"

dark <- "#333333"

theme_void_pub <- theme_void() +
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    plot.title = element_text(
      face = "bold",
      size = 11.5
    ),
    plot.margin = margin(
      8, 8, 8, 8
    )
  )


# ============================================================
# 4. SAVING FUNCTIONS
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
      paste0(stem, ".pdf")
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
    pC,
    pD,
    width = 13,
    height = 10) {

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
            c(0.93, 1.07),
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
    "Figure_1_experimental_design_and_resource_overview.pdf"
  )

  grDevices::cairo_pdf(
    pdf_file,
    width = width,
    height = height
  )

  draw()
  dev.off()


  tif_file <- file.path(
    outdir,
    paste0(
      "Figure_1_experimental_design_and_",
      "resource_overview_600dpi.tiff"
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
# 5. PANEL A
# 2026 PRIMARY EXPERIMENTAL DESIGN
# ============================================================

boxes_A <- data.frame(
  xmin = c(
    0.35,
    0.04,
    0.38,
    0.70
  ),
  xmax = c(
    0.65,
    0.30,
    0.64,
    0.96
  ),
  ymin = c(
    0.72,
    0.40,
    0.40,
    0.40
  ),
  ymax = c(
    0.88,
    0.60,
    0.60,
    0.60
  ),
  fill = c(
    "#F1F1F1",
    light_blue,
    light_orange,
    light_green
  ),
  label = c(
    "Experimental animals\nGroup allocation",

    paste0(
      "Baseline controls\n",
      "Uninjured\n",
      "n = ",
      n26_baseline
    ),

    paste0(
      "Injury cohort\n",
      "70-cm gravity-hammer\n",
      "skeletal muscle contusion"
    ),

    paste0(
      "Post-injury sampling\n",
      "0–72 h\n",
      "n = ",
      n26_injury
    )
  ),
  stringsAsFactors = FALSE
)

arrows_A <- data.frame(
  x = c(
    0.44,
    0.56,
    0.64
  ),
  y = c(
    0.72,
    0.72,
    0.50
  ),
  xend = c(
    0.17,
    0.51,
    0.70
  ),
  yend = c(
    0.60,
    0.60,
    0.50
  )
)

p1a <- ggplot() +

  geom_rect(
    data = boxes_A,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    colour = "grey35",
    linewidth = 0.45,
    show.legend = FALSE
  ) +

  scale_fill_identity() +

  geom_segment(
    data = arrows_A,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    arrow = arrow(
      length = unit(
        0.17,
        "cm"
      )
    ),
    linewidth = 0.7,
    colour = dark
  ) +

  annotate(
    "text",
    x = (
      boxes_A$xmin +
      boxes_A$xmax
    ) / 2,
    y = (
      boxes_A$ymin +
      boxes_A$ymax
    ) / 2,
    label = boxes_A$label,
    size = 3.25,
    lineheight = 1.06
  ) +

  annotate(
    "rect",
    xmin = 0.26,
    xmax = 0.74,
    ymin = 0.17,
    ymax = 0.28,
    fill = "#F5F5F5",
    colour = "grey65",
    linewidth = 0.35
  ) +

  annotate(
    "text",
    x = 0.50,
    y = 0.225,
    label = paste0(
      "2026 primary data resource  •  ",
      "77 RNA-seq samples"
    ),
    size = 3.4,
    fontface = "bold"
  ) +

  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0.12, 0.94),
    clip = "off"
  ) +

  labs(
    tag = "A",
    title = "Experimental design of the 2026 primary cohort"
  ) +

  theme_void_pub


# ============================================================
# 6. PANEL B
# EQUALLY SPACED SAMPLING-STAGE DESIGN
# ============================================================

timeline26 <- sum26

timeline26$x <- seq_len(
  nrow(timeline26)
)

timeline26$Type <- ifelse(
  timeline26$Group == "Baseline",
  "Baseline",
  "Post-injury"
)

timeline26$Label <- paste0(
  timeline26$Stage,
  "\n",
  "n=",
  timeline26$N
)

line_B <- data.frame(
  x = 1,
  xend = nrow(timeline26),
  y = 0,
  yend = 0
)

p1b <- ggplot(
  timeline26,
  aes(
    x = x,
    y = 0
  )
) +

  geom_segment(
    data = line_B,
    aes(
      x = x,
      xend = xend,
      y = y,
      yend = yend
    ),
    inherit.aes = FALSE,
    linewidth = 0.75,
    colour = "grey40"
  ) +

  geom_point(
    aes(
      fill = Type
    ),
    shape = 21,
    size = 5.4,
    stroke = 0.5,
    colour = "grey30"
  ) +

  geom_text(
    aes(
      y = 0.19,
      label = Label
    ),
    size = 3.0,
    lineheight = 1.05
  ) +

  geom_segment(
    x = 1.5,
    xend = 1.5,
    y = -0.17,
    yend = 0.18,
    linetype = "dashed",
    linewidth = 0.55,
    colour = orange
  ) +

  annotate(
    "text",
    x = 1.5,
    y = -0.26,
    label = "Injury",
    size = 3.0,
    fontface = "bold"
  ) +

  annotate(
    "text",
    x = 6.2,
    y = -0.27,
    label = paste0(
      "10 post-injury stages • ",
      "70 injury samples"
    ),
    size = 3.0
  ) +

  annotate(
    "text",
    x = 6.2,
    y = -0.39,
    label = paste0(
      "Sampling-stage positions are ",
      "equally spaced for readability"
    ),
    size = 2.7,
    fontface = "italic",
    colour = "grey40"
  ) +

  scale_fill_manual(
    values = c(
      "Baseline" = "#F2C94C",
      "Post-injury" = blue
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    limits = c(
      0.55,
      11.45
    ),
    expand = c(0, 0)
  ) +

  coord_cartesian(
    ylim = c(
      -0.46,
      0.38
    ),
    clip = "off"
  ) +

  labs(
    tag = "B",
    title = "2026 primary-cohort sampling design"
  ) +

  theme_void_pub


# ============================================================
# 7. PANEL C
# ANALYSIS RESOURCE WORKFLOW
# ============================================================

boxes_C <- data.frame(
  xmin = c(
    0.28,
    0.28,
    0.02,
    0.35,
    0.68,
    0.02,
    0.52
  ),

  xmax = c(
    0.72,
    0.72,
    0.31,
    0.64,
    0.98,
    0.44,
    0.98
  ),

  ymin = c(
    0.84,
    0.66,
    0.39,
    0.39,
    0.39,
    0.12,
    0.12
  ),

  ymax = c(
    0.96,
    0.78,
    0.56,
    0.56,
    0.56,
    0.27,
    0.27
  ),

  fill = c(
    light_green,
    light_blue,
    light_blue,
    light_orange,
    light_purple,
    "#F1F1F1",
    light_orange
  ),

  label = c(
    "2026 primary RNA-seq resource\n77 samples",

    paste0(
      "Raw counts + TMM.TPM\n",
      "+ sample metadata"
    ),

    paste0(
      "Count filtering + DESeq2 VST\n",
      "PCA • correlation\n",
      "PERMANOVA / PERMDISP"
    ),

    paste0(
      "Filtered raw counts\n",
      "timepoint-vs-Baseline DESeq2\n",
      "+ injury-only spline LRT"
    ),

    paste0(
      "Filtered TMM.TPM\n",
      "group medians\n",
      "Tau • SPM • peak stage"
    ),

    paste0(
      "2021 historical-reference cohort\n",
      "independent reanalysis"
    ),

    paste0(
      "Baseline-relative cross-cohort concordance\n",
      "+ sample-exclusion sensitivity"
    )
  ),

  stringsAsFactors = FALSE
)


# arrows: top -> input
# input -> 3 branches
# historical reference -> cross-cohort
# DE branch -> cross-cohort

arrows_C <- data.frame(
  x = c(
    0.50,
    0.50,
    0.50,
    0.50,
    0.23,
    0.64
  ),

  y = c(
    0.84,
    0.66,
    0.66,
    0.66,
    0.195,
    0.39
  ),

  xend = c(
    0.50,
    0.165,
    0.495,
    0.83,
    0.52,
    0.75
  ),

  yend = c(
    0.78,
    0.56,
    0.56,
    0.56,
    0.195,
    0.27
  )
)

p1c <- ggplot() +

  geom_rect(
    data = boxes_C,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = fill
    ),
    colour = "grey40",
    linewidth = 0.4,
    show.legend = FALSE
  ) +

  scale_fill_identity() +

  geom_segment(
    data = arrows_C,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    arrow = arrow(
      length = unit(
        0.14,
        "cm"
      )
    ),
    colour = "grey35",
    linewidth = 0.6
  ) +

  annotate(
    "text",
    x = (
      boxes_C$xmin +
      boxes_C$xmax
    ) / 2,
    y = (
      boxes_C$ymin +
      boxes_C$ymax
    ) / 2,
    label = boxes_C$label,
    size = 2.75,
    lineheight = 1.03
  ) +

  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(
      0.06,
      1
    ),
    clip = "off"
  ) +

  labs(
    tag = "C",
    title = "Data processing and validation framework"
  ) +

  theme_void_pub


# ============================================================
# 8. PANEL D
# COHORT DESIGN COMPARISON TABLE
# ============================================================

comparison <- data.frame(
  Feature = c(
    "Role",
    "Injury paradigm",
    "Total samples",
    "Baseline samples",
    "Post-injury samples",
    "Post-injury stages",
    "Temporal coverage"
  ),

  Cohort2026 = c(
    "Primary data resource",
    "70-cm gravity-hammer",
    "77",
    as.character(n26_baseline),
    as.character(n26_injury),
    "0, 1, 2, 6, 10, 14, 18, 36, 60, 72 h",
    "0–72 h"
  ),

  Cohort2021 = c(
    "Previously published\nhistorical reference",
    "50-cm injury paradigm",
    "60",
    as.character(n21_baseline),
    as.character(n21_injury),
    "4, 8, 12, 16, 20, 24, 48 h",
    "4–48 h"
  ),

  stringsAsFactors = FALSE
)

comparison$y <- rev(
  seq_len(
    nrow(comparison)
  )
)

row_fill <- rep(
  c(
    "#FFFFFF",
    "#F7F7F7"
  ),
  length.out = nrow(comparison)
)

p1d <- ggplot() +

  # row backgrounds
  geom_rect(
    data = comparison,
    aes(
      xmin = 0,
      xmax = 1,
      ymin = y - 0.43,
      ymax = y + 0.43
    ),
    fill = row_fill,
    colour = NA
  ) +

  # header boxes
  annotate(
    "rect",
    xmin = 0.00,
    xmax = 0.28,
    ymin = 7.60,
    ymax = 8.35,
    fill = "#EFEFEF",
    colour = "grey65"
  ) +

  annotate(
    "rect",
    xmin = 0.28,
    xmax = 0.64,
    ymin = 7.60,
    ymax = 8.35,
    fill = light_blue,
    colour = "grey65"
  ) +

  annotate(
    "rect",
    xmin = 0.64,
    xmax = 1.00,
    ymin = 7.60,
    ymax = 8.35,
    fill = light_orange,
    colour = "grey65"
  ) +

  annotate(
    "text",
    x = 0.14,
    y = 7.98,
    label = "Feature",
    fontface = "bold",
    size = 3.4
  ) +

  annotate(
    "text",
    x = 0.46,
    y = 7.98,
    label = "2026 primary cohort",
    fontface = "bold",
    size = 3.4
  ) +

  annotate(
    "text",
    x = 0.82,
    y = 7.98,
    label = "2021 historical-reference cohort",
    fontface = "bold",
    size = 3.2
  ) +

  geom_text(
    data = comparison,
    aes(
      x = 0.02,
      y = y,
      label = Feature
    ),
    hjust = 0,
    size = 3.0,
    fontface = "bold"
  ) +

  geom_text(
    data = comparison,
    aes(
      x = 0.46,
      y = y,
      label = Cohort2026
    ),
    size = 2.85,
    lineheight = 1.0
  ) +

  geom_text(
    data = comparison,
    aes(
      x = 0.82,
      y = y,
      label = Cohort2021
    ),
    size = 2.85,
    lineheight = 1.0
  ) +

  geom_segment(
    x = 0.28,
    xend = 0.28,
    y = 0.55,
    yend = 8.35,
    colour = "grey78",
    linewidth = 0.35
  ) +

  geom_segment(
    x = 0.64,
    xend = 0.64,
    y = 0.55,
    yend = 8.35,
    colour = "grey78",
    linewidth = 0.35
  ) +

  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(
      0.5,
      8.4
    ),
    clip = "off"
  ) +

  labs(
    tag = "D",
    title = "Primary and historical-reference cohort design"
  ) +

  theme_void_pub


# ============================================================
# 9. SAVE INDIVIDUAL PANELS
# ============================================================

save_panel(
  p1a,
  "Figure_1A_experimental_design_v2",
  width = 7,
  height = 4.6
)

save_panel(
  p1b,
  "Figure_1B_sampling_design_v2",
  width = 7,
  height = 4.6
)

save_panel(
  p1c,
  "Figure_1C_analysis_framework_v2",
  width = 7,
  height = 5.3
)

save_panel(
  p1d,
  "Figure_1D_cohort_comparison_v2",
  width = 7,
  height = 5.3
)


# ============================================================
# 10. SAVE FINAL FIGURE 1
# ============================================================

save_combined(
  p1a,
  p1b,
  p1c,
  p1d,
  width = 13,
  height = 10
)


# ============================================================
# 11. SAVE OBJECTS
# ============================================================

saveRDS(
  list(
    Figure1A = p1a,
    Figure1B = p1b,
    Figure1C = p1c,
    Figure1D = p1d
  ),
  file = file.path(
    outdir,
    "Figure_1_all_plot_objects_v2.rds"
  )
)

write.table(
  sum26,
  file = file.path(
    panel_dir,
    "Figure_1_2026_sample_summary_v2.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  sum21,
  file = file.path(
    panel_dir,
    "Figure_1_2021_sample_summary_v2.tsv"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# 12. FINAL REPORT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 1 v2 COMPLETED\n")
cat("============================================================\n")

cat(
  "2026 primary cohort:\n"
)

cat(
  "  Baseline       : ",
  n26_baseline,
  "\n",
  sep = ""
)

cat(
  "  Injury samples : ",
  n26_injury,
  "\n",
  sep = ""
)

cat(
  "  Total          : ",
  sum(sum26$N),
  "\n\n",
  sep = ""
)

cat(
  "2021 historical-reference cohort:\n"
)

cat(
  "  Baseline       : ",
  n21_baseline,
  "\n",
  sep = ""
)

cat(
  "  Injury samples : ",
  n21_injury,
  "\n",
  sep = ""
)

cat(
  "  Total          : ",
  sum(sum21$N),
  "\n\n",
  sep = ""
)

cat(
  "2021 is NOT represented as deriving from ",
  "the 70-cm 2026 injury experiment.\n"
)

cat(
  "Panel B uses equally spaced sampling-stage positions ",
  "for readability.\n"
)

cat(
  "No proportional elapsed-time interpretation ",
  "should be applied to Panel B spacing.\n\n"
)

cat(
  "Output:\n",
  file.path(
    outdir,
    "Figure_1_experimental_design_and_resource_overview.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  file.path(
    outdir,
    paste0(
      "Figure_1_experimental_design_and_",
      "resource_overview_600dpi.tiff"
    )
  ),
  "\n\n",
  sep = ""
)

cat("PASS\n")
