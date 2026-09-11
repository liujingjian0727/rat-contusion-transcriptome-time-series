options(stringsAsFactors = FALSE)

infile <- "06_2026_continuous_time_spline_LRT_result/06_2026_dynamic_genes_FDR005_with_VST_metrics.tsv"
outdir <- "06_2026_continuous_time_spline_LRT_result"

df <- read.delim(infile, sep="\t", check.names=FALSE)

# ================================
# thresholds
# ================================
RANGE_CUTOFF <- 1

df$High_confidence_dynamic <- df$Mean_VST_temporal_range >= RANGE_CUTOFF

high_df <- df[df$High_confidence_dynamic, ]

write.table(
  high_df,
  file=file.path(outdir, "06B_high_confidence_dynamic_genes.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE
)

write.table(
  data.frame(GeneID=high_df$GeneID),
  file=file.path(outdir, "06B_high_confidence_dynamic_gene_IDs.tsv"),
  sep="\t", row.names=FALSE, quote=FALSE
)

cat("========================================\n")
cat("High-confidence dynamic genes\n")
cat("Range cutoff:", RANGE_CUTOFF, "\n")
cat("Total dynamic genes (LRT):", nrow(df), "\n")
cat("High-confidence dynamic:", nrow(high_df), "\n")
cat("Percent retained:", round(100*nrow(high_df)/nrow(df),2), "%\n")
cat("========================================\n")
