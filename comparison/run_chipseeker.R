#!/usr/bin/env Rscript
# Run ChIPseeker on the bundled comparison GTF/BED.
#
# Usage:
#   Rscript run_chipseeker.R <gtf> <bed> <out_prefix> [level] [tssRegion] [overlap]
#
# Example:
#   Rscript run_chipseeker.R example.gtf example_peaks.bed r_default transcript "-3000,3000" TSS
#
# Outputs:
#   <out_prefix>.tsv       one row per peak (ChIPseeker as.data.frame columns)
#   <out_prefix>.annoStat.tsv
#   <out_prefix>.detail.tsv
#   <out_prefix>.sessionInfo.txt

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    library(txdbmaker)
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run_chipseeker.R <gtf> <bed> <out_prefix> [level] [tssRegion] [overlap]")
}

gtf <- args[[1]]
bed <- args[[2]]
out_prefix <- args[[3]]
level <- if (length(args) >= 4) args[[4]] else "transcript"
tss_text <- if (length(args) >= 5) args[[5]] else "-3000,3000"
overlap <- if (length(args) >= 6) args[[6]] else "TSS"

tss_region <- as.integer(strsplit(tss_text, ",", fixed = TRUE)[[1]])
if (length(tss_region) != 2L) {
  stop("tssRegion must be two comma-separated integers")
}

# ChIPseeker's package startup option (300 bp) is the actual default; make it
# explicit so an interactive R profile cannot change the comparison.
options(ChIPseeker.downstreamDistance = 300)

message("ChIPseeker ", as.character(packageVersion("ChIPseeker")))
message("GTF: ", gtf)
message("BED: ", bed)
message("level=", level, " tssRegion=", paste(tss_region, collapse = ","),
        " overlap=", overlap)

txdb <- makeTxDbFromGFF(
  gtf,
  format = "gtf",
  dataSource = "pychipseeker comparison",
  organism = "Homo sapiens"
)

peaks <- readPeakFile(bed)
message("peaks: ", length(peaks))

anno <- annotatePeak(
  peaks,
  TxDb = txdb,
  tssRegion = tss_region,
  level = level,
  overlap = overlap,
  verbose = FALSE
)

df <- as.data.frame(anno)
# ChIPseeker 1.48.0 converts the GRanges factor columns geneChr/geneStrand to
# integer codes in as.data.frame().  Restore their labels for comparison.
if (is.factor(anno@anno$geneChr)) {
  df$geneChr <- as.character(anno@anno$geneChr)
}
if (is.factor(anno@anno$geneStrand)) {
  df$geneStrand <- as.character(anno@anno$geneStrand)
}
df$peakId <- df$V4
df$peakStart0 <- df$start - 1L
df$geneStart0 <- df$geneStart - 1L
write.table(
  df,
  paste0(out_prefix, ".tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

stat <- as.data.frame(anno@annoStat)
write.table(
  stat,
  paste0(out_prefix, ".annoStat.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

detail <- as.data.frame(anno@detailGenomicAnnotation)
detail <- cbind(peakId = df$peakId, detail)
write.table(
  detail,
  paste0(out_prefix, ".detail.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

writeLines(
  capture.output(sessionInfo()),
  paste0(out_prefix, ".sessionInfo.txt")
)

message("wrote ", paste0(out_prefix, ".tsv"))
