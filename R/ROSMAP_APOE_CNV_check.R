################################################################################
# APOE Locus CNV Overlap Check
################################################################################
#
# Purpose: Verify that no CNVs in the ROSMAP dataset overlap the APOE gene or
# its immediate flanking region. APOE genotype (ε2/ε3/ε4) is the strongest
# known genetic risk factor for late-onset AD and is driven by two coding SNPs
# (rs429358, rs7412) rather than copy number variation. Nevertheless, we
# explicitly checked whether any CNVs in this dataset overlap the APOE locus
# to rule out potential confounding.
#
# APOE genomic coordinates (GRCh38/hg38):
#   Chromosome 19: 44,905,754 – 44,909,393
#   Ensembl gene ID: ENSG00000130203
#   NCBI Gene ID: 348  |  HGNC: 613
#
# References for APOE locus coordinates:
#   Ensembl Release 112 (GRCh38.p14):
#     https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130203
#   NCBI Gene entry 348:
#     https://www.ncbi.nlm.nih.gov/gene/348
#
# Prerequisites:
#   data/raw/CNVcalls.tsv (produced upstream; CHROM is numeric, no "chr" prefix)
#
################################################################################

# Set working directory to project root
if (!file.exists("data/raw/CNVcalls.tsv")) {
  if (file.exists("../data/raw/CNVcalls.tsv")) {
    setwd("..")
  } else {
    stop("Cannot find data/raw/CNVcalls.tsv. Run from project root or R/ directory.")
  }
}


# ---------------------------------------------------------------------------
# APOE locus definition (GRCh38/hg38)
#   Ensembl ENSG00000130203 | NCBI Gene 348
# ---------------------------------------------------------------------------
APOE_CHR   <- 19L
APOE_START <- 44905754L
APOE_END   <- 44909393L
APOE_GENE  <- "APOE"

# Flanking window used for a secondary "wide" check (500 kb each side)
FLANK_BP <- 500000L

# ---------------------------------------------------------------------------
# Load CNV calls
# ---------------------------------------------------------------------------
message("Loading CNV calls...")
cnv <- read.csv("data/raw/CNVcalls.tsv", sep = "\t", header = TRUE,
                stringsAsFactors = FALSE)

# Columns: CHROM START END NAME SVTYPE QUAL <sample1> <sample2> ...
stopifnot(all(c("CHROM", "START", "END", "NAME", "SVTYPE") %in% names(cnv)))

meta_cols  <- c("CHROM", "START", "END", "NAME", "SVTYPE", "QUAL")
sample_ids <- setdiff(names(cnv), meta_cols)
n_cnvs     <- nrow(cnv)
n_samples  <- length(sample_ids)

message(sprintf("  %d CNV loci across %d samples", n_cnvs, n_samples))

# ---------------------------------------------------------------------------
# Overlap helper
# Two intervals [a1,a2] and [b1,b2] overlap if a1 <= b2 AND a2 >= b1
# ---------------------------------------------------------------------------
overlaps_apoe <- function(chrom, start, end, flank = 0L) {
  chrom == APOE_CHR &
    start <= (APOE_END   + flank) &
    end   >= (APOE_START - flank)
}

# ---------------------------------------------------------------------------
# Check 1: exact overlap with APOE gene body
# ---------------------------------------------------------------------------
cnv_apoe_exact <- cnv[overlaps_apoe(cnv$CHROM, cnv$START, cnv$END, flank = 0L), ]
n_exact <- nrow(cnv_apoe_exact)

# ---------------------------------------------------------------------------
# Check 2: overlap within 500 kb flanking window
# ---------------------------------------------------------------------------
cnv_apoe_wide <- cnv[overlaps_apoe(cnv$CHROM, cnv$START, cnv$END, flank = FLANK_BP), ]
n_wide <- nrow(cnv_apoe_wide)

# ---------------------------------------------------------------------------
# For any overlapping CNVs, count how many samples carry a non-reference CN
# (copy number != 2)
# ---------------------------------------------------------------------------
count_carriers <- function(cnv_subset) {
  if (nrow(cnv_subset) == 0L) return(data.frame())
  sample_dat <- cnv_subset[, sample_ids, drop = FALSE]
  carriers <- apply(sample_dat != 2, 1, sum, na.rm = TRUE)
  cbind(cnv_subset[, meta_cols], n_carriers = carriers)
}

carriers_exact <- count_carriers(cnv_apoe_exact)
carriers_wide  <- count_carriers(cnv_apoe_wide)

# ---------------------------------------------------------------------------
# Nearest CNV on chr19 (for context)
# ---------------------------------------------------------------------------
chr19 <- cnv[cnv$CHROM == APOE_CHR, ]
if (nrow(chr19) > 0) {
  chr19$dist_to_apoe <- pmax(0L,
    pmax(chr19$START - APOE_END, APOE_START - chr19$END))
  chr19_nearest <- chr19[order(chr19$dist_to_apoe), ]
  nearest_5 <- chr19_nearest[1:min(5, nrow(chr19_nearest)),
                              c(meta_cols, "dist_to_apoe")]
} else {
  nearest_5 <- data.frame(message = "No CNVs on chromosome 19")
}

# ---------------------------------------------------------------------------
# Summary report
# ---------------------------------------------------------------------------
cat("================================================================================\n")
cat("APOE Locus CNV Overlap Check — ROSMAP Dataset\n")
cat("================================================================================\n\n")

cat("APOE GENE COORDINATES (GRCh38/hg38)\n")
cat("-------------------------------------\n")
cat(sprintf("  Gene:        %s\n", APOE_GENE))
cat(sprintf("  Chromosome:  %d\n", APOE_CHR))
cat(sprintf("  Start:       %s\n", format(APOE_START, big.mark = ",")))
cat(sprintf("  End:         %s\n", format(APOE_END,   big.mark = ",")))
cat(sprintf("  Length:      %s bp\n",
            format(APOE_END - APOE_START + 1L, big.mark = ",")))
cat("  Source:      Ensembl Release 112 (ENSG00000130203); NCBI Gene 348\n\n")

cat("DATASET SUMMARY\n")
cat("---------------\n")
cat(sprintf("  Total CNV loci:  %d\n", n_cnvs))
cat(sprintf("  Total samples:   %d\n", n_samples))
cat(sprintf("  Chr19 loci:      %d\n", nrow(chr19)), "\n")

cat("OVERLAP RESULTS\n")
cat("---------------\n")
cat(sprintf("  CNVs overlapping APOE gene body (exact):           %d\n", n_exact))
cat(sprintf("  CNVs within +/- %s bp of APOE:  %d\n",
            format(FLANK_BP, big.mark = ","), n_wide))

if (n_exact == 0L && n_wide == 0L) {
  cat("\n  RESULT: No CNVs overlap the APOE gene or its 500 kb flanking region.\n")
  cat("  The nearest CNV on chromosome 19 is located approximately\n")
  cat(sprintf("  %s bp (~%.1f Mb) from the APOE locus.\n",
              format(nearest_5$dist_to_apoe[1], big.mark = ","),
              nearest_5$dist_to_apoe[1] / 1e6))
  cat("\n  CONCLUSION: CNV-based confounding of APOE genotype effects is not\n")
  cat("  a concern in this dataset. No individual in the cohort carries a CNV\n")
  cat("  overlapping or proximal to the APOE locus.\n")
} else if (n_exact > 0L) {
  cat(sprintf("\n  WARNING: %d CNV(s) overlap the APOE gene body.\n", n_exact))
  print(carriers_exact)
} else {
  cat(sprintf("\n  NOTE: No CNVs overlap APOE directly, but %d CNV(s) fall within\n", n_wide))
  cat("  the 500 kb flanking window.\n")
  print(carriers_wide)
}

cat("\nNEAREST CNVs ON CHROMOSOME 19\n")
cat("------------------------------\n")
print(nearest_5)

cat("\nREFERENCES\n")
cat("----------\n")
cat("  APOE coordinates: Ensembl Release 112, GRCh38.p14 (ENSG00000130203)\n")
cat("    https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000130203\n")
cat("  APOE gene entry:  NCBI Gene ID 348\n")
cat("    https://www.ncbi.nlm.nih.gov/gene/348\n")
cat("  APOE & AD risk:   Corder et al. (1993) Science 261:921-923.\n")
cat("                    DOI: 10.1126/science.8346443\n")
cat("  APOE SNP basis:   rs429358 (Cys112Arg) and rs7412 (Arg158Cys)\n")
cat("    determine ε2/ε3/ε4 alleles — point mutations, not CNVs.\n")
cat("================================================================================\n")
