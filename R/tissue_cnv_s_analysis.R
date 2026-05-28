# CNV-S by sequencing source tissue: reproduces Supp Table tissue effects using
# Pass QC rows in WGS_sample_QC_info (WGS_id matches CNVcalls column names).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

root <- "/Users/u1082736/Documents/Personal/cnvrs-study"
setwd(root)

if (!file.exists("data/pLI_LOEUF_data.rds")) {
  source("R/CNVRS_preprocess.R")
}

pli_loeuf <- readRDS("data/pLI_LOEUF_data.rds")
pli_loeuf$oe_lof_upper <- pli_loeuf$oe_lof_upper^(-1)

haplo_triplo <- readRDS("data/pHaplo_pTriplo_data.rds")
colnames(haplo_triplo) <- c("gene", "pHaplo", "pTriplo")
pHaploLst <- haplo_triplo
pHaploLst$pTriplo <- NULL
pTriploLst <- haplo_triplo
pTriploLst$pHaplo <- NULL

thresh1 <- 0.86
pHaplo_thresh <- pHaploLst
pHaplo_thresh[which(pHaplo_thresh$pHaplo < thresh1), 2] <- 0
pHaplo_thresh[which(pHaplo_thresh$pHaplo >= thresh1), 2] <- 1
thresh2 <- 0.94
pTriplo_thresh <- pTriploLst
pTriplo_thresh[which(pTriplo_thresh$pTriplo < thresh2), 2] <- 0
pTriplo_thresh[which(pTriplo_thresh$pTriplo >= thresh2), 2] <- 1

trnsfm <- function(count) count

CNVcalls <- read.delim(
  "data/raw/CNVcalls.tsv",
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
)

VEP1 <- read.delim("data/raw/01_04_2023_VEP_output_filtered.txt", stringsAsFactors = FALSE,
                  comment.char = "#")
VEP <- VEP1[, c(2, 3, 6)]
VEP <- VEP[!duplicated(VEP), ]
colnames(VEP) <- c("location", "SVTYPE", "gene")
VEP <- VEP[!grepl("-", VEP$gene), ]
VEP2 <- tidyr::separate(VEP, col = location, into = c("CHROM", "loc"), sep = ":", remove = FALSE)
VEP2 <- tidyr::separate(VEP2, col = loc, into = c("START", "END"), sep = "-")
VEP2[VEP2 == "deletion"] <- "DEL"
VEP2[VEP2 == "duplication"] <- "DUP"
VEP2$CHROM <- as.integer(VEP2$CHROM)
VEP2$START <- as.integer(VEP2$START)
VEP2$END <- as.integer(VEP2$END)
VEP <- VEP2[, c("CHROM", "START", "END", "SVTYPE", "gene")]

geneContentTable <- crossing(CNVcalls[, c(1, 2, 3, 5)], VEP, .name_repair = "minimal")
colnames(geneContentTable) <- c(
  "CHROM.x", "START.x", "END.x", "SVTYPE.x",
  "CHROM.y", "START.y", "END.y", "SVTYPE.y", "gene"
)
geneContentTable <- filter(geneContentTable, SVTYPE.x == SVTYPE.y, CHROM.x == CHROM.y,
                          START.x - 1 <= START.y, END.x + 1 >= END.y)
geneContentTable <- geneContentTable[, c("CHROM.x", "START.x", "END.x", "SVTYPE.x", "gene")]
colnames(geneContentTable) <- c("CHROM", "START", "END", "SVTYPE", "gene")
geneContentTable <- unique(geneContentTable)

wgs_qc <- read.csv("data/cnv_meta_data/WGS_sample_QC_info.csv", stringsAsFactors = FALSE)
wgs_qc <- wgs_qc[wgs_qc$QC == "Pass", ]

read_depth <- read.delim("data/cnv_meta_data/TableS2.ReadDepth.GenomeMean.txt",
                        stringsAsFactors = FALSE)
colnames(read_depth) <- c("sample_id", "read_depth")

# Match manuscript-style coarse labels (supplementary table caption).
coarse_tissue <- function(src) {
  s <- tolower(trimws(as.character(src)))
  if (is.na(s) || !nzchar(s)) return(NA_character_)
  if (grepl("lymphocytes|ebv", s)) return("unspecified")
  if (grepl("whole blood|blood-pbmc|^blood$|blood-cerebellum", s)) return("blood")
  if (grepl("brain|cingulate|dlpfc|cortex|caudate|cerebellum|ba[0-9]|unknown", s)) {
    return("brain")
  }
  return("other")
}

lookup_rd <- function(wgs_id) {
  hit <- read_depth[tolower(read_depth$sample_id) == tolower(wgs_id), "read_depth"]
  if (length(hit)) return(as.numeric(hit[1]))
  NA_real_
}

one_sample_scores <- function(mat_row_indices, count_vals) {
  lstCNV <- data.frame(
    CHROM = CNVcalls$CHROM[mat_row_indices],
    START = CNVcalls$START[mat_row_indices],
    END = CNVcalls$END[mat_row_indices],
    SVTYPE = CNVcalls$SVTYPE[mat_row_indices],
    count = count_vals,
    stringsAsFactors = FALSE
  )
  lstCNV <- lstCNV[lstCNV$SVTYPE != "mCNV", , drop = FALSE]
  if (!nrow(lstCNV)) {
    return(rep(NA_real_, 8))
  }

  lstCNV[lstCNV$SVTYPE == "DUP", "count"] <- lstCNV[lstCNV$SVTYPE == "DUP", "count"] - 2
  lstCNV[lstCNV$SVTYPE == "DEL", "count"] <- 2 - lstCNV[lstCNV$SVTYPE == "DEL", "count"]
  lstCNV$CNVsizes <- lstCNV$END - lstCNV$START
  lstCNV$count <- trnsfm(lstCNV$count)

  genes <- merge(lstCNV, geneContentTable, by = c("CHROM", "START", "END", "SVTYPE"))
  sum_gcount_del <- sum(genes[genes$SVTYPE == "DEL", "count"])
  sum_gcount_dup <- sum(genes[genes$SVTYPE == "DUP", "count"])

  pli <- unique(merge(genes, pli_loeuf, by = "gene"))
  pli$w_pli <- pli$count * pli$pLI
  pli$w_loeuf <- pli$count * pli$oe_lof_upper
  sum_pli_del <- sum(pli[pli$SVTYPE == "DEL", "w_pli"], na.rm = TRUE)
  sum_pli_dup <- sum(pli[pli$SVTYPE == "DUP", "w_pli"], na.rm = TRUE)
  sum_loeuf_del <- sum(pli[pli$SVTYPE == "DEL", "w_loeuf"], na.rm = TRUE)
  sum_loeuf_dup <- sum(pli[pli$SVTYPE == "DUP", "w_loeuf"], na.rm = TRUE)

  phi <- unique(merge(genes, pHaploLst, by = "gene"))
  phi <- phi[phi$SVTYPE == "DEL", ]
  pts <- unique(merge(genes, pTriploLst, by = "gene"))
  pts <- pts[pts$SVTYPE == "DUP", ]
  sum_phi <- sum(phi$count * phi$pHaplo)
  sum_pts <- sum(pts$count * pts$pTriplo)

  phi_thresh <- unique(merge(genes, pHaplo_thresh, by = "gene"))
  pts_thresh <- unique(merge(genes, pTriplo_thresh, by = "gene"))
  phi_thresh <- phi_thresh[phi_thresh$SVTYPE == "DEL", ]
  pts_thresh <- pts_thresh[pts_thresh$SVTYPE == "DUP", ]
  sum_phi_thresh <- sum(phi_thresh$count * phi_thresh$pHaplo)
  sum_pts_thresh <- sum(pts_thresh$count * pts_thresh$pTriplo)

  c(pli_del = sum_pli_del, pli_dup = sum_pli_dup,
    loeuf_del = sum_loeuf_del, loeuf_dup = sum_loeuf_dup,
    pHI = sum_phi, pTS = sum_pts,
    pHI_bin = sum_phi_thresh, pTS_bin = sum_pts_thresh)
}

cn <- colnames(CNVcalls)
sample_cols <- 7:ncol(CNVcalls)

del_mask <- CNVcalls$SVTYPE == "DEL"
dup_mask <- CNVcalls$SVTYPE == "DUP"
ndel <- sapply(sample_cols, function(i) sum(2 - CNVcalls[del_mask, i], na.rm = TRUE))
ndup <- sapply(sample_cols, function(i) sum(CNVcalls[dup_mask, i] - 2, na.rm = TRUE))
bad <- ndel > 6000 | ndup > 1000

results <- list()
for (k in seq_along(sample_cols)) {
  if (bad[k]) next
  i <- sample_cols[k]
  smid <- cn[i]
  wgs_hyphen <- if (grepl("^MAP|^ROS", smid)) smid else gsub(".", "-", smid, fixed = TRUE)
  tr <- wgs_qc[tolower(wgs_qc$WGS_id) == tolower(wgs_hyphen), , drop = FALSE]
  if (!nrow(tr)) next
  grp <- coarse_tissue(tr$Source.Tissue.Type[1])
  rd <- lookup_rd(wgs_hyphen)
  if (is.na(rd)) rd <- as.numeric(tr$GQN[1])

  ok <- CNVcalls$SVTYPE != "mCNV"
  cnvec <- unlist(CNVcalls[ok, i, drop = TRUE])
  sc <- one_sample_scores(which(ok), cnvec)

  results[[length(results) + 1]] <- data.frame(
    wgs_id = wgs_hyphen,
    projid = tr$projid[1],
    source_tissue = grp,
    read_depth = rd,
    t(sc),
    stringsAsFactors = FALSE
  )
}

df <- bind_rows(results)
df <- df[df$source_tissue %in% c("brain", "blood", "unspecified"), ]

score_cols <- c(
  "pli_del", "pli_dup", "loeuf_del", "loeuf_dup",
  "pHI", "pTS", "pHI_bin", "pTS_bin"
)

# --- Model A (matches narrative): brain vs blood only ---
df_bb <- df[df$source_tissue %in% c("brain", "blood"), ]
df_bb$tissue_f <- factor(df_bb$source_tissue, levels = c("blood", "brain"))

# --- Model B: 3-level factor like supplement (blood ref) ---
df$source_3 <- factor(df$source_tissue,
                      levels = c("blood", "brain", "unspecified"))

out_path <- file.path(root, "output_tables/tissue_cnv_s_lm_results.txt")
sink(out_path)
cat("Generated by R/tissue_cnv_s_analysis.R\n")
cat("=== Sample counts ===\n")
print(table(df$source_tissue, useNA = "ifany"))
cat("\nBrain vs blood only: n =", nrow(df_bb),
    "(brain:", sum(df_bb$source_tissue == "brain"),
    ", blood:", sum(df_bb$source_tissue == "blood"), ")\n\n")

cat("=== lm(CNV-S ~ brain[vs blood] + read_depth), brain vs blood sample set ===\n")
for (y in score_cols) {
  fit <- lm(as.formula(paste(y, "~ tissue_f + read_depth")), data = df_bb)
  cat("\n", y, ":\n", sep = "")
  print(coef(summary(fit)))
}

cat("\n=== lm(CNV-S ~ source_3 + read_depth) 3-level (blood ref); includes unspecified ===\n")
for (y in score_cols) {
  fit <- lm(as.formula(paste(y, "~ source_3 + read_depth")), data = df)
  cat("\n", y, ":\n", sep = "")
  print(coef(summary(fit)))
}

cat("\n=== Mean (SD) raw CNV-S by coarse tissue ===\n")
agg <- df %>%
  group_by(source_tissue) %>%
  summarise(across(all_of(score_cols), list(m = mean, s = sd), .names = "{.col}_{.fn}"),
            .groups = "drop")
print(agg, width = 200)
sink()

message("Wrote ", out_path)
