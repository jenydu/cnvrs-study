################################################################################
# ROSMAP PGS-LOAD vs CNV Score Correlations
################################################################################
#
# Explores correlations between PGS-LOAD (polygenic score for late-onset
# Alzheimer's disease) and the CNV scores (CNV-S). SNPs in PGS-LOAD can be in LD
# with CNVs in the CNV scores; quantifying these correlations is relevant for
# interpreting independent vs overlapping genetic contributions.
#
# Outputs:
#   - output_tables/PGS_CNV_correlations.csv
#   - output_figs/PGS_CNV_correlation_heatmap.png
#   - output_figs/PGS_CNV_correlation_heatmap.svg
#   - output_figs/PGS_CNV_scatter_plots.png
#
# Prerequisites: Run ROSMAP_CNVRS_pheno_preprocess.R first to create
# data/scoreWithPheno.rds
#
################################################################################

# Set working directory to project root
if (file.exists("data/scoreWithPheno.rds")) {
  # Already in project root
} else if (file.exists("../data/scoreWithPheno.rds")) {
  setwd("..")
} else {
  stop("Cannot find data/scoreWithPheno.rds. Run from project root or R/ directory.")
}

# Load packages
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
# ppcor for partial correlations (optional)
has_ppcor <- requireNamespace("ppcor", quietly = TRUE)
if (!has_ppcor) {
  message("Package 'ppcor' not installed. Partial correlations will be skipped.")
  message("Install with: install.packages('ppcor')")
}

# Create output directories
dir.create("output_tables", showWarnings = FALSE)
dir.create("output_figs", showWarnings = FALSE)

################################################################################
# Load data and merge PGS-LOAD
################################################################################

prs_file <- "data/raw/raw_score_kunkle_full.tsv"
if (!file.exists(prs_file)) {
  stop("PGS-LOAD file not found: ", prs_file)
}

scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

# Load and merge PRS (PGS-LOAD)
PRS <- read.csv(prs_file, sep = '\t', header = TRUE)
ROSmaster <- readRDS("data/ROSmaster.rds")
ROSmaster_keylist <- ROSmaster[c('projid', 'study', 'IID')]
rm(ROSmaster)

ROSmaster_keylist$study <- as.character(ROSmaster_keylist$study)
ROSmaster_keylist$projid <- as.integer(ROSmaster_keylist$projid)
ROSmaster_keylist$study <- gsub('\\s+', '', ROSmaster_keylist$study)

AD <- PRS[which(grepl('AD', PRS$IID)), ]
temp <- merge(AD, ROSmaster_keylist, by.x = "IID", by.y = "IID")
temp$IID <- NULL

ROS <- PRS[-which(grepl('AD', PRS$IID)), ]
study <- substring(ROS$IID, 1, 3)
projid <- as.numeric(substring(ROS$IID, 4))
temp2 <- data.frame(ROS$PRS, projid, study)
colnames(temp2) <- c('PRS', 'projid', 'study')
PRS <- rbind(temp, temp2)

scoreWithPheno <- merge(scoreWithPheno, PRS,
                        by.x = c('projid', 'study.x'),
                        by.y = c('projid', 'study'))
rm(temp, temp2, AD, ROS, ROSmaster_keylist)

# CNV scores (CNV-S)
lst_risk_scores <- c('pli_del', 'pli_dup', 'loeuf_del', 'loeuf_dup',
                     'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

# Display names for plots/tables
score_labels <- c(
  'pli_del'     = 'pLI (DEL)',
  'pli_dup'     = 'pLI (DUP)',
  'loeuf_del'   = 'LOEUF (DEL)',
  'loeuf_dup'   = 'LOEUF (DUP)',
  'pHI'         = 'pHI',
  'pTS'         = 'pTS',
  'pHI_thresh'  = 'pHI (binarized)',
  'pTS_thresh'  = 'pTS (binarized)'
)

################################################################################
# Compute correlations
################################################################################

cor_results <- data.frame(
  CNV_score = character(length(lst_risk_scores)),
  Pearson_r = numeric(length(lst_risk_scores)),
  Pearson_p = character(length(lst_risk_scores)),
  Pearson_p_val = numeric(length(lst_risk_scores)),  # raw p for significance
  Spearman_rho = numeric(length(lst_risk_scores)),
  Spearman_p = character(length(lst_risk_scores)),
  Partial_r = numeric(length(lst_risk_scores)),
  Partial_p = character(length(lst_risk_scores)),
  N = integer(length(lst_risk_scores)),
  stringsAsFactors = FALSE
)

# Covariates for partial correlation (age, sex, 10 eigenvectors)
covar_cols <- c('msex', 'age_death', paste0('egvec', 1:10))
covar_cols <- covar_cols[covar_cols %in% colnames(scoreWithPheno)]

for (i in seq_along(lst_risk_scores)) {
  score <- lst_risk_scores[i]
  valid <- complete.cases(scoreWithPheno[, c('PRS', score)])
  if (length(covar_cols) > 0) {
    valid <- valid & complete.cases(scoreWithPheno[, covar_cols])
  }

  x <- scoreWithPheno$PRS[valid]
  y <- scoreWithPheno[[score]][valid]
  n <- sum(valid)

  # Pearson
  pearson <- cor.test(x, y, method = "pearson", exact = FALSE)
  # Spearman
  spearman <- cor.test(x, y, method = "spearman", exact = FALSE)

  # Partial correlation (PGS-LOAD vs CNV score, adjusting for covariates)
  partial_r <- NA
  partial_p <- NA
  if (has_ppcor && length(covar_cols) > 0 && n > length(covar_cols) + 2) {
    tryCatch({
      Z <- as.matrix(scoreWithPheno[valid, covar_cols])
      pcor_res <- ppcor::pcor.test(x, y, Z, method = "pearson")
      partial_r <- pcor_res$estimate
      partial_p <- pcor_res$p.value
    }, error = function(e) NULL)
  }

  cor_results[i, ] <- c(
    score_labels[score],
    round(pearson$estimate, 4),
    format(pearson$p.value, scientific = TRUE, digits = 3),
    pearson$p.value,
    round(spearman$estimate, 4),
    format(spearman$p.value, scientific = TRUE, digits = 3),
    if (is.na(partial_r)) NA else round(partial_r, 4),
    if (is.na(partial_p)) NA else format(partial_p, scientific = TRUE, digits = 3),
    n
  )
}

# FDR adjustment for Pearson p-values (8 tests)
cor_results$Pearson_p_fdr <- p.adjust(cor_results$Pearson_p_val, method = "fdr")
cor_results$Pearson_p_fdr_fmt <- format(cor_results$Pearson_p_fdr, scientific = TRUE, digits = 3)

# Save correlation table (include raw and FDR-adjusted p, exclude internal columns)
cor_results_out <- cor_results[, c("CNV_score", "Pearson_r", "Pearson_p", "Pearson_p_fdr_fmt",
                                   "Spearman_rho", "Spearman_p", "Partial_r", "Partial_p", "N")]
colnames(cor_results_out)[colnames(cor_results_out) == "Pearson_p_fdr_fmt"] <- "Pearson_p_FDR"
write.csv(cor_results_out, 'output_tables/PGS_CNV_correlations.csv', row.names = FALSE)
message("Saved: output_tables/PGS_CNV_correlations.csv")

################################################################################
# Correlation matrix (PGS-LOAD + all CNV scores)
################################################################################

cor_vars <- c('PRS', lst_risk_scores)
cor_data <- scoreWithPheno[, cor_vars]
cor_data <- cor_data[complete.cases(cor_data), ]

cor_mat <- cor(cor_data, method = "pearson")
n_obs <- nrow(cor_data)

# P-values for Pearson correlations: t = r * sqrt(n-2) / sqrt(1-r^2)
cor_pval <- function(r, n) {
  if (n <= 2 || abs(r) >= 1) return(NA)
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  2 * pt(-abs(t_stat), df = n - 2)
}

# Reshape for ggplot
cor_long <- as.data.frame(cor_mat)
cor_long$Var1 <- rownames(cor_long)
cor_long <- pivot_longer(cor_long, -Var1, names_to = "Var2", values_to = "correlation")
cor_long$pval <- sapply(cor_long$correlation, function(r) cor_pval(r, n_obs))
cor_long$significant <- cor_long$pval < 0.05 & !is.na(cor_long$pval)
cor_long$Var1 <- factor(cor_long$Var1, levels = c('PRS', lst_risk_scores))
cor_long$Var2 <- factor(cor_long$Var2, levels = c('PRS', lst_risk_scores))
levels(cor_long$Var1) <- c('PGS-LOAD', score_labels[lst_risk_scores])
levels(cor_long$Var2) <- c('PGS-LOAD', score_labels[lst_risk_scores])

p_heatmap <- ggplot(cor_long, aes(Var1, Var2, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", correlation), fontface = ifelse(significant, "bold", "plain")), size = 2.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(x = NULL, y = NULL, fill = "Pearson r",
       title = "Correlation: PGS-LOAD and CNV Scores") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('output_figs/PGS_CNV_correlation_heatmap.png', p_heatmap,
       width = 10, height = 8, dpi = 300)
ggsave('output_figs/PGS_CNV_correlation_heatmap.svg', p_heatmap,
       width = 10, height = 8)
message("Saved: output_figs/PGS_CNV_correlation_heatmap.png/.svg")

################################################################################
# Scatter plots: PGS-LOAD vs each CNV score
################################################################################

plot_list <- list()
for (i in seq_along(lst_risk_scores)) {
  score <- lst_risk_scores[i]
  label <- score_labels[score]
  pearson_r <- cor_results$Pearson_r[i]
  pearson_p <- cor_results$Pearson_p[i]
  pearson_p_fdr <- cor_results$Pearson_p_fdr_fmt[i]
  is_sig_fdr <- cor_results$Pearson_p_fdr[i] < 0.05

  p <- ggplot(scoreWithPheno, aes(x = PRS, y = .data[[score]])) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "#B2182B", fill = "#F4A582", alpha = 0.3) +
    labs(x = "PGS-LOAD", y = label,
         title = paste0(label, " vs PGS-LOAD"),
         subtitle = sprintf("r = %.3f, p = %s (FDR = %s)", as.numeric(pearson_r), pearson_p, pearson_p_fdr)) +
    theme_minimal(base_size = 10) +
    theme(plot.subtitle = element_text(face = if (is_sig_fdr) "bold" else "plain"))
  plot_list[[i]] <- p
}

p_combined <- ggarrange(plotlist = plot_list, ncol = 2, nrow = 4)
p_combined <- annotate_figure(p_combined,
  top = text_grob("PGS-LOAD vs CNV Score (CNV-S)", face = "bold", size = 14))

ggsave('output_figs/PGS_CNV_scatter_plots.png', p_combined,
       width = 10, height = 14, dpi = 300)
message("Saved: output_figs/PGS_CNV_scatter_plots.png")

################################################################################
# Done
################################################################################

message("
=============================================================================
PGS-LOAD vs CNV score correlation analysis complete.
Outputs saved to:
  - output_tables/PGS_CNV_correlations.csv
  - output_figs/PGS_CNV_correlation_heatmap.png/.svg
  - output_figs/PGS_CNV_scatter_plots.png
=============================================================================
")
