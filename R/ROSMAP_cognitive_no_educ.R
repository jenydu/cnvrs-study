################################################################################
# ROSMAP Cognitive Analysis - Without Years of Education Correction
################################################################################
#
# Re-runs cognitive phenotype association models WITHOUT adjusting for years of
# education (educ). Years of education is genetically correlated with cognitive
# function (https://pmc.ncbi.nlm.nih.gov/articles/PMC6344370) and adjusting for
# it would lead to over-correction.
#
# Cognitive outcomes: Global Cognitive Function (cogn_global), Final Consensus
# Cognitive Diagnosis (cogdx). Neuropathological models are NOT run.
#
# Cognitive models adjusted for: sex (msex), age at cognitive evaluation
# (age_at_visit), and 10 eigenvectors.
#
# Outputs: output_tables_no_educ/ and output_figs_no_educ/
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

# Create output directories
dir.create("output_tables_no_educ", showWarnings = FALSE)
dir.create("output_figs_no_educ", showWarnings = FALSE)

# Load packages
library(stats)
library(reshape2)
library(dplyr)
library(scales)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(ggsci)
library(gridExtra)
library(AICcmodavg)
library(grid)
library(forestploter)
library(lm.beta)
library(pscl)
library(lmtest)

################################################################################
# Load data and define cognitive covariates (NO educ)
################################################################################

scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

# cogdx_binom: cogdx = {1,2} -> 0; cogdx = {4} -> 1; {3,5} -> NA
scoreWithPheno$cogdx_binom <- scoreWithPheno$cogdx
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 2)] <- 1
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 3)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 5)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 1)] <- 0
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 4)] <- 1

str_egvec <- '+ egvec1 + egvec2 + egvec3 + egvec4 + egvec5 + egvec6 + egvec7 + egvec8 + egvec9 + egvec10'
str_egvec_stdz <- '+ scale(egvec1) + scale(egvec2) + scale(egvec3) + scale(egvec4) + scale(egvec5) + scale(egvec6) + scale(egvec7) + scale(egvec8) + scale(egvec9) + scale(egvec10)'

# Cognitive covariates: sex + age_at_visit (NO educ - genetically correlated with cognition)
str_covar_cogn <- '+ msex + age_at_visit'
str_covar_cogn_stdz <- '+ scale(msex) + scale(age_at_visit)'

str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')
str_scores_stdz <- c('~ scale(pli_del) + scale(pli_dup)',
                     '~ scale(loeuf_del) + scale(loeuf_dup)',
                     '~ scale(pHI) + scale(pTS)',
                     '~ scale(pHI_thresh) + scale(pTS_thresh)')

lst_pheno_cog <- c('cogn_global')
lst_risk_scores <- c('pli_del', 'pli_dup', 'loeuf_del', 'loeuf_dup',
                    'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

################################################################################
# 1. MAIN PHENOTYPE ASSOCIATION (cognitive only)
################################################################################

p_t_val_calculation_cogn <- function(scoreWithPheno, lst_risk_scores, str_scores,
                                    str_scores_stdz, str_covar_cogn, str_covar_cogn_stdz,
                                    str_egvec, str_egvec_stdz) {
  cnames <- c('phenotype', 'risk_score', 'standardized_beta', 'p_val', 'ci_lw', 'ci_up', 't_val')
  p_t_val <- matrix(nrow = 0, ncol = 7)
  colnames(p_t_val) <- cnames

  # cogn_global (continuous)
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogn_global', str_scores[j], str_covar_cogn, str_egvec)
    stdz_formula <- paste0('scale(cogn_global)', str_scores_stdz[j],
                           str_covar_cogn_stdz, str_egvec_stdz)
    fit <- lm(formula, data = scoreWithPheno)
    fit_stdz <- lm(stdz_formula, data = scoreWithPheno)
    confint1 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j-1], ')'), level = 0.95)
    confint2 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j], ')'), level = 0.95)
    standardized_betas <- lm.beta(fit_stdz)
    p_t_val <- rbind(p_t_val,
      c('cogn_global', lst_risk_scores[2*j-1],
        unname(standardized_betas$coefficients[2]),
        summary(fit)$coefficients[,4][2],
        confint1[1], confint1[2],
        summary(fit)$coefficients[,3][2]),
      c('cogn_global', lst_risk_scores[2*j],
        unname(standardized_betas$coefficients[3]),
        summary(fit)$coefficients[,4][3],
        confint2[1], confint2[2],
        summary(fit)$coefficients[,3][3]))
  }

  # cogdx (categorical)
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec)
    stdz_formula <- paste0('scale(cogdx_binom)', str_scores_stdz[j],
                           str_covar_cogn_stdz, str_egvec_stdz)
    fit <- glm(formula, data = scoreWithPheno, family = 'binomial')
    fit_stdz <- lm(stdz_formula, data = scoreWithPheno)
    confint1 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j-1], ')'), level = 0.95)
    confint2 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j], ')'), level = 0.95)
    standardized_betas <- lm.beta(fit_stdz)
    p_t_val <- rbind(p_t_val,
      c('cogdx', lst_risk_scores[2*j-1],
        unname(standardized_betas$coefficients[2]),
        summary(fit)$coefficients[,4][2],
        confint1[1], confint1[2],
        summary(fit)$coefficients[,3][2]),
      c('cogdx', lst_risk_scores[2*j],
        unname(standardized_betas$coefficients[3]),
        summary(fit)$coefficients[,4][3],
        confint2[1], confint2[2],
        summary(fit)$coefficients[,3][3]))
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val[c("standardized_beta", "p_val", "t_val", "ci_lw", "ci_up")] <-
    lapply(p_t_val[c("standardized_beta", "p_val", "t_val", "ci_lw", "ci_up")], as.double)
  p_t_val$p_val_adj <- p.adjust(p_t_val$p_val, 'fdr')
  return(p_t_val)
}

p_t_val <- p_t_val_calculation_cogn(scoreWithPheno, lst_risk_scores, str_scores,
                                    str_scores_stdz, str_covar_cogn, str_covar_cogn_stdz,
                                    str_egvec, str_egvec_stdz)

# Rename for display
rename_mappings <- c(
  "pTS_thresh" = "pTS (binarized)",
  "pHI_thresh" = "pHI (binarized)",
  "pli_del" = "pLI (DEL)",
  "pli_dup" = "pLI (DUP)",
  "loeuf_del" = "LOEUF (DEL)",
  "loeuf_dup" = "LOEUF (DUP)",
  "cogn_global" = "Global Cognitive Function\n(n=1011)",
  "cogdx" = "Final Consensus\nCognitive Diagnosis\n(n=946)"
)
for (old_name in names(rename_mappings)) {
  p_t_val[p_t_val == old_name] <- rename_mappings[old_name]
}

write.csv(p_t_val, file = "output_tables_no_educ/p_t_val.csv", row.names = FALSE)
message("Saved: output_tables_no_educ/p_t_val.csv")

# Heatmap
rs_order <- c('pLI (DEL)', 'pLI (DUP)', 'LOEUF (DEL)', 'LOEUF (DUP)',
              'pHI', 'pTS', 'pHI (binarized)', 'pTS (binarized)')
p_t_val$phenotype <- factor(p_t_val$phenotype)
p_t_val$risk_score <- factor(p_t_val$risk_score, levels = rs_order)

pheno_heatmap <- ggplot(data = p_t_val, aes(x = risk_score, y = phenotype, fill = t_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-6, 6), space = "Lab", name = "t-statistic") +
  theme_minimal() + xlab('CNV-S') + ylab('Phenotype') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05, paste0(sprintf("%.3f", t_val), "*"), ' ')),
            color = "grey30", size = 3.5) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05, paste0(sprintf("%.3f", t_val), "**"), ' '), fontface = 'bold'),
            color = "black", size = 3.5)
ggsave("output_figs_no_educ/cognitive_heatmap.png", pheno_heatmap, width = 14, height = 6, units = "cm")
ggsave("output_figs_no_educ/cognitive_heatmap.svg", pheno_heatmap, width = 14, height = 6, units = "cm")

# Forest plot
p_t_val$DEL_DUP <- rep(c("DEL", "DUP"), length.out = nrow(p_t_val))
p_t_val <- p_t_val %>%
  mutate(Annotation = case_when(p_val < 0.05 & p_val_adj >= 0.05 ~ '*',
                                p_val < 0.05 & p_val_adj < 0.05 ~ '**',
                                TRUE ~ ''))
lst_shapes <- c(19, 15)

forest_plot <- ggplot(data = p_t_val, aes(y = phenotype, x = standardized_beta,
                                          xmin = ci_lw, xmax = ci_up, col = risk_score, fill = risk_score)) +
  geom_linerange(size = 0.7, position = position_dodge(width = 0.5), alpha = ifelse(p_t_val$p_val >= 0.05, 0.5, 1)) +
  geom_vline(xintercept = 0, lty = 2) +
  geom_point(size = 3, shape = ifelse(p_t_val$DEL_DUP == 'DEL', lst_shapes[1], lst_shapes[2]),
             stroke = 0.5, position = position_dodge(width = 0.5),
             alpha = ifelse(p_t_val$p_val >= 0.05, 0.3, 1)) +
  scale_x_continuous(limits = c(-0.15, 0.25), breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
                    name = 'Standardized Effect (Beta Coefficient)') +
  scale_y_discrete(name = "Phenotypic Outcomes") +
  theme_minimal() +
  scale_color_discrete(name = "CNV-S") +
  scale_fill_discrete(name = "CNV-S") +
  theme(legend.position = "bottom") +
  geom_text(aes(x = ci_up, label = Annotation), hjust = -1, vjust = 0.7,
            position = position_dodge(width = 0.5), size = 5, color = 'black') +
  guides(color = guide_legend(override.aes = list(shape = rep(lst_shapes, length.out = 8))))
ggsave("output_figs_no_educ/main_forest_plot.png", forest_plot, width = 20, height = 12, units = "cm")
ggsave("output_figs_no_educ/main_forest_plot.svg", forest_plot, width = 20, height = 12, units = "cm")
message("Saved: output_figs_no_educ/main_forest_plot.png, .svg")

################################################################################
# 2. SEX STRATIFIED ANALYSIS (cognitive only)
################################################################################

str_covar_cogn_strat <- '+ age_at_visit'
str_covar_cogn_stdz_strat <- '+ scale(age_at_visit)'

scoreWithPheno_f <- scoreWithPheno[which(scoreWithPheno$msex == 0), ]
scoreWithPheno_m <- scoreWithPheno[which(scoreWithPheno$msex == 1), ]

p_t_val_f <- p_t_val_calculation_cogn(scoreWithPheno_f, lst_risk_scores, str_scores,
                                      str_scores_stdz, str_covar_cogn_strat, str_covar_cogn_stdz_strat,
                                      str_egvec, str_egvec_stdz)
p_t_val_m <- p_t_val_calculation_cogn(scoreWithPheno_m, lst_risk_scores, str_scores,
                                      str_scores_stdz, str_covar_cogn_strat, str_covar_cogn_stdz_strat,
                                      str_egvec, str_egvec_stdz)

# Rename for display
p_t_val_f[p_t_val_f == "cogn_global"] <- "Global Cognitive Function\n(n=663)"
p_t_val_f[p_t_val_f == "cogdx"] <- "Final Consensus\nCognitive Diagnosis\n(n=663)"
p_t_val_m[p_t_val_m == "cogn_global"] <- "Global Cognitive Function\n(n=348)"
p_t_val_m[p_t_val_m == "cogdx"] <- "Final Consensus\nCognitive Diagnosis\n(n=348)"
for (old_name in c("pTS_thresh", "pHI_thresh", "pli_del", "pli_dup", "loeuf_del", "loeuf_dup")) {
  new_name <- rename_mappings[old_name]
  if (!is.null(new_name)) {
    p_t_val_f[p_t_val_f == old_name] <- new_name
    p_t_val_m[p_t_val_m == old_name] <- new_name
  }
}

# FDR correction across both strata
n_cog <- nrow(p_t_val_f) + nrow(p_t_val_m)
p_vals <- c(p_t_val_f$p_val, p_t_val_m$p_val)
fdr_p_vals <- p.adjust(p_vals, 'fdr')
p_t_val_f$p_val_adj <- fdr_p_vals[1:nrow(p_t_val_f)]
p_t_val_m$p_val_adj <- fdr_p_vals[(nrow(p_t_val_f) + 1):n_cog]

# Heatmaps
plot_m <- ggplot(data = p_t_val_m, aes(x = factor(risk_score, level = rs_order), y = phenotype, fill = t_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-4, 4), space = "Lab", name = "t-statistic") +
  theme_minimal() + xlab('CNV-S') + ylab('Outcomes') + ggtitle('Male') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05, paste0(sprintf("%.2f", t_val), "*"), ' ')),
            color = "black", size = 3.5) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05, paste0(sprintf("%.2f", t_val), "**"), ' '), fontface = 'bold'),
            color = "black", size = 4)

plot_f <- ggplot(data = p_t_val_f, aes(x = factor(risk_score, level = rs_order), y = phenotype, fill = t_val)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-4, 4), space = "Lab", name = "t-statistic") +
  theme_minimal() + xlab('CNV-S') + ylab('Outcomes') + ggtitle('Female') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05, paste0(sprintf("%.2f", t_val), "*"), ' ')),
            color = "black", size = 3.5) +
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05, paste0(sprintf("%.2f", t_val), "**"), ' '), fontface = 'bold'),
            color = "black", size = 4)

g <- arrangeGrob(plot_m, plot_f)
ggsave("output_figs_no_educ/sex_stratified.png", g, width = 22, height = 10, units = "cm")
ggsave("output_figs_no_educ/sex_stratified.svg", g, width = 22, height = 10, units = "cm")
message("Saved: output_figs_no_educ/sex_stratified.png, .svg")

################################################################################
# 3. SEX INTERACTION ANOVA (cognitive only)
################################################################################

str_scores_interaction <- c('~ pli_del * msex + pli_dup * msex ',
                            '~ loeuf_del * msex + loeuf_dup * msex ',
                            '~ pHI * msex + pTS * msex ',
                            '~ pHI_thresh * msex + pTS_thresh * msex ')
str_covar_cogn_int <- '+ age_at_visit'

sex_interactions_anova_cogn <- function(scoreWithPheno, lst_risk_scores,
                                        str_scores_interaction, str_scores,
                                        str_covar_cogn_int) {
  lst_pheno <- c('cogn_global', 'cogdx_binom')
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 6)
  colnames(p_t_val) <- c('phenotype', 'risk_score', 'F', 'Pr(>F)', 'p', 'p_int_msex')

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      formula <- paste0(lst_pheno[i], str_scores[j], '+msex', str_covar_cogn_int, str_egvec)
      formula_int <- paste0(lst_pheno[i], str_scores_interaction[j], str_covar_cogn_int, str_egvec)
      is_glm <- lst_pheno[i] == 'cogdx_binom'

      fit1 <- if (is_glm) glm(formula, data = scoreWithPheno, family = 'binomial') else lm(formula, data = scoreWithPheno)
      fit2 <- if (is_glm) glm(formula_int, data = scoreWithPheno, family = 'binomial') else lm(formula_int, data = scoreWithPheno)
      fit <- anova(fit1, fit2, test = if (is_glm) "Chisq" else "F")

      anova_stat <- if (is_glm) fit$Deviance[2] else fit$F[2]
      anova_p <- if (is_glm) fit$`Pr(>Chi)`[2] else fit$`Pr(>F)`[2]

      # Interaction terms pli_del:msex (coef 4), pli_dup:msex (coef 6) in expanded formula
      idx_int1 <- 4
      idx_int2 <- 6

      a <- (i - 1) * 8 + 2 * j - 1
      p_t_val[a, ] <- c(lst_pheno[i], lst_risk_scores[2*j-1],
                        anova_stat, anova_p,
                        summary(fit1)$coefficients[,4][2],
                        summary(fit2)$coefficients[,4][idx_int1])
      p_t_val[a+1, ] <- c(lst_pheno[i], lst_risk_scores[2*j],
                          anova_stat, anova_p,
                          summary(fit1)$coefficients[,4][3],
                          summary(fit2)$coefficients[,4][idx_int2])
    }
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
  p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
  p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
  p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
  p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
  p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"
  p_t_val[p_t_val == "cogn_global"] <- "Global Cognitive Function"
  p_t_val[p_t_val == "cogdx_binom"] <- "Final Consensus Cognitive Diagnosis"
  return(p_t_val)
}

sex_int <- sex_interactions_anova_cogn(scoreWithPheno, lst_risk_scores,
                                      str_scores_interaction, str_scores,
                                      str_covar_cogn_int)
write.csv(sex_int, 'output_tables_no_educ/sex_stratified_outputs.csv', row.names = FALSE)
message("Saved: output_tables_no_educ/sex_stratified_outputs.csv")

################################################################################
# 4. PGS ASSOCIATION (CNV-S vs AD PGS) - cognitive only
################################################################################

prs_file <- "data/raw/raw_score_kunkle_full.tsv"
if (!file.exists(prs_file)) {
  message("Skipping PGS analysis: ", prs_file, " not found")
} else {
  scoreWithPheno_prs <- readRDS('data/scoreWithPheno.rds')
  scoreWithPheno_prs$cogdx_binom <- scoreWithPheno_prs$cogdx
  scoreWithPheno_prs$cogdx_binom[which(scoreWithPheno_prs$cogdx_binom == 2)] <- 1
  scoreWithPheno_prs$cogdx_binom[which(scoreWithPheno_prs$cogdx_binom == 3)] <- NA
  scoreWithPheno_prs$cogdx_binom[which(scoreWithPheno_prs$cogdx_binom == 5)] <- NA
  scoreWithPheno_prs$cogdx_binom[which(scoreWithPheno_prs$cogdx_binom == 1)] <- 0
  scoreWithPheno_prs$cogdx_binom[which(scoreWithPheno_prs$cogdx_binom == 4)] <- 1

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

  scoreWithPheno_prs <- merge(scoreWithPheno_prs, PRS, by.x = c('projid', 'study.x'), by.y = c('projid', 'study'))
  rm(temp, temp2, AD, ROS, ROSmaster_keylist)

  str_covar_cogn_prs <- '+ msex + age_at_visit'
  str_scores_prs <- c('+ pli_del + pli_dup ',
                      '+ loeuf_del + loeuf_dup ',
                      '+ pHI + pTS ',
                      '+ pHI_thresh + pTS_thresh ')

  scitf_note <- function(num, digit) formatC(num, format = "E", digits = digit)

  p_t_val_prs_calculation <- function(scoreWithPheno, lst_risk_scores, str_scores) {
    p_t_val <- matrix(nrow = 0, ncol = 9)
    colnames(p_t_val) <- c('phenotype', 'risk_score', 'PRS_t', 'PRS_p', 'adj.r', 'adj.r_rs',
                           'anova_p', 'AIC', 'AIC_w_RS')

    for (j in 1:length(str_scores)) {
      formula <- paste0('cogn_global', '~ PRS', str_covar_cogn_prs, str_egvec)
      formula_rs <- paste0('cogn_global', '~ PRS', str_scores[j], str_covar_cogn_prs, str_egvec)
      fit <- lm(formula, data = scoreWithPheno)
      fit_rs <- lm(formula_rs, data = scoreWithPheno)
      anova_res <- anova(fit, fit_rs)
      aic <- aictab(cand.set = list(fit, fit_rs), modnames = c('fit', 'fit_rs'))
      fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
      p_t_val <- rbind(p_t_val,
        c('cogn_global', lst_risk_scores[2*j-1],
          fit_coeff[,3][2], fit_coeff[,4][2],
          scitf_note(summary(fit)$adj.r.squared, 2),
          scitf_note(summary(fit_rs)$adj.r.squared, 2),
          scitf_note(anova_res$`Pr(>F)`[2], 2),
          aic$AICc[1], aic$AICc[2]),
        c('cogn_global', lst_risk_scores[2*j], NA, NA, NA, NA, NA, NA, NA))
    }

    for (j in 1:length(str_scores)) {
      formula <- paste0('cogdx_binom', '~ PRS', str_covar_cogn_prs, str_egvec)
      formula_rs <- paste0('cogdx_binom', '~ PRS', str_scores[j], str_covar_cogn_prs, str_egvec)
      fit <- glm(formula, data = scoreWithPheno, family = 'binomial')
      fit_rs <- glm(formula_rs, data = scoreWithPheno, family = 'binomial')
      fit_coeff <- scitf_note(summary(fit)$coefficients, 2)
      model_comparison <- lrtest(fit, fit_rs)
      p_t_val <- rbind(p_t_val,
        c('cogdx', lst_risk_scores[2*j-1],
          fit_coeff[,3][2], fit_coeff[,4][2],
          scitf_note(pR2(fit, type = "mcfadden")['r2ML'], 2),
          scitf_note(pR2(fit_rs, type = "mcfadden")['r2ML'], 2),
          scitf_note(model_comparison$`Pr(>Chisq)`[2], 20),
          AIC(fit), AIC(fit_rs)),
        c('cogdx', lst_risk_scores[2*j], NA, NA, NA, NA, NA, NA, NA))
    }

    p_t_val <- as.data.frame(p_t_val)
    p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
    p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
    p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
    p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
    p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
    p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"
    return(p_t_val)
  }

  p_t_val_prs <- p_t_val_prs_calculation(scoreWithPheno_prs, lst_risk_scores, str_scores_prs)
  p_t_val_prs$PRS_p_fdr <- p.adjust(as.numeric(as.character(p_t_val_prs$PRS_p)), method = "fdr")

  p_t_val_prs_filtered <- p_t_val_prs %>%
    filter(!is.na(anova_p)) %>%
    select(-c(AIC, AIC_w_RS))

  write.csv(p_t_val_prs_filtered, 'output_tables_no_educ/CNVRS_ANOVA_updated.csv', row.names = FALSE)
  message("Saved: output_tables_no_educ/CNVRS_ANOVA_updated.csv")
}

################################################################################
# DONE
################################################################################

message("
=============================================================================
Analysis complete. Cognitive models run WITHOUT years of education adjustment.
Outputs saved to:
  - output_tables_no_educ/
  - output_figs_no_educ/
=============================================================================
")
