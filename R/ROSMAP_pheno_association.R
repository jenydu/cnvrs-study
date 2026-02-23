# load packages
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

# load output from `ROSMAP_CNVRS_pheno_preprocess.R`
scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

# list of 10 eigenvectors
str_egvec <- '+ egvec1 + egvec2 + egvec3 + egvec4 + egvec5 + egvec6 + egvec7 + egvec8 + egvec9 + egvec10'
str_egvec_stdz <- '+ scale(egvec1) + scale(egvec2) + scale(egvec3) + scale(egvec4) + scale(egvec5) + scale(egvec6) + scale(egvec7) + scale(egvec8) + scale(egvec9) + scale(egvec10)'

# other covariates (neurological)
str_covar_autop <- '+ msex + age_death + pmi'
str_covar_autop_stdz <- '+ scale(msex) + scale(age_death) + scale(pmi)'

# other covariates (cognitive)
str_covar_cogn <- '+ msex + educ + age_at_visit'
str_covar_cogn_stdz <- '+ scale(msex) + scale(educ) + scale(age_at_visit)'

# 4 models, 2 scores per model
str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')
str_scores_stdz <- c('~ scale(pli_del) + scale(pli_dup)',
                     '~ scale(loeuf_del) + scale(loeuf_dup)',
                     '~ scale(pHI) + scale(pTS)',
                     '~ scale(pHI_thresh) + scale(pTS_thresh)')

# continuous variables
lst_pheno_autop <- c('sqrt(tangles)', 'sqrt(amyloid)', 'arteriol_scler', 'cvda_4gp2')
lst_pheno_cog <- c('cogn_global')

lst_risk_scores <- c('pli_del', 'pli_dup', 'loeuf_del', 'loeuf_dup',
                     'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

#######################
p_t_val_calculation <- function(scoreWithPheno, lst_risk_scores, str_scores,
                                lst_pheno_autop, lst_pheno_cog) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 7)
  colnames(p_t_val) <- c('phenotype', 'risk_score', 'standardized_beta',
                         'p_val', 'ci_lw', 'ci_up', 't_val')

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_autop, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j],
                               str_covar_autop_stdz, str_egvec_stdz)
      } else {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_cogn, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j],
                               str_covar_cogn_stdz, str_egvec_stdz)
      }

      fit <- lm(formula, data=scoreWithPheno)
      fit_stdz <- lm(stdz_formula, data=scoreWithPheno)

      confint1 <- confint(fit_stdz,
                          paste0('scale(', lst_risk_scores[2*j-1], ')'),
                          level=0.95)
      confint2 <- confint(fit_stdz,
                          paste0('scale(', lst_risk_scores[2*j], ')'),
                          level=0.95)
      standardized_betas <- lm.beta(fit_stdz)

      curr_row <- (i-1)*8+2*j-1
      p_t_val[curr_row, ] <- c(lst_pheno[i], lst_risk_scores[2*j-1],
                               unname(standardized_betas$coefficients[2]),
                               summary(fit)$coefficients[,4][2],
                               confint1[1], confint1[2],
                               summary(fit)$coefficients[,3][2])
      p_t_val[curr_row+1, ] <- c(lst_pheno[i], lst_risk_scores[2*j],
                                 unname(standardized_betas$coefficients[3]),
                                 summary(fit)$coefficients[,4][3],
                                 confint2[1], confint2[2],
                                 summary(fit)$coefficients[,3][3])
    }
  }

  ### categorical
  # cogdx
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec)
    stdz_formula <- paste0('scale(cogdx_binom)', str_scores_stdz[j],
                           str_covar_cogn_stdz, str_egvec_stdz)

    fit <- glm(formula, data=scoreWithPheno, family='binomial')
    fit_stdz <- lm(stdz_formula, data=scoreWithPheno)

    confint1 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j-1], ')'), level=0.95)
    confint2 <- confint(fit_stdz, paste0('scale(', lst_risk_scores[2*j], ')'), level=0.95)
    standardized_betas <- lm.beta(fit_stdz)

    new_rows <- matrix(nrow=2, ncol=7)
    new_rows[1, ] <- c('cogdx',
                       lst_risk_scores[2*j-1],
                       unname(standardized_betas$coefficients[2]),
                       summary(fit)$coefficients[,4][2],
                       confint1[1],
                       confint1[2],
                       summary(fit)$coefficients[,3][2])
    new_rows[2, ] <- c('cogdx',
                       lst_risk_scores[2*j],
                       unname(standardized_betas$coefficients[3]),
                       summary(fit)$coefficients[,4][3],
                       confint2[1],
                       confint2[2],
                       summary(fit)$coefficients[,3][3])

    p_t_val <- rbind(p_t_val, new_rows)
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val[c("standardized_beta",
            "p_val", "t_val",
            "ci_lw", "ci_up")] <- lapply(p_t_val[c("standardized_beta",
                                                   "p_val", "t_val",
                                                   "ci_lw", "ci_up")], as.double)

  # multiple testing correction
  p_t_val$p_val_adj <- p.adjust(p_t_val$p_val, 'fdr')
  return(p_t_val)
}

p_t_val <- p_t_val_calculation(scoreWithPheno, lst_risk_scores, str_scores,
                               lst_pheno_autop, lst_pheno_cog)

rename_mappings <- c(
  "pTS_thresh" = "pTS (binarized)",
  "pHI_thresh" = "pHI (binarized)",
  "pli_del" = "pLI (DEL)",
  "pli_dup" = "pLI (DUP)",
  "loeuf_del" = "LOEUF (DEL)",
  "loeuf_dup" = "LOEUF (DUP)",

  "sqrt(tangles)" = "PHF-Tau\n(n=1001)",
  "sqrt(amyloid)" = "Overall Amyloid Load\n(n=1005)",
  "cvda_4gp2" = "Cerebral Atherosclerosis\nRating\n(n=1004)",
  "cogn_global" = "Global Cognitive Function\n(n=1011)",
  "cogdx" = "Final Consensus\nCognitive Diagnosis\n(n=946)",
  "arteriol_scler" = "Arteriolosclerosis\n(n=1003)"
)

# Rename columns
for (old_name in names(rename_mappings)) {
  p_t_val[p_t_val == old_name] <- rename_mappings[old_name]
}

write.csv(p_t_val, file = "output_tables/p_t_val.csv", row.names = FALSE)

# PLOTS

# p_val < 0.05 & p_val_adj >= 0.05 -> '*'
# p_val < 0.05 & p_val_adj < 0.05 -> '**',

rs_order <- c('pLI (DEL)', 'pLI (DUP)', 'LOEUF (DEL)', 'LOEUF (DUP)',
              'pHI', 'pTS', 'pHI (binarized)', 'pTS (binarized)')

# heatmap
p_t_val$phenotype <- factor(p_t_val$phenotype)
p_t_val$risk_score <- factor(p_t_val$risk_score, levels = rs_order)

pheno_heatmap <- ggplot(data = p_t_val, aes(x=risk_score, y=phenotype, fill=t_val)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-6,6), space = "Lab", name="t-statistic") +
  theme_minimal() + xlab('Risk Score') + ylab('Phenotype') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) +
  # '*'
  geom_text(aes(risk_score, phenotype, label = if_else(p_val < 0.05 & p_val_adj >= 0.05, paste0(sprintf(sprintf("%.3f", t_val)), "*"), ' ')),
            color = "grey30", size = 3.5) + theme(axis.text.x = element_text(size = 8)) +
  # '**'
  geom_text(aes(risk_score, phenotype, label = if_else(p_val_adj < 0.05, paste0(as.character(sprintf("%.3f", t_val)), "**"), ' '), fontface = 'bold'),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8))
pheno_heatmap

# forest plot
p_t_val$DEL_DUP <- rep(c("DEL", "DUP"), length.out = nrow(p_t_val))
p_t_val <- p_t_val %>%
  mutate(Annotation = case_when(p_val < 0.05 & p_val_adj >= 0.05 ~ '*',
                                p_val < 0.05 & p_val_adj < 0.05 ~ '**',
                                TRUE ~ ''))
lst_shapes <- c(19, 15)

forest_plot <- ggplot(data=p_t_val, aes(y=phenotype, x=standardized_beta,
                                        xmin=ci_lw, xmax=ci_up, col=risk_score, fill=risk_score)) +
  geom_linerange(size=0.7, position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val >= 0.05, 0.5, 1)) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(size=3, shape=ifelse(p_t_val$DEL_DUP == 'DEL', lst_shapes[1], lst_shapes[2]),
             stroke = 0.5, position=position_dodge(width = 0.5),
             alpha=ifelse(p_t_val$p_val >= 0.05, 0.3, 1)) +
  scale_x_continuous(limits = c(-0.15, 0.25), breaks = c(-0.2,-0.1,0, 0.1, 0.2),
                     name = 'Standardized Effect (Beta Coefficient)') +
  scale_y_discrete(name = "Phenotypic Outcomes") +
  theme_minimal() +
  scale_color_discrete(name = "CNV-RS") +
  scale_fill_discrete(name = "CNV-RS") +
  theme(legend.position = "bottom") +
  # adding the stars '*'/'**'
  geom_text(aes(x = ci_up, label = Annotation), hjust = -1, vjust=0.7,
            position = position_dodge(width = 0.5), size = 5, color='black') +
  guides(color = guide_legend(
    override.aes=list(shape = rep(lst_shapes, length.out = 8))))
forest_plot

ggsave("output_figs/main_forest_plot.png", forest_plot, width = 20, height = 25, units = "cm")
ggsave(file="output_figs/main_forest_plot.svg", forest_plot, width = 20, height = 25, units = "cm")

###################
# cvda plots (fig.1)
generate_plot <- function(data, y_var, y_label) {
  ggplot(data = subset(data, !is.na(cvda_4gp2)),
         aes(x = as.factor(cvda_4gp2), y = !!sym(y_var), color = as.factor(cvda_4gp2))) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(cex = 0.7, size = 1) +
    theme_bw() +
    xlab(' ') +
    ylab(y_label) +
    scale_y_log10() +
    labs(color = 'Cerebral Atherosclerosis Rating') +
    scale_color_jco() +
    stat_summary(fun.y = "mean", color = 'black', shape = 18)
}

# Variables and labels
variables <- c("pli_del", "loeuf_del", "pHI", "pHI_thresh")
labels <- c("pLI Score (DEL)", "LOEUF Score (DEL)", "pHI", "pHI (binarized)")

# Create and arrange plots
cvda_plots <- lapply(seq_along(variables), function(i) {
  generate_plot(scoreWithPheno, variables[i], labels[i])
})

arranged_plots <- ggarrange(plotlist = cvda_plots, common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(left = text_grob("CNV-RS", rot = 90, x = unit(0, "npc"), y = unit(0.6, "npc")),
                  top = NULL, bottom = NULL, right = NULL) +
  theme(plot.margin = margin(0, 0, 0, 0.5, "cm"))

arranged_plots

ggsave("output_figs/cvda_plot.png", arranged_plots, width = 20, height = 16, units = "cm")
ggsave(file="output_figs/cvda_plot.svg", arranged_plots, width = 20, height = 16, units = "cm")

###################
# more plots
p1 <- ggplot(scoreWithPheno, aes(x = loeuf_dup, y = sqrt(tangles))) +
  geom_point(size = 1, color = '#000000BB', stroke = NA) +
  geom_smooth(method = lm) + theme_bw() +
  xlab('LOEUF Score (DUP)') +
  ylab('PHF-Tau')

p2 <- ggplot(scoreWithPheno, aes(x = pTS_thresh, y = sqrt(tangles)))  +
  geom_smooth(method = lm) + theme_bw() +
  geom_quasirandom(cex = 1, color = '#00000088', stroke = NA) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab('pTS Score (binarized)') +
  ylab(paste0('Tangles Density'))

p3 <- ggplot(scoreWithPheno, aes(x = pTS_thresh, y = sqrt(amyloid))) +
  geom_quasirandom(cex = 1, color = '#00000088', stroke = NA) +
  geom_smooth(method = lm) + theme_bw() +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  xlab('pTS Score (binarized)') +
  ylab(paste0('Overall Amyloid Load'))

p4 <- ggplot(scoreWithPheno, aes(x = pTS, y = cogn_global)) +
  geom_point(size = 1, color = '#00000088', stroke = NA) +
  geom_smooth(method = lm) + theme_bw() +
  xlab('pTS Score') +
  ylab(paste0('Global Cognitive Function'))

p5 <- ggplot(scoreWithPheno, aes(x = pTS_thresh, y = cogn_global)) +
  geom_quasirandom(cex = 1, size = 0.7, color = '#00000088', stroke = NA) +
  geom_smooth(method = lm) + theme_bw() +
  xlab('pTS Score (binarized)') +
  ylab(paste0('Global Cognitive Function')) +
  scale_x_continuous(breaks = seq(0, 10, by = 2))

p6 <- ggplot(data = subset(scoreWithPheno, !is.na(cogdx_binom)),
             aes(x = as.factor(cogdx_binom), y = pTS, color = as.factor(cogdx_binom))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex = 0.7, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Final Consensus\nCognitive Diagnosis') +
  ylab('pTS') + scale_y_log10() +
  scale_color_jco()

p7 <- ggplot(data = subset(scoreWithPheno,!is.na(cogdx_binom)),
             aes(x = as.factor(cogdx_binom), y = pTS_thresh, color = as.factor(cogdx_binom))) +
  geom_boxplot()+
  geom_beeswarm(cex = 0.3, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Final Consensus\nCognitive Diagnosis') +
  ylab('pTS (binarized)') +
  scale_color_jco()

p8 <- ggplot(data = subset(scoreWithPheno, !is.na(arteriol_scler)),
             aes(x = as.factor(arteriol_scler), y = pli_del, color = as.factor(arteriol_scler))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex = 0.7, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Arteriolosclerosis') +
  ylab('pLI (DEL)') +
  scale_y_log10() +
  scale_color_jco()

p9 <- ggplot(data = subset(scoreWithPheno, !is.na(arteriol_scler)),
             aes(x = as.factor(arteriol_scler), y = loeuf_del, color = as.factor(arteriol_scler))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex = 0.7, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Arteriolosclerosis') +
  ylab('LOEUF (DEL)') +
  scale_y_log10() +
  scale_color_jco()

p10 <- ggplot(data = subset(scoreWithPheno, !is.na(arteriol_scler)),
              aes(x = as.factor(arteriol_scler), y = pHI, color = as.factor(arteriol_scler))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex = 0.7, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Arteriolosclerosis') +
  ylab('pHI') +
  scale_y_log10() +
  scale_color_jco()

p11 <- ggplot(data = subset(scoreWithPheno, !is.na(arteriol_scler)),
              aes(x = as.factor(arteriol_scler), y = pHI_thresh, color = as.factor(arteriol_scler))) +
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(cex = 0.7, size = 1) +
  theme_bw() +
  xlab(' ') + labs(color = 'Arteriolosclerosis') +
  ylab('pHI (binarized)') +
  scale_y_log10() +
  scale_color_jco()

supp_plotA <- ggarrange(p1, p4, p2, p3, p5, ncol = 2, nrow = 3) +
  theme(plot.margin = margin(0.5, 0, 0.5, 0, "cm"))
supp_plotB <- ggarrange(p8, p9, p10, p11, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.5, 0, 2, 0, "cm"))
supp_plotC <- ggarrange(p6, p7, common.legend = TRUE, legend = "bottom") +
  theme(plot.margin = margin(0.5, 0, 0, 0, "cm"))

# Plot arrangement
supp_plot <- ggarrange(supp_plotA, supp_plotB, supp_plotC,
                       ncol = 2, nrow = 2,
                       labels = c("A", "B", "C"),
                       heights = c(2, 1))

ggsave("output_figs/supp_plot.png", supp_plot, width = 30, height = 25, units = "cm")
ggsave(file="output_figs/supp_plot.svg", supp_plot, width = 30, height = 25, units = "cm")
