##########
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
##############

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

lst_risk_scores <- c('pli_del','pli_dup','loeuf_del','loeuf_dup',
                     'pHI','pTS', 'pHI_thresh','pTS_thresh')

# make a new column named cogdx_binom: cogdx = {1,2} -> 0; cogdx = {4} -> 1;
scoreWithPheno$cogdx_binom <- scoreWithPheno$cogdx
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==2)] <- 1
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==3)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==5)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==1)] <- 0
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom==4)] <- 1

# saveRDS(scoreWithPheno, 'scoreWithPheno.rds')

p_t_val_calculation <- function(scoreWithPheno, lst_risk_scores, str_scores,
                                lst_pheno_autop, lst_pheno_cog) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)
  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 7)
  cnames <- c('phenotype', 'risk_score', 'standardized_beta','p_val', 'ci_lw', 'ci_up', 't_val')
  colnames(p_t_val) <- cnames

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_autop, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j], str_covar_autop_stdz, str_egvec_stdz)
      } else {
        formula <- paste0(lst_pheno[i], str_scores[j], str_covar_cogn, str_egvec)
        stdz_formula <- paste0('scale(', lst_pheno[i], ')', str_scores_stdz[j], str_covar_cogn_stdz, str_egvec_stdz)
      }
      fit <- lm(formula, data=scoreWithPheno)
      fit_stdz <- lm(stdz_formula, data=scoreWithPheno)

      confint1 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j-1],')'),level=0.95)
      confint2 <- confint(fit_stdz, paste0('scale(',lst_risk_scores[2*j],')'), level=0.95)

      standardized_betas <- lm.beta(fit_stdz)

      a <-(i-1)*8+2*j-1
      p_t_val[a, ] <- c(lst_pheno[i],
                        lst_risk_scores[2*j-1],
                        unname(standardized_betas$coefficients[2]),
                        summary(fit)$coefficients[,4][2],
                        confint1[1],
                        confint1[2],
                        summary(fit)$coefficients[,3][2])
      p_t_val[a+1, ] <- c(lst_pheno[i],
                          lst_risk_scores[2*j],
                          unname(standardized_betas$coefficients[3]),
                          summary(fit)$coefficients[,4][3],
                          confint2[1],
                          confint2[2],
                          summary(fit)$coefficients[,3][3])
    }
  }
  ### categorical
  # cogdx
  for (j in 1:length(str_scores)) {
    formula <- paste0('cogdx_binom', str_scores[j], str_covar_cogn, str_egvec)

    fit <- glm(formula, data=scoreWithPheno, family='binomial')

    new_rows <- matrix(nrow=2, ncol=7)
    colnames(new_rows) <- cnames

    confint1 <- confint(fit, lst_risk_scores[2*j-1], level=0.95)
    confint2 <- confint(fit, lst_risk_scores[2*j], level=0.95)

    new_rows[1, ] <- c('cogdx',
                       lst_risk_scores[2*j-1],
                       summary(fit)$coefficients[,1][2],
                       summary(fit)$coefficients[,4][2],
                       confint1[1],
                       confint1[2],
                       summary(fit)$coefficients[,3][2])
    new_rows[2, ] <- c('cogdx',
                       lst_risk_scores[2*j],
                       summary(fit)$coefficients[,1][3],
                       summary(fit)$coefficients[,4][3],
                       confint2[1],
                       confint2[2],
                       summary(fit)$coefficients[,3][3])
    p_t_val <- rbind(p_t_val, new_rows)
  }

  p_t_val <- as.data.frame(p_t_val)
  p_t_val$standardized_beta <- as.double(p_t_val$standardized_beta)
  p_t_val$p_val <- as.double(p_t_val$p_val)
  p_t_val$t_val <- as.double(p_t_val$t_val)
  p_t_val$ci_lw <- as.double(p_t_val$ci_lw)
  p_t_val$ci_up <- as.double(p_t_val$ci_up)

  # multiple testing correction
  p_t_val$p_val_adj <- p.adjust(p_t_val$p_val, 'fdr')
  p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
  p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
  p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
  p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
  p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
  p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"

  p_t_val[p_t_val == "sqrt(tangles)"] <- 'PHF-Tau'
  p_t_val[p_t_val == "sqrt(amyloid)"] <- 'Overall Amyloid Load'
  p_t_val[p_t_val == "cvda_4gp2"] <- "Cerebral\nAtherosclerosis Rating"
  p_t_val[p_t_val == "cogn_global"] <- "Global Cognitive Function"
  p_t_val[p_t_val == "cogdx"] <- 'Final Consensus\nCognitive Diagnosis'
  p_t_val[p_t_val == "arteriol_scler"] <- 'Arteriolosclerosis'
  return(p_t_val)
}

##################
# sex stratified #
##################
# other covariates (neurological)
str_covar_autop <- '+ age_death + pmi'
str_covar_autop_stdz <- '+ scale(age_death) + scale(pmi)'

# other covariates (cognitive)
str_covar_cogn <- '+ educ + age_at_visit'
str_covar_cogn_stdz <- '+ scale(educ) + scale(age_at_visit)'

p_t_val <- p_t_val_calculation(scoreWithPheno, lst_risk_scores, str_scores,
                               lst_pheno_autop, lst_pheno_cog)
scoreWithPheno_f <- scoreWithPheno[which(scoreWithPheno$msex==0),]
p_t_val_f <- p_t_val_calculation(scoreWithPheno_f, lst_risk_scores, str_scores,
                                 lst_pheno_autop, lst_pheno_cog)

p_t_val_f[p_t_val_f == 'PHF-Tau'] <- 'PHF-Tau\n(n=654)'
p_t_val_f[p_t_val_f == 'Overall Amyloid Load'] <- 'Overall Amyloid Load\n(n=658)'
p_t_val_f[p_t_val_f == "Cerebral\nAtherosclerosis Rating"] <- "Cerebral Atherosclerosis\nRating\n(n=656)"
p_t_val_f[p_t_val_f == "Global Cognitive Function"] <- "Global Cognitive Function\n(n=663)"
p_t_val_f[p_t_val_f == 'Final Consensus\nCognitive Diagnosis'] <- 'Final Consensus\nCognitive Diagnosis\n(n=663)'
p_t_val_f[p_t_val_f == 'Arteriolosclerosis'] <- 'Arteriolosclerosis\n(n=659)'

scoreWithPheno_m <- scoreWithPheno[which(scoreWithPheno$msex==1),]
p_t_val_m <- p_t_val_calculation(scoreWithPheno_m, lst_risk_scores, str_scores,
                                 lst_pheno_autop, lst_pheno_cog)

p_t_val_m[p_t_val_m == 'PHF-Tau'] <- 'PHF-Tau\n(n=347)'
p_t_val_m[p_t_val_m == 'Overall Amyloid Load'] <- 'Overall Amyloid Load\n(n=347)'
p_t_val_m[p_t_val_m == "Cerebral\nAtherosclerosis Rating"] <- "Cerebral Atherosclerosis\nRating\n(n=348)"
p_t_val_m[p_t_val_m == "Global Cognitive Function"] <- "Global Cognitive Function\n(n=348)"
p_t_val_m[p_t_val_m == 'Final Consensus\nCognitive Diagnosis'] <- 'Final Consensus\nCognitive Diagnosis\n(n=348)'
p_t_val_m[p_t_val_m == 'Arteriolosclerosis'] <- 'Arteriolosclerosis\n(n=344)'

# fdr correction

p_vals <- c(p_t_val_f$p_val, p_t_val_m$p_val)
fdr_p_vals <- p.adjust(p_vals, 'fdr')
p_t_val_f$p_val_adj <- fdr_p_vals[1:48]
p_t_val_m$p_val_adj <- fdr_p_vals[49:96]
#########

rs_order <- c('pLI (DEL)', 'pLI (DUP)',
              'LOEUF (DEL)', 'LOEUF (DUP)',
              'pHI', 'pTS',
              'pHI (binarized)', 'pTS (binarized)')
plot_m <- ggplot(data = p_t_val_m, aes(x=factor(risk_score, level = rs_order), y=phenotype, fill=t_val)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-4,4), space = "Lab", name="t-statistic") +
  theme_minimal()+
  xlab('CNV-RS')+
  ylab('Outcomes')+
  ggtitle('Male')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  geom_text(aes(risk_score, phenotype,label = if_else(p_val < 0.05 & p_val_adj >= 0.05 ,  paste0(as.character(sprintf("%.2f", t_val)), "*"), ' ')),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8))+
  geom_text(aes(risk_score, phenotype,label = if_else(p_val_adj < 0.05, paste0(as.character(sprintf("%.2f", t_val)), "**"), ' '), fontface = 'bold'),
            color = "black", size = 4) + theme(axis.text.x = element_text(size = 8))

plot_f <- ggplot(data = p_t_val_f, aes(x=factor(risk_score, level = rs_order), y=phenotype, fill=t_val)) +
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-4,4), space = "Lab", name="t-statistic") +
  theme_minimal()+
  xlab('CNV-RS')+
  ylab('Outcomes')+
  ggtitle('Female')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  geom_text(aes(risk_score, phenotype,label = if_else(p_val < 0.05 & p_val_adj >= 0.05 , paste0(as.character(sprintf("%.2f", t_val)), "*"), ' ')),
            color = "black", size = 3.5) + theme(axis.text.x = element_text(size = 8))+
  geom_text(aes(risk_score, phenotype,label = if_else(p_val_adj < 0.05, paste0(as.character(sprintf("%.2f", t_val)), "**"), ' '), fontface = 'bold'),
            color = "black", size = 4) + theme(axis.text.x = element_text(size = 8))

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

grid.arrange(plot_m, plot_f)
g <- arrangeGrob(plot_m, plot_f)
ggsave("output_figs/sex_stratified.png", g, width = 22, height = 22, units = "cm")
ggsave("output_figs/sex_stratified.svg", g, width = 22, height = 22, units = "cm")
rm(plot_m,plot_f)

p_t_val_f$phenotype = factor (p_t_val_f$phenotype)
forest_plot <- ggplot(data=p_t_val_f, aes(y=phenotype, x=standardized_beta, xmin=ci_lw, xmax=ci_up, col=risk_score,fill=risk_score)) +
  geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val_adj >= 0.05, 0.5, 1)) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(size=3, shape=20, stroke = 0.5,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val_adj >= 0.05, 0.3, 1)) +
  scale_x_continuous(limits = c(-0.15, 0.25), breaks = c(-0.2,-0.1,0, 0.1, 0.2), name = 'Standardized Effect (Beta Coefficient)') +
  scale_y_discrete(name="Phenotype") +
  theme_minimal()

p_t_val_m$phenotype = factor (p_t_val_m$phenotype)
forest_plot <- ggplot(data=p_t_val_m, aes(y=phenotype, x=standardized_beta, xmin=ci_lw, xmax=ci_up, col=risk_score,fill=risk_score)) +
  geom_linerange(size=0.7,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val_adj >= 0.05, 0.5, 1)) +
  geom_vline(xintercept=0, lty=2) +
  geom_point(size=3, shape=20, stroke = 0.5,position=position_dodge(width = 0.5), alpha=ifelse(p_t_val$p_val_adj >= 0.05, 0.3, 1)) +
  scale_x_continuous(limits = c(-0.15, 0.25), breaks = c(-0.2,-0.1,0, 0.1, 0.2), name = 'Standardized Effect (Beta Coefficient)') +
  scale_y_discrete(name="Phenotype") +
  theme_minimal()

#####################
# interaction terms #
#####################

lst_pheno_autop <- c('sqrt(tangles)', 'sqrt(amyloid)', 'arteriol_scler', 'cvda_4gp2')
lst_pheno_cog <- c('cogn_global')

# scoreWithPheno$cvda_4gp2 <- as.integer(scoreWithPheno$cvda_4gp2)

str_scores_interaction <- c('~ pli_del * msex + pli_dup * msex ',
                            '~ loeuf_del * msex + loeuf_dup * msex ',
                            '~ pHI * msex + pTS * msex ',
                            '~ pHI_thresh * msex + pTS_thresh * msex ')

str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')

# other covariates (neurological)
str_covar_autop <- '+ age_death + pmi'

# other covariates (cognitive)
str_covar_cogn <- '+ educ + age_at_visit'

sex_interactions_anova <- function(scoreWithPheno, lst_risk_scores,
                                   str_scores_interaction, str_scores,
                                   lst_pheno_autop, lst_pheno_cog) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)

  p_t_val <- matrix(nrow = length(lst_pheno) * length(lst_risk_scores), ncol = 6)
  colnames(p_t_val) <- c('phenotype', 'risk_score', 'F','Pr(>F)', 'p', 'p_int_msex')

  for (i in 1:length(lst_pheno)) {
    for (j in 1:length(str_scores)) {
      if (i <= length(lst_pheno_autop)) {
        formula <- paste0(lst_pheno[i], str_scores[j], '+msex', str_covar_autop, str_egvec)
        formula_int <- paste0(lst_pheno[i], str_scores_interaction[j], str_covar_autop, str_egvec)
      } else {
        formula <- paste0(lst_pheno[i], str_scores[j], '+msex', str_covar_cogn, str_egvec)
        formula_int <- paste0(lst_pheno[i], str_scores_interaction[j], str_covar_cogn, str_egvec)
      }

      fit1 <- lm(formula, data=scoreWithPheno)
      fit2 <- lm(formula_int, data=scoreWithPheno)

      fit <- anova(fit1, fit2)
      # print(fit)
      a <-(i-1)*8+2*j-1
      p_t_val[a, ] <- c(lst_pheno[i],
                        lst_risk_scores[2*j-1],
                        fit$F[2],
                        fit$`Pr(>F)`[2],
                        summary(fit1)$coefficients[,4][2],
                        summary(fit2)$coefficients[,4][17])
      p_t_val[a+1, ] <- c(lst_pheno[i],
                          lst_risk_scores[2*j],
                          fit$F[2],
                          fit$`Pr(>F)`[2],
                          summary(fit1)$coefficients[,4][3],
                          summary(fit2)$coefficients[,4][18])
    }
  }

  p_t_val <- as.data.frame(p_t_val)

  p_t_val[p_t_val == "pTS_thresh"] <- "pTS (binarized)"
  p_t_val[p_t_val == "pHI_thresh"] <- "pHI (binarized)"
  p_t_val[p_t_val == "pli_del"] <- "pLI (DEL)"
  p_t_val[p_t_val == "pli_dup"] <- "pLI (DUP)"
  p_t_val[p_t_val == "loeuf_del"] <- "LOEUF (DEL)"
  p_t_val[p_t_val == "loeuf_dup"] <- "LOEUF (DUP)"

  p_t_val[p_t_val == "sqrt(tangles)"] <- 'PHF-Tau'
  p_t_val[p_t_val == "sqrt(amyloid)"] <- 'Overall Amyloid Load'
  p_t_val[p_t_val == "cvda_4gp2"] <- "Cerebral\nAtherosclerosis Rating"
  p_t_val[p_t_val == "cogn_global"] <- "Global Cognitive Function"
  # p_t_val[p_t_val == "cogdx"] <- 'Final Consensus\nCognitive Diagnosis'
  p_t_val[p_t_val == "arteriol_scler"] <- 'Arteriolosclerosis'

  return(p_t_val)
}

sex_int <- sex_interactions_anova(scoreWithPheno, lst_risk_scores,
                                  str_scores_interaction, str_scores,
                                  lst_pheno_autop, lst_pheno_cog)
write.csv(sex_int, 'output_tables/sex_stratified_outputs.csv', row.names = FALSE)

# https://stats.stackexchange.com/questions/331244/how-to-test-if-an-interaction-is-significant-interaction-terms-or-model-compari
# https://bookdown.org/ndphillips/YaRrr/comparing-regression-models-with-anova.html

wo_interaction_formula <- paste0('sqrt(amyloid)', "~ pHI_thresh + pTS_thresh + msex", str_covar_autop, str_egvec)
w_interaction_formula <- paste0('sqrt(amyloid)', "~ pHI_thresh *msex + pTS_thresh* msex", str_covar_autop, str_egvec)
anova(lm(formula=wo_interaction_formula, data=scoreWithPheno),
      lm(formula=w_interaction_formula, data=scoreWithPheno))
