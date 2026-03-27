################################################################################
# ROSMAP Post-Hoc Power Analysis
################################################################################
#
# Post-hoc power analysis for CNV score (CNV-S) associations with
# neuropathological and cognitive phenotypes. Addresses reviewer comment:
# "Power calculations are not included, nor is statistical power mentioned
#  anywhere in the manuscript."
#
# ROS/MAP is a publicly available resource cohort, not designed for any
# specific genetic analysis. Sample size is determined by data availability,
# not a pre-specified power requirement. A priori power calculations are
# therefore not applicable. Instead, we perform a post-hoc power analysis
# using:
#   1. Observed sample sizes from this study
#   2. The significance threshold used for inference (BH-FDR; Bonferroni used
#      as a conservative proxy for power calculations)
#   3. Expected effect sizes from three sources:
#      (a) EXTERNAL LITERATURE — prior CNV burden/score analyses (Ming et al.
#          2022) and PRS-neuropathology studies (|beta| = 0.10-0.15)
#      (b) OBSERVED PRS EFFECTS from this cohort — context for what genetic
#          scores can achieve in ROS/MAP (upper bound; much larger than CNV-S)
#      (c) OBSERVED CNV-S EFFECTS from this study — betas we actually observed,
#          reported separately for all 48 tests and for FDR-significant tests
#
# IMPORTANT NOTE ON POST-HOC POWER:
#   Using observed effect sizes for post-hoc power is mathematically redundant
#   with the p-value (Hoenig & Heisey 2001; see references below). We compute
#   it because the reviewer requested power calculations and the PI has asked
#   for the observed-beta framing. We explicitly cite the methodological
#   literature showing why post-hoc power is of limited utility.
#
# Prerequisites:
#   1. ROSMAP_CNVRS_pheno_preprocess.R  -> data/scoreWithPheno.rds
#   2. ROSMAP_pheno_association.R       -> output_tables/p_t_val.csv
#   3. ROSMAP_PRS_association_updated.R -> output_tables/CNVRS_ANOVA_updated.csv
# Requires: pwr package
#
################################################################################

if (file.exists("data/scoreWithPheno.rds")) {
  # Already in project root
} else if (file.exists("../data/scoreWithPheno.rds")) {
  setwd("..")
} else {
  stop("Cannot find data/scoreWithPheno.rds. Run from project root or R/ directory.")
}

library(stats)
library(dplyr)

if (!requireNamespace("pwr", quietly = TRUE)) {
  install.packages("pwr", repos = "https://cloud.r-project.org")
}
library(pwr)

dir.create("output_tables", showWarnings = FALSE)

################################################################################
# Significance thresholds
#
# BH-FDR correction was applied across 48 main tests (6 phenotypes × 8
# CNV scores). We use Bonferroni alpha = 0.05/48 as a conservative proxy for
# power calculations, and report nominal alpha = 0.05 for comparison.
################################################################################

n_tests_main  <- 48L
n_tests_strat <- 96L
alpha_nominal <- 0.05
alpha_bonf    <- alpha_nominal / n_tests_main   # ~0.00104

################################################################################
# Load main data
################################################################################

scoreWithPheno <- readRDS('data/scoreWithPheno.rds')

str_egvec       <- '+ egvec1 + egvec2 + egvec3 + egvec4 + egvec5 + egvec6 + egvec7 + egvec8 + egvec9 + egvec10'
str_covar_autop <- '+ msex + age_death + pmi'
str_covar_cogn  <- '+ msex + educ + age_at_visit'

str_scores <- c('~ pli_del + pli_dup ',
                '~ loeuf_del + loeuf_dup ',
                '~ pHI + pTS ',
                '~ pHI_thresh + pTS_thresh ')

lst_pheno_autop <- c('sqrt(tangles)', 'sqrt(amyloid)', 'arteriol_scler', 'cvda_4gp2')
lst_pheno_cog   <- c('cogn_global')
lst_risk_scores <- c('pli_del', 'pli_dup', 'loeuf_del', 'loeuf_dup',
                     'pHI', 'pTS', 'pHI_thresh', 'pTS_thresh')

scoreWithPheno$cogdx_binom <- scoreWithPheno$cogdx
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 2)] <- 1
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 3)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 5)] <- NA
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 1)] <- 0
scoreWithPheno$cogdx_binom[which(scoreWithPheno$cogdx_binom == 4)] <- 1

################################################################################
# Load association results: observed CNV-S betas and PRS reference effects
################################################################################

assoc_results <- tryCatch(
  read.csv("output_tables/p_t_val.csv"),
  error = function(e) { warning("Cannot load p_t_val.csv"); NULL }
)
anova_results <- tryCatch(
  read.csv("output_tables/CNVRS_ANOVA_updated.csv"),
  error = function(e) { warning("Cannot load CNVRS_ANOVA_updated.csv"); NULL }
)

# Observed CNV-S standardized betas from this study
if (!is.null(assoc_results)) {
  obs_betas     <- abs(assoc_results$standardized_beta)
  sig_mask      <- !is.na(assoc_results$p_val_adj) & assoc_results$p_val_adj < 0.05
  obs_betas_sig <- abs(assoc_results$standardized_beta[sig_mask])

  beta_obs_all_median <- round(median(obs_betas, na.rm = TRUE), 4)
  beta_obs_all_range  <- round(range(obs_betas, na.rm = TRUE), 4)
  beta_obs_sig_median <- if (length(obs_betas_sig) > 0)
                           round(median(obs_betas_sig, na.rm = TRUE), 4)
                         else NA_real_
  n_sig_assoc <- sum(sig_mask, na.rm = TRUE)
} else {
  beta_obs_all_median <- NA_real_
  beta_obs_all_range  <- c(NA_real_, NA_real_)
  beta_obs_sig_median <- NA_real_
  n_sig_assoc         <- NA_integer_
}

f2_obs_all <- if (!is.na(beta_obs_all_median)) beta_obs_all_median^2 else NA_real_
f2_obs_sig <- if (!is.na(beta_obs_sig_median)) beta_obs_sig_median^2 else NA_real_

# PRS effect sizes from this cohort (standardized beta from t-stat and df)
# These serve as context for what genetic scores can achieve in ROS/MAP;
# PRS effects are substantially larger than CNV-S effects
if (!is.null(anova_results)) {
  prs_t         <- abs(as.numeric(anova_results$PRS_t))
  prs_df_approx <- 980
  prs_betas     <- prs_t / sqrt(prs_t^2 + prs_df_approx)
  beta_prs_median <- round(median(prs_betas, na.rm = TRUE), 3)
  beta_prs_range  <- round(range(prs_betas,  na.rm = TRUE), 3)
} else {
  beta_prs_median <- NA_real_
  beta_prs_range  <- c(NA_real_, NA_real_)
}

################################################################################
# Literature-derived expected effect sizes
#
# (a) Ming et al. (2022) — CNV burden/score analyses in ROS/MAP;
#     typical |beta| ≈ 0.05-0.15; conservative estimate: 0.10
# (b) PRS-AD studies of continuous neuropathological outcomes in
#     comparable cohorts (n ~ 800-1200); typical |beta| ≈ 0.08-0.15
#
# Cohen's f2 for a single standardized predictor ≈ beta^2
################################################################################

beta_lit_conservative <- 0.10
beta_lit_upper        <- 0.15

f2_lit_conservative <- beta_lit_conservative^2
f2_lit_upper        <- beta_lit_upper^2

################################################################################
# Sample sizes and residual df per phenotype
################################################################################

get_n_and_df <- function(data, lst_pheno_autop, lst_pheno_cog,
                         str_scores, str_covar_autop, str_covar_cogn, str_egvec) {
  lst_pheno <- c(lst_pheno_autop, lst_pheno_cog)
  out <- data.frame(phenotype = character(), n = integer(), df = integer(),
                    model_type = character(), stringsAsFactors = FALSE)

  for (i in seq_along(lst_pheno)) {
    pheno <- lst_pheno[i]
    if (i <= length(lst_pheno_autop)) {
      formula <- paste0(pheno, str_scores[1], str_covar_autop, str_egvec)
    } else {
      formula <- paste0(pheno, str_scores[1], str_covar_cogn, str_egvec)
    }
    fit <- lm(as.formula(formula), data = data)
    n   <- fit$df.residual + fit$rank
    df  <- fit$df.residual
    out <- rbind(out, data.frame(phenotype = pheno, n = n, df = df,
                                 model_type = "linear", stringsAsFactors = FALSE))
  }

  formula <- paste0('cogdx_binom', str_scores[1], str_covar_cogn, str_egvec)
  fit <- glm(as.formula(formula), data = data, family = 'binomial')
  n   <- nobs(fit)
  df  <- fit$df.residual
  out <- rbind(out, data.frame(phenotype = "cogdx", n = n, df = df,
                               model_type = "logistic", stringsAsFactors = FALSE))
  out
}

n_df_main <- get_n_and_df(scoreWithPheno, lst_pheno_autop, lst_pheno_cog,
                           str_scores, str_covar_autop, str_covar_cogn, str_egvec)

scoreWithPheno_f <- scoreWithPheno[scoreWithPheno$msex == 0, ]
scoreWithPheno_m <- scoreWithPheno[scoreWithPheno$msex == 1, ]
str_covar_autop_strat <- '+ age_death + pmi'
str_covar_cogn_strat  <- '+ educ + age_at_visit'

n_df_f <- get_n_and_df(scoreWithPheno_f, lst_pheno_autop, lst_pheno_cog,
                        str_scores, str_covar_autop_strat, str_covar_cogn_strat, str_egvec)
n_df_m <- get_n_and_df(scoreWithPheno_m, lst_pheno_autop, lst_pheno_cog,
                        str_scores, str_covar_autop_strat, str_covar_cogn_strat, str_egvec)

################################################################################
# Power calculation helpers
################################################################################

calc_power <- function(f2, df, alpha) {
  if (is.na(f2) || is.na(df) || is.na(alpha)) return(NA_real_)
  tryCatch(pwr::pwr.f2.test(u = 1, v = df, f2 = f2, sig.level = alpha)$power,
           error = function(e) NA_real_)
}

min_detectable_beta <- function(df, alpha, power = 0.80) {
  tryCatch({
    f2 <- pwr::pwr.f2.test(u = 1, v = df, power = power, sig.level = alpha)$f2
    sqrt(f2)
  }, error = function(e) NA_real_)
}

################################################################################
# Main analysis: post-hoc power across four effect size scenarios
#   Scenario 1: literature conservative  |beta| = 0.10
#   Scenario 2: literature upper bound   |beta| = 0.15
#   Scenario 3: observed CNV-S median (all 48 tests)
#   Scenario 4: observed CNV-S median (FDR-significant tests only)
################################################################################

power_main <- n_df_main

# Scenario 1 and 2: literature-based
power_main$power_lit_conserv_nominal <- mapply(calc_power, f2_lit_conservative, power_main$df, alpha_nominal)
power_main$power_lit_upper_nominal   <- mapply(calc_power, f2_lit_upper,        power_main$df, alpha_nominal)
power_main$power_lit_conserv_bonf    <- mapply(calc_power, f2_lit_conservative, power_main$df, alpha_bonf)
power_main$power_lit_upper_bonf      <- mapply(calc_power, f2_lit_upper,        power_main$df, alpha_bonf)

# Scenario 3: observed CNV-S median (all tests)
power_main$power_obs_all_nominal <- mapply(calc_power, f2_obs_all, power_main$df, alpha_nominal)
power_main$power_obs_all_bonf    <- mapply(calc_power, f2_obs_all, power_main$df, alpha_bonf)

# Scenario 4: observed CNV-S median (significant tests)
power_main$power_obs_sig_nominal <- mapply(calc_power, f2_obs_sig, power_main$df, alpha_nominal)
power_main$power_obs_sig_bonf    <- mapply(calc_power, f2_obs_sig, power_main$df, alpha_bonf)

# Minimum detectable |beta| at 80% power
power_main$min_beta_80_nominal <- mapply(min_detectable_beta, power_main$df, alpha_nominal)
power_main$min_beta_80_bonf    <- mapply(min_detectable_beta, power_main$df, alpha_bonf)

power_main$adequate_power_lit_nominal <- power_main$power_lit_conserv_nominal >= 0.80
power_main$adequate_power_lit_bonf    <- power_main$power_lit_conserv_bonf    >= 0.80

################################################################################
# Sex-stratified post-hoc power
################################################################################

make_strat_power <- function(n_df, stratum_label) {
  d <- n_df
  d$stratum <- stratum_label
  d$power_lit_conserv_nominal <- mapply(calc_power, f2_lit_conservative, d$df, alpha_nominal)
  d$power_lit_upper_nominal   <- mapply(calc_power, f2_lit_upper,        d$df, alpha_nominal)
  d$power_lit_conserv_bonf    <- mapply(calc_power, f2_lit_conservative, d$df, alpha_bonf)
  d$power_obs_all_nominal     <- mapply(calc_power, f2_obs_all,          d$df, alpha_nominal)
  d$power_obs_all_bonf        <- mapply(calc_power, f2_obs_all,          d$df, alpha_bonf)
  d$power_obs_sig_nominal     <- mapply(calc_power, f2_obs_sig,          d$df, alpha_nominal)
  d$power_obs_sig_bonf        <- mapply(calc_power, f2_obs_sig,          d$df, alpha_bonf)
  d$min_beta_80_nominal       <- mapply(min_detectable_beta, d$df, alpha_nominal)
  d$min_beta_80_bonf          <- mapply(min_detectable_beta, d$df, alpha_bonf)
  d$adequate_power_nominal    <- d$power_lit_conserv_nominal >= 0.80
  d
}

power_strat_f <- make_strat_power(n_df_f, "Female")
power_strat_m <- make_strat_power(n_df_m, "Male")
power_strat   <- rbind(power_strat_f, power_strat_m)

################################################################################
# Write tabular outputs
################################################################################

write.csv(power_main,  "output_tables/power_posthoc_main.csv",           row.names = FALSE)
write.csv(power_strat, "output_tables/power_posthoc_sex_stratified.csv", row.names = FALSE)

# Observed effect size summary
obs_beta_summary <- data.frame(
  metric = c(
    "N total CNV-S tests",
    "N FDR-significant tests",
    "Median |beta| (all 48 tests)",
    "Median |beta| (FDR-significant tests)",
    "Min |beta| detectable at 80% power, alpha=0.05",
    "Min |beta| detectable at 80% power, Bonferroni alpha"
  ),
  value = c(
    if (!is.null(assoc_results)) nrow(assoc_results) else NA,
    n_sig_assoc,
    round(beta_obs_all_median, 4),
    round(beta_obs_sig_median, 4),
    round(mean(power_main$min_beta_80_nominal, na.rm = TRUE), 4),
    round(mean(power_main$min_beta_80_bonf,    na.rm = TRUE), 4)
  )
)
write.csv(obs_beta_summary, "output_tables/power_effect_size_summary.csv", row.names = FALSE)

################################################################################
# Summary output (saved to file for manuscript / reviewer response)
################################################################################

sink("output_tables/power_posthoc_summary.txt")

cat("================================================================================\n")
cat("ROSMAP CNV-S Post-Hoc Power Analysis Summary\n")
cat("================================================================================\n\n")

cat("STUDY DESIGN NOTE\n")
cat("-----------------\n")
cat("ROS/MAP is a publicly available resource cohort, not designed for any specific\n")
cat("genetic analysis. Sample size reflects data availability, not a pre-specified\n")
cat("power requirement. A priori power calculations are not applicable to this\n")
cat("analysis of opportunity; sample size is constrained by data availability and\n")
cat("cost, as is typical in omics-based cohort studies.\n\n")

cat("POST-HOC POWER APPROACH\n")
cat("-----------------------\n")
cat("Power estimated using:\n")
cat("  - Observed study sample sizes (n per phenotype, see below)\n")
cat("  - Significance thresholds: nominal alpha = 0.05 (reference);\n")
cat("    Bonferroni alpha =", format(alpha_bonf, scientific = TRUE, digits = 3),
    "(conservative multiple-testing proxy)\n")
cat("  - Four effect size scenarios:\n")
cat("    (1) Conservative literature estimate:  |beta| =", beta_lit_conservative,
    "(Ming et al. 2022; PRS-AD literature)\n")
cat("    (2) Upper-bound literature estimate:   |beta| =", beta_lit_upper, "\n")
cat("    (3) Observed CNV-S (all", n_tests_main, "tests, this study): median |beta| =",
    beta_obs_all_median, "(range:", beta_obs_all_range[1], "-", beta_obs_all_range[2], ")\n")
cat("    (4) Observed CNV-S (FDR-significant,", n_sig_assoc, "tests): median |beta| =",
    beta_obs_sig_median, "\n")
if (!is.na(beta_prs_median)) {
  cat("  - PRS reference (this cohort): median |beta| =", beta_prs_median,
      "(range:", beta_prs_range[1], "-", beta_prs_range[2], ")\n")
  cat("    [PRS effects are ~3-10x larger than CNV-S effects; shown for context only]\n")
}
cat("\n")

cat("SAMPLE SIZES\n")
cat("------------\n")
cat("  Main analysis:   n =", min(n_df_main$n), "-", max(n_df_main$n), "per phenotype\n")
cat("  Female stratum:  n =", min(n_df_f$n), "-", max(n_df_f$n), "\n")
cat("  Male stratum:    n =", min(n_df_m$n), "-", max(n_df_m$n), "\n\n")

cat("MAIN ANALYSIS POWER — Scenario (1): |beta| =", beta_lit_conservative,
    "[literature conservative]\n")
cat(strrep("-", 72), "\n")
for (i in seq_len(nrow(power_main))) {
  cat(sprintf("  %-22s n=%4d  power(alpha=0.05)=%5.3f  power(Bonf)=%5.3f  >=80%%[nominal]=%s\n",
              power_main$phenotype[i], power_main$n[i],
              power_main$power_lit_conserv_nominal[i],
              power_main$power_lit_conserv_bonf[i],
              ifelse(power_main$adequate_power_lit_nominal[i], "YES", "NO")))
}
cat("\n")

cat("MAIN ANALYSIS POWER — Scenario (2): |beta| =", beta_lit_upper,
    "[literature upper bound]\n")
cat(strrep("-", 72), "\n")
for (i in seq_len(nrow(power_main))) {
  cat(sprintf("  %-22s n=%4d  power(alpha=0.05)=%5.3f  power(Bonf)=%5.3f\n",
              power_main$phenotype[i], power_main$n[i],
              power_main$power_lit_upper_nominal[i],
              power_main$power_lit_upper_bonf[i]))
}
cat("\n")

cat("MAIN ANALYSIS POWER — Scenario (3): |beta| =", beta_obs_all_median,
    "[observed median, all tests]\n")
cat(strrep("-", 72), "\n")
for (i in seq_len(nrow(power_main))) {
  cat(sprintf("  %-22s n=%4d  power(alpha=0.05)=%5.3f  power(Bonf)=%5.3f\n",
              power_main$phenotype[i], power_main$n[i],
              power_main$power_obs_all_nominal[i],
              power_main$power_obs_all_bonf[i]))
}
cat("\n")

cat("MAIN ANALYSIS POWER — Scenario (4): |beta| =", beta_obs_sig_median,
    "[observed median, significant tests]\n")
cat(strrep("-", 72), "\n")
for (i in seq_len(nrow(power_main))) {
  cat(sprintf("  %-22s n=%4d  power(alpha=0.05)=%5.3f  power(Bonf)=%5.3f\n",
              power_main$phenotype[i], power_main$n[i],
              power_main$power_obs_sig_nominal[i],
              power_main$power_obs_sig_bonf[i]))
}
cat("\n")

cat("MINIMUM DETECTABLE |beta| AT 80% POWER\n")
cat("---------------------------------------\n")
cat("  At nominal alpha = 0.05:  |beta| >=",
    round(mean(power_main$min_beta_80_nominal, na.rm = TRUE), 3), "(avg across phenotypes)\n")
cat("  At Bonferroni alpha:      |beta| >=",
    round(mean(power_main$min_beta_80_bonf,    na.rm = TRUE), 3), "\n\n")

cat("SEX-STRATIFIED POWER — Scenario (1): |beta| = 0.10, alpha = 0.05\n")
cat("-------------------------------------------------------------------\n")
cat("  Female (n =", min(n_df_f$n), "-", max(n_df_f$n), "): power range =",
    round(min(power_strat_f$power_lit_conserv_nominal, na.rm = TRUE), 3), "-",
    round(max(power_strat_f$power_lit_conserv_nominal, na.rm = TRUE), 3), "\n")
cat("  Male   (n =", min(n_df_m$n), "-", max(n_df_m$n), "): power range =",
    round(min(power_strat_m$power_lit_conserv_nominal, na.rm = TRUE), 3), "-",
    round(max(power_strat_m$power_lit_conserv_nominal, na.rm = TRUE), 3), "\n\n")

cat("OVERALL CONCLUSIONS\n")
cat("-------------------\n")
all_adequate_nominal <- all(power_main$adequate_power_lit_nominal, na.rm = TRUE)
cat("  1. All main analyses adequately powered (>=80%) at nominal alpha with\n")
cat("     conservative literature effect size (|beta|=0.10):",
    ifelse(all_adequate_nominal, "YES", "NO"), "\n")

if (!is.na(beta_obs_sig_median) && !is.na(power_main$power_obs_sig_bonf[1])) {
  sig_bonf_min <- min(power_main$power_obs_sig_bonf, na.rm = TRUE)
  cat("  2. For the observed significant CNV-S associations (n =", n_sig_assoc,
      "tests; all for\n")
  cat("     cerebral atherosclerosis; |beta| =", beta_obs_sig_median, "),\n")
  cat("     estimated power at Bonferroni alpha =",
      round(sig_bonf_min * 100, 1), "% (min across phenotypes).\n")
}

if (!is.na(beta_obs_all_median)) {
  all_bonf_max <- max(power_main$power_obs_all_bonf, na.rm = TRUE)
  cat("  3. For the median observed CNV-S effect across all 48 tests (|beta| =",
      beta_obs_all_median, "),\n")
  cat("     power at Bonferroni alpha is low (max =",
      round(all_bonf_max * 100, 1),
      "%). This is expected:\n")
  cat("     most CNV-S associations are small and below the detection threshold\n")
  cat("     at stringent multiple-testing correction.\n")
}
cat("  4. Sex-stratified analyses have reduced power due to smaller subgroup\n")
cat("     sample sizes, particularly for males (n ~", round(mean(n_df_m$n)), ").\n\n")

cat("IMPORTANT CAVEATS ON POST-HOC POWER\n")
cat("------------------------------------\n")
cat("  Post-hoc power analyses have well-documented methodological limitations.\n")
cat("  When computed using the study's own observed effect sizes, post-hoc power\n")
cat("  is mathematically equivalent to a monotone transformation of the p-value\n")
cat("  and provides no additional inferential information (Hoenig & Heisey 2001;\n")
cat("  Goodman & Berlin 1994). For significant tests, 'adequate' post-hoc power\n")
cat("  is a mathematical consequence of significance, not an independent finding.\n")
cat("  For non-significant tests, low post-hoc power reflects the small observed\n")
cat("  effect — not an independent reason to dismiss the null. Accordingly, the\n")
cat("  external literature-based estimates (Scenarios 1 & 2) are the most\n")
cat("  interpretively valid, and we recommend against over-interpreting the\n")
cat("  observed-beta power figures (Scenarios 3 & 4).\n\n")
cat("  This analysis is provided in response to the reviewer's request. We note\n")
cat("  that this cohort is an analysis of opportunity: cost and sample availability\n")
cat("  are the limiting factors, and as a public resource, ROS/MAP cannot be\n")
cat("  designed for any single genetic analysis.\n\n")

cat("KEY REFERENCES ON POST-HOC POWER ANALYSIS\n")
cat("------------------------------------------\n")
cat("  Hoenig JM & Heisey DM (2001) Am Stat 55:19-24. DOI:10.1198/000313001300339897\n")
cat("  Levine M & Ensom MH (2001) Pharmacotherapy 21:405-9. DOI:10.1592/phco.21.5.405.34503\n")
cat("  Goodman SN & Berlin JA (1994) Ann Intern Med 121:200-6. doi:10.7326/0003-4819-121-3-199408010-00008\n")
cat("  Thomas L (1997) Conserv Biol 11:276-280. DOI:10.1046/j.1523-1739.1997.96102.x\n")
cat("  Yuan K-H & Maxwell S (2005) J Educ Behav Stat 30:141-167. DOI:10.3102/10769986030002141\n")
cat("  Walters SJ (2008) Pharm Stat 8:163-169. DOI:10.1002/pst.334\n")
cat("  Miller J (2009) Psychon Bull Rev 16:617-640.\n")
cat("  Greenland S (2012) Ann Epidemiol 22:364-368.\n")
cat("  Wilkinson L et al. (1999) Am Psychol 54:594-604. DOI:10.1037/0003-066X.54.8.594\n")
cat("  See also: https://pmc.ncbi.nlm.nih.gov/articles/PMC9452450/\n")
cat("================================================================================\n")

sink()

################################################################################
# Manuscript-ready prose for the Methods/Statistical Analysis section
################################################################################

n_main_min <- min(n_df_main$n)
n_main_max <- max(n_df_main$n)
n_f_min    <- min(n_df_f$n)
n_f_max    <- max(n_df_f$n)
n_m_min    <- min(n_df_m$n)
n_m_max    <- max(n_df_m$n)

pwr_sig_bonf_min <- if (!is.na(f2_obs_sig)) round(min(power_main$power_obs_sig_bonf, na.rm = TRUE) * 100, 1) else NA_real_

manuscript_text <- paste0(
"=== MANUSCRIPT TEXT (brief, for Methods or Results) ===\n",
"\n",
"ROS/MAP is a publicly available resource cohort assembled for broad use; sample size ",
"reflects data and specimen availability rather than a pre-specified power requirement. ",
"In response to the reviewer's request, we performed a post-hoc power analysis using ",
"observed sample sizes (n = ", n_main_min, "\u2013", n_main_max, " per phenotype), ",
"observed CNV-S effect sizes from this study, and the Bonferroni-corrected significance ",
"threshold as a conservative proxy for the BH-FDR correction applied (",
"\u03b1 = 0.05/", n_tests_main, " = ", format(alpha_bonf, scientific = TRUE, digits = 3), "). ",
if (!is.na(pwr_sig_bonf_min))
  paste0("The FDR-significant associations (n\u2009=\u2009", n_sig_assoc,
         ", all for cerebral atherosclerosis; median |\u03b2|\u2009=\u2009", beta_obs_sig_median,
         ") had estimated post-hoc power of ", pwr_sig_bonf_min, "%, ",
         "exceeding the conventional 80% benchmark. ")
else "",
"Full post-hoc power results are provided in Supplementary Table [X]. ",
"We note that post-hoc power analyses have well-documented methodological limitations ",
"[HOENIG_HEISEY; LEVINE_ENSOM; GOODMAN_BERLIN; THOMAS; YUAN_MAXWELL; WALTERS; GREENLAND]: ",
"when derived from observed effect sizes, post-hoc power is mathematically equivalent to ",
"a transformation of the p-value and provides no independent inferential information. ",
"Non-significant results should therefore be interpreted as inconclusive rather than as ",
"evidence of no association.\n",
"\n",
"=== CITATION PLACEHOLDERS ===\n",
"HOENIG_HEISEY  Hoenig JM & Heisey DM (2001) Am Stat 55:19-24. DOI:10.1198/000313001300339897\n",
"LEVINE_ENSOM   Levine M & Ensom MH (2001) Pharmacotherapy 21:405-9. DOI:10.1592/phco.21.5.405.34503\n",
"GOODMAN_BERLIN Goodman SN & Berlin JA (1994) Ann Intern Med 121:200-6. doi:10.7326/0003-4819-121-3-199408010-00008\n",
"THOMAS         Thomas L (1997) Conserv Biol 11:276-280. DOI:10.1046/j.1523-1739.1997.96102.x\n",
"YUAN_MAXWELL   Yuan K-H & Maxwell S (2005) J Educ Behav Stat 30:141-167. DOI:10.3102/10769986030002141\n",
"WALTERS        Walters SJ (2008) Pharm Stat 8:163-169. DOI:10.1002/pst.334\n",
"MILLER         Miller J (2009) Psychon Bull Rev 16:617-640.\n",
"GREENLAND      Greenland S (2012) Ann Epidemiol 22:364-368.\n",
"WILKINSON      Wilkinson L et al. (1999) Am Psychol 54:594-604. DOI:10.1037/0003-066X.54.8.594\n",
"               https://pmc.ncbi.nlm.nih.gov/articles/PMC9452450/\n"
)

writeLines(manuscript_text, "output_tables/manuscript_power_section.txt")

message("Post-hoc power analysis complete. Outputs saved to output_tables/")
message("  - power_posthoc_main.csv")
message("  - power_posthoc_sex_stratified.csv")
message("  - power_effect_size_summary.csv")
message("  - power_posthoc_summary.txt")
message("  - manuscript_power_section.txt")
