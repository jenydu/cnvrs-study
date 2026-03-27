################################################################################
# Comparison Table: Cognitive Models With vs Without Education Adjustment
################################################################################
#
# Side-by-side comparison of cognitive association results:
#   - WITH years of education (educ) covariate  [original]
#   - WITHOUT years of education (educ) covariate [reviewer revision]
#
# Also checks whether removing educ changes sample sizes (missing data check).
#
# Outputs:
#   output_tables/cognitive_educ_comparison.csv
#   output_tables/cognitive_educ_comparison_n.csv
################################################################################

if (file.exists("data/scoreWithPheno.rds")) {
  # already in project root
} else if (file.exists("../data/scoreWithPheno.rds")) {
  setwd("..")
} else {
  stop("Cannot find data/scoreWithPheno.rds. Run from project root or R/ directory.")
}

library(stats)
library(dplyr)
library(lm.beta)
library(pscl)

dir.create("output_tables", showWarnings = FALSE)

################################################################################
# Load data
################################################################################

scoreWithPheno <- readRDS("data/scoreWithPheno.rds")

scoreWithPheno$cogdx_binom <- scoreWithPheno$cogdx
scoreWithPheno$cogdx_binom[scoreWithPheno$cogdx_binom == 2] <- 1
scoreWithPheno$cogdx_binom[scoreWithPheno$cogdx_binom == 3] <- NA
scoreWithPheno$cogdx_binom[scoreWithPheno$cogdx_binom == 5] <- NA
scoreWithPheno$cogdx_binom[scoreWithPheno$cogdx_binom == 1] <- 0
scoreWithPheno$cogdx_binom[scoreWithPheno$cogdx_binom == 4] <- 1

str_egvec <- paste(paste0("egvec", 1:10), collapse = " + ")
str_egvec <- paste0("+ ", str_egvec)
str_egvec_stdz <- paste(paste0("scale(egvec", 1:10, ")"), collapse = " + ")
str_egvec_stdz <- paste0("+ ", str_egvec_stdz)

lst_risk_scores <- c("pli_del", "pli_dup", "loeuf_del", "loeuf_dup",
                     "pHI", "pTS", "pHI_thresh", "pTS_thresh")

str_scores <- c("~ pli_del + pli_dup ",
                "~ loeuf_del + loeuf_dup ",
                "~ pHI + pTS ",
                "~ pHI_thresh + pTS_thresh ")

str_scores_stdz <- c("~ scale(pli_del) + scale(pli_dup)",
                     "~ scale(loeuf_del) + scale(loeuf_dup)",
                     "~ scale(pHI) + scale(pTS)",
                     "~ scale(pHI_thresh) + scale(pTS_thresh)")

################################################################################
# Core function: run models and extract stats + n
################################################################################

run_cogn_models <- function(data, covar_str, covar_stdz_str) {
  results <- list()

  for (j in seq_along(str_scores)) {
    # --- cogn_global (continuous) ---
    formula_lm <- paste0("cogn_global", str_scores[j], covar_str, str_egvec)
    formula_stdz <- paste0("scale(cogn_global)", str_scores_stdz[j], covar_stdz_str, str_egvec_stdz)

    fit_lm   <- lm(formula_lm,   data = data)
    fit_stdz <- lm(formula_stdz, data = data)

    ci1 <- confint(fit_stdz, paste0("scale(", lst_risk_scores[2*j-1], ")"), level = 0.95)
    ci2 <- confint(fit_stdz, paste0("scale(", lst_risk_scores[2*j]   , ")"), level = 0.95)
    sb  <- lm.beta(fit_stdz)
    n   <- nobs(fit_lm)

    results[[length(results)+1]] <- data.frame(
      phenotype   = "cogn_global",
      risk_score  = lst_risk_scores[2*j-1],
      beta        = unname(sb$coefficients[2]),
      ci_lw       = ci1[1],
      ci_up       = ci1[2],
      p_val       = summary(fit_lm)$coefficients[, 4][2],
      t_val       = summary(fit_lm)$coefficients[, 3][2],
      n           = n,
      stringsAsFactors = FALSE
    )
    results[[length(results)+1]] <- data.frame(
      phenotype   = "cogn_global",
      risk_score  = lst_risk_scores[2*j],
      beta        = unname(sb$coefficients[3]),
      ci_lw       = ci2[1],
      ci_up       = ci2[2],
      p_val       = summary(fit_lm)$coefficients[, 4][3],
      t_val       = summary(fit_lm)$coefficients[, 3][3],
      n           = n,
      stringsAsFactors = FALSE
    )

    # --- cogdx_binom (logistic) ---
    formula_glm  <- paste0("cogdx_binom", str_scores[j], covar_str, str_egvec)
    formula_stdz2 <- paste0("scale(cogdx_binom)", str_scores_stdz[j], covar_stdz_str, str_egvec_stdz)

    fit_glm   <- glm(formula_glm, data = data, family = "binomial")
    fit_stdz2 <- lm(formula_stdz2, data = data)

    ci3 <- confint(fit_stdz2, paste0("scale(", lst_risk_scores[2*j-1], ")"), level = 0.95)
    ci4 <- confint(fit_stdz2, paste0("scale(", lst_risk_scores[2*j]   , ")"), level = 0.95)
    sb2 <- lm.beta(fit_stdz2)
    n2  <- nobs(fit_glm)

    results[[length(results)+1]] <- data.frame(
      phenotype   = "cogdx",
      risk_score  = lst_risk_scores[2*j-1],
      beta        = unname(sb2$coefficients[2]),
      ci_lw       = ci3[1],
      ci_up       = ci3[2],
      p_val       = summary(fit_glm)$coefficients[, 4][2],
      t_val       = summary(fit_glm)$coefficients[, 3][2],
      n           = n2,
      stringsAsFactors = FALSE
    )
    results[[length(results)+1]] <- data.frame(
      phenotype   = "cogdx",
      risk_score  = lst_risk_scores[2*j],
      beta        = unname(sb2$coefficients[3]),
      ci_lw       = ci4[1],
      ci_up       = ci4[2],
      p_val       = summary(fit_glm)$coefficients[, 4][3],
      t_val       = summary(fit_glm)$coefficients[, 3][3],
      n           = n2,
      stringsAsFactors = FALSE
    )
  }

  out <- do.call(rbind, results)
  out$p_val_adj <- p.adjust(out$p_val, "fdr")
  out
}

################################################################################
# Run both model sets
################################################################################

# With education
res_educ <- run_cogn_models(
  data         = scoreWithPheno,
  covar_str    = "+ msex + educ + age_at_visit",
  covar_stdz_str = "+ scale(msex) + scale(educ) + scale(age_at_visit)"
)

# Without education
res_no_educ <- run_cogn_models(
  data         = scoreWithPheno,
  covar_str    = "+ msex + age_at_visit",
  covar_stdz_str = "+ scale(msex) + scale(age_at_visit)"
)

################################################################################
# Check sample size differences
################################################################################

n_check <- res_educ %>%
  select(phenotype, risk_score, n_educ = n) %>%
  left_join(res_no_educ %>% select(phenotype, risk_score, n_no_educ = n),
            by = c("phenotype", "risk_score")) %>%
  mutate(n_diff = n_no_educ - n_educ)

message("\n--- Sample size check ---")
print(n_check)

################################################################################
# Build side-by-side comparison table
################################################################################

fmt_beta_ci <- function(beta, ci_lw, ci_up) {
  sprintf("%.2f [%.2f, %.2f]", beta, ci_lw, ci_up)
}

fmt_p <- function(p) {
  ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.2f", p))
}

sig_label <- function(p_raw, p_adj) {
  case_when(
    p_adj <  0.05 ~ "** (FDR sig.)",
    p_raw <  0.05 ~ "*  (nom. sig.)",
    TRUE          ~ ""
  )
}

comparison <- res_educ %>%
  rename_with(~ paste0(.x, "_educ"), -c(phenotype, risk_score)) %>%
  left_join(
    res_no_educ %>% rename_with(~ paste0(.x, "_no_educ"), -c(phenotype, risk_score)),
    by = c("phenotype", "risk_score")
  ) %>%
  mutate(
    Phenotype   = recode(phenotype,
      "cogn_global" = "Global Cognitive Function",
      "cogdx"       = "Final Consensus Cognitive Diagnosis"),
    CNV_Score   = recode(risk_score,
      "pli_del"    = "pLI (DEL)",
      "pli_dup"    = "pLI (DUP)",
      "loeuf_del"  = "LOEUF (DEL)",
      "loeuf_dup"  = "LOEUF (DUP)",
      "pHI_thresh" = "pHI (binarized)",
      "pTS_thresh" = "pTS (binarized)"),
    `n` = n_educ,
    `beta [95% CI] — with educ`  = fmt_beta_ci(beta_educ,    ci_lw_educ,    ci_up_educ),
    `p — with educ`              = fmt_p(p_val_educ),
    `p_adj — with educ`          = fmt_p(p_val_adj_educ),
    `Sig — with educ`            = sig_label(p_val_educ, p_val_adj_educ),
    `beta [95% CI] — no educ`    = fmt_beta_ci(beta_no_educ, ci_lw_no_educ, ci_up_no_educ),
    `p — no educ`                = fmt_p(p_val_no_educ),
    `p_adj — no educ`            = fmt_p(p_val_adj_no_educ),
    `Sig — no educ`              = sig_label(p_val_no_educ, p_val_adj_no_educ)
  ) %>%
  select(Phenotype, CNV_Score,
         `n`,
         `beta [95% CI] — with educ`, `p — with educ`, `p_adj — with educ`, `Sig — with educ`,
         `beta [95% CI] — no educ`,   `p — no educ`,   `p_adj — no educ`,   `Sig — no educ`)

write.csv(comparison, "output_tables/cognitive_educ_comparison.csv", row.names = FALSE)
message("Saved: output_tables/cognitive_educ_comparison.csv")

################################################################################
# Print summary to console
################################################################################

message("\n=== COMPARISON TABLE: With educ vs Without educ ===\n")
message("* = nominally significant (p < 0.05); ** = FDR-significant (p_adj < 0.05)\n")
print(as.data.frame(comparison), row.names = FALSE)

message("\n=== ORIGINAL (with educ): Significant results ===")
sig_orig <- res_educ %>%
  filter(p_val < 0.05) %>%
  mutate(fdr_sig = p_val_adj < 0.05)
if (nrow(sig_orig) == 0) {
  message("  No nominally significant results.")
} else {
  print(sig_orig %>% select(phenotype, risk_score, beta, ci_lw, ci_up, p_val, p_val_adj, n, fdr_sig))
}

message("\n=== REVISED (no educ): Significant results ===")
sig_new <- res_no_educ %>%
  filter(p_val < 0.05) %>%
  mutate(fdr_sig = p_val_adj < 0.05)
if (nrow(sig_new) == 0) {
  message("  No nominally significant results.")
} else {
  print(sig_new %>% select(phenotype, risk_score, beta, ci_lw, ci_up, p_val, p_val_adj, n, fdr_sig))
}

message("\n=== DONE ===")
