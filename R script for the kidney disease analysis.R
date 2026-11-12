# ==============================================================================
# Manuscript Pipeline: Heterogeneous Machine Learning and Deep Neural Network
# Ensembles for Screening Renal Dysfunction from Routine Serum Biochemistry:
# A Nationally Representative NHANES Study (2021-2023 Cycle: _L)
# ==============================================================================

# 1. SETUP & DIRECTORY CONFIGURATION -------------------------------------------
setwd("C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\Kidney disease")

dir.create("Tables", showWarnings = FALSE)
dir.create("Plots", showWarnings = FALSE)

required_packages <- c(
  # Data I/O & Wrangling
  "haven", "tidyverse", "janitor", "scales", "xtable", "gt", "gtsummary",
  # Classical Biostatistics, Survey & Econometrics
  "survey", "rstatix", "moments", "quantreg", "fitdistrplus", "survival",
  "splines", "broom",
  # Visualization
  "patchwork", "hexbin", "pROC", "PRROC",
  # Heterogeneous Machine Learning, Deep Learning & Ensembles
  "caret", "ranger", "xgboost", "e1071", "glmnet", "nnet", "mclust"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)
invisible(lapply(required_packages, library, character.only = TRUE))

# Prevent MASS::select() from masking dplyr::select()
select <- dplyr::select

# Strictly NO titles/subtitles on figures for clean LaTeX \caption{} integration
theme_set(
  theme_classic(base_size = 11.5) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 9.5),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 9.5),
      legend.text = element_text(size = 9),
      strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.6),
      strip.text = element_text(face = "bold", size = 10)
    )
)

# Helper to export both CSV and booktabs LaTeX (.tex) tables safely
export_table <- function(df, base_name, caption_text = "", label_text = "") {
  clean_df <- df %>%
    mutate(across(where(is.character), ~ iconv(.x, from = "", to = "UTF-8", sub = "u")))
  write_csv(clean_df, file.path("Tables", paste0(base_name, ".csv")))
  tex_out <- print(
    xtable(as.data.frame(clean_df), caption = caption_text, label = label_text),
    include.rownames = FALSE,
    booktabs = TRUE,
    print.results = FALSE
  )
  writeLines(tex_out, file.path("Tables", paste0(base_name, ".tex")))
}

# Helper to sanitize Windows-1252 / Latin-1 NHANES SAS labels (e.g., <b5>g/dL)
clean_sas_labels <- function(df) {
  for (col in names(df)) {
    lbl <- attr(df[[col]], "label")
    if (!is.null(lbl)) {
      raw_bytes <- charToRaw(lbl)
      raw_bytes[raw_bytes == as.raw(0xb5)] <- as.raw(0x75)
      lbl_clean <- rawToChar(raw_bytes)
      lbl_clean <- iconv(lbl_clean, from = "latin1", to = "UTF-8", sub = "")
      attr(df[[col]], "label") <- lbl_clean
    }
  }
  return(df)
}

# 2. MULTI-FILE NHANES MERGE (ALB_CR_L + BIOPRO_L + DEMO_L) --------------------
find_nhanes_file <- function(prefix) {
  list.files(pattern = paste0("^", prefix), ignore.case = TRUE)[1]
}

f_alb  <- find_nhanes_file("ALB_CR_L")
f_bio  <- find_nhanes_file("BIOPRO_L")
f_demo <- find_nhanes_file("DEMO_L")

if (is.na(f_alb)) stop("ALB_CR_L.xpt not found in working directory.")
df_alb <- read_xpt(f_alb) %>% clean_sas_labels()

has_biopro <- !is.na(f_bio)
has_demo   <- !is.na(f_demo)

merged_df <- df_alb
if (has_demo) {
  df_demo <- read_xpt(f_demo) %>%
    clean_sas_labels() %>%
    dplyr::select(any_of(c(
      "SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH3", "DMDEDUC2", "INDFMPIR",
      "WTINT2YR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA"
    )))
  merged_df <- left_join(merged_df, df_demo, by = "SEQN")
}

if (has_biopro) {
  df_bio <- read_xpt(f_bio) %>% clean_sas_labels()
  merged_df <- left_join(merged_df, df_bio, by = "SEQN")
}

# 3. CLINICAL ENDPOINT DERIVATION (KDIGO ACR + CKD-EPI 2021 eGFR) --------------
calc_egfr_2021 <- function(scr, age, sex) {
  kappa <- ifelse(sex == 2, 0.7, 0.9)
  alpha <- ifelse(sex == 2, -0.241, -0.302)
  sex_mult <- ifelse(sex == 2, 1.012, 1.000)
  142 * (pmin(scr / kappa, 1)^alpha) * (pmax(scr / kappa, 1)^(-1.200)) * (0.9938^age) * sex_mult
}

kidney_df <- merged_df %>%
  mutate(
    Albumin_LLOD = factor(URDUMALC, levels = c(0, 1), labels = c("Quantifiable (>=LLOD)", "Left-Censored (<LLOD)")),
    Creatinine_LLOD = factor(URDUCRLC, levels = c(0, 1), labels = c("Quantifiable (>=LLOD)", "Left-Censored (<LLOD)")),
    Albuminuria_Stage = case_when(
      is.na(URDACT) ~ NA_character_,
      URDACT < 30 ~ "A1: Normal/Mild (<30 mg/g)",
      URDACT >= 30 & URDACT <= 300 ~ "A2: Moderate (30-300 mg/g)",
      URDACT > 300 ~ "A3: Severe (>300 mg/g)"
    ),
    Albuminuria_Stage = factor(
      Albuminuria_Stage,
      levels = c("A1: Normal/Mild (<30 mg/g)", "A2: Moderate (30-300 mg/g)", "A3: Severe (>300 mg/g)")
    ),
    Elevated_ACR = factor(ifelse(URDACT >= 30, "Yes", "No"), levels = c("No", "Yes")),
    log10_URXUMA = log10(URXUMA),
    log10_URXUCR = log10(URXUCR),
    log10_URDACT = log10(URDACT),
    ln_URXUMA    = log(URXUMA),
    ln_URXUCR    = log(URXUCR)
  )

if (has_demo) {
  kidney_df <- kidney_df %>%
    mutate(
      Sex = factor(RIAGENDR, levels = c(1, 2), labels = c("Male", "Female")),
      Age_Group = cut(RIDAGEYR, breaks = c(18, 39, 59, 120), labels = c("18-39 yrs", "40-59 yrs", "60+ yrs"), include.lowest = TRUE),
      Race_Ethnicity = factor(
        RIDRETH3,
        levels = c(1, 2, 3, 4, 6, 7),
        labels = c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian", "Other/Multi-Racial")
      )
    )
}

if (has_biopro && has_demo && ("LBXSCR" %in% names(kidney_df))) {
  kidney_df <- kidney_df %>%
    mutate(
      eGFR_2021 = calc_egfr_2021(LBXSCR, RIDAGEYR, RIAGENDR),
      Reduced_eGFR = factor(ifelse(eGFR_2021 < 60, "Yes", "No"), levels = c("No", "Yes")),
      Renal_Dysfunction = factor(
        ifelse(eGFR_2021 < 60 | URDACT >= 30, "Dysfunction", "Normal"),
        levels = c("Normal", "Dysfunction")
      ),
      BUN_Creatinine_Ratio   = LBXSBU / LBXSCR,
      UricAcid_eGFR_Ratio    = LBXSUA / eGFR_2021,
      Albumin_Globulin_Ratio = LBXSAL / LBXSGB,
      Anion_Gap              = LBXSNASI - (LBXSCLSI + LBXSC3SI)
    )
}

analytic_df <- kidney_df %>%
  filter(!is.na(URXUMA) & !is.na(URXUCR) & !is.na(URDACT))
if (has_demo) {
  analytic_df <- analytic_df %>% filter(RIDAGEYR >= 18)
}
analytic_df <- clean_sas_labels(analytic_df)

# 4. CLASSICAL BIOSTATISTICS, SURVEY WEIGHTING & ECONOMETRIC MODELS ------------

# Table 1: Baseline Cohort & Routine Serum Biochemistry Characteristics Stratified by Outcome
if (has_biopro && has_demo) {
  serum_vars <- intersect(
    c("RIDAGEYR", "Sex", "Race_Ethnicity", "URDACT", "Albuminuria_Stage", "eGFR_2021",
      "LBXSCR", "LBXSBU", "LBXSUA", "LBXSAL", "LBXSGL", "LBXSCH", "LBXSTR",
      "LBXSNASI", "LBXSKSI", "LBXSCLSI", "LBXSC3SI", "LBXSCA", "LBXSPH", "LBXSIR",
      "BUN_Creatinine_Ratio", "Anion_Gap"),
    names(analytic_df)
  )
  
  table1_baseline <- analytic_df %>%
    filter(!is.na(Renal_Dysfunction)) %>%
    dplyr::select(Renal_Dysfunction, all_of(serum_vars)) %>%
    tbl_summary(
      by = Renal_Dysfunction,
      statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p}%)"),
      digits = all_continuous() ~ 2,
      missing = "no"
    ) %>%
    add_overall() %>%
    add_p() %>%
    bold_labels()
  
  t1_df <- as_tibble(table1_baseline)
  export_table(t1_df, "Table1_Baseline_Serum_Characteristics",
               "Baseline demographic and routine serum biochemistry characteristics stratified by renal dysfunction.", "tab:baseline")
  gt::gtsave(as_gt(table1_baseline), filename = "Tables/Table1_Baseline_Serum_Characteristics.html")
} else {
  table1_baseline <- analytic_df %>%
    dplyr::select(Albuminuria_Stage, URXUMA, URXUCR, URDACT, Albumin_LLOD) %>%
    tbl_summary(by = Albuminuria_Stage, statistic = list(all_continuous() ~ "{median} ({p25}, {p75})")) %>%
    add_overall() %>%
    add_p() %>%
    bold_labels()
  
  t1_df <- as_tibble(table1_baseline)
  export_table(t1_df, "Table1_Baseline_Albuminuria",
               "Baseline urinary biomarker characteristics stratified by KDIGO albuminuria stage.", "tab:baseline")
  gt::gtsave(as_gt(table1_baseline), filename = "Tables/Table1_Baseline_Albuminuria.html")
}

# Table 2 & 3: Complex Survey-Weighted Nationally Representative Prevalence & Econometric Models
if (has_demo && all(c("WTMEC2YR", "SDMVPSU", "SDMVSTRA") %in% names(analytic_df))) {
  options(survey.lonely.psu = "adjust")
  nhanes_design <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC2YR,
    nest = TRUE,
    data = analytic_df %>% filter(!is.na(WTMEC2YR) & WTMEC2YR > 0)
  )
  
  svy_prev <- svymean(~Albuminuria_Stage, design = nhanes_design, na.rm = TRUE)
  svy_ci   <- confint(svy_prev)
  table2_svy_prev <- tibble(
    KDIGO_Stage = gsub("Albuminuria_Stage", "", names(coef(svy_prev))),
    Weighted_Prevalence_Pct = round(as.numeric(coef(svy_prev)) * 100, 2),
    SE_Pct = round(as.numeric(SE(svy_prev)) * 100, 2),
    CI_95_Lower = round(svy_ci[, 1] * 100, 2),
    CI_95_Upper = round(svy_ci[, 2] * 100, 2)
  )
  export_table(table2_svy_prev, "Table2_National_Survey_Weighted_Prevalence",
               "Nationally representative survey-weighted prevalence of KDIGO albuminuria stages.", "tab:svyprev")
  
  if (has_biopro) {
    svy_logit <- svyglm(
      I(Renal_Dysfunction == "Dysfunction") ~ RIDAGEYR + Sex + LBXSCR + LBXSBU + LBXSUA + LBXSAL + LBXSGL + LBXSKSI + LBXSC3SI,
      design = nhanes_design,
      family = quasibinomial()
    )
    table3_svy_logit <- tidy(svy_logit, conf.int = TRUE, exponentiate = TRUE) %>%
      rename(Adjusted_OR = estimate) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
    export_table(table3_svy_logit, "Table3_Survey_Weighted_Multivariable_Logit",
                 "Survey-weighted multivariable logistic regression for renal dysfunction.", "tab:svylogit")
  }
}

# Table 4: Econometric Left-Censored Tobit MLE vs. OLS
analytic_df <- analytic_df %>% mutate(cens_status = ifelse(URDUMALC == 0, 1, 0))
ols_fit   <- lm(ln_URXUMA ~ ln_URXUCR, data = analytic_df)
tobit_fit <- survreg(Surv(ln_URXUMA, cens_status, type = "left") ~ ln_URXUCR, data = analytic_df, dist = "gaussian")

table4_tobit <- bind_rows(
  tidy(ols_fit, conf.int = TRUE) %>% mutate(Estimator = "OLS (LLOD/sqrt(2) Substitution)"),
  tidy(tobit_fit, conf.int = TRUE) %>% filter(term != "Log(scale)") %>% mutate(Estimator = "Left-Censored Tobit MLE")
) %>%
  dplyr::select(Estimator, term, estimate, std.error, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
export_table(table4_tobit, "Table4_Econometric_Tobit_vs_OLS",
             "Econometric comparison of OLS vs. Left-Censored Tobit MLE for urinary albumin censoring.", "tab:tobit")

# Table 5: Quantile Regression Across Conditional Tail Percentiles
taus <- c(0.25, 0.50, 0.75, 0.90, 0.95)
table5_qr <- map_dfr(taus, function(tau_val) {
  qr_fit <- rq(log10_URXUMA ~ log10_URXUCR, tau = tau_val, data = analytic_df)
  tidy(qr_fit, se.type = "nid", conf.int = TRUE) %>% mutate(Quantile_Tau = tau_val)
}) %>%
  dplyr::select(Quantile_Tau, term, estimate, std.error, conf.low, conf.high, p.value) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
export_table(table5_qr, "Table5_Quantile_Regression_Tail_Effects",
             "Quantile regression slope heterogeneity across conditional percentiles.", "tab:qr")

# Table 6: Parametric Distribution Goodness-of-Fit
acr_vec <- analytic_df$URDACT[analytic_df$URDACT > 0]
fit_lnorm   <- fitdist(acr_vec, "lnorm")
fit_gamma   <- fitdist(acr_vec, "gamma", method = "mme")
fit_weibull <- fitdist(acr_vec, "weibull")

table6_gof <- tibble(
  Distribution = c("Log-Normal", "Gamma (MME)", "Weibull"),
  Param1 = round(c(fit_lnorm$estimate[1], fit_gamma$estimate[1], fit_weibull$estimate[1]), 4),
  Param2 = round(c(fit_lnorm$estimate[2], fit_gamma$estimate[2], fit_weibull$estimate[2]), 4),
  LogLikelihood = round(c(fit_lnorm$loglik, fit_gamma$loglik, fit_weibull$loglik), 2),
  AIC = round(c(fit_lnorm$aic, fit_gamma$aic, fit_weibull$aic), 2),
  BIC = round(c(fit_lnorm$bic, fit_gamma$bic, fit_weibull$bic), 2)
) %>% arrange(AIC)
export_table(table6_gof, "Table6_Parametric_Distribution_GOF",
             "Goodness-of-fit comparison (AIC/BIC) for continuous ACR distributions.", "tab:gof")


# 5. HETEROGENEOUS ML, DEEP NEURAL NETWORKS & ENSEMBLE SCREENING BENCHMARK ----
set.seed(2026)

if (has_biopro && has_demo) {
  candidate_predictors <- intersect(
    c("RIDAGEYR", "RIAGENDR", "LBXSCR", "LBXSBU", "LBXSUA", "LBXSAL", "LBXSGB",
      "LBXSGL", "LBXSCH", "LBXSTR", "LBXSNASI", "LBXSKSI", "LBXSCLSI", "LBXSC3SI",
      "LBXSCA", "LBXSPH", "LBXSIR", "LBXSTB", "LBXSASSI", "LBXSATSI", "LBXSAPSI",
      "LBXSGTSI", "LBXSLDS", "LBXSOSSI", "BUN_Creatinine_Ratio", "Anion_Gap"),
    names(analytic_df)
  )
  ml_df <- analytic_df %>%
    mutate(Target = factor(ifelse(URDACT >= 30, "Case", "Control"), levels = c("Control", "Case"))) %>%
    dplyr::select(Target, all_of(candidate_predictors)) %>%
    drop_na()
} else {
  ml_df <- analytic_df %>%
    mutate(
      Target = factor(ifelse(URDACT >= 30, "Case", "Control"), levels = c("Control", "Case")),
      URDUMALC_num = as.numeric(URDUMALC),
      log10_URXUMA_sq = log10_URXUMA^2,
      log10_URXUMA_cu = log10_URXUMA^3
    ) %>%
    dplyr::select(Target, URXUMA, log10_URXUMA, log10_URXUMA_sq, log10_URXUMA_cu, URDUMALC_num) %>%
    drop_na()
}

train_idx  <- createDataPartition(ml_df$Target, p = 0.80, list = FALSE)
train_data <- ml_df[train_idx, ]
test_data  <- ml_df[-train_idx, ]

# Create explicit 5-fold stratified index shared across all learners (including native XGBoost)
folds <- createFolds(train_data$Target, k = 5, returnTrain = TRUE)

cv_ctrl <- trainControl(
  method = "cv",
  number = 5,
  index = folds,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# Family 1: Classical Parametric Benchmark — Logistic Regression (GLM)
mod_glm <- train(Target ~ ., data = train_data, method = "glm", family = "binomial",
                 metric = "ROC", preProcess = c("center", "scale"), trControl = cv_ctrl)

# Family 2: Regularized Linear ML — Elastic Net (L1/L2 Penalized)
mod_enet <- train(Target ~ ., data = train_data, method = "glmnet",
                  metric = "ROC", preProcess = c("center", "scale"), tuneLength = 5, trControl = cv_ctrl)

# Family 3: Kernel Machine Learning — Support Vector Machine (Radial Basis Kernel)
mod_svm <- train(Target ~ ., data = train_data, method = "svmRadial",
                 metric = "ROC", preProcess = c("center", "scale"), tuneLength = 4, trControl = cv_ctrl)

# Family 4: Bagging Tree Ensemble — Random Forest (ranger)
p_dim <- ncol(train_data) - 1
mod_rf <- train(
  Target ~ ., data = train_data, method = "ranger", metric = "ROC",
  importance = "permutation",
  tuneGrid = expand.grid(mtry = unique(pmax(1, round(c(sqrt(p_dim), p_dim / 3)))),
                         splitrule = "gini", min.node.size = c(15, 30)),
  num.trees = 400, trControl = cv_ctrl
)

# Family 5: Gradient Boosting Ensemble — Native XGBoost (5-Fold Out-of-Fold + Full Refit)
# Avoids caret::xgbTree compatibility issues with modern xgboost (>= 1.7)
X_train_mat <- as.matrix(train_data %>% dplyr::select(-Target))
y_train_num <- ifelse(train_data$Target == "Case", 1, 0)
X_test_mat  <- as.matrix(test_data %>% dplyr::select(-Target))

xgb_params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 3,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 5
)

xgb_oof_preds <- numeric(nrow(train_data))
for (k in seq_along(folds)) {
  tr_i <- folds[[k]]
  va_i <- setdiff(seq_len(nrow(train_data)), tr_i)
  dtr <- xgb.DMatrix(data = X_train_mat[tr_i, , drop = FALSE], label = y_train_num[tr_i])
  dva <- xgb.DMatrix(data = X_train_mat[va_i, , drop = FALSE], label = y_train_num[va_i])
  bst_k <- xgb.train(params = xgb_params, data = dtr, nrounds = 150, verbose = 0)
  xgb_oof_preds[va_i] <- predict(bst_k, dva)
}

dfull_xgb <- xgb.DMatrix(data = X_train_mat, label = y_train_num)
mod_xgb_full <- xgb.train(params = xgb_params, data = dfull_xgb, nrounds = 150, verbose = 0)
xgb_test_preds <- predict(mod_xgb_full, xgb.DMatrix(data = X_test_mat))

# Family 6: Deep Learning — Multi-Layer Perceptron Neural Network (MLP)
mod_mlp <- train(
  Target ~ ., data = train_data, method = "nnet", metric = "ROC",
  preProcess = c("center", "scale"),
  tuneGrid = expand.grid(size = c(12, 24), decay = c(0.01, 0.1)),
  maxit = 400, MaxNWts = 5000, trace = FALSE, trControl = cv_ctrl
)

# Family 7: Heterogeneous Super-Learner Ensemble (Stacking RF + XGBoost + SVM + Deep MLP)
oof_df <- tibble(
  RF     = mod_rf$pred %>% arrange(rowIndex) %>% pull(Case),
  XGB    = xgb_oof_preds,
  SVM    = mod_svm$pred %>% arrange(rowIndex) %>% pull(Case),
  MLP    = mod_mlp$pred %>% arrange(rowIndex) %>% pull(Case),
  Target = train_data$Target
)
meta_stack <- glm(Target ~ RF + XGB + SVM + MLP, data = oof_df, family = "binomial")

# Held-Out Test Set Predictions
test_probs <- tibble(
  `Logistic Regression (Classical)`   = predict(mod_glm, test_data, type = "prob")$Case,
  `Elastic Net (Regularized ML)`      = predict(mod_enet, test_data, type = "prob")$Case,
  `Radial SVM (Kernel ML)`            = predict(mod_svm, test_data, type = "prob")$Case,
  `Random Forest (Bagging Ensemble)`  = predict(mod_rf, test_data, type = "prob")$Case,
  `XGBoost (Boosting Ensemble)`       = xgb_test_preds,
  `Deep Neural Network (MLP)`         = predict(mod_mlp, test_data, type = "prob")$Case
)

test_probs$`Heterogeneous Super-Ensemble (RF+XGB+SVM+MLP)` <- predict(
  meta_stack,
  newdata = tibble(
    RF  = test_probs$`Random Forest (Bagging Ensemble)`,
    XGB = test_probs$`XGBoost (Boosting Ensemble)`,
    SVM = test_probs$`Radial SVM (Kernel ML)`,
    MLP = test_probs$`Deep Neural Network (MLP)`
  ),
  type = "response"
)

# Table 7: Comprehensive Test-Set Discrimination & Calibration Benchmark
y_true_bin <- ifelse(test_data$Target == "Case", 1, 0)

table7_benchmark <- map_dfr(names(test_probs), function(m_name) {
  p_hat <- test_probs[[m_name]]
  roc_o <- roc(response = test_data$Target, predictor = p_hat, levels = c("Control", "Case"), quiet = TRUE)
  ci_o  <- ci.auc(roc_o)
  pr_o  <- pr.curve(scores.class0 = p_hat[y_true_bin == 1], scores.class1 = p_hat[y_true_bin == 0])
  opt_c <- coords(roc_o, "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity", "ppv", "npv"))
  pred_lbl <- factor(ifelse(p_hat >= opt_c$threshold[1], "Case", "Control"), levels = c("Control", "Case"))
  cm <- confusionMatrix(pred_lbl, test_data$Target, positive = "Case")
  brier <- mean((p_hat - y_true_bin)^2)
  
  tibble(
    Model_Architecture = m_name,
    ROC_AUC = round(as.numeric(auc(roc_o)), 4),
    ROC_AUC_95CI = paste0(round(ci_o[1], 4), "-", round(ci_o[3], 4)),
    PR_AUC = round(pr_o$auc.integral, 4),
    Sensitivity_Pct = round(opt_c$sensitivity[1] * 100, 2),
    Specificity_Pct = round(opt_c$specificity[1] * 100, 2),
    Balanced_Acc_Pct = round(cm$byClass["Balanced Accuracy"] * 100, 2),
    F1_Score = round(cm$byClass["F1"], 4),
    Brier_Score = round(brier, 4)
  )
}) %>% arrange(desc(ROC_AUC))

export_table(table7_benchmark, "Table7_Heterogeneous_ML_DL_Ensemble_Benchmark",
             "Out-of-sample discrimination (ROC-AUC, PR-AUC) and calibration (Brier) across model families.", "tab:mlbench")

# Table 8: Top Serum Biochemistry Feature Importance
rf_imp <- varImp(mod_rf)$importance %>%
  rownames_to_column("Biomarker_Feature") %>%
  arrange(desc(Overall)) %>%
  mutate(Rank = row_number(), Relative_Importance = round(Overall, 2)) %>%
  dplyr::select(Rank, Biomarker_Feature, Relative_Importance)
export_table(rf_imp, "Table8_Serum_Biomarker_Feature_Importance",
             "Permutation feature importance ranking of routine serum biochemistry predictors.", "tab:featimp")

# Table 9: Unsupervised Gaussian Mixture Model (GMM) Latent Renal Phenotyping
gmm_fit <- Mclust(analytic_df %>% dplyr::select(log10_URXUMA, log10_URXUCR), G = 2:4, verbose = FALSE)
analytic_df$Latent_Phenotype <- factor(paste0("Phenotype ", gmm_fit$classification))

table9_gmm <- analytic_df %>%
  group_by(Latent_Phenotype) %>%
  summarise(
    N = n(),
    Pct_Cohort = round((n() / nrow(analytic_df)) * 100, 2),
    Median_URXUMA = round(median(URXUMA), 2),
    Median_URXUCR = round(median(URXUCR), 2),
    Median_URDACT = round(median(URDACT), 2),
    Pct_Elevated_ACR = round(mean(URDACT >= 30) * 100, 2),
    .groups = "drop"
  )
export_table(table9_gmm, "Table9_Unsupervised_GMM_Phenotypes",
             "Unsupervised Gaussian Mixture Model (GMM) latent renal phenotypes.", "tab:gmm")


# 6. PUBLICATION-GRADE PLOTS (STRICTLY NO TITLES FOR LATEX \CAPTION) -----------

# Figure 1: Multi-Panel Urinary Biomarker Distributions (Tagged A, B, C)
p1_a <- ggplot(analytic_df, aes(x = URXUMA)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#3182bd", color = "white", alpha = 0.85) +
  geom_density(color = "#08519c", linewidth = 0.9) +
  scale_x_log10(labels = comma) +
  labs(x = "Urine Albumin, ug/mL (log10 scale)", y = "Density")

p1_b <- ggplot(analytic_df, aes(x = URXUCR)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#31a354", color = "white", alpha = 0.85) +
  geom_density(color = "#006d2c", linewidth = 0.9) +
  scale_x_log10(labels = comma) +
  labs(x = "Urine Creatinine, mg/dL (log10 scale)", y = "Density")

p1_c <- ggplot(analytic_df, aes(x = URDACT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 45, fill = "#756bb1", color = "white", alpha = 0.85) +
  geom_density(color = "#54278f", linewidth = 0.9) +
  geom_vline(xintercept = c(30, 300), linetype = "dashed", color = "#d95f02", linewidth = 0.9) +
  annotate("text", x = 34, y = 0.85, label = "30 mg/g", hjust = 0, size = 3.4, fontface = "bold") +
  annotate("text", x = 340, y = 0.85, label = "300 mg/g", hjust = 0, size = 3.4, fontface = "bold") +
  scale_x_log10(labels = comma) +
  labs(x = "Albumin-to-Creatinine Ratio, mg/g (log10 scale)", y = "Density")

fig1 <- (p1_a | p1_b) / p1_c + plot_annotation(tag_levels = "A")
ggsave("Plots/Fig1_Biomarker_Distributions.pdf", fig1, width = 10, height = 7.5, dpi = 300)
ggsave("Plots/Fig1_Biomarker_Distributions.png", fig1, width = 10, height = 7.5, dpi = 300)

# Figure 2: Joint 2D Hexbin Density with KDIGO Iso-Ratio Boundaries
fig2 <- ggplot(analytic_df, aes(x = URXUCR, y = URXUMA)) +
  geom_hex(bins = 45, alpha = 0.9) +
  scale_fill_viridis_c(option = "mako", trans = "log10", name = "Participant Count") +
  geom_abline(intercept = log10(0.3), slope = 1, linetype = "dashed", color = "#fdb863", linewidth = 1.1) +
  geom_abline(intercept = log10(3.0), slope = 1, linetype = "solid", color = "#e34a33", linewidth = 1.1) +
  annotate("text", x = 12, y = 4.8, label = "ACR = 30 mg/g (A2)", angle = 38, color = "#d95f02", fontface = "bold", size = 3.5) +
  annotate("text", x = 12, y = 48, label = "ACR = 300 mg/g (A3)", angle = 38, color = "#e34a33", fontface = "bold", size = 3.5) +
  scale_x_log10(labels = comma, limits = c(5, 600)) +
  scale_y_log10(labels = comma, limits = c(0.1, 10000)) +
  labs(x = "Urine Creatinine, mg/dL (log10 scale)", y = "Urine Albumin, ug/mL (log10 scale)")
ggsave("Plots/Fig2_Joint_Hexbin_KDIGO_Boundaries.pdf", fig2, width = 8.5, height = 6.5, dpi = 300)
ggsave("Plots/Fig2_Joint_Hexbin_KDIGO_Boundaries.png", fig2, width = 8.5, height = 6.5, dpi = 300)

# Figure 3: Parametric Distribution Goodness-of-Fit Overlay
x_grid <- seq(min(acr_vec), quantile(acr_vec, 0.98), length.out = 500)
gof_curves <- tibble(
  x = rep(x_grid, 3),
  density = c(
    dlnorm(x_grid, meanlog = fit_lnorm$estimate[1], sdlog = fit_lnorm$estimate[2]),
    dgamma(x_grid, shape = fit_gamma$estimate[1], rate = fit_gamma$estimate[2]),
    dweibull(x_grid, shape = fit_weibull$estimate[1], scale = fit_weibull$estimate[2])
  ),
  Distribution = rep(c("Log-Normal", "Gamma", "Weibull"), each = length(x_grid))
)

fig3 <- ggplot(analytic_df %>% filter(URDACT <= quantile(URDACT, 0.98)), aes(x = URDACT)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "grey82", color = "white") +
  geom_line(data = gof_curves, aes(x = x, y = density, color = Distribution, linetype = Distribution), linewidth = 1.1) +
  scale_color_manual(values = c("Log-Normal" = "#08519c", "Gamma" = "#e34a33", "Weibull" = "#238b45")) +
  labs(x = "Albumin-to-Creatinine Ratio, mg/g (Truncated at 98th Percentile)", y = "Probability Density")
ggsave("Plots/Fig3_Parametric_Goodness_of_Fit.pdf", fig3, width = 8.5, height = 5.5, dpi = 300)
ggsave("Plots/Fig3_Parametric_Goodness_of_Fit.png", fig3, width = 8.5, height = 5.5, dpi = 300)

# Figure 4: Quantile Regression Process Across Conditional Percentiles
qr_process <- map_dfr(seq(0.10, 0.95, by = 0.05), function(tau_val) {
  fit <- rq(log10_URXUMA ~ log10_URXUCR, tau = tau_val, data = analytic_df)
  tidy(fit, se.type = "nid", conf.int = TRUE) %>%
    filter(term == "log10_URXUCR") %>%
    mutate(tau = tau_val)
})
ols_slope <- coef(lm(log10_URXUMA ~ log10_URXUCR, data = analytic_df))["log10_URXUCR"]

fig4 <- ggplot(qr_process, aes(x = tau, y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#2b8cbe", alpha = 0.25) +
  geom_line(color = "#08519c", linewidth = 1.1) +
  geom_point(color = "#08519c", size = 2) +
  geom_hline(yintercept = ols_slope, linetype = "dashed", color = "#e34a33", linewidth = 1) +
  annotate("text", x = 0.25, y = ols_slope + 0.03, label = paste0("OLS Mean Slope = ", round(ols_slope, 3)), color = "#e34a33", fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = seq(0.1, 0.9, by = 0.1)) +
  labs(x = "Conditional Quantile (tau) of log10 Urine Albumin", y = "Regression Slope (Beta_1)")
ggsave("Plots/Fig4_Quantile_Regression_Process.pdf", fig4, width = 8.5, height = 5.5, dpi = 300)
ggsave("Plots/Fig4_Quantile_Regression_Process.png", fig4, width = 8.5, height = 5.5, dpi = 300)

# Figure 5: Multi-Model ROC-AUC Comparison
roc_plot_df <- map_dfr(names(test_probs), function(m_name) {
  r <- roc(response = test_data$Target, predictor = test_probs[[m_name]], levels = c("Control", "Case"), quiet = TRUE)
  tibble(
    Model = paste0(m_name, " (AUC = ", sprintf("%.3f", as.numeric(auc(r))), ")"),
    FPR = 1 - r$specificities,
    TPR = r$sensitivities
  )
})

fig5 <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(linewidth = 0.95) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(x = "1 - Specificity (False Positive Rate)", y = "Sensitivity (True Positive Rate)", color = "Model Architecture") +
  guides(color = guide_legend(ncol = 2))
ggsave("Plots/Fig5_Heterogeneous_Ensemble_ROC_Curves.pdf", fig5, width = 9.5, height = 6.5, dpi = 300)
ggsave("Plots/Fig5_Heterogeneous_Ensemble_ROC_Curves.png", fig5, width = 9.5, height = 6.5, dpi = 300)

# Figure 6: Probability Calibration Curves Across Model Families
cal_df <- test_probs %>%
  mutate(Observed = y_true_bin) %>%
  pivot_longer(-Observed, names_to = "Model", values_to = "Pred_Prob") %>%
  mutate(Prob_Bin = cut(Pred_Prob, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) %>%
  group_by(Model, Prob_Bin) %>%
  summarise(Mean_Pred = mean(Pred_Prob), Obs_Freq = mean(Observed), N = n(), .groups = "drop") %>%
  filter(N >= 5)

fig6 <- ggplot(cal_df, aes(x = Mean_Pred, y = Obs_Freq, color = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format()) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
  labs(x = "Mean Predicted Risk Probability", y = "Observed Empirical Proportion", color = "Model Architecture") +
  guides(color = guide_legend(ncol = 2))
ggsave("Plots/Fig6_Calibration_Curves.pdf", fig6, width = 9.5, height = 6.5, dpi = 300)
ggsave("Plots/Fig6_Calibration_Curves.png", fig6, width = 9.5, height = 6.5, dpi = 300)

# Figure 7: Top Serum Biochemistry Feature Importance Bar Chart
top_imp <- rf_imp %>% slice_head(n = 15)
fig7 <- ggplot(top_imp, aes(x = reorder(Biomarker_Feature, Relative_Importance), y = Relative_Importance)) +
  geom_col(fill = "#08519c", width = 0.7, alpha = 0.9) +
  coord_flip() +
  labs(x = "Biomarker / Clinical Predictor", y = "Relative Permutation Importance (0-100 Scale)")
ggsave("Plots/Fig7_Feature_Importance_Ranking.pdf", fig7, width = 8.5, height = 6, dpi = 300)
ggsave("Plots/Fig7_Feature_Importance_Ranking.png", fig7, width = 8.5, height = 6, dpi = 300)

# Figure 8: Unsupervised Gaussian Mixture Model (GMM) Latent Biomarker Clusters
fig8 <- ggplot(analytic_df, aes(x = URXUCR, y = URXUMA, color = Latent_Phenotype)) +
  geom_point(alpha = 0.35, size = 1.5) +
  stat_ellipse(level = 0.95, linewidth = 1.1) +
  geom_abline(intercept = log10(0.3), slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_color_brewer(palette = "Dark2", name = "GMM Latent Phenotype") +
  labs(x = "Urine Creatinine, mg/dL (log10 scale)", y = "Urine Albumin, ug/mL (log10 scale)")
ggsave("Plots/Fig8_Unsupervised_GMM_Phenotypes.pdf", fig8, width = 8.5, height = 6, dpi = 300)
ggsave("Plots/Fig8_Unsupervised_GMM_Phenotypes.png", fig8, width = 8.5, height = 6, dpi = 300)

cat("\nPipeline complete! All .tex tables and title-free .pdf/.png figures generated.\n")