# Heterogeneous Machine Learning and Deep Neural Network Ensembles for Screening Renal Dysfunction from Routine Serum Biochemistry: A Nationally Representative NHANES Study

This repository contains the complete, reproducible R analysis pipeline accompanying the manuscript **"Heterogeneous Machine Learning and Deep Neural Network Ensembles for Screening Renal Dysfunction from Routine Serum Biochemistry: A Nationally Representative NHANES Study (2021–2023 Cycle)"**.

The pipeline integrates complex survey-weighted epidemiology, censored and tail-sensitive econometrics, unsupervised Gaussian Mixture Modeling (GMM), and a seven-architecture machine learning and deep learning benchmark to screen for elevated albuminuria ($\text{ACR} \ge 30\text{ mg/g}$) using strictly routine outpatient serum biochemistry panels without urinary inputs.

---

## Key Analytical Components

1. **Multi-File NHANES 2021–2023 (`_L`) Integration & Clinical Phenotyping**
   - Merges urinary albumin/creatinine (`ALB_CR_L`), standard serum biochemistry (`BIOPRO_L`), and demographics/survey weights (`DEMO_L`) by `SEQN`.
   - Computes race-free **CKD-EPI 2021 eGFR**, **KDIGO albuminuria stages (A1–A3)**, composite renal dysfunction ($\text{eGFR}_{2021} < 60\text{ mL/min/1.73 m}^2$ or $\text{ACR} \ge 30\text{ mg/g}$), and derived clinical ratios (BUN-to-creatinine ratio, uric acid-to-eGFR ratio, albumin-to-globulin ratio, and serum anion gap).

2. **Complex Survey-Weighted Biostatistics & Econometrics**
   - **Survey-Weighted Inference (`survey`):** Estimates nationally representative KDIGO albuminuria prevalence and multivariable quasibinomial logistic regression models accounting for `WTMEC2YR`, `SDMVPSU`, and `SDMVSTRA`.
   - **Left-Censored Tobit MLE (`survival`):** Corrects attenuation bias caused by lower-limit-of-detection ($\text{LLOD}/\sqrt{2}$) substitution in urinary albumin excretion[cite: 9].
   - **Conditional Quantile Regression (`quantreg`):** Evaluates albumin-to-creatinine scaling heterogeneity across conditional percentiles ($\tau \in [0.10, 0.95]$)[cite: 9].
   - **Parametric Goodness-of-Fit (`fitdistrplus`):** Benchmarks Log-Normal, Weibull, and Gamma distributions via AIC/BIC[cite: 9].
   - **Unsupervised Latent Phenotyping (`mclust`):** Identifies bivariate $(\log_{10}\text{URXUMA}, \log_{10}\text{URXUCR})$ clusters via BIC-guided Expectation-Maximization[cite: 9].

3. **Heterogeneous ML, Deep Neural Network & Stacked Super-Ensemble Benchmark**
   - Evaluates seven distinct model families using stratified 5-fold cross-validation (`seed = 2026`) and an 80/20 train-test split[cite: 9]:
     1. Classical Logistic Regression (`glm`)
     2. Regularized Elastic Net (`glmnet`)
     3. Support Vector Machine with Radial Basis Kernel (`svmRadial`)
     4. Random Forest Bagging Ensemble (`ranger`) with permutation feature importance
     5. Gradient Boosted Decision Trees (`xgboost` native out-of-fold + full refit)
     6. Deep Multi-Layer Perceptron Neural Network (`nnet`)
     7. **Heterogeneous Stacked Super-Ensemble** combining out-of-fold probability predictions from Random Forest, XGBoost, Radial SVM, and Deep MLP[cite: 9].

---

## Repository Structure

```text
├── R script for the kidney disease analysis.R   # End-to-end R pipeline
├── README.md                                      # Project documentation
├── Tables/                                        # Auto-generated CSV, HTML, and booktabs LaTeX tables
│   ├── Table1_Baseline_Serum_Characteristics.*
│   ├── Table2_National_Survey_Weighted_Prevalence.*
│   ├── Table3_Survey_Weighted_Multivariable_Logit.*
│   ├── Table4_Econometric_Tobit_vs_OLS.*
│   ├── Table5_Quantile_Regression_Tail_Effects.*
│   ├── Table6_Parametric_Distribution_GOF.*
│   ├── Table7_Heterogeneous_ML_DL_Ensemble_Benchmark.*
│   ├── Table8_Serum_Biomarker_Feature_Importance.*
│   └── Table9_Unsupervised_GMM_Phenotypes.*
└── Plots/                                         # Auto-generated publication-ready PDF & PNG figures (300 DPI)
    ├── Fig1_Biomarker_Distributions.(pdf|png)
    ├── Fig2_Joint_Hexbin_KDIGO_Boundaries.(pdf|png)
    ├── Fig3_Parametric_Goodness_of_Fit.(pdf|png)
    ├── Fig4_Quantile_Regression_Process.(pdf|png)
    ├── Fig5_Heterogeneous_Ensemble_ROC_Curves.(pdf|png)
    ├── Fig6_Calibration_Curves.(pdf|png)
    ├── Fig7_Feature_Importance_Ranking.(pdf|png)
    └── Fig8_Unsupervised_GMM_Phenotypes.(pdf|png)
