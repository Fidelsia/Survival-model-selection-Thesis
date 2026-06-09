# Advanced Survival Model Selection Techniques for Health Economic Evaluations

# Overview
This repository contains the R code and data for my Master's thesis project investigating two approaches to improve survival model selection in health economic evaluations, a field where researchers must predict long-term patient outcomes that extend well beyond what clinical trial data can directly show.

Specifically, it examines whether k-fold cross-validation can reduce overfitting in traditional model selection, and whether inverse probability of censoring weighting (IPCW) can improve extrapolation accuracy when a large proportion of patients are lost to follow-up.

# Research Questions
## RQ1 — Cross-Validation in Cure Models
Can k-fold cross-validation improve model selection and extrapolation accuracy, particularly in datasets showing long-term survival plateaus?
Scripts and data are in the `Cure_model_Analysis` folder

## RQ2 — IPCW Under High Censoring
Does Inverse Probability of Censoring Weighting (IPCW) improve extrapolation accuracy in datasets with high censoring rates?
Scripts and data are in the `IPCW_Analysis` folder. Run the four scripts in the following order:

1. `cleaned_IPCW_datasets.R` – Prepares and cleans the datasets for analysis.
2. `RMST_groundTruth_calc.R` – Calculates the restricted mean survival times (RMST) as the ground truth reference benchmark.
3. `artificial_censoring_IPCW.R` – Introduces artificial censoring at 40%, 60%, and 80% levels to simulate high censoring scenarios. 
4. `IPCW_final_model_working.R` – Fits the final survival models with and without IPCW weighting and compares results against the ground truth.

# Key Findings

- Cross-validation improved extrapolation accuracy in approximately half of the datasets tested, though benefits varied considerably by context
- Despite clinical evidence supporting survival plateaus, flexible spline models were consistently selected over mixture cure models
- IPCW improved prediction accuracy in approximately 91% of model comparisons, with the strongest benefit observed at around 60% censoring

# Practical Implications
The findings have direct relevance for health technology assessments, particularly when working with immature trial data or patient populations where long-term cure is clinically plausible. Better model selection methods lead to more accurate cost-effectiveness decisions and ultimately better healthcare resource allocation.

# Software & Tools
R (version 4.5.0) Main packages include:

- survival
- flexsurv
- flexsurvcure
- dplyr
- ggplot2



