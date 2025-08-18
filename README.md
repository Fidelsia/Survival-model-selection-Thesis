# Survival Model Selection & IPCW Analysis

This repository contains the R code and data for my thesis on Advanced Survival Model Selection Techniques for Health Economic Evaluations.

## The contents are: 
Two main folders for two research questions:
- `Cure_model_Analysis`: Scripts and data for evaluating model selection using k-Fold Cross Validation in cure models.
- `IPCW_Analysis`: Scripts and data for testing whether Inverse Probability of Censoring Weighting (IPCW) improves model fit.

### Running the IPCW Analysis

The IPCW analysis consists of four R scripts that should be executed sequentially as given below: 

1. `cleaned_IPCW_datasets.R` – Prepares and cleans the datasets for analysis.
2. `RMST_groundTruth_calc.R` – Calculates the restricted mean survival times (RMST) as the ground truth reference.
3. `artificial_censoring_IPCW.R` – Introduces artificial censoring to simulate different censoring scenarios.
4. `IPCW_final_model_working.R` – Fits the final IPCW models using the prepared and censored datasets.
