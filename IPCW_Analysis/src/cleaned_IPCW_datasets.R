# Step 1: Dataset Loading and Preparation

library(survival)
library(flexsurv)
library(flexsurvcure)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("C:/Users/frifi/OneDrive/Desktop/School Mat UHasselt/Fourth Semester/Master's Thesis/SurvivalFitCrossValidation-main/Cure_data_melanoma")

# Create results directory
if (!dir.exists("ipcw_extrapolation_results")) {
  dir.create("ipcw_extrapolation_results")
}

set.seed(12345)

cat("COMPREHENSIVE IPCW EXTRAPOLATION ANALYSIS\n")
cat("Design: 3 datasets × 3 censoring levels × All models (Traditional vs IPCW)\n")
cat("IPCW Method: Kaplan-Meier censoring survival function\n")
cat("Analysis: Model-by-model comparison\n")
cat("Reference: KM RMST from full dataset\n")
cat("Metric: RMST extrapolation accuracy\n\n")

# Dataset loading function
load_selected_datasets <- function() {
  
  datasets <- list()
  
  # SEER Pancreas (Registry data)
  if (file.exists('seer_selected.csv')) {
    df_seer_selected <- read.csv('seer_selected.csv')
    seer_pancreas <- df_seer_selected %>% 
      subset(Site.recode.ICD.O.3.WHO.2008 == 'Pancreas') %>%
      rename(time = Survival.months, status = status) %>%
      filter(time > 0, !is.na(time), !is.na(status)) %>%
      mutate(time = as.numeric(time), status = as.integer(status))
    
    if (nrow(seer_pancreas) > 100) {
      datasets[["SEER_Pancreas"]] <- list(
        data = seer_pancreas, 
        type = "registry",
        description = "Pancreatic cancer registry data"
      )
      cat("SEER_Pancreas loaded successfully\n")
    } else {
      cat("SEER_Pancreas: Insufficient data\n")
    }
  } else {
    cat("SEER_Pancreas: seer_selected.csv not found\n")
  }
  
  # NSCLC Immuno therapy (Clinical trial)
  if (file.exists("IPDfromKM_immuno_with_or_without_chemo.csv")) {
    nsclc_immuno <- read.csv("IPDfromKM_immuno_with_or_without_chemo.csv") %>%
      rename(time = 1, status = 2) %>% 
      filter(time > 0) %>%
      mutate(time = as.numeric(time), status = as.integer(status))
    
    if (nrow(nsclc_immuno) > 100) {
      datasets[["NSCLC_Immunotherapy"]] <- list(
        data = nsclc_immuno, 
        type = "immunotherapy",
        description = "NSCLC immunotherapy trial"
      )
      cat("NSCLC_Immunotherapy loaded successfully\n")
    } else {
      cat("NSCLC_Immunotherapy: Insufficient data\n")
    }
  } else {
    cat("NSCLC_Immunotherapy: IPDfromKM_immuno_with_or_without_chemo.csv not found\n")
  }
  
  # SEER SCLC (Registry data)
  if (file.exists('seer_selected.csv') && exists("df_seer_selected")) {
    df_seer_lung <- df_seer_selected %>% subset(Site.recode.ICD.O.3.WHO.2008 == 'Lung and Bronchus')
    seer_sclc <- df_seer_lung %>% 
      subset(Histologic.Type.ICD.O.3 %in% c("8041", "8042", "8043", "8044", "8045")) %>%
      rename(time = Survival.months, status = status) %>%
      filter(time > 0, !is.na(time), !is.na(status)) %>%
      mutate(time = as.numeric(time), status = as.integer(status))
    
    if (nrow(seer_sclc) > 100) {
      datasets[["SEER_SCLC"]] <- list(
        data = seer_sclc, 
        type = "registry",
        description = "Small cell lung cancer registry data"
      )
      cat("SEER_SCLC loaded successfully\n")
    } else {
      cat("SEER_SCLC: Insufficient data\n")
    }
  } else {
    cat("SEER_SCLC: Data not available\n")
  }
  
  return(datasets)
}

# Load datasets
datasets <- load_selected_datasets()

# Display dataset information
cat(sprintf("\nDataset Summary: %d datasets loaded\n\n", length(datasets)))

if (length(datasets) > 0) {
  for (name in names(datasets)) {
    ds <- datasets[[name]]
    n_total <- nrow(ds$data)
    n_events <- sum(ds$data$status)
    censoring_rate <- round((1 - mean(ds$data$status)) * 100, 1)
    
    cat(sprintf("Dataset: %s\n", name))
    cat(sprintf("Description: %s\n", ds$description))
    cat(sprintf("Type: %s\n", ds$type))
    cat(sprintf("Sample size: %d patients\n", n_total))
    cat(sprintf("Events: %d (%.1f%%)\n", n_events, (n_events/n_total)*100))
    cat(sprintf("Censoring: %.1f%%\n", censoring_rate))
    cat(sprintf("Time range: %.1f - %.1f months\n", min(ds$data$time), max(ds$data$time)))
    cat(sprintf("Median time: %.1f months\n\n", median(ds$data$time)))
  }
} else {
  stop("No datasets loaded successfully. Check file paths and data availability.")
}

# Save datasets to CSV
cat("Saving datasets to CSV files:\n")

for (name in names(datasets)) {
  filename <- paste0("ipcw_extrapolation_results/", name, "_cleaned.csv")
  write.csv(datasets[[name]]$data, filename, row.names = FALSE)
  cat(sprintf("%s saved to: %s\n", name, filename))
}

# Create dataset summary
dataset_summary <- data.frame()
for (name in names(datasets)) {
  ds <- datasets[[name]]
  n_total <- nrow(ds$data)
  n_events <- sum(ds$data$status)
  censoring_rate <- round((1 - mean(ds$data$status)) * 100, 1)
  
  dataset_summary <- rbind(dataset_summary, data.frame(
    dataset_name = name,
    description = ds$description,
    type = ds$type,
    n_total = n_total,
    n_events = n_events,
    censoring_rate = censoring_rate,
    min_time = min(ds$data$time),
    max_time = max(ds$data$time),
    median_time = median(ds$data$time),
    stringsAsFactors = FALSE
  ))
}

write.csv(dataset_summary, "ipcw_extrapolation_results/dataset_summary.csv", row.names = FALSE)
cat("Dataset summary saved to: ipcw_extrapolation_results/dataset_summary.csv\n")

cat("\nFinal Dataset Summary:\n")
print(dataset_summary)

# Data quality checks
cat("\nData Quality Checks:\n")

perform_quality_checks <- function(datasets) {
  all_checks_passed <- TRUE
  
  for (name in names(datasets)) {
    cat(sprintf("\nChecking %s:\n", name))
    data <- datasets[[name]]$data
    
    required_cols <- c("time", "status")
    missing_cols <- setdiff(required_cols, colnames(data))
    if (length(missing_cols) == 0) {
      cat("  Required columns present\n")
    } else {
      cat(sprintf("  Missing columns: %s\n", paste(missing_cols, collapse = ", ")))
      all_checks_passed <- FALSE
    }
    
    if (is.numeric(data$time) && is.numeric(data$status)) {
      cat("  Correct data types\n")
    } else {
      cat("  Incorrect data types\n")
      all_checks_passed <- FALSE
    }
    
    unique_status <- unique(data$status)
    if (all(unique_status %in% c(0, 1))) {
      cat("  Status values valid\n")
    } else {
      cat(sprintf("  Invalid status values: %s\n", paste(unique_status, collapse = ", ")))
      all_checks_passed <- FALSE
    }
    
    if (all(data$time > 0) && !any(is.na(data$time))) {
      cat("  Time values valid\n")
    } else {
      cat("  Invalid time values\n")
      all_checks_passed <- FALSE
    }
    
    if (nrow(data) >= 100) {
      cat(sprintf("  Adequate sample size (n = %d)\n", nrow(data)))
    } else {
      cat(sprintf("  Small sample size (n = %d)\n", nrow(data)))
      all_checks_passed <- FALSE
    }
    
    event_rate <- mean(data$status)
    if (event_rate >= 0.1 && event_rate <= 0.95) {
      cat(sprintf("  Reasonable event rate (%.1f%%)\n", event_rate * 100))
    } else {
      cat(sprintf("  Unusual event rate (%.1f%%)\n", event_rate * 100))
    }
  }
  
  return(all_checks_passed)
}

quality_passed <- perform_quality_checks(datasets)

if (quality_passed) {
  cat("\nAll quality checks passed. Data ready for analysis.\n")
} else {
  cat("\nSome quality issues detected. Review before proceeding.\n")
}

cat("\nStep 1 Completed: Dataset Loading and Preparation\n")
cat("What was accomplished:\n")
cat("• Loaded and cleaned datasets from source files\n")
cat("• Standardized variable names (time, status)\n")
cat("• Filtered invalid observations\n")
cat("• Converted data types appropriately\n")
cat("• Saved cleaned datasets to CSV files\n")
cat("• Generated dataset summary statistics\n\n")

cat("Files created:\n")
for (name in names(datasets)) {
  cat(sprintf("• %s_cleaned.csv\n", name))
}
cat("• dataset_summary.csv\n\n")

cat("Ready for Step 2: RMST Ground Truth Calculation\n")