## Comprehensive Spline-Cure Analysis - Full Dataset

# ===============================================================================
#                    COMPLETE ALL-IN-ONE SPLINE-CURE ANALYSIS
#                         WITH COMPREHENSIVE RESULTS EXPORT
# ===============================================================================

library(survival)
library(flexsurvcure)
library(flexsurv)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(openxlsx)
library(jsonlite)
library(gridExtra)


# Set working directory
setwd("C:/Users/frifi/OneDrive/Desktop/School Mat UHasselt/Fourth Semester/Master's Thesis/SurvivalFitCrossValidation-main/Cure_data_melanoma")

# Create results directory
if (!dir.exists("analysis_results")) {
  dir.create("analysis_results")
}

cat("===============================================================================\n")
cat("                    COMPREHENSIVE SPLINE-CURE MODEL ANALYSIS\n")
cat("                          FULL SPLINE CONFIGURATION\n")
cat("===============================================================================\n\n")

# ===============================================================================
#                              LOAD DATASETS
# ===============================================================================

cat("Loading datasets...\n")


# Load built-in datasets
#data(gbsg, package = "survival")

#data(gbsg)
#head(gbsg)



# Load your custom datasets
df_seer_selected <- read.csv('seer_selected.csv')
df_seer_breast <- df_seer_selected %>%
  subset(Site.recode.ICD.O.3.WHO.2008 == 'Breast') %>%
  rename(time = Survival.months, status = status)

data.ovarian <- read.csv("dataOvarian1.csv") %>%
  filter(time > 0) %>%
  rename(time = time, status = event)

# Digitized KM data
ipd_ipilimumab <- read.csv("Os_green_ipili_Data (IPD) From Kaplan-Meier Survival Curve.csv") %>%
  rename(time = 1, status = 2) %>% filter(time > 0)

ipd_antiPD1 <- read.csv("From Kaplan-Meier Survival Curve_blue_AntiPD.csv") %>%
  rename(time = 1, status = 2) %>% filter(time > 0)

ipd_ipinivo <- read.csv("OS_orange_IpilliNivoData (IPD) From Kaplan-Meier.csv") %>%
  rename(time = 1, status = 2) %>% filter(time > 0)

ipd_nsclc_immuno <- read.csv("IPDfromKM_immuno_with_or_without_chemo.csv") %>%
  rename(time = 1, status = 2) %>% filter(time > 0)

ipd_nsclc_chemo <- read.csv("IPDfromKM _chemo.csv") %>%
  rename(time = 1, status = 2) %>% filter(time > 0)

# Prepare GBSG dataset
gbsg_cleaned <- gbsg %>%
  rename(time = rfstime, status = status) %>%
  filter(time > 0) %>%
  mutate(status = as.integer(status))

head(gbsg_cleaned)

cat("Datasets loaded successfully!\n\n")

# ===============================================================================
#                         MODEL FITTING FUNCTIONS
# ===============================================================================

fit_enhanced_cure_models <- function(data, time_col = 'time', status_col = 'status') {
  
  # Create survival formula
  surv_formula <- as.formula(paste("Surv(", time_col, ",", status_col, ") ~ 1"))
  
  # Initialize storage
  fitted_models <- list()
  model_names <- character()
  cure_fractions <- numeric()
  
  #cat("  Fitting models...")
  
  # =============================
  # PARAMETRIC CURE MODELS
  # =============================
  parametric_dists <- list(
    "Cure_Weibull" = "weibull",
    "Cure_Lognormal" = "lnorm", 
    "Cure_Loglogistic" = "llogis",
    "Cure_Exponential" = "exp",
    "Cure_Gompertz" = "gompertz"
  )
  
  for (dist_name in names(parametric_dists)) {
    flex_dist <- parametric_dists[[dist_name]]
    
    model_fit <- tryCatch({
      suppressWarnings(suppressMessages(
        flexsurvcure(
          formula = surv_formula,
          data = data,
          dist = flex_dist,
          cureform = ~1
        )
      ))
    }, error = function(e) NULL)
    
    if (!is.null(model_fit)) {
      fitted_models[[dist_name]] <- model_fit
      model_names <- c(model_names, dist_name)
      
      # Extract cure fraction
      if ("theta" %in% names(model_fit$coefficients)) {
        cure_fractions[dist_name] <- plogis(model_fit$coefficients[["theta"]])
      } else {
        cure_fractions[dist_name] <- NA
      }
    }
  }
  
  # =============================
  # COMPREHENSIVE SPLINE MODELS
  # =============================
  # Spline configurations: k=1,2,3,4 with hazard, odds, and normal scales
  spline_configs <- expand.grid(
    k = 1:4,
    scale = c("hazard", "odds", "normal"),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:nrow(spline_configs)) {
    k_val <- spline_configs$k[i]
    scale_val <- spline_configs$scale[i]
    model_name <- paste0("Spline_k", k_val, "_", scale_val)
    
    spline_model <- tryCatch({
      suppressWarnings(suppressMessages(
        flexsurvspline(
          formula = surv_formula,
          data = data,
          k = k_val,
          scale = scale_val
        )
      ))
    }, error = function(e) NULL)
    
    if (!is.null(spline_model)) {
      fitted_models[[model_name]] <- spline_model
      model_names <- c(model_names, model_name)
      cure_fractions[model_name] <- 0  # Standard models assume no cure
    }
  }
  
  
  # Calculate AIC/BIC
  if (length(fitted_models) == 0) {
    return(list(
      models = list(),
      model_names = character(0),
      performance_metrics = data.frame(aic = numeric(0), bic = numeric(0)),
      cure_fractions = numeric(0)
    ))
  }
  
  aic_values <- sapply(fitted_models, function(m) {
    tryCatch(AIC(m), error = function(e) Inf)
  })
  
  bic_values <- sapply(fitted_models, function(m) {
    tryCatch(BIC(m), error = function(e) Inf)
  })
  
  return(list(
    models = fitted_models,
    model_names = model_names,
    performance_metrics = data.frame(aic = aic_values, bic = bic_values, row.names = model_names),
    cure_fractions = cure_fractions
  ))
}

# ===============================================================================
#                           VALIDATION AND CV FUNCTIONS
# ===============================================================================

assess_models_on_validation <- function(fitted_models, validation_data, 
                                        time_col = 'time', status_col = 'status') {
  
  n_models <- length(fitted_models)
  metrics <- data.frame(
    aic = rep(Inf, n_models),
    bic = rep(Inf, n_models),
    row.names = names(fitted_models)
  )
  
  for (i in seq_along(fitted_models)) {
    model <- fitted_models[[i]]
    
    if (is.null(model)) next
    
    tryCatch({
      log_lik <- logLik(model, newdata = validation_data)
      n_params <- length(model$coefficients)
      n_obs <- nrow(validation_data)
      
      if (is.finite(log_lik) && is.finite(n_params)) {
        metrics$aic[i] <- -2 * log_lik + 2 * n_params
        metrics$bic[i] <- -2 * log_lik + log(n_obs) * n_params
      }
    }, error = function(e) {
      # Keep Inf values
    })
  }
  
  return(metrics)
}

cross_validate_models <- function(data, k_folds = 10, time_col = 'time', status_col = 'status') {
  
  # Shuffle and create folds
  shuffled_data <- data[sample(nrow(data)), ]
  fold_size <- floor(nrow(shuffled_data) / k_folds)
  
  # Get model template from first fit
  temp_size <- min(nrow(shuffled_data) - fold_size, 200)
  temp_data <- shuffled_data[1:temp_size, ]
  
  initial_fit <- suppressWarnings(suppressMessages(
    fit_enhanced_cure_models(temp_data, time_col, status_col)
  ))
  
  template_names <- initial_fit$model_names
  
  if (length(template_names) == 0) {
    return(data.frame(aic = numeric(0), bic = numeric(0)))
  }
  
  # Storage for fold metrics
  fold_metrics <- list()
  
  cat("    CV:")
  
  # Perform CV
  for (fold in 1:k_folds) {
    
    # Create train/test split
    test_start <- (fold - 1) * fold_size + 1
    test_end <- min(fold * fold_size, nrow(shuffled_data))
    test_indices <- test_start:test_end
    
    train_data <- shuffled_data[-test_indices, ]
    test_data <- shuffled_data[test_indices, ]
    
    # Fit models on training data
    train_fits <- suppressWarnings(suppressMessages(
      fit_enhanced_cure_models(train_data, time_col, status_col)
    ))
    
    # Initialize fold metrics
    fold_metric <- data.frame(
      aic = rep(Inf, length(template_names)),
      bic = rep(Inf, length(template_names)),
      row.names = template_names
    )
    
    # Evaluate on test data
    if (length(train_fits$models) > 0) {
      test_metrics <- assess_models_on_validation(
        train_fits$models, test_data, time_col, status_col
      )
      
      # Update metrics for common models
      common_models <- intersect(rownames(fold_metric), rownames(test_metrics))
      if (length(common_models) > 0) {
        fold_metric[common_models, ] <- test_metrics[common_models, ]
      }
    }
    
    fold_metrics[[fold]] <- fold_metric
    
    if (fold %% 2 == 0) cat(".")
  }
  
  cat(" completed\n")
  
  # Average across folds
  if (length(fold_metrics) > 0) {
    valid_folds <- fold_metrics[sapply(fold_metrics, function(x) {
      is.data.frame(x) && nrow(x) == length(template_names)
    })]
    
    if (length(valid_folds) > 0) {
      averaged_metrics <- Reduce("+", valid_folds) / length(valid_folds)
      return(averaged_metrics)
    }
  }
  
  # Return empty if failed
  return(data.frame(
    aic = rep(Inf, length(template_names)),
    bic = rep(Inf, length(template_names)),
    row.names = template_names
  ))
}

# ===============================================================================
#                              HELPER FUNCTIONS
# ===============================================================================

calculate_model_rmst <- function(model, time_horizon) {
  if (is.null(model)) return(NA)
  
  rmst_value <- tryCatch({
    summary_result <- summary(model, t = time_horizon, type = "rmst")
    
    if (is.list(summary_result) && length(summary_result) > 0 && "est" %in% names(summary_result[[1]])) {
      summary_result[[1]]$est
    } else if (is.data.frame(summary_result) && "est" %in% colnames(summary_result)) {
      summary_result$est[1]
    } else {
      NA
    }
  }, error = function(e) NA)
  
  return(rmst_value)
}

calculate_rmst_km <- function(km_fit, tau) {
  if (!inherits(km_fit, "survfit")) stop("km_fit must be a survfit object.")
  if (length(km_fit$time) == 0) return(0)
  
  plot_times <- unique(c(0, km_fit$time[km_fit$time <= tau], tau))
  plot_times <- plot_times[order(plot_times)]
  
  s_summary <- summary(km_fit, times = plot_times, extend = TRUE)
  plot_surv <- s_summary$surv
  
  if(plot_times[1] == 0 && (length(plot_surv) < length(plot_times) || plot_surv[1] !=1) ){
    temp_times <- c(0, km_fit$time)
    temp_surv <- c(1, km_fit$surv)
    plot_surv_approx <- approx(temp_times, temp_surv, xout = plot_times, method="constant", f=0, rule=2)$y
    plot_surv <- plot_surv_approx
  }
  
  if (length(plot_times) < 2) return(0)
  
  idx <- 1:(length(plot_times) - 1)
  rmst <- sum( (plot_surv[idx] + plot_surv[idx+1]) / 2 * diff(plot_times) )
  return(rmst)
}

apply_artificial_censoring <- function(data, time_col = 'time', status_col = 'status', 
                                       sample_size = 250, survival_cutoff = 0.5) {
  
  data_work <- data
  data_work$time_work <- data_work[[time_col]]
  data_work$status_work <- data_work[[status_col]]
  
  if (nrow(data_work) < 1) return(data_work[0,])
  
  n_sample <- min(nrow(data_work), sample_size)
  if (n_sample < 1) return(data_work[0,])
  
  sampled_data <- data_work[sample(nrow(data_work), n_sample), ]
  if (nrow(sampled_data) == 0) return(sampled_data)
  
  km_fit <- suppressWarnings(
    survfit(Surv(time_work, status_work) ~ 1, data = sampled_data)
  )
  
  censor_time <- Inf
  if (any(km_fit$surv < survival_cutoff)) {
    cutoff_indices <- which(km_fit$surv < survival_cutoff)
    if (length(cutoff_indices) > 0) {
      censor_time <- km_fit$time[min(cutoff_indices)]
    }
  }
  
  if (is.finite(censor_time)) {
    sampled_data$status_work <- ifelse(sampled_data$time_work > censor_time, 0, sampled_data$status_work)
    sampled_data$time_work <- pmin(sampled_data$time_work, censor_time)
  }
  
  sampled_data[[time_col]] <- sampled_data$time_work
  sampled_data[[status_col]] <- sampled_data$status_work
  sampled_data <- sampled_data[, !names(sampled_data) %in% c("time_work", "status_work")]
  
  return(sampled_data)
}

# ===============================================================================
#                              MAIN ANALYSIS FUNCTION
# ===============================================================================

run_dataset_analysis <- function(dataset, dataset_name, time_col = 'time', status_col = 'status',
                                 sample_size = 250, n_simulations = 100, cv_folds = 10) {
  
  cat("\n=== ANALYZING DATASET:", dataset_name, "===\n")
  cat("Original size:", nrow(dataset), "observations\n")
  cat("Configuration: Sample size =", sample_size, "| Simulations =", n_simulations, "| CV folds =", cv_folds, "\n")
  
  # Calculate gold standard
  max_time <- max(dataset[[time_col]], na.rm = TRUE)
  full_km_fit <- suppressWarnings(survfit(as.formula(paste("Surv(", time_col, ",", status_col, ") ~ 1")), 
                                          data = dataset))
  gold_rmst <- calculate_rmst_km(full_km_fit, max_time)
  
  cat("Gold standard RMST:", round(gold_rmst, 2), "\n")
  
  # Storage for results
  results <- list(
    simulation_id = 1:n_simulations,
    rmst_traditional_aic = numeric(n_simulations),
    rmst_cv_aic = numeric(n_simulations),
    rmst_traditional_bic = numeric(n_simulations),
    rmst_cv_bic = numeric(n_simulations),
    pi_traditional_aic = numeric(n_simulations),
    pi_cv_aic = numeric(n_simulations),
    pi_traditional_bic = numeric(n_simulations),
    pi_cv_bic = numeric(n_simulations),
    model_traditional_aic = character(n_simulations),
    model_cv_aic = character(n_simulations),
    model_traditional_bic = character(n_simulations),
    model_cv_bic = character(n_simulations)
  )
  
  # Initialize with defaults
  results$rmst_traditional_aic[] <- NA
  results$rmst_cv_aic[] <- NA
  results$rmst_traditional_bic[] <- NA
  results$rmst_cv_bic[] <- NA
  results$pi_traditional_aic[] <- NA
  results$pi_cv_aic[] <- NA
  results$pi_traditional_bic[] <- NA
  results$pi_cv_bic[] <- NA
  results$model_traditional_aic[] <- "NOT_RUN"
  results$model_cv_aic[] <- "NOT_RUN"
  results$model_traditional_bic[] <- "NOT_RUN"
  results$model_cv_bic[] <- "NOT_RUN"
  
  # Progress tracking
  progress_marks <- seq(10, n_simulations, by = 10)
  
  # Run simulations
  for (sim in 1:n_simulations) {
    cat("Simulation", sim, "of", n_simulations, "...")
    
    # Generate sample
    sample_data <- apply_artificial_censoring(dataset, time_col, status_col, sample_size, 0.5)
    
    if (nrow(sample_data) < cv_folds * 2) {
      cat(" skipped (sample too small)\n")
      next
    }
    
    # Fit models
    sample_fits <- suppressWarnings(suppressMessages(
      fit_enhanced_cure_models(sample_data, time_col, status_col)
    ))
    
    if (length(sample_fits$models) == 0) {
      cat(" skipped (no models converged)\n")
      next
    }
    
    # Traditional selection
    traditional_metrics <- sample_fits$performance_metrics
    
    if (nrow(traditional_metrics) > 0 && any(is.finite(traditional_metrics$aic))) {
      best_aic_idx <- which.min(traditional_metrics$aic)
      best_aic_model <- sample_fits$model_names[best_aic_idx]
      results$model_traditional_aic[sim] <- best_aic_model
      results$rmst_traditional_aic[sim] <- calculate_model_rmst(sample_fits$models[[best_aic_model]], max_time)
      results$pi_traditional_aic[sim] <- sample_fits$cure_fractions[best_aic_model]
    }
    
    if (nrow(traditional_metrics) > 0 && any(is.finite(traditional_metrics$bic))) {
      best_bic_idx <- which.min(traditional_metrics$bic)
      best_bic_model <- sample_fits$model_names[best_bic_idx]
      results$model_traditional_bic[sim] <- best_bic_model
      results$rmst_traditional_bic[sim] <- calculate_model_rmst(sample_fits$models[[best_bic_model]], max_time)
      results$pi_traditional_bic[sim] <- sample_fits$cure_fractions[best_bic_model]
    }
    
    # Cross-validation selection
    cv_metrics <- suppressWarnings(suppressMessages(
      cross_validate_models(sample_data, cv_folds, time_col, status_col)
    ))
    
    if (nrow(cv_metrics) > 0 && any(is.finite(cv_metrics$aic))) {
      best_cv_aic_idx <- which.min(cv_metrics$aic)
      best_cv_aic_model <- rownames(cv_metrics)[best_cv_aic_idx]
      
      if (best_cv_aic_model %in% names(sample_fits$models)) {
        results$model_cv_aic[sim] <- best_cv_aic_model
        results$rmst_cv_aic[sim] <- calculate_model_rmst(sample_fits$models[[best_cv_aic_model]], max_time)
        results$pi_cv_aic[sim] <- sample_fits$cure_fractions[best_cv_aic_model]
      }
    }
    
    if (nrow(cv_metrics) > 0 && any(is.finite(cv_metrics$bic))) {
      best_cv_bic_idx <- which.min(cv_metrics$bic)
      best_cv_bic_model <- rownames(cv_metrics)[best_cv_bic_idx]
      
      if (best_cv_bic_model %in% names(sample_fits$models)) {
        results$model_cv_bic[sim] <- best_cv_bic_model
        results$rmst_cv_bic[sim] <- calculate_model_rmst(sample_fits$models[[best_cv_bic_model]], max_time)
        results$pi_cv_bic[sim] <- sample_fits$cure_fractions[best_cv_bic_model]
      }
    }
    
    # Progress indication
    if (sim %in% progress_marks) {
      cat(" [", round(sim/n_simulations*100), "% complete]")
    }
    cat("\n")
  }
  
  # Calculate summary
  rmst_error_trad_aic <- mean(abs(results$rmst_traditional_aic - gold_rmst), na.rm = TRUE)
  rmst_error_cv_aic <- mean(abs(results$rmst_cv_aic - gold_rmst), na.rm = TRUE)
  rmst_error_trad_bic <- mean(abs(results$rmst_traditional_bic - gold_rmst), na.rm = TRUE)
  rmst_error_cv_bic <- mean(abs(results$rmst_cv_bic - gold_rmst), na.rm = TRUE)
  
  cat("\n=== RESULTS FOR", dataset_name, "===\n")
  cat("RMST Errors:\n")
  cat("Traditional AIC:", round(rmst_error_trad_aic, 2), "\n")
  cat("CV AIC:         ", round(rmst_error_cv_aic, 2), "\n")
  cat("Traditional BIC:", round(rmst_error_trad_bic, 2), "\n")
  cat("CV BIC:         ", round(rmst_error_cv_bic, 2), "\n")
  
  return(list(
    dataset = dataset_name,
    gold_rmst = gold_rmst,
    rmst_error_trad_aic = rmst_error_trad_aic,
    rmst_error_cv_aic = rmst_error_cv_aic,
    rmst_error_trad_bic = rmst_error_trad_bic,
    rmst_error_cv_bic = rmst_error_cv_bic,
    detailed_results = as.data.frame(results)
  ))
}

# ===============================================================================
#                              RUN MAIN ANALYSIS
# ===============================================================================

# Define your datasets
datasets_list <- list(
  list(name = "OS_Ipilimumab", data = ipd_ipilimumab, time_var = "time", status_var = "status"),
  list(name = "OS_AntiPD1", data = ipd_antiPD1, time_var = "time", status_var = "status"),
  list(name = "OS_IpiNivo", data = ipd_ipinivo, time_var = "time", status_var = "status"),
  list(name = "NSCLC_Immunotherapy", data = ipd_nsclc_immuno, time_var = "time", status_var = "status"),
  list(name = "NSCLC_Chemotherapy", data = ipd_nsclc_chemo, time_var = "time", status_var = "status"),
  list(name = "SEER_Breast", data = df_seer_breast, time_var = "time", status_var = "status"),
  list(name = "Ovarian", data = data.ovarian, time_var = "time", status_var = "status"),
  list(name = "GBSG", data = gbsg_cleaned, time_var = "time", status_var = "status")
)


# Run analysis on all datasets
all_results <- list()

for (ds_info in datasets_list) {
  # Clean data
  cleaned_data <- ds_info$data %>%
    filter(!is.na(.data[[ds_info$time_var]]) & !is.na(.data[[ds_info$status_var]])) %>%
    filter(.data[[ds_info$time_var]] > 0) %>%
    mutate(status = as.integer(.data[[ds_info$status_var]]))
  
  if (nrow(cleaned_data) < 100) {
    cat("Skipping", ds_info$name, "- too small (", nrow(cleaned_data), "rows)\n")
    next
  }
  
  # Run comprehensive analysis
  result <- run_dataset_analysis(
    dataset = cleaned_data,
    dataset_name = ds_info$name,
    time_col = ds_info$time_var,
    status_col = "status",
    sample_size = 250,
    n_simulations = 100,  # Increased from 10 to 100
    cv_folds = 10         # Increased from 5 to 10
  )
  
  all_results[[ds_info$name]] <- result
}



# ===============================================================================
#                              COMPREHENSIVE RESULTS EXPORT
# ===============================================================================

cat("\n===============================================================================\n")
cat("                              EXPORTING RESULTS\n")
cat("===============================================================================\n")

# 1. SAVE MASTER RDS FILE
cat("1. Saving master RDS file...\n")
saveRDS(all_results, "analysis_results/master_results.rds")

# 2. SUMMARY DATA FRAME
cat("2. Creating summary tables...\n")
summary_df <- data.frame(
  Dataset = character(),
  Gold_Standard_RMST = numeric(),
  Traditional_AIC_Error = numeric(),
  CV_AIC_Error = numeric(),
  AIC_Improvement_Percent = numeric(),
  Traditional_BIC_Error = numeric(),
  CV_BIC_Error = numeric(),
  BIC_Improvement_Percent = numeric(),
  stringsAsFactors = FALSE
)

for (result in all_results) {
  aic_improvement <- round((result$rmst_error_trad_aic - result$rmst_error_cv_aic) / result$rmst_error_trad_aic * 100, 1)
  bic_improvement <- round((result$rmst_error_trad_bic - result$rmst_error_cv_bic) / result$rmst_error_trad_bic * 100, 1)
  
  summary_df <- rbind(summary_df, data.frame(
    Dataset = result$dataset,
    Gold_Standard_RMST = round(result$gold_rmst, 2),
    Traditional_AIC_Error = round(result$rmst_error_trad_aic, 2),
    CV_AIC_Error = round(result$rmst_error_cv_aic, 2),
    AIC_Improvement_Percent = aic_improvement,
    Traditional_BIC_Error = round(result$rmst_error_trad_bic, 2),
    CV_BIC_Error = round(result$rmst_error_cv_bic, 2),
    BIC_Improvement_Percent = bic_improvement,
    stringsAsFactors = FALSE
  ))
}

# 3. SAVE CSV FILES
cat("3. Saving CSV files...\n")
write.csv(summary_df, "analysis_results/summary_table.csv", row.names = FALSE)

# Detailed results saved for each dataset
for (name in names(all_results)) {
  if (!is.null(all_results[[name]]$detailed_results)) {
    detailed_df <- all_results[[name]]$detailed_results
    write.csv(detailed_df, paste0("analysis_results/detailed_", name, ".csv"), row.names = FALSE)
  }
}

# 4. EXCEL WORKBOOK
cat("4. Creating Excel workbook...\n")
wb <- createWorkbook()

# Summary sheet
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_df)

# Summary sheet
headerStyle <- createStyle(fontSize = 12, fontColour = "white", fgFill = "#4472C4", 
                           halign = "center", valign = "center", textDecoration = "bold")
addStyle(wb, "Summary", headerStyle, rows = 1, cols = 1:ncol(summary_df), gridExpand = TRUE)

# Detailed sheets for each dataset
for (name in names(all_results)) {
  if (!is.null(all_results[[name]]$detailed_results)) {
    addWorksheet(wb, name)
    writeData(wb, name, all_results[[name]]$detailed_results)
    
    # Style headers
    addStyle(wb, name, headerStyle, rows = 1, cols = 1:ncol(all_results[[name]]$detailed_results), gridExpand = TRUE)
  }
}

saveWorkbook(wb, "analysis_results/complete_results.xlsx", overwrite = TRUE)

# 5. SAVE JSON FILE 
cat("5. Saving JSON file...\n")
write_json(all_results, "analysis_results/results.json", pretty = TRUE)

# 6. CREATE VISUALIZATIONS
cat("6. Creating visualizations...\n")

# Custom theme
clean_theme <- theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_line(color = "gray95", size = 0.3),
    text = element_text(color = "black", family = "Arial"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 11, face = "bold"),
    plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "gray30", size = 11, hjust = 0.5),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray95", color = "gray80"),
    strip.text = element_text(color = "black", face = "bold")
  )

# Plot 1: RMST Error Comparison (AIC)
p1 <- ggplot(summary_df, aes(x = reorder(Dataset, Traditional_AIC_Error))) +
  geom_col(aes(y = Traditional_AIC_Error, fill = "Traditional AIC"), 
           alpha = 0.8, width = 0.4, position = position_nudge(x = -0.2)) +
  geom_col(aes(y = CV_AIC_Error, fill = "CV AIC"), 
           alpha = 0.8, width = 0.4, position = position_nudge(x = 0.2)) +
  scale_fill_manual(values = c("Traditional AIC" = "#D32F2F", "CV AIC" = "#1976D2")) +
  labs(title = "RMST Prediction Error: Traditional vs Cross-Validation (AIC)",
       subtitle = "Lower is better - CV compared to traditional selection",
       x = "Dataset", 
       y = "Mean Absolute RMST Error", 
       fill = "Method") +
  clean_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "bottom"
  ) +
  coord_flip() +
  # Add value labels on bars
  geom_text(aes(y = Traditional_AIC_Error, label = round(Traditional_AIC_Error, 1)), 
            position = position_nudge(x = -0.2), hjust = -0.1, size = 3, color = "black") +
  geom_text(aes(y = CV_AIC_Error, label = round(CV_AIC_Error, 1)), 
            position = position_nudge(x = 0.2), hjust = -0.1, size = 3, color = "black")

ggsave("analysis_results/rmst_error_comparison_aic.png", p1, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 2: RMST Error Comparison (BIC)
p2 <- ggplot(summary_df, aes(x = reorder(Dataset, Traditional_BIC_Error))) +
  geom_col(aes(y = Traditional_BIC_Error, fill = "Traditional BIC"), 
           alpha = 0.8, width = 0.4, position = position_nudge(x = -0.2)) +
  geom_col(aes(y = CV_BIC_Error, fill = "CV BIC"), 
           alpha = 0.8, width = 0.4, position = position_nudge(x = 0.2)) +
  scale_fill_manual(values = c("Traditional BIC" = "#F57C00", "CV BIC" = "#388E3C")) +
  labs(title = "RMST Prediction Error: Traditional vs Cross-Validation (BIC)",
       subtitle = "Lower is better - CV performance varies more with BIC",
       x = "Dataset", 
       y = "Mean Absolute RMST Error", 
       fill = "Method") +
  clean_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "bottom"
  ) +
  coord_flip() +
  # Add value labels on bars
  geom_text(aes(y = Traditional_BIC_Error, label = round(Traditional_BIC_Error, 1)), 
            position = position_nudge(x = -0.2), hjust = -0.1, size = 3, color = "black") +
  geom_text(aes(y = CV_BIC_Error, label = round(CV_BIC_Error, 1)), 
            position = position_nudge(x = 0.2), hjust = -0.1, size = 3, color = "black")

ggsave("analysis_results/rmst_error_comparison_bic.png", p2, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 3: Improvement Percentage
improvement_df <- summary_df %>%
  select(Dataset, AIC_Improvement_Percent, BIC_Improvement_Percent) %>%
  pivot_longer(cols = c(AIC_Improvement_Percent, BIC_Improvement_Percent), 
               names_to = "Criterion", values_to = "Improvement") %>%
  mutate(Criterion = ifelse(Criterion == "AIC_Improvement_Percent", "AIC", "BIC"))

p3 <- ggplot(improvement_df, aes(x = reorder(Dataset, Improvement), y = Improvement, fill = Criterion)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40", size = 1) +
  scale_fill_manual(values = c("AIC" = "#1976D2", "BIC" = "#388E3C")) +
  labs(title = "Cross-Validation Improvement Over Traditional Selection",
       subtitle = "Positive values indicate CV outperforms traditional methods",
       x = "Dataset", 
       y = "Improvement (%)", 
       fill = "Criterion") +
  clean_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "bottom"
  ) +
  coord_flip() +
  # Add value labels on bars
  geom_text(aes(label = paste0(round(Improvement, 1), "%")), 
            position = position_dodge(width = 0.8), hjust = ifelse(improvement_df$Improvement >= 0, -0.1, 1.1), 
            size = 3, color = "black")

ggsave("analysis_results/improvement_percentage.png", p3, 
       width = 12, height = 8, dpi = 300, bg = "white")

# Plot 4: Model Selection Frequency
# Create model selection frequency analysis
model_freq_data <- data.frame()

for (name in names(all_results)) {
  detailed <- all_results[[name]]$detailed_results
  
  if (!is.null(detailed) && nrow(detailed) > 0) {
    # Count model selections
    trad_aic_freq <- table(detailed$model_traditional_aic)
    cv_aic_freq <- table(detailed$model_cv_aic)
    
    # Get most frequent models
    if (length(trad_aic_freq) > 0) {
      most_freq_trad <- names(sort(trad_aic_freq, decreasing = TRUE))[1]
      freq_trad <- as.numeric(trad_aic_freq[most_freq_trad])
    } else {
      most_freq_trad <- "None"
      freq_trad <- 0
    }
    
    if (length(cv_aic_freq) > 0) {
      most_freq_cv <- names(sort(cv_aic_freq, decreasing = TRUE))[1]
      freq_cv <- as.numeric(cv_aic_freq[most_freq_cv])
    } else {
      most_freq_cv <- "None"
      freq_cv <- 0
    }
    
    model_freq_data <- rbind(model_freq_data, data.frame(
      Dataset = name,
      Traditional_Model = most_freq_trad,
      Traditional_Freq = freq_trad,
      CV_Model = most_freq_cv,
      CV_Freq = freq_cv,
      Total_Sims = nrow(detailed)
    ))
  }
}

# Save model frequency analysis
write.csv(model_freq_data, "analysis_results/model_selection_frequency.csv", row.names = FALSE)

# ===============================================================================
#                        DISPLAY MODEL SELECTION FREQUENCY TABLE
# ===============================================================================

cat("\n===============================================================================\n")
cat("                        MODEL SELECTION FREQUENCY ANALYSIS\n")
cat("===============================================================================\n\n")

cat("This table shows which models were most frequently selected by each method:\n\n")

# Display the model frequency table in a formatted way
if (nrow(model_freq_data) > 0) {
  
  # Calculate percentages for better interpretation
  model_freq_display <- model_freq_data %>%
    mutate(
      Traditional_Percent = round((Traditional_Freq / Total_Sims) * 100, 1),
      CV_Percent = round((CV_Freq / Total_Sims) * 100, 1)
    ) %>%
    select(Dataset, Traditional_Model, Traditional_Percent, CV_Model, CV_Percent, Total_Sims)
  
  # Print formatted table
  cat(sprintf("%-20s %-20s %-8s %-20s %-8s %-8s\n", 
              "Dataset", "Traditional Model", "Trad%", "CV Model", "CV%", "N_Sims"))
  cat(paste(rep("=", 85), collapse = ""), "\n")
  
  for (i in 1:nrow(model_freq_display)) {
    row <- model_freq_display[i, ]
    cat(sprintf("%-20s %-20s %-8s %-20s %-8s %-8s\n",
                substr(row$Dataset, 1, 19),           # Truncate if needed
                substr(row$Traditional_Model, 1, 19), # Truncate if needed
                paste0(row$Traditional_Percent, "%"),
                substr(row$CV_Model, 1, 19),          # Truncate if needed
                paste0(row$CV_Percent, "%"),
                row$Total_Sims))
  }
  
  cat("\n")
  
  # Summary statistics
  cat("SUMMARY OF MODEL SELECTION PATTERNS:\n")
  cat("------------------------------------\n")
  
  # Count cure vs non-cure models
  cure_models <- c("Cure_Weibull", "Cure_Lognormal", "Cure_Loglogistic", "Cure_Exponential", "Cure_Gompertz")
  spline_models <- c("Spline_k1_hazard", "Spline_k1_odds", "Spline_k1_normal",
                     "Spline_k2_hazard", "Spline_k2_odds", "Spline_k2_normal",
                     "Spline_k3_hazard", "Spline_k3_odds", "Spline_k3_normal",
                     "Spline_k4_hazard", "Spline_k4_odds", "Spline_k4_normal")
  
  # Traditional method preferences
  trad_cure_selections <- sum(model_freq_display$Traditional_Model %in% cure_models, na.rm = TRUE)
  trad_spline_selections <- sum(model_freq_display$Traditional_Model %in% spline_models, na.rm = TRUE)
  
  # CV method preferences  
  cv_cure_selections <- sum(model_freq_display$CV_Model %in% cure_models, na.rm = TRUE)
  cv_spline_selections <- sum(model_freq_display$CV_Model %in% spline_models, na.rm = TRUE)
  
  cat("Traditional Method:\n")
  cat("  - Datasets where cure models were most frequent:", trad_cure_selections, "/", nrow(model_freq_display), "\n")
  cat("  - Datasets where spline models were most frequent:", trad_spline_selections, "/", nrow(model_freq_display), "\n")
  
  cat("\nCross-Validation Method:\n")
  cat("  - Datasets where cure models were most frequent:", cv_cure_selections, "/", nrow(model_freq_display), "\n")
  cat("  - Datasets where spline models were most frequent:", cv_spline_selections, "/", nrow(model_freq_display), "\n")
  
  # Agreement analysis
  agreement_count <- sum(model_freq_display$Traditional_Model == model_freq_display$CV_Model, na.rm = TRUE)
  cat("\nMethod Agreement:\n")
  cat("  - Datasets where both methods selected the same model type:", agreement_count, "/", nrow(model_freq_display), 
      "(", round(agreement_count/nrow(model_freq_display)*100, 1), "%)\n")
  
  # Most popular models overall
  all_trad_models <- model_freq_display$Traditional_Model[!is.na(model_freq_display$Traditional_Model) & 
                                                            model_freq_display$Traditional_Model != "NOT_RUN"]
  all_cv_models <- model_freq_display$CV_Model[!is.na(model_freq_display$CV_Model) & 
                                                 model_freq_display$CV_Model != "NOT_RUN"]
  
  if (length(all_trad_models) > 0) {
    most_popular_trad <- names(sort(table(all_trad_models), decreasing = TRUE))[1]
    cat("\nMost frequently selected by Traditional method:", most_popular_trad, "\n")
  }
  
  if (length(all_cv_models) > 0) {
    most_popular_cv <- names(sort(table(all_cv_models), decreasing = TRUE))[1]
    cat("Most frequently selected by Cross-Validation:", most_popular_cv, "\n")
  }
  
  cat("\n")
  
} else {
  cat("No model frequency data available.\n\n")
}

cat("Note: Percentages show how often each model was selected across simulations.\n")
cat("Higher percentages indicate more consistent model selection.\n\n")



# ===============================================================================
#                              FINAL SUMMARY AND INSTRUCTIONS
# ===============================================================================

cat("\n===============================================================================\n")
cat("                              ANALYSIS COMPLETE!\n")
cat("===============================================================================\n\n")

# Display final summary
print(summary_df)

cat("\n=== FILES CREATED ===\n")
cat(" analysis_results/\n")
cat("    master_results.rds          - Main results file (for R)\n")
cat("    summary_table.csv           - Summary table \n")
cat("    complete_results.xlsx       - Excel workbook with all results\n")
cat("    results.json               - Human-readable JSON format\n")
cat("    rmst_error_comparison_aic.png - AIC comparison plot\n")
cat("    rmst_error_comparison_bic.png - BIC comparison plot\n")
cat("    improvement_percentage.png   - Improvement visualization\n")
cat("    model_selection_frequency.csv - Model selection patterns\n")
cat("    view_results.R              - Quick results viewer function\n")
cat("    detailed_[dataset].csv      - Individual dataset details\n\n")

cat("=== NEXT STEPS ===\n")
cat("results saved in multiple formats\n")
cat("visualizations \n")
cat("summary table \n")
cat("All raw data preserved for further analysis\n\n")

# Calculate and display overall improvement
all_trad_aic_errors <- sapply(all_results, function(x) x$rmst_error_trad_aic)
all_cv_aic_errors <- sapply(all_results, function(x) x$rmst_error_cv_aic)
overall_improvement <- round((mean(all_trad_aic_errors) - mean(all_cv_aic_errors)) / mean(all_trad_aic_errors) * 100, 1)

cat("KEY FINDING: Cross-validation improved prediction accuracy by", overall_improvement, "% on average\n")

cat("===============================================================================\n")

