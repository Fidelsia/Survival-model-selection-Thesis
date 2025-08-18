# IPCW Bootstrap and Model Fitting

library(survival)
library(flexsurv)
library(flexsurvcure)
library(dplyr)

# Configuration
SAVE_RESULTS <- TRUE

# Load data sets and references
master_summary <- read.csv("ipcw_extrapolation_results/master_datasets_summary_CORRECTED.csv", stringsAsFactors = FALSE)
reference_standards <- read.csv("ipcw_extrapolation_results/reference_standards.csv", stringsAsFactors = FALSE)

# Load all test data sets
test_datasets <- list()
for (i in 1:nrow(master_summary)) {
  master_id <- master_summary$master_id[i]
  filename <- paste0("ipcw_extrapolation_results/", master_id, "_master_CORRECTED.csv")
  
  if (file.exists(filename)) {
    master_data <- read.csv(filename, stringsAsFactors = FALSE)
    test_datasets[[master_id]] <- list(
      data = master_data,
      original_dataset = master_summary$original_dataset[i],
      target_censoring = master_summary$target_censoring[i],
      achieved_censoring = master_summary$achieved_censoring[i],
      type = master_summary$original_type[i]
    )
  }
}

# RMST extraction function
safe_extract_rmst <- function(model_fit, tau, model_name = "", verbose = FALSE) {
  
  if (inherits(model_fit, "try-error") || is.null(model_fit)) return(NA)
  
  tryCatch({
    
    if (inherits(model_fit, "survreg")) {
      times <- seq(0, tau, length.out = 50)
      surv_probs <- numeric(length(times))
      linear_pred <- predict(model_fit, type = "linear")[1]
      
      for (i in 1:length(times)) {
        if (times[i] == 0) {
          surv_probs[i] <- 1
        } else {
          if (model_fit$dist == "weibull") {
            z <- (log(times[i]) - linear_pred) / model_fit$scale
            surv_probs[i] <- 1 / (1 + exp(z))
          } else if (model_fit$dist == "lognormal") {
            z <- (log(times[i]) - linear_pred) / model_fit$scale
            surv_probs[i] <- 1 - pnorm(z)
          }
        }
      }
      
      rmst_val <- sum(diff(times) * (head(surv_probs, -1) + tail(surv_probs, -1)) / 2)
      
    } else if (inherits(model_fit, "flexsurvcure")) {
      rmst_result <- summary(model_fit, t = tau, type = "rmst")
      rmst_val <- rmst_result[[1]]$est
      
    } else if (inherits(model_fit, c("flexsurvspline", "flexsurvreg"))) {
      rmst_val <- tryCatch({
        rmst_result <- summary(model_fit, t = tau, type = "rmst")
        if (is.list(rmst_result) && length(rmst_result) > 0) {
          rmst_result[[1]]$est
        } else {
          rmst_result$est
        }
      }, error = function(e1) {
        tryCatch({
          times <- seq(0, tau, length.out = 100)
          surv_probs <- summary(model_fit, t = times, type = "survival")[[1]]$est
          sum(diff(times) * (head(surv_probs, -1) + tail(surv_probs, -1)) / 2)
        }, error = function(e2) NA)
      })
      
    } else {
      return(NA)
    }
    
    if (is.null(rmst_val) || any(!is.finite(rmst_val))) return(NA)
    return(rmst_val)
    
  }, error = function(e) NA)
}

# Model fitting function
fit_survival_models_detailed <- function(data, dataset_type, use_weights = FALSE, weights = NULL, verbose = FALSE) {
  
  if (use_weights && is.null(weights)) {
    return(list(models = list(), results = data.frame()))
  }
  
  if (use_weights) {
    data$ipcw_weights <- weights
  }
  
  if (dataset_type == "immunotherapy") {
    model_specs <- list(
      list(name = "weibull", type = "survreg", dist = "weibull"),
      list(name = "lognormal", type = "survreg", dist = "lognormal"),
      list(name = "cure_weibull", type = "flexsurvcure", dist = "weibull"),
      list(name = "spline", type = "flexsurvspline", k = 1, scale = "hazard")
    )
  } else {
    model_specs <- list(
      list(name = "weibull", type = "survreg", dist = "weibull"),
      list(name = "lognormal", type = "survreg", dist = "lognormal"),
      list(name = "spline", type = "flexsurvspline", k = 2, scale = "hazard")
    )
  }
  
  fitted_models <- list()
  model_results <- data.frame()
  
  for (spec in model_specs) {
    
    model_fit <- tryCatch({
      
      if (spec$type == "survreg") {
        if (use_weights) {
          survreg(Surv(time, status) ~ 1, data = data, dist = spec$dist, weights = data$ipcw_weights)
        } else {
          survreg(Surv(time, status) ~ 1, data = data, dist = spec$dist)
        }
        
      } else if (spec$type == "flexsurvcure") {
        if (use_weights) {
          flexsurvcure(Surv(time, status) ~ 1, data = data, 
                       dist = spec$dist, weights = ipcw_weights)
        } else {
          flexsurvcure(Surv(time, status) ~ 1, data = data, dist = spec$dist)
        }
        
      } else if (spec$type == "flexsurvspline") {
        fallback_configs <- list(
          list(k = spec$k, scale = spec$scale),
          list(k = 2, scale = "hazard"),
          list(k = 1, scale = "hazard"),
          list(k = 1, scale = "odds"),
          list(k = 1, scale = "normal"),
          list(k = 0, scale = "hazard")
        )
        
        spline_fit <- NULL
        
        for (config in fallback_configs) {
          spline_fit <- tryCatch({
            if (use_weights) {
              flexsurvspline(Surv(time, status) ~ 1, data = data, 
                             k = config$k, scale = config$scale, weights = weights)
            } else {
              flexsurvspline(Surv(time, status) ~ 1, data = data, 
                             k = config$k, scale = config$scale)
            }
          }, error = function(e) NULL)
          
          if (!is.null(spline_fit)) break
        }
        
        spline_fit
      }
      
    }, error = function(e) NULL)
    
    if (!is.null(model_fit)) {
      fitted_models[[spec$name]] <- model_fit
      
      cure_fraction <- NA
      if (spec$type == "flexsurvcure") {
        cure_fraction <- tryCatch({
          theta <- model_fit$coefficients[["theta"]]
          plogis(theta)
        }, error = function(e) NA)
      }
      
      model_results <- rbind(model_results, data.frame(
        model_name = spec$name,
        model_type = spec$type,
        converged = TRUE,
        aic = AIC(model_fit),
        cure_fraction = cure_fraction,
        stringsAsFactors = FALSE
      ))
      
    } else {
      model_results <- rbind(model_results, data.frame(
        model_name = spec$name,
        model_type = spec$type,
        converged = FALSE,
        aic = NA,
        cure_fraction = NA,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(list(models = fitted_models, results = model_results))
}

# IPCW weight calculation
calculate_ipcw_weights_detailed <- function(data, stabilization_percentile = 0.95, verbose = FALSE) {
  n_events <- sum(data$status)
  n_censored <- sum(1 - data$status)
  
  if (n_events < 10 || n_censored < 5) {
    return(list(weights = rep(1, nrow(data)), success = FALSE))
  }
  
  tryCatch({
    censoring_indicator <- 1 - data$status
    censoring_surv <- survfit(Surv(data$time, censoring_indicator) ~ 1)
    event_times <- data$time
    surv_summary <- summary(censoring_surv, times = event_times, extend = TRUE)
    surv_probs <- surv_summary$surv
    
    surv_probs[is.na(surv_probs)] <- 1
    surv_probs[surv_probs <= 0.001] <- 0.001
    
    weights <- rep(1, nrow(data))
    event_indices <- which(data$status == 1)
    weights[event_indices] <- 1 / surv_probs[event_indices]
    
    if (length(weights[event_indices]) > 1) {
      weight_cap <- quantile(weights[event_indices], stabilization_percentile, na.rm = TRUE)
      weights[event_indices] <- pmin(weights[event_indices], weight_cap)
    }
    
    weights[is.na(weights) | is.infinite(weights)] <- 1
    weights[weights <= 0] <- 1
    
    return(list(weights = weights, success = TRUE))
    
  }, error = function(e) {
    return(list(weights = rep(1, nrow(data)), success = FALSE))
  })
}

# Main analysis
model_test_results <- data.frame()

for (master_id in names(test_datasets)) {
  
  master_info <- test_datasets[[master_id]]
  master_data <- master_info$data
  original_dataset <- master_info$original_dataset
  
  ref_row <- reference_standards[reference_standards$dataset_name == original_dataset, ]
  if (nrow(ref_row) == 0) next
  
  reference_rmst <- ref_row$reference_rmst_90th
  time_horizon <- ref_row$time_horizon_90th
  
  cat(sprintf("--- Testing: %s ---\n", master_id))
  cat(sprintf("Original dataset: %s (%s)\n", original_dataset, master_info$type))
  cat(sprintf("Reference RMST: %.2f (time horizon: %.1f)\n", reference_rmst, time_horizon))
  cat(sprintf("Target censoring: %.1f%%, Achieved: %.1f%%\n", 
              master_info$target_censoring, master_info$achieved_censoring))
  
  if (nrow(master_data) < 500) {
    sample_size <- min(200, as.integer(nrow(master_data) * 0.6))
  } else {
    sample_size <- 300
  }
  
  for (test_iter in 1:2) {
    
    cat(sprintf("\n  Simulation %d (n=%d):\n", test_iter, sample_size))
    
    set.seed(12345 + test_iter)
    bootstrap_indices <- sample(nrow(master_data), sample_size, replace = TRUE)
    bootstrap_sample <- master_data[bootstrap_indices, ]
    
    sample_events <- sum(bootstrap_sample$status)
    sample_censoring <- round((1 - mean(bootstrap_sample$status)) * 100, 1)
    
    cat(sprintf("    Sample: n=%d, events=%d, censoring=%.1f%%\n", 
                sample_size, sample_events, sample_censoring))
    
    if (sample_events < 20) next
    
    weight_result <- calculate_ipcw_weights_detailed(bootstrap_sample, verbose = FALSE)
    if (!weight_result$success) next
    
    ipcw_weights <- weight_result$weights
    
    traditional_models <- fit_survival_models_detailed(
      bootstrap_sample, master_info$type, use_weights = FALSE, verbose = FALSE
    )
    
    ipcw_models <- fit_survival_models_detailed(
      bootstrap_sample, master_info$type, use_weights = TRUE, 
      weights = ipcw_weights, verbose = FALSE
    )
    
    for (model_name in names(traditional_models$models)) {
      if (model_name %in% names(ipcw_models$models)) {
        
        rmst_traditional <- safe_extract_rmst(
          traditional_models$models[[model_name]], time_horizon, 
          paste(model_name, "traditional"), verbose = FALSE
        )
        
        rmst_ipcw <- safe_extract_rmst(
          ipcw_models$models[[model_name]], time_horizon, 
          paste(model_name, "IPCW"), verbose = FALSE
        )
        
        if (!is.na(rmst_traditional) && !is.na(rmst_ipcw)) {
          error_traditional <- abs(rmst_traditional - reference_rmst)
          error_ipcw <- abs(rmst_ipcw - reference_rmst)
          improvement <- error_traditional - error_ipcw
          
          trad_result <- traditional_models$results[traditional_models$results$model_name == model_name, ]
          ipcw_result <- ipcw_models$results[ipcw_models$results$model_name == model_name, ]
          
          event_weights <- ipcw_weights[bootstrap_sample$status == 1]
          
          model_test_results <- rbind(model_test_results, data.frame(
            master_dataset = master_id,
            original_dataset = original_dataset,
            dataset_type = master_info$type,
            target_censoring = master_info$target_censoring,
            test_iteration = test_iter,
            model_name = model_name,
            sample_size = sample_size,
            sample_events = sample_events,
            sample_censoring = sample_censoring,
            reference_rmst = reference_rmst,
            time_horizon = time_horizon,
            rmst_traditional = rmst_traditional,
            rmst_ipcw = rmst_ipcw,
            error_traditional = error_traditional,
            error_ipcw = error_ipcw,
            improvement = improvement,
            aic_traditional = if(nrow(trad_result) > 0) trad_result$aic else NA,
            aic_ipcw = if(nrow(ipcw_result) > 0) ipcw_result$aic else NA,
            cure_fraction_traditional = if(nrow(trad_result) > 0) trad_result$cure_fraction else NA,
            cure_fraction_ipcw = if(nrow(ipcw_result) > 0) ipcw_result$cure_fraction else NA,
            weight_min = min(event_weights),
            weight_max = max(event_weights),
            weight_mean = mean(event_weights),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  cat("\n")  # Add space after each dataset
}

# Results summary
if (nrow(model_test_results) > 0) {
  
  if (SAVE_RESULTS) {
    # Save raw results as RDS
    test_results_file <- "ipcw_extrapolation_results/model_fitting_final_results.rds"
    saveRDS(model_test_results, test_results_file)
    cat(sprintf("Detailed results saved to: %s\n", test_results_file))
  }
  
  summary_results <- model_test_results %>%
    group_by(original_dataset, target_censoring, model_name) %>%
    summarise(
      n_bootstrap_runs = n(),
      mean_rmst_traditional = round(mean(rmst_traditional), 2),
      mean_rmst_ipcw = round(mean(rmst_ipcw), 2),
      mean_improvement = round(mean(improvement), 3),
      ci_lower = round(quantile(improvement, 0.025), 3),
      ci_upper = round(quantile(improvement, 0.975), 3),
      pct_beneficial = round(sum(improvement > 0) / n() * 100, 1),
      .groups = 'drop'
    ) %>%
    mutate(
      improvement_ci = paste0(mean_improvement, " (", ci_lower, ", ", ci_upper, ")")
    ) %>%
    arrange(original_dataset, target_censoring, model_name)
  
  final_summary <- summary_results[, c("original_dataset", "target_censoring", "model_name", 
                                       "mean_rmst_traditional", "mean_rmst_ipcw", "improvement_ci")]
  colnames(final_summary) <- c("Dataset", "Censoring_%", "Model", "RMST_Trad", "RMST_IPCW", "Improvement_95CI")
  
  cat("\n\nAGGREGATED RESULTS:\n")
  print(final_summary, n = Inf)
  
  if (SAVE_RESULTS) {
    # Save aggregated results as RDS
    aggregated_file <- "ipcw_extrapolation_results/aggregated_final_results.rds"
    saveRDS(summary_results, aggregated_file)
    cat(sprintf("\n\n Aggregated results saved to: %s\n", aggregated_file))
    
    # Also save the formatted summary table
    final_summary_file <- "ipcw_extrapolation_results/final_summary_table.rds"
    saveRDS(final_summary, final_summary_file)
    cat(sprintf(" Summary table saved to: %s\n", final_summary_file))
  }
  
  overall_improvement <- mean(model_test_results$improvement)
  pct_beneficial <- sum(model_test_results$improvement > 0) / nrow(model_test_results) * 100
  
  cat(sprintf("\nTOTAL MODEL COMPARISONS: %d\n", nrow(model_test_results)))
  cat(sprintf("OVERALL IPCW PERFORMANCE:\n"))
  cat(sprintf("  Average improvement: %.3f RMST units\n", overall_improvement))
  cat(sprintf("  IPCW beneficial in: %.1f%% of tests\n", pct_beneficial))
  
  model_summary <- model_test_results %>%
    group_by(model_name) %>%
    summarise(
      n_tests = n(),
      mean_improvement = round(mean(improvement), 3),
      pct_beneficial = round(sum(improvement > 0) / n() * 100, 1),
      best_improvement = round(max(improvement), 3),
      worst_improvement = round(min(improvement), 3),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_improvement))
  
  cat(sprintf("\nPERFORMANCE BY MODEL TYPE:\n"))
  print(model_summary)
  
  censoring_summary <- model_test_results %>%
    group_by(target_censoring) %>%
    summarise(
      n_tests = n(),
      mean_improvement = round(mean(improvement), 3),
      pct_beneficial = round(sum(improvement > 0) / n() * 100, 1),
      .groups = 'drop'
    ) %>%
    arrange(target_censoring)
  
  cat(sprintf("\nPERFORMANCE BY CENSORING LEVEL:\n"))
  print(censoring_summary)
  
  cat("\n✓ Analysis completed successfully\n")
}
