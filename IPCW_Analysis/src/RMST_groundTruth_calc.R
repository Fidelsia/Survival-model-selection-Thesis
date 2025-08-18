# Step 2: Reference Standard Establishment

cat("STEP 2: REFERENCE STANDARD Calculation\n")
cat("Loading cleaned datasets and calculating reference RMST values\n")
cat("Time horizon: 90th percentile of observed times\n")
cat("Method: Kaplan-Meier non-parametric estimation\n\n")

# Load cleaned datasets
load_cleaned_datasets <- function() {
  datasets <- list()
  
  if (!dir.exists("ipcw_extrapolation_results")) {
    stop("Results directory not found. Please run Step 1 first.")
  }
  
  summary_file <- "ipcw_extrapolation_results/dataset_summary.csv"
  if (!file.exists(summary_file)) {
    stop("Dataset summary not found. Please run Step 1 first.")
  }
  
  dataset_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
  cat(sprintf("Found %d datasets from Step 1:\n", nrow(dataset_summary)))
  
  for (i in 1:nrow(dataset_summary)) {
    dataset_name <- dataset_summary$dataset_name[i]
    filename <- paste0("ipcw_extrapolation_results/", dataset_name, "_cleaned.csv")
    
    if (file.exists(filename)) {
      cleaned_data <- read.csv(filename, stringsAsFactors = FALSE)
      
      if (all(c("time", "status") %in% colnames(cleaned_data))) {
        datasets[[dataset_name]] <- list(
          data = cleaned_data,
          type = dataset_summary$type[i],
          description = dataset_summary$description[i],
          original_n = dataset_summary$n_total[i],
          original_events = dataset_summary$n_events[i],
          original_censoring = dataset_summary$censoring_rate[i]
        )
        cat(sprintf("%s: n=%d, events=%d\n", 
                    dataset_name, nrow(cleaned_data), sum(cleaned_data$status)))
      } else {
        cat(sprintf("%s: Invalid data structure\n", dataset_name))
      }
    } else {
      cat(sprintf("%s: File not found\n", dataset_name))
    }
  }
  
  return(datasets)
}

datasets <- load_cleaned_datasets()

if (length(datasets) == 0) {
  stop("No datasets loaded successfully. Please check Step 1 results.")
}

# RMST calculation functions
calculate_km_rmst <- function(data, max_time = NULL, time_col = 'time', status_col = 'status') {
  
  if (is.null(max_time)) {
    max_time <- max(data[[time_col]])
  }
  
  tryCatch({
    km_fit <- survfit(Surv(data[[time_col]], data[[status_col]]) ~ 1)
    
    times <- c(0, km_fit$time[km_fit$time <= max_time], max_time)
    surv_probs <- c(1, km_fit$surv[km_fit$time <= max_time])
    
    if (length(times) != length(surv_probs)) {
      surv_probs <- c(surv_probs, tail(surv_probs, 1))
    }
    
    rmst <- sum(diff(times) * (head(surv_probs, -1) + tail(surv_probs, -1)) / 2)
    
    return(rmst)
    
  }, error = function(e) {
    cat(sprintf("Error calculating RMST: %s\n", e$message))
    return(NA)
  })
}

calculate_model_rmst <- function(model, max_time) {
  
  if (is.null(model)) return(NA)
  
  tryCatch({
    rmst_result <- summary(model, type = "rmst", t = max_time)
    if (inherits(model, "flexsurvcure")) {
      return(rmst_result[[1]]$est)
    } else {
      return(rmst_result$est)
    }
  }, error = function(e) {
    return(NA)
  })
}

# Calculate reference standards
cat("\nCalculating reference standards:\n")

reference_standards <- data.frame()

for (dataset_name in names(datasets)) {
  cat(sprintf("\nProcessing %s:\n", dataset_name))
  
  dataset_info <- datasets[[dataset_name]]
  full_data <- dataset_info$data
  
  cat(sprintf("  Dataset: %s\n", dataset_info$description))
  cat(sprintf("  Sample size: %d patients\n", nrow(full_data)))
  cat(sprintf("  Events: %d (%.1f%%)\n", sum(full_data$status), mean(full_data$status)*100))
  cat(sprintf("  Censoring: %.1f%%\n", (1-mean(full_data$status))*100))
  cat(sprintf("  Time range: %.1f - %.1f months\n", min(full_data$time), max(full_data$time)))
  
  time_horizon_90th <- quantile(full_data$time, 0.9)
  cat(sprintf("  Time horizon (90th percentile): %.1f months\n", time_horizon_90th))
  
  km_fit <- tryCatch({
    survfit(Surv(full_data$time, full_data$status) ~ 1)
  }, error = function(e) {
    cat(sprintf("  Error fitting KM curve: %s\n", e$message))
    return(NULL)
  })
  
  if (!is.null(km_fit)) {
    cat("  Kaplan-Meier curve fitted successfully\n")
    
    km_summary <- summary(km_fit)
    median_survival <- ifelse(is.na(km_summary$table["median"]), 
                              "Not reached", 
                              sprintf("%.1f months", km_summary$table["median"]))
    cat(sprintf("  Median survival: %s\n", median_survival))
    
    reference_rmst <- calculate_km_rmst(full_data, time_horizon_90th)
    
    if (!is.na(reference_rmst)) {
      cat(sprintf("  Reference RMST (t=%.1f): %.2f months\n", time_horizon_90th, reference_rmst))
      
      time_horizon_median <- median(full_data$time)
      time_horizon_85th <- quantile(full_data$time, 0.85)
      time_horizon_95th <- quantile(full_data$time, 0.95)
      
      rmst_median <- calculate_km_rmst(full_data, time_horizon_median)
      rmst_85th <- calculate_km_rmst(full_data, time_horizon_85th)
      rmst_95th <- calculate_km_rmst(full_data, time_horizon_95th)
      
      cat(sprintf("  Additional RMST values:\n"))
      cat(sprintf("    RMST (median time, %.1f): %.2f\n", time_horizon_median, rmst_median))
      cat(sprintf("    RMST (85th percentile, %.1f): %.2f\n", time_horizon_85th, rmst_85th))
      cat(sprintf("    RMST (95th percentile, %.1f): %.2f\n", time_horizon_95th, rmst_95th))
      
      reference_standards <- rbind(reference_standards, data.frame(
        dataset_name = dataset_name,
        description = dataset_info$description,
        type = dataset_info$type,
        n_patients = nrow(full_data),
        n_events = sum(full_data$status),
        censoring_rate = round((1-mean(full_data$status))*100, 1),
        min_time = min(full_data$time),
        max_time = max(full_data$time),
        median_time = median(full_data$time),
        time_horizon_90th = time_horizon_90th,
        reference_rmst_90th = reference_rmst,
        rmst_median_time = rmst_median,
        rmst_85th = rmst_85th,
        rmst_95th = rmst_95th,
        median_survival = km_summary$table["median"],
        stringsAsFactors = FALSE
      ))
      
    } else {
      cat("  Failed to calculate reference RMST\n")
    }
  }
}

# Save reference standards
cat("\nSaving reference standards:\n")

if (nrow(reference_standards) > 0) {
  reference_file <- "ipcw_extrapolation_results/reference_standards.csv"
  write.csv(reference_standards, reference_file, row.names = FALSE)
  cat(sprintf("Reference standards saved to: %s\n", reference_file))
  
  cat("\nReference Standards Summary:\n")
  display_table <- reference_standards[, c("dataset_name", "n_patients", "n_events", 
                                           "censoring_rate", "time_horizon_90th", "reference_rmst_90th")]
  colnames(display_table) <- c("Dataset", "N", "Events", "Censoring%", "Time_Horizon", "Reference_RMST")
  print(display_table)
  
  # Create KM curves visualization
  cat("\nCreating Kaplan-Meier visualization:\n")
  
  km_plot_data <- data.frame()
  
  for (dataset_name in names(datasets)) {
    if (dataset_name %in% reference_standards$dataset_name) {
      full_data <- datasets[[dataset_name]]$data
      km_fit <- survfit(Surv(full_data$time, full_data$status) ~ 1)
      
      km_times <- km_fit$time
      km_surv <- km_fit$surv
      
      plot_data <- data.frame(
        time = km_times,
        survival = km_surv,
        dataset = dataset_name
      )
      km_plot_data <- rbind(km_plot_data, plot_data)
    }
  }
  
  if (nrow(km_plot_data) > 0) {
    library(ggplot2)
    
    km_plot <- ggplot(km_plot_data, aes(x = time, y = survival, color = dataset)) +
      geom_step(size = 1) +
      labs(
        title = "Kaplan-Meier Survival Curves - Reference Standards",
        subtitle = "Vertical lines show 90th percentile time horizons",
        x = "Time (months)",
        y = "Survival Probability",
        color = "Dataset"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom") +
      ylim(0, 1)
    
    for (i in 1:nrow(reference_standards)) {
      km_plot <- km_plot + 
        geom_vline(xintercept = reference_standards$time_horizon_90th[i], 
                   linetype = "dashed", alpha = 0.7)
    }
    
    plot_file <- "ipcw_extrapolation_results/reference_km_curves.png"
    ggsave(plot_file, km_plot, width = 12, height = 8, dpi = 300)
    cat(sprintf("KM curves plot saved to: %s\n", plot_file))
    
    print(km_plot)
  }
  
} else {
  stop("No reference standards calculated successfully.")
}

cat("\nStep 2 Completed: Reference Standard Calculation\n")
cat("Files created:\n")
cat("• reference_standards.csv\n")
cat("• reference_km_curves.png\n\n")

cat("Reference Standards Summary:\n")
for (i in 1:nrow(reference_standards)) {
  row <- reference_standards[i, ]
  cat(sprintf("• %s: RMST = %.2f months (t = %.1f)\n", 
              row$dataset_name, row$reference_rmst_90th, row$time_horizon_90th))
}

cat("\nReady for Step 3: Artificial Censoring Simulation\n")