# Step 3: Master Artificially Censored Data sets

cat("STEP 3: MASTER ARTIFICIALLY CENSORED DATASETS\n")
cat("Creating master datasets with controlled high censoring levels\n")
cat("Approach: Full dataset censoring before sampling\n")
cat("Target censoring levels: 40%, 60%, 80%\n\n")

# Load previous results
cat("Loading cleaned datasets and reference standards:\n")

if (!dir.exists("ipcw_extrapolation_results")) {
  stop("Results directory not found. Please run Steps 1 and 2 first.")
}

dataset_summary <- read.csv("ipcw_extrapolation_results/dataset_summary.csv", stringsAsFactors = FALSE)
reference_standards <- read.csv("ipcw_extrapolation_results/reference_standards.csv", stringsAsFactors = FALSE)

datasets <- list()
for (i in 1:nrow(dataset_summary)) {
  dataset_name <- dataset_summary$dataset_name[i]
  filename <- paste0("ipcw_extrapolation_results/", dataset_name, "_cleaned.csv")
  
  if (file.exists(filename)) {
    cleaned_data <- read.csv(filename, stringsAsFactors = FALSE)
    datasets[[dataset_name]] <- list(
      data = cleaned_data,
      type = dataset_summary$type[i],
      description = dataset_summary$description[i]
    )
  }
}

cat(sprintf("Loaded %d datasets from previous steps\n", length(datasets)))
cat(sprintf("Loaded reference standards for %d datasets\n", nrow(reference_standards)))

# Artificial censoring function
apply_artificial_censoring_master_fixed <- function(df, target_censoring_rate, dataset_name, time_col = 'time', status_col = 'status') {
  
  cat(sprintf("    Target censoring rate: %.0f%%\n", target_censoring_rate * 100))
  cat(sprintf("    Original dataset size: n = %d\n", nrow(df)))
  
  current_censoring_rate <- 1 - mean(df[[status_col]])
  cat(sprintf("    Current censoring: %.1f%%\n", current_censoring_rate * 100))
  
  if (current_censoring_rate >= target_censoring_rate) {
    cat(sprintf("    Target already achieved (current %.1f%% >= target %.1f%%)\n", 
                current_censoring_rate * 100, target_censoring_rate * 100))
    return(list(
      data = df,
      cutoff_time = max(df[[time_col]]),
      target_censoring = target_censoring_rate * 100,
      achieved_censoring = current_censoring_rate * 100,
      original_n = nrow(df),
      final_n = nrow(df),
      final_events = sum(df[[status_col]])
    ))
  }
  
  additional_censoring_needed = target_censoring_rate - current_censoring_rate
  cat(sprintf("    Additional censoring needed: %.1f%%\n", additional_censoring_needed * 100))
  
  current_event_rate <- mean(df[[status_col]])
  fraction_of_events_to_censor <- additional_censoring_needed / current_event_rate
  
  cat(sprintf("    Need to censor %.1f%% of current events\n", fraction_of_events_to_censor * 100))
  
  event_times <- df[[time_col]][df[[status_col]] == 1]
  
  if (fraction_of_events_to_censor >= 1) {
    cutoff_time <- min(event_times) - 0.01
    cat(sprintf("    Censoring ALL events (cutoff: %.2f)\n", cutoff_time))
  } else {
    cutoff_quantile <- 1 - fraction_of_events_to_censor
    cutoff_time <- quantile(event_times, cutoff_quantile, type = 1)
    cat(sprintf("    Cutoff time (%.1f%% quantile of events): %.1f months\n", 
                cutoff_quantile * 100, cutoff_time))
  }
  
  df_censored <- df
  
  censor_mask <- (df_censored[[status_col]] == 1) & (df_censored[[time_col]] > cutoff_time)
  
  df_censored[[status_col]][censor_mask] <- 0
  df_censored[[time_col]][censor_mask] <- cutoff_time
  
  achieved_censoring_rate <- 1 - mean(df_censored[[status_col]])
  achieved_censoring_pct <- achieved_censoring_rate * 100
  
  cat(sprintf("    Achieved censoring: %.1f%%\n", achieved_censoring_pct))
  cat(sprintf("    Final events: %d\n", sum(df_censored[[status_col]])))
  
  censoring_error <- abs(achieved_censoring_pct - (target_censoring_rate * 100))
  cat(sprintf("    Censoring accuracy: ±%.1f%% from target\n", censoring_error))
  
  return(list(
    data = df_censored,
    cutoff_time = cutoff_time,
    target_censoring = target_censoring_rate * 100,
    achieved_censoring = achieved_censoring_pct,
    censoring_error = censoring_error,
    original_n = nrow(df),
    final_n = nrow(df_censored),
    final_events = sum(df_censored[[status_col]]),
    events_censored = sum(censor_mask)
  ))
}

# Create master datasets with corrected logic
cat("\nCreating master datasets with corrected logic:\n")

if (!exists("datasets") || length(datasets) == 0) {
  dataset_summary <- read.csv("ipcw_extrapolation_results/dataset_summary.csv", stringsAsFactors = FALSE)
  datasets <- list()
  for (i in 1:nrow(dataset_summary)) {
    dataset_name <- dataset_summary$dataset_name[i]
    filename <- paste0("ipcw_extrapolation_results/", dataset_name, "_cleaned.csv")
    if (file.exists(filename)) {
      cleaned_data <- read.csv(filename, stringsAsFactors = FALSE)
      datasets[[dataset_name]] <- list(
        data = cleaned_data,
        type = dataset_summary$type[i],
        description = dataset_summary$description[i]
      )
    }
  }
}

target_censoring_levels <- c(0.4, 0.6, 0.8)
master_datasets_corrected <- list()
master_summary_corrected <- data.frame()

for (dataset_name in names(datasets)) {
  cat(sprintf("\nProcessing %s:\n", dataset_name))
  
  original_dataset <- datasets[[dataset_name]]
  original_data <- original_dataset$data
  
  cat(sprintf("Dataset: %s\n", original_dataset$description))
  
  for (target_rate in target_censoring_levels) {
    
    cat(sprintf("  Creating master dataset with %.0f%% target censoring:\n", target_rate * 100))
    
    master_result <- apply_artificial_censoring_master_fixed(
      original_data, 
      target_rate, 
      dataset_name
    )
    
    if (!is.null(master_result)) {
      master_id <- sprintf("%s_Censored_%dpct", dataset_name, target_rate * 100)
      
      master_datasets_corrected[[master_id]] <- list(
        data = master_result$data,
        original_dataset = dataset_name,
        target_censoring = target_rate * 100,
        achieved_censoring = master_result$achieved_censoring,
        censoring_error = master_result$censoring_error,
        cutoff_time = master_result$cutoff_time,
        events_censored = master_result$events_censored,
        type = original_dataset$type,
        description = sprintf("%s (%.0f%% artificial censoring)", 
                              original_dataset$description, target_rate * 100)
      )
      
      master_summary_corrected <- rbind(master_summary_corrected, data.frame(
        master_id = master_id,
        original_dataset = dataset_name,
        original_type = original_dataset$type,
        target_censoring = target_rate * 100,
        achieved_censoring = round(master_result$achieved_censoring, 1),
        censoring_error = round(master_result$censoring_error, 1),
        cutoff_time = round(master_result$cutoff_time, 1),
        original_n = master_result$original_n,
        final_n = master_result$final_n,
        final_events = master_result$final_events,
        events_censored = master_result$events_censored,
        stringsAsFactors = FALSE
      ))
      
      cat(sprintf("    Master dataset created: %s\n", master_id))
      
    } else {
      cat(sprintf("    Failed to create master dataset for %.0f%% censoring\n", target_rate * 100))
    }
  }
}

# Save corrected master datasets
cat("\nSaving corrected master datasets:\n")

if (length(master_datasets_corrected) > 0) {
  
  for (master_id in names(master_datasets_corrected)) {
    filename <- paste0("ipcw_extrapolation_results/", master_id, "_master_CORRECTED.csv")
    write.csv(master_datasets_corrected[[master_id]]$data, filename, row.names = FALSE)
    cat(sprintf("%s saved\n", master_id))
  }
  
  summary_file <- "ipcw_extrapolation_results/master_datasets_summary_CORRECTED.csv"
  write.csv(master_summary_corrected, summary_file, row.names = FALSE)
  cat(sprintf("Summary saved to: %s\n", summary_file))
  
  cat("\nMaster Datasets Summary:\n")
  display_summary_corrected <- master_summary_corrected[, c("master_id", "target_censoring", "achieved_censoring", 
                                                            "censoring_error", "final_n", "final_events")]
  colnames(display_summary_corrected) <- c("Master_Dataset", "Target_%", "Achieved_%", "Error_%", "N", "Events")
  print(display_summary_corrected)
  
  cat("\nCensoring Accuracy Assessment:\n")
  accuracy_summary <- master_summary_corrected %>%
    group_by(target_censoring) %>%
    summarise(
      n_datasets = n(),
      mean_achieved = round(mean(achieved_censoring), 1),
      mean_error = round(mean(censoring_error), 1),
      max_error = round(max(censoring_error), 1),
      .groups = 'drop'
    )
  
  cat("Accuracy by target censoring level:\n")
  print(accuracy_summary)
  
  acceptable_accuracy <- all(master_summary_corrected$censoring_error <= 5)
  if (acceptable_accuracy) {
    cat("\nAll censoring targets achieved within ±5% accuracy\n")
  } else {
    high_error_datasets <- master_summary_corrected$master_id[master_summary_corrected$censoring_error > 5]
    cat(sprintf("\nHigh censoring error (>5%%) in: %s\n", paste(high_error_datasets, collapse = ", ")))
  }
  
} else {
  stop("No corrected master datasets created successfully.")
}

cat("\nStep 3 Completed: Artificial Censoring\n")
cat("Ready for simulation sampling\n")