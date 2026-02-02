### All functions needed for step 2C

# Function to stop all clusters
stopAllClusters <- function() {
  if (exists("cl") && inherits(cl, "cluster") && !is.null(cl)) {
    cat("Stopping existing cluster...\n")
    tryCatch({
      stopCluster(cl)
      registerDoSEQ() # Reset to sequential processing
    }, error = function(e) {
      cat("Error stopping cluster: ", e$message, "\n")
    })
  } else {
    cat("No active clusters to stop.\n")
  }
}

# Function to ensure all predictors are numeric
ensure_numeric <- function(df, predictors) {
  for (col in predictors) {
    df[[col]] <- as.numeric(as.character(df[[col]]))
  }
  return(df)
}

# Custom MAPE metric function
mape_summary <- function(data, lev = NULL, model = NULL) {
  mape <- mean(abs((data$obs - data$pred) / data$obs), na.rm = TRUE)
  c(MAPE = mape)
}

# Custom summary function for MSE
mse_summary <- function(data, lev = NULL, model = NULL) {
  out <- mean((data$obs - data$pred)^2)
  names(out) <- "MSE"
  out
}

# Function process_data: for data of low intensity, its highest Mn signal's true intensity is predicted 
# through the previously trained model, so the predicted intensity could be considered as the true intensity to
# predict the true intenstiy of other Mns based on the theoretical ratios
process_data <- function(Matrix_data, info_data, clf_rf) {
  # Initialize lists
  Matrix_m0 <- list()
  Matrix_m1 <- list()
  Matrix_m2 <- list()
  info_m0 <- list()
  info_m1 <- list()
  info_m2 <- list()
  
  # Create a row with all NA values for both data.tables
  na_row_dt1 <- as.list(setNames(rep(NA, ncol(Matrix_data)), names(Matrix_data)))
  na_row_dt2 <- as.list(setNames(rep(NA, ncol(info_data)), names(info_data)))
  
  # Combined function to process rows
  process_combined_rows <- function(dt1, dt2, m0_1, m1_1, m2_1, m0_2, m1_2, m2_2, zero_row1, zero_row2) {
    for (i in 1:nrow(dt1)) {
      if (dt1[i, mn] == 0 && i == 1) {
        m0_1 <- append(m0_1, list(dt1[i, ]))
        m0_2 <- append(m0_2, list(dt2[i, ]))
      } else if (dt1[i, mn] == 0 && i > 1 && dt1[i-1, mn] != 1) {
        m0_1 <- append(m0_1, list(dt1[i, ]))
        m0_2 <- append(m0_2, list(dt2[i, ]))
      } else if (i == nrow(dt1) && dt1[i, mn] == 1) {
        m1_1 <- append(m1_1, list(dt1[i, ]))
        m1_2 <- append(m1_2, list(dt2[i, ]))
        m2_1 <- append(m2_1, list(zero_row1))
        m2_2 <- append(m2_2, list(zero_row2))
      } else if (i > 1 && dt1[i-1, mn] == 0 && dt1[i, mn] == 1) {
        m1_1 <- append(m1_1, list(dt1[i, ]))
        m1_2 <- append(m1_2, list(dt2[i, ]))
      } else if (i > 1 && dt1[i-1, mn] == 1) {
        if (dt1[i, mn] == 2) {
          m2_1 <- append(m2_1, list(dt1[i, ]))
          m2_2 <- append(m2_2, list(dt2[i, ]))
        } else if (dt1[i, mn] == 0) {
          # Insert row with all 0 values
          m2_1 <- append(m2_1, list(zero_row1))
          m2_2 <- append(m2_2, list(zero_row2))
          # Copy row to m0
          m0_1 <- append(m0_1, list(dt1[i, ]))
          m0_2 <- append(m0_2, list(dt2[i, ]))
        }
      }
    }
    return(list(
      m0_1 = rbindlist(m0_1, fill = TRUE), 
      m1_1 = rbindlist(m1_1, fill = TRUE), 
      m2_1 = rbindlist(m2_1, fill = TRUE), 
      m0_2 = rbindlist(m0_2, fill = TRUE), 
      m1_2 = rbindlist(m1_2, fill = TRUE), 
      m2_2 = rbindlist(m2_2, fill = TRUE)
    ))
  }
  
  # Convert the initial data.tables to as.data.table if not already
  Matrix_data <- as.data.table(Matrix_data)
  info_data <- as.data.table(info_data)
  
  # Process the combined data
  result_combined <- process_combined_rows(Matrix_data, info_data, 
                                           Matrix_m0, Matrix_m1, Matrix_m2,
                                           info_m0, info_m1, info_m2,
                                           na_row_dt1, na_row_dt2)
  
  Matrix_m0 <- as.matrix(result_combined$m0_1)
  Matrix_m1 <- as.matrix(result_combined$m1_1)
  Matrix_m2 <- as.matrix(result_combined$m2_1)
  info_m0 <- result_combined$m0_2
  info_m1 <- result_combined$m1_2
  info_m2 <- result_combined$m2_2
  
  n <- nrow(info_m0)
  
  # Initialize a vector to store the results
  highest_ratio_source <- integer(n)
  
  # Iterate through each row from i to n
  for (i in 1:n) {
    # Extract the "ratio" values from the current row, treating NA as -Inf
    ratio_m0 <- ifelse(is.na(info_m0[i, "ratio"]), -Inf, info_m0[i, "ratio"])
    ratio_m1 <- ifelse(is.na(info_m1[i, "ratio"]), -Inf, info_m1[i, "ratio"])
    ratio_m2 <- ifelse(is.na(info_m2[i, "ratio"]), -Inf, info_m2[i, "ratio"])
    
    # Determine which ratio is the highest with priority handling
    if (as.numeric(ratio_m0) >= as.numeric(ratio_m1) & as.numeric(ratio_m0) >= as.numeric(ratio_m2)) {
      highest_ratio_source[i] <- 0
    } else if (as.numeric(ratio_m1) >= as.numeric(ratio_m0) & as.numeric(ratio_m1) >= as.numeric(ratio_m2)) {
      highest_ratio_source[i] <- 1
    } else if (as.numeric(ratio_m2) > as.numeric(ratio_m0) & as.numeric(ratio_m2) > as.numeric(ratio_m1)) {
      highest_ratio_source[i] <- 2
    }
  }
  
  # Set -Inf values back to NA in the original data frames
  info_m0[info_m0 == -Inf] <- NA
  info_m1[info_m1 == -Inf] <- NA
  info_m2[info_m2 == -Inf] <- NA
  
  # Optionally, convert to a data frame for better readability
  result_df <- data.frame(row_index = 1:n, highest_ratio_source)
  
  # Identify the last row with non-NA data in Matrix_m2
  last_non_na_row_m2 <- max(which(rowSums(is.na(Matrix_m2)) == 0))
  
  # Precompute predictions for each model
  predictions_m0_full <- info_m0$intensity_minus_bg_smoothed / predict(clf_rf, newdata = Matrix_m0)
  predictions_m1_full <- info_m1$intensity_minus_bg_smoothed / predict(clf_rf, newdata = Matrix_m1)
  
  # For Matrix_m2, predict only up to the last row with non-NA data
  predictions_m2_full <- rep(NA, n)
  predictions_m2_full[1:last_non_na_row_m2] <- info_m2$intensity_minus_bg_smoothed[1:last_non_na_row_m2] / predict(clf_rf, newdata = Matrix_m2[1:last_non_na_row_m2, , drop = FALSE])
  
  # Initialize predictions vectors
  predictions_m0 <- numeric(n)
  predictions_m1 <- numeric(n)
  predictions_m2 <- numeric(n)
  
  # Iterate through each row from i to n
  for (i in 1:n) {
    if (result_df$highest_ratio_source[i] == 0) {
      predictions_m0[i] <- predictions_m0_full[i]
      predictions_m1[i] <- predictions_m0[i] * info_m1$ratio[i] / info_m0$ratio[i]
      predictions_m2[i] <- predictions_m0[i] * info_m2$ratio[i] / info_m0$ratio[i]
    } else if (result_df$highest_ratio_source[i] == 1) {
      predictions_m1[i] <- predictions_m1_full[i]
      predictions_m0[i] <- predictions_m1[i] * info_m0$ratio[i] / info_m1$ratio[i]
      predictions_m2[i] <- predictions_m1[i] * info_m2$ratio[i] / info_m1$ratio[i]
    } else if (result_df$highest_ratio_source[i] == 2) {
      predictions_m2[i] <- predictions_m2_full[i]
      predictions_m0[i] <- predictions_m2[i] * info_m0$ratio[i] / info_m2$ratio[i]
      predictions_m1[i] <- predictions_m2[i] * info_m1$ratio[i] / info_m2$ratio[i]
    }
  }
  
  info_m0$intensity_smoothed_predicted <- predictions_m0
  info_m1$intensity_smoothed_predicted <- predictions_m1
  info_m2$intensity_smoothed_predicted <- predictions_m2
  
  # Initialize an empty data frame for info_pseudoL
  info_pseudoL <- data.frame()
  
  # Loop through each row
  for (i in 1:nrow(info_m0)) {
    # Copy row i from info_m0
    info_pseudoL <- rbind(info_pseudoL, info_m0[i, ])
    
    # Copy row i from info_m1
    info_pseudoL <- rbind(info_pseudoL, info_m1[i, ])
    
    # Check if row i from info_m2 is all NA
    if (!all(is.na(info_m2[i, ]))) {
      # Copy row i from info_m2 if not all NA
      info_pseudoL <- rbind(info_pseudoL, info_m2[i, ])
    }
  }
  
  list(
    Matrix_m0 = Matrix_m0,
    Matrix_m1 = Matrix_m1,
    Matrix_m2 = Matrix_m2,
    info_m0 = info_m0,
    info_m1 = info_m1,
    info_m2 = info_m2,
    predictions_m0 = predictions_m0,
    predictions_m1 = predictions_m1,
    predictions_m2 = predictions_m2,
    info_pseudoL = info_pseudoL
  )
}

# Function to validate predictions using the ratio
validate_predictions <- function(predictions, actuals, ratios) {
  groups <- split(actuals, interaction(actuals$compound, actuals$sample, actuals$slice))
  total_diff <- 0
  total_count <- 0
  
  for (group in groups) {
    predicted_group <- predictions[rownames(group)]
    actual_ratios <- group$ratio
    sum_predicted <- sum(predicted_group)
    
    if (sum_predicted > 0) {
      predicted_ratios <- predicted_group / sum_predicted
      total_diff <- total_diff + sum(abs(predicted_ratios - actual_ratios))
      total_count <- total_count + length(predicted_group)
    }
  }
  
  mean_ratio_diff <- total_diff / total_count
  return(mean_ratio_diff)
}

# Function needed to do actual corrects
process_and_export_corrected_df <- function(df_by_mn_input, df_by_slice_input, df_input, clf_model, seed_i, output_prefix) {
  set.seed(seed_i)
  
  df_by_mn_tc <- as.data.table(df_by_mn_input)
  df_by_mn_tc_info <- df_by_mn_tc[, .(compound, sample, formula, slice,
                                      intensity_predicted, intensity_smoothed_predicted,
                                      ratio, intensity_minus_bg, intensity_minus_bg_smoothed)]

  # Clean matrix to get same structure as model required
  df_by_mn_tc[, c("compound", "sample", "formula", "ratio",
                  "intensity_predicted", "intensity_smoothed_predicted",
                  "intensity_minus_bg", "intensity_minus_bg_smoothed",
                  "rt_left", "rt_peak", "rt_right", "slice",
                  "total_slices", "slice_peak_height", "mplus1_intensity_smoothed",
                  "highest_original_intensity", "rt_width",
                  "slice_diff_to_m0peak", "Si") := NULL]
  
  df_by_mn_tc <- df_by_mn_tc[, lapply(.SD, as.numeric)]
  
  # Make prediction
  Matrix_tc <- as.matrix(df_by_mn_tc)
  predictions <- df_by_mn_tc_info$intensity_minus_bg_smoothed / predict(clf_model, newdata = Matrix_tc)
  
  # Append back compound/sample/slice
  df_by_mn_tc$intensity_corrected <- predictions
  df_by_mn_tc$compound <- df_by_mn_tc_info$compound
  df_by_mn_tc$sample <- df_by_mn_tc_info$sample
  df_by_mn_tc$slice <- df_by_mn_tc_info$slice
  
  # Merge with by slice data
  df_by_slice_tc <- df_by_slice_input %>%
    mutate(across(c(compound, sample, slice), as.character)) %>%
    mutate(m0_intensity_corrected = NA_real_,
           m1_intensity_corrected = NA_real_,
           m2_intensity_corrected = NA_real_)
  
  df_by_mn_tc <- df_by_mn_tc %>%
    mutate(across(c(compound, sample, slice), as.character))
  
  merged_df <- df_by_slice_tc %>%
    left_join(df_by_mn_tc, by = c("compound", "sample", "slice"))
  
  merged_df <- merged_df %>%
    mutate(
      m0_intensity_corrected = ifelse(mn == 0, intensity_corrected, m0_intensity_corrected),
      m1_intensity_corrected = ifelse(mn == 1, intensity_corrected, m1_intensity_corrected),
      m2_intensity_corrected = ifelse(mn == 2, intensity_corrected, m2_intensity_corrected)
    ) %>%
    select(-intensity_corrected, -mn)
  
  # Summarize corrected intensities by slice
  result_df <- merged_df %>%
    group_by(compound, sample, slice) %>%
    summarize(
      m0_intensity_corrected = ifelse(all(is.na(m0_intensity_corrected)), NA_real_, max(m0_intensity_corrected, na.rm = TRUE)),
      m1_intensity_corrected = ifelse(all(is.na(m1_intensity_corrected)), NA_real_, max(m1_intensity_corrected, na.rm = TRUE)),
      m2_intensity_corrected = ifelse(all(is.na(m2_intensity_corrected)), NA_real_, max(m2_intensity_corrected, na.rm = TRUE)),
      .groups = 'drop'
    )
  
  df_by_slice_corrected <- df_by_slice_tc %>%
    left_join(result_df, by = c("compound", "sample", "slice")) %>%
    mutate(
      m0_intensity_corrected = coalesce(m0_intensity_corrected.x, m0_intensity_corrected.y),
      m1_intensity_corrected = coalesce(m1_intensity_corrected.x, m1_intensity_corrected.y),
      m2_intensity_corrected = coalesce(m2_intensity_corrected.x, m2_intensity_corrected.y)
    ) %>%
    select(-m0_intensity_corrected.x, -m1_intensity_corrected.x, -m2_intensity_corrected.x,
           -m0_intensity_corrected.y, -m1_intensity_corrected.y, -m2_intensity_corrected.y)
  
  # Merge back to by metabolite/by sample 
  summarized_df <- df_by_slice_corrected %>%
    group_by(compound, sample) %>%
    summarize(
      m0_area_corrected = sum(m0_intensity_corrected, na.rm = TRUE),
      m1_area_corrected = sum(m1_intensity_corrected, na.rm = TRUE),
      m2_area_corrected = sum(m2_intensity_corrected, na.rm = TRUE),
      .groups = 'drop'
    )
  
  df_corrected <- df_input %>%
    left_join(summarized_df, by = c("compound", "sample")) %>%
    mutate(
      m0_observed_corrected = ifelse(integrated_to_area == 0,
                                     m0_area_corrected / (m0_area_corrected + m1_area_corrected),
                                     ifelse(integrated_to_area == 1,
                                            m0_area_corrected / (m0_area_corrected + m1_area_corrected + m2_area_corrected),
                                            NA_real_)),
      m1_observed_corrected = ifelse(integrated_to_area == 0,
                                     m1_area_corrected / (m0_area_corrected + m1_area_corrected),
                                     ifelse(integrated_to_area == 1,
                                            m1_area_corrected / (m0_area_corrected + m1_area_corrected + m2_area_corrected),
                                            NA_real_)),
      m2_observed_corrected = ifelse(integrated_to_area == 1,
                                     m2_area_corrected / (m0_area_corrected + m1_area_corrected + m2_area_corrected),
                                     NA_real_),
      delta_m0_corrected = m0_observed_corrected - m0_predicted,
      delta_m1_corrected = m1_observed_corrected - m1_predicted,
      delta_m2_corrected = ifelse(integrated_to_area == 1, m2_observed_corrected - m2_predicted, NA_real_),
      delta_m0_percentage_corrected = delta_m0_corrected / m0_predicted,
      delta_m1_percentage_corrected = delta_m1_corrected / m1_predicted,
      delta_m2_percentage_corrected = ifelse(integrated_to_area == 1, delta_m2_corrected / m2_predicted, NA_real_)
    )
  
  # Save both by sample and by slice data sheets
  output_list <- list(
    Corrected_Compound_Level = df_corrected,
    Corrected_Slice_Level = df_by_slice_corrected
  )
  
  write_xlsx(output_list, path = paste0(output_prefix, "_seed", seed_i, "_3Sig_sup1.xlsx"))
}