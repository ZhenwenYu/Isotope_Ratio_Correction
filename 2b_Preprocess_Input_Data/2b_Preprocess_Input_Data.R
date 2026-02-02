###################################

# File: 2b_Preprocess_Input_Data
# Input: "All_1c_Combined_2a.xlsx" is the output of step 2a, which is all combined data of peak and by slice infor extracted in step 1.
# Function: Filter the data to remove any extreme outliers (which is likely due to misidentifidation);
#         Calculaten some more features that will be used in step 2c model training;
#         Inspect the input data structure to determine an intensity cutoff to know which signals could be considered
# Output: Two Rdata files: "df_by_mn_high_2b.RData" and "df_by_mn_low_2b.RData" , which contains data above and below the determined intensity cutoff

# Usage: Before running, copy the 2a.xlsx file to the same folder as this 2b codes. During running, keep an eye on the plot and final summary
#       to check if the results make sense
# Next to Do: Copy the four Rdata "df_by_mn_high_2b.RData" and "df_by_mn_low_2b.RData" 
#           "df_2b.RData" and "df_by_slice_2b.RData" to 2cRF and 2cEN folders' "input"  folder to do actual model training
#
# Disclaimer: Chatgpt is used to debug this code.

####################################
#####################################

library(readxl) #load xlsx
library(purrr) #load xlsx by iterations if input files are not combined through step 2a
library(openxlsx) #write xlsx
library(magrittr) #pipe
library(tidyr) #df cleaning
library(plyr)  #df splitting and combining
library(dplyr) #df data manipulation
library(caret) #training
library(randomForest) #rf related plot
library(doParallel) # parallel computation
library(data.table) # handle large data without using too much RAM
library(ggplot2) # for plotting

# Indicate the outlier_cutoff of data points. By default does not need to be modified by user, 3 means the data point's percentage error of
# MID (mass isotopomer distribution) should be within the +-3 sd of percentage error of all data points
outlier_cutoff <- 3
# Indicate the percentage_error_cutoff. By default, it does not need to be modified by the user. A value of 0.05 means that if any two of a single slice's M0, M1, and M2
# are within ±5% error of the true intensity, the code will find the minimum intensity such that the slice's data errors remain within
# the cutoff range and do not need to be corrected. In our data, we find that if any two of a slice's M0, M1, or M2 signal intensities are ≥ 1e7,
# that slice's MID will be accurate enough and will not need correction. Therefore, we can use this ground-truth Mn intensity to calculate
# the true intensity of the remaining Mn signals, as the predicted ratio of M0, M1, and M2 is known.
percentage_error_cutoff <- 0.05
# epsilon is a very small number to be added to 0 in case when the noise or background = 0 makes the 
# ratio of Signal/Noise or Signal/Baseline to cause mathematical error. By default, it should not be modifed by the user.
epsilon <- 1e-9

####################################
#####################################

# Get and check the current working directory
directory <- getwd()
cat("Current working directory:\n")
print(directory)
# List all 1c XLSX files in the directory
files <- list.files(directory, pattern = "\\_2a.xlsx$", full.names = TRUE)
cat("2a xlsx result files found:\n")
print(files)
# Check if files were found
if (length(files) == 0) {
  stop("No 2a xlsx result files found in the current directory.")
}
# Find all 2a files in case multiple 2a files are inputted at a time
file_paths <- list.files(path = directory, pattern = "\\_2a.xlsx$", full.names = TRUE, recursive = TRUE)

# Function read_first_sheet and read_second_sheet to read sheets seperately
read_first_sheet <- function(file_path) {
  suppressWarnings(read_excel(file_path, sheet = 1))
}
read_second_sheet <- function(file_path) {
  suppressWarnings(read_excel(file_path, sheet = 2))
}
# loop for each 2a files
data_sheet_1 <- map(file_paths, read_first_sheet)
data_sheet_2 <- map(file_paths, read_second_sheet)
# Store all sheet1 integrated area data and sheet2 by slice data seperately
df_raw <- janitor::clean_names(bind_rows(data_sheet_1))
df_by_slice_raw <- janitor::clean_names(bind_rows(data_sheet_2))

#######################################################################

# Identify rows with any area equal to 0
rows_with_null_area <- with(df_raw, is.na(m0_area) & is.na(m1_area) & is.na(m2_area))

# Identify rows with super outliers by first
# calculate the mean and standard deviation for delta_m1_percentage and delta_m2_percentage
# Combine the columns into a single vector, removing NA values
combined_values <- c(df_raw$delta_m0_percentage_sub_b, df_raw$delta_m1_percentage_sub_b, df_raw$delta_m2_percentage_sub_b)
combined_values <- combined_values[!is.na(combined_values)]
mean_combined <- mean(combined_values)
sd_combined <- sd(combined_values)
upper_bound_combined <- mean_combined + outlier_cutoff * sd_combined
lower_bound_combined <- mean_combined - outlier_cutoff * sd_combined
rows_with_super_outliers <- with(df_raw, 
                                 (!is.na(delta_m0_percentage_sub_b) & 
                                  (delta_m0_percentage_sub_b > upper_bound_combined | delta_m0_percentage_sub_b < lower_bound_combined)) |
                                 (!is.na(delta_m1_percentage_sub_b) & 
                                  (delta_m1_percentage_sub_b > upper_bound_combined | delta_m1_percentage_sub_b < lower_bound_combined)) |
                                 (!is.na(delta_m2_percentage_sub_b) & 
                                  (delta_m2_percentage_sub_b > upper_bound_combined | delta_m2_percentage_sub_b < lower_bound_combined)))


# Combined info of rows with any area equal to 0 and outlier for row removal
rows_to_remove <- rows_with_null_area | rows_with_super_outliers
# Subset the dataframe to remove specified rows from the sheet 1 integrated area df and keep them in a separate dataframe
df_removed <- df_raw[rows_to_remove, ]
df <- df_raw[!rows_to_remove, ]

# Reorder sheet 1 df by the first column after filtering and realign the sheet 2 by slice data with
df <- df[order(df[[1]]), ]
df_M0M1M2 <- subset(df, integrated_to_area == 1)
df_M0M1 <- subset(df, integrated_to_area_only_m0m1 == 1 & integrated_to_area != 1)
temp_df_by_slice <- semi_join(df_by_slice_raw, df, by = c("compound", "sample"))
temp_df_by_slice_M0M1M2 <- semi_join(df_by_slice_raw, df_M0M1M2, by = c("compound", "sample"))
temp_df_by_slice_M0M1 <- semi_join(df_by_slice_raw, df_M0M1, by = c("compound", "sample"))
df_by_slice_M0M1M2 <- subset(temp_df_by_slice_M0M1M2, integrated_to_area == 1)
df_by_slice_M0M1 <- subset(temp_df_by_slice_M0M1, integrated_to_area_only_m0m1 == 1 & integrated_to_area != 1)

# Identify sheet 2 rows (that have all M0M1M2 data) with outliers
M0M1M2_rows_with_super_outliers <- with(df_by_slice_M0M1M2, 
                                 (!is.na(delta_m0_percentage_sub_b) & 
                                  (delta_m0_percentage_sub_b > upper_bound_combined | delta_m0_percentage_sub_b < lower_bound_combined)) |
                                 (!is.na(delta_m1_percentage_sub_b) & 
                                  (delta_m1_percentage_sub_b > upper_bound_combined | delta_m1_percentage_sub_b < lower_bound_combined)) |
                                 (!is.na(delta_m2_percentage_sub_b) & 
                                  (delta_m2_percentage_sub_b > upper_bound_combined | delta_m2_percentage_sub_b < lower_bound_combined)))
# Subset the sheet 2 df to remove outlier rows and keep them in a separate df
df_by_slice_M0M1M2_removed <- df_by_slice_M0M1M2[M0M1M2_rows_with_super_outliers, ]
df_by_slice_M0M1M2 <- df_by_slice_M0M1M2[!M0M1M2_rows_with_super_outliers, ]
# Repeat to identify rows with outliers in slices where only M0 and M1 are recoreded
M0M1_rows_with_super_outliers <- with(df_by_slice_M0M1,
                                (!is.na(delta_m0_percentage_sub_b) &
                                  (delta_m0_percentage_sub_b > upper_bound_combined | delta_m0_percentage_sub_b < lower_bound_combined)) |
                                (!is.na(delta_m1_percentage_sub_b) &
                                  (delta_m1_percentage_sub_b > upper_bound_combined | delta_m1_percentage_sub_b < lower_bound_combined)))
# Subset the dataframe to remove specified rows and keep them in a separate dataframe
df_by_slice_M0M1_removed <- df_by_slice_M0M1[M0M1_rows_with_super_outliers, ]
df_by_slice_M0M1 <- df_by_slice_M0M1[!M0M1_rows_with_super_outliers, ]

# Delete temp files to release RAM usage
rm(data_sheet_1)
rm(data_sheet_2)
rm(temp_df_by_slice)
rm(temp_df_by_slice_M0M1)
rm(temp_df_by_slice_M0M1M2)
rm(M0M1_rows_with_super_outliers)
rm(M0M1M2_rows_with_super_outliers)
rm(read_first_sheet)
rm(read_second_sheet)

#############################################

# # Optional: show a summary which compounds are often being filtered due to low data quality
# # Function to normalize compound names
# normalize_compound_name <- function(name) {
#   gsub("\\s*(possible|Possible|peak|Peak|minor|Minor)\\s*\\d*", "", name, ignore.case = TRUE)
# }
# # Apply the normalization to the compound names in the original and removed dataframes
# df_raw$compound_normalized <- sapply(df_raw$compound, normalize_compound_name)
# df_removed$compound_normalized <- sapply(df_removed$compound, normalize_compound_name)
# # Summary for 'compound_normalized' in the original and removed dataframes
# compound_summary_original <- table(df_raw$compound_normalized)
# compound_summary_removed <- table(df_removed$compound_normalized)
# # For each normalized compound, calculate how many rows are in the original and how many are removed
# compound_analysis <- merge(data.frame(compound_summary_original), data.frame(compound_summary_removed), by = "Var1", all.x = TRUE)
# names(compound_analysis) <- c("Compound", "OriginalCount", "RemovedCount")
# compound_analysis$RemovedCount[is.na(compound_analysis$RemovedCount)] <- 0
# # Add a column for the percentage of removed rows
# compound_analysis$RemovedPercentage <- (compound_analysis$RemovedCount / compound_analysis$OriginalCount) * 100
# compound_analysis$RemovedPercentage <- round(compound_analysis$RemovedPercentage, 2)

#############################################

# Find out above which intensity the M1/M0 ratio is usually accurate (< 5% percentage error)
# Therefore if a metabolite has at least one M0, M1, or M2 signal with intensity higher than that threshold,
# That intensity would be considered to be the ground true intensity signal and other Mn (usually M1 abd M2) with
# use theorotical MID and this found ground true intensity to calculate what the true intensity should be.
# Create a modified copy of df_by_slice_M0M1M2
df_by_slice_M0M1M2_modified <- df_by_slice_M0M1M2 %>%
  mutate(log_m0_intensity_sub_b = log10(m0_intensity_sub_b),
         log_m1_intensity_sub_b = log10(m1_intensity_sub_b),
         min_log_m0m1_intensity_sub_b = pmin(log_m0_intensity_sub_b, log_m1_intensity_sub_b),
         source = "M0M1M2")
# Create a modified copy of df_by_slice_M0M1
df_by_slice_M0M1_modified <- df_by_slice_M0M1 %>%
  mutate(log_m0_intensity_sub_b = log10(m0_intensity_sub_b),
         log_m1_intensity_sub_b = log10(m1_intensity_sub_b),
         min_log_m0m1_intensity_sub_b = pmin(log_m0_intensity_sub_b, log_m1_intensity_sub_b),
         source = "M0M1")
# Combine the two modified data frames
df_by_slice_combined <- bind_rows(df_by_slice_M0M1M2_modified, df_by_slice_M0M1_modified)

# Define the expression to plot
expression <- ((df_by_slice_combined$m1_intensity_sub_b / df_by_slice_combined$m0_intensity_sub_b) - 
               (df_by_slice_combined$m1_predicted / df_by_slice_combined$m0_predicted)) / 
              (df_by_slice_combined$m1_predicted / df_by_slice_combined$m0_predicted)
# Add the expression as a new column in the combined data frame
df_by_slice_combined <- df_by_slice_combined %>%
  mutate(expression = expression)
# Define the ranges
ranges <- seq(2, 11, by = 0.5)
# Function to calculate medians for a given intensity_sub_b column and expression
calculate_medians <- function(df, log_intensity_col, expression_col, ranges) {
  sapply(1:(length(ranges) - 1), function(i) {
    lower_bound <- ranges[i]
    upper_bound <- ranges[i + 1]
    subset <- df %>%
      filter(!!sym(log_intensity_col) >= lower_bound & !!sym(log_intensity_col) < upper_bound)
    median(subset[[expression_col]], na.rm = TRUE)
  })
}
# Calculate medians for min_log_intensity_sub_b and the expression
median_expression_min_log <- calculate_medians(df_by_slice_combined, "min_log_m0m1_intensity_sub_b", "expression", ranges)
# Create a data frame with the results
percentage_error_vs_intensity_results <- data.frame(
  range_start = ranges[1:(length(ranges) - 1)],
  range_end = ranges[2:length(ranges)],
  median_expression = median_expression_min_log
)
percentage_error_vs_intensity_results <- percentage_error_vs_intensity_results[!is.na(percentage_error_vs_intensity_results$median_expression), ]
# Print the results
print(percentage_error_vs_intensity_results)

#############################################
### Based on the printed results, when both M0 and M1 intensity are >= 1e7, their 
### ratio in a slice is usually accurate, here is some additional confirmation
#############################################

range_larger7 <- c(7.0, 11.0)
median_larger7 <- calculate_medians(df_by_slice_combined, "min_log_m0m1_intensity_sub_b", "expression", range_larger7)


# Function to calculate sdfor a given intensity_sub_b column and expression
calculate_sd <- function(df, log_intensity_col, expression_col, ranges) {
  sapply(1:(length(ranges) - 1), function(i) {
    lower_bound <- ranges[i]
    upper_bound <- ranges[i + 1]
    subset <- df %>%
      filter(!!sym(log_intensity_col) >= lower_bound & !!sym(log_intensity_col) < upper_bound)
    sd(subset[[expression_col]], na.rm = TRUE)
  })
}
sd_larger7 <- calculate_sd(df_by_slice_combined, "min_log_m0m1_intensity_sub_b", "expression", range_larger7)

print(paste0("We can observe that the > E7 M1 measured intensity is very similar to the true intensity, error = ", median_larger7))
print(paste0("and sd = ", sd_larger7))
# The 0.27 sd may looks bad in the first glance. However, a typical peak usually has
# 50 slices, so the final integrated ration will have sd = 0.04 (0.27*50^0.5)

# A futher look on historgram to find out the distribution is 
ratios <- suppressWarnings(as.numeric(unlist(expression)))
larger7v <- df_by_slice_combined$min_log_m0m1_intensity_sub_b >= 7
larger7v[is.na(larger7v)] <- FALSE

plot <- unlist(expression)[df_by_slice_combined$min_log_m0m1_intensity_sub_b >= 7]
plot <- plot[is.finite(plot) & plot >= -1 & plot <= 2] 

hist(plot, breaks = seq(-1, 2, length.out = 201),
     main = "Expression (min_log_m0m1_intensity_sub_b ≥ 7)",
     xlab = "Expression")
abline(v = 0, col = "red", lwd = 2, lty = 2)

#############################################################
# Next we we want to confirm if we we choose 1e7 to divide low and high intensity slices data, how many slices will we have

# Identify the minimal intensity that makes ratio of two Mn of a slice to be accurate
target_rows <- which(abs(percentage_error_vs_intensity_results$median_expression) <= percentage_error_cutoff)
if (length(target_rows) > 0) {
  rle_diff <- rle(c(1, diff(target_rows)) == 1)
  lengths <- rle_diff$lengths
  values <- rle_diff$values
  run_ends <- cumsum(lengths)
  run_starts <- c(1, head(run_ends + 1, -1))
  adjacent_runs <- mapply(function(start, len) {
    if (values[start]) {
      target_rows[run_starts[start]:(run_ends[start])]
    } else {
      NULL
    }
  }, seq_along(values), lengths, SIMPLIFY = FALSE)
  adjacent_runs <- Filter(Negate(is.null), adjacent_runs)
  adjacent_runs <- Filter(function(x) length(x) > 0, adjacent_runs)
  if (length(adjacent_runs) > 0) {
    longest_run <- adjacent_runs[[which.max(sapply(adjacent_runs, length))]]
    subset_df <- percentage_error_vs_intensity_results[longest_run, ]
    min_range_row <- subset_df[which.min(subset_df$range_start), ]
    good_int_cutoff <- 10^min_range_row$range_start
  } else {
    good_int_cutoff <- NA
  }
} else {
  good_int_cutoff <- NA  # No rows with median_expression <= 0.05
}

# Next to plot for each slice what is the max intensity among m0m1m2 to check if 
# Function to calculate the maximum mn intensity, considering the source
# calculate_max_mn_intensity <- function(row) {
#   if (row$source == "M0M1M2") {
#     return(max(c(row$m0_intensity, row$m1_intensity, row$m2_intensity), na.rm = TRUE))
#   } else {
#     return(max(c(row$m0_intensity, row$m1_intensity), na.rm = TRUE))
#   }
# }
# Apply the function to each row to calculate the max mn intensity

df_by_slice_combined <- df_by_slice_combined %>%
  rowwise() %>%
  #mutate(max_mn_intensity = calculate_max_mn_intensity(cur_data())) %>%
  #mutate(max_mn_intensity = calculate_max_mn_intensity(row = pick(everything()))) %>%
  mutate(max_mn_intensity = if_else(source == "M0M1M2",
                                    max(c(m0_intensity, m1_intensity, m2_intensity), na.rm = TRUE),
                                    max(c(m0_intensity, m1_intensity), na.rm = TRUE))) %>%
  ungroup()

# Calculate the maximum intensity for each row in df_by_slice_combined




max_intensity_sub_b <- apply(df_by_slice_combined[, c("m0_intensity_sub_b", "m1_intensity_sub_b", "m2_intensity_sub_b")], 1, max, na.rm = TRUE)
# Apply log10 transformation to the maximum intensity values
log_max_intensity_sub_b <- log10(max_intensity_sub_b)
# Sort the log-transformed maximum intensity values
sorted_log_max_intensity_sub_b <- sort(log_max_intensity_sub_b)
# Create the final combined data frame df_by_slice_combined with necessary columns
df_by_slice_combined <- df_by_slice_combined %>%
  select(-source)  # Remove the source column if not needed

log_cutoff <- log10(good_int_cutoff)
# Create data frame for plotting and count and calculate statistics
df_ggplot <- data.frame(
  Index = seq_along(sorted_log_max_intensity_sub_b),
  logMaxIntensity = sorted_log_max_intensity_sub_b
)

num_above_cutoff <- sum(df_ggplot$logMaxIntensity >= log_cutoff)
total_points <- nrow(df_ggplot)
per_above_int_cutoff <- num_above_cutoff / total_points * 100

# Print both number and percentage
cat(sprintf(
  "Number of points above log10(good_int_cutoff) = %.2f: %d out of %d (%.2f%%)\n",
  log_cutoff,
  num_above_cutoff,
  total_points,
  per_above_int_cutoff
))

# Create the ggplot
ggplot(df_ggplot, aes(x = Index, y = logMaxIntensity)) +
  geom_line() +
  geom_hline(yintercept = log_cutoff, color = "red", linetype = "dashed", size = 1) +
  labs(
    x = "Index",
    y = "log10(Max intensity_sub_b)",
    title = "Spread of log10(Max intensity_sub_b)",
    subtitle = sprintf("Red dashed line: log10(good_int_cutoff) = %.2f", log_cutoff)
  ) +
  theme_minimal()

# rm all temp files to release RAM
rm(df_by_slice_M0M1M2_modified)
rm(df_by_slice_M0M1_modified)
rm(expression)
rm(log_max_intensity_sub_b)
rm(max_intensity_sub_b)
rm(calculate_medians)
# rm(calculate_max_mn_intensity)
rm(df_ggplot)

#######################################

# Additional check on by slice data row that the observed M0M1M2 intensity if very far away from true ratio: 
# If the highest Mn ratio is at least 0.1 higher than other Mn ratios but the observed intensity is not the max among all mn intensity, 
# that row of by slice data is discarded before training as the true intensity cannot be determined
# Function to check the condition
check_highest_predicted <- function(row) {
  max_predicted <- which.max(c(row$m0_predicted, row$m1_predicted, row$m2_predicted))
  predicted_values <- c(row$m0_predicted, row$m1_predicted, row$m2_predicted)
  sorted_indices <- order(predicted_values, decreasing = TRUE)
  max_index <- sorted_indices[1]
  second_max_index <- sorted_indices[2]
  
  if (row$integrated_to_area == 1) {
    if (predicted_values[max_index] - predicted_values[second_max_index] < 0.1) {
      if ((max_index == 1 && (row$m0_intensity >= row$m1_intensity & row$m0_intensity >= row$m2_intensity)) ||
          (max_index == 2 && (row$m1_intensity >= row$m0_intensity & row$m1_intensity >= row$m2_intensity)) ||
          (max_index == 3 && (row$m2_intensity >= row$m0_intensity & row$m2_intensity >= row$m1_intensity)) ||
          (second_max_index == 1 && (row$m0_intensity >= row$m1_intensity & row$m0_intensity >= row$m2_intensity)) ||
          (second_max_index == 2 && (row$m1_intensity >= row$m0_intensity & row$m1_intensity >= row$m2_intensity)) ||
          (second_max_index == 3 && (row$m2_intensity >= row$m0_intensity & row$m2_intensity >= row$m1_intensity))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      if ((max_index == 1 && (row$m0_intensity >= row$m1_intensity & row$m0_intensity >= row$m2_intensity)) ||
          (max_index == 2 && (row$m1_intensity >= row$m0_intensity & row$m1_intensity >= row$m2_intensity)) ||
          (max_index == 3 && (row$m2_intensity >= row$m0_intensity & row$m2_intensity >= row$m1_intensity))) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  } else {
    if (predicted_values[max_index] - predicted_values[second_max_index] < 0.1) {
      if ((max_index == 1 && row$m0_intensity >= row$m1_intensity) ||
          (max_index == 2 && row$m1_intensity >= row$m0_intensity) ||
          (second_max_index == 1 && row$m0_intensity >= row$m1_intensity) ||
          (second_max_index == 2 && row$m1_intensity >= row$m0_intensity)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      if ((max_index == 1 && row$m0_intensity >= row$m1_intensity) ||
          (max_index == 2 && row$m1_intensity >= row$m0_intensity)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
  }
}

# Create a new dataframe with the additional check_result column
df_by_slice_combined_checked <- df_by_slice_combined %>%
  rowwise() %>%
  #mutate(check_result = check_highest_predicted(cur_data())) %>%
  #mutate(check_result = check_highest_predicted(pick(everything())))

  mutate(check_result = {
    predicted_values <- c(m0_predicted, m1_predicted, m2_predicted)
    sorted_indices <- order(predicted_values, decreasing = TRUE)
    max_index <- sorted_indices[1]
    second_max_index <- sorted_indices[2]
    
    case_when(
      # integrated_to_area == 1 and difference < 0.1
      integrated_to_area == 1 & (predicted_values[max_index] - predicted_values[second_max_index] < 0.1) ~
        (max_index == 1 && (m0_intensity >= m1_intensity & m0_intensity >= m2_intensity)) ||
        (max_index == 2 && (m1_intensity >= m0_intensity & m1_intensity >= m2_intensity)) ||
        (max_index == 3 && (m2_intensity >= m0_intensity & m2_intensity >= m1_intensity)) ||
        (second_max_index == 1 && (m0_intensity >= m1_intensity & m0_intensity >= m2_intensity)) ||
        (second_max_index == 2 && (m1_intensity >= m0_intensity & m1_intensity >= m2_intensity)) ||
        (second_max_index == 3 && (m2_intensity >= m0_intensity & m2_intensity >= m1_intensity)),
      
      # integrated_to_area == 1 and difference >= 0.1
      integrated_to_area == 1 ~
        (max_index == 1 && (m0_intensity >= m1_intensity & m0_intensity >= m2_intensity)) ||
        (max_index == 2 && (m1_intensity >= m0_intensity & m1_intensity >= m2_intensity)) ||
        (max_index == 3 && (m2_intensity >= m0_intensity & m2_intensity >= m1_intensity)),
      
      # integrated_to_area != 1 and difference < 0.1
      predicted_values[max_index] - predicted_values[second_max_index] < 0.1 ~
        (max_index == 1 && m0_intensity >= m1_intensity) ||
        (max_index == 2 && m1_intensity >= m0_intensity) ||
        (second_max_index == 1 && m0_intensity >= m1_intensity) ||
        (second_max_index == 2 && m1_intensity >= m0_intensity),
      
      # integrated_to_area != 1 and difference >= 0.1
      TRUE ~
        (max_index == 1 && m0_intensity >= m1_intensity) ||
        (max_index == 2 && m1_intensity >= m0_intensity)
    )
  }) %>%
  ungroup()



# Filter the rows where check_result is FALSE and save in a separate dataframe
df_by_slice_combined_false <- df_by_slice_combined_checked %>%
  filter(check_result == FALSE)
df_by_slice_combined <- df_by_slice_combined_checked %>%
  filter(check_result == TRUE) %>%
  select(-check_result)

rm(df_by_slice_combined_checked)
rm(check_highest_predicted)

##############################################################################################################

# Create df_by_mn_high and df_by_mn_low based on if its signal is above or below the determined good_int_cutoff which is used for rf training
# For each slice of df_by_mn_high, the highest mn signal will be used to calculate the true intensity of other Mn signals
# For each slice of df_by_mn_low, they will not be used in the first model training of step 2c as its highest intensity is still off,
# Only after the intial traing, these data will be added for further training, after its true inetnsity be predicted by the inital model
# Function to generate new rows with the highest intensity parameter
generate_new_rows <- function(compound, sample, formula, rt_left, rt_peak, rt_right, slice_rt, slice, total_slices, slice_peak_height, reference_mz, precursor_mz, m1_mz, m2_mz, total_ion, total_ion_smoothed, m0_intensity, m1_intensity, m2_intensity = NA, m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed = NA, m0_baseline, m1_baseline, m2_baseline = NA, m0_noise, m1_noise, m2_noise = NA, mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity = NA, m3_intensity = NA, m3_intensity_smoothed = NA, m0_intensity_sub_b, m1_intensity_sub_b, m2_intensity_sub_b = NA, m0_intensity_smoothed_sub_b, m1_intensity_smoothed_sub_b, m2_intensity_smoothed_sub_b = NA, m0_predicted, m1_predicted, m2_predicted = NA, integrated_to_area, source) {
  highest_original_intensity <- max(m0_intensity, m1_intensity, m2_intensity, na.rm = TRUE)
  
  if (source == "M0M1M2") {
    if (highest_original_intensity >= good_int_cutoff) {
      if (m0_intensity == highest_original_intensity) {
        new_rows <- data.frame(
          compound = rep(compound, 3),
          sample = rep(sample, 3),
          formula = rep(formula, 3),
          mn = c(0, 1, 2),
          rt_left = rep(rt_left, 3),
          rt_peak = rep(rt_peak, 3),
          rt_right = rep(rt_right, 3),
          slice_rt = rep(slice_rt, 3),
          slice = rep(slice, 3),
          total_slices = rep(total_slices, 3),
          slice_peak_height = rep(slice_peak_height, 3),
          reference_mz = rep(reference_mz, 3),
          mn_mz = c(precursor_mz, m1_mz, m2_mz),
          total_ion = rep(total_ion, 3),
          total_ion_smoothed = rep(total_ion_smoothed, 3),
          mn_intensity = c(m0_intensity, m1_intensity, m2_intensity),
          mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed),
          mn_baseline = c(m0_baseline, m1_baseline, m2_baseline),
          mn_noise = c(m0_noise, m1_noise, m2_noise),
          mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity),
          mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity),
          mplus1_intensity = c(m1_intensity, m2_intensity, m3_intensity),
          mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed, m3_intensity_smoothed),
          intensity_predicted = c(m0_intensity_sub_b, (m0_intensity_sub_b * m1_predicted) / m0_predicted, (m0_intensity_sub_b * m2_predicted) / m0_predicted),
          intensity_smoothed_predicted = c(m0_intensity_smoothed_sub_b, (m0_intensity_smoothed_sub_b * m1_predicted) / m0_predicted, (m0_intensity_smoothed_sub_b * m2_predicted) / m0_predicted),
          ratio = c(NA, NA, NA),
          highest_original_intensity = highest_original_intensity
        )
      } else if (m1_intensity == highest_original_intensity) {
        new_rows <- data.frame(
          compound = rep(compound, 3),
          sample = rep(sample, 3),
          formula = rep(formula, 3),
          mn = c(0, 1, 2),
          rt_left = rep(rt_left, 3),
          rt_peak = rep(rt_peak, 3),
          rt_right = rep(rt_right, 3),
          slice_rt = rep(slice_rt, 3),
          slice = rep(slice, 3),
          total_slices = rep(total_slices, 3),
          slice_peak_height = rep(slice_peak_height, 3),
          reference_mz = rep(reference_mz, 3),
          mn_mz = c(precursor_mz, m1_mz, m2_mz),
          total_ion = rep(total_ion, 3),
          total_ion_smoothed = rep(total_ion_smoothed, 3),
          mn_intensity = c(m0_intensity, m1_intensity, m2_intensity),
          mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed),
          mn_baseline = c(m0_baseline, m1_baseline, m2_baseline),
          mn_noise = c(m0_noise, m1_noise, m2_noise),
          mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity),
          mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity),
          mplus1_intensity = c(m1_intensity, m2_intensity, m3_intensity),
          mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed, m3_intensity_smoothed),
          intensity_predicted = c((m1_intensity_sub_b * m0_predicted) / m1_predicted, m1_intensity_sub_b, (m1_intensity_sub_b * m2_predicted) / m1_predicted),
          intensity_smoothed_predicted = c((m1_intensity_smoothed_sub_b * m0_predicted) / m1_predicted, m1_intensity_smoothed_sub_b, (m1_intensity_smoothed_sub_b * m2_predicted) / m1_predicted),
          ratio = c(NA, NA, NA),
          highest_original_intensity = highest_original_intensity
        )
      } else if (m2_intensity == highest_original_intensity) {
        new_rows <- data.frame(
          compound = rep(compound, 3),
          sample = rep(sample, 3),
          formula = rep(formula, 3),
          mn = c(0, 1, 2),
          rt_left = rep(rt_left, 3),
          rt_peak = rep(rt_peak, 3),
          rt_right = rep(rt_right, 3),
          slice_rt = rep(slice_rt, 3),
          slice = rep(slice, 3),
          total_slices = rep(total_slices, 3),
          slice_peak_height = rep(slice_peak_height, 3),
          reference_mz = rep(reference_mz, 3),
          mn_mz = c(precursor_mz, m1_mz, m2_mz),
          total_ion = rep(total_ion, 3),
          total_ion_smoothed = rep(total_ion_smoothed, 3),
          mn_intensity = c(m0_intensity, m1_intensity, m2_intensity),
          mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed),
          mn_baseline = c(m0_baseline, m1_baseline, m2_baseline),
          mn_noise = c(m0_noise, m1_noise, m2_noise),
          mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity),
          mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity),
          mplus1_intensity = c(m1_intensity, m2_intensity, m3_intensity),
          mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed, m3_intensity_smoothed),
          intensity_predicted = c((m2_intensity_sub_b * m0_predicted) / m2_predicted, (m2_intensity_sub_b * m1_predicted) / m2_predicted, m2_intensity_sub_b),
          intensity_smoothed_predicted = c((m2_intensity_smoothed_sub_b * m0_predicted) / m2_predicted, (m2_intensity_smoothed_sub_b * m1_predicted) / m2_predicted, m2_intensity_smoothed_sub_b),
          ratio = c(NA, NA, NA),
          highest_original_intensity = highest_original_intensity
        )
      }
    } else {
      new_rows <- data.frame(
        compound = rep(compound, 3),
        sample = rep(sample, 3),
        formula = rep(formula, 3),
        mn = c(0, 1, 2),
        rt_left = rep(rt_left, 3),
        rt_peak = rep(rt_peak, 3),
        rt_right = rep(rt_right, 3),
        slice_rt = rep(slice_rt, 3),
        slice = rep(slice, 3),
        total_slices = rep(total_slices, 3),
        slice_peak_height = rep(slice_peak_height, 3),
        reference_mz = rep(reference_mz, 3),
        mn_mz = c(precursor_mz, m1_mz, m2_mz),
        total_ion = rep(total_ion, 3),
        total_ion_smoothed = rep(total_ion_smoothed, 3),
        mn_intensity = c(m0_intensity, m1_intensity, m2_intensity),
        mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed),
        mn_baseline = c(m0_baseline, m1_baseline, m2_baseline),
        mn_noise = c(m0_noise, m1_noise, m2_noise),
        mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity),
        mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity),
        mplus1_intensity = c(m1_intensity, m2_intensity, m3_intensity),
        mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed, m3_intensity_smoothed),
        intensity_predicted = c(NA, NA, NA),
        intensity_smoothed_predicted = c(NA, NA, NA),
        ratio = c(m0_predicted, m1_predicted, m2_predicted),
        highest_original_intensity = highest_original_intensity
      )
    }
  } else if (source == "M0M1") {
    highest_original_intensity <- max(m0_intensity, m1_intensity, na.rm = TRUE)
    if (highest_original_intensity >= good_int_cutoff) {
      if (m0_intensity >= m1_intensity) {
        new_rows <- data.frame(
          compound = rep(compound, 2),
          sample = rep(sample, 2),
          formula = rep(formula, 2),
          mn = c(0, 1),
          rt_left = rep(rt_left, 2),
          rt_peak = rep(rt_peak, 2),
          rt_right = rep(rt_right, 2),
          slice_rt = rep(slice_rt, 2),
          slice = rep(slice, 2),
          total_slices = rep(total_slices, 2),
          slice_peak_height = rep(slice_peak_height, 2),
          reference_mz = rep(reference_mz, 2),
          mn_mz = c(precursor_mz, m1_mz),
          total_ion = rep(total_ion, 2),
          total_ion_smoothed = rep(total_ion_smoothed, 2),
          mn_intensity = c(m0_intensity, m1_intensity),
          mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed),
          mn_baseline = c(m0_baseline, m1_baseline),
          mn_noise = c(m0_noise, m1_noise),
          mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity),
          mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity),
          mplus1_intensity = c(m1_intensity, m2_intensity),
          mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed),
          intensity_predicted = c(m0_intensity_sub_b, (m0_intensity_sub_b * m1_predicted) / m0_predicted),
          intensity_smoothed_predicted = c(m0_intensity_smoothed_sub_b, (m0_intensity_smoothed_sub_b * m1_predicted) / m0_predicted),
          ratio = c(NA, NA),
          highest_original_intensity = highest_original_intensity
        )
      } else if (m1_intensity >= m0_intensity) {
        new_rows <- data.frame(
          compound = rep(compound, 2),
          sample = rep(sample, 2),
          formula = rep(formula, 2),
          mn = c(0, 1),
          rt_left = rep(rt_left, 2),
          rt_peak = rep(rt_peak, 2),
          rt_right = rep(rt_right, 2),
          slice_rt = rep(slice_rt, 2),
          slice = rep(slice, 2),
          total_slices = rep(total_slices, 2),
          slice_peak_height = rep(slice_peak_height, 2),
          reference_mz = rep(reference_mz, 2),
          mn_mz = c(precursor_mz, m1_mz),
          total_ion = rep(total_ion, 2),
          total_ion_smoothed = rep(total_ion_smoothed, 2),
          mn_intensity = c(m0_intensity, m1_intensity),
          mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed),
          mn_baseline = c(m0_baseline, m1_baseline),
          mn_noise = c(m0_noise, m1_noise),
          mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity),
          mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity),
          mplus1_intensity = c(m1_intensity, m2_intensity),
          mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed),
          intensity_predicted = c((m1_intensity_sub_b * m0_predicted) / m1_predicted, m1_intensity_sub_b),
          intensity_smoothed_predicted = c((m1_intensity_smoothed_sub_b * m0_predicted) / m1_predicted, m1_intensity_smoothed_sub_b),
          ratio = c(NA, NA),
          highest_original_intensity = highest_original_intensity
        )
      }
    } else {
      new_rows <- data.frame(
        compound = rep(compound, 2),
        sample = rep(sample, 2),
        formula = rep(formula, 2),
        mn = c(0, 1),
        rt_left = rep(rt_left, 2),
        rt_peak = rep(rt_peak, 2),
        rt_right = rep(rt_right, 2),
        slice_rt = rep(slice_rt, 2),
        slice = rep(slice, 2),
        total_slices = rep(total_slices, 2),
        slice_peak_height = rep(slice_peak_height, 2),
        reference_mz = rep(reference_mz, 2),
        mn_mz = c(precursor_mz, m1_mz),
        total_ion = rep(total_ion, 2),
        total_ion_smoothed = rep(total_ion_smoothed, 2),
        mn_intensity = c(m0_intensity, m1_intensity),
        mn_intensity_smoothed = c(m0_intensity_smoothed, m1_intensity_smoothed),
        mn_baseline = c(m0_baseline, m1_baseline),
        mn_noise = c(m0_noise, m1_noise),
        mnminus1_to_mn_intensity = c(mminus1_to_m0_intensity, m0_to_m1_intensity),
        mn_to_mplus1_intensity = c(m0_to_m1_intensity, m1_to_m2_intensity),
        mplus1_intensity = c(m1_intensity, m2_intensity),
        mplus1_intensity_smoothed = c(m1_intensity_smoothed, m2_intensity_smoothed),
        intensity_predicted = c(NA, NA),
        intensity_smoothed_predicted = c(NA, NA),
        ratio = c(m0_predicted / (m0_predicted + m1_predicted), m1_predicted / (m0_predicted + m1_predicted)),
        highest_original_intensity = highest_original_intensity
      )
    }
  }
  return(new_rows)
}

# Add a source column to tell the function generate_new_rows how to process
df_by_slice_M0M1M2 <- df_by_slice_M0M1M2 %>%
  mutate(source = "M0M1M2")
df_by_slice_M0M1 <- df_by_slice_M0M1 %>%
  mutate(source = "M0M1")
# Select the relevant columns from the original data frames including the new parameters
df_selected_M0M1M2 <- df_by_slice_M0M1M2 %>%
  select(compound, sample, formula, rt_left, rt_peak, rt_right, slice_rt, slice, total_slices, slice_peak_height, reference_mz, precursor_mz, m1_mz, m2_mz, total_ion, total_ion_smoothed,
         m0_intensity, m1_intensity, m2_intensity, m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed, m0_baseline, m1_baseline, m2_baseline, m0_noise, m1_noise, m2_noise,
         mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity, m3_intensity, m3_intensity_smoothed, m0_intensity_sub_b, m1_intensity_sub_b, m2_intensity_sub_b,
         m0_intensity_smoothed_sub_b, m1_intensity_smoothed_sub_b, m2_intensity_smoothed_sub_b, m0_predicted, m1_predicted, m2_predicted, integrated_to_area, source)
df_selected_M0M1 <- df_by_slice_M0M1 %>%
  select(compound, sample, formula, rt_left, rt_peak, rt_right, slice_rt, slice, total_slices, slice_peak_height, reference_mz, precursor_mz, m1_mz, total_ion, total_ion_smoothed,
         m0_intensity, m1_intensity, m2_intensity, m0_intensity_smoothed, m1_intensity_smoothed, m2_intensity_smoothed, m0_baseline, m1_baseline, m0_noise, m1_noise,
         mminus1_to_m0_intensity, m0_to_m1_intensity, m1_to_m2_intensity, m2_to_m3_intensity, m0_intensity_sub_b, m1_intensity_sub_b, m0_intensity_smoothed_sub_b, m1_intensity_smoothed_sub_b, m0_predicted, m1_predicted, integrated_to_area, source)
# Apply the function generate_new_rows to each row of the selected data frames using pmap
df_by_mn_M0M1M2 <- pmap_dfr(df_selected_M0M1M2, generate_new_rows)
df_by_mn_M0M1 <- pmap_dfr(df_selected_M0M1, generate_new_rows)
# Combine the new rows from both data frames
df_by_mn <- bind_rows(df_by_mn_M0M1M2, df_by_mn_M0M1)
# Separate the new rows into high and low intensity data frames
df_by_mn_high <- df_by_mn %>% filter(highest_original_intensity >= good_int_cutoff)
df_by_mn_low <- df_by_mn %>% filter(highest_original_intensity < good_int_cutoff)

rm(df_selected_M0M1)
rm(df_selected_M0M1M2)
rm(generate_new_rows)

#######################################

# # Optional Check: Identify rows where M0_predicted is less than M1_predicted or less than M2_predicted
# rows_M0_less_than_M1_or_M2 <- which(df_by_slice_M0M1M2$m0_predicted < df_by_slice_M0M1M2$m1_predicted |
#                                     df_by_slice_M0M1M2$m0_predicted < df_by_slice_M0M1M2$m2_predicted)
# 
# cat("Rows where M0_predicted is less than M1_predicted or less than M2_predicted:\n")
# print(rows_M0_less_than_M1_or_M2)
# cat("Details of rows where M0_predicted is less than M1_predicted or less than M2_predicted:\n")
# print(df_by_slice_M0M1M2[rows_M0_less_than_M1_or_M2, ])

###########################################

# Should check which cols have NA values and replace them!
high_columns_with_na <- names(df_by_mn_high)[sapply(df_by_mn_high, function(x) sum(is.na(x)) > 0)]
low_columns_with_na <- names(df_by_mn_low)[sapply(df_by_mn_low, function(x) sum(is.na(x)) > 0)]
# Calculate the median of mn_noise excluding NA values
median_mn_noise <- median(df_by_mn$mn_noise, na.rm = TRUE)
# Replace NA values in mn_noise with the calculated median
df_by_mn_high$mn_noise[is.na(df_by_mn_high$mn_noise)] <- median_mn_noise
df_by_mn_low$mn_noise[is.na(df_by_mn_low$mn_noise)] <- median_mn_noise
df_by_mn_high$rt_width <- df_by_mn_high$rt_right - df_by_mn_high$rt_left
df_by_mn_low$rt_width <- df_by_mn_low$rt_right - df_by_mn_low$rt_left
df_by_mn_high$slice_diff_to_m0peak <- df_by_mn_high$slice - df_by_mn_high$slice_peak_height
df_by_mn_low$slice_diff_to_m0peak <- df_by_mn_low$slice - df_by_mn_low$slice_peak_height
df_by_mn_high$mz_to_center <- abs(df_by_mn_high$mn_mz - (69+360)/2)
df_by_mn_low$mz_to_center <- abs(df_by_mn_low$mn_mz - (69+360)/2)
df_by_mn_high$intensity_over_TIC  <- df_by_mn_high$mn_intensity/df_by_mn_high$total_ion
df_by_mn_low$intensity_over_TIC  <- df_by_mn_low$mn_intensity/df_by_mn_low$total_ion
df_by_mn_high$intensity_over_TIC_smoothed  <- df_by_mn_high$mn_intensity_smoothed/df_by_mn_high$total_ion_smoothed
df_by_mn_low$intensity_over_TIC_smoothed  <- df_by_mn_low$mn_intensity_smoothed/df_by_mn_low$total_ion_smoothed
df_by_mn_high$intensity_minus_bg  <- df_by_mn_high$mn_intensity - df_by_mn_high$mn_baseline
df_by_mn_low$intensity_minus_bg  <- df_by_mn_low$mn_intensity - df_by_mn_low$mn_baseline
df_by_mn_high$intensity_minus_bg_smoothed  <- df_by_mn_high$mn_intensity_smoothed - df_by_mn_high$mn_baseline
df_by_mn_low$intensity_minus_bg_smoothed  <- df_by_mn_low$mn_intensity_smoothed - df_by_mn_low$mn_baseline
df_by_mn_high$intensity_over_bg  <- df_by_mn_high$mn_intensity / (df_by_mn_high$mn_baseline + epsilon)
df_by_mn_low$intensity_over_bg  <- df_by_mn_low$mn_intensity / (df_by_mn_low$mn_baseline + epsilon)
df_by_mn_high$intensity_over_ns  <- df_by_mn_high$mn_intensity / (df_by_mn_high$mn_noise + epsilon)
df_by_mn_low$intensity_over_ns  <- df_by_mn_low$mn_intensity / (df_by_mn_low$mn_noise + epsilon)

# Function to count elements number from formula, be aware that it has different name and output format as 1c functions
countelements <- function(inputformula) {
  elements <- c("H", "C", "N", "O", "S", "P", "Si")
  counts <- setNames(numeric(length(elements)), elements)
  formula <- gsub("Si", "Xx", inputformula)
  
  for (element in elements) {
    search_el <- if (element == "Si") "Xx" else element
    pattern <- paste0(search_el, "(\\d*)")
    raw_matches <- regmatches(formula, gregexpr(pattern, formula, perl = TRUE))[[1]]
    
    if (length(raw_matches) > 0) {
      # Extract numbers; implied '1' if empty eg: CO2 means C1O2
      nums <- vapply(raw_matches, function(x) {
        val <- sub(search_el, "", x)
        if (val == "") 1 else as.numeric(val)
      }, numeric(1))
      counts[element] <- sum(nums)
    }
  }
  return(counts)
}

# Function to add element count cols
add_element_counts <- function(df) {
  required_columns <- c("H", "C", "N", "O", "S", "P", "Si")
  df <- df[, !(names(df) %in% required_columns), drop = FALSE]
  element_counts <- t(sapply(df$formula, countelements))
  element_counts <- as.data.frame(element_counts)
  df <- cbind(df, element_counts)
  df$atom_number <- rowSums(df[required_columns], na.rm = TRUE)
  return(df)
}

df_by_mn_high <- add_element_counts(df_by_mn_high)
df_by_mn_low <- add_element_counts(df_by_mn_low)

rm(add_element_counts)
rm(countelements)

##########################################################

# Now we will save and summary the preprocessed data
save(df, file = "df_2b.RData")
save(df_by_mn_high, file = "df_by_mn_high_2b.RData")
save(df_by_mn_low, file = "df_by_mn_low_2b.RData")
df_by_slice <- df_by_slice_combined
save(df_by_slice, file = "df_by_slice_2b.RData")

# save space for github upload
rm(df_by_slice_combined) 
rm(df_raw)
rm(df_removed)
rm(df_by_slice_raw)

# Print summary
red <- function(x) paste0("\033[31m", x, "\033[0m")
good_cutoff_str <- red(sprintf("%.2e", good_int_cutoff))
error_cutoff_str <- red(sprintf("%.2f", percentage_error_cutoff))
log_cutoff_str <- red(sprintf("%.2f", log_cutoff))
num_above_str <- red(sprintf("%d", num_above_cutoff))
total_str <- red(sprintf("%d", total_points))
percent_str <- red(sprintf("%.2f%%", per_above_int_cutoff))
cat(
  "DATA STATS SUMMARY:\n",
  sprintf("• In this input data set, If two Mn signals of a slice are both >= %s,", good_cutoff_str), "\n",
  sprintf("  their median percentage relative ratio error will be within +-%s.", error_cutoff_str), "\n",
  sprintf("• Therefore, %s is chosen as the cutoff to decide if a Mn signal intesntiy is accurate.", good_cutoff_str), "\n",
  "  Accurate signals need no bias correction and can help compute true intensity\n",
  "  of other Mn signals based on predicted MID (Mass Isotopomer Distribution).\n",
  sprintf("• %s out of %s slices (before ratio checks) have Mn signals above this cutoff (%s).\n", num_above_str, total_str, percent_str),
  sprintf("• %s intensity cutoff based df_by_mn_high and df_by_mn_low data \n", good_cutoff_str),
  "  along with the filtered df and df by slice has been saved,\n",
  "  please continue either step 2cRF or 2cEN bias correction model training.\n"
)


