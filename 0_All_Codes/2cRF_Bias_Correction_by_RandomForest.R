###################################

# File: 2cRF_Bias_Correction_by_RandomForest
# Input: Four files are needed  from step 2b outputs "df_2b.RData" contains the compound level data 
#         "df_by_slice_2b.RData" contains the scan level data, "df_by_mn_high_2b.RData" and "df_by_mn_low_2b.RData" 
#          are mn level data (every scan level data are separated into two or three mn level data based on if only M0 and M1 or 
#          all M0, M1, and M2 stages could be identified, and then sperated again by intensity) 
# Function: The randomforest model training have three iterations: i0 use only high intensity data that at least 
#         one of M0, M1, and M2's measured intensity is known to be accurate so the accurate M1 and M2 intensity could be calculated
#         based on known ratio and used as an answer key to evaluate the i0 model training.
#         The trained i0 model is then used to predict the isotopologue intensity of the low-intensity group data, as none of M0, M1, and M2 in
#        low-intensity data group would likely have high enough intensity to be measured accurately. Based on this predicted intensity, the true intensity of
#         all M0, M1, and M2 intenisty thus could be calculated, which is why we call it as pseudo-labeled. The psedudo-labeled low intenisty group
#        data and the highly accurate high-intensity data are used to train model i1.
#         Model i1 will be used to predict and pseudo-label the low-intensity group data again, and the pseudo-labeled low-intensity data will be used again
#         with the high-intenisty group data to train the final model i2
# Output: Three model files "clf_rf_high_i0_seedxxxx.rds", "clf_rf_all_i1_seedxxxx.rds" , and "clf_rf_all_i2_seedxxxx.rds"
#         Two Rdata files: "df_by_mn_high_2b.RData" and "df_by_mn_low_2b.RData" , which contains data above and below the determined intensity cutoff
#
# Disclaimer: Chatgpt is used to debug this code.

####################################
#####################################

library(readxl) #load xlsx
library(purrr) #load xlsx by iterations
library(openxlsx) #write xlsx
library(magrittr) #pipe
library(tidyr) #df cleaning
library(plyr)  #df splitting and combining
library(dplyr) #df data manipulation
library(caret) #training
library(randomForest) #rf related plot
library(data.table) # handle large data without using too much RAM
library(tidyverse) # ggplot2, dplyr and more
library(ggplot2)
library(writexl)
library(gridExtra) # enable displaying multiple ggplot2 in one page
directory <- getwd()
source("2cRF_Functions.R") # custom functions

# seed_i is the random seed used to separate 80% training and 20% test data by sample
seed_i <- 2024
# total number of model training reruns to check if the result is robust, in each rerun, the used seed will +1 from the previous
num_rerun <- 1
# inner_cv is the within training cross-validation fold, 3 is chose to speed up training
inner_cv <- 3
# Recommended ntree for this data size
ntree <- 2000  # Adjust as necessary, 2000 is a good one for this study but be careful to not over RAM if your computer is not good. You may try differnt values for optimization
mtry_fixed <- 5 # A common starting point for mtry(tree depth) in RF

####################################
#####################################

input_dir <- file.path(directory, "input")
if (!dir.exists(input_dir)) {
  dir.create(input_dir)
}
file_path_df <- file.path(input_dir, "df_2b.RData")
file_path_high <- file.path(input_dir, "df_by_mn_high_2b.RData")
file_path_low  <- file.path(input_dir, "df_by_mn_low_2b.RData")
file_path_slice  <- file.path(input_dir, "df_by_slice_2b.RData")
load(file = file_path_df)
load(file = file_path_high)
load(file = file_path_low)
load(file = file_path_slice)

# Define log file path and start logging
current_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file_name <- sprintf("output_log_seed%d_rerun%d_time%s.txt", seed_i, num_rerun, current_time)
write("", file = log_file_name)
sink(log_file_name, append = TRUE)


#######################################
######################################

for (i in 1:num_rerun) 
  {
    cat("------------------------------------------------------\n")
    cat("Start running trial", i, "out of total trials", num_rerun, "by using random seed", seed_i, "\n")
    cat("The used inner cross-validation fold is", inner_cv, "\n")
    cat("The used number of trees in random forest training is", ntree, "\n")
    cat("The used depth of trees in random forest training is", mtry_fixed, "\n")
    cat("------------------------------------------------------\n")
    
    # Randomize samples into 80% training and 20% test based on sampple name and random seed
    df$sample <- as.character(df$sample)
    set.seed(seed_i)
    unique_values <- unique(df$sample)
    index <- seq(5, length(unique_values), by=5)    # this is a fixed stratified sampling to make unlabeled standard, labeled standard, and unlabeled samples all been represented to show model average performance
    # index <- sample(seq_along(unique_values), length(unique_values) %/% 5)   #  use this line instead of the above line for simple random sampling on training and test set instead of only randomize the model training step
    selected_values <- unique_values[index]
    training_df <- df %>% filter(!sample %in% selected_values)
    test_df <- df %>% filter(sample %in% selected_values)
    training_df_by_mn_high <- df_by_mn_high %>% filter(!sample %in% selected_values)
    training_df_by_mn_low <- df_by_mn_low %>% filter(!sample %in% selected_values)
    test_df_by_mn_high <- df_by_mn_high %>% filter(sample %in% selected_values)
    test_df_by_mn_low <- df_by_mn_low %>% filter(sample %in% selected_values)
    
    # Calculate and print the statistics for M1 M2 RIA
    cat("Uncorrected Mean Absolute Error of M1 RIA (%):", mean(abs(100 * test_df$delta_m1_percentage_sub_b), na.rm = TRUE), "\n")
    cat("Uncorrected Median Absolute Error of M1 RIA (%):", median(abs(100 * test_df$delta_m1_percentage_sub_b), na.rm = TRUE), "\n")
    cat("Uncorrected sd of M1 RIA (%):", sd(100 * test_df$delta_m1_percentage_sub_b, na.rm = TRUE), "\n")
    cat("Uncorrected Mean Absolute Error of M2 RIA (%):", mean(abs(100 * test_df$delta_m2_percentage_sub_b), na.rm = TRUE), "\n")
    cat("Uncorrected Median Absolute Error of M2 RIA (%):", median(abs(100 * test_df$delta_m2_percentage_sub_b), na.rm = TRUE), "\n")
    cat("Uncorrected sd of M2 RIA (%):", sd(100 * test_df$delta_m2_percentage_sub_b, na.rm = TRUE), "\n")
    
    # First convert the format to subset, then subset the needed non-numerical columns into new data.tables
    training_df_by_mn_high <- as.data.table(training_df_by_mn_high)
    training_df_by_mn_low <- as.data.table(training_df_by_mn_low)
    test_df_by_mn_high <- as.data.table(test_df_by_mn_high)
    test_df_by_mn_low <- as.data.table(test_df_by_mn_low)
    training_high_info <- training_df_by_mn_high[, .(compound, sample, formula, slice, intensity_predicted, intensity_smoothed_predicted, ratio, intensity_minus_bg, intensity_minus_bg_smoothed)]
    training_low_info <- training_df_by_mn_low[, .(compound, sample, formula, slice, intensity_predicted, intensity_smoothed_predicted, ratio, intensity_minus_bg, intensity_minus_bg_smoothed)]
    test_high_info <- test_df_by_mn_high[, .(compound, sample, formula, slice, intensity_predicted, intensity_smoothed_predicted, ratio, intensity_minus_bg, intensity_minus_bg_smoothed)]
    test_low_info <- test_df_by_mn_low[, .(compound, sample, formula, slice, intensity_predicted, intensity_smoothed_predicted, ratio, intensity_minus_bg, intensity_minus_bg_smoothed)]
    
    # Remove unused non-numerical columns and force set rest everytying are all numericals
    training_df_by_mn_high[, `:=`(compound = NULL, sample = NULL, formula = NULL, ratio = NULL, intensity_predicted = NULL, intensity_smoothed_predicted = NULL,  intensity_minus_bg = NULL, intensity_minus_bg_smoothed = NULL)]
    training_df_by_mn_low[, `:=`(compound = NULL, sample = NULL, formula = NULL, ratio = NULL, intensity_predicted = NULL, intensity_smoothed_predicted = NULL,  intensity_minus_bg = NULL, intensity_minus_bg_smoothed = NULL)]
    test_df_by_mn_high[, `:=`(compound = NULL, sample = NULL, formula = NULL, ratio = NULL, intensity_predicted = NULL, intensity_smoothed_predicted = NULL,  intensity_minus_bg = NULL, intensity_minus_bg_smoothed = NULL)]
    test_df_by_mn_low[, `:=`(compound = NULL, sample = NULL, formula = NULL, ratio = NULL, intensity_predicted = NULL, intensity_smoothed_predicted = NULL,  intensity_minus_bg = NULL, intensity_minus_bg_smoothed = NULL)]
    training_df_by_mn_high[, `:=`(rt_left = NULL, rt_peak = NULL, rt_right = NULL, slice = NULL, total_slices = NULL, slice_peak_height = NULL, mplus1_intensity_smoothed = NULL, highest_original_intensity = NULL, rt_width = NULL)]
    training_df_by_mn_low[, `:=`(rt_left = NULL, rt_peak = NULL, rt_right = NULL, slice = NULL, total_slices = NULL, slice_peak_height = NULL, mplus1_intensity_smoothed = NULL, highest_original_intensity = NULL, rt_width = NULL)]
    test_df_by_mn_high[, `:=`(rt_left = NULL, rt_peak = NULL, rt_right = NULL, slice = NULL, total_slices = NULL, slice_peak_height = NULL, mplus1_intensity_smoothed = NULL, highest_original_intensity = NULL, rt_width = NULL)]
    test_df_by_mn_low[, `:=`(rt_left = NULL, rt_peak = NULL, rt_right = NULL, slice = NULL, total_slices = NULL, slice_peak_height = NULL, mplus1_intensity_smoothed = NULL, highest_original_intensity = NULL, rt_width = NULL)]
    training_df_by_mn_high <- training_df_by_mn_high[, lapply(.SD, as.numeric)]
    training_df_by_mn_low <- training_df_by_mn_low[, lapply(.SD, as.numeric)]
    test_df_by_mn_high <- test_df_by_mn_high[, lapply(.SD, as.numeric)]
    test_df_by_mn_low <- test_df_by_mn_low[, lapply(.SD, as.numeric)]
    
    #########################################
    # For this study, we also remove the columns "slice_diff_to_m0peak" and "Si" as after preliminary inspectation, slice_diff have nothing to do with bias and no metabolites contain "Si"
    training_df_by_mn_high[, c("slice_diff_to_m0peak", "Si") := NULL]
    training_df_by_mn_low[, c("slice_diff_to_m0peak", "Si") := NULL]
    test_df_by_mn_high[, c("slice_diff_to_m0peak", "Si") := NULL]
    test_df_by_mn_low[, c("slice_diff_to_m0peak", "Si") := NULL]
    ########################################
    
    # Convert data.tables to matrix for latering inputing into training
    Matrix_training_high <- as.matrix(training_df_by_mn_high)
    Matrix_training_low <- as.matrix(training_df_by_mn_low)
    Matrix_test_high <- as.matrix(test_df_by_mn_high)
    Matrix_test_low <- as.matrix(test_df_by_mn_low)
    # Extract intensity_predicted column from training_high_info
    DataInt_training_high <- training_high_info$intensity_minus_bg_smoothed / training_high_info$intensity_smoothed_predicted
    DataInt_test_high <- test_high_info$intensity_minus_bg_smoothed / test_high_info$intensity_smoothed_predicted
    
    # Inspect data types before training
    cat("Inspecting data types before training...\n")
    str(Matrix_training_high)
    

    # Set up cross-validation control
    rfControl <- trainControl(method = "cv", number = inner_cv, verboseIter = TRUE, summaryFunction = mse_summary)
    # Initial training on high intensity data
    options(expressions = 100000)  # Set to a higher value than default to allow for more iterations
    cat("Starting initial training on high-intensity data...\n")
    set.seed(seed_i)
    clf_rf_high <- train(x = Matrix_training_high,
                         y = DataInt_training_high,
                         method = "rf",
                         metric = "MSE",
                         trControl = rfControl,
                         ntree = ntree,
                         tuneGrid = data.frame(mtry = mtry_fixed))
    cat("i0 training completed.\n")
    # Save the initial model
    saveRDS(clf_rf_high, file = paste0("clf_rf_high_i0_seed", seed_i, ".rds"))
    # Extract and print importance
    importance_high <- as.data.frame(importance(clf_rf_high$finalModel))
    print(importance_high)

    
    # Test initial high intensity model accuracy
    predictions_high_training <- predict(clf_rf_high, newdata = Matrix_training_high)
    predictions_high_test <- predict(clf_rf_high, newdata = Matrix_test_high)
    
    cat("Training Performance (high-intensity data) of the iteration 0 model trained only by high-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    p_high_training_Int <- training_high_info$intensity_minus_bg_smoothed/predictions_high_training
    p_high_test_Int <- test_high_info$intensity_minus_bg_smoothed/predictions_high_test
    
    # Identify and print "super off" values (e.g., APE > 100%)
    training_ape <- abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    # Identify and print "super off" values (e.g., APE > 100%)
    test_ape <- abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")

    # For low intensity data below the cutoff, the initial clf_rf_high modell is used to predict their true intensity
    result_training <- process_data(Matrix_training_low, training_low_info, clf_rf_high)
    Matrix_training_low_m0 <- result_training$Matrix_m0
    Matrix_training_low_m1 <- result_training$Matrix_m1
    Matrix_training_low_m2 <- result_training$Matrix_m2
    training_low_info_m0 <- result_training$info_m0
    training_low_info_m1 <- result_training$info_m1
    training_low_info_m2 <- result_training$info_m2
    training_low_info_pseudoL <- result_training$info_pseudoL
    result_test <- process_data(Matrix_test_low, test_low_info, clf_rf_high)
    Matrix_test_low_m0 <- result_test$Matrix_m0
    Matrix_test_low_m1 <- result_test$Matrix_m1
    Matrix_test_low_m2 <- result_test$Matrix_m2
    test_low_info_m0 <- result_test$info_m0
    test_low_info_m1 <- result_test$info_m1
    test_low_info_m2 <- result_test$info_m2
    test_low_info_pseudoL <- result_test$info_pseudoL
    
    # Testing clf_rf_high on low data
    predictions_low_training <- predict(clf_rf_high, newdata = Matrix_training_low)
    predictions_low_test <- predict(clf_rf_high, newdata = Matrix_test_low)
    cat("Training Performance (low-intensity data) of the iteration 0 model trained only by high-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    p_low_training_Int <- training_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_training
    p_low_test_Int <- test_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_test
    
    training_ape <- abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    
    test_ape <- abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    # Then, the pesudo labled low intenisty data will be added to the high-intensity data for retraining the model
    Matrix_training_all <- rbind(Matrix_training_high, Matrix_training_low)
    training_all_info <- rbind(training_high_info, training_low_info_pseudoL)
    DataInt_training_all <- training_all_info$intensity_minus_bg_smoothed / training_all_info$intensity_smoothed_predicted
    Matrix_training_all <- as.matrix(Matrix_training_all)
    Matrix_test_all <- rbind(Matrix_test_high, Matrix_test_low)
    test_all_info <- rbind(test_high_info, test_low_info_pseudoL)
    DataInt_test_all <- test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted
    Matrix_test_all <- as.matrix(Matrix_test_all)
    
    predictions_all_training <- predict(clf_rf_high, newdata = Matrix_training_all)
    predictions_all_test <- predict(clf_rf_high, newdata = Matrix_test_all)
    cat("Training Performance (high and low intensity data pooled) of the iteration 0 model trained only by high-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    training_ape <- abs(100 * ((training_all_info$intensity_minus_bg_smoothed / predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted)
    
    # Identify and print "super off" values (e.g., APE > 100%)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    test_ape <- abs(100 * ((test_all_info$intensity_minus_bg_smoothed / predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    rm(clf_rf_high) # remove clf_rf_high to save RAM, if needed can use readRDS to load back manually
    
    # Inspect data types before training
    cat("Inspecting data types before training...\n")
    str(Matrix_training_all)
    
    # Set up cross-validation control
    rfControl <- trainControl(method = "cv", number = inner_cv, verboseIter = TRUE, summaryFunction = mse_summary)
    # Initial training on high intensity data
    options(expressions = 100000)  # Set to a higher value than default to allow for more iterations
    cat("Starting initial training on all-intensity data...\n")
    set.seed(seed_i)
    clf_rf_all <- train(x = Matrix_training_all, 
                         y = DataInt_training_all, 
                         method = "rf", 
                         metric = "MSE", 
                         trControl = rfControl, 
                         ntree = ntree, 
                         tuneGrid = data.frame(mtry = mtry_fixed))
    cat("i1 training completed.\n")
    # Save the initial model
    saveRDS(clf_rf_all, file = paste0("clf_rf_all_i1_seed", seed_i, ".rds"))
    # Extract and print importance
    importance_all <- as.data.frame(importance(clf_rf_all$finalModel))
    print(importance_all)
    
    # Test iteration 1 all intensity model accuracy
    predictions_high_training <- predict(clf_rf_all, newdata = Matrix_training_high)
    predictions_high_test <- predict(clf_rf_all, newdata = Matrix_test_high)
    
    cat("Training Performance (high-intensity data) of the iteration 1 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    p_high_training_Int <- training_high_info$intensity_minus_bg_smoothed/predictions_high_training
    p_high_test_Int <- test_high_info$intensity_minus_bg_smoothed/predictions_high_test
    
    # Identify and print "super off" values (e.g., APE > 100%)
    training_ape <- abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")

    # Identify and print "super off" values (e.g., APE > 100%)
    test_ape <- abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    # For low intensity data below the cutoff, the initial clf_rf_all modell is used to predict their true intensity
    result_training <- process_data(Matrix_training_low, training_low_info, clf_rf_all)
    Matrix_training_low_m0 <- result_training$Matrix_m0
    Matrix_training_low_m1 <- result_training$Matrix_m1
    Matrix_training_low_m2 <- result_training$Matrix_m2
    training_low_info_m0 <- result_training$info_m0
    training_low_info_m1 <- result_training$info_m1
    training_low_info_m2 <- result_training$info_m2
    training_low_info_pseudoL <- result_training$info_pseudoL
    result_test <- process_data(Matrix_test_low, test_low_info, clf_rf_all)
    Matrix_test_low_m0 <- result_test$Matrix_m0
    Matrix_test_low_m1 <- result_test$Matrix_m1
    Matrix_test_low_m2 <- result_test$Matrix_m2
    test_low_info_m0 <- result_test$info_m0
    test_low_info_m1 <- result_test$info_m1
    test_low_info_m2 <- result_test$info_m2
    test_low_info_pseudoL <- result_test$info_pseudoL
    
    # Testing clf_rf_all on low data
    predictions_low_training <- predict(clf_rf_all, newdata = Matrix_training_low)
    predictions_low_test <- predict(clf_rf_all, newdata = Matrix_test_low)
    cat("Training Performance (low-intensity data) of the iteration 1 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    p_low_training_Int <- training_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_training
    p_low_test_Int <- test_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_test
    
    training_ape <- abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    test_ape <- abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    # Then, the pesudo labled low intenisty data will be added to the high-intensity data for retraining the model
    Matrix_training_all <- rbind(Matrix_training_high, Matrix_training_low) # Matrix_training_all already exists but I want to keep the copy and past the same in iterations
    training_all_info <- rbind(training_high_info, training_low_info_pseudoL)
    DataInt_training_all <- training_all_info$intensity_minus_bg_smoothed / training_all_info$intensity_smoothed_predicted
    Matrix_training_all <- as.matrix(Matrix_training_all)
    Matrix_test_all <- rbind(Matrix_test_high, Matrix_test_low)
    test_all_info <- rbind(test_high_info, test_low_info_pseudoL)
    DataInt_test_all <- test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted
    Matrix_test_all <- as.matrix(Matrix_test_all)
    
    predictions_all_training <- predict(clf_rf_all, newdata = Matrix_training_all)
    predictions_all_test <- predict(clf_rf_all, newdata = Matrix_test_all)
    cat("Training Performance (high and low intensity data pooled) of the iteration 1 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    training_ape <- abs(100 * ((training_all_info$intensity_minus_bg_smoothed / predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted)
    
    # Identify and print "super off" values (e.g., APE > 100%)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    test_ape <- abs(100 * ((test_all_info$intensity_minus_bg_smoothed / predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    rm(clf_rf_all)
    
    # Inspect data types before training
    cat("Inspecting data types before training...\n")
    str(Matrix_training_all)

    # Set up cross-validation control
    rfControl <- trainControl(method = "cv", number = inner_cv, verboseIter = TRUE, summaryFunction = mse_summary)
    # Initial training on high intensity data
    options(expressions = 100000)  # Set to a higher value than default to allow for more iterations
    cat("Starting iteration 2 training on all-intensity data...\n")
    set.seed(seed_i)
    clf_rf_all_i2 <- train(x = Matrix_training_all, 
                        y = DataInt_training_all, 
                        method = "rf", 
                        metric = "MSE", 
                        trControl = rfControl, 
                        ntree = ntree, 
                        tuneGrid = data.frame(mtry = mtry_fixed))
    cat("i2 training completed.\n")
    # Save the initial model
    saveRDS(clf_rf_all_i2, file = paste0("clf_rf_all_i2_seed", seed_i, "_2cRF.rds"))
    # Extract and print importance
    importance_all <- as.data.frame(importance(clf_rf_all_i2$finalModel))
    print(importance_all)

    
    # Test iteration 1 all intensity model accuracy
    predictions_high_training <- predict(clf_rf_all_i2, newdata = Matrix_training_high)
    predictions_high_test <- predict(clf_rf_all_i2, newdata = Matrix_test_high)
    
    cat("Training Performance (high-intensity data) of the iteration 1 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    p_high_training_Int <- training_high_info$intensity_minus_bg_smoothed/predictions_high_training
    p_high_test_Int <- test_high_info$intensity_minus_bg_smoothed/predictions_high_test
    
    # Identify and print "super off" values (e.g., APE > 100%)
    training_ape <- abs(100 * ((training_high_info$intensity_minus_bg_smoothed / predictions_high_training) - training_high_info$intensity_smoothed_predicted) / training_high_info$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    # Identify and print "super off" values (e.g., APE > 100%)
    test_ape <- abs(100 * ((test_high_info$intensity_minus_bg_smoothed / predictions_high_test) - test_high_info$intensity_smoothed_predicted) / test_high_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    # For low intensity data below the cutoff, the initial clf_rf_high modell is used to predict their true intensity
    result_training <- process_data(Matrix_training_low, training_low_info, clf_rf_all_i2)
    Matrix_training_low_m0 <- result_training$Matrix_m0
    Matrix_training_low_m1 <- result_training$Matrix_m1
    Matrix_training_low_m2 <- result_training$Matrix_m2
    training_low_info_m0 <- result_training$info_m0
    training_low_info_m1 <- result_training$info_m1
    training_low_info_m2 <- result_training$info_m2
    training_low_info_pseudoL <- result_training$info_pseudoL
    result_test <- process_data(Matrix_test_low, test_low_info, clf_rf_all_i2)
    Matrix_test_low_m0 <- result_test$Matrix_m0
    Matrix_test_low_m1 <- result_test$Matrix_m1
    Matrix_test_low_m2 <- result_test$Matrix_m2
    test_low_info_m0 <- result_test$info_m0
    test_low_info_m1 <- result_test$info_m1
    test_low_info_m2 <- result_test$info_m2
    test_low_info_pseudoL <- result_test$info_pseudoL
    
    # Testing clf_rf_all_i2 on low data
    predictions_low_training <- predict(clf_rf_all_i2, newdata = Matrix_training_low)
    predictions_low_test <- predict(clf_rf_all_i2, newdata = Matrix_test_low)
    cat("Training Performance (low-intensity data) of the iteration 2 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    p_low_training_Int <- training_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_training
    p_low_test_Int <- test_low_info_pseudoL$intensity_minus_bg_smoothed/predictions_low_test
    
    training_ape <- abs(100 * ((training_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_training) - training_low_info_pseudoL$intensity_smoothed_predicted) / training_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    test_ape <- abs(100 * ((test_low_info_pseudoL$intensity_minus_bg_smoothed / predictions_low_test) - test_low_info_pseudoL$intensity_smoothed_predicted) / test_low_info_pseudoL$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    # Then, the pesudo labled low intenisty data will be added to the high-intensity data for retraining the model
    Matrix_training_all <- rbind(Matrix_training_high, Matrix_training_low)
    training_all_info <- rbind(training_high_info, training_low_info_pseudoL)
    DataInt_training_all <- training_all_info$intensity_minus_bg_smoothed / training_all_info$intensity_smoothed_predicted
    Matrix_training_all <- as.matrix(Matrix_training_all)
    Matrix_test_all <- rbind(Matrix_test_high, Matrix_test_low)
    test_all_info <- rbind(test_high_info, test_low_info_pseudoL)
    DataInt_test_all <- test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted
    Matrix_test_all <- as.matrix(Matrix_test_all)
    
    predictions_all_training <- predict(clf_rf_all_i2, newdata = Matrix_training_all)
    predictions_all_test <- predict(clf_rf_all_i2, newdata = Matrix_test_all)
    cat("Training Performance (high and low intensity data pooled) of the iteration 2 model trained by both high-intensity and pseudo-labeled low-intensity data:\n")
    cat("Training Mean APE Error:", mean(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Training Median APE Error:", median(abs(100 * ((training_all_info$intensity_minus_bg_smoothed/predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Mean APE Error:", mean(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    cat("Test Median APE Error:", median(abs(100 * ((test_all_info$intensity_minus_bg_smoothed/predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted), na.rm = TRUE), "\n")
    
    training_ape <- abs(100 * ((training_all_info$intensity_minus_bg_smoothed / predictions_all_training) - training_all_info$intensity_smoothed_predicted) / training_all_info$intensity_smoothed_predicted)
    
    # Identify and print "super off" values (e.g., APE > 100%)
    super_off_training <- which(training_ape > 50)
    cat("Super Off Training Values (APE > 50%):\n")
    cat(length(super_off_training), "out of", length(training_ape), "training Mn Signals are > 50% off.\n")
    test_ape <- abs(100 * ((test_all_info$intensity_minus_bg_smoothed / predictions_all_test) - test_all_info$intensity_smoothed_predicted) / test_all_info$intensity_smoothed_predicted)
    super_off_test <- which(test_ape > 50)
    cat("Super Off Test Values (APE > 50%):\n")
    cat(length(super_off_test), "out of", length(test_ape), "test Mn Signals are > 50% off.\n")
    
    ###############################################################
    # Next, we will calculate permutation importance of clf_rf_all_i2
    # Be careful that if a feature has both unsmoothed and smoothed version, they will be permutated together in group
    # Define the groups of features that need to permute together (smoothed and unsmoothed)
    cat("Start Calculating Permutation Importance:\n")
    feature_groups <- list(
      "intensity_over_TIC_group" = c("intensity_over_TIC_smoothed", "intensity_over_TIC"),
      "mn_intensity_group" = c("mn_intensity", "mn_intensity_smoothed"),
      "total_ion_group" = c("total_ion", "total_ion_smoothed")
    )
    grouped_features <- unlist(feature_groups)
    independent_features <- setdiff(colnames(Matrix_test_all), grouped_features)
    
    # Step 4: Calculate individual permutation importance
    set.seed(seed_i)
    importance_scores_ind <- numeric(length(independent_features))
    names(importance_scores_ind) <- independent_features
    
    original_predictions <- predict(clf_rf_all_i2$finalModel, Matrix_test_all) # only use test data to calculate permutation score
    original_mse <- mean((original_predictions - test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted)^2)
    
    for (feature in independent_features) {
      cat("Processing independent feature:", feature, "\n")
      permuted_data <- Matrix_test_all
      permuted_data[, feature] <- sample(permuted_data[, feature])
      permuted_predictions <- predict(clf_rf_all_i2$finalModel, permuted_data)
      permuted_mse <- mean((permuted_predictions - test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted)^2)
      importance_scores_ind[feature] <- permuted_mse - original_mse
    }
    
    importance_ind_df <- data.frame(
      Name = names(importance_scores_ind),
      PermutationImportance = importance_scores_ind,
      Type = "Feature"
    )
    
    importance_scores_grouped <- numeric(length(feature_groups))
    names(importance_scores_grouped) <- names(feature_groups)
    for (group_name in names(feature_groups)) {
      cat("Processing feature group:", group_name, "\n")
      permuted_data <- Matrix_test_all
      group_indices <- which(colnames(permuted_data) %in% feature_groups[[group_name]])
      permuted_rows <- sample(nrow(permuted_data))
      permuted_data[, group_indices] <- permuted_data[permuted_rows, group_indices]
      permuted_predictions <- predict(clf_rf_all_i2$finalModel, permuted_data)
      permuted_mse <- mean((permuted_predictions - test_all_info$intensity_minus_bg_smoothed / test_all_info$intensity_smoothed_predicted)^2)
      importance_scores_grouped[group_name] <- permuted_mse - original_mse
    }
    
    importance_grouped_df <- data.frame(
      Name = names(importance_scores_grouped),
      PermutationImportance = importance_scores_grouped,
      Type = "Group"
    )

    # Save importance data as xlsx in aphabatical order
    importance_combined_df <- rbind(importance_ind_df, importance_grouped_df)
    sheet1 <- importance_combined_df
    sheet2 <- importance_combined_df
    total_importance <- sum(sheet2$PermutationImportance)
    sheet2$PermutationImportance <- round(100 * sheet2$PermutationImportance / total_importance, 2)
    colnames(sheet2)[which(colnames(sheet2) == "PermutationImportance")] <- "ImportanceWeightPercent"
    importance_sheets <- list(
      RawImportance = sheet1,
      RelativeWeight = sheet2
    )
    write_xlsx(importance_sheets, path = paste0("clf_rf_all_i2_seed", seed_i, "_2cRF_importance_rank_fig3.xlsx"))
    
    # Plot importance
    ggplot(importance_combined_df, aes(x = reorder(Name, PermutationImportance), y = PermutationImportance, fill = Type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = "Combined Permutation Importance", x = "Feature or Group", y = "Importance (Increase in MSE)") +
      theme_minimal() +
      scale_fill_manual(values = c("Feature" = "steelblue", "Group" = "tomato"))
    
    #############################################################
    # After the final iteration 2 model is trained, we will next apply this final model to each slice of each metabolites
    # to generate single data points to plot figure 2, 3, 4 ...
    # Convert the data frames to data.tables (if they are not already)
    set.seed(seed_i)
    df_by_mn <- rbind(df_by_mn_high, df_by_mn_low)
    training_df_by_slice <- df_by_slice %>% filter(!sample %in% selected_values)
    test_df_by_slice <- df_by_slice %>% filter(sample %in% selected_values)
    training_df_by_mn <- df_by_mn %>% filter(!sample %in% selected_values)
    test_df_by_mn <- df_by_mn %>% filter(sample %in% selected_values)
    
    # Apply the function to the training data
    process_and_export_corrected_df(
      df_by_mn_input = training_df_by_mn,
      df_by_slice_input = training_df_by_slice,
      df_input = training_df,
      clf_model = clf_rf_all_i2,
      seed_i = seed_i,
      output_prefix = "Corrected_training"
    )
    
    # Apply the function to the test data
    process_and_export_corrected_df(
      df_by_mn_input = test_df_by_mn,
      df_by_slice_input = test_df_by_slice,
      df_input = test_df,
      clf_model = clf_rf_all_i2,
      seed_i = seed_i,
      output_prefix = "Corrected_test"
    )

    ##############################################################
    
    rm(clf_rf_all_i2)
    seed_i <- seed_i + 1 # plus 1 to seed and start next rerun trials
  }

sink() # end logging 



