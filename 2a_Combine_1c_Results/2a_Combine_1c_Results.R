###################################

# File: 2a_Combine_1c_Results
# Input: "Sample_XXX_1c.xlsx" of all experiments should be moved to the same folder of this 2a codes.
# Function: Step 2a codes will combine xlsx files of all experiments into one 2a xlsx file.
# Output: "All_1c_Combined_2a.xlsx" have two combined sheets of all inputted 1c: sheet1 contains the integrated information, and sheet2 contains data of each individual slices
# Usage: Put this 2a codes and all 1c results files in one folder, run the codes and be careful of if any bugs happened.
# Next to Do: Quick eyeballing on the "All_1c_Combined_2a.xlsx" to see if it contains the correct data, next move this
#             output to the "2b_Preprocess_Input_Data" folder to preprocess the data before bias correction model training.

####################################
#####################################

# Load necessary libraries
library(readxl) # handle multiple sheets xlsx files
library(openxlsx) # handle xlsx in general
library(dplyr) # to bind df

# Get and check the current working directory
directory <- getwd()
cat("Current working directory:\n")
print(directory)
# List all 1c XLSX files in the directory
files <- list.files(directory, pattern = "\\_1c.xlsx$", full.names = TRUE)
cat("1c xlsx result files found:\n")
print(files)
# Check if files were found
if (length(files) == 0) {
  stop("No 1c xlsx result files found in the current directory.")
}

# Determine the number of sheets in the first file
first_file_sheets <- excel_sheets(files[1])
num_sheets <- length(first_file_sheets)
# Print the sheet names found in the first file
cat("Sheet names in the first file:\n")
print(first_file_sheets)

# Create a workbook to save combined sheets
wb <- createWorkbook()
# Function combine_sheets to read and combine all sheets of different 1c xlsx
combine_sheets <- function(sheet_name, files) {
  combined_data <- NULL
  
  for (file in files) {
    if (file.exists(file)) {
      actual_sheets <- excel_sheets(file)
      cat("Available sheets in", file, ":\n")
      print(actual_sheets)
      
      if (sheet_name %in% actual_sheets) {
        sheet_data <- read_excel(file, sheet = sheet_name)
        # Only for the first xlsx, keep the col names
        if (is.null(combined_data)) {
          combined_data <- sheet_data
        } else {
          combined_data <- bind_rows(combined_data, sheet_data)
        }
      } else {
        warning(paste("Sheet", sheet_name, "not found in file:", file))
      }
    } else {
      warning(paste("File does not exist:", file))
    }
  }
  return(combined_data)
}

# Start to Combine data for each sheet and add to the workbook
for (sheet_num in 1:num_sheets) {
  sheet_name <- first_file_sheets[sheet_num]
  cat("Processing sheet:", sheet_name, "\n") # Show Current Progress
  combined_data <- combine_sheets(sheet_name, files)
  # Additional check if combined_data is not empty before adding to workbook
  if (!is.null(combined_data)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, combined_data)
  } else {
    warning(paste("No data to combine for sheet:", sheet_name))
  }
}

# Save the combined workbook
saveWorkbook(wb, file = "All_1c_Combined_2a.xlsx", overwrite = TRUE)
cat("Combined workbook saved as 'All_1c_Combined_2a.xlsx'.")



