###################################

# File: 1b_Extract_Metabolites_Identification_Info.R
# Input: "PeakID_DateXXX_1a.txt" is the alignment file exported from MS-DIAL ver 4.8 and has undergo step 1a processing,
#         This "1a" alignment file should have been filtered through codes "1a", and been manually corrected
#         if necessary to correct wrong alignment information given by MS-DIAL
#         "Sample_XXX.txt"  is the sample peak information file exported directly from MS-DIAL ver 4.8
# Function: Based on the correct "PeakID_DateXXX_1a.txt" alignment file, this codes filter and extract metabolites peak
#           identification info from all sample files.
# Output: "Sample_XXX_1b.txt" for each sample that contains the correwct identifciation information such as formula, RT, m/z.

# Usage: Shall change input and output file names, mzml file m/z range as needed. Could change m/z ppm window.
# Next to Do: Manually check the output "Sample_XXX_1b.txt" with excel and do any necessary modifications.
#             For example, if MS-DIAL identify one metabolite peak to be two separated peak, we can manually
#             combine the two peaks row data into one row by updating the left, peak height, and right retention time info
#             we will use this "Sample_XXX_1b.txt" file to do "1c_Extract_Metabolites_M0M1M2_using_Identification.R"

####################################
######################################

#IMPORTANT: change the alignment file name in line 26 and samples name from line 28 to 33

# install.packages("openxlsx")
library(openxlsx)  # For reading xlsx files

xlsx_name <- "PeakID_0_202431195"

Peak_list <- c(
  "std_e1",
  "std_e2",
  "std_e3",
  "std_e4"
)

# Sometimes, the MS-DIAL alignment file identify a metabolite with confidence in most samples 
# but not in the specific sample this iteration processes. In this condition, we will use "MS1_mz_tolerance"
# to do an additional check to confirm this samples' corresponding peak have similar MS1 m/z with the alignment one.
# "0.005" should be a robust value for Orbitrap
MS1_mz_tolerance <- 0.005

######################################
########################################

# Set working directory to be the current folder
current_path <- getwd()

# Load the xlsx/txt file
xlsx_1a_name <- paste0(xlsx_name, "_1a.txt")
xlsx_data <- read.csv(xlsx_1a_name, header = TRUE, sep = "\t")

for (peak_list in Peak_list){
  # Find the alignment file columns that match the samples and get the peak ID
  ID_list <- xlsx_data[[peak_list]]
  
  # Find sample files of the same folder
  txt_file_name <- paste0(peak_list, ".txt")
  input_file <- read.csv(txt_file_name, header = TRUE, sep = "\t")
  
  # Filter out sample entries with IDs in alignment peak ID
  txt_data <- input_file[input_file$PeakID %in% ID_list, ]
  unknown_data <- txt_data[txt_data$Title == "Unknown", ] # Some data of the sample may not be identified with confidence, extract them to do one more check later
  txt_data <- txt_data[txt_data$Title != "Unknown", ]
  
  # Iterate through each data row of the txt data
  if (nrow(txt_data) > 0) {  
    for (i in 1:nrow(txt_data)) {
      peak_id <- txt_data[i, "PeakID"]
      
      # Find the row in xlsx_data with matching Peak ID
      matching_row <- xlsx_data[xlsx_data[[peak_list]] == peak_id, ]
      
      if (nrow(matching_row) > 0) {  # If a match is found
        new_name <- matching_row$`Metabolite.name`
        txt_data[i, "Title"] <- new_name # Use Alignment metabolite name to replace sample metabolite name for consistency
        new_formula <- matching_row$`Formula`
        txt_data[i, "Formula"] <- new_formula
      }
    }
  }  
  if (nrow(unknown_data) > 0) { 
    for (n in 1:nrow(unknown_data)){
      peak_id <- unknown_data[n, "PeakID"]
      matching_row <- xlsx_data[xlsx_data[[peak_list]] == peak_id, ]
      
      # Additional check to handle the condition when the alignment file identify a metabolite with confidence in most samples but not the specific sample this iteration processes
      if (nrow(matching_row) > 0) { 
        if (abs(unknown_data$Precursor.m.z[n] - matching_row$Reference.m.z) <= MS1_mz_tolerance){
          new_name <- matching_row$`Metabolite.name`
          unknown_data[n, "Title"] <- new_name
          new_formula <- matching_row$`Formula`
          unknown_data[n, "Formula"] <- new_formula
          txt_data <- rbind(txt_data, unknown_data[n,])
        }
      }
    }
  }
  
  # Save the updated txt data to a new txt file
  txt_file_1a_name <- paste0(peak_list, "_1b.txt")
  write.table(txt_data, file = txt_file_1a_name, sep = "\t", row.names = FALSE)
}

# Next to Do: Manually check the output "Sample_XXX_1b.txt" with excel and do any necessary modifications.
#             For example, if MS-DIAL identify one metabolite peak to be two separated peak, we can manually
#             combine the two peaks row data into one row by updating the left, peak height, and right retention time info
#             we will use this "Sample_XXX_1b.txt" file to do "1c_Extract_Metabolites_M0M1M2_using_Identification.R"