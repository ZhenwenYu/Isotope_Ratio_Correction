###################################

# File: 1a_Filter_Alignment_PeakID_Info.R
# Input: "PeakID_DateXXX.txt" is the alignment file exported from MS-DIAL ver 4.8,
#         which contains the PeakID information indicating which peaks are the metabolites of interest.
#         If MS-DIAL assigned the wrong peak during group alignment, we can manually correct 
#         this text file using excel before or after importing.
# Function: Corerct the alignment file format and keep only only identified metabolites in the alignment
# Output: "PeakID_DateXXX_1a.txt" the alignment file that kept only identified metabolites in the correct file format 

# Usage: Shall change the alignment file name in line 21
# Next to Do: Manually check the output "PeakID_DateXXX_1a.txt" with excel and do any necessary modifications, 
#             we will use this file to do "1b_Extract_Metabolites_Identification_Info.R"

####################################
#####################################

#IMPORTANT: change the alignment file name in line 23

# install.packages("openxlsx")
library(openxlsx)  # For reading xlsx files

xlsx_name <- "PeakID_0_20246231817"

######################################
#####################################

# Set working directory to be the current folder
current_path <- getwd()

xlsx_path <- paste0(xlsx_name, ".txt")
xlsx_1a_name <- paste0(xlsx_name, "_1a.txt")

# Load the xlsx/txt file
xlsx_data <- read.csv(xlsx_path, header = FALSE, sep = "\t")

# Set row 5 data as column name and remove the first 4 empty rows. Be Cautious that MS-DIAL versions other than
# ver 4.8 may have a different export file format and requires modification of these codes.
colnames(xlsx_data) <- as.character(xlsx_data[5, ])
xlsx_data <- xlsx_data[-c(1:5), ]

# Remove rows where "Metabolite name" contains "Unknown" or "w/o MS2:" so the alignment will only keeps "A:Confidence"
# and "B:Unsettled" marked metabolites
xlsx_data <- xlsx_data[!grepl("Unknown|w/o MS2:", xlsx_data[["Metabolite name"]]), ]

# Save the new alignment file
write.table(xlsx_data, file = xlsx_1a_name, sep = "\t", row.names = FALSE)

# Next to Do: After saving the 1a file, it is always a good idea to check the filtered alignment file manually to
#             correct any wrong alignment annotation given by MS-DIAL and manually modify the "B:Unsettled" metabolites
