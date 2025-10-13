###################################

# File: 1c_Extract_Metabolites_M0M1M2_using_Identification
# Input: "Sample_XXX_1b.txt" is the output of step 1b, which contains the peak identification information of 
#         each metabolite of interest of each sample.
#         "Sample_XXX.mzML" is the open-source mzML format of mass spectrometery experients result files, which could be easily
#         converted from each company's result files through MSConvert (part of ProteoWizard)
# Function: Extract the M0, M1, and M2 information plus any information that may contribute to measurement bias
#         from the mzML result files based on the metabolite peak identification information given in the 1b .txt files
# Output: "Sample_XXX_1c.xlsx" have two sheets: sheet1 contains the integrated information, and sheet2 contains data of each individual slices

# Usage: For first time usage, please install all library properly. Must tell the program
#         the sample name and mzML result sample names (if different), and the method mz range
# Next to Do: Manually check the output "Sample_XXX_1c.xlsx" with excel and see if the calculated M0, M1 and M2 ratio is close
#             to the theory predictions. If not, this indicates the Orbitrap measurement of MID indeed has a bias and would need step 2
#             coding to train and apply the bias correction models.

####################################
#####################################

# IMPORTANT: for 1st time usage, please intall the libarary properly
# install.packages("BiocManager")
# BiocManager::install("MSnbase", force = TRUE)
# install.packages("openxlsx")
# install.packages("reticulate")
# install.packages("signal")
# install.packages("BiocParallel")

library(MSnbase) # To process mzml files
library(openxlsx) # For reading xlsx files
library(reticulate) # For using python files "mida_alg.py" and "numeric_functions.py"
library(signal)  # For Savitzky-Golay smoothing
library(BiocParallel)
register(SerialParam()) # library(BiocParallel) restarts all parallel workers in case some used 
                        # by library(MSnbase) are not been properly turned off in force stop


current_path <- getwd()
source("1c_Functions.R") # Seperate files to contain all additional functions, some functions 
                         # require paramaters given in this 1c file thus may not be used independently 

# IMPORTANT: Python needs to be already installed in the computer separately, and then run the following 
# three lines with library(reticulate) loaded at the 1st time. Answer "No" for creating default reticulate enviorment
# np <- import("numpy")
# py_install("scipy")
# py_install("pandas")

source_python("numeric_functions.py") 
source_python("mida_alg.py") # Two python files from https://github.com/naveedziari/PyMIDA/tree/main/pymida 
                            # but with one line commented out to do Mass Isotopomer Distribution Analysis in R

# IMPORTANT: change the peak and mzML file name in following lines accordingly
Peak_list <- c(
  "std_e1",
  "std_e2",
  "std_e3",
  "std_e4"
)
mzML_list <- c(
  "std_e1",
  "std_e2",
  "std_e3",
  "std_e4"
)

# IMPORTANT: Tell codes the mass window of the mass spec method in Dalton
mz_total_region <- c(69,360)

####################################################

# The following 4 parameters are not recommenced be modified, but some mass spec/methods may need a different value
ppm = 10                # The m/z ppm tolerance window to extract M0M1M2, some mass spec may have poorer resolution thus needs a higher ppm tolerance
slice_int_cutoff = 10000 # will be used in function pick_edges, which is the minimum Mn intensity a slice need to have to let that slice be
                        # included to do M0M1M2 calculations. Empirically, if a slice's M1 or M2 have very low intensity, that slice's M1 or M2 
S_N_BL_ratio_cutoff = 3 # set the Signal/Noise and Signal/Baseline ratio cutoff. 3 is field-recognized default
                        # Baseline: background signal of the area. For example, peak baseline is largely based on 
                        # the intensity when the peak starts and ends 
                        # Noise: calculated using the RMS method by checking the variation of the background baseline before the metabolite peak
correlation_cutoff = 0.800 # the correlation score cutoff when comparing M0, M1 and M2 peaks to do additional check on if the M1 or M2 peak are been
                          # disturbed by other peaks. Either Pearson score or Spearman score pass the cutoff threshold is sufficient.

# Physical constants, by default should not be modified
massdiffH1H2=2.01410177812-1.00782503223
massdiffC12C13=13.00335483507-12
massdiffN14N15=15.00010889888-14.00307400443
massdiffO16O17=16.99913175650-15.99491461957
massdiffO16O18=17.99915961286-15.99491461957
massdiffS32S33=32.9714587098-31.9720711744
massdiffS32S34=33.967867004-31.9720711744

################################################################
#################################################################

# Start c_Extract_Metabolites_M0M1M2_using_Identification
for (peak_list in Peak_list)
{
  
  # Specify the file name, need to be cautious the name format are correct here
  anno_name <- peak_list
  annotation_name <- paste0(anno_name, "_1b.txt")
  plposition <- match(peak_list, Peak_list)
  sample_name <- mzML_list[plposition]
  MS_name <- paste0(sample_name, ".mzML")
  output_name <- paste0(anno_name, "_1c.xlsx")
  
  # Create the file path by combining the current path and file name
  annotation_path <- file.path(current_path, annotation_name)
  MS_path <- file.path(current_path, MS_name)
  
  # Open the file connection
  annotation <- read.csv(annotation_path, header = TRUE, sep = "\t")
  MS1file <- readMSData(MS_path, msLevel = 1)
  
  # Initialize two new empty tables to store the extracted information: integrated area and individual slices
  extracted_table <- data.frame()
  extracted_table2 <- data.frame()
  
  # Loop through each row in the compound table
  for (i in 1:nrow(annotation)) 
    tryCatch(
      {
        # Initialize one new empty table to temporally carry the extracted individual slices info
        extracted_slices <- data.frame()
        
        # Extract the compound and its corresponding values from "1b" files
        compound <- annotation$Title[i]
        heightnumber <- annotation$Height[i]
        reference_mz <- annotation$Reference.m.z[i]
        mz_M0 <- annotation$`Precursor.m.z`[i]
        formula <- annotation$Formula[i]
        elementscount <- countElements(formula)
        warningmessage <- ""
        # Additional check on if the observed mz is similar to reference libarary
        if (abs(reference_mz - mz_M0) > 0.003) {
          warningmessage <- paste(warningmessage, "mz_M0 may be far from reference_mz")
        }
        
        # Get the left and right retention time values and create the rt_range with left and right values
        left_value <- annotation$'RT.left.min.'[i] * 60
        right_value <- annotation$'RT.right..min.'[i] * 60
        middle_value <-annotation$'RT..min.'[i] * 60
        rt_range <- c(max(0, left_value - 5),  min(900, right_value + 5))
        rt_range_180s <- c(max(0, left_value - 180), min(900, right_value + 180))
      
        # Calculate the m/z values for M0, M1, and M2 using the physical constants stated at the beginning, 
        # using H1H2 and N14N15 mass diff should cover all the isotopologues variants.
        mz_M1_low <- mz_M0 + massdiffN14N15
        mz_M1_high <- mz_M0 + massdiffH1H2
        mz_M2_low <- mz_M0 + massdiffN14N15*2
        mz_M2_high <- mz_M0 + massdiffH1H2*2
        mz_M3_low <- mz_M0 + massdiffN14N15*3
        mz_M3_high <- mz_M0 + massdiffH1H2*3
        mz_M4_low <- mz_M0 + massdiffN14N15*4
        mz_Mminus1_high <- mz_M0 - massdiffH1H2
        
        # Specify the m/z regions
        mz_M0_region <- c(mz_M0*(1-(ppm)/(10^6)), mz_M0*(1+(ppm)/(10^6)))
        mz_M1_region <- c(mz_M1_low*(1-(ppm)/(10^6)), mz_M1_high*(1+(ppm)/(10^6)))
        mz_M2_region <- c(mz_M2_low*(1-(ppm)/(10^6)), mz_M2_high*(1+(ppm)/(10^6)))
        mz_M3_region <- c(mz_M3_low*(1-(ppm)/(10^6)), mz_M3_high*(1+(ppm)/(10^6)))
        mz_Mminus1_to_M0_region <- c(mz_Mminus1_high, mz_M0*(1-(ppm)/(10^6)))  
        mz_M0_to_M1_region <- c(mz_M0*(1+(ppm)/(10^6)), mz_M1_low*(1-(ppm)/(10^6)))   
        mz_M1_to_M2_region <- c(mz_M1_high*(1+(ppm)/(10^6)), mz_M2_low*(1-(ppm)/(10^6)))   
        mz_M2_to_M3_region <- c(mz_M2_high*(1+(ppm)/(10^6)), mz_M3_low*(1-(ppm)/(10^6)))   
        mz_M3_to_M4_region <- c(mz_M3_high*(1+(ppm)/(10^6)), mz_M4_low*(1-(ppm)/(10^6))) 
        
        # Extract the chromatographic peaks within the RT range of the peak, and also within the +- 180s range
        rt_filtered_ms1 <- filterRt(MS1file, rt = rt_range)
        rt_180s_filtered_ms1 <- filterRt(MS1file, rt = rt_range_180s)
        
        suppressWarnings(
        M0_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M0_region)
        )
        peaksM0 <- chromatogram(M0_mz_filtered_ms1, rt = rt_range, mz = mz_M0_region)
        suppressWarnings(
        M1_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M1_region)
        )
        peaksM1 <- chromatogram(M1_mz_filtered_ms1, rt = rt_range, mz = mz_M1_region)
        suppressWarnings(
        M2_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M2_region)
        )
        peaksM2 <- chromatogram(M2_mz_filtered_ms1, rt = rt_range, mz = mz_M2_region)
        
        suppressWarnings(
          M0_180s_mz_filtered_ms1 <- filterMz(rt_180s_filtered_ms1, mz = mz_M0_region)
        )
        peaksM0_180s <- chromatogram(M0_180s_mz_filtered_ms1, rt = rt_range_180s, mz = mz_M0_region)
        suppressWarnings(
          M1_180s_mz_filtered_ms1 <- filterMz(rt_180s_filtered_ms1, mz = mz_M1_region)
        )
        peaksM1_180s <- chromatogram(M1_180s_mz_filtered_ms1, rt = rt_range_180s, mz = mz_M1_region)
        suppressWarnings(
          M2_180s_mz_filtered_ms1 <- filterMz(rt_180s_filtered_ms1, mz = mz_M2_region)
        )
        peaksM2_180s <- chromatogram(M2_180s_mz_filtered_ms1, rt = rt_range_180s, mz = mz_M2_region)
        
        # Info of the M3 mz region, total ion mz range region, and intervals between each Mn are also stored for later research on if they contribute to the bias
        suppressWarnings(
        M3_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M3_region)
        )
        peaksM3 <- chromatogram(M3_mz_filtered_ms1, rt = rt_range, mz = mz_M3_region)
        suppressWarnings(
          M3_180s_mz_filtered_ms1 <- filterMz(rt_180s_filtered_ms1, mz = mz_M3_region)
        )
        peaksM3_180s <- chromatogram(M3_180s_mz_filtered_ms1, rt = rt_range_180s, mz = mz_M3_region)
        
        suppressWarnings(
        total_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_total_region)
        )
        peakstotal <- chromatogram(total_mz_filtered_ms1, rt = rt_range, mz = mz_total_region)
        suppressWarnings(
        Mminus1_to_M0_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_Mminus1_to_M0_region)
        )
        peaksMminus1_to_M0 <- chromatogram(Mminus1_to_M0_mz_filtered_ms1, rt = rt_range, mz = mz_Mminus1_to_M0_region)
        suppressWarnings(
        M0_to_M1_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M0_to_M1_region)
        )
        peaksM0_to_M1 <- chromatogram(M0_to_M1_mz_filtered_ms1, rt = rt_range, mz = mz_M0_to_M1_region)
        suppressWarnings(
        M1_to_M2_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M1_to_M2_region)
        )
        peaksM1_to_M2 <- chromatogram(M1_to_M2_mz_filtered_ms1, rt = rt_range, mz = mz_M1_to_M2_region)
        suppressWarnings(
        M2_to_M3_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M2_to_M3_region)  
        )
        peaksM2_to_M3 <- chromatogram(M2_to_M3_mz_filtered_ms1, rt = rt_range, mz = mz_M2_to_M3_region)
        suppressWarnings(
        M3_to_M4_mz_filtered_ms1 <- filterMz(rt_filtered_ms1, mz = mz_M2_to_M3_region)  
        )
        peaksM3_to_M4 <- chromatogram(MS1file, rt = rt_range, mz = mz_M3_to_M4_region) 
      
         
        # Extract the height and area
        intensitiesM0 <- as.vector(intensity(peaksM0[1]))
        intensitiesM0[is.na(intensitiesM0)] <- 0
        total_slices <- length(intensitiesM0)
        slices_rt <- as.vector(rtime(peaksM0[1]))
        
        intensitiesM1 <- as.vector(intensity(peaksM1[1]))
        intensitiesM1[is.na(intensitiesM1)] <- 0
      
        intensitiesM2 <- as.vector(intensity(peaksM2[1]))
        intensitiesM2[is.na(intensitiesM2)] <- 0
      
        intensitiesM0_180s <- as.vector(intensity(peaksM0_180s[1]))
        intensitiesM0_180s[is.na(intensitiesM0_180s)] <- 0
        total_slices_180s <- length(intensitiesM0_180s)
        slices_rt_180s <- as.vector(rtime(peaksM0_180s[1]))
        
        intensitiesM1_180s <- as.vector(intensity(peaksM1_180s[1]))
        intensitiesM1_180s[is.na(intensitiesM1_180s)] <- 0
        
        intensitiesM2_180s <- as.vector(intensity(peaksM2_180s[1]))
        intensitiesM2_180s[is.na(intensitiesM2_180s)] <- 0
        
        # Determine which slices have good quality and can be integrated to area
        # First initialize lists to tell which slice are selected
        integrated_to_area_list <- rep(0, total_slices)
        integrated_to_area_list_M0 <- rep(0, total_slices)
        integrated_to_area_list_M1 <- rep(0, total_slices)
        integrated_to_area_list_M2 <- rep(0, total_slices)
        integrated_to_area_list_all1 <- rep(1, total_slices)
        
        # SG smoothing to get the smoothed intensity vales
        intensitiesM0_sm <- sgolayfilt(intensitiesM0, p = 3, n = 7)
        intensitiesM1_sm <- sgolayfilt(intensitiesM1, p = 3, n = 7)
        intensitiesM2_sm <- sgolayfilt(intensitiesM2, p = 3, n = 7)
        intensitiesM0_180s_sm <- sgolayfilt(intensitiesM0_180s, p = 3, n = 7)
        intensitiesM1_180s_sm <- sgolayfilt(intensitiesM1_180s, p = 3, n = 7)
        intensitiesM2_180s_sm <- sgolayfilt(intensitiesM2_180s, p = 3, n = 7)
        
        # Find which slice is the peak height
        positionM0 <- which.max(intensitiesM0_sm)
        positionM1 <- which.max(intensitiesM1_sm)
        # Check if M1 and M2 peak height are different from M0, if it happens, it could mean M1 and M2 are contaminated or have high noises
        if (abs(positionM1-positionM0) > 5) {
          warntemp <- paste0("M1_peak at slice ", as.character(positionM1), " but M0_peak at slice ", as.character(positionM0))
          warningmessage <- paste(warningmessage, warntemp)
        }
        # If M1 or M2's peak height are found to be at least 5 slices away from M0 peak height, they are likely not true peak heights,
        # thus another slice within 5 slices of M0 peak height will be set to be the actual M1`or M2 peak height`
        startRange <- max(1, positionM0 - 5)
        endRange <- min(length(intensitiesM1_sm), positionM0 + 5)
        positionM1 <- startRange + which.max(intensitiesM1_sm[startRange:endRange]) - 1
        positionM2 <- which.max(intensitiesM2_sm)
        if (abs(positionM2-positionM0) > 5) {
          warntemp <- paste0("M2_peak at slice ", as.character(positionM2), " but M0_peak at slice ", as.character(positionM0))
          warningmessage <- paste(warningmessage, warntemp)
        }
        startRange <- max(1, positionM0 - 5)
        endRange <- min(length(intensitiesM2_sm), positionM0 + 5)
        positionM2 <- startRange + which.max(intensitiesM2_sm[startRange:endRange]) - 1
        
        # Use arPLS method to calculate peak baseline of M0, M1, and M2. Anything above the baseline are considered to be actual peak
        M0_baseline_180s <- arPLS_baseline_pos(intensitiesM0_180s)
        M0_baseline_180s[M0_baseline_180s > -1 & M0_baseline_180s < 1] <- 0
        M0_baseline <- rv_convert_baseline(M0_baseline_180s, slices_rt, slices_rt_180s)
        M0_edgeList <- pick_edges(intensitiesM0, intensitiesM0_sm, positionM0, M0_baseline)
        M0_edgeList_180s <- convert_List(M0_edgeList, slices_rt, slices_rt_180s)
        if (all(M0_edgeList_180s == 0)) {
          warningmessage <- paste(warningmessage, "There is no peak in M0")
        }
        M0_noise_180s <- estimate_noise(M0_edgeList_180s, intensitiesM0_180s, intensitiesM0_180s_sm, slices_rt_180s, M0_baseline_180s)
        
        M1_baseline_180s <- arPLS_baseline_pos(intensitiesM1_180s)
        M1_baseline_180s[M1_baseline_180s > -1 & M1_baseline_180s < 1] <- 0
        M1_baseline <- rv_convert_baseline(M1_baseline_180s, slices_rt, slices_rt_180s)
        M1_edgeList <- pick_edges(intensitiesM1, intensitiesM1_sm, positionM1, M1_baseline)
        M1_edgeList_180s <- convert_List(M1_edgeList, slices_rt, slices_rt_180s)
        if (all(M1_edgeList_180s == 0)) {
          warningmessage <- paste(warningmessage, "There is no peak in M1")
        }
        M1_noise_180s <- estimate_noise(M1_edgeList_180s, intensitiesM1_180s, intensitiesM1_180s_sm, slices_rt_180s, M1_baseline_180s)
        
        M2_baseline_180s <- arPLS_baseline_pos(intensitiesM2_180s)
        M2_baseline_180s[M2_baseline_180s > -1 & M2_baseline_180s < 1] <- 0
        M2_baseline <- rv_convert_baseline(M2_baseline_180s, slices_rt, slices_rt_180s)
        M2_edgeList <- pick_edges(intensitiesM2, intensitiesM2_sm, positionM2, M2_baseline)
        M2_edgeList_180s <- convert_List(M2_edgeList, slices_rt, slices_rt_180s)
        if (all(M2_edgeList_180s == 0)) {
          warningmessage <- paste(warningmessage, "There is no peak in M2")
        }
        M2_noise_180s <- estimate_noise(M2_edgeList_180s, intensitiesM2_180s, intensitiesM2_180s_sm, slices_rt_180s, M2_baseline_180s)
        
        # Check if the noises could be calculated correctly. If not, set a warning
        if (is.na(M0_noise_180s)) {
          warningmessage <- paste(warningmessage, "M0 noise is NA")
        }
        if (is.na(M1_noise_180s)) {
          warningmessage <- paste(warningmessage, "M1 noise is NA")
        }
        if (is.na(M2_noise_180s)) {
          warningmessage <- paste(warningmessage, "M2 noise is NA")
        }

        # Use function pick_slices to pick high quality slices
        integrated_to_area_list_M0 <- pick_slices(M0_edgeList, intensitiesM0, intensitiesM0_sm, positionM0, M0_baseline, M0_noise_180s, integrated_to_area_list_all1, intensitiesM0_sm, M0_baseline)
        integrated_to_area_list_M1 <- pick_slices(M1_edgeList, intensitiesM1, intensitiesM1_sm, positionM1, M1_baseline, M1_noise_180s, integrated_to_area_list_M0, intensitiesM0_sm, M0_baseline)
        integrated_to_area_list_M2 <- pick_slices(M2_edgeList, intensitiesM2, intensitiesM2_sm, positionM2, M2_baseline, M2_noise_180s, integrated_to_area_list_M0, intensitiesM0_sm, M0_baseline)
        # Check if both M0 and M1 slices are high quality or if all M0, M1, and M2 are all high quality to know what to be kept
        integrated_to_area_list <- as.integer(integrated_to_area_list_M0 & integrated_to_area_list_M1 & integrated_to_area_list_M2)
        integrated_to_area_list_onlyM0M1 <- as.integer(integrated_to_area_list_M0 & integrated_to_area_list_M1)
        
        
        # Calculating both unsmoothed and smoothed minus baseline intensity of each Mn, 
        intensitiesM0_subB <- intensitiesM0 - M0_baseline
        intensitiesM0_subB[intensitiesM0_subB < 0] <- 0
        intensitiesM1_subB <- intensitiesM1 - M1_baseline
        intensitiesM1_subB[intensitiesM1_subB < 0] <- 0
        intensitiesM2_subB <- intensitiesM2 - M2_baseline
        intensitiesM2_subB[intensitiesM2_subB < 0] <- 0
        intensitiesM0_subB_sm <- intensitiesM0_sm - M0_baseline
        intensitiesM0_subB_sm[intensitiesM0_subB_sm < 0] <- 0
        intensitiesM1_subB_sm <- intensitiesM1_sm - M1_baseline
        intensitiesM1_subB_sm[intensitiesM1_subB_sm < 0] <- 0
        intensitiesM2_subB_sm <- intensitiesM2_sm - M2_baseline
        intensitiesM2_subB_sm[intensitiesM2_subB_sm < 0] <- 0
        
        # Continue to repeat on other less important regions
        intensitiesM3 <- as.vector(intensity(peaksM3[1]))
        intensitiesM3[is.na(intensitiesM3)] <- 0
        intensitiesM3_180s <- as.vector(intensity(peaksM3_180s[1]))
        intensitiesM3_180s[is.na(intensitiesM3_180s)] <- 0
        intensitiesM3_sm <- sgolayfilt(intensitiesM3, p = 3, n = 7)
        intensitiesM3_180s_sm <- sgolayfilt(intensitiesM3_180s, p = 3, n = 7)
        M3_baseline_180s <- arPLS_baseline_pos(intensitiesM3_180s)
        M3_baseline_180s[M3_baseline_180s > -1 & M3_baseline_180s < 1] <- 0
        M3_baseline <- rv_convert_baseline(M3_baseline_180s, slices_rt, slices_rt_180s)
        
        startRange <- max(1, positionM0 - 5)
        endRange <- min(length(intensitiesM3_sm), positionM0 + 5)
        positionM3 <- startRange + which.max(intensitiesM3_sm[startRange:endRange]) - 1
        M3_edgeList <- pick_edges(intensitiesM3, intensitiesM3_sm, positionM3, M3_baseline)
        M3_edgeList_180s <- convert_List(M3_edgeList, slices_rt, slices_rt_180s)
        if (all(M3_edgeList_180s == 0)) {
          warningmessage <- paste(warningmessage, "There is no peak in M3")
        }
        M3_noise_180s <- estimate_noise(M3_edgeList_180s, intensitiesM3_180s, intensitiesM3_180s_sm, slices_rt_180s, M3_baseline_180s)
        
        intensitiesM3_subB <- intensitiesM3 - M3_baseline
        intensitiesM3_subB[intensitiesM3_subB < 0] <- 0
        intensitiesM3_subB_sm <- intensitiesM3_sm - M3_baseline
        intensitiesM3_subB_sm[intensitiesM3_subB_sm < 0] <- 0
        
        intensitiestotal <- as.vector(intensity(peakstotal[1]))
        intensitiestotal[is.na(intensitiestotal)] <- 0
        intensitiestotal_sm <- sgolayfilt(intensitiestotal, p = 3, n = 7)
        intensitiestotal_sm[intensitiestotal_sm < 0] <- 0
        
        intensitiesM0_sm[intensitiesM0_sm < 0] <- 0
        intensitiesM1_sm[intensitiesM1_sm < 0] <- 0
        intensitiesM2_sm[intensitiesM2_sm < 0] <- 0
        intensitiesM3_sm[intensitiesM3_sm < 0] <- 0
      
        intensitiesMminus1_to_M0 <- as.vector(intensity(peaksMminus1_to_M0[1]))
        intensitiesMminus1_to_M0[is.na(intensitiesMminus1_to_M0)] <- 0
        intensitiesM0_to_M1 <- as.vector(intensity(peaksM0_to_M1[1]))
        intensitiesM0_to_M1[is.na(intensitiesM0_to_M1)] <- 0
        intensitiesM1_to_M2 <- as.vector(intensity(peaksM1_to_M2[1]))
        intensitiesM1_to_M2[is.na(intensitiesM1_to_M2)] <- 0
        intensitiesM2_to_M3 <- as.vector(intensity(peaksM2_to_M3[1]))
        intensitiesM2_to_M3[is.na(intensitiesM2_to_M3)] <- 0
        intensitiesM3_to_M4 <- as.vector(intensity(peaksM3_to_M4[1]))
        intensitiesM3_to_M4[is.na(intensitiesM3_to_M4)] <- 0
        
        # Check if formula is illegal, which could often be an issue due to uncompleted reference library, 
        # si that one missing formula metabolite won't cause all other metabolites to be re-extracted.
        if (formula != "" & formula != "null") 
          {
            # Calculate Percentage
            passing_formula <- elementscount$Count
            predicted_M012 <- Abundance$get_MID(Abundance(passing_formula))
            predicted_M012Sum <- predicted_M012[1] + predicted_M012[2] + predicted_M012[3]
            predicted_M0 <- predicted_M012[1]/predicted_M012Sum
            predicted_M1 <- predicted_M012[2]/predicted_M012Sum
            predicted_M2 <- predicted_M012[3]/predicted_M012Sum
            
            for (s in 1:total_slices) 
              {
                # Calculate values for each slices
                M012_slice <- intensitiesM0[s] + intensitiesM1[s] + intensitiesM2[s]
                M0_slice <- intensitiesM0[s]/M012_slice
                M1_slice <- intensitiesM1[s]/M012_slice
                M2_slice <- intensitiesM2[s]/M012_slice
                dM0 <- M0_slice - predicted_M0
                dM1 <- M1_slice - predicted_M1
                dM2 <- M2_slice - predicted_M2  
                
                M012_slice_subB <- intensitiesM0_subB[s] + intensitiesM1_subB[s] + intensitiesM2_subB[s]
                M0_slice_subB <- intensitiesM0_subB[s]/M012_slice_subB
                M1_slice_subB <- intensitiesM1_subB[s]/M012_slice_subB
                M2_slice_subB <- intensitiesM2_subB[s]/M012_slice_subB
                dM0_subB <- M0_slice_subB - predicted_M0
                dM1_subB <- M1_slice_subB - predicted_M1
                dM2_subB <- M2_slice_subB - predicted_M2  
                
                M012_slice_sm <- intensitiesM0_sm[s] + intensitiesM1_sm[s] + intensitiesM2_sm[s]
                M0_slice_sm <- intensitiesM0_sm[s]/M012_slice_sm
                M1_slice_sm <- intensitiesM1_sm[s]/M012_slice_sm
                M2_slice_sm <- intensitiesM2_sm[s]/M012_slice_sm
                dM0_sm <- M0_slice_sm - predicted_M0
                dM1_sm <- M1_slice_sm - predicted_M1
                dM2_sm <- M2_slice_sm - predicted_M2  
                
                M012_slice_subB_sm <- intensitiesM0_subB_sm[s] + intensitiesM1_subB_sm[s] + intensitiesM2_subB_sm[s]
                M0_slice_subB_sm <- intensitiesM0_subB_sm[s]/M012_slice_subB_sm
                M1_slice_subB_sm <- intensitiesM1_subB_sm[s]/M012_slice_subB_sm
                M2_slice_subB_sm <- intensitiesM2_subB_sm[s]/M012_slice_subB_sm
                dM0_subB_sm <- M0_slice_subB_sm - predicted_M0
                dM1_subB_sm <- M1_slice_subB_sm - predicted_M1
                dM2_subB_sm <- M2_slice_subB_sm - predicted_M2  
                
                extracted_slices <- rbind(extracted_slices, 
                                          data.frame(
                                                      Compound = compound,
                                                      Sample = sample_name,
                                                      Formula = formula,
                                                      Compound_Abundance = heightnumber,
                                                      RT_left = left_value/60,
                                                      RT_peak = middle_value/60,
                                                      RT_right = right_value/60,
                                                      Slice_rt = slices_rt[s]/60,
                                                      Slice = s,
                                                      Total_Slices = total_slices,
                                                      Slice_Peak_Height = positionM0,
                                                      Integrated_To_Area = integrated_to_area_list[s],
                                                      Integrated_To_Area_onlyM0M1 = integrated_to_area_list_onlyM0M1[s],
                                                      Integrated_To_Area_M0 = integrated_to_area_list_M0[s],
                                                      Integrated_To_Area_M1 = integrated_to_area_list_M1[s],
                                                      Integrated_To_Area_M2 = integrated_to_area_list_M2[s],
                                                      Reference_mz = reference_mz,
                                                      Precursor_mz = mz_M0,
                                                      M1_mz = (mz_M1_low + mz_M1_high)/2,
                                                      M2_mz = (mz_M2_low + mz_M2_high)/2,  
                                                      M0_Predicted = predicted_M0,
                                                      M1_Predicted = predicted_M1,
                                                      M2_Predicted = predicted_M2,
                                                      M0_Observed = M0_slice,
                                                      M1_Observed = M1_slice,
                                                      M2_Observed = M2_slice,
                                                      DeltaM0	= dM0,
                                                      DeltaM1	= dM1,
                                                      DeltaM2 = dM2,
                                                      DeltaM0_Percentage	= dM0/predicted_M0,
                                                      DeltaM1_Percentage	= dM1/predicted_M1,
                                                      DeltaM2_Percentage  = dM2/predicted_M2,
                                                      Total_Ion = intensitiestotal[s],
                                                      Total_Ion_Smoothed = intensitiestotal_sm[s],
                                                      M0_Intensity = intensitiesM0[s], 
                                                      M1_Intensity = intensitiesM1[s], 
                                                      M2_Intensity = intensitiesM2[s],
                                                      M3_Intensity = intensitiesM3[s],
                                                      M0_Baseline = M0_baseline[s],
                                                      M1_Baseline = M1_baseline[s], 
                                                      M2_Baseline = M2_baseline[s],
                                                      M3_Baseline = M3_baseline[s],
                                                      M0_Noise = M0_noise_180s, 
                                                      M1_Noise = M1_noise_180s, 
                                                      M2_Noise = M2_noise_180s,
                                                      M3_Noise = M3_noise_180s,
                                                      Mminus1_to_M0_Intensity = intensitiesMminus1_to_M0[s], 
                                                      M0_to_M1_Intensity = intensitiesM0_to_M1[s], 
                                                      M1_to_M2_Intensity = intensitiesM1_to_M2[s],
                                                      M2_to_M3_Intensity = intensitiesM2_to_M3[s],
                                                      M3_to_M4_Intensity = intensitiesM3_to_M4[s],        
                                                      M0_Observed_subB = M0_slice_subB,
                                                      M1_Observed_subB = M1_slice_subB,
                                                      M2_Observed_subB = M2_slice_subB,
                                                      DeltaM0_subB	= dM0_subB,
                                                      DeltaM1_subB	= dM1_subB,
                                                      DeltaM2_subB = dM2_subB,
                                                      DeltaM0_Percentage_subB	= dM0_subB/predicted_M0,
                                                      DeltaM1_Percentage_subB	= dM1_subB/predicted_M1,
                                                      DeltaM2_Percentage_subB  = dM2_subB/predicted_M2,
                                                      M0_Intensity_subB = intensitiesM0_subB[s], 
                                                      M1_Intensity_subB = intensitiesM1_subB[s], 
                                                      M2_Intensity_subB = intensitiesM2_subB[s],
                                                      M3_Intensity_subB = intensitiesM3_subB[s],
                                                      M0_Smoothed_Observed = M0_slice_sm,
                                                      M1_Smoothed_Observed = M1_slice_sm,
                                                      M2_Smoothed_Observed = M2_slice_sm,
                                                      DeltaM0_Smoothed	= dM0_sm,
                                                      DeltaM1_Smoothed	= dM1_sm,
                                                      DeltaM2_Smoothed = dM2_sm,
                                                      DeltaM0_Smoothed_Percentage	= dM0_sm/predicted_M0,
                                                      DeltaM1_Smoothed_Percentage	= dM1_sm/predicted_M1,
                                                      DeltaM2_Smoothed_Percentage  = dM2_sm/predicted_M2,
                                                      M0_Intensity_Smoothed = intensitiesM0_sm[s], 
                                                      M1_Intensity_Smoothed = intensitiesM1_sm[s], 
                                                      M2_Intensity_Smoothed = intensitiesM2_sm[s],
                                                      M3_Intensity_Smoothed = intensitiesM3_sm[s],
                                                      M0_Smoothed_Observed_subB = M0_slice_subB_sm,
                                                      M1_Smoothed_Observed_subB = M1_slice_subB_sm,
                                                      M2_Smoothed_Observed_subB = M2_slice_subB_sm,
                                                      DeltaM0_Smoothed_subB	= dM0_subB_sm,
                                                      DeltaM1_Smoothed_subB	= dM1_subB_sm,
                                                      DeltaM2_Smoothed_subB = dM2_subB_sm,
                                                      DeltaM0_Smoothed_Percentage_subB	= dM0_subB_sm/predicted_M0,
                                                      DeltaM1_Smoothed_Percentage_subB	= dM1_subB_sm/predicted_M1,
                                                      DeltaM2_Smoothed_Percentage_subB  = dM2_subB_sm/predicted_M2,
                                                      M0_Intensity_Smoothed_subB = intensitiesM0_subB_sm[s], 
                                                      M1_Intensity_Smoothed_subB = intensitiesM1_subB_sm[s], 
                                                      M2_Intensity_Smoothed_subB = intensitiesM2_subB_sm[s],
                                                      M3_Intensity_Smoothed_subB = intensitiesM3_subB_sm[s],
                                                      Warning = warningmessage
                                                    )
                                          )
              }
            
            
        
            
            
            # Calculate Integrated Area Delta, and only do that if there are at least 5 high-quality slices 
            if (sum(sapply(integrated_to_area_list, function(x) x == 1)) >= 5)
              {
                # Get correlation score between M0 and M1, and M0 and M2 of the full peak area              
                all_edge_list <-   M0_edgeList & M1_edgeList & M2_edgeList
                left_slice_p <- which(all_edge_list == 1)[1]
                right_slice_p <- tail(which(all_edge_list == 1), n = 1)
                score_M1_edge <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "spearman"))
                score_M2_edge <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM2_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM2_sm[left_slice_p: right_slice_p], method = "spearman"))
                
                # Get correlation score between M0 and M1, and M0 and M2 of only the area with slices considered high quality    
                ITA_list <- extracted_slices$Integrated_To_Area
                left_slice_p <- which(ITA_list == 1)[1]
                right_slice_p <- tail(which(ITA_list == 1), n = 1)
                score_M1 <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "spearman"))
                score_M2 <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM2_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM2_sm[left_slice_p: right_slice_p], method = "spearman"))
                
                # If Peak M1 and M2 has very different shape from M0 based on correlation_cutoff, they will not be included to do M0M1M2 ratio calculations
                if (min(score_M1, score_M2) < correlation_cutoff) 
                  {
                  integrated_to_area_list <- rep(0, total_slices)
                  integrated_to_area_list_M2 <- rep(0, total_slices)
                  extracted_slices$Integrated_To_Area <- rep(0, total_slices)
                  extracted_slices$Integrated_To_Area_M2 <- rep(0, total_slices)
                  } 
                  else 
                    {
                      # Additional check to see if the selected region of Peak M1 and M2's correlation is different from full peak area similarity
                      if (abs(min(score_M1_edge, score_M2_edge) - min(score_M1, score_M2)) > 0.1) 
                        {
                        warningmessage <- paste(warningmessage, "Selected slices' max(pearson, spearman) score is >0.1 different than the full peak max(pearson, spearman) score")
                        }
                      # If everything looks good, start to integrate data from all high-quality slices.
                      areasM0 <- sum(extracted_slices$M0_Intensity[extracted_slices$Integrated_To_Area == 1])
                      areasM1 <- sum(extracted_slices$M1_Intensity[extracted_slices$Integrated_To_Area == 1])
                      areasM2 <- sum(extracted_slices$M2_Intensity[extracted_slices$Integrated_To_Area == 1])
                      areasM3 <- sum(extracted_slices$M3_Intensity[extracted_slices$Integrated_To_Area == 1])
                      deltaM0 <- areasM0/(areasM0+areasM1+areasM2) - predicted_M0
                      deltaM1 <- areasM1/(areasM0+areasM1+areasM2) - predicted_M1
                      deltaM2 <- areasM2/(areasM0+areasM1+areasM2) - predicted_M2
                      
                      areastotal <- sum(extracted_slices$Total_Ion[extracted_slices$Integrated_To_Area == 1])
                      areastotal_sm <- sum(extracted_slices$Total_Ion_Smoothed[extracted_slices$Integrated_To_Area == 1])
                      areasMminus1_to_M0 <- sum(extracted_slices$Mminus1_to_M0_Intensity[extracted_slices$Integrated_To_Area == 1])      
                      areasM0_to_M1 <- sum(extracted_slices$M0_to_M1_Intensity[extracted_slices$Integrated_To_Area == 1])
                      areasM1_to_M2 <- sum(extracted_slices$M1_to_M2_Intensity[extracted_slices$Integrated_To_Area == 1])
                      areasM2_to_M3 <- sum(extracted_slices$M2_to_M3_Intensity[extracted_slices$Integrated_To_Area == 1])      
                      areasM3_to_M4 <- sum(extracted_slices$M3_to_M4_Intensity[extracted_slices$Integrated_To_Area == 1])  
                      
                      areasM0_subB <- sum(extracted_slices$M0_Intensity_subB[extracted_slices$Integrated_To_Area == 1] )
                      areasM1_subB <- sum(extracted_slices$M1_Intensity_subB[extracted_slices$Integrated_To_Area == 1])
                      areasM2_subB <- sum(extracted_slices$M2_Intensity_subB[extracted_slices$Integrated_To_Area == 1])
                      deltaM0_subB <- areasM0_subB/(areasM0_subB+areasM1_subB+areasM2_subB) - predicted_M0
                      deltaM1_subB <- areasM1_subB/(areasM0_subB+areasM1_subB+areasM2_subB) - predicted_M1
                      deltaM2_subB <- areasM2_subB/(areasM0_subB+areasM1_subB+areasM2_subB) - predicted_M2
                      
                      areasM0_sm <- sum(extracted_slices$M0_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1])
                      areasM1_sm <- sum(extracted_slices$M1_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1])
                      areasM2_sm <- sum(extracted_slices$M2_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1])
                      deltaM0_sm <- areasM0_sm/(areasM0_sm+areasM1_sm+areasM2_sm) - predicted_M0
                      deltaM1_sm <- areasM1_sm/(areasM0_sm+areasM1_sm+areasM2_sm) - predicted_M1
                      deltaM2_sm <- areasM2_sm/(areasM0_sm+areasM1_sm+areasM2_sm) - predicted_M2
                      
                      areasM0_subB_sm <- sum(extracted_slices$M0_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1] )
                      areasM1_subB_sm <- sum(extracted_slices$M1_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1])
                      areasM2_subB_sm <- sum(extracted_slices$M2_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1])
                      deltaM0_subB_sm <- areasM0_subB_sm/(areasM0_subB_sm+areasM1_subB_sm+areasM2_subB_sm) - predicted_M0
                      deltaM1_subB_sm <- areasM1_subB_sm/(areasM0_subB_sm+areasM1_subB_sm+areasM2_subB_sm) - predicted_M1
                      deltaM2_subB_sm <- areasM2_subB_sm/(areasM0_subB_sm+areasM1_subB_sm+areasM2_subB_sm) - predicted_M2
                      
                      # Append the extracted information to the table
                      extracted_table <- rbind(extracted_table, 
                                               data.frame(
                                                            Compound = compound,
                                                            Sample = sample_name,
                                                            Formula = formula,
                                                            Compound_Abundance = heightnumber,
                                                            RT_left = left_value/60,
                                                            RT_peak = middle_value/60,
                                                            RT_right = right_value/60,
                                                            Reference_mz = reference_mz,
                                                            Precursor_mz = mz_M0,
                                                            M1_mz = (mz_M1_low + mz_M1_high)/2,
                                                            M2_mz = (mz_M2_low + mz_M2_high)/2,
                                                            Integrated_To_Area = 1,
                                                            Integrated_To_Area_onlyM0M1 = 1,
                                                            M0_Observed = areasM0/(areasM0+areasM1+areasM2),
                                                            M1_Observed = areasM1/(areasM0+areasM1+areasM2),
                                                            M2_Observed = areasM2/(areasM0+areasM1+areasM2),
                                                            M0_Predicted = predicted_M0,
                                                            M1_Predicted = predicted_M1,
                                                            M2_Predicted = predicted_M2,
                                                            DeltaM0	= deltaM0,
                                                            DeltaM1	= deltaM1,
                                                            DeltaM2 = deltaM2,
                                                            DeltaM0_Percentage	= deltaM0/predicted_M0,
                                                            DeltaM1_Percentage	= deltaM1/predicted_M1,
                                                            DeltaM2_Percentage  = deltaM2/predicted_M2,
                                                            Total_Ion = areastotal,
                                                            Total_Ion_Smoothed = areastotal_sm,
                                                            M0_Area = areasM0, 
                                                            M1_Area = areasM1, 
                                                            M2_Area = areasM2,
                                                            M3_Area = areasM3,
                                                            M0_Height = max(extracted_slices$M0_Intensity[extracted_slices$Integrated_To_Area == 1]), 
                                                            M1_Height = max(extracted_slices$M1_Intensity[extracted_slices$Integrated_To_Area == 1]),  
                                                            M2_Height = max(extracted_slices$M2_Intensity[extracted_slices$Integrated_To_Area == 1]), 
                                                            M3_Height = max(extracted_slices$M3_Intensity[extracted_slices$Integrated_To_Area == 1]), 
                                                            M0_Baseline = mean(extracted_slices$M0_Baseline[extracted_slices$Integrated_To_Area == 1]),
                                                            M1_Baseline = mean(extracted_slices$M1_Baseline[extracted_slices$Integrated_To_Area == 1]),
                                                            M2_Baseline = mean(extracted_slices$M2_Baseline[extracted_slices$Integrated_To_Area == 1]),
                                                            M0_Noise = M0_noise_180s, 
                                                            M1_Noise = M1_noise_180s, 
                                                            M2_Noise = M2_noise_180s,
                                                            Mminus1_to_M0_Area = areasMminus1_to_M0, 
                                                            M0_to_M1_Area = areasM0_to_M1,
                                                            M1_to_M2_Area = areasM1_to_M2,
                                                            M2_to_M3_Area = areasM2_to_M3,
                                                            M3_to_M4_Area = areasM3_to_M4,   
                                                            M0_Observed_subB = areasM0_subB/(areasM0_subB+areasM1_subB+areasM2_subB), 
                                                            M1_Observed_subB = areasM1_subB/(areasM0_subB+areasM1_subB+areasM2_subB), 
                                                            M2_Observed_subB = areasM2_subB/(areasM0_subB+areasM1_subB+areasM2_subB), 
                                                            DeltaM0_subB	= deltaM0_subB,
                                                            DeltaM1_subB	= deltaM1_subB,
                                                            DeltaM2_subB = deltaM2_subB,
                                                            DeltaM0_Percentage_subB	= deltaM0_subB/predicted_M0,
                                                            DeltaM1_Percentage_subB	= deltaM1_subB/predicted_M1,
                                                            DeltaM2_Percentage_subB  = deltaM2_subB/predicted_M2,
                                                            M0_Area_subB = areasM0_subB, 
                                                            M1_Area_subB = areasM1_subB, 
                                                            M2_Area_subB = areasM2_subB,
                                                            M0_Height_subB = max(extracted_slices$M0_Intensity_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M1_Height_subB = max(extracted_slices$M1_Intensity_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M2_Height_subB = max(extracted_slices$M2_Intensity_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M0_Smoothed_Observed = areasM0_sm/(areasM0_sm+areasM1_sm+areasM2_sm),
                                                            M1_Smoothed_Observed = areasM1_sm/(areasM0_sm+areasM1_sm+areasM2_sm),
                                                            M2_Smoothed_Observed = areasM2_sm/(areasM0_sm+areasM1_sm+areasM2_sm),
                                                            DeltaM0_Smoothed	= deltaM0_sm,
                                                            DeltaM1_Smoothed	= deltaM1_sm,
                                                            DeltaM2_Smoothed = deltaM2_sm,
                                                            DeltaM0_Smoothed_Percentage	= deltaM0_sm/predicted_M0,
                                                            DeltaM1_Smoothed_Percentage	= deltaM1_sm/predicted_M1,
                                                            DeltaM2_Smoothed_Percentage  = deltaM2_sm/predicted_M2,
                                                            M0_Area_Smoothed = areasM0_sm, 
                                                            M1_Area_Smoothed = areasM1_sm, 
                                                            M2_Area_Smoothed = areasM2_sm,
                                                            M0_Height_Smoothed = max(extracted_slices$M0_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1]), 
                                                            M1_Height_Smoothed = max(extracted_slices$M1_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1]), 
                                                            M2_Height_Smoothed = max(extracted_slices$M2_Intensity_Smoothed[extracted_slices$Integrated_To_Area == 1]), 
                                                            DeltaM0_Smoothed_subB	= deltaM0_subB_sm,
                                                            DeltaM1_Smoothed_subB	= deltaM1_subB_sm,
                                                            DeltaM2_Smoothed_subB = deltaM2_subB_sm,
                                                            DeltaM0_Smoothed_Percentage_subB	= deltaM0_subB_sm/predicted_M0,
                                                            DeltaM1_Smoothed_Percentage_subB	= deltaM1_subB_sm/predicted_M1,
                                                            DeltaM2_Smoothed_Percentage_subB  = deltaM2_subB_sm/predicted_M2,
                                                            M0_Area_Smoothed_subB = areasM0_subB_sm, 
                                                            M1_Area_Smoothed_subB = areasM1_subB_sm, 
                                                            M2_Area_Smoothed_subB = areasM2_subB_sm,
                                                            M0_Height_Smoothed_subB = max(extracted_slices$M0_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M1_Height_Smoothed_subB = max(extracted_slices$M1_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M2_Height_Smoothed_subB = max(extracted_slices$M2_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area == 1]),
                                                            M1_FullPeak_max_pea_spe_score = score_M1_edge,
                                                            M1_Slices_max_pea_spe_score = score_M1,
                                                            M2_FullPeak_max_pea_spe_score = score_M2_edge,
                                                            M2_Slices_max_pea_spe_score = score_M2,
                                                            Warning = warningmessage
                                                          )
                                                )
                    } 
              }
            
              # Ending the condition of M0M1M2 all have enough high-quality slices, and stard dealing with only if M0 and M1 have enough high-quality slices            
              if (sum(sapply(integrated_to_area_list, function(x) x == 1)) < 5 & sum(sapply(integrated_to_area_list_onlyM0M1, function(x) x == 1)) >= 5)
                {
                  pM0M1sum <- predicted_M0 + predicted_M1
                  predicted_M0 <- predicted_M0/pM0M1sum
                  predicted_M1 <- predicted_M1/pM0M1sum
                  predicted_M2 <- NA
          
                  areasM0 <- sum(extracted_slices$M0_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM1 <- sum(extracted_slices$M1_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM2 <- sum(extracted_slices$M2_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM3 <- sum(extracted_slices$M3_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  deltaM0 <- areasM0/(areasM0+areasM1) - predicted_M0
                  deltaM1 <- areasM1/(areasM0+areasM1) - predicted_M1
                  deltaM2 <- NA
                  
                  areastotal <- sum(extracted_slices$Total_Ion[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areastotal_sm <- sum(extracted_slices$Total_Ion_Smoothed[extracted_slices$Integrated_To_Area == 1])
                  areasMminus1_to_M0 <- sum(extracted_slices$Mminus1_to_M0_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])      
                  areasM0_to_M1 <- sum(extracted_slices$M0_to_M1_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM1_to_M2 <- sum(extracted_slices$M1_to_M2_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM2_to_M3 <- sum(extracted_slices$M2_to_M3_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])      
                  areasM3_to_M4 <- sum(extracted_slices$M3_to_M4_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])  
                  
                  areasM0_subB <- sum(extracted_slices$M0_Intensity_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1] )
                  areasM1_subB <- sum(extracted_slices$M1_Intensity_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM2_subB <- NA
                  deltaM0_subB <- areasM0_subB/(areasM0_subB+areasM1_subB) - predicted_M0
                  deltaM1_subB <- areasM1_subB/(areasM0_subB+areasM1_subB) - predicted_M1
                  deltaM2_subB <- NA
                  
                  all_edge_list <- M0_edgeList & M1_edgeList
                  left_slice_p <- which(all_edge_list == 1)[1]
                  right_slice_p <- tail(which(all_edge_list == 1), n = 1)
                  score_M1_edge <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "spearman"))
                  
                  ITA_list <- extracted_slices$Integrated_To_Area_onlyM0M1
                  left_slice_p <- which(ITA_list == 1)[1]
                  right_slice_p <- tail(which(ITA_list == 1), n = 1)
                  score_M1 <- max(cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "pearson"), cor(intensitiesM0_sm[left_slice_p: right_slice_p], intensitiesM1_sm[left_slice_p: right_slice_p], method = "spearman"))
                  if (abs(score_M1_edge - score_M1) > 0.1) {
                    warningmessage <- paste(warningmessage, "Selected slices' max(pearson, spearman) score is >0.1 different than the full peak max(pearson, spearman) score")
                  }
                  
                  areasM0_sm <- sum(extracted_slices$M0_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM1_sm <- sum(extracted_slices$M1_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM2_sm <- sum(extracted_slices$M2_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  deltaM0_sm <- areasM0_sm/(areasM0_sm+areasM1_sm) - predicted_M0
                  deltaM1_sm <- areasM1_sm/(areasM0_sm+areasM1_sm) - predicted_M1
                  deltaM2_sm <- NA
                  
                  areasM0_subB_sm <- sum(extracted_slices$M0_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1] )
                  areasM1_subB_sm <- sum(extracted_slices$M1_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1])
                  areasM2_subB_sm <- NA
                  deltaM0_subB_sm <- areasM0_subB_sm/(areasM0_subB_sm+areasM1_subB_sm) - predicted_M0
                  deltaM1_subB_sm <- areasM1_subB_sm/(areasM0_subB_sm+areasM1_subB_sm) - predicted_M1
                  deltaM2_subB_sm <- NA
                  
                  # Append the extracted information to the table
                  extracted_table <- rbind(extracted_table, 
                                           data.frame(
                                                        Compound = compound,
                                                        Sample = sample_name,
                                                        Formula = formula,
                                                        Compound_Abundance = heightnumber,
                                                        RT_left = left_value/60,
                                                        RT_peak = middle_value/60,
                                                        RT_right = right_value/60,
                                                        Reference_mz = reference_mz,
                                                        Precursor_mz = mz_M0,
                                                        M1_mz = (mz_M1_low + mz_M1_high)/2,
                                                        M2_mz = (mz_M2_low + mz_M2_high)/2,
                                                        Integrated_To_Area = 0,
                                                        Integrated_To_Area_onlyM0M1 = 1,
                                                        M0_Predicted = predicted_M0,
                                                        M1_Predicted = predicted_M1,
                                                        M2_Predicted = NA,
                                                        M0_Observed = areasM0/(areasM0+areasM1),
                                                        M1_Observed = areasM1/(areasM0+areasM1),
                                                        M2_Observed = NA,
                                                        DeltaM0	= deltaM0,
                                                        DeltaM1	= deltaM1,
                                                        DeltaM2 = NA,
                                                        DeltaM0_Percentage	= deltaM0/predicted_M0,
                                                        DeltaM1_Percentage	= deltaM1/predicted_M1,
                                                        DeltaM2_Percentage  = NA,
                                                        Total_Ion = areastotal,
                                                        Total_Ion_Smoothed = areastotal_sm,
                                                        M0_Area = areasM0, 
                                                        M1_Area = areasM1, 
                                                        M2_Area = areasM2,
                                                        M3_Area = areasM3,
                                                        M0_Height = max(extracted_slices$M0_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        M1_Height = max(extracted_slices$M1_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),  
                                                        M2_Height = max(extracted_slices$M2_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        M3_Height = max(extracted_slices$M3_Intensity[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        M0_Baseline = mean(extracted_slices$M0_Baseline[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M1_Baseline = mean(extracted_slices$M1_Baseline[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M2_Baseline = NA,
                                                        M0_Noise = M0_noise_180s, 
                                                        M1_Noise = M1_noise_180s, 
                                                        M2_Noise = NA,
                                                        Mminus1_to_M0_Area = areasMminus1_to_M0, 
                                                        M0_to_M1_Area = areasM0_to_M1,
                                                        M1_to_M2_Area = areasM1_to_M2,
                                                        M2_to_M3_Area = areasM2_to_M3,
                                                        M3_to_M4_Area = areasM3_to_M4,   
                                                        M0_Observed_subB = areasM0_subB/(areasM0_subB+areasM1_subB), 
                                                        M1_Observed_subB = areasM1_subB/(areasM0_subB+areasM1_subB), 
                                                        M2_Observed_subB = NA, 
                                                        DeltaM0_subB	= deltaM0_subB,
                                                        DeltaM1_subB	= deltaM1_subB,
                                                        DeltaM2_subB = NA,
                                                        DeltaM0_Percentage_subB	= deltaM0_subB/predicted_M0,
                                                        DeltaM1_Percentage_subB	= deltaM1_subB/predicted_M1,
                                                        DeltaM2_Percentage_subB  = NA,
                                                        M0_Area_subB = areasM0_subB, 
                                                        M1_Area_subB = areasM1_subB, 
                                                        M2_Area_subB = NA,
                                                        M0_Height_subB = max(extracted_slices$M0_Intensity_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M1_Height_subB = max(extracted_slices$M1_Intensity_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M2_Height_subB = NA,
                                                        M0_Smoothed_Observed = areasM0_sm/(areasM0_sm+areasM1_sm),
                                                        M1_Smoothed_Observed = areasM1_sm/(areasM0_sm+areasM1_sm),
                                                        M2_Smoothed_Observed = NA,
                                                        DeltaM0_Smoothed	= deltaM0_sm,
                                                        DeltaM1_Smoothed	= deltaM1_sm,
                                                        DeltaM2_Smoothed = NA,
                                                        DeltaM0_Smoothed_Percentage	= deltaM0_sm/predicted_M0,
                                                        DeltaM1_Smoothed_Percentage	= deltaM1_sm/predicted_M1,
                                                        DeltaM2_Smoothed_Percentage  = NA,
                                                        M0_Area_Smoothed = areasM0_sm, 
                                                        M1_Area_Smoothed = areasM1_sm, 
                                                        M2_Area_Smoothed = areasM2_sm,
                                                        M0_Height_Smoothed = max(extracted_slices$M0_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        M1_Height_Smoothed = max(extracted_slices$M1_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        M2_Height_Smoothed = max(extracted_slices$M2_Intensity_Smoothed[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]), 
                                                        DeltaM0_Smoothed_subB	= deltaM0_subB_sm,
                                                        DeltaM1_Smoothed_subB	= deltaM1_subB_sm,
                                                        DeltaM2_Smoothed_subB = NA,
                                                        DeltaM0_Smoothed_Percentage_subB	= deltaM0_subB_sm/predicted_M0,
                                                        DeltaM1_Smoothed_Percentage_subB	= deltaM1_subB_sm/predicted_M1,
                                                        DeltaM2_Smoothed_Percentage_subB  = NA,
                                                        M0_Area_Smoothed_subB = areasM0_subB_sm, 
                                                        M1_Area_Smoothed_subB = areasM1_subB_sm, 
                                                        M2_Area_Smoothed_subB = NA,
                                                        M0_Height_Smoothed_subB = max(extracted_slices$M0_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M1_Height_Smoothed_subB = max(extracted_slices$M1_Intensity_Smoothed_subB[extracted_slices$Integrated_To_Area_onlyM0M1 == 1]),
                                                        M2_Height_Smoothed_subB = NA,
                                                        M1_FullPeak_max_pea_spe_score = score_M1_edge,
                                                        M1_Slices_max_pea_spe_score = score_M1,
                                                        M2_FullPeak_max_pea_spe_score = NA,
                                                        M2_Slices_max_pea_spe_score = NA,
                                                        Warning = warningmessage
                                                    )
                  )
                } else if(sum(sapply(integrated_to_area_list, function(x) x == 1)) < 5 & sum(sapply(integrated_to_area_list_onlyM0M1, function(x) x == 1)) < 5)
                    { 
                      warningmessage <- paste(warningmessage, "Fewer than three five have M0M1 all above noises")
                      extracted_table <- rbind(extracted_table, 
                                               data.frame(
                                                            Compound = compound,
                                                            Sample = sample_name,
                                                            Formula = formula,
                                                            Compound_Abundance = heightnumber,
                                                            RT_left = left_value/60,
                                                            RT_peak = middle_value/60,
                                                            RT_right = right_value/60,
                                                            Reference_mz = reference_mz,
                                                            Precursor_mz = mz_M0,
                                                            M1_mz = (mz_M1_low + mz_M1_high)/2,
                                                            M2_mz = (mz_M2_low + mz_M2_high)/2,    
                                                            Integrated_To_Area = 0,
                                                            Integrated_To_Area_onlyM0M1 = 0,
                                                            M0_Predicted = predicted_M0,
                                                            M1_Predicted = predicted_M1,
                                                            M2_Predicted = predicted_M2,
                                                            M0_Observed = NA,
                                                            M1_Observed = NA,
                                                            M2_Observed = NA,
                                                            DeltaM0	= NA,
                                                            DeltaM1	= NA,
                                                            DeltaM2 = NA,
                                                            DeltaM0_Percentage	= NA,
                                                            DeltaM1_Percentage	= NA,
                                                            DeltaM2_Percentage  = NA,      
                                                            Total_Ion = NA,
                                                            Total_Ion_Smoothed = NA,
                                                            M0_Area = NA, 
                                                            M1_Area = NA, 
                                                            M2_Area = NA,
                                                            M3_Area = NA,
                                                            M0_Height = NA, 
                                                            M1_Height = NA, 
                                                            M2_Height = NA,
                                                            M3_Height = NA,
                                                            M0_Baseline = NA,
                                                            M1_Baseline = NA,
                                                            M2_Baseline = NA,
                                                            M0_Noise = NA, 
                                                            M1_Noise = NA, 
                                                            M2_Noise = NA,
                                                            Mminus1_to_M0_Area = NA, 
                                                            M0_to_M1_Area = NA,
                                                            M1_to_M2_Area = NA,
                                                            M2_to_M3_Area = NA,
                                                            M3_to_M4_Area = NA,
                                                            M0_Observed_subB = NA,
                                                            M1_Observed_subB = NA,
                                                            M2_Observed_subB = NA,
                                                            DeltaM0_subB	= NA,
                                                            DeltaM1_subB	= NA,
                                                            DeltaM2_subB = NA,
                                                            DeltaM0_Percentage_subB	= NA,
                                                            DeltaM1_Percentage_subB	= NA,
                                                            DeltaM2_Percentage_subB  = NA,
                                                            M0_Area_subB = NA, 
                                                            M1_Area_subB = NA, 
                                                            M2_Area_subB = NA,
                                                            M0_Height_subB = NA, 
                                                            M1_Height_subB = NA, 
                                                            M2_Height_subB = NA,
                                                            M0_Smoothed_Observed = NA,
                                                            M1_Smoothed_Observed = NA,
                                                            M2_Smoothed_Observed = NA,
                                                            DeltaM0_Smoothed	= NA,
                                                            DeltaM1_Smoothed	= NA,
                                                            DeltaM2_Smoothed = NA,
                                                            DeltaM0_Smoothed_Percentage	= NA,
                                                            DeltaM1_Smoothed_Percentage	= NA,
                                                            DeltaM2_Smoothed_Percentage  = NA,
                                                            M0_Area_Smoothed = NA, 
                                                            M1_Area_Smoothed = NA, 
                                                            M2_Area_Smoothed = NA,
                                                            M0_Height_Smoothed = NA, 
                                                            M1_Height_Smoothed = NA, 
                                                            M2_Height_Smoothed = NA,
                                                            DeltaM0_Smoothed_subB	= NA,
                                                            DeltaM1_Smoothed_subB	= NA,
                                                            DeltaM2_Smoothed_subB = NA,
                                                            DeltaM0_Smoothed_Percentage_subB	= NA,
                                                            DeltaM1_Smoothed_Percentage_subB	= NA,
                                                            DeltaM2_Smoothed_Percentage_subB  = NA,
                                                            M0_Area_Smoothed_subB = NA, 
                                                            M1_Area_Smoothed_subB = NA, 
                                                            M2_Area_Smoothed_subB = NA,
                                                            M0_Height_Smoothed_subB = NA, 
                                                            M1_Height_Smoothed_subB = NA, 
                                                            M2_Height_Smoothed_subB = NA,
                                                            M1_FullPeak_max_pea_spe_score = NA,
                                                            M1_Slices_max_pea_spe_score = NA,
                                                            M2_FullPeak_max_pea_spe_score = NA,
                                                            M2_Slices_max_pea_spe_score = NA,
                                                            Warning = warningmessage
                                                          )
                                                )
                    }
        } else # This "else" follows around line 398 to handle the specific issue when the formula is missing
          {
            warningmessage <- paste(warningmessage, "Missing formula to calculate predicted M0M1M2")
            extracted_table <- rbind(extracted_table, 
                                     data.frame(
                                                  Compound = compound,
                                                  Sample = sample_name,
                                                  Formula = formula,
                                                  Compound_Abundance = heightnumber,
                                                  RT_left = left_value/60,
                                                  RT_peak = middle_value/60,
                                                  RT_right = right_value/60,
                                                  Reference_mz = reference_mz,
                                                  Precursor_mz = mz_M0,
                                                  M1_mz = (mz_M1_low + mz_M1_high)/2,
                                                  M2_mz = (mz_M2_low + mz_M2_high)/2,    
                                                  Integrated_To_Area = 0,
                                                  Integrated_To_Area_onlyM0M1 = 0,
                                                  M0_Predicted = NA,
                                                  M1_Predicted = NA,
                                                  M2_Predicted = NA,
                                                  M0_Observed = NA,
                                                  M1_Observed = NA,
                                                  M2_Observed = NA,
                                                  DeltaM0	= NA,
                                                  DeltaM1	= NA,
                                                  DeltaM2 = NA,
                                                  DeltaM0_Percentage	= NA,
                                                  DeltaM1_Percentage	= NA,
                                                  DeltaM2_Percentage  = NA,      
                                                  Total_Ion = NA,
                                                  Total_Ion_Smoothed = NA,
                                                  M0_Area = NA, 
                                                  M1_Area = NA, 
                                                  M2_Area = NA,
                                                  M3_Area = NA,
                                                  M0_Height = NA, 
                                                  M1_Height = NA, 
                                                  M2_Height = NA,
                                                  M3_Height = NA,
                                                  M0_Baseline = NA,
                                                  M1_Baseline = NA,
                                                  M2_Baseline = NA,
                                                  M0_Noise = NA, 
                                                  M1_Noise = NA, 
                                                  M2_Noise = NA,
                                                  Mminus1_to_M0_Area = NA, 
                                                  M0_to_M1_Area = NA,
                                                  M1_to_M2_Area = NA,
                                                  M2_to_M3_Area = NA,
                                                  M3_to_M4_Area = NA,
                                                  M0_Observed_subB = NA,
                                                  M1_Observed_subB = NA,
                                                  M2_Observed_subB = NA,
                                                  DeltaM0_subB	= NA,
                                                  DeltaM1_subB	= NA,
                                                  DeltaM2_subB = NA,
                                                  DeltaM0_Percentage_subB	= NA,
                                                  DeltaM1_Percentage_subB	= NA,
                                                  DeltaM2_Percentage_subB  = NA,
                                                  M0_Area_subB = NA, 
                                                  M1_Area_subB = NA, 
                                                  M2_Area_subB = NA,
                                                  M0_Height_subB = NA, 
                                                  M1_Height_subB = NA, 
                                                  M2_Height_subB = NA,
                                                  M0_Smoothed_Observed = NA,
                                                  M1_Smoothed_Observed = NA,
                                                  M2_Smoothed_Observed = NA,
                                                  DeltaM0_Smoothed	= NA,
                                                  DeltaM1_Smoothed	= NA,
                                                  DeltaM2_Smoothed = NA,
                                                  DeltaM0_Smoothed_Percentage	= NA,
                                                  DeltaM1_Smoothed_Percentage	= NA,
                                                  DeltaM2_Smoothed_Percentage  = NA,
                                                  M0_Area_Smoothed = NA, 
                                                  M1_Area_Smoothed = NA, 
                                                  M2_Area_Smoothed = NA,
                                                  M0_Height_Smoothed = NA, 
                                                  M1_Height_Smoothed = NA, 
                                                  M2_Height_Smoothed = NA,
                                                  DeltaM0_Smoothed_subB	= NA,
                                                  DeltaM1_Smoothed_subB	= NA,
                                                  DeltaM2_Smoothed_subB = NA,
                                                  DeltaM0_Smoothed_Percentage_subB	= NA,
                                                  DeltaM1_Smoothed_Percentage_subB	= NA,
                                                  DeltaM2_Smoothed_Percentage_subB  = NA,
                                                  M0_Area_Smoothed_subB = NA, 
                                                  M1_Area_Smoothed_subB = NA, 
                                                  M2_Area_Smoothed_subB = NA,
                                                  M0_Height_Smoothed_subB = NA, 
                                                  M1_Height_Smoothed_subB = NA, 
                                                  M2_Height_Smoothed_subB = NA,
                                                  M1_FullPeak_max_pea_spe_score = NA,
                                                  M1_Slices_max_pea_spe_score = NA,
                                                  M2_FullPeak_max_pea_spe_score = NA,
                                                  M2_Slices_max_pea_spe_score = NA,
                                                  Warning = warningmessage
                                              )
                                    )
        }
        
        extracted_table2 <- rbind(extracted_table2, extracted_slices) # if there is no more issues, slices information will be updated to table2
        
      }, # This "}" follows the around line 110 trycatch function so that if any unexpected error happened, the program
        # will continue running other metabolites instead of stopping there      
        error=function(e) {
          message('An Error Occurred')
          cat("Error Compound: ", compound)
          warningmessage <<- paste(warningmessage, "An error happened")
          extracted_table <<- rbind(extracted_table, 
                                    data.frame(
                                                Compound = compound,
                                                Sample = sample_name,
                                                Formula = formula,
                                                Compound_Abundance = heightnumber,
                                                RT_left = left_value/60,
                                                RT_peak = middle_value/60,
                                                RT_right = right_value/60,
                                                Reference_mz = reference_mz,
                                                Precursor_mz = mz_M0,
                                                M1_mz = NA,
                                                M2_mz = NA,    
                                                Integrated_To_Area = 0,
                                                Integrated_To_Area_onlyM0M1 = 0,
                                                M0_Predicted = NA,
                                                M1_Predicted = NA,
                                                M2_Predicted = NA,
                                                M0_Observed = NA,
                                                M1_Observed = NA,
                                                M2_Observed = NA,
                                                DeltaM0	= NA,
                                                DeltaM1	= NA,
                                                DeltaM2 = NA,
                                                DeltaM0_Percentage	= NA,
                                                DeltaM1_Percentage	= NA,
                                                DeltaM2_Percentage  = NA,      
                                                Total_Ion = NA,
                                                Total_Ion_Smoothed = NA,
                                                M0_Area = NA, 
                                                M1_Area = NA, 
                                                M2_Area = NA,
                                                M3_Area = NA,
                                                M0_Height = NA, 
                                                M1_Height = NA, 
                                                M2_Height = NA,
                                                M3_Height = NA,
                                                M0_Baseline = NA,
                                                M1_Baseline = NA,
                                                M2_Baseline = NA,
                                                M0_Noise = NA, 
                                                M1_Noise = NA, 
                                                M2_Noise = NA,
                                                Mminus1_to_M0_Area = NA, 
                                                M0_to_M1_Area = NA,
                                                M1_to_M2_Area = NA,
                                                M2_to_M3_Area = NA,
                                                M3_to_M4_Area = NA,
                                                M0_Observed_subB = NA,
                                                M1_Observed_subB = NA,
                                                M2_Observed_subB = NA,
                                                DeltaM0_subB	= NA,
                                                DeltaM1_subB	= NA,
                                                DeltaM2_subB = NA,
                                                DeltaM0_Percentage_subB	= NA,
                                                DeltaM1_Percentage_subB	= NA,
                                                DeltaM2_Percentage_subB  = NA,
                                                M0_Area_subB = NA, 
                                                M1_Area_subB = NA, 
                                                M2_Area_subB = NA,
                                                M0_Height_subB = NA, 
                                                M1_Height_subB = NA, 
                                                M2_Height_subB = NA,
                                                M0_Smoothed_Observed = NA,
                                                M1_Smoothed_Observed = NA,
                                                M2_Smoothed_Observed = NA,
                                                DeltaM0_Smoothed	= NA,
                                                DeltaM1_Smoothed	= NA,
                                                DeltaM2_Smoothed = NA,
                                                DeltaM0_Smoothed_Percentage	= NA,
                                                DeltaM1_Smoothed_Percentage	= NA,
                                                DeltaM2_Smoothed_Percentage  = NA,
                                                M0_Area_Smoothed = NA, 
                                                M1_Area_Smoothed = NA, 
                                                M2_Area_Smoothed = NA,
                                                M0_Height_Smoothed = NA, 
                                                M1_Height_Smoothed = NA, 
                                                M2_Height_Smoothed = NA,
                                                DeltaM0_Smoothed_subB	= NA,
                                                DeltaM1_Smoothed_subB	= NA,
                                                DeltaM2_Smoothed_subB = NA,
                                                DeltaM0_Smoothed_Percentage_subB	= NA,
                                                DeltaM1_Smoothed_Percentage_subB	= NA,
                                                DeltaM2_Smoothed_Percentage_subB  = NA,
                                                M0_Area_Smoothed_subB = NA, 
                                                M1_Area_Smoothed_subB = NA, 
                                                M2_Area_Smoothed_subB = NA,
                                                M0_Height_Smoothed_subB = NA, 
                                                M1_Height_Smoothed_subB = NA, 
                                                M2_Height_Smoothed_subB = NA,
                                                M1_FullPeak_max_pea_spe_score = NA,
                                                M1_Slices_max_pea_spe_score = NA,
                                                M2_FullPeak_max_pea_spe_score = NA,
                                                M2_Slices_max_pea_spe_score = NA,,
                                                Warning = warningmessage
                                              )
                                  )
          message('Error handling of this compound finished')
      }
    ) # This is the end of the loop of one compound, then it will loop again to extract info of the next compound
 
  # After extracting M0M1M2 of all compounds and Create a new workbook
  wb <- createWorkbook()
  # Add sheets with data
  addWorksheet(wb, "Sheet1")
  writeData(wb, "Sheet1", extracted_table) # integrated data
  addWorksheet(wb, "Sheet2")
  writeData(wb, "Sheet2", extracted_table2) # individual slices data
  # Save the workbook
  saveWorkbook(wb, file = output_name, overwrite = TRUE)
} # This is the end of the loop of every sample file, then it will loop again to the next sample 1b file until it processed all sample files in the list