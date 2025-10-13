# Isotope_Ratio_Correction
Isotope_Ratio_Correction are R codes to train correction models that correct isotopic ratio measurement bias in Orbitrap mass spectrometers based on the input of all reported bias causing factors

# Package Requirements and Software Versions
All required packages for each “.R” coding file can be found at their headers. Most of the packages can be easily installed using the install.packages() function, while a few require the installation of BiocManager first, followed by BiocManager::install(). The Step 1c codes also require Python installation for R usage and can be executed by running the included comments once.
The R version used is 4.3.3, and at least some earlier 4.x.x versions should also work.
The MS-DIAL metabolite annotation software version used is 4.8. We have observed that other versions of MS-DIAL may produce output files with different formats and naming conventions, which could affect step 1's slice-level data retrieval. In addition, the provided annotation parameters and metabolite database format are specifically suited for version 4.8, so we recommend sticking to MS-DIAL version 4.8.

# Index
This respiratory contain several folders: 0_All_Codes, 1_raw_and_MS-DIAL, 2a_Combine_1c_Results, 2b_Preprocess_Input_Data, 2cEN_Bias_Correction_by_ElasticNet, 2cRF_Bias_Correction_by_RandomForest, and example_calculations.

# Step by Step Tutorials
The folder 0_All_Codes contains all “.R” coding files from the later folders as a backup and for inspection.
To replicate our work, proceed through the folders in numerical order and run all “.R” coding files according to their respective numbers. For example: 1a → 1b → 1c → 2a → 2b → (2cEN) → 2cRF. Each “.R” coding file’s header also includes step-specific tutorials describing the required input files and the parameters that may need inspection and modification.
The folder “1_raw_and_MS-DIAL” contains steps for extracting all slice-level metabolite isotope ratio information and measurement bias related factors.
Due to GitHub's size limitations, this folder does not include the full MS-DIAL project folder or the “.mzML” raw metabolomics result files to run the codes. (Different company’s’ raw metabolomics result files can all be converted to the open-source “.mzML” format by using software ProteoWizard's Msconvert) However, these raw data will be available through the Metabolomics Workbench or upon direct request. Still, all samples’ middle-step files and final slice-level info data are included. This folder also includes our MS-DIAL annotation parameters and used metabolite annotation database file.  MS-DIAL metabolite annotation tutorial could be found at https://systemsomicslab.github.io/compms/msdial/main.html
Once step 1 data retrieval is completed, the output files need to be copied into the folder 2a_Combine_1c_Results to combine all results into a single file.
The output file from step 2a should then be copied to the folder 2b_Preprocess_Input_Data, where we will inspect which portions of the slice-level data contain accurate isotopic ratios suitable for training the initial correction model, and which portions contain inaccurate isotopic ratios.
These inaccurate portions can still be used to train later iteration models, as the initial model can be applied to correct this inaccurate data and assess whether the corrected values follow the true isotopic ratios.
After all data has been clustered and pre-processed, the step 2b output files need to be copied to both 2cEN_Bias_Correction_by_ElasticNet and 2cRF_Bias_Correction_by_RandomForest folders.
We tested both a linear Elastic Net model and a nonlinear Random Forest model to correct isotopic measurement bias. The Elastic Net method is faster but performs worse than the Random Forest method. Therefore, Step 2cEN can be skipped, and only Step 2cRF (Random Forest training) needs to be performed. A full Random Forest training session may take up to 12 hours and requires a minimum of 16 GB of RAM. However, once the final model is trained, it can be saved and used to correct future data without retraining. Due to GitHub size limitations, we include only the final model generated using a randomness seed of 2024 as an example.
Finally, folder example_calculations contains examples showing how the fully trained Random Forest correction model could correct isotopic ratio measurement bias and a real application in MIDA.

