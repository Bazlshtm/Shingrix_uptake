# Shingrix_uptake

The following codelists and code were used in the pulbication of "Disparities in uptake of Shingrix® vaccine in immunosuppressed individuals in England: a population-based cohort study"

The following files are contained:
- Scripts - R scripts used, with numerical order in which they were run
- Codelists - codelists for the study compatible with CPRD, under build March 2024. 

## DEPENDENCIES
- **R version:** ≥ 4.2.0  
- **Packages needed:** Listed in `00_Shingrix_set_up_directory_GITHUB.R`

The preprint can be found here: (TBC)

## DATA SOURCES
- Source: CPRD, March 2024 build

---

## USAGE INSTRUCTIONS
The scripts are designed to be run **sequentially**, with script 13 being the **output** script. 

---

## SUMMARY SCRIPTS

- **00_Shingrix_Set_up_directory_GITHUB** – Files/directories of the data used in this study, personalized to the user.  

- **01_Shingrix_parquet_transformation.R** – Transforms raw CPRD files flagged with relevant immunosuppression codes and returns all those with relevant codes only.  

- **02_Shingrix_filtered_patients** – Subsets patients into those who are CPRD acceptable and within the years of interest.  

- **03_Shingrix_filtered_med.R** – Filters all patients who have relevant medical codes.  

- **04_Shingrix_steroid_prescriptions.R** – Processes steroid prescription data for eligible patients, converting doses to daily prednisolone-equivalents, standardizing durations, and resolving duplicate or overlapping prescriptions. Identifies periods of immunosuppression risk based on predefined dose–duration thresholds, and saves both the cleaned prescription dataset and calculated risk periods.  

- **05_Shingrix_standard_prescriptions.R** – Processes standard prescription records for eligible patients, calculates standardized daily doses and durations, spreads multiple injection doses evenly across their course, and saves the cleaned dataset.  

- **06_Shingrix_risk_periods.R** – Finds the times at risk and allocates patients to the relevant risk categories.  

- **07_Shingrix_BMI.R** – Calculates the BMI for all patients.  

- **08_Shingrix_Ethnicity.R** – Finds the ethnicity of relevant patients.  

- **09_Shingrix_CKD.R** – Finds chronic kidney disease for relevant patients.  

- **10_Shingrix_covariates.R** – Assigns relevant covariates of interest.  

- **11_Shingrix_zoster_uptake.R** – Finds the first and second vaccination date for each person and removes people if they have vaccination prior to entry.  

- **12_Shingrix_flu.R** – Finds flu records for relevant patients.  

- **13_Output_code.R** – Outputs all code seen in the article.  

---

## OUTPUT

### Main article:
- **Figure 1:** Sample size prior to immunosuppression status is derived at the end of script 02. The sample size of those who are immunosuppressed over the study period is derived at the end of script 06. The final sample size is derived at the end of script 11.  
- **Table 2:** Derived in lines 1–171 of `13_Output_code.R`.  
- **Table 2 (summary):** Derived in lines 179–214 of `13_Output_code.R`.  
- **Figure 2:** Produced using the ORs derived in lines 179–214 of `13_Output_code.R`. The plotting code can be found in `14_Output_visualisation`.  

**Flu uptake:** Additional code in lines 280–304 of `13_Output_code.R` finds the time between vaccinations.  

### Supplementary material:

- **Supplementary Table S3** – Descriptive coding for `05_Shingrix_standard_prescriptions.R`  
- **Supplementary Table S3** – Descriptive coding for `04_Shingrix_steroid_prescriptions.R`  
- **Supplementary Table S4** – Describes the hierarchical modelling strategy (lines 179–185 of `13_Output_code.R`).  

- **Supplementary Table S5/8:**  
    - *Static cohort* – Cohort reduced to those in the study on 1 Sept 2021. Command run from console, analysis repeated from `13_Output_code.R`.  
    - *Time-updated covariates* – Run through scripts 1-12. Script `10b_Shingrix_covariates_timevarying.R` can be used in conjection with `07b_Shingrix_BMI_timevary` and `09b_Shingrix_CKD_timevary.R`, then script 13 re-run with the updated dataset from 10b.  
    - *Restricted to 13 months* – See line 45 in `13_Output_code.R` (restricts to this follow-up time).  
    - *Restricted to age 70 only* – Cohort reduced to those aged 70 on entry. Command run from console, analysis repeated from `13_Output_code.R`.  
    - *Sixteen categories of ethnicity* – Adjusted in `10_Shingrix_covariates.R` and `13_Output_code.R`.  
    - *Five immunosuppression risk categories* – Formulated in `06_Shingrix_risk_periods.R` and adjusted in `13_Output_code.R`.  

 - **Supplementary Table S11** – Lines 225–274 of `13_Output_cority ethnic backgrounds. 
