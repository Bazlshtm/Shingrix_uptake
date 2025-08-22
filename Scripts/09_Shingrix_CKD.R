################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 29/11/2024
# Version: R 4.3.0
# File name: 09_Shingrix_CKD.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID-filtered.rds
## filtered_med.parquet
## combined_med
## codelist_ckd_aurum
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
## df_CKD
# Description of file: This has been altered from Patrick's algorithm
# Actions:
## Link necessary files
## Reduce to valid records
## Apply algorithm to find egfr code which allows to find CKD risk thresholds
################################################################################

# This has been altered from Dr Patrick Bidulka's code from Stata to R

############################# Load in medcode data #############################
    
# Generate paths to 25 combined_med parquet files
file_paths <- sprintf(med_combined, 1:25)

# Read and combine all medical data
combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Load CKD-related medical code codelist (Aurum format)
codelist_ckd_aurum <- read_delim(
  paste0(codelists, "Medcodes/codelist_ckd_aurum.txt"),
  delim = "\t", 
  escape_double = FALSE, 
  col_types = cols(
    medcodeid = col_character(), 
    originalreadcode = col_character(), 
    snomedctconceptid = col_character(), 
    snomedctdescriptionid = col_character(), 
    release = col_skip(), 
    emiscodecategoryid = col_skip()
  ), 
  trim_ws = TRUE
)

# Load patient risk window and demographics
Shingrix_risktimes <- read_rds(paste0(processed_data, "Shingrix_Risks.rds")) %>%
  select(patid, gender, dob, min_entry, max_end)

# Join medical codes with CKD codelist and patient risk times
obsfile <- combined_med %>%
  inner_join(codelist_ckd_aurum, by = "medcodeid") %>%
  inner_join(Shingrix_risktimes, by = "patid")

############################# Clean and prepare CKD data #############################

# Extract and clean serum creatinine (SCr) data
df_ckd <- obsfile %>%
  mutate(SCr = as.numeric(value)) %>%
  rename(
    unit = numunitid,
    rangeFrom = numrangelow,
    rangeTo = numrangehigh
  ) %>%
  select(patid, SCr, unit, rangeFrom, rangeTo, obsdate, enterdate, gender, dob, min_entry, max_end) %>%
  filter(!is.na(SCr) & SCr != 0)

# Remove duplicates and implausible records
df_ckd <- df_ckd %>%
  distinct() %>%  # Drop exact duplicates
  filter(
    obsdate >= dob,                           # Must occur after birth
    obsdate <= max_end,                      # Must occur before end of follow-up
    obsdate >= as.Date("1942-01-01"),        # Drop very old values
    SCr >= 20 & SCr <= 3000,                 # Filter out implausible creatinine levels
    !is.na(obsdate)
  ) %>%
  mutate(age_at_mmt = as.numeric(difftime(as.Date(obsdate), dob, units = "days")) / 365.25)

# Keep only one observation per day per patient
df_ckd <- df_ckd %>%
  distinct(patid, obsdate, .keep_all = TRUE)

############################# Calculate eGFR using CKD-EPI formula #############################

# Convert SCr to mg/dL equivalent
df_ckd$SCr_adj <- (df_ckd$SCr * 0.95) / 88.4

# Compute 'min' term of eGFR formula based on sex
df_ckd$min <- ifelse(
  df_ckd$gender == 2,
  (df_ckd$SCr_adj / 0.7)^(-0.329),
  (df_ckd$SCr_adj / 0.9)^(-0.411)
)
df_ckd$min <- ifelse(df_ckd$min < 1, 1, df_ckd$min)

# Compute 'max' term of eGFR formula
df_ckd$max <- ifelse(df_ckd$gender == 2, df_ckd$SCr_adj / 0.7, df_ckd$SCr_adj / 0.9)
df_ckd$max <- df_ckd$max^(-1.209)
df_ckd$max[df_ckd$max > 1] <- 1

# Final CKD-EPI eGFR calculation
df_ckd$egfr <- df_ckd$min * df_ckd$max * 141
df_ckd$egfr <- df_ckd$egfr * (0.993^df_ckd$age_at_mmt)
df_ckd$egfr[df_ckd$gender == 2] <- df_ckd$egfr[df_ckd$gender == 2] * 1.018

############################# Categorise CKD Stage #############################

# Classify eGFR into CKD stages
df_ckd$egfr_cat <- cut(
  df_ckd$egfr, 
  breaks = c(0, 15, 30, 45, 60, 5000), 
  labels = c("stage 5", "stage 4", "stage 3b", "stage 3a", "no CKD")
)

# Create binary CKD indicator
df_ckd$CKD <- recode(
  df_ckd$egfr_cat, 
  "stage 5" = 1, 
  "stage 4" = 1, 
  "stage 3b" = 1, 
  "stage 3a" = 1, 
  "no CKD" = 0
)

# Retain only the most recent observation for each patient
df_ckd <- df_ckd %>%
  group_by(patid) %>%
  # drop patients if their ckd_date is later than their min_entry
  filter(!(obsdate > min_entry)) %>%
  # then keep only the row with the max obsdate
  filter(obsdate == max(obsdate)) %>%
  ungroup()

# Save final CKD classification (latest per patient)
write_rds(df_ckd[c("patid", "CKD")], paste0(processed_data, "CKD.rds"))