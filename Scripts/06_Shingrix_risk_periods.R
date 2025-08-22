################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 31/10/2024
# Version: R 4.3.0
# File name: 06_Shingrix_risk_periods.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID-filtered.rds
## filtered_med.parquet
## combined_med
# R scripts needed: 
## 00_Set_up_directory.R
## 02_Shingrix_filtered_patients.R
## 03_Shingrix_filtered_med.R
## 04_Shingrix_steroid_prescriptions.R
## 05_Shingrix_standard_prescriptions.R
# Data sets created: 
# Description of file: The following code establishes which patients are at risk over what time period, and will create two datasets to consider hierarchy of risks. One will be for the main analysis, where the category at start of follow-up dictates the overall which will have the category that is the max over the time under study.
# Actions performed: 
## Read in libraries
## Find risk under study for all patients
## Find overall risk times 
################################################################################



########################## Find start/end dates for categories #################
combined_med_immuno=read_parquet(raw_processed_med)
PATID_filtered <- read_rds(paste0(processed_data, "PATID_filtered.rds")) 
PATID_filtered <- subset(PATID_filtered, PATID_filtered$gender!=3)

# Add minimum and maximum follow-up to immuno records
PATID_fup = PATID_filtered %>% select(patid, min_cohort_entry, max_cohort_exit)
combined_med_immuno = left_join(combined_med_immuno, PATID_fup, by="patid")

################### 1A: haematologic malignancy ########################
immuno_cat_1A = combined_med_immuno %>% 
  filter(immuno_cat=="haema") %>%
  filter(obsdate >= as.Date("2019-09-01") & obsdate < as.Date("2023-08-31")) 
# filter to 2 years preceding potential study entry
nrow(immuno_cat_1A) # 86,122 records
length(unique(immuno_cat_1A$patid))  # 18,606 unique patids

# Summarise by patient
immuno_cat_1A_pat = immuno_cat_1A %>% 
  group_by(patid) %>%
  summarise(
    start_fup_1A = max(min(obsdate), min_cohort_entry),
    end_fup_1A = min(max(obsdate) %m+% months(24), max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_1A_pat = immuno_cat_1A_pat %>% 
  filter(start_fup_1A < end_fup_1A) 
nrow(immuno_cat_1A_pat) # 17,732

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_1A_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_1A_pat$patid, 5)
immuno_cat_1A %>% filter(patid %in% subsample_patids[1]) %>% 
  select(obsdate, min_cohort_entry, max_cohort_exit) %>% arrange(obsdate)
PATID_filtered %>% filter(patid %in% subsample_patids[1]) %>% select(start_fup_1A, end_fup_1A)
# Follow-up correctly assigned in all instances 


######################## 1B: Stem cell transplant ##############################
immuno_cat_1B = combined_med_immuno %>% 
  filter(immuno_cat=="stemcell") %>%
  filter(obsdate < as.Date("2023-08-31"))
nrow(immuno_cat_1B) # 1212 records
length(unique(immuno_cat_1B$patid))  # 667 unique patids

# Summarise by patient
immuno_cat_1B_pat = immuno_cat_1B %>% 
  group_by(patid) %>%
  summarise(
    start_fup_1B = max(min(obsdate), min_cohort_entry),
    end_fup_1B = max_cohort_exit[1] 
  ) 

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_1B_pat = immuno_cat_1B_pat %>% 
  filter(start_fup_1B < end_fup_1B) 
nrow(immuno_cat_1B_pat) # 667 

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_1B_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_1B_pat$patid, 5)
immuno_cat_1B %>% filter(patid %in% subsample_patids[5]) %>% 
  select(obsdate, min_cohort_entry, max_cohort_exit) %>% arrange(obsdate)
PATID_filtered %>% filter(patid %in% subsample_patids[5]) %>% 
  select(start_fup_1B, end_fup_1B)
# Follow-up correctly assigned in all instances 

# Assign combined follow-up for category 1
PATID_filtered = PATID_filtered %>%
  mutate(
    start_fup_1_combined = pmin(start_fup_1A, start_fup_1B, na.rm=TRUE),
    end_fup_1_combined = pmax(end_fup_1A, end_fup_1B, na.rm=TRUE),
  ) 

# Write output
write_rds(PATID_filtered, paste0(processed_data, "PATID_filtered_1.rds" ))

############################### 2: Solid Organ Transplant ###############################

immuno_cat_2 = combined_med_immuno %>% 
  filter(immuno_cat=="SOT") %>%
  filter(obsdate < as.Date("2023-08-31"))
nrow(immuno_cat_2) # 28,762 records
length(unique(immuno_cat_2$patid))  # 2,668 unique patids

# Summarise by patient
immuno_cat_2_pat = immuno_cat_2 %>% 
  group_by(patid) %>%
  summarise(
    start_fup_2 = max(min(obsdate), min_cohort_entry),
    end_fup_2 = max_cohort_exit[1] 
  ) 

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_2_pat = immuno_cat_2_pat %>% 
  filter(start_fup_2 < end_fup_2) 
nrow(immuno_cat_2_pat) # 2,665

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_2_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_2_pat$patid, 5)
immuno_cat_2 %>% filter(patid %in% subsample_patids[5]) %>% 
  select(obsdate, min_cohort_entry, max_cohort_exit) %>% arrange(obsdate)
PATID_filtered %>% filter(patid %in% subsample_patids[5]) %>% 
  select(start_fup_2, end_fup_2)
# Follow-up correctly assigned in all instances 

# Write output
write_rds(PATID_filtered, paste0(processed_data, "PATID_filtered_2.rds" ))

############################### 3: Cancer ###############################

immuno_cat_3 = combined_med_immuno %>% 
  filter(immuno_cat=="cancer") %>%
  filter(obsdate >= as.Date("2020-09-01") & obsdate < as.Date("2023-08-31")) 
# filter to 1 years preceding potential study entry
nrow(immuno_cat_3) # 253,925 records
length(unique(immuno_cat_3$patid))  # 67,837 unique patids

# Summarise by patient
immuno_cat_3_pat = immuno_cat_3 %>% 
  group_by(patid) %>%
  summarise(
    start_fup_3 = max(min(obsdate), min_cohort_entry),
    end_fup_3 = min(max(obsdate) %m+% months(12), max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_3_pat = immuno_cat_3_pat %>% 
  filter(start_fup_3 < end_fup_3) 
nrow(immuno_cat_3_pat) # 62,854

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_3_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_3_pat$patid, 5)
immuno_cat_3 %>% filter(patid %in% subsample_patids[1]) %>% 
  select(obsdate, min_cohort_entry, max_cohort_exit) %>% arrange(obsdate)
PATID_filtered %>% filter(patid %in% subsample_patids[1]) %>% 
  select(start_fup_3, end_fup_3)
# Follow-up correctly assigned in all instances 

############################### Load in CPRD Aurum drug parquet ################

# Generate a vector of file paths for x between 1 and 20
file_paths <- sprintf(drug_combined, 1:20)

# Read and combine all the parquet files
combined_drug <- file_paths %>%
  map_dfr(~ read_parquet(.x))
nrow(combined_drug) # 13,011,576 records

############################### 4A: Targeted ##################################

# Load in standard codelist
codelist_targeted <- read_delim(codelist_path_targeted,    
                                delim = "\t", escape_double = FALSE, 
                                col_types = cols(
                                  prodcodeid = col_character(),  
                                  dmdid = col_character()), 
                                trim_ws = TRUE)

Targeted_1 = combined_drug %>%
  filter(prodcodeid %in% codelist_targeted$prodcodeid)
nrow(Targeted_1) # 16,779 records

rm(combined_drug)

immuno_cat_4A = Targeted_1 %>% 
  filter(issuedate >= as.Date("2020-09-01") & issuedate < as.Date("2023-08-31")) 
# filter to 1 years preceding potential study entry
nrow(immuno_cat_4A) # 6,683 records
length(unique(immuno_cat_4A$patid))  # 5,499 unique patids

immuno_cat_4A=left_join(immuno_cat_4A, PATID_fup, by="patid")

# Summarise by patient
immuno_cat_4A_pat = immuno_cat_4A %>% 
  group_by(patid) %>%
  summarise(
    start_fup_4A = max(min(issuedate), min_cohort_entry),
    end_fup_4A = min(max(issuedate) %m+% months(12), max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_4A_pat = immuno_cat_4A_pat %>% 
  filter(start_fup_4A < end_fup_4A) 
nrow(immuno_cat_4A_pat) # 4,672

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_4A_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_4A_pat$patid, 5)
immuno_cat_4A %>% filter(patid %in% subsample_patids[1]) %>% 
  select(issuedate, min_cohort_entry, max_cohort_exit) %>% arrange(issuedate)
PATID_filtered %>% filter(patid %in% subsample_patids[1]) %>% 
  select(start_fup_4A, end_fup_4A)
# Follow-up correctly assigned in all instances

############################### 4B: Standard ###############################

immuno_cat_4B=read_parquet(paste0(parquet_processed, "/filtered_standard.parquet")) 
#706711
immuno_cat_4B$end_risk_add_3_months <- immuno_cat_4B$issuedate + lubridate::days(as.integer(immuno_cat_4B$calc_duration_days))+ months(3)

immuno_cat_4B = immuno_cat_4B %>% 
  filter(issuedate >= as.Date("2020-09-01") & issuedate < as.Date("2023-08-31")) 
# filter to 1 years preceding potential study entry
nrow(immuno_cat_4B) # 893,402 records
length(unique(immuno_cat_4B$patid))  # 27,410 unique patids

immuno_cat_4B=left_join(immuno_cat_4B, PATID_fup, by="patid")

# Summarise by patient
immuno_cat_4B_pat = immuno_cat_4B %>% 
  group_by(patid) %>%
  summarise(
    start_fup_4B = max(min(issuedate), min_cohort_entry),
    end_fup_4B = min(max(end_risk_add_3_months ), max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_4B_pat = immuno_cat_4B_pat %>% 
  filter(start_fup_4B < end_fup_4B) 
nrow(immuno_cat_4B_pat) # 15,959

PATID_filtered = left_join(PATID_filtered, immuno_cat_4B_pat, by="patid") 

############################### 4C: Steroid ###############################

immuno_cat_4C=read_parquet(paste0(parquet_processed, "/risk_steroids.parquet")) 
#102,818


immuno_cat_4C=left_join(immuno_cat_4C, PATID_fup, by="patid")

# Summarise by patient
immuno_cat_4C_pat = immuno_cat_4C %>% 
  group_by(patid) %>%
  summarise(
    start_fup_4C = max(min(risk_start), min_cohort_entry),
    end_fup_4C = min(max(end_risk_add_3_months), max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_4C_pat = immuno_cat_4C_pat %>% 
  filter(start_fup_4C < end_fup_4C) 
nrow(immuno_cat_4C_pat) # 53,423

PATID_filtered = left_join(PATID_filtered, immuno_cat_4C_pat, by="patid") 

# Assign combined follow-up for category 1
PATID_filtered = PATID_filtered %>%
  mutate(start_fup_4_combined = pmin(start_fup_4A, start_fup_4B,start_fup_4C, na.rm=TRUE),
  end_fup_4_combined = pmax(end_fup_4A, end_fup_4B,  end_fup_4C, na.rm=TRUE)) 

############################### 5A: HIV ###############################

immuno_cat_5A = combined_med_immuno %>% 
  filter(immuno_cat=="hiv") 
nrow(immuno_cat_5A) # 7621 records
length(unique(immuno_cat_5A$patid))  # 1,182 unique patids

# Summarise by patient
immuno_cat_5A_pat = immuno_cat_5A %>% 
  group_by(patid) %>%
  summarise(
    start_fup_5A = max(min(obsdate), min_cohort_entry),
    end_fup_5A = min(max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_5A_pat = immuno_cat_5A_pat %>% 
  filter(start_fup_5A < end_fup_5A) 
nrow(immuno_cat_5A_pat) # 1,176

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_5A_pat, by="patid") 

# Double check data for random subset of 5 patients
subsample_patids = sample(immuno_cat_5A_pat$patid, 5)
immuno_cat_5A %>% filter(patid %in% subsample_patids[1]) %>% 
  select(obsdate, min_cohort_entry, max_cohort_exit) %>% arrange(obsdate)
PATID_filtered %>% filter(patid %in% subsample_patids[1]) %>% select(start_fup_5A, end_fup_5A)
# Follow-up correctly assigned in all instances 

############################### 5B: Other conditions ###############################

immuno_cat_5B = combined_med_immuno %>% 
  filter(immuno_cat=="other_cond") 
nrow(immuno_cat_5B) # 661,560 records
length(unique(immuno_cat_5B$patid))  # 50,206 unique patids

# Summarise by patient
immuno_cat_5B_pat = immuno_cat_5B %>% 
  group_by(patid) %>%
  summarise(
    start_fup_5B = max(min(obsdate), min_cohort_entry),
    end_fup_5B = min(max_cohort_exit)
  ) %>%
  ungroup()

# Filter to individuals with at least 1 valid day of follow-up in category
immuno_cat_5B_pat = immuno_cat_5B_pat %>% 
  filter(start_fup_5B < end_fup_5B) 
nrow(immuno_cat_5B_pat) # 49,211

# Merge with patient line list
PATID_filtered = left_join(PATID_filtered, immuno_cat_5B_pat, by="patid") 

PATID_filtered = PATID_filtered %>%
  mutate(start_fup_5_combined = pmin(start_fup_5A, start_fup_5B, na.rm=TRUE),
         end_fup_5_combined = pmax(end_fup_5A, end_fup_5B,  na.rm=TRUE)  ) 

############################################## All together, with high/low and categories ####################

## Remove patients who have no date of immunosuppression 
PATID_filtered <- PATID_filtered %>%
  filter(!if_all(starts_with("start_fu"), is.na)) # 180,687

## Find the overall start and end date for each patient


PATID_filtered_2 <- PATID_filtered %>%
  rowwise() %>%
  mutate(min_entry = min(c_across(starts_with("start_fup")), na.rm = TRUE), 
         max_end = max(c_across(starts_with("end_fup")), na.rm = TRUE))

# Set code to high and low depending on what the patient enters as 
PATID_filtered_2$risk_condition="low"
PATID_filtered_2$risk_condition[PATID_filtered_2$min_entry==PATID_filtered_2$start_fup_1_combined |  
                                PATID_filtered_2$min_entry==PATID_filtered_2$start_fup_2 | 
                                PATID_filtered_2$min_entry==PATID_filtered_2$start_fup_3]="high"

PATID_filtered_3 <- PATID_filtered_2 %>%
  mutate(risk_category = case_when(
    !is.na(min_entry) & min_entry == start_fup_1_combined ~ 1,
    !is.na(min_entry) & min_entry == start_fup_2 ~ 2,
    !is.na(min_entry) & min_entry == start_fup_3 ~ 3,
    !is.na(min_entry) & min_entry == start_fup_4_combined ~ 4,
    !is.na(min_entry) & min_entry == start_fup_5_combined ~ 5,
    ))


## Use plot for overall idea of time under study

# hist(PATID_filtered_2$max_entry- PATID_filtered_2$min_entry)
# 
# ggplot(PATID_filtered_2[1:100,], aes(x = min_entry, xend = max_end, y = factor(patid), yend = factor(patid))) +
#   geom_segment(size = 2, color = "blue") +  # Create lines for each patient
#   scale_x_date(date_labels = "%b %d, %Y", date_breaks = "1 month") +  # Format the x-axis with dates
#   labs(x = "Study Dates", y = "Patient ID") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Angle the x-axis labels for readability

write_rds(PATID_filtered_3, paste0(processed_data, "Shingrix_Risks.rds" ))

