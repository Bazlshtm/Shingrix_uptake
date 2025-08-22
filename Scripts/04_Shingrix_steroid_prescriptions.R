################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 29/10/2024
# Version: R 4.3.0
# File name: 03_Shingrix_steroid_prescriptions.R
# Status: Checked
# CPRD version: March 2024
# Data sets used:
## Processed_data/PATID_filtered_immunocategories.rds
## drug_combined
## common_dosages_sep2024_path
# R scripts needed: 
## 00_Set_up_directory.R
## 02_Shingrix_filtered_patients.R
# Data sets created: None 
# Actions performed: 
## Read in steroid prescriptions
## Link to common dosages
## Derive daily prednisolone-equivalent dose and duration
## Deal with multiple patient records
## 
################################################################################


############################## 1 Read in data ##################################
# Read in patient data
PATID <- read_rds(paste0(processed_data, "PATID_filtered.rds")) 
nrow(PATID) # 1,523,450

# Generate a vector of file paths for x between 1 and 20
file_paths <- sprintf(drug_combined, 1:20)

# Read and combine all the parquet files
combined_drug <- file_paths %>%
  map_dfr(~ read_parquet(.x))
nrow(combined_drug) # 12,504,786 records
 
# Load in steroid codelist
codelist_steroids <- read_delim(codelist_path_steroids, 
                                             delim = "\t", escape_double = FALSE, 
                                             col_types = cols(prodcodeid = col_character(), 
                                             dmdid = col_character()), trim_ws = TRUE)

# Select subset of prescription records that relate to steroids 
steroids_1 = combined_drug %>%
  filter(prodcodeid %in% codelist_steroids$prodcodeid)
nrow(steroids_1) # 7,088,174 records

# Limit to records in the year preceding study start date
steroids_2 = steroids_1 %>%
  filter(issuedate >= as.Date("2020-09-01") & issuedate<"2023-09-01") 
nrow(steroids_2) # 1,763,768 records

# Merge to add codelist info
steroids_2 = left_join(steroids_2, codelist_steroids, by="prodcodeid")

# Limit to records from potentially eligible patients
steroids_3 = steroids_2 %>%
  filter(patid %in% PATID$patid)
nrow(steroids_3) # 1,633,361 records

# Drop injections 
steroids_3 = steroids_3 %>%
  filter(
    !(steroids_3$quantunitid==95 | steroids_3$quantunitid==62 | steroids_3$quantunitid==3)
  ) 
# 1597708/1633361 (2.2%)

############################ 2 Link common dosage ##############################

# Read in common dosages for Sep 2024 build and subset to steroids
common_dosages_sep2024 <- read_delim(common_dosages_sep2024_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(dosageid %in% steroids_3$dosageid)
nrow(common_dosages_sep2024) #3,094

# Merge common dosages with prescriptions
steroids_4 = left_join(steroids_3, common_dosages_sep2024, by="dosageid")

# Remove packs prescribed 'as needed' or 'as required' 
check = steroids_4 %>% 
  filter(grepl("REQUIRED", dosage_text)| grepl("NEEDED", dosage_text)) %>% 
  count(dosage_text, sort=T)
# All should be removed - cannot reliably determine daily dose taken
steroids_4 <- steroids_4 %>% 
  filter(!grepl("REQUIRED", dosage_text) & !grepl("NEEDED", dosage_text)) # 1000 removed

# Recode numerics
steroids_4$duration = as.numeric(steroids_4$duration)
steroids_4$quantity = as.numeric(steroids_4$quantity)

# Remove unused data
rm(list=ls(pattern="common_dosages_m"))
rm(combined_drug, steroids_1, steroids_2, steroids_3)

# Export drug substances requiring prednisolone equivalent dose
ped_list = steroids_4 %>% count(drugsubstancename)
#write.csv(ped_list, paste0(processed_data,"PED_unassigned.csv"))

# Read in assigned prednisolone equivalent dose
ped_assigned = read.csv(paste0(processed_data,"PED_assigned.csv")) %>%
  select(drugsubstancename, pred_equiv_dose)
ped_assigned$pred_equiv_dose = as.numeric(ped_assigned$pred_equiv_dose)

# Merge with steroids prescriptions and filter non-active substances
steroids_5 = left_join(steroids_4, ped_assigned, by = "drugsubstancename")

# Filter prescriptions with no activity
steroids_6 = steroids_5 %>%
  filter(pred_equiv_dose>0)
nrow(steroids_6) # 1,477,241 records

# Create final variables to populate for each prescription
steroids_6 = steroids_6 %>%
  mutate(
    calc_duration_days = NA,
    calc_strength_per_unit_mg = NA,
    calc_units_per_day = NA
  )

##### Assign strength per unit in mg
# Turn off scientific notation
options(scipen=999)

# Convert all strengths to mg
steroids_6 = steroids_6 %>%
  mutate(
    calc_strength_per_unit_mg = case_when(
      str_detect(substancestrength, "mg$") ~ as.numeric(str_extract(substancestrength, "\\d+(\\.\\d+)?")) * 1,
      str_detect(substancestrength, "microgram$") ~ as.numeric(str_extract(substancestrength, "\\d+(\\.\\d+)?")) / 1000,
      str_detect(substancestrength, "microgram/1\\.000ml$") ~ as.numeric(str_extract(substancestrength, "\\d+(\\.\\d+)?")) / 1000,
      str_detect(substancestrength, "gram$") ~ as.numeric(str_extract(substancestrength, "\\d+(\\.\\d+)?")) * 1000,
      str_detect(substancestrength, "mg/1\\.000ml$") ~ as.numeric(str_extract(substancestrength, "\\d+(\\.\\d+)?")),
      TRUE ~ NA_real_  # For any other format, assign NA
    ))


##### Clean daily_dose to ensure that always expressed as N units per day
# Check potential exclusions
check = steroids_6 %>%
  filter(
    (grepl("MG", dosage_text) | grepl("ML", dosage_text) | grepl("MCG", dosage_text)) & daily_dose>1) %>%
  count(duration, quantity, substancestrength, dosage_text, daily_dose, dose_duration, sort=T)

# All doses in this format are expressed as dosage rather than N units, but possible to derive based on substancestrength
steroids_6 = steroids_6 %>%
  mutate(
    daily_dose = if_else(
      ((grepl("MG", dosage_text) & !grepl("tablet", dosage_text, ignore.case = T)) |  grepl("MCG", dosage_text)) & daily_dose>1, daily_dose/calc_strength_per_unit_mg, daily_dose
    )
  )

# However, ML is problematic and cannot be guaranteed, so remove
steroids_6 = steroids_6 %>%
  mutate(
    daily_dose = if_else(
      ( grepl("ML", dosage_text) | (grepl("MG", dosage_text) & grepl("tablet", dosage_text, ignore.case = T)) ) & daily_dose>1, NA, daily_dose
    )
  )

# Check derived values
check = steroids_6 %>%
  filter(
    (grepl("MG", dosage_text) | grepl("ML", dosage_text) | grepl("MCG", dosage_text)) ) %>%
  count(duration, quantity, substancestrength, calc_strength_per_unit_mg, dosage_text, daily_dose, dose_duration, sort=T)

# Check remainder of daily doses
check = steroids_6 %>% filter(daily_dose>0) %>% count(dosage_text, daily_dose, sort=T) # 1477241




############################### 3 Assign duration ##############################

# Use duration from look-up where available
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(dose_duration>0, dose_duration, calc_duration_days))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 2.5% assigned

# If daily_dose available from look-up and quantity in prescription, set duration as quantity/daily_dose
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & daily_dose>0 & quantity>0, quantity/daily_dose, calc_duration_days))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 34.0% assigned

# If quantity <=1, set duration and quantity to 1
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & quantity>0 & quantity<=1, 1, calc_duration_days),
         quantity = if_else(is.na(calc_duration_days) & quantity>0 & quantity<=1, 1, quantity))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 34.0% assigned

# Use recorded prescription duration for remainder
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & duration>0, duration, calc_duration_days))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 98.3% assigned

# Check median duration
median(steroids_6$calc_duration_days, na.rm=T)

# Check profiles if duration 0
check = steroids_6 %>% 
  filter(duration==0 & is.na(calc_duration_days)) %>% 
  count(duration, quantity, dosage_text, calc_strength_per_unit_mg, daily_dose, sort=T)

# Quantities and strengths predominantly compatible with a duration of 28 days, so reassign
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & duration==0, 28, calc_duration_days))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 100% assigned

# Check profiles if duration >=365
check = steroids_6 %>% 
  filter(calc_duration_days >= 365) %>% 
  count(duration, quantity, dosage_text, calc_strength_per_unit_mg, daily_dose, sort=T) #34

# Quantities and strengths predominantly compatible with a duration of 28 days, so reassign
steroids_6 = steroids_6 %>%
  mutate(calc_duration_days = if_else(calc_duration_days >= 365, 28, calc_duration_days))
sum(!is.na(steroids_6$calc_duration_days))/nrow(steroids_6)*100 # 100% assigned


############################## Assign units per day ############################
# Use daily_dose from look-up where available
steroids_6 = steroids_6 %>%
  mutate(calc_units_per_day = if_else(daily_dose>0, daily_dose, calc_units_per_day))
sum(!is.na(steroids_6$calc_units_per_day))/nrow(steroids_6)*100 # 34.0% assigned

# If quantity>0, use quantity/calc_duration_days 
steroids_6 = steroids_6 %>%
  mutate(calc_units_per_day = if_else(is.na(calc_units_per_day) & !is.na(calc_duration_days) & quantity>0, quantity/calc_duration_days, calc_units_per_day))
sum(!is.na(steroids_6$calc_units_per_day))/nrow(steroids_6)*100 # 99.9% assigned

# Remove incomplete records
steroids_6 <- steroids_6 %>%
  mutate(
    steroid_complete = if_else(
      !is.na(pred_equiv_dose) & 
        !is.na(calc_duration_days) & 
        !is.na(calc_strength_per_unit_mg) & 
        !is.na(calc_units_per_day), 1, 0))  

sum(steroids_6$steroid_complete == 1) / nrow(steroids_6) * 100 #99.9

# Limit to prescriptions with complete data
steroids_7 = steroids_6 %>% filter(steroid_complete==1)
nrow(steroids_7) # 1,476,031

# Calculate daily prednisolone equivalent dose for valid prescriptions
steroids_7 = steroids_7 %>% mutate(daily_pred_equiv_dose = pred_equiv_dose*calc_strength_per_unit_mg*calc_units_per_day)

############################# Multiple patient records #########################

# 1. Count duplicates to assess the size of the problem
# Check duplicate entries based on `patid` and `eventdate`
steroids_8 <- steroids_7 %>%
  mutate(dupe = duplicated(select(., patid, issuedate)))
sum(steroids_8$dupe)*100/nrow(steroids_8) # 8.8%

# Drop the `dupe` column after analysis
steroids_8 <- steroids_8 %>% select(-dupe)

# 2. Count duplicate prescriptions with same `patid`, `issuedate`, and `ped`
steroids_8 <- steroids_8 %>%
  mutate(dupe = duplicated(select(., patid, issuedate, daily_pred_equiv_dose)))
sum(steroids_8$dupe)*100/nrow(steroids_8) 

steroids_8 <- steroids_8 %>% select(-dupe)

# Two records on the same day with the same strength
# Combine and increase prescription duration
steroids_8 <- steroids_8 %>%
  arrange(patid,issuedate, daily_pred_equiv_dose, desc(calc_duration_days)) %>% 
  group_by(patid, issuedate, daily_pred_equiv_dose) %>%
  mutate(count = row_number(), # which number for patid is this (i.e. first, second)
         COUNT = n(), # How many duplicates of patid
         totaldays = sum(calc_duration_days)) %>%
  ungroup()

# Cap `totaldays` at 365 if COUNT > 1 and `totaldays > 365` (not needed)
steroids_8 <- steroids_8 %>%
  mutate(totaldays = ifelse(COUNT > 1 & totaldays > 365, 365, totaldays))

# Drop duplicate prescriptions on the same day
steroids_8  <- steroids_8  %>%
  distinct(patid, issuedate, daily_pred_equiv_dose, .keep_all = TRUE)

# Update duration with `totaldays`
steroids_8  <- steroids_8  %>%
  mutate(calc_duration_days = ifelse(calc_duration_days != totaldays, totaldays, calc_duration_days)) %>%
  select(-totaldays, -count, -COUNT)

# Two records on the same day with DIFFERENT strengths
# Combine and increase prescription dose

steroids_9 <- steroids_8  %>%
  group_by(patid, issuedate) %>%
  summarise(
    calc_duration_days = max(calc_duration_days, na.rm = TRUE), 
    daily_pred_equiv_dose = sum(daily_pred_equiv_dose, na.rm = TRUE),  
    .groups = "drop"  
  )

steroids_9 <- steroids_9 %>%
  mutate(enddate = issuedate + (calc_duration_days-1))

# Double check data for random subset of 5 patients
subsample_patids = sample(steroids_9$patid, 5)
subset_7 <- steroids_7 %>% filter(patid %in% subsample_patids[1])
subset_9 <- steroids_9 %>% filter(patid %in% subsample_patids[1])

# If patient's have overlapping prescription, later prescriptions should 
# override earlier ones in terms of daily prednisolone equivalent dose; the 
# duration will run until the end of the later prescriptions
steroids_10 <- steroids_9 %>%
  arrange(patid, issuedate) %>%  # Sort by patient ID and issue date
  group_by(patid) %>%  # Group by patient to compare within patient
  mutate(
    next_dose = lead(daily_pred_equiv_dose),  # Get the dose of the next prescription
    next_issuedate = lead(issuedate), # Get the issue date of next prescription
    next_enddate= lead(enddate)
    ) %>%
  # Adjust duration if prescript overlap
  mutate(enddate = if_else(!is.na(next_issuedate) & next_issuedate <= enddate & next_enddate>=enddate, next_issuedate-1, enddate)
  )

# Investigate distribution
quantile(steroids_10$daily_pred_equiv_dose, prob=c(0,0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999,1))

# Modify the dataset
steroids_11 <- steroids_10 %>%
  filter(daily_pred_equiv_dose  <= 100) 

## Remove patients whose whole dose is under lowest limit per day
filtered_data <- steroids_11  %>%
  group_by(patid) %>%
  mutate(sum_days = sum(calc_duration_days)) %>%
  filter(sum_days >= 7) %>%  # Keep only groups with a dose 
  ungroup()
length(unique(filtered_data$patid ))/length(unique(steroids_11$patid )) # now down to 83.1% retained

# Save processed steroid prescriptions to secure drive
write_parquet(filtered_data, paste0(parquet_processed, "filtered_steroids.parquet"))



######################## Find risk periods depending on intake #################

################################################################################
# The following code takes each patient and finds their immunosuppression times. 
# This does this for # each of the three criteria (long, moderate, short dose). 
# It takes the time the patient is under study # and finds the dose per day. 
# Doses are then flagged (when they meet a sufficient threshold in lookback  
# period. These TRUE times are then condensed when they are still taking the 
# amount that puts then at risk (dose_mg >= condition$mg_level). The earliest 
# and latest of these dates are taken as drug starting and ending. Three months 
# is then added as the end. 
################################################################################


get_risk_periods <- function(df) {
  
  # Define risk conditions
  risk_conditions <- list(
    list(duration_days = 28, window_days = 90, mg_level = 10),
    list(duration_days = 10, window_days = 30, mg_level = 20),
    list(duration_days = 7, window_days = 7, mg_level = 40)
  )
  
  # Initialize dataframe to store risk periods
  risk_periods <- data.frame()
  
  # Loop over each patient
  unique_patids <- unique(df$patid)
  
  for (i in 1:length(unique_patids)) {
    patient_data <- df %>% filter(patid == unique_patids[i])
    
    # Progress tracker
    print(paste0("Progress: ", round((i/length(unique_patids))*100, 2), "%"))
    
    # Define treatment period sequence
    treatment_dates <- seq(min(patient_data$issuedate), 
                           max(patient_data$enddate), 
                           by = "day")
    
    # Calculate cumulative daily doses
    daily_doses <- sapply(treatment_dates, function(date) {
      # Extract the relevant doses based on the issuedate and enddate conditions
      doses <- patient_data$daily_pred_equiv_dose[patient_data$issuedate <= date & patient_data$enddate >= date]
      
      # If there are no doses, return 0. Otherwise, return the maximum dose
      if (length(doses) == 0) {
        return(0)
      } else {
        return(max(doses))
      }
    })
    
    daily_dose_df <- data.frame(date = treatment_dates, dose_mg = daily_doses)
    
    # Loop over each risk condition
    for (condition_id in 1:length(risk_conditions)) {
      condition <- risk_conditions[[condition_id]]
      
      # Identify dates that meet risk condition's dose threshold
      daily_dose_df <- daily_dose_df %>%
        mutate(at_risk_flag = ifelse(dose_mg >= condition$mg_level, 1, 0)) %>%
        mutate(risk_count = slide_int(at_risk_flag, sum, .before = condition$window_days, .complete = FALSE),
               at_risk = risk_count >= condition$duration_days)
      
      # Check and store risk periods
      if (any(daily_dose_df$at_risk)) {
        
        at_risk_dates <- daily_dose_df %>%
          filter(at_risk == TRUE & dose_mg >= condition$mg_level)
        
        # Find the first and last dates where at_risk is TRUE
        first_risk_date <- min(at_risk_dates$date, na.rm = TRUE)
        last_risk_date <- max(at_risk_dates$date, na.rm = TRUE)
        
        # Append to final dataframe
        risk_periods <- rbind(risk_periods, data.frame(
          patid = unique_patids[i],                    
          risk_condition = condition_id,               
          risk_start = first_risk_date, 
          risk_end = last_risk_date
        ))
      }
    }
  }
  
  return(risk_periods)
}

# Apply the function
final_risk_periods <- get_risk_periods(filtered_data)
final_risk_periods$end_risk_add_3_months=final_risk_periods$risk_end %m+% months(3)   

write_parquet(final_risk_periods, paste0(parquet_processed, "/risk_steroids.parquet"))

# Double check data for random subset of 5 patients
subsample_patids = sample(final_risk_periods$patid, 5)
subset2=final_risk_periods %>% filter(patid %in% subsample_patids[1]) 
subset1=steroids_9 %>% filter(patid %in% subsample_patids[1])
# Follow-up correctly assigned in all instances 



