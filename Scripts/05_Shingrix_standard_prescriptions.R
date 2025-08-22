################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 29/10/2024
# Version: R 4.3.0
# File name: 05_Shingrix_standard_prescriptions.R
# Status: Checked
# CPRD version: March 2024
# Data sets used:
## Processed_data/PATID_filtered_immunocategories.rds
## combined_drug
## common_dosages_sep2024_path
# R scripts needed: 
## 00_Set_up_directory.R
## 02_Shingrix_filtered_patients.R
# Data sets created: "filtered_standard.parquet"
# Actions performed: 
## Read in standard prescriptions
## Link to common dosages
## Derive daily dose and duration
## Sort injections so they are spread over the duration calculated
## Save
################################################################################

# Read in patient data
PATID <- read_rds(paste0(processed_data, "PATID_filtered.rds")) 
nrow(PATID) # 1,523,450

# Generate a vector of file paths for x between 1 and 20
file_paths <- sprintf(drug_combined, 1:20)

# Read and combine all the parquet files
combined_drug <- file_paths %>%
  map_dfr(~ read_parquet(.x))
nrow(combined_drug) # 12,504,786 records
 
# Load in standard codelist
codelist_standard <- read_delim(codelist_path_standard,    
                                             delim = "\t", escape_double = FALSE, 
                                             col_types = cols(
                                               prodcodeid = col_character(),  
                                               dmdid = col_character()), 
                                               trim_ws = TRUE)

# Select subset of prescription records that relate to standard 
standard_1 = combined_drug %>%
  filter(prodcodeid %in% codelist_standard$prodcodeid)
nrow(standard_1) # 3,612,392 records

# Limit to records in the year preceding study start date
standard_2 = standard_1 %>%
  filter(issuedate >= as.Date("2020-09-01"))
nrow(standard_2) # 767,114 records

# Merge to add codelist metadata
standard_2 = left_join(standard_2, codelist_standard, by="prodcodeid")

# Limit to records from potentially eligible patients
standard_3 = standard_2 %>%
  filter(patid %in% PATID$patid)
nrow(standard_3) # 730,487 records

# Read in common dosages for Sep 2024 build and subset to standard
common_dosages_sep2024 <- read_delim(common_dosages_sep2024_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(dosageid %in% standard_3$dosageid)
nrow(common_dosages_sep2024) #1,056

# Merge common dosages with prescriptions
standard_4 = left_join(standard_3, common_dosages_sep2024, by="dosageid")

# Remove 'Required/needed packs' 
standard_5 = standard_4 %>%
  filter(
    !(grepl("REQUIRED", dosage_text)| grepl("NEEDED", dosage_text))
  ) 
nrow(standard_5) # 73,0478 (9 removed)

# Recode numerics
standard_5$duration = as.numeric(standard_5$duration)
standard_5$quantity = as.numeric(standard_5$quantity)

# Remove unused data
rm(combined_drug)
rm(standard_1, standard_2, standard_3)

# Create final variables to populate for each prescription
standard_6 = standard_5 %>%
  mutate(
    calc_duration_days = NA,
    calc_units_per_day = NA
  )

##### Clean daily_dose to ensure that always expressed as N units per day
# Check potential exclusions
check = standard_6 %>%
  filter(
    (grepl("MG", dosage_text) | grepl("ML", dosage_text) | grepl("MCG", dosage_text)) & daily_dose>1 
  ) %>%
  count(dosage_text, daily_dose, sort=T) 

# ML is not efficient, MG could be is, no MCG
# The number of doses of this format is 6988, which is 0.4% of all doses and therefore not worth making mistakes over, so not uding daily dose
standard_6 = standard_6 %>%
  mutate(
    daily_dose = if_else(
      (grepl("MG", dosage_text) | grepl("ML", dosage_text) | grepl("MCG", dosage_text)) , NA, daily_dose
    )
  )

# Check remainder of daily doses
check = standard_6 %>% filter(daily_dose>0) %>% count(dosage_text, daily_dose, sort=T)
# All seem to reflect N units per day

##### Assign duration
# Use duration from look-up where available
standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(dose_duration>0, dose_duration, calc_duration_days))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # <0.1% assigned

# If daily_dose available from look-up and quantity in prescription, set duration as quantity/daily_dose
standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & daily_dose>0 & quantity>0, quantity/daily_dose, calc_duration_days))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # 13.7% assigned

# If quantity <=1, set duration and quantity to 1
standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & quantity>0 & quantity<=1, 1, calc_duration_days),
         quantity = if_else(is.na(calc_duration_days) & quantity>0 & quantity<=1, 1, quantity))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # 14.1% assigned

# Use recorded prescription duration for remainder
standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & duration>0, duration, calc_duration_days))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # 96.8% assigned

# Check median duration
median(standard_6$calc_duration_days, na.rm=T)

# Check profiles if duration 0
check = standard_6 %>% 
  filter(duration==0 & is.na(calc_duration_days)) %>% 
  count(duration, quantity, dosage_text,  daily_dose, sort=T)

standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(is.na(calc_duration_days) & duration==0, 28, calc_duration_days))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # 100% assigned

check = standard_6 %>% 
  filter(calc_duration_days >= 365) %>% 
  count(duration, quantity, dosage_text,  daily_dose, sort=T) # 

# Quantities and strengths predominantly compatible with a duration of 28 days, so reassign
standard_6 = standard_6 %>%
  mutate(calc_duration_days = if_else(calc_duration_days >= 365, 28, calc_duration_days))
sum(!is.na(standard_6$calc_duration_days))/nrow(standard_6)*100 # 100% assigned

##### Assign strength per unit in mg
# Turn off scientific notation
options(scipen=999)

# If quantity>0, use quantity/calc_duration_days 

# Remove incomplete records
standard_6 <- standard_6 %>%
  mutate(
    steroid_complete = if_else(
        !is.na(calc_duration_days), 1, 0))  

sum(standard_6$steroid_complete == 1) / nrow(standard_6) * 100 #100

# Limit to prescriptions with complete data
standard_6 = standard_6 %>% filter(steroid_complete==1)
nrow(standard_6) # 730478



############ Change multiple injections spaced between duration ################

standard_notinject <- standard_6 %>%
  filter((quantunitid!=95 & quantunitid!=62 & quantunitid!=3) | quantity<=1 | calc_duration_days==0)

standard_inject <- standard_6 %>%
  filter((quantunitid==95 | quantunitid==62 | quantunitid==3) & quantity>1 & calc_duration_days!=0 ) 

standard_inject$calc_duration_days=as.numeric(standard_inject$calc_duration_days)
standard_inject$quantity=as.numeric(standard_inject$quantity)

standard_inject_seperate <- data.frame()

# Loop through each row of the dataset
for (i in 1:nrow(standard_inject)) {
  print(100*i/nrow(standard_inject))
  # Get the current row data
  row <- standard_inject[i, ]
  
  # Calculate dose interval
  dose_interval <- row$calc_duration_days / (row$quantity )
  
  # Generate sequence of dates where doses will be given
  dose_dates <- sapply(1:row$quantity, function(j) {
    row$issuedate + (j - 1) * dose_interval
  })
  
  # Create new rows with quantity set to 1 for each dose
  new_rows <- row %>%
    # Repeat the row data for the number of doses
    slice(rep(1, row$quantity)) %>%
    # Update the issuedate column for each dose
    mutate(issuedate = dose_dates, quantity = 1, calc_duration_days=1)  # Set quantity to 1 for each dose
  
  
  # Append to the result dataframe
  standard_inject_seperate  <- rbind(standard_inject_seperate , new_rows)
  
}
standard_inject_seperate$issuedate=as.Date(standard_inject_seperate$issuedate)
standard_6=rbind(standard_notinject, standard_inject_seperate )
rm(standard_inject_seperate, standard_notinject, standard_inject)


############################# Remove missing data ##############################


# Assign variable for complete records
standard_6 = standard_6 %>%
  mutate(standard_complete = if_else(!is.na(calc_duration_days) , 1, 0))  
sum(standard_6$standard_complete==1)/nrow(standard_6)*100 

# Limit to prescriptions with complete data
standard_7 = standard_6 %>% filter(standard_complete==1)
nrow(standard_7) # 996,528

# Save processed standard prescriptions to secure drive
write_parquet(standard_7, paste0(parquet_processed, "/filtered_standard.parquet"))




