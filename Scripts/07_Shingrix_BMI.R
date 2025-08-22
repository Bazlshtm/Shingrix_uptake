################################################################################
# Author: Eleanor Barry
# Date: 29/11/2024
# Version: R 4.3.0
# File name: 07_Shingrix_BMI.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## Shingrix_risktimes.rds
## codelist_bmi_aurum.txt
## combined_med
## Numunit.txt
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
# Description of file: This has been altered from Dr Angel Wong's code (https://github.com/ehr-lshtm/cprd-algorithms/blob/main/CPRD_AURUM/pr_getallbmirecords_Aurum.do) 
# Actions:
## Read in libraries
## Finds age at measurement and reduces those who are <18 at measurement
## Reduce to those with suitable convertible units to find height, weight, BMI
## Reduce implausible values
## Reduce conflicting records on the same day
## Find latest measurements 
## Save file as binary overweight/obese (1) or not (0) flag
################################################################################


##### Read in medical codes to find BMI codes ##################################


file_paths <- sprintf(med_combined, 1:25)

combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Load in relevant patients

Shingrix_risktimes=read_rds(paste0(processed_data, "Shingrix_Risks.rds" ))
Shingrix_risktimes = Shingrix_risktimes %>%
  select(patid, yob, dob, min_entry)
obsfile=inner_join(combined_med, Shingrix_risktimes, by="patid" )

###  BMI codes

codelist_bmi_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_bmi_aurum.txt"), 
                                 delim = "\t", escape_double = FALSE, 
                                 col_types = cols(medcodeid = col_character(), 
                                 observations = col_character(), originalreadcode = col_character(), 
                                 snomedctconceptid = col_character(), 
                                 snomedctdescriptionid = col_character(), 
                                 release = col_skip(), emiscodecategoryid = col_skip(), 
                                 weight = col_character(), height = col_character(), 
                                 bmi = col_character()), trim_ws = TRUE)

# Join relavent patients to bmi medcodes
obs_data=inner_join(obsfile, codelist_bmi_aurum, by=c("medcodeid"))

# Look up units of measurements
numunit_data=read_delim(paste0(data_dir_lookup, "2024_03/NumUnit.txt"), 
                        delim = "\t", escape_double = FALSE, 
                        col_types = cols(numunitid = col_character()), 
                        trim_ws = TRUE)
  
# Merge dataset and measurement units by patid
  df_BMI <- obs_data %>%
    left_join(numunit_data, by = "numunitid")
  
# Convert to numeric
  df_BMI <- df_BMI %>%
    mutate(yob = as.numeric(yob),
           value = as.numeric(value))

# Drop rows with missing observation date or missing value, restrict to last 10 years
  df_BMI2 <- df_BMI %>%
    filter(!is.na(obsdate) & !is.na(value) & obsdate>=as.Date("01-01-2011", format="%d-%m-%Y"))
    nrow(df_BMI2)/nrow(df_BMI)
    length(unique(df_BMI2$patid))/length(unique(df_BMI$patid)) 
  
# Create enttype variable to classify measurement type
  df_BMI2 <- df_BMI2 %>%
    mutate(enttype = case_when(
      weight==9 ~ "weight",
      height==9 ~ "height",
      TRUE ~ NA_character_
    ))

# Height
  
  # Make sure units are viable, currently there are a number with units feet, kg, 
  # kg/m2, kgs, mm/Hg, mmHg and O/E -height within 10% average
  # Convert height from cm to meters
  df_BMI3 <- df_BMI2 %>%
    filter(
      (enttype != "height") |  # Keep all rows where enttype is not weight
      (enttype == "height" & Description %in% c("cm", "cms", "m", "metres"))  
    )
  nrow(df_BMI3)/nrow(df_BMI2) 
  
  # Turn cm into m
  df_BMI3 <- df_BMI3 %>%
    mutate(value = ifelse(enttype == "height" & (Description == "cm" | Description == "cms"), value/100, value))

# Weight
  
  # Remove all weights that do not relate to Kg or Stone
  df_BMI4 <- df_BMI3 %>%
    filter(
      (enttype != "weight") |  # Keep all rows where enttype is not weight
        (enttype == "weight" & Description %in% c("/kg(body wt)", "kg", "kg.", "Kgs", "kilograms", "Kilos", "Weight (kgs)", "Weight in Kg", "WEIGHT IN KILOS", "st", "stone", "Stones"))  
    )
  nrow(df_BMI4)/nrow(df_BMI3) 
  
df_BMI4 <- df_BMI4 %>%
    mutate(value = ifelse(enttype == "weight" & (Description == "stone" | Description == "st" | Description == "Stones"), value*6.35, value))
  

  # Drop implausible weights and heights
  df_BMI6 <- df_BMI4 %>%
    filter(!(enttype == "weight" & value < 20),
           !(enttype == "weight" & value >= 400),
           !(enttype == "height" & value < 1.2),
           !(enttype == "height" & value >= 3))
  nrow(df_BMI6)/nrow(df_BMI4) 
  length(unique(df_BMI6$patid))/length(unique(df_BMI4$patid))
  
  # Remove patient df_BMI if it is later than min entry date (would need to be change for time-varying covariates)
  df_BMI6  <- df_BMI6  %>%
    filter(min_entry >= obsdate)  
  
##################### Find patient's whose BMI is latest #######################
  # Step 1: Organise data from long to wide
  df_height <- df_BMI6 %>%
    filter(enttype == "height") %>%
    select(patid, value, obsdate) %>%
    rename(height_value = value) %>%
    group_by(patid, obsdate) %>%
    filter(diff(range(height_value)) <= 0.5) # If differences between same day records >0.5m, reject
  
  df_height <- df_height %>% distinct()
    
  df_weight <- df_BMI6 %>%
    filter(enttype == "weight") %>%
    select(patid, value, obsdate) %>%
    rename(weight_value = value) %>%
    group_by(patid, obsdate) %>%
    filter(diff(range(weight_value)) <= 5) # If differences between same day records >5kg, reject
  
  df_weight <- df_weight %>% distinct()

  # Step 2: Keep only the latest observation date for each patient
  df_height_latest <- df_height %>%
    group_by(patid) %>%
    filter(obsdate == max(obsdate)) %>%
    summarise(height = mean(height_value, na.rm = TRUE), height_obsdate = max(obsdate), .groups = "drop")
  
  df_weight_latest <- df_weight %>%
    group_by(patid) %>%
    filter(obsdate == max(obsdate)) %>%
    summarise(weight = mean(weight_value, na.rm = TRUE), weight_obsdate = max(obsdate), .groups = "drop")
  

  # Step 3: Join df_BMI frames by patientID 
  df_long <- df_height_latest %>%
    merge(df_weight_latest, by = c("patid")) %>%
    mutate(BMI =  weight / (height ^ 2))  # Use bmi if bmi_date is more recent
  
  df_long$BMI=ifelse(df_long$BMI>=25, 1,0)
  
  write_rds(df_long[c("patid", "BMI")], paste0(processed_data, "BMI.rds" ))
  

 
