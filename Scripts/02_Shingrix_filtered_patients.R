################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 01/09/2024
# Version: R 4.3.0
# File name: 02_Shingrix_filtered_patients.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID.dta
# R scripts needed: 
## 00_Shingrix_set_up_directory.R
# Data sets created: 
## PATID_filtered.rds
# Description of file: 
## Takes CPRD's patient file and works out study start and end for each patient 
## valid for our study
# Actions performed: 
## Load in CPRD Aurum data
## Create cohort study start and end dates
## Restrict based on those who are not in the study for at least a day
################################################################################



############################# Load in Patient ID  ##############################

# Load in patient data from CPRD acceptable patients who are born between 
# 1942 and 1953
PATID <- read_dta(paste0(codelists, "PATID.dta"))

# Specify dates
PATID$regstartdate=as.Date(PATID$regstartdate, "%d/%m/%Y")
PATID$regenddate=as.Date(PATID$regenddate, "%d/%m/%Y")
PATID$lcd=as.Date(PATID$lcd, "%d/%m/%Y")
PATID$cprd_ddate=as.Date(PATID$cprd_ddate, "%d/%m/%Y")
PATID$study_startdate = as.Date("01/09/2021", "%d/%m/%Y")
PATID$study_enddate = as.Date("31/08/2023", "%d/%m/%Y")

# Set birth date to 01/07 of year of birth
PATID$dob = as.Date(paste0("01/07/",PATID$yob), "%d/%m/%Y")
PATID$dob_plus70 = PATID$dob %m+% years(70)
PATID$dob_plus80 = PATID$dob %m+% years(80)

# Set minimum cohort entry (prior to immunosuppression assignment) based on 
# study start, age, start of registration 
PATID$min_cohort_entry = pmax(PATID$study_startdate, 
                              PATID$dob_plus70, 
                              PATID$regstartdate %m+% years(1), na.rm=TRUE)

# Set maximum cohort exit (prior to immunosuppression assignment) based on on 
# study end, age, death, last collection, end of registration
PATID$max_cohort_exit = pmin(PATID$study_enddate, 
                             PATID$dob_plus80, 
                             PATID$cprd_ddate, 
                             PATID$lcd, 
                             PATID$regenddate, na.rm=TRUE)

############################# Initial filtering  ###############################

# Count number of potentially eligible individuals
nrow(PATID) # 3,669,708

# Check all born between 1942 and 1953 as per extract
all(PATID$yob>=1942 & PATID$yob<=1953) # TRUE

# Filter to individuals with alive, registered and at collecting practices after 
# study_startdate
PATID_filtered = PATID %>%
  filter(min_cohort_entry <= max_cohort_exit) %>%
  filter(region!=12) %>%
  filter(gender==1 | gender==2)
nrow(PATID_filtered) # 1,519,720


# Write temporary file of eligible patients
write_rds(PATID_filtered, paste0(processed_data, "PATID_filtered.rds" ))
# filtered PATIDs to run script from here

# Remove unused data
rm(PATID)

##End file## 


