################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 29/10/2024
# Version: R 4.3.0
# File name: 00_Set_up_directory.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: None
# R scripts needed: None
# Data sets created: None 
# Description of file: 
## Read in libraries
## Set up directories for inputs and outputs
################################################################################

# Libraries
library(readr)
library(haven)
library(dplyr)
library(lubridate)
library(slider)
library(stringr)
library(tidyr)
library(purrr)
library(arrow)
library(ggplot2)
library(kableExtra)
library(data.table)
library(janitor)
library(forestplot)
library(forester)
library(multiwayvcov)
library(lmtest)

# File paths
project <- # File location of the Shingrix project
codelists <- # File location of the Shingrix codelists
scripts <- # File location of the Shingrix R scripts
processed_data <- # File location of the processed data
med_combined <-  # Extracted CPRD medical files with relevant files for covariates, made in script 01
drug_combined <- # Extracted CPRD product files with relevant files for covariates, made in script 01
raw_processed_med <- # Extracted from CPRD based on immunosuppressed definitions 
data_dir_lookup <- # File location for CPRD data lookups
imd_linkage_file <- # File location for IMD linkage 

# Codelists
codelist_path_steroids <- # File location of the Shingrix steroids codelists
codelist_path_standard <- # File location of the Shingrix standard drugs codelists
codelist_path_targeted <- # File location of the Shingrix targeted drugs codelists
codelist_path_medical_immuno <-  # File location of the Shingrix medical immunosuppression codelists

# Lookup
product_lookup_path <- # Directions to 202403_EMISProductDictionary.txt
common_dosages_sep2024_path <- # Directions to 2024_09/common_dosages.txt
numunitid <- # Directions to 2024_09/numunitid.txt
imd_linkage <- # Directions to IMD linkage

# Extract files link
med_zip_file <- # Direction to CPRD medical files
med_extract_text_start <- # Direction to CPRD medical files
med_extract_text_mid <- # Direction to CPRD medical files
raw_data_meds <- # Direction to CPRD medical files
raw_data_drugs <- # Direction to CPRD product files
med_combined_file <- # Direction to CPRD medical files
drug_extract_text_start <- # Direction to CPRD product files
drug_extract_text_mid <- # Direction to CPRD product files
drug_combined_file <- # Direction to CPRD product files
imd_linked <- # Directions to CPRD IMD linkage
 
