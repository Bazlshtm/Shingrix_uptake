################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 01/09/2024
# Version: R 4.3.0
# File name: 03_Shingrix_filtered_med.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID.rds
## med_combined 
## codelist_immuno
# R scripts needed: None
## 00_Set_up_directory.R
# Data sets created: 
## filtered_med.parquet
## Description of file: As extract file to simplify the medcodes needed for 
## medical study
## Create possible study start and ends 
## Restrict based on those who are not in the study for at least a day
## Save file
################################################################################



############################ Load in med data ##################################

# Generate a vector of file paths for x between 1 and 25
file_paths <- sprintf(med_combined, 1:25)

# Read and combine all the parquet files
combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))
nrow(combined_med)/1e6 # 90.8 million records

# Load in codelist #
codelist_immuno <- read_delim(paste0(codelists, "Medcodes/codelist_immunosuppresion_medical_aurum.txt"), 
                                                delim = "\t", escape_double = FALSE, 
                                                col_types = cols(medcodeid = col_character(), 
                                                snomedctconceptid = col_character(), 
                                                snomedctdescriptionid = col_character(), 
                                                emiscodecategoryid = col_character()), 
                                                trim_ws = TRUE)

# Load in patient list
PATID_filtered <- read_rds(paste0(processed_data, "PATID_filtered.rds")) 

# Subset medical data to eligible codes and patients
combined_med_immuno = combined_med %>%
  filter(medcodeid %in% codelist_immuno$medcodeid) %>%
  filter(patid %in% PATID_filtered$patid )
nrow(combined_med_immuno) #1,811,429 relevant records retained 

# Subset to eligible entries in valid date range
combined_med_immuno = combined_med_immuno %>%
  filter(obsdate>=as.Date("1942-01-01") | enterdate>=as.Date("1942-01-01"))
nrow(combined_med_immuno) #1,811,093 relevant records retained 

# Remove invalid obsdate to allow later filling in by enterdate
combined_med_immuno <- combined_med_immuno %>%
  mutate(obsdate = if_else(obsdate <= as.Date("1942-01-01"), NA_Date_, obsdate))

# Remove unused data
rm(combined_med)

# Merge with codelist metadata
combined_med_immuno=left_join(combined_med_immuno, 
                              codelist_immuno, 
                              by="medcodeid")

# Quantify entries with missing obsdate
sum(is.na(combined_med_immuno$obsdate)) #1,496
sum(is.na(combined_med_immuno$obsdate))/nrow(combined_med_immuno)*100 #0.08%

# Check association between enterdate and obsdate where observed
combined_med_immuno$enterdate_diff = as.numeric(combined_med_immuno$enterdate-combined_med_immuno$obsdate)
sum(combined_med_immuno$enterdate_diff %in% -60:60, na.rm=T)/nrow(combined_med_immuno) # >87% of enter dates within 60d of obs date where both observed

# Plot consistency across categories
ggplot(combined_med_immuno, aes(enterdate_diff)) + 
  geom_density() + xlim(-365,365) +
  facet_grid(immuno_cat~.)

## Handling of missing obsdate
# Overall, >99.9% of records have a recorded obsdate
# In the remaining <0.1% records, we estimate that 87% can be assigned an obsdate to within 60 days based on differences between enterdate and obsdate where both observed
# Consistency between enterdate and obsdate apparent across subgroups
# On this basis, we used enterdate where obsdate unavailable
combined_med_immuno$obsdate[is.na(combined_med_immuno$obsdate)]=combined_med_immuno$enterdate[is.na(combined_med_immuno$obsdate)]

# Write temporary file of eligible records

write_parquet(combined_med_immuno, raw_processed_med)
# combined_med_immuno <- read_parquet(raw_processed_med) # Read in filtered parquet to run script from here



