################################################################################
# Author: Eleanor Barry
# Date: 29/11/2024
# Version: R 4.3.0
# File name: 08_Shingrix_Ethnicity.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID-filtered.rds
## filtered_med.parquet
## combined_med
## codelist_ethnicity.txt
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
# Description of file: From Dr Helen Strongman's codelist to Ethnicity coding. Discrepancies have been classed as "Mixed"
# Actions:
## Read in libraries
## Read in ethnicity
## Find overall risk times 
################################################################################

################################  Load in medcodes #############################

# Generate paths to 25 combined_med parquet files
file_paths <- sprintf(med_combined, 1:25)

# Read and combine the parquet files into one dataframe
combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Load the ethnicity codelist (Aurum format)
codelist_ethnicity_aurum <- read_delim(
  paste0(codelists, "Medcodes/codelist_ethnicity_aurum.txt"), 
  delim = "\t", 
  escape_double = FALSE, 
  col_types = cols(
    medcodeid = col_character(), 
    snomedctconceptid = col_character(), 
    snomedctdescriptionid = col_character()
  ), 
  trim_ws = TRUE
)

################################ Ethnicity5 analysis ################################

# Load the Shingrix risk observation times
Shingrix_risktimes <- read_rds(paste0(processed_data, "Shingrix_Risks.rds")) %>%
  select(patid)

# Join ethnicity codes to medical codes, then join to risk times
obsfile <- combined_med %>%
  inner_join(codelist_ethnicity_aurum, by = "medcodeid") %>%
  inner_join(Shingrix_risktimes, by = "patid")

# Tabulate Ethnicity5 and Ethnicity16 counts for checking
janitor::tabyl(obsfile$eth5)
janitor::tabyl(obsfile$eth16)

# Remove duplicates and records with "Not stated" ethnicity
obsfile <- obsfile %>%
  arrange(patid) %>%
  mutate(sysyear = year(obsdate)) %>%
  group_by(patid, sysyear, medcodeid) %>%
  mutate(duplicates = n()) %>%
  group_by(patid, sysyear) %>%
  mutate(
    first_entry = if_else(duplicates > 1, min(obsdate), NA_Date_),
    keep_row = (duplicates == 1) | (obsdate == first_entry)
  ) %>%
  ungroup() %>%
  filter(keep_row, eth5 != "5. Not Stated")

# Count the number of distinct Ethnicity5 codes per patient
df_count <- obsfile %>%
  group_by(patid) %>%
  summarise(eth5_count = n_distinct(eth5))

table(df_count$eth5_count)

# Collapse multiple ethnicities to "Other/Mixed" where needed
Ethnicity <- obsfile %>%
  filter(eth5 != "5. Not Stated") %>%
  group_by(patid) %>%
  select(patid, Ethnicity = eth5) %>%
  mutate(
    Ethnicity = case_when(
      Ethnicity == "3. Other" ~ "3. Other/Mixed",
      Ethnicity == "4. Mixed" ~ "3. Other/Mixed",
      TRUE ~ Ethnicity
    )
  ) %>%
  mutate(
    Ethnicity = if_else(
      n_distinct(Ethnicity) > 1,   # multiple ethnicities recorded
      "3. Other/Mixed",            # collapse to Other/Mixed
      first(Ethnicity)             # otherwise use the first
    )
  ) %>%
  distinct() %>%
  ungroup()

# Save Ethnicity5-coded data
write_rds(Ethnicity, paste0(processed_data, "Ethnicity.rds"))

################################ Ethnicity16 analysis ################################

# Reload the Shingrix risk times to ensure clean environment
Shingrix_risktimes <- read_rds(paste0(processed_data, "Shingrix_Risks.rds")) %>%
  select(patid)

# Join ethnicity codes to medical codes, then join to risk times
obsfile <- combined_med %>%
  inner_join(codelist_ethnicity_aurum, by = "medcodeid") %>%
  inner_join(Shingrix_risktimes, by = "patid")

# Tabulate Ethnicity5 and Ethnicity16 codes for checking
janitor::tabyl(obsfile$eth5)
janitor::tabyl(obsfile$eth16)

# Remove duplicates and records with "Not stated" ethnicity
obsfile <- obsfile %>%
  arrange(patid) %>%
  mutate(sysyear = year(obsdate)) %>%
  group_by(patid, sysyear, medcodeid) %>%
  mutate(duplicates = n()) %>%
  group_by(patid, sysyear) %>%
  mutate(
    first_entry = if_else(duplicates > 1, min(obsdate), NA_Date_),
    keep_row = (duplicates == 1) | (obsdate == first_entry)
  ) %>%
  ungroup() %>%
  filter(keep_row, eth16 != "17. Not Stated")

# Count distinct Ethnicity16 codes per patient
df_count <- obsfile %>%
  group_by(patid) %>%
  summarise(eth16_count = n_distinct(eth16))

table(df_count$eth16_count)

# Collapse Ethnicity16 codes into consistent groupings:
#   - resolve conflicts if multiple ethnicities recorded
#   - prioritise by category, otherwise collapse to Other/Mixed
Ethnicity16 <- obsfile %>%
  group_by(patid) %>%
  mutate(
    ethnicity = case_when(
      n_distinct(eth16) == 1 ~ first(eth16),
      
      any(eth16 %in% c("1. British", "2. Irish", "3. Other White")) &
        !any(eth16 %in% c("4. White and Black Caribbean", "5. White and Black African",
                          "6. White and Asian", "7. Other Mixed", "8. Indian", 
                          "9. Pakistani", "10. Bangladeshi", "11. Other Asian", 
                          "12. Caribbean", "13. African", "14. Other Black", 
                          "15. Chinese", "16. Other ethnic group")) ~ "3. Other White",
      
      any(eth16 %in% c("8. Indian", "9. Pakistani", "10. Bangladeshi", "11. Other Asian")) &
        !any(eth16 %in% c("1. British", "2. Irish", "3. Other White", "4. White and Black Caribbean",
                          "5. White and Black African", "6. White and Asian", "7. Other Mixed", 
                          "12. Caribbean", "13. African", "14. Other Black", "15. Chinese", 
                          "16. Other ethnic group")) ~ "11. Other Asian",
      
      any(eth16 %in% c("12. Caribbean", "13. African", "14. Other Black")) &
        !any(eth16 %in% c("1. British", "2. Irish", "3. Other White", "4. White and Black Caribbean",
                          "5. White and Black African", "6. White and Asian", "7. Other Mixed", 
                          "8. Indian", "9. Pakistani", "10. Bangladeshi", "11. Other Asian", 
                          "15. Chinese", "16. Other ethnic group")) ~ "14. Other Black",
      
      TRUE ~ "7. Other Mixed"
    )
  ) %>%
  select(patid, ethnicity) %>%
  distinct() %>%
  ungroup()

# Save Ethnicity16-coded data
write_rds(Ethnicity16, paste0(processed_data, "Ethnicity16.rds")

                                       
