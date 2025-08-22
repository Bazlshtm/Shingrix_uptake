################################################################################
# Author: Eleanor Barry
# Date: 06/11/2024
# Version: R 4.3.0
# File name: 12_Shingrix_Zoster_Uptake.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## Shingrix_Cov
## Combined_med
## Combined_prod
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
## Full_data
# Description of file: Finds the first and second vacc date for each person. Removes people if they have vaccination prior to entry
# Actions:
## Finds relevant zoster vacc med and prod codes for each patient
## Removed patients who had vaccine before 01-09-2021 and if date not given
## Finds first date of vaccine, and type of vaccine up to a week after that, assigns categories based on the week of data
## Remove patients whose vaccine was before min entry
##
################################################################################


############################## Load in data  ###############################

# Generate file paths for combined_med data
file_paths <- sprintf(med_combined, 1:25)

# Read in and combine all medical observation files
combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Generate file paths for combined_drug data
file_paths <- sprintf(drug_combined, 1:20)

# Read in and combine all drug prescription files
combined_drug <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Check total number of drug records
nrow(combined_drug) # 12,504,786 records

# Load in main cohort file
Shingrix_df = read_rds(paste0(processed_data, "Shingrix_Cov.rds"))

# Load medical zoster vaccine codelist
codelist_zoster_vacc_medical_aurum <- read_delim(
  "J:/EHR-Working/Baz/Shingrix/Codelists/Medcodes/codelist_zoster_vacc_medical_aurum.txt", 
  delim = "\t", escape_double = FALSE,
  col_types = cols(
    medcodeid = col_character(), 
    Zostavax = col_character(), 
    Shingrix = col_character(), 
    Neutral = col_character()
  ), trim_ws = TRUE
)

# Load product zoster vaccine codelist
codelist_zoster_vacc_product_aurum <- read_delim(
  "J:/EHR-Working/Baz/Shingrix/Codelists/Prodcodes/codelist_zoster_vacc_product_aurum.txt", 
  delim = "\t", escape_double = FALSE,
  col_types = cols(
    prodcodeid = col_character(), 
    dmdid = col_character(), 
    Zostavax = col_character(), 
    Shingrix = col_character(), 
    Neutral = col_character()
  ), trim_ws = TRUE
)

# Extract and join medical vaccine records for cohort
zoster_med <- combined_med %>%
  filter(medcodeid %in% codelist_zoster_vacc_medical_aurum$medcodeid & patid %in% Shingrix_df$patid) %>%
  select(patid, medcodeid, obsdate) %>%
  left_join(Shingrix_df %>% select(patid), by = "patid") %>%
  filter(obsdate >= as.Date("1942-01-01")) %>%
  rename(vacc_date = obsdate)

# Attach vaccine type labels
zoster_med = left_join(zoster_med, codelist_zoster_vacc_medical_aurum, by = "medcodeid")

# Extract and join product vaccine records for cohort
zoster_prod <- combined_drug %>%
  filter(prodcodeid %in% codelist_zoster_vacc_product_aurum$prodcodeid & patid %in% Shingrix_df$patid) %>%
  select(patid, prodcodeid, issuedate) %>%
  left_join(Shingrix_df %>% select(patid), by = "patid") %>%
  filter(issuedate >= as.Date("1942-01-01")) %>%
  rename(vacc_date = issuedate)

# Attach product descriptions
zoster_prod = left_join(zoster_prod, codelist_zoster_vacc_product_aurum, by = "prodcodeid")    
zoster_prod = rename(zoster_prod, term = termfromemis)

# Combine medical and product vaccine records
zoster_vacc = bind_rows(
  zoster_med[c("patid", "Shingrix", "Neutral", "Zostavax", "term", "vacc_date", "medcodeid")],
  zoster_prod[c("patid", "Shingrix", "Neutral", "Zostavax", "term", "vacc_date")]
)

################## Remove any patient who has an existing vaccine ##########

# Join vaccine data with cohort and remove those vaccinated before entry
Shingrix_df_2 = left_join(Shingrix_df, zoster_vacc, by = "patid")
Shingrix_df_2 <- Shingrix_df_2 %>%
  group_by(patid) %>%
  mutate(before_entry = any(vacc_date < min_entry, na.rm = TRUE)) %>%
  filter(!before_entry | is.na(vacc_date)) %>%
  select(-before_entry)

# Number of patients remaining in the cohort
length(unique(Shingrix_df_2$patid)) # 99,647

# Store original, cleaned cohort
Shingrix_df_original = Shingrix_df[Shingrix_df$patid %in% Shingrix_df_2$patid,]

##################### Find first vaccine ###################################

# Determine first vaccine type within 7-day window
Shingrix_df_3 <- Shingrix_df_2 %>%
  group_by(patid) %>%
  mutate(
    first_vacc_date = min(vacc_date),
    diff = as.numeric(difftime(vacc_date, first_vacc_date, units = "days")),
    within_seven_days = diff >= 0 & diff <= 7
  ) %>%
  filter(within_seven_days) %>%
  filter(term != "Administration of second dose of Varicella-zoster vaccine for shingles" &
           term != "Administration of second dose of vaccine product containing only Human alphaherpesvirus 3 antigen for shingles") %>%
  mutate(
    first_vacc_type = case_when(
      sum(Shingrix == 1, na.rm=TRUE) > 0 & sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) == 0 & sum(Zostavax == 1, na.rm=TRUE) == 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) > 0 ~ "Shingrix",
      sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Zostavax",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(patid, first_vacc_date, first_vacc_type) %>%
  distinct()

# Add first vaccine info to cohort
Shingrix_df_4 = left_join(Shingrix_df_original, Shingrix_df_3, by = "patid")

# Remove vaccine records after follow-up end
Shingrix_df_4 <- Shingrix_df_4 %>%
  mutate(
    condition = first_vacc_date > max_end,
    first_vacc_type = if_else(condition, NA_character_, first_vacc_type),
    first_vacc_date = if_else(condition, NA_Date_, first_vacc_date)
  ) %>%
  select(-condition)

# Join mid-stage vaccine data
Shingrix_df_midway = left_join(Shingrix_df_2, Shingrix_df_4[c("patid", "first_vacc_type", "first_vacc_date")], by = "patid")

##################### Find second vaccine ###################################

# Identify valid second doses (7 weeks to 13 months)
Shingrix_df_5 <- Shingrix_df_midway %>%
  filter(
    vacc_date >= first_vacc_date + weeks(7) &
      vacc_date <= first_vacc_date + months(13)
  ) %>%
  group_by(patid) %>%
  mutate(
    second_vacc_date = min(vacc_date),
    diff = as.numeric(difftime(vacc_date, second_vacc_date, units = "days")),
    within_seven_days = diff >= 0 & diff <= 7
  ) %>%
  filter(within_seven_days) %>%
  filter(term != "Administration of first dose of Varicella-zoster vaccine for shingles" &
           term != "Administration of first dose of vaccine product containing only Human alphaherpesvirus 3 antigen for shingles") %>%
  mutate(
    second_vacc_type = case_when(
      sum(Shingrix == 1, na.rm=TRUE) > 0 & sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) == 0 & sum(Zostavax == 1, na.rm=TRUE) == 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) > 0 ~ "Shingrix",
      sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Zostavax",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(patid, second_vacc_date, second_vacc_type) %>%
  distinct()

################## Count other vaccine doses (exploratory) ##################

Shingrix_df_other <- Shingrix_df_midway %>%
  filter(
    vacc_date >= first_vacc_date + weeks(1) &
      vacc_date <= first_vacc_date + months(24)
  ) %>%
  group_by(patid) %>%
  mutate(
    other_vacc_date = min(vacc_date),
    diff = as.numeric(difftime(vacc_date, other_vacc_date, units = "days")),
    within_seven_days = diff >= 0 & diff <= 7
  ) %>%
  filter(within_seven_days) %>%
  filter(term != "Administration of first dose of Varicella-zoster vaccine for shingles" &
           term != "Administration of first dose of vaccine product containing only Human alphaherpesvirus 3 antigen for shingles") %>%
  mutate(
    other_vacc_type = case_when(
      sum(Shingrix == 1, na.rm=TRUE) > 0 & sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) == 0 & sum(Zostavax == 1, na.rm=TRUE) == 0 ~ "Neutral",
      sum(Shingrix == 1, na.rm=TRUE) > 0 ~ "Shingrix",
      sum(Zostavax == 1, na.rm=TRUE) > 0 ~ "Zostavax",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  select(patid, other_vacc_date) %>%
  distinct()

######################## Merge final dataset ################################

Shingrix_df_6 = left_join(Shingrix_df_4, Shingrix_df_5, by = "patid")

# Optional: Remove second dose after max_end (currently commented out)
# Shingrix_df_7 <- Shingrix_df_6 %>%
#   mutate(
#     condition = second_vacc_date > max_end,
#     second_vacc_type = if_else(condition, NA_character_, second_vacc_type),
#     second_vacc_date = if_else(condition, NA_Date_, second_vacc_date)
#   ) %>%
#   select(-condition)

# Save final output
write_rds(Shingrix_df_6, paste0(processed_data, "Shingrix_zoster.rds"))

    
