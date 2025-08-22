################################################################################
# Author: Eleanor Barry
# Date: 06/11/2024
# Version: R 4.3.0
# File name: 12_Shingrix_flu.R
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
# Description of file: Reconstructed from Dr Anne Suffel's code
# Suffel, A.M., Walker, J.L., Campbell, C., Carreira, H., Warren‐Gash, C. and McDonald, H.I., 2024. 
# A New Validated Approach for Identifying Childhood Immunizations in Electronic Health Records in the United Kingdom. 
# Pharmacoepidemiology and drug safety, 33(8), p.e5848.
# Imports and combines medical and drug data  
# Identifies flu vaccination records using codelists  
# Classifies patients (vaccinated, delayed, declined, unvaccinated, conflict)  
# Extracts first flu vaccination dates for 2021/22 and 2022/23 seasons  
# Compares flu and Shingrix vaccination dates (same day / within 1–2 weeks)  
# Saves updated cohort dataset for further analyses
################################################################################


############################## Load in data  ###############################

# Create file paths for medical data parquet files (assumes med_combined is a format string)
file_paths <- sprintf(med_combined, 1:25)

# Read and combine all medical parquet files into one dataframe
combined_med <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Create file paths for drug data parquet files
file_paths <- sprintf(drug_combined, 1:20)

# Read and combine all drug parquet files into one dataframe
combined_drug <- file_paths %>%
  map_dfr(~ read_parquet(.x))

# Check number of records in drug data
nrow(combined_drug) # 12,504,786 records

# Load flu vaccine codelists: medcodes
codelist_flu_medical_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_flu_medical_aurum.txt"), 
                                         delim = "\t", escape_double = FALSE, 
                                         col_types = cols(medcodeid = col_character()), 
                                         trim_ws = TRUE)

# Remove duplicates and codes marked for removal
codelist_flu_medical_aurum <- unique(codelist_flu_medical_aurum)
codelist_flu_medical_aurum <- codelist_flu_medical_aurum[is.na(codelist_flu_medical_aurum$remove),]

# Load flu vaccine codelists: prodcodes
codelist_flu_product_aurum <- read_delim(paste0(codelists, "Prodcodes/codelist_flu_product_aurum.txt"), 
                                         col_types = cols(prodcodeid = col_character()), 
                                         delim = " ")

# Remove duplicates and incomplete records
codelist_flu_product_aurum <- unique(codelist_flu_product_aurum)
codelist_flu_product_aurum <- codelist_flu_product_aurum[complete.cases(codelist_flu_product_aurum),]

# Load main cohort data
Shingrix_flu = read_rds(paste0(processed_data, "Shingrix_zoster.rds"))

# Extract flu-related medical records for patients in the cohort from August 2021 onward
Flu_med <- combined_med %>%
  filter(medcodeid %in% codelist_flu_medical_aurum$medcodeid & patid %in% Shingrix_flu$patid) %>%
  left_join(Shingrix_flu %>% select(patid), by = "patid") %>%
  filter(obsdate >= as.Date("2021-08-01")) %>%
  rename(flu_date = obsdate) %>%
  left_join(codelist_flu_medical_aurum, by = "medcodeid")

# Extract flu-related drug records for patients in the cohort from August 2021 onward
Flu_prod <- combined_drug %>%
  filter(prodcodeid %in% codelist_flu_product_aurum$prodcodeid & patid %in% Shingrix_flu$patid) %>%
  select(patid, prodcodeid, issuedate) %>%
  left_join(Shingrix_flu %>% select(patid), by = "patid") %>%
  filter(issuedate >= as.Date("2021-08-01")) %>%
  rename(flu_date = issuedate)

# Add an indicator for product-based records
Flu_prod$product = 1

# Combine medical and drug-based flu records into one dataset
Flu = bind_rows(Flu_med, Flu_prod)

# Replace NA values in classification flags with 0
Flu <- Flu %>%
  mutate(across(c("given", "given_lag", "neutral", "adv_reac", "declined", "remove", "product"), ~replace(., is.na(.), 0)))

# Convert to data.table and aggregate codes per patient and flu date
setDT(Flu)
tmp_filtered <- Flu[, by = c("patid", "flu_date"), lapply(.SD, sum),
                    .SDcols = c("given", "given_lag", "neutral", "declined", "product")]
tmp_filtered <- tmp_filtered[order(patid)]

# Assign flu vaccination status based on coding combinations
tmp_filtered[declined == 0 & product == 0 & given > 0 & neutral == 0, status := "vaccinated"]
tmp_filtered[declined == 0 & product > 0 & given == 0 & neutral == 0, status := "vaccinated"]
tmp_filtered[declined == 0 & product == 0 & given == 0 & neutral > 0, status := "vaccinated"]
tmp_filtered[neutral > 0 & given > 0, status := "vaccinated"]
tmp_filtered[neutral > 0 & product > 0, status := "vaccinated"]
tmp_filtered[given > 0 & product > 0, status := "vaccinated"]

# Flag declined cases
tmp_filtered[neutral > 0 & declined > 0, status := "declined"]
tmp_filtered[declined > 0 & product == 0 & given == 0 & neutral == 0, status := "declined"]

# Flag conflicting records (e.g., both given and declined codes)
tmp_filtered[declined > 0 & given > 0, status := "conflict"]
tmp_filtered[declined > 0 & product > 0, status := "conflict"]

# Flag delayed flu vaccinations
tmp_filtered[declined == 0 & product == 0 & given == 0 & given_lag > 0, status := "delayed"]

# Flag patients with no evidence of flu vaccination
tmp_filtered[declined == 0 & product == 0 & given == 0 & given_lag == 0, status := "unvaccinated"]

# Check if any status was not assigned
if (any(is.na(tmp_filtered$status))) {
  n <- nrow(tmp_filtered[is.na(status)])
  print(paste0("Error - ", n, " unassigned vaccination statuses"))
}

# Keep only vaccinated or delayed statuses
tmp_filtered <- tmp_filtered[status == "vaccinated" | status == "delayed"]

# First flu season: Sep 2021 to Mar 2022
tmp_filtered_21 = subset(tmp_filtered, tmp_filtered$flu_date >= "2021-09-01" & tmp_filtered$flu_date <= "2022-03-31")
tmp_filtered_21[status == "vaccinated", dose_records := seq_len(.N), by = "patid"]
tmp_filtered_21[status == "delayed", dose_records_del := seq_len(.N), by = "patid"]
tmp_filtered_21 <- tmp_filtered_21[status == "delayed" | dose_records == 1]

# Select earliest record for vaccinated or delayed status
tmp1 <- tmp_filtered_21[status == "vaccinated", list(patid, flu_date, status)]
tmp2 <- tmp_filtered_21[status == "delayed", list(patid, flu_date, status)]
both <- merge(tmp1, tmp2, by = "patid", all = TRUE)

# Determine earliest flu date
both[!is.na(flu_date.x) & is.na(flu_date.y), first_flu_date := flu_date.x]
both[is.na(flu_date.x) & !is.na(flu_date.y), first_flu_date := flu_date.y]
both[!is.na(flu_date.x) & !is.na(flu_date.y), first_flu_date := min(flu_date.x, flu_date.y)]

# Handle duplicates (only keep earliest per patient)
duplicates <- both[which(duplicated(both$patid))]
dupes <- both[patid %chin% duplicates$patid]
dupes <- dupes[order(first_flu_date)]
dupes[, num := seq_len(.N), by = "patid"]
dupes[, num := rowid(patid)]
dupes <- dupes[num == 1]
no_dupes <- both[!patid %chin% duplicates$patid]
no_dupes <- no_dupes[, list(patid, first_flu_date)]
dupes <- data.table(dupes)[, list(patid, first_flu_date)]
all_flu_vac <- rbind(no_dupes, dupes)

# Merge with main cohort
Shingrix_flu = left_join(Shingrix_flu, all_flu_vac, by = "patid")

# Second flu season: Sep 2022 to Mar 2023
tmp_filtered_22 = subset(tmp_filtered, tmp_filtered$flu_date >= "2022-09-01" & tmp_filtered$flu_date <= "2023-03-31")
tmp_filtered_22[status == "vaccinated", dose_records := seq_len(.N), by = "patid"]
tmp_filtered_22[status == "delayed", dose_records_del := seq_len(.N), by = "patid"]
tmp_filtered_22 <- tmp_filtered_22[status == "delayed" | dose_records == 1]

# Select earliest record for second season
tmp1 <- tmp_filtered_22[status == "vaccinated", list(patid, flu_date, status)]
tmp2 <- tmp_filtered_22[status == "delayed", list(patid, flu_date, status)]
both <- merge(tmp1, tmp2, by = "patid", all = TRUE)

# Determine earliest flu date
both[!is.na(flu_date.x) & is.na(flu_date.y), second_flu_date := flu_date.x]
both[is.na(flu_date.x) & !is.na(flu_date.y), second_flu_date := flu_date.y]
both[!is.na(flu_date.x) & !is.na(flu_date.y), second_flu_date := min(flu_date.x, flu_date.y)]

# Handle duplicates for second season
duplicates <- both[which(duplicated(both$patid))]
dupes <- both[patid %chin% duplicates$patid]
dupes <- dupes[order(second_flu_date)]
dupes[, num := seq_len(.N), by = "patid"]
dupes[, num := rowid(patid)]
dupes <- dupes[num == 1]
no_dupes <- both[!patid %chin% duplicates$patid]
no_dupes <- no_dupes[, list(patid, second_flu_date)]
dupes <- data.table(dupes)[, list(patid, second_flu_date)]
all_flu_vac <- rbind(no_dupes, dupes)

# Merge second season dates into main cohort
Shingrix_flu = left_join(Shingrix_flu, all_flu_vac, by = "patid")

# Calculate time difference (in days) between flu and Shingrix doses
Shingrix_flu$first_flu_vacc_diff = as.numeric(pmin(abs(Shingrix_flu$first_flu_date - Shingrix_flu$first_vacc_date), abs(Shingrix_flu$first_flu_date - Shingrix_flu$second_vacc_date)))
Shingrix_flu$second_flu_vacc_diff = as.numeric(pmin(abs(Shingrix_flu$second_flu_date - Shingrix_flu$first_vacc_date), abs(Shingrix_flu$second_flu_date - Shingrix_flu$second_vacc_date)))

# Initialise columns
Shingrix_flu$first_flu_vacc_type = NA
Shingrix_flu$second_flu_vacc_type = NA

# Classify first flu vaccination timing
for (i in 1:nrow(Shingrix_flu)) {
  if (!is.na(Shingrix_flu$first_flu_vacc_diff[i])) {
    if (Shingrix_flu$first_flu_vacc_diff[i] == 0) {
      Shingrix_flu$first_flu_vacc_type[i] <- "Same day"
    } else if (Shingrix_flu$first_flu_vacc_diff[i] >= 1 && Shingrix_flu$first_flu_vacc_diff[i] <= 6) {
      Shingrix_flu$first_flu_vacc_type[i] <- "Within a week"
    } else if (Shingrix_flu$first_flu_vacc_diff[i] >= 7 && Shingrix_flu$first_flu_vacc_diff[i] <= 14) {
      Shingrix_flu$first_flu_vacc_type[i] <- "1-2 weeks"
    } else {
      Shingrix_flu$first_flu_vacc_type[i] <- NA  
    }
  }
}

# Classify second flu vaccination timing
for (i in 1:nrow(Shingrix_flu)) {
  if (!is.na(Shingrix_flu$second_flu_vacc_diff[i])) {
    if (Shingrix_flu$second_flu_vacc_diff[i] == 0) {
      Shingrix_flu$second_flu_vacc_type[i] <- "Same day"
    } else if (Shingrix_flu$second_flu_vacc_diff[i] >= 1 && Shingrix_flu$second_flu_vacc_diff[i] <= 6) {
      Shingrix_flu$second_flu_vacc_type[i] <- "Within a week"
    } else if (Shingrix_flu$second_flu_vacc_diff[i] >= 7 && Shingrix_flu$second_flu_vacc_diff[i] <= 14) {
      Shingrix_flu$second_flu_vacc_type[i] <- "1-2 weeks"
    } else {
      Shingrix_flu$second_flu_vacc_type[i] <- NA  
    }
  }
}

# Save the updated dataset
write_rds(Shingrix_flu, paste0(parquet_processed, "/Shingrix_postFlu.rds"))


    


