################################################################################
# Author: Eleanor Barry
# Date: 06/11/2024
# Version: R 4.3.0
# File name: 10b_Shingrix_Covariates.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## PATID-filtered.rds
## filtered_med.parquet
## combined_med
## [almost all codelists]
## Ethnicity16
## CKD
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
## Full_data
# Description of file: Finds all patients covariates when we consider covariates to be time-varying
# Actions:
## Subset medical/drug records to only be past 1st Jan 2011
## Link relevant patient ids to their medical records
################################################################################


    ############################## Load in data  ###############################
    
    file_paths <- sprintf(med_combined, 1:25)
    
    combined_med <- file_paths %>%
      map_dfr(~ read_parquet(.x))
    combined_med <- combined_med %>%
      filter(obsdate>="2011-01-01")
    
    file_paths <- sprintf(drug_combined, 1:20)
    combined_drug <- file_paths %>%
      map_dfr(~ read_parquet(.x))
    combined_drug <- combined_drug %>%
      filter(issuedate>="2011-01-01")
    
    
    Shingrix_full=read_rds(paste0(parquet_processed, "/Shingrix_postFlu.rds"))
    
    
    
    ########################## Binary outcomes #################################
    
    
    #01 Alcohol
    
    codelist_alcohol_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_alcohol_aurum.txt"), 
                                         delim = "\t", escape_double = FALSE, 
                                         col_types = cols(medcodeid = col_character(), 
                                                          snomedctconceptid = col_character(), 
                                                          snomedctdescriptionid = col_character()), 
                                         trim_ws = TRUE)
    
    codelist_alcohol_product_aurum <- read_delim(paste0(codelists, "Prodcodes/codelist_alcohol_product_aurum.txt"), 
                                                 delim = "\t", escape_double = FALSE, 
                                                 col_types = cols(prodcodeid = col_character(), 
                                                                  dmdid = col_character(), termfromemis = col_character(), 
                                                                  formulation = col_character()),  
                                                 trim_ws = TRUE)
    
    Alcohol <- combined_med %>%
      filter(medcodeid %in% codelist_alcohol_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end)
    
    Alcohol_prod <- combined_drug %>%
      filter(prodcodeid %in% codelist_alcohol_product_aurum$prodcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, prodcodeid, issuedate) %>%
      left_join( Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(issuedate <= max_end)
    
    Shingrix_full$Alcohol=ifelse(Shingrix_full$patid %in% Alcohol$patid | Shingrix_full$patid %in% Alcohol_prod$patid, 1,0)
    rm(codelist_alcohol_aurum, codelist_alcohol_product_aurum, Alcohol, Alcohol_prod)
    
    #02 Smoking 
    names(Shingrix_full)[names( Shingrix_full) == "Smoking"] <- "Smoking_old"
    codelist_smoking_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_smoking_aurum.txt"), 
                                         delim = "\t", escape_double = FALSE, 
                                         col_types = cols(medcodeid = col_character(), 
                                                          smokstatus = col_character(), status = col_character()), 
                                         trim_ws = TRUE)
    Smoking=subset(combined_med, combined_med$medcodeid %in% codelist_smoking_aurum$medcodeid & combined_med$patid %in% Shingrix_full$patid)
    Smoking$obsdate[is.na(Smoking$obsdate)]=Smoking$enterdate[is.na(Smoking$obsdate)]

    Smoking2 <- Shingrix_full[c("patid", "max_end")] %>%
      left_join(Smoking, by = "patid") %>%
      filter(obsdate <= max_end) %>%
      left_join(codelist_smoking_aurum, by = "medcodeid") %>%
      mutate(smoke = ifelse(status=='smoker or ex-smoker', 1, 0)) %>%
      group_by(patid) %>%
      mutate(smoke = max(smoke)) %>%
      ungroup() %>%
      select(patid, smoke) %>%  # Keep only the 'patid' and 'smoke' columns
      distinct()  # Remove duplicates
    
    Shingrix_full <- Shingrix_full %>%
      left_join(Smoking2, by="patid") %>%
      rename(Smoking = smoke)
       
    rm(Smoking, Smoking2, codelist_smoking_aurum)
    
    # 03 Region/Deprivation
    
    # NA for now, to be filled in
    
    # 04 CMD
    
    codelist_cmd_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_cmd_aurum.txt"), 
                                     delim = "\t", escape_double = FALSE, 
                                     col_types = cols(medcodeid = col_character(), 
                                                      observations = col_character(), snomedctconceptid = col_character(), 
                                                      snomedctdescriptionid = col_character()), 
                                     trim_ws = TRUE)
    
    cmd <- combined_med %>%
      filter(medcodeid %in% codelist_cmd_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1942-01-01") )
    
    Shingrix_full$CMD=ifelse(Shingrix_full$patid %in% cmd$patid, 1,0)
    rm(codelist_cmd_aurum, cmd)

    # 05 SMI
    
    codelist_smi_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_smi_aurum.txt"), 
                                     delim = "\t", escape_double = FALSE, 
                                     col_types = cols(medcodeid = col_character(), 
                                                      snomedctconceptid = col_character(), 
                                                      snomedctdescriptionid = col_character()), 
                                     trim_ws = TRUE)
    
    smi <- combined_med %>%
      filter(medcodeid %in% codelist_smi_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1941-01-01") )
    
    Shingrix_full$SMI=ifelse(Shingrix_full$patid %in% smi$patid , 1,0)
    rm(codelist_smi_aurum, smi)
    
    # Living alone
    
    codelist_livingalone_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_livingalone_aurum.txt"), 
                                             delim = "\t", escape_double = FALSE, 
                                             col_types = cols(medcodeid = col_character(), 
                                                              snomedctconceptid = col_character(), 
                                                              snomedctdescriptionid = col_character(), 
                                                              emiscodecategoryid = col_character()), 
                                             trim_ws = TRUE)
    
    Live_alone <- combined_med %>%
      filter(medcodeid %in% codelist_livingalone_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1960-01-01") ) # Since 18
    
    Shingrix_full$LivingAlone=ifelse(Shingrix_full$patid %in% Live_alone $patid , 1,0)
    rm(codelist_livingalone_aurum, Live_alone)
    
    #Carehome
    codelist_carehome_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_carehome_aurum.txt"), 
                                          delim = "\t", escape_double = FALSE, 
                                          col_types = cols(snomedctconceptid = col_character(), 
                                                           snomedctdescriptionid = col_character(), 
                                                           emiscodecategoryid = col_character(), 
                                                           medcodeid = col_character()), trim_ws = TRUE)
    
      Care_home <- combined_med %>%
      filter(medcodeid %in% codelist_carehome_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1960-01-01") ) #Since 18
    
    Shingrix_full$CareHome=ifelse(Shingrix_full$patid %in% Care_home$patid , 1,0)
    rm(codelist_carehome_aurum, Care_home)
    
    
    # Dementia
    codelist_dementia_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_dementia_aurum.txt"), 
                                             delim = "\t", escape_double = FALSE, 
                                             col_types = cols(medcodeid = col_character(), 
                                                              snomedctconceptid = col_character(), 
                                                              snomedctdescriptionid = col_character(), 
                                                              emiscodecategoryid = col_character()), 
                                             trim_ws = TRUE)
    
      Dementia <- combined_med %>%
      filter(medcodeid %in% codelist_dementia_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1960-01-01") ) #Since 18
    
    Shingrix_full$Dementia=ifelse(Shingrix_full$patid %in% Dementia$patid , 1,0)
    rm(codelist_dementia_aurum, Dementia)
    
    #DM
    
    codelist_dm_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_dm_aurum.txt"), 
                                          delim = "\t", escape_double = FALSE, 
                                          col_types = cols(medcodeid = col_character(), 
                                                           snomedctconceptid = col_character(), 
                                                           snomedctdescriptionid = col_character(), 
                                                           emiscodecategoryid = col_character()), 
                                          trim_ws = TRUE)
    
    dm <- combined_med %>%
      filter(medcodeid %in% codelist_dm_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1941-01-01") ) # Since 18
    
    Shingrix_full$DM=ifelse(Shingrix_full$patid %in% dm$patid , 1,0)
    rm(codelist_dm_aurum, dm)
    
    # COPD
    codelist_copd_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_copd_aurum.txt"), 
                                      delim = "\t", escape_double = FALSE, 
                                      col_types = cols(medcodeid = col_character(), 
                                                       incident = col_character(), prevalent = col_character(), 
                                                       snomedctconceptid = col_character(), 
                                                       snomedctdescriptionid = col_character(), 
                                                       emiscodecategoryid = col_character()), 
                                      trim_ws = TRUE)
    
    copd <- combined_med %>%
      filter(medcodeid %in% codelist_copd_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, max_end), by = "patid") %>%
      filter(obsdate <= max_end) %>%
      filter(obsdate >= as.Date("1942-01-01") )
    
    Shingrix_full$COPD=ifelse(Shingrix_full$patid %in% copd$patid , 1,0)
    rm(codelist_copd_aurum, copd)
    

    # CKD
    names(Shingrix_full)[names( Shingrix_full) == "CKD"] <- "CKD_old"
    ckd <- read_rds(paste0(processed_data, "CKD_timevary.rds"))
    Shingrix_full=left_join(Shingrix_full, ckd, by="patid")
    rm(ckd)
    Shingrix_full$CKD[is.na(Shingrix_full$CKD)]=0

    #BMI
    names(Shingrix_full)[names( Shingrix_full) == "BMI"] <- "_old"
    bmi <- read_rds(paste0(processed_data, "BMI_timevary.rds"))
    Shingrix_full=left_join(Shingrix_full, bmi, by="patid")
    rm(bmi)
    
    # Ethnicity
    names(Shingrix_full)[names( Shingrix_full) == "ethnicity5"] <- "_ethnicity5_old"
    names(Shingrix_full)[names( Shingrix_full) == "ethnicity16"] <- "_ethnicity16_old"
    ethnicity5 <- read_rds(paste0(processed_data, "Ethnicity.rds"))
    colnames(ethnicity5)[colnames(ethnicity5) == "Ethnicity"] <- "ethnicity5"
    ethnicity16 <- read_rds(paste0(processed_data, "Ethnicity16.rds"))
    colnames(ethnicity16)[colnames(ethnicity16) == "Ethnicity"] <- "ethnicity16"
    Shingrix_full=left_join(Shingrix_full, ethnicity5, by="patid")
    Shingrix_full=left_join(Shingrix_full, ethnicity16, by="patid")
    rm(ethnicity5, ethnicity16 )
    
    #Zoster 
    #We will consequently adjust for recent shingles in our regression analyses of potential
    #factors associated with vaccine uptake to account for recent shingles affecting the choice to
    #vaccinate or not (e.g., because a clinician does not think it is worthwhile due to the natural
                      #boosting effect of shingles, or because an individual was offered vaccination while still
                      #experiencing the rash of shingles) and the potential for previous shingles to be related to the
    #various potential factors associated with vaccine uptake. We will define recent shingles as any record indicating shingles in the year before individuals become eligible for Shingrix® vaccination
    #(i.e., the latest of: study start [1st September 2021], when Shingrix® became available in the UK;
     # 70th birthday; or the date of meeting our Shingrix®-eligible immunosuppression definition).
    
    codelist_zoster_aurum <- read_delim(paste0(codelists, "Medcodes/codelist_zoster_aurum.txt"), 
                                        delim = "\t", escape_double = FALSE, 
                                        col_types = cols(medcodeid = col_character(), 
                                                         marker = col_character(), site = col_character(), 
                                                         hosp = col_character(), originalreadcode = col_character(), 
                                                         snomedctconceptid = col_character(), 
                                                         snomedctdescriptionid = col_character(), 
                                                         emiscodecategoryid = col_character(), 
                                                         term = col_character()), trim_ws = TRUE)
    
    zoster <- combined_med %>%
      filter(medcodeid %in% codelist_zoster_aurum$medcodeid & patid %in% Shingrix_full$patid) %>%
      select(patid, medcodeid, obsdate) %>%
      left_join(Shingrix_full %>% select(patid, min_entry), by = "patid") %>%
      filter(obsdate <= min_entry) %>%
      filter(obsdate >= as.Date("1941-01-01") ) %>%
      filter(obsdate >= as.Date(min_entry-365.25) &  obsdate <= as.Date(min_entry))
   
     Shingrix_full$Zoster=ifelse(Shingrix_full$patid %in% zoster$patid , 1,0)
    
     # IMD
     patient_imd <- read_delim(paste0(imd_linked , "/patient_2019_imd_Shingrix.txt"), 
                                              delim = "\t", escape_double = FALSE, 
                                              col_types = cols(patid = col_character(), 
                                              pracid = col_character(), e2019_imd_5 = col_character()), 
                                              trim_ws = TRUE)
     
     practice_imd <- read_delim(paste0(imd_linked, "/practice_imd_Shingrix.txt"), 
                                         delim = "\t", escape_double = FALSE, 
                                col_types = cols(pracid = col_character()),
                                         trim_ws = TRUE)
     names(Shingrix_full)[names( Shingrix_full) == "pracid.y"] <- "pracid"
     names(Shingrix_full)[names( Shingrix_full) == "IMD"] <- "IMD_old"
     Shingrix_IMD = left_join(Shingrix_full[c("patid", "pracid")], patient_imd, by=c("patid", "pracid")) # 142421/174918 (81.4%)
     names(Shingrix_IMD)[names(Shingrix_IMD) == "e2019_imd_5"] <- "PatientIMD"
     Shingrix_IMD = left_join(Shingrix_IMD, practice_imd[c("pracid", "e2019_imd_5")], by=c("pracid")) # 142421/174918 (81.4%)
     names(Shingrix_IMD)[names(Shingrix_IMD) == "e2019_imd_5"] <- "PracticeIMD"
     Shingrix_IMD$IMD=Shingrix_IMD$PatientIMD
     Shingrix_IMD$IMD[is.na(Shingrix_IMD$IMD)]=Shingrix_IMD$PracticeIMD[is.na(Shingrix_IMD$IMD)]
     Shingrix_full=left_join(Shingrix_full, Shingrix_IMD, by="patid")
     
    write_rds(Shingrix_full, paste0(processed_data, "Shingrix_TimeVar.rds" ))

    


