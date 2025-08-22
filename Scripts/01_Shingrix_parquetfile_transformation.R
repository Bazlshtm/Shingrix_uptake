################################################################################
# Author: Eleanor Barry, Edward Parker
# Date: 01/09/2024
# Version: R 4.3.0
# File name: 01_Shingrix_parquetfile_transformation.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: Extract files from CPRD Aurum
# R scripts needed: 00_Shingrix_set_up_directory.R
# Data sets created: 
## combined_med 
## combined_drug
# Description of file: 
##As extract files come with all patient information, this 
## subsets it to relevant medcodes/prodcodes needed for the study
# Actions performed: 
## Compile all codelists used, subsetting by relevant med/prodcodes 
################################################################################



####################### Compile observation files###############################

file_list <- list.files(path = file.path(codelists, "Medcodes"), 
                        pattern = "codelist_.*_aurum.txt")
Major_codelist= read_delim(paste0(codelists, "Medcodes/", file_list[1]),  
                           delim = "\t", 
                           escape_double = FALSE, 
                           col_types = cols(.default = "c"), 
                           trim_ws = TRUE)

for (i in 2:length(file_list)){
  Major_codelist=bind_rows(Major_codelist, 
                    read_delim(paste0(codelists, "Medcodes/", file_list[i]),  
                               delim = "\t", escape_double = FALSE, 
                               col_types = cols(.default = "c"), 
                               trim_ws = TRUE))
}


# Filter chunk
filter_chunk <- function(chunk, pos) {
 filtered_chunk <- chunk %>%
    filter((medcodeid %in% Major_codelist$medcodeid) )
 return(filtered_chunk)
}

# Define the path to your directory


# Loop over the 24 zip files
for (zip_num in 1:24) {
  # Initialize an empty list to store data frames
  data_list <- list()
  
  # Loop through the 5 text files in each zip archive
  for (txt_num in 0:4) {
    # Calculate the file number
    file_num <- sprintf("%03d", (zip_num - 1) * 5 + txt_num + 1)
    
    # Construct the file paths
    zip_file <- paste0(raw_data_meds, med_extract_text_start, zip_num, ".zip")
    txt_file <- paste0(med_extract_text_start, zip_num, med_extrat_text_mid, file_num, ".txt")
    
    # Read the data
    d_filtered <- read_delim_chunked(
      unz(zip_file, txt_file), 
      delim = "\t", 
      escape_double = FALSE, 
      col_types = cols(patid = col_character(), consid = col_skip(), pracid = col_character(), obsid = col_skip(), obsdate = col_date(format = "%d/%m/%Y"), enterdate = col_date(format = "%d/%m/%Y"), staffid = col_skip(), parentobsid = col_skip(), medcodeid = col_character(), value = col_character(), numunitid = col_character(), obstypeid = col_skip(), numrangelow = col_character(), numrangehigh =col_character(), probobsid = col_skip()), trim_ws = TRUE,
      chunk_size = 10000,
      callback = DataFrameCallback$new(filter_chunk)
    )
    
    # Remove duplicates
    d_filtered <- d_filtered[!duplicated(d_filtered), ]
    
    # Add the filtered data to the list
    data_list[[length(data_list) + 1]] <- d_filtered
  }

  # Bind the data frames together
  combined_data <- bind_rows(data_list)
  
  # Save the combined data as a Parquet file
  write_parquet(combined_data, paste0(raw_data_meds, med_combined_file, zip_num, ".parquet"))
}

for (zip_num in 25:25) {
  # Initialize an empty list to store data frames
  data_list <- list()
  
  # Loop through the 1 text files in each zip archive
  for (txt_num in 0:0) {
    # Calculate the file number
    file_num <- sprintf("%03d", (zip_num - 1) * 5 + txt_num + 1)
    
    # Construct the file paths
    zip_file <- paste0(raw_data_meds, med_extract_text_start, zip_num, ".zip")
    txt_file <- paste0(med_extract_text_start, zip_num, med_extrat_text_mid, file_num, ".txt")
    
    # Read the data
    d_filtered <- read_delim_chunked(
      unz(zip_file, txt_file), 
      delim = "\t", 
      escape_double = FALSE, 
      col_types = cols(patid = col_character(), consid = col_skip(), pracid = col_character(), obsid = col_character(), obsdate = col_date(format = "%d/%m/%Y"), enterdate = col_date(format = "%d/%m/%Y"), staffid = col_skip(), parentobsid = col_skip(), medcodeid = col_character(), value = col_character(), numunitid = col_character(), obstypeid = col_skip(), numrangelow = col_skip(), numrangehigh = col_skip(), probobsid = col_skip()), trim_ws = TRUE,
      chunk_size = 10000,
      callback = DataFrameCallback$new(filter_chunk)
    )
    
    # Remove duplicates
    d_filtered <- d_filtered[!duplicated(d_filtered), ]
    
    # Add the filtered data to the list
    data_list[[length(data_list) + 1]] <- d_filtered
  }
  
  # Bind the data frames together
  combined_data <- bind_rows(data_list)
  
  # Save the combined data as a Parquet file
  write_parquet(combined_data, paste0(raw_data_meds, med_combined_file, zip_num, ".parquet"))
}

############################## Compile drugs files #############################

file_list <- list.files(path = file.path(codelists, "Prodcodes"), 
                        pattern = "codelist_.*_aurum.txt")
Major_codelist_prod= read_delim(paste0(codelists, "Prodcodes/", file_list[1]),  
                                delim = "\t", 
                                escape_double = FALSE, 
                                col_types = cols(.default = "c"), 
                                trim_ws = TRUE)

for (i in 2:length(file_list)){
  Major_codelist_prod=bind_rows(Major_codelist_prod,
                                read_delim(paste0(codelists, "Prodcodes/", file_list[i]),  
                                delim = "\t", 
                                escape_double = FALSE, 
                                col_types = cols(.default = "c"), 
                                trim_ws = TRUE))
}


filter_chunk <- function(chunk, pos) {
 filtered_chunk <- chunk %>%
    filter((prodcodeid %in% Major_codelist_prod$prodcodeid)  )
 return(filtered_chunk)
}
# Define the path to your directory

# Loop over the 20 zip files
for (zip_num in 1:19) {
  # Initialize an empty list to store data frames
  data_list <- list()
  # Loop through the 5 text files in each zip archive
  for (txt_num in 0:4) {
    # Calculate the file number
    file_num <- sprintf("%03d", (zip_num - 1) * 5 + txt_num + 1)
    
    # Construct the file paths
    zip_file <- paste0(raw_data_drugs, drug_extract_text_start, zip_num, ".zip")
    txt_file <- paste0(drug_extract_text_start, 
                       zip_num, 
                       drug_extract_text_mid, 
                       file_num, ".txt")
    
    # Read the data
    d_filtered <- read_delim_chunked(
      unz(zip_file, txt_file),
      delim = "\t", 
      escape_double = FALSE, 
      col_types = cols(
        patid = col_character(), 
        issueid = col_skip(), 
        pracid = col_character(), 
        probobsid = col_skip(), 
        drugrecid = col_skip(), 
        issuedate = col_date(format = "%d/%m/%Y"), 
        enterdate = col_date(format = "%d/%m/%Y"), 
        staffid = col_skip(), 
        prodcodeid = col_character(), 
        quantity = col_number(), 
        quantunitid = col_character(), 
        duration = col_character(), 
        estnhscost = col_skip()
      ), 
      trim_ws = TRUE,
      chunk_size = 10000,
      callback = DataFrameCallback$new(filter_chunk)
    )
    
    # Remove duplicates
    d_filtered <- d_filtered[!duplicated(d_filtered), ]
    
    # Add the filtered data to the list
    data_list[[length(data_list) + 1]] <- d_filtered
  }

  # Bind the data frames together
  combined_data <- bind_rows(data_list)
  
  # Save the combined data as a Parquet file
  output_file <- paste0(raw_data_drugs, drug_combined_file, zip_num, ".parquet")
  write_parquet(combined_data, output_file)
}

# For last file
for (zip_num in 20:20) {
  # Initialize an empty list to store data frames
  data_list <- list()
  
  # Loop through the 4 text files in each zip archive
  for (txt_num in 0:3) {
    file_num <- sprintf("%03d", (zip_num - 1) * 5 + txt_num + 1)
    
    # Construct the file paths
    zip_file <- paste0(raw_data_drugs, drug_extract_text_start, zip_num, ".zip")
    txt_file <- paste0(drug_extract_text_start, 
                       zip_num, 
                       drug_extract_text_mid, 
                       file_num, ".txt")
    
    # Read the data
    d_filtered <- read_delim_chunked(
      unz(zip_file, txt_file),
      delim = "\t", 
      escape_double = FALSE, 
      col_types = cols(
        patid = col_character(), 
        issueid = col_skip(), 
        pracid = col_character(), 
        probobsid = col_skip(), 
        drugrecid = col_skip(), 
        issuedate = col_date(format = "%d/%m/%Y"), 
        enterdate = col_date(format = "%d/%m/%Y"), 
        staffid = col_skip(), 
        prodcodeid = col_character(), 
        quantity = col_number(), 
        quantunitid = col_character(), 
        duration = col_character(), 
        estnhscost = col_skip()
      ), 
      trim_ws = TRUE,
      chunk_size = 10000,
      callback = DataFrameCallback$new(filter_chunk)
    )
    
    # Remove duplicates
    d_filtered <- d_filtered[!duplicated(d_filtered), ]
    
    # Add the filtered data to the list
    data_list[[length(data_list) + 1]] <- d_filtered
  }
  
  # Bind the data frames together
  combined_data <- bind_rows(data_list)
  
  # Save the combined data as a Parquet file
  output_file <- paste0(raw_data_drugs, drug_combined_file, zip_num, ".parquet")
  write_parquet(combined_data, output_file)
}
