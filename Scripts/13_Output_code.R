################################################################################
# Author: Eleanor Barry
# Date: 06/01/2024
# Version: R 4.3.0
# File name: 13_Output_code.R
# Status: Complete
# CPRD version: March 2024
# Data sets used: 
## Shingrix_postFlu
# R scripts needed: 
## 00_Set_up_directory.R
# Data sets created: 
## Full_data
# Description of file: This code is what was used for the output data for this study.
# It includes the following: 
# 1. Baseline summary statistics, 
# 2. OR + 95% CI, 
# 3. Entry plot (Supp Fig 1)
# 4. Risk changes
# 5. IMD-ethnicity interaction
# 6. Flu uptake
## For sensitivity analysis, the following will need to be changed in the following ways
# Static cohort – Cohort reduced to those in the study on 1 Sept 2021. Command run from console, analysis repeated from 13_Output_code.R, line 65.
# Time-updated covariates – Time-varying covariates established in 10b_Shingrix_covariates_timevarying.R (alternative to 10_Shingrix_covariates.R), then subsequent scripts re-run.
# Restricted to 13 months – See line 45 in 13_Output_code.R (restricts to this follow-up time).
# Restricted to age 70 only – Cohort reduced to those aged 70 on entry. Command run from console, analysis repeated from 13_Output_code.R, line 65.
# Sixteen categories of ethnicity – Adjusted in 10_Shingrix_covariates.R and 13_Output_code.R, line 50.
# Five immunosuppression risk categories – Formulated in 06_Shingrix_risk_periods.R and adjusted in 13_Output_code.R as "risk_category" rather than "risk level".
################################################################################



#################### Sort data #################################################

df= read_rds(paste0(parquet_processed, "/Shingrix_postFlu.rds"))
df$pracid.x=df$pracid
df$Sex=ifelse(df$gender==1, "Male", "Female")
df$`Year of birth`=as.character(df$yob)
Region <- read_delim("[filepath/region.txt]", 
    delim = "\t", escape_double = FALSE, 
    trim_ws = TRUE)
Region= Region %>% 
  rename(Region=Description, region =regionid)
df=left_join(df, Region, by=c("region"))
df$`Risk level`=ifelse(df$risk_condition=="high", "High", "Low")
df$`Risk category`=df$risk_category
df$Age=floor(as.numeric(difftime(df$min_entry, df$dob, units = "days")) / 365.25)
df$Age[df$Age==69]=70 # Due to /365.25, some have age 69.998 which will round to 70, this corrects this 
df$Age=as.character(df$Age)
df$Ethnicity=df$ethnicity5
df$Vaccine=ifelse( !is.na(df$first_vacc_type) & (df$first_vacc_type %in% c("Shingrix", "Neutral")),1,0)
df$Vaccine2=ifelse(!is.na(df$second_vacc_type) & (df$second_vacc_type %in% c("Shingrix", "Neutral")),1,0)
df$Region=relevel(factor(df$Region), ref="London")
df$IMD=df$PatientIMD
df$time_in_study=as.numeric(df$max_end-df$min_entry)
df2 <- subset(df, interval(df$first_vacc_date, df$max_cohort_exit) %/% months(1) >= 13 & df$Vaccine==1 )

### The following "create summary" has two options, either (1) counts and percentages or (2) median time and IQR.
### You will need to comment/uncomment which you want. 

## 1
create_summary <- function(data, group_var, filter_var = NULL, filter_val = NULL) {
  if (!is.null(filter_var)) {
    data <- data %>% filter(get(filter_var) == filter_val)
  }

  data %>%
    group_by(!!sym(group_var)) %>%  # Use !!sym() to handle variable names correctly
    summarise(
      count = n(),
      proportion = count / nrow(data) * 100
    ) %>%
    mutate(
      # Change item names (categories) to "Yes" for 1 and "No" for 0
      item = case_when(
        !!sym(group_var) == 1 ~ "Yes",
        !!sym(group_var) == 0 ~ "No",
        TRUE ~ as.character(!!sym(group_var))  # Keep other values unchanged
      ),
      count_proportion = paste0(count, " (", format(round(proportion, 1), nsmall = 1), "%)")
    ) %>%
    select(item, count_proportion)  # Select the necessary columns
}

## 2 
create_summary <- function(data, group_var, filter_var = NULL, filter_val = NULL) {
  if (!is.null(filter_var)) {
    data <- data %>% filter(get(filter_var) == filter_val)
  }
  
  # Ensure that time_in_study is numeric (and remove any non-numeric values if they exist)
  # Group by the selected grouping variable (e.g., "Age")
  data %>%
    group_by(!!sym(group_var)) %>%  # Group by the variable (e.g., "Age")
    summarise(
      median_time = median(time_in_study, na.rm = TRUE) ,  # Median of "time_in_study"
      lower_IQR = quantile(time_in_study, 0.25, na.rm = TRUE),  # Lower IQR (25th percentile)
      upper_IQR = quantile(time_in_study, 0.75, na.rm = TRUE)  # Upper IQR (75th percentile)
    ) %>%
    mutate(
      # Create the summary format as "median [lower IQR - upper IQR]"
      time_summary = paste0(
        as.numeric(format(round(median_time, 0))), " [", 
        as.numeric(format(round(lower_IQR, 0))), " - ", 
        as.numeric(format(round(upper_IQR, 0))), "]"
      ),
      item = as.character(!!sym(group_var))  # Preserve the group value as "item" (e.g., "70" for Age)
    ) %>%
    select(item, time_summary)  # Select the necessary columns: "item" (group) and the summary
}

# Overall summary
sex_summary <- create_summary(df, "Sex")
age_summary <- create_summary(df, "Age")
prior_zoster <- create_summary(df, "Zoster")
region_summary <- create_summary(df, "Region")
imd_summary <- create_summary(df, "IMD")
#ethnicity_summary <- create_summary(df, "ethnicity")
zoster_summary <- create_summary(df, "Zoster")
alcohol_summary <- create_summary(df, "Alcohol")
smoking_summary <- create_summary(df, "Smoking")
bmi_summary <- create_summary(df, "BMI")
anxdep_summary <- create_summary(df, "CMD")
smi_summary <- create_summary(df, "SMI")
dementia_summary <- create_summary(df, "Dementia")
dm_summary <- create_summary(df, "DM")
copd_summary <- create_summary(df, "COPD")
ckd_summary <- create_summary(df, "CKD")
livingalone_summary <- create_summary(df, "LivingAlone")
carehome_summary <- create_summary(df, "CareHome")
vaccine_summary <- create_summary(df, "Vaccine")
ethnicity_summary <- create_summary(df, "Ethnicity")
risklevel_summary <- create_summary(df, "Risk level")
riskcategory_summary <- create_summary(df, "Risk category")


# Combine summaries
combine_summaries <- function(overall, subset1) {
  overall %>%
    left_join(subset1, by = "item", suffix = c("", "_vaccine"))
}


final_table <- bind_rows(age_summary,sex_summary,  zoster_summary, region_summary, imd_summary,  ethnicity_summary,  risklevel_summary, riskcategory_summary, anxdep_summary, smi_summary, dementia_summary, dm_summary, copd_summary, ckd_summary,  alcohol_summary, smoking_summary, bmi_summary, livingalone_summary, carehome_summary) %>%
  mutate(category = rep(c("Age","Sex",  "Zoster", "Region", "IMD", "Ethnicity",  "Risk level",  "Risk category", "CMD", "SMI", "Dementia","DM", "COPD", "CKD", "Alcohol", "Smoking", "BMI", "LivingAlone" , "CareHome"), c(nrow(age_summary),nrow(sex_summary),  nrow(zoster_summary), nrow(region_summary),  nrow(imd_summary) , nrow(ethnicity_summary) , nrow(risklevel_summary) ,nrow(riskcategory_summary) , nrow(anxdep_summary), nrow(smi_summary), nrow(dementia_summary), nrow(dm_summary), nrow(copd_summary), nrow(ckd_summary),nrow(alcohol_summary), nrow(smoking_summary), nrow(bmi_summary), nrow(livingalone_summary), nrow(carehome_summary)  )))

final_table=as.data.table(final_table)

final_table=final_table[, category := ifelse(duplicated(category), "", category)]

final_table=as.data.frame(final_table)

# Generate the table (select first for counts, second for IQR)
## 1
final_table %>%
  select(category, item, count_proportion) %>%
  kable(col.names = c("", "", "Primary analysis N = 86197 from 1668 practices,
Median duration (days) of registration (IQR): 348 (IQR 104–729)")) %>%
  kable_styling(full_width = FALSE) %>%
  collapse_rows(columns = 1, valign = "top")

#df=df_backup
#df <- subset(df, interval(df$min_entry, df$max_cohort_exit) %/% months(1) >= 13 )
#df=subset(df, df$Vaccine==0 & df$Vaccine2==0)

## 2
final_table %>%
  select(category, item, time_summary) %>%
  kable(col.names = c("", "", "Primary analysis N = 86197 from 1668 practices,
Median duration (days) of registration (IQR): 348 (IQR 104–729)")) %>%
  kable_styling(full_width = FALSE) %>%
  collapse_rows(columns = 1, valign = "top")

nrow(df)
# Practices 
length(unique(df$pracid))
#

median(df$max_end-df$min_entry)
quantile(df$max_end-df$min_entry)
# All 346 [104-729]
# Vaccinated 

############## 2. ORs and size of data for each model #############################





model1 = glm(Vaccine ~ Age + Sex, family="binomial", data=df)
model2 = glm(Vaccine ~ Age + Sex + Zoster , family="binomial", data=df)
model3 = glm(Vaccine ~ Age + Sex + Zoster + Region + IMD + Ethnicity, family="binomial", data=df)
model4 = glm(Vaccine ~ Age + Sex + Zoster + Region + IMD + Ethnicity + as.factor(`Risk level`), family="binomial", data=df)
model5 = glm(Vaccine ~ Age + Sex + Zoster + Region + IMD + Ethnicity + as.factor(`Risk level`) + CMD + SMI + Dementia + DM + COPD + CKD, family="binomial", data=df)
model6 = glm(Vaccine ~ Age + Sex + Zoster + Region + IMD + Ethnicity + as.factor(`Risk level`) + CMD + SMI + Dementia + DM + COPD + CKD + Alcohol + Smoking + BMI, family="binomial", data=df)
model7 = glm(Vaccine ~ Age + Sex + Zoster + Region + IMD + Ethnicity + as.factor(`Risk level`) + CMD + SMI + Dementia + DM + COPD + CKD + Alcohol + Smoking + BMI + LivingAlone + CareHome, family="binomial", data=df)

# Change for model
cluster_se <- cluster.vcov(model6, cluster = df$pracid.y)
#ci <- confint(model6, vcov. = cluster_se)
coefs <- coef(model6)
robust_se <- sqrt(diag(cluster_se))  # Robust standard errors

# Confidence intervals (CI) calculation
ci_lower <- exp(coefs - 1.96 * robust_se)  # Lower bound of 95% CI
ci_upper <- exp(coefs + 1.96 * robust_se)  # Upper bound of 95% CI

# Prepare the data for forest plot
forest_data <- data.frame(
  variable = names(coefs)[-1],  # Extract variable names
  OR = exp(coefs)[-1], 
  lower = ci_lower[-1], 
  upper = ci_upper[-1], 
  se = robust_se[-1]  # Robust standard errors
)


############################## 3.   Plot for entry dates ########################
data=df
data$month_year <- format(data$min_entry, "%Y-%m")
# OR
# 

# Count the occurrences of each month_year
data_summary <- data %>%
  group_by(month_year) %>%
  summarise(count = n())

# Create the ggplot histogram showing the count per month
ggplot(data_summary, aes(x = month_year, y = count)) +
  geom_bar(stat = "identity", col="blue") +
  labs(x = "Date", y = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##################### 4. Changing risk from max overall to at entry ###############

df$risk_condition_timevary="low"
df$risk_condition_timevary[!is.na(df$start_fup_1_combined) | !is.na(df$start_fup_2) | !is.na(df$start_fup_3)]="high"

df <- df %>%
  mutate(risk_category_timevary = case_when(
    !is.na(start_fup_1_combined) ~ 1,
    !is.na(start_fup_2) ~ 2,
    !is.na(start_fup_3) ~ 3,
    !is.na(start_fup_4_combined) ~ 4,
    !is.na(start_fup_5_combined) ~ 5,
  ))

################ 5. Interaction for IMD and ethnicity##############################

IMD_eth_table <- df %>%
  group_by(IMD, Ethnicity) %>%
  summarise(
    count = n(),
    uptake_rate = mean(Vaccine, na.rm = TRUE)  # Proportion of vaccinated individuals
  ) %>%
  arrange(IMD, Ethnicity)

IMD_eth_table<- IMD_eth_table%>%
  mutate(uptake_rate = sprintf("%.2f%%", uptake_rate))  # Format as percentage with 2dp

# Display the table in RStudio Viewer with enhanced styling
IMD_eth_table %>%
  kable(col.names = c("IMD Quintile", "Ethnicity", "Count", "Vaccine Uptake (%)"),
        caption = "Vaccine Uptake by Ethnicity within IMD Quintiles") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "responsive"), full_width = FALSE)

ggplot(heatmap_data, aes(x = factor(IMD), y = Ethnicity, fill = uptake_rate)) +
  geom_tile(color = "white") +  # Add white grid lines
  geom_text(aes(label = label), color = "black", size = 3) +  # Smaller text inside boxes
  scale_fill_gradient(low = "lightblue", high = "red", name = "Uptake (%)") +  # Color gradient
  labs( 
    x = "IMD Quintile", 
    y = "Ethnicity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
                
model6_interaction =glm(Vaccine ~ Age + Sex + Zoster + Region + IMD*Ethnicity + as.factor(`Risk level`) + CMD + SMI + Dementia + DM + COPD + CKD + Alcohol + Smoking + BMI, family="binomial", data=df)


# Perform Likelihood Ratio Test
lrt_result <- lrtest(model6, model6_interaction)

print(lrt_result)

## Find size of compete case of data

df2=subset(df, !is.na(df$Ethnicity) & !is.na(df$BMI) & !is.na(df$IMD) & !is.na(df$Smoking) )

################## 6. Flu uptake  #############################

# Flu - time betweeen vaccines

df_flu <- df %>%
  # Ensure the date columns are in Date format
  mutate(
    first_vacc_date = as.Date(first_vacc_date),
    first_flu_vacc_date = as.Date(first_flu_date),
    vacc_month = month(first_vacc_date)
  ) %>%
  # Filter for September (9) to March (3)
  filter(vacc_month %in% c(9, 10, 11, 12, 1, 2, 3)) %>%
  # Calculate the difference in days between dates
  mutate(
    date_diff = pmin(abs(as.numeric(difftime(first_flu_date, first_vacc_date, units = "days"))),abs(as.numeric(difftime(second_flu_date, first_vacc_date, units = "days"))), na.rm=TRUE),
    flu_timing = case_when(
      abs(date_diff)==0 ~ "Same day",
      abs(date_diff)!= 0 &  abs(date_diff) < 7 ~ "under a week",
      abs(date_diff) >= 7 & abs(date_diff) <= 14 ~ "7-14 days",
      abs(date_diff) > 14 ~ ">14 days",
      TRUE ~ NA_character_
    )
  )







