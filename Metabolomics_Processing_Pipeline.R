library(tidyverse)
library(dplyr)
library(readxl)

## Read in raw data from excel file, select sheet with 'raw' data
peak_area_data = read_excel("/home/bneely/shulmanRotation/metabolomicsWGSIntegration/BAYL-06-22MD+ DATA TABLES.XLSX",
                            sheet = "Peak Area Data")

## Batch normalize data and convert to dataframe
# The mutate_if function will operate on each column that is of the numeric type
# i.e. the metabolite columns
batch_normalized = peak_area_data %>%
  mutate_if(is.numeric, function(x) x / median(x, na.rm = T))

## Drop subjects/proteins with too many missing values
# Drop sample name for next step
batch_normalized_no_IDs = subset(batch_normalized, select = -c(PARENT_SAMPLE_NAME))
# Identify and remove proteins with >=30% missing; subjects with >=40% missing. So...
# Keep columns (proteins) with >70% of values not missing
batch_normalized_no_IDs_clean_1 = batch_normalized_no_IDs[,which(colMeans(!is.na(batch_normalized_no_IDs)) > 0.7)]
# Keep rows (subjects) with >60% of values not missing
batch_normalized_no_IDs_cleaned = batch_normalized_no_IDs_clean_1[which(rowMeans(!is.na(batch_normalized_no_IDs_clean_1)) > 0.6),]
batch_normalized_no_IDs_cleaned$PARENT_SAMPLE_NAME = batch_normalized$PARENT_SAMPLE_NAME

## Impute missing values using minimum
# Calculate the minimum value across each metabolite to be used for imputing missing values
batch_normalized_cleaned_imputed = batch_normalized_no_IDs_cleaned %>% 
  mutate_if(is.numeric, function(x) ifelse(is.na(x), min(x, na.rm = T), x))

## Log transform the data, but don't do anything to col 1227 which is the PARENT_SAMPLE_NAME col
batch_normalized_cleaned_imputed_log = batch_normalized_cleaned_imputed
batch_normalized_cleaned_imputed_log[-1227] = log(batch_normalized_cleaned_imputed[-1227])
# Rearrange column order so ID is first
batch_normalized_cleaned_imputed_log = batch_normalized_cleaned_imputed_log[,c(1227, 1:1226)]

## Perform z-scoring per metabolite - indicates no. of sd away from mean
batch_normalized_cleaned_imputed_log_zscored = batch_normalized_cleaned_imputed_log %>% 
  mutate_if(is.numeric, function(x) (x - mean(x)) / sd(x))
