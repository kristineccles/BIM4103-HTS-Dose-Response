##########################################################################################
# Sample code for assessing Toxcast/Tox21 Concentration-Response curves
# Example - CCTE_Simmons_AUR_TPO_dn
# Retrieved from: https://ice.ntp.niehs.nih.gov/Tools?tool=curvesurfer
# By: Kristin Eccles
# Written February 4th, 2022
##########################################################################################

#load libraries
install.packages("stringr") #only need to install once
install.packages("drc")

library(stringr) #must be loaded every time
library(drc)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)

# load data
thyroid_df <- read.csv("tox21_thyroid.csv")
# keep only the needed columns in r
thyroid_df <- 

##########################################################################################
#### Prepare Data ####
# parse data for excercise - e.g. chemical with the Bisphenol A
thyroid_subset <- subset(thyroid_df, chemName == "Bisphenol A")

#isolate the dose-response
thyroid_dr <-  thyroid_subset$Concentration.Response..Index.Concentration.Response.

# Remove parentheses and split by commas
ranges <- strsplit(gsub("\\(|\\)", "", thyroid_dr), ",")

# Initialize lists to store results
table_data <- list()

# Loop through each range and extract values
for (range in ranges[[1]]) {
  values <- as.numeric(unlist(strsplit(range, ":")))
  table_data <- c(table_data, values)
}

# Convert the list to a matrix
num_cols <- length(values)
table_matrix <- matrix(unlist(table_data), ncol = num_cols, byrow = TRUE)
# Convert the matrix to a data frame
result_df <- as.data.frame(table_matrix)
# Set column names
colnames(result_df) <- c("rep", "concentration", "response")
result_df

##########################################################################################
##### Dose Response Modelling ####
# fit for mixtures modeling
individual_model<- drm(response~concentration, data=result_df,
                       fct=LL.4(fixed=c(NA, NA , NA, NA), 
                                names = c("Slope", "Lower Limit", "Upper Limit", "ED50")))
summary(individual_model)
# Quick Plot of curve fits
plot(individual_model)

# get coefficients
individual_model_wCI <- tidy(individual_model, conf.int = TRUE)
# add 95% CI
individual_model_wCI$CI95 <- individual_model_wCI$std.error *1.96
individual_model_wCI


