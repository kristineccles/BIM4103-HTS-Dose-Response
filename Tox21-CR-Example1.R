##########################################################################################
# Sample code for assessing Toxcast/Tox21 Concentration-Response curves
# Example - CCTE_Simmons_AUR_TPO_dn
# Retrieved from: https://ice.ntp.niehs.nih.gov/Tools?tool=curvesurfer
# By: Kristin Eccles
# Written February 4th, 2024
##########################################################################################

#load libraries
install.packages("stringr") #only need to install once
install.packages("drc")
install.packages("broom")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")

library(stringr) #must be loaded every time
library(drc)
library(broom)
library(dplyr)
library(tidyr)
library(ggplot2)

# load data
thyroid_df <- read.csv("tox21_thyroid.csv")

##########################################################################################
#### Data Exploration ####
# Explore the TPO assays data - plot example
p1 <- ggplot(data = subset(thyroid_df, call == "Active") , 
             aes(x=chemName, y=log10(ac50), color = assay))+
  geom_point()+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_minimal()+
  theme(text=element_text(size=8))+
  labs(y= "Log10 AC50", x="Chemical Name")
p1

#### Example 1 ####
#### Prepare Data ####
# parse data for exercise - e.g. chemical with the Bisphenol B
thyroid_AUR_BPB <- subset(thyroid_df, chemName == "Bisphenol B" & 
                            assay == "CCTE_Simmons_AUR_TPO_dn")

#isolate the dose-response
thyroid_dr_AUR_BPB <-  thyroid_AUR_BPB$Concentration.Response..Index.Concentration.Response.

# Remove parentheses and split by commas
ranges_AUR_BPB <- strsplit(gsub("\\(|\\)", "", thyroid_dr_AUR_BPB), ",")

# Initialize lists to store results
table_data_AUR_BPB <- list()

# Loop through each range and extract values
for (range in ranges_AUR_BPB[[1]]) {
  values_AUR_BPB <- as.numeric(unlist(strsplit(range, ":")))
  table_data_AUR_BPB <- c(table_data_AUR_BPB, values_AUR_BPB)
}

# Convert the list to a matrix
num_cols_AUR_BPB <- length(values_AUR_BPB)
table_matrix_AUR_BPB <- matrix(unlist(table_data_AUR_BPB), ncol = num_cols_AUR_BPB, byrow = TRUE)
# Convert the matrix to a data frame
result_AUR_BPB <- as.data.frame(table_matrix_AUR_BPB)
# Set column names
colnames(result_AUR_BPB) <- c("rep", "concentration", "response")
result_AUR_BPB

##### Dose Response Modelling ####
# fit for mixtures modeling
model_AUR_BPB<- drm(response~concentration, data=result_AUR_BPB,
                       fct=LL.4(fixed=c(NA, 0 , NA, NA), 
                                names = c("Slope", "Lower Limit", "Upper Limit", "AC50")))
summary(model_AUR_BPB)
# Quick Plot of curve fits
plot(model_AUR_BPB)

# get coefficients and 95% CI
model_AUR_BPB_model_wCI <- tidy(model_AUR_BPB, conf.int = TRUE)

# Question 1: Compare your drc predicted AC50 value to the tcpl derived AC50 value? 
# Are the values the same? How can you tell?

##########################################################################################
#### Example 2 ####

# Estimate the AC50 value for Bisphenol B in the CCTE_Simmons_GUA_TPO_dn using sample code in the example above




# Question 2: Compare your drc predicted AC50 value to the tcpl derived AC50 value? 
# Are the values the same? How can you tell?

# Question 3: Compare the AC50 measures between the CCTE_Simmons_GUA_TPO_dn
# and the CCTE_Simmons_GUA_TPO_dn assays. Both assays inform on similar the same molecular target (you can find more information
# about each assay on the US EPA CompTox Dashboard under the Assay/Gene tab: https://comptox.epa.gov/dashboard/)
# Which AC50 value would you use for conducting a risk assessment and why?

