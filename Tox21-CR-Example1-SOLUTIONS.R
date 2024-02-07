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
  labs(y= "AC50", x="Chemical Name")
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

# get coefficients
model_AUR_BPB_model_wCI <- tidy(model_AUR_BPB, conf.int = TRUE)
# add 95% CI
model_AUR_BPB_model_wCI$CI95 <- model_AUR_BPB_model_wCI$std.error *1.96
model_AUR_BPB_model_wCI

# Question 1: Compare your drc predicted AC50 value to the tcpl derived AC50 value? 
# Are the values the same? How can you tell?

#Answer: While the two AC50 values are different(tcpl:1.47, drm:1.33 ), the tcpl derived AC50 is within the 
#95% confidence interval (CI:0.884-1.78) of our drm derived AC50 so we can say they are statistically different from one another 

##########################################################################################
#### Example 2 ####

# Estimate the AC50 value for Bisphenol B in the CCTE_Simmons_GUA_TPO_dn using sample code in the example above
# parse data for exercise - e.g. chemical with the Bisphenol B
thyroid_GUA_BPB <- subset(thyroid_df, chemName == "Bisphenol B" & 
                            assay == "CCTE_Simmons_GUA_TPO_dn")

#isolate the dose-response
thyroid_dr_GUA_BPB <-  thyroid_GUA_BPB$Concentration.Response..Index.Concentration.Response.

# Remove parentheses and split by commas
ranges_GUA_BPB <- strsplit(gsub("\\(|\\)", "", thyroid_dr_GUA_BPB), ",")

# Initialize lists to store results
table_data_GUA_BPB <- list()

# Loop through each range and extract values
for (range in ranges_GUA_BPB[[1]]) {
  values_GUA_BPB <- as.numeric(unlist(strsplit(range, ":")))
  table_data_GUA_BPB <- c(table_data_GUA_BPB, values_GUA_BPB)
}

# Convert the list to a matrix
num_cols_GUA_BPB <- length(values_GUA_BPB)
table_matrix_GUA_BPB <- matrix(unlist(table_data_GUA_BPB), ncol = num_cols_GUA_BPB, byrow = TRUE)
# Convert the matrix to a data frame
result_GUA_BPB <- as.data.frame(table_matrix_GUA_BPB)
# Set column names
colnames(result_GUA_BPB) <- c("rep", "concentration", "response")
result_GUA_BPB

##### Dose Response Modelling ####
# fit for mixtures modeling
model_GUA_BPB<- drm(response~concentration, data=result_GUA_BPB,
                    fct=LL.4(fixed=c(NA, 0,  NA, NA), 
                             names = c("Slope", "Lower Limit", "Upper Limit", "AC50")))
summary(model_GUA_BPB)
# Quick Plot of curve fits
plot(model_GUA_BPB)

# get coefficients
model_GUA_BPB_model_wCI <- tidy(model_GUA_BPB, conf.int = TRUE)
# add 95% CI
model_GUA_BPB_model_wCI$CI95 <- model_GUA_BPB_model_wCI$std.error *1.96
model_GUA_BPB_model_wCI

# Question 2: Compare your drc predicted AC50 value to the tcpl derived AC50 value? 
# Are the values the same? How can you tell?

#Answer: this model is not a good fit according to DRM (high p-value on the model paramters) - looks more exponential than log-logistic
# the DRM AC50  is different the tcpl AC50 - tcpl has a couple pre-processing steps that combining the 
# replicates which many explain the difference we observe

# Question 3: Compare the AC50 measures between the CCTE_Simmons_GUA_TPO_dn
# and the CCTE_Simmons_GUA_TPO_dn assays. Both assays inform on similar the same molecular target.
# Which AC50 value would you use for conducting a risk assessment and why?

# Answer: Use the lowest AC50, it is the most potent, and therefore the most conservative (health protective) estimate
# of chemical hazard

