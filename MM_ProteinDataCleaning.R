### Cleaning the Protein Data Set ###

# libraries
library(tidyverse)
library(DataExplorer)
library(zoo)
library(gridExtra)
library(car)
library(missForest)

# Load in data and create one tibble
protein.train <- read_csv("ProteinTrain.csv")
protein.test <- read_csv("ProteinTest.csv")

# Merge the training and test data
protein <- bind_rows(train = protein.train, test = protein.test, .id = "Set")

# Make the Amino.Acid and Response variables factors
protein <- protein %>% mutate_at(vars(Amino.Acid, Response), factor)

# Summary
summary(protein)

# Plot missing values
plot_missing(protein)
# Consensus and PSSM have NAs

# Table of Amino Acid and Response
addmargins(table(protein$Amino.Acid, protein$Response))

# Use random forests to fill NA values
filled_nas <- protein %>% 
  select(-c(Response, Set)) %>% 
  as.data.frame() %>% 
  missForest()

# Merge Response and Set back with other variables
protein.clean <- filled_nas$ximp %>%  cbind(., protein %>% select(Set, Response)) %>% select(Set, SiteNum, Response, Amino.Acid, Iupred.score, ANN, PSSM, SVM, Consensus, normalization)

# Correlation matrix
plot_correlation(protein.clean, 
                 type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

# Confirming that all NAs have been dealt with
plot_missing(protein.clean)

# Create a .csv file of the cleaned data
protein.clean %>% write_csv("MM_ProteinCleaned.csv")
