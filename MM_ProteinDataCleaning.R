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

# Correlation matrix
plot_correlation(protein, 
                 type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

PSSM.lm <- lm(PSSM ~ Iupred.score + ANN + SVM + normalization,
              data = protein)

summary(PSSM.lm)

PSSM.preds <- predict(PSSM.lm, 
                      newdata = (protein %>% filter(is.na(PSSM)))) + rnorm(sum(is.na(protein$PSSM)), 0, sigma(PSSM.lm))

protein <- protein %>%
  mutate(PSSM = replace(PSSM, is.na(PSSM), PSSM.preds))

protein <- protein %>% rowwise %>% mutate(Consensus = replace(Consensus, is.na(Consensus), mean(SVM, PSSM, ANN)))

# Use random forests to fill NA values
# filled_nas <- protein %>% 
#  select(-c(Response, Set)) %>% 
#  as.data.frame() %>% 
#  missForest()

# Merge Response and Set back with other variables
# protein.clean <- filled_nas$ximp %>%  cbind(., protein %>% select(Set, Response)) %>% select(Set, SiteNum, Response, Amino.Acid, Iupred.score, ANN, PSSM, SVM, Consensus, normalization)

# Confirming that all NAs have been dealt with
plot_missing(protein.clean)

# Create a .csv file of the cleaned data
protein.clean %>% write_csv("MM_ProteinCleaned.csv")
