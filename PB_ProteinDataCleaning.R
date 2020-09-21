### Cleaning the Protein Data Set

# libraries
library(tidyverse)
library(DataExplorer)
library(zoo)
library(gridExtra)
library(car)
library(missForest)

### Load in data and create one tibble
protein.train <- read_csv("ProteinTrain.csv")
protein.test <- read_csv("ProteinTest.csv")
protein <- bind_rows(train = protein.train, test = protein.test, .id = "Set")
# make Amino Acid a factor
protein <- protein %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid)) %>% 
  arrange(SiteNum)

# Summary
summary(protein)


### Missing values
plot_missing(protein)

# Correlation matrix
plot_correlation(protein.c, type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

# Use random forest to fill PSSM NA values
filled_nas <- protein %>% 
  select(-c(Response, Set, SiteNum, Consensus)) %>% 
  as.data.frame() %>% 
  missForest()

filled_nas <- filled_nas$ximp

protein <- protein %>%  
  mutate(PSSM = filled_nas$PSSM,
         Response = as.logical(Response))

# Fill Consensus variable with average of SVM, PSSM, and ANN
protein.c <- protein.c %>%
  rowwise() %>% 
  mutate(Consensus = replace_na(Consensus, mean(c(SVM, PSSM, ANN))))

# table of Amino Acid and Response
addmargins(table(protein$Amino.Acid, protein$Response))

### Visualization

source("HelpfulFunctions.R")
# Scatterplots with smoothers (Using a function I made found in HelpfulFunctions.R)
sp1 <- protein %>% scp_smooth(SVM, Response, col = Amino.Acid)
sp2 <- protein %>% scp_smooth(normalization, Response, col = Amino.Acid)
sp3 <- protein %>% scp_smooth(Iupred.score, Response, col = Amino.Acid)

grid.arrange(sp1, sp2, sp3, ncol = 2)

# Histograms
protein.c %>% ggplot(mapping = aes(SVM)) + 
  geom_histogram(bins = 100)


# Check that there are no missing values
plot_missing(protein.c)

# Write cleaned data out to file
protein.c %>% write_csv("ProteinCleaned.csv")




