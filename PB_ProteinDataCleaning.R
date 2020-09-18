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
  mutate(Amino.Acid = as_factor(Amino.Acid))

# Summary
summary(protein)


### Missing values
plot_missing(protein)

# Correlation matrix
plot_correlation(protein.c, type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

# Use random forests to fill NA values
filled_nas <- protein %>% 
  select(-c(Response, Set)) %>% 
  as.data.frame() %>% 
  missForest()

# merge Response and Set back with other variables
protein.c <- filled_nas$ximp %>%  
  merge(protein)

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


protein.c %>% write_csv("ProteinCleaned.csv")




