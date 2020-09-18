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

# Use random forests to fill NA values
filled_nas <- protein %>% 
  select(-c(Response, Set)) %>% 
  as.data.frame() %>% 
  missForest()

<<<<<<< Updated upstream
# protein.c %>% 
#  bind_rows(filled_nas %>% )
# forest$ximp
=======
# merge Response and Set back with other variables
protein.c <- filled_nas$ximp %>%  
  merge(protein)
>>>>>>> Stashed changes

# Correlation matrix
plot_correlation(protein.c, type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

# Replace missing "Consensus" values with the mean of it's SVM and ANN values
#protein.c <- protein %>%
 # rowwise() %>% 
#  mutate(Consensus = replace_na(Consensus, mean(c(ANN, SVM)))) %>% 
#  arrange(SiteNum)

# Make linear model to impute missing PSSM values; regress on SVM & Consensus variables

# model
<<<<<<< Updated upstream
imp_lm <- lm(PSSM ~ SVM + Consensus, data = protein.c)

# Impute predicted values into our data
protein.c <- protein.c %>%
  rowwise() %>% 
  mutate(PSSM = replace_na(PSSM, predict(imp_lm, newdata = tibble(SVM, Consensus))))

protein.test <- aregImpute(~ Consensus + SVM + normalization, data = protein.c)

protein.test$imputed

=======
#imp_lm <- lm(PSSM ~ SVM+Consensus, data = protein.c)

# Impute predicted values into our data
#protein.c <- protein.c %>%
#  rowwise() %>% 
#  mutate(PSSM = replace_na(PSSM, predict(imp_lm, newdata = tibble(SVM, Consensus))))
>>>>>>> Stashed changes
source("HelpfulFunctions.R")
# Scatterplots with smoothers (Using a function I made found in HelpfulFunctions.R)
sp1 <- protein %>% scp_smooth(SVM, Response, col = Amino.Acid)
sp2 <- protein %>% scp_smooth(normalization, Response, col = Amino.Acid)
sp3 <- protein %>% scp_smooth(Iupred.score, Response, col = Amino.Acid)

grid.arrange(sp1, sp2, sp3, ncol = 2)

# table of Amino Acid and Response
addmargins(table(protein$Amino.Acid, protein$Response))

### Visualization

plot_missing(protein.c)

# Histograms
protein.c %>% ggplot(mapping = aes(SVM)) + 
  geom_histogram(bins = 100)

protein.c %>% write_csv("ProteinCleaned.csv")




