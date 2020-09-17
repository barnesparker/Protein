### Cleaning the Protein Data Set

# libraries
library(tidyverse)
library(DataExplorer)
library(zoo)
library(gridExtra)
library(car)
### Load in data and create one tibble
protein.train <- read_csv("ProteinTrain.csv")
protein.test <- read_csv("ProteinTest.csv")
protein <- bind_rows(train=protein.train, test=protein.test, .id="Set")
# make Amino Acid a factor
protein <- protein %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid))

# Summary
summary(protein)


### Missing values
plot_missing(protein.c)

# Correlation matrix
plot_correlation(protein.c, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))

# Replace missing "Consensus" values with the mean of it's SVM and ANN values
protein.c <- protein %>%
  rowwise() %>% 
  mutate(Consensus = mean(c(ANN, SVM)))

# Make linear model to impute missing PSSM values; regress on SVM & Consensus variables

# model
imp_lm <- lm(PSSM ~ SVM+Consensus, data = protein.c)
# Get the rows with missing PSSM values
missing_PSSM <- protein.c %>%
  filter(is.na(PSSM))
# predict those values with model
new_PSSMs <- predict(imp_lm, missing_PSSM)

# Impute predicted values into our data
protein.c <- protein.c %>%
  sapply()

protein.c <- protein.c %>%
  mutate(PSSM =)


source("HelpfulFunctions.R")
# Jittered scatterplot of Consensus and Response, with color as Amino Acid
# (Using a function I made found in HelpfulFunctions.R)
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




