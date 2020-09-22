### Cleaning the Protein Data Set

# libraries
library(tidyverse)
library(DataExplorer)
library(zoo)
library(gridExtra)
library(car)
library(missForest)
library(caret)

### Load in data and create one tibble
protein.train <- read_csv("ProteinTrain.csv")
protein.test <- read_csv("ProteinTest.csv")
protein <- bind_rows(train = protein.train, test = protein.test, .id = "Set")
# make Amino Acid a factor and Response a boolean
protein <- protein %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid)) %>% 
  mutate(Response = as.logical(Response))

### Missing values
plot_missing(protein)


# Correlation matrix
plot_correlation(protein.c, type = "continuous", 
                 cor_args = list(use = "pairwise.complete.obs"))

### Machine Learning Imputation for PSSM
control <- trainControl(method="repeatedcv",
                        number=10,
                        repeats=3
)
alg = 'gbm'
t.grid <- expand.grid(shrinkage=.1, interaction.depth=3, n.trees=100, n.minobsinnode=10)
gbm <- train(form=PSSM~.,
               data=protein %>% select(-c(Set, Response, SiteNum, Consensus)) %>% filter(!is.na(PSSM)),
               method=alg,
               trControl=control
               ,tuneGrid=t.grid
)

protein <- protein %>% 
  rowwise %>% 
  mutate(PSSM = replace_na(PSSM, predict(gbm, data.frame(SVM, ANN, Iupred.score, Amino.Acid, normalization))))

# Fill Consensus variable with average of SVM, PSSM, and ANN
protein.c <- protein %>%
  rowwise() %>% 
  mutate(Consensus = replace_na(Consensus, mean(c(SVM, PSSM, ANN))))

### Visualization

source("HelpfulFunctions.R")
# Scatterplots with smoothers (Using a function I made found in HelpfulFunctions.R)
sp1 <- protein %>% scp_smooth(SVM, Response, col = Amino.Acid)
sp2 <- protein %>% scp_smooth(normalization, Response, col = Amino.Acid)
sp3 <- protein %>% scp_smooth(Iupred.score, Response, col = Amino.Acid)

grid.arrange(sp1, sp2, sp3, ncol = 2)

# Check that there are no missing values
plot_missing(protein.c)

# Write cleaned data out to file
protein.c %>% write_csv("ProteinCleaned.csv")




