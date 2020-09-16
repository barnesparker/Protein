### Cleaning the Protein Data Set

# libraries
library(tidyverse)
library(DataExplorer)

protein.train <- read_csv("ProteinTrain.csv", stringsasfactors=T)
protein.train$Amino.Acid <- as.factor(protein.train$Amino.Acid)
unique(protein.train$Amino.Acid)
