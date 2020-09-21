# libraries
library(bestglm)
library(caret)
library(tidyverse)
library(car)
library(beepr)

# Read in data
filename <- "ProteinCleaned.csv"

protein <- read_csv(filename) %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid),
         Response = as_factor(Response))
     
# separate by train and test    
protein.train <- protein %>% 
  filter(!is.na(Response)) %>% 
  select(-c(Set, SiteNum))
protein.test <- protein %>% filter(is.na(Response))

source("HelpfulFunctions.R")

control <- trainControl(method="repeatedcv",
                       number=10,
                       repeats=3,
                       summaryFunction = f1)
### Machine Learning Algorithms ###
alg = 'bagFDA'
t.grid <- expand.grid(degree=1, nprune=13)
model <- train(form=Response~.,
            data=protein.train,
            method=alg,
            metric="F1",
            preProcess="pca",
            trControl=control
            ,tuneGrid=t.grid
)
beep()
plot(model)
model$bestTune
model$results

# Predict test set and write results to file
preds <- data.frame(Id=protein.test$SiteNum, Predicted=predict(model, newdata=protein.test))
write_csv(preds, paste0("PB_", toupper(alg), "_Preds.csv"))




