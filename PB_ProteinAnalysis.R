# libraries
library(caret)
library(tidyverse)
library(beepr)
library(MLmetrics)
# Read in data
filename <- "ProteinCleaned.csv"

protein.01 <- filename %>% 
  read_csv() %>% 
  mutate_at(vars(Amino.Acid, Response), factor) %>% 
  select(-c(SiteNum, Response)) %>% 
  preProcess(method="range", rangeBounds=c(0,1)) %>% 
  predict(protein.clean)

# separate by train and test    
protein.train <- protein.01 %>% 
  filter(!is.na(Response)) %>% 
  select(-c(Set, SiteNum))
protein.test <- protein.01 %>% filter(is.na(Response))

source("HelpfulFunctions.R")

control <- trainControl(method="repeatedcv",
                       number=10,
                       repeats=3,
                       summaryFunction = f1)

### Machine Learning Algorithms ###
alg = 'svmLinearWeights2'
t.grid <- expand.grid(degree=1,
                      nprune=12:15)

model <- train(form=Response~.,
            data=protein.train,
            method=alg,
            metric="F1",
            trControl=control
            #,tuneGrid=t.grid
)
beep()
plot(model)
model$bestTune
model$results
# Predict test set and write results to file
preds <- data.frame(Id=protein.test$SiteNum, Predicted=predict(model, newdata=protein.test))
write_csv(preds, paste0("PB_", toupper(alg), "_Preds.csv"))




