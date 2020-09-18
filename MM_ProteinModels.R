# libraries
library(bestglm)
library(caret)
library(tidyverse)
library(car)
library(beepr)

protein.c <- read_csv("ProteinCleaned.csv")
protein.c <- protein.c %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid)) %>% 
  mutate(Response = as_factor(Response))

# trans01 <- preProcess(x = protein.c %>% select(-c(SiteNum, Response)), 
#                      method = "range", 
#                      rangeBounds = c(0,1))
# protein.01 <- predict(trans01, newdata = protein.c)

# PCA Transformation of the data
# pcTrans <- preProcess(x = protein.c %>% select(-Response), 
#                       method = "pca")
# protein.pca <- predict(pcTrans, newdata = protein.c)

######################
## Machine Learning ##
######################

protein_train <- protein.c %>% filter(!is.na(Response))
protein_test <- protein.c %>% filter(is.na(Response))

rpart <- train(Response ~ ., 
               data = protein.train %>% select(-Set, -SiteNum),
               method = "rpart",
               metric = "F",
               trControl = trainControl(summaryFunction = prSummary, 
                                        classProbs = TRUE,
                                        method = 'cv',
                                        number = 10),
               tuneLength = 5)