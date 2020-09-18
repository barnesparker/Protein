# libraries
library(bestglm)
library(caret)
library(tidyverse)
library(car)
library(beepr)
library(MLmetrics)

protein.clean <- read_csv("ProteinCleaned.csv")

# Make the Amino.Acid and Response variables factors
protein.clean <- protein.clean %>% mutate_at(vars(Amino.Acid, Response), factor)

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

protein_train <- protein.clean %>% filter(Set == 'train') %>% 
  select(-Set, -SiteNum) %>% 
  mutate(Response = ifelse(Response == 1, 'Yes', 'No'))

protein_test <- protein.clean %>% filter(Set == 'test') %>% select(-Set) 

gbm.model <- train(Response ~ .,
                   data = protein_train,
                   method = "gbm",
                   metric = "F",
                   trControl = trainControl(summaryFunction = prSummary, 
                                            classProbs = TRUE,
                                            method = 'cv',
                                            number = 10),
                   preProcess = 'corr',
                   tuneGrid = expand.grid(n.trees = 100,
                                          interaction.depth = 1,
                                          shrinkage = .1,
                                          n.minobsinnode = 10
                                          )
                   )
gbm.model

gbm.preds <- predict(gbm.model, protein_test)

gbm.submission <- data.frame(Id = protein_test %>% select(SiteNum), 
                             Predicted = gbm.preds) %>% 
  mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(gbm.submission, "MM-gbm-preds.csv", row.names = FALSE)
