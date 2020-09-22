# libraries
library(caret)
library(tidyverse)
library(car)
library(beepr)
library(MLmetrics)

protein.clean <- read_csv("ProteinCleaned.csv")

# Make the Amino.Acid and Response variables factors
protein.clean <- protein.clean %>% mutate_at(vars(Amino.Acid, Response), factor)

# Trying to filter out highly correlated variables
# corrTrans <- preProcess(x = protein.clean %>% select(-c(SiteNum, Response)), 
#                        method = "corr")

# protein.corr <- predict(corrTrans, newdata = protein.clean)

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
  select(-Set, -SiteNum) # %>% 
  # mutate(Response = ifelse(Response == 1, 'Yes', 'No'))

protein_test <- protein.clean %>% filter(Set == 'test') %>% select(-Set) 

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = TRUE)
  c(F1 = f1_val)
}

xgbTree.model <- train(Response ~ .,
                       data = protein_train,
                       method = "xgbTree",
                       metric = 'F1',
                       trControl = trainControl(summaryFunction = f1, 
                                                method = 'cv',
                                                number = 10),
                       tuneGrid = expand.grid(nrounds = c(50, 250),
                                              max_depth = 1,
                                              eta = c(.35, .355, .36, .365, .37, .375),
                                              gamma = 0,
                                              colsample_bytree = .6,
                                              min_child_weight = 1,
                                              subsample = .833333)
                       )

xgbTree.model

xgbTree.preds <- predict(xgbTree.model, protein_test)

xgbTree.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                                 Predicted = xgbTree.preds) # %>% 
# mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(xgbTree.submission, "MM-xgbTree-preds.csv", row.names = FALSE)

gbm.model <- train(Response ~ .,
                   data = protein_train,
                   method = "gbm",
                   metric = 'F1',
                   trControl = trainControl(summaryFunction = f1, 
                                            method = 'cv',
                                            number = 10),
                   tuneGrid = expand.grid(n.trees = 50,
                                          interaction.depth = c(1, 2, 3),
                                          shrinkage = .1,
                                          n.minobsinnode = 10)
)

gbm.model

gbm.preds <- predict(gbm.model, protein_test)

gbm.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                             Predicted = gbm.preds) # %>% 
# mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(gbm.submission, "MM-gbm-preds.csv", row.names = FALSE)

rf.model <- train(Response ~ .,
                  data = protein_train,
                  method = "rf",
                  metric = 'F1',
                  trControl = trainControl(summaryFunction = f1,
                                           method = 'cv',
                                           number = 10),
                  tuneGrid = expand.grid(mtry = c(2, 3, 4, 5))
                  )

rf.model

rf.preds <- predict(rf.model, protein_test)

rf.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                             Predicted = rf.preds) # %>% 
# mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(rf.submission, "MM-rf-preds.csv", row.names = FALSE)

hdda.model <- train(Response ~ .,
                    data = protein_train,
                    method = "hdda",
                    metric = 'F1',
                    trControl = trainControl(summaryFunction = f1,
                                             method = 'cv',
                                             number = 10),
                    tuneGrid = expand.grid(threshold = c(.1, .133333, .15),
                                           model = 'all')
                    )

hdda.model

hdda.preds <- predict(hdda.model, protein_test)

hdda.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                             Predicted = hdda.preds) # %>% 
# mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(hdda.submission, "MM-hdda-preds.csv", row.names = FALSE)

pda2.model <- train(Response ~ .,
                    data = protein_train,
                    method = "pda2",
                    metric = 'F1',
                    trControl = trainControl(summaryFunction = f1,
                                             method = 'cv',
                                             number = 10),
                    tuneGrid = expand.grid(df = c(4, 5, 6, 7))
                    )

pda2.model

pda2.preds <- predict(pda2.model, protein_test)

pda2.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                              Predicted = pda2.preds) # %>% 
# mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(pda2.submission, "MM-pda2-preds.csv", row.names = FALSE)

models.train.df <- data.frame(Response = protein_train %>% select(Response),
                              P1 = predict(xgbTree.model, protein_train),
                              P2 = predict(gbm.model, protein_train),
                              P3 = predict(rf.model, protein_train),
                              P4 = predict(hdda.model, protein_train),
                              P5 = predict(pda2.model, protein_train))

models.test.df <- data.frame(P1 = predict(xgbTree.model, protein_test),
                             P2 = predict(gbm.model, protein_test),
                             P3 = predict(rf.model, protein_test),
                             P4 = predict(hdda.model, protein_test),
                             P5 = predict(pda2.model, protein_test))

models.df <- bind_rows(train = models.train.df, test = models.test.df, .id = "Set")

ensemble.rf <- train(Response ~ .,
                     data = models.train.df,
                     method = 'rf',
                     metric = 'F1',
                     trControl = trainControl(summaryFunction = f1,
                                              method = 'cv',
                                              number = 10),
                     tuneLength = 3)

ensemble.rf

ensemble.preds <- predict(ensemble.rf, models.test.df)

ensemble.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                                  Predicted = ensemble.preds) # %>% 
  # mutate(Predicted = ifelse(Predicted == 'Yes', 1, 0))

write.csv(ensemble.submission, "MM-ensemble-preds.csv", row.names = FALSE)

ensemble.ada <- train(Response ~ .,
                     data = models.train.df,
                     method = 'ada',
                     metric = 'F1',
                     trControl = trainControl(summaryFunction = f1,
                                              method = 'cv',
                                              number = 10),
                     tuneLength = 3)

ensemble.ada

voting.preds <- models.test.df %>% mutate(P1 = ifelse(P1 == TRUE, 1, 0),
                                          P2 = ifelse(P2 == TRUE, 1, 0),
                                          P3 = ifelse(P3 == TRUE, 1, 0),
                                          P4 = ifelse(P4 == TRUE, 1, 0),
                                          P5 = ifelse(P5 == TRUE, 1, 0),
                                          Predicted = ifelse(P1 + P2 + P3 + P4 + P5 == 2, TRUE, FALSE)
                                          )

voting.submission <- data.frame(Id = protein_test %>% pull(SiteNum),
                                Predicted = voting.preds %>% select(Predicted))

write.csv(voting.submission, "MM-voting-preds.csv", row.names = FALSE)
