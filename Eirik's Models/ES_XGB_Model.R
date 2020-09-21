# Clear workspace
rm(list=ls())

# libraries
library(caret)
library(tidyverse)
library(car)
library(MLmetrics)

set.seed(777)

protein.clean <- read_csv("../ProteinCleaned.csv")

# Make the Amino.Acid and Response variables factors
protein.clean <- protein.clean %>% mutate_at(vars(Amino.Acid, Response), factor)

trans01 <- preProcess(x=protein.clean %>% select(-c(SiteNum, Response)), method="range",
                      rangeBounds=c(0,1))
protein.01 <- predict(trans01, newdata=protein.clean)

pcTrans <- preProcess(x=protein.clean %>% select(-Response), method="pca")
protein.pca <- predict(pcTrans, newdata=protein.clean)

protein.train <- protein.01 %>% filter(!is.na(Response)) %>%
  mutate(Response = as_factor(Response))

protein.test <- protein.01 %>% filter(is.na(Response))
protein.con <- protein.train %>% select(-c(SVM, ANN, PSSM))


######################
## Machine Learning ##
######################

## Function to compute F1 values
f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = TRUE)
  c(F1=f1_val)
}

## Baseline model using caret's XGBTree
grid_default <- expand.grid(
  nrounds = 100,
  max_depth = 6,
  eta = 0.3,
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

train_control <- caret::trainControl(
  method = "none",
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_base <- train(form=Response~.,
                  data=(protein.train %>% select(-Set)),
                  method="xgbTree",
                  trControl=train_control,
                  tuneGrid=grid_default,
                  verbose=TRUE
)

## Next, start tuning hyperparameters
nrounds <- 1000
tune_grid <- expand.grid(
  nrounds = seq(from = 200, to = nrounds, by = 50),
  eta = c(0.025, 0.05, 0.1, 0.3),
  max_depth = c(2, 3, 4, 5, 6),
  gamma = 0,
  colsample_bytree = 1,
  min_child_weight = 1,
  subsample = 1
)

tune_control <- caret::trainControl(
  method = "cv", # cross-validation
  number = 3, # with n folds 
  #index = createFolds(tr_treated$Id_clean), # fix the folds
  verboseIter = FALSE, # no training log
  allowParallel = TRUE # FALSE for reproducible results 
)

xgb_tune <-train(form=Response~.,
                 data=(protein.train %>% select(-Set)),
                 method="xgbTree",
                 metric=,
                 trControl=tune_control,
                 tuneGrid=tune_grid,
                 verbose=TRUE
)

# helper function for the plots
tuneplot <- function(x, probs = .90) {
  ggplot(x) +
    coord_cartesian(ylim = c(min(x$results$Accuracy), max(x$results$Accuracy))) +
    theme_bw()
}

tuneplot(xgb_tune)
xgb_tune$bestTune

## Next round of tuning
tune_grid2 <- expand.grid(nrounds = seq(from = 50, to = nrounds, by = 50),
                          eta = xgb_tune$bestTune$eta,
                          max_depth = ifelse(xgb_tune$bestTune$max_depth == 2,
                                             c(xgb_tune$bestTune$max_depth:4),
                                             xgb_tune$bestTune$max_depth - 1:xgb_tune$bestTune$max_depth + 1),
                          gamma = 0,
                          colsample_bytree = 1,
                          min_child_weight = c(1, 2, 3),
                          subsample = 1
)

xgb_tune2 <- caret::train(
  form=Response~.,
  data=(protein.train %>% select(-Set)),
  method="xgbTree",
  trControl=tune_control,
  tuneGrid=tune_grid2,
  verbose=TRUE
)

tuneplot(xgb_tune2)
xgb_tune2$bestTune
max(xgb_tune$results$Accuracy)
max(xgb_tune2$results$Accuracy)

## Next tuning round
tune_grid3 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 50),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = 0,
  colsample_bytree = c(0.4, 0.6, 0.8, 1.0),
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = c(0.5, 0.75, 1.0)
)

xgb_tune3 <- caret::train(
  form=Response~.,
  data=(protein.train %>% select(-Set)),
  method="xgbTree",
  trControl=tune_control,
  tuneGrid=tune_grid3,
  verbose=TRUE
)

tuneplot(xgb_tune3)
xgb_tune3$bestTune
max(xgb_tune$results$Accuracy)
max(xgb_tune2$results$Accuracy)
max(xgb_tune3$results$Accuracy)

## Tuning the Gamma
tune_grid4 <- expand.grid(
  nrounds = seq(from = 50, to = nrounds, by = 50),
  eta = xgb_tune$bestTune$eta,
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = c(0, 0.05, 0.1, 0.5, 0.7, 0.9, 1.0),
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune4 <- caret::train(
  form=Response~.,
  data=(protein.train %>% select(-Set)),
  method="xgbTree",
  trControl=tune_control,
  tuneGrid=tune_grid4,
  verbose=TRUE
)

tuneplot(xgb_tune4)
xgb_tune4$bestTune
max(xgb_tune$results$Accuracy)
max(xgb_tune2$results$Accuracy)
max(xgb_tune3$results$Accuracy)
max(xgb_tune4$results$Accuracy)

## Reduce learning rate
tune_grid5 <- expand.grid(
  nrounds = seq(from = 100, to = 5000, by = 50),
  eta = c(0.01, 0.015, 0.025, 0.05, 0.1),
  max_depth = xgb_tune2$bestTune$max_depth,
  gamma = xgb_tune4$bestTune$gamma,
  colsample_bytree = xgb_tune3$bestTune$colsample_bytree,
  min_child_weight = xgb_tune2$bestTune$min_child_weight,
  subsample = xgb_tune3$bestTune$subsample
)

xgb_tune5 <- caret::train(
  form=Response~.,
  data=(protein.train %>% select(-Set)),
  method="xgbTree",
  trControl=tune_control,
  tuneGrid=tune_grid5,
  verbose=TRUE
)

tuneplot(xgb_tune5)
xgb_tune5$bestTune
max(xgb_tune$results$Accuracy)
max(xgb_tune2$results$Accuracy)
max(xgb_tune3$results$Accuracy)
max(xgb_tune4$results$Accuracy)
max(xgb_tune5$results$Accuracy)

## Fit the model and predict
final_grid <- expand.grid(
  nrounds = xgb_tune5$bestTune$nrounds,
  eta = xgb_tune5$bestTune$eta,
  max_depth = xgb_tune5$bestTune$max_depth,
  gamma = xgb_tune5$bestTune$gamma,
  colsample_bytree = xgb_tune5$bestTune$colsample_bytree,
  min_child_weight = xgb_tune5$bestTune$min_child_weight,
  subsample = xgb_tune5$bestTune$subsample
)

xgb_model <- caret::train(
  form=Response~.,
  data=(protein.train %>% select(-Set)),
  method="xgbTree",
  trControl=tune_control,
  tuneGrid=final_grid,
  verbose=TRUE
)

xgb.preds <- data.frame(Id=protein.test$SiteNum, Predicted=predict(xgb_model, newdata=protein.test))
head(xgb.preds, 25)
# Test the F1 score
F1_Score(y_pred = predict(xgb_model, newdata=protein.train),
                             y_true = protein.train$Response, positive = TRUE)
write_csv(x=xgb.preds, path="./ES_XGBPredictions.csv")
