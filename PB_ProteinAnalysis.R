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

trans01 <- preProcess(x=protein.c %>% select(-c(SiteNum, Response)), method="range",
                      rangeBounds=c(0,1))
protein.01 <- predict(trans01, newdata=protein.c)

pcTrans <- preProcess(x=protein.c %>% select(-Response), method="pca")
protein.pca <- predict(pcTrans, newdata=protein.c)

### Machine Learning

protein.train <- protein.01 %>% filter(!is.na(Response)) %>%
  mutate(Response = as_factor(Response))

protein.test <- protein.01 %>% filter(is.na(Response))
protein.con <- protein.train %>% select(-c(SVM, ANN, PSSM))
# Random Forest

rf <- train(form=Response~.,
            data=(protein.train %>% select(-Set, -SiteNum)),
            method="rf",
            metric="F",
            trControl=trainControl(method="repeatedcv",
                                   number=10, #Number of pieces of your data
                                   repeats=1,
                                   summaryFunction = prSummary), #repeats=1 = "cv"
            tuneGrid=expand.grid(mtry = 2)
)
beep()
plot(rf)
rf$bestTune

### Support Vector Machine

svm <- train(form=Response~.,
             data=(protein.train %>% select(-Set, -SiteNum)),
             method='svmSpectrumString',
             metric = "F",
             trControl=trainControl(method="repeatedcv",
                                    number=10, #Number of pieces of your data
                                    repeats=3,
                                    #classProbs = TRUE,
                                    summaryFunction = prSummary)
             #tuneGrid=expand.grid(degree = 3, scale = c(.25, .5, .75), C = .25)
)
beep()
plot(svm)
svm$bestTune
svm$results
svm.preds <- data.frame(Id=protein.test$SiteNum, Predicted=predict(svm, newdata=protein.test))
write_csv(x=svm.preds, path="./PB_SVM_POLY_Preds.csv")

### Naive Bayes

nb <- train(form=Response~.,
            data=(protein.train %>% select(-Set, -SiteNum)),
            method='nb',
            metric = 'F',
            trControl=trainControl(method="repeatedcv",
                                   summaryFunction = prSummary,
                                   number=10, #Number of pieces of your data
                                   repeats=3) #repeats=1 = "cv"
            #tuneGrid=enet.grid
)
beep()
plot(nb)
nb$results


### Gradient Boosting
g.grid <- expand.grid(n.trees = 60, interaction.depth=4, shrinkage=c(.8,.1, .1), n.minobsinnode = c(9, 10, 11))
gbm <- train(form=Response~.,
             data=(protein.train %>% select(-Set, -SiteNum)),
             method='gbm',
             metric="F",
             trControl=trainControl(method="repeatedcv",
                                    #classProbs = TRUE,
                                    number=10, #Number of pieces of your data
                                    repeats=3,
                                    summaryFunction = prSummary), #repeats=1 = "cv"
             tuneGrid=g.grid
)
beep()
plot(gbm)
gbm$results
gbm$bestTune

g.grid <- expand.grid(n.trees = 60, interaction.depth=4, shrinkage=c(.8,.1, .1), n.minobsinnode = c(9, 10, 11))
qda <- train(form=Response~.,
             data=(protein.train %>% select(-Set, -SiteNum)),
             method='ordinalNet',
             metric="F",
             trControl=trainControl(method="repeatedcv",
                                    #classProbs = TRUE,
                                    number=10, #Number of pieces of your data
                                    repeats=3,
                                    summaryFunction = prSummary) #repeats=1 = "cv"
             #tuneGrid=g.grid
)
beep()
plot(qda)
qda$results
qda$bestTune
