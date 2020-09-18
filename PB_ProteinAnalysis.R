# libraries
library(bestglm)
library(caret)
library(tidyverse)
library(car)
library(beepr)
protein.c <- read_csv("ProteinCleaned.csv")
protein.c <- protein.c %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid))

trans01 <- preProcess(x=protein.c %>% select(-c(SiteNum, Response)), method="range",
                      rangeBounds=c(0,1))
protein.01 <- predict(trans01, newdata=protein.c)

pcTrans <- preProcess(x=protein.c %>% select(-Response), method="pca")
protein.pca <- predict(pcTrans, newdata=protein.c)

### Machine Learning

protein.train <- protein.01 %>% filter(!is.na(Response)) %>% 
  mutate(Response = as_factor(Response))

protein.test <- protein.01 %>% filter(is.na(Response))


# Random Forest

rf <- train(form=Response~.,
            data=(protein.train %>% select(-Set, -SiteNum)),
            method="rf",
            trControl=trainControl(method="repeatedcv",
                                   number=10, #Number of pieces of your data
                                   repeats=3) #repeats=1 = "cv"
            #tuneGrid=enet.grid
)
beep()
plot(rf)
rf$bestTune

### Support Vector Machine

svm <- train(form=Response~.,
             data=(protein.train %>% select(-Set, -SiteNum)),
             method='svmLinear',
             trControl=trainControl(method="repeatedcv",
                                    number=10, #Number of pieces of your data
                                    repeats=3), #repeats=1 = "cv"
             tuneGrid=expand.grid(C = 1:3)
)
beep()
plot(svm)
svm$bestTune
svm$results
confusionMatrix(data = test_set$pred, reference = test_set$obs, mode = "prec_recall")
svm.preds <- data.frame(Id=protein.test$SiteNum, Predicted=predict(svm, newdata=protein.test))
write_csv(x=svm.preds, path="./PB_SVM_Preds.csv")

### Naive Bayes

nb <- train(form=Response~.,
            data=(protein.train %>% select(-Set, -SiteNum)),
            method='nb',
            trControl=trainControl(method="repeatedcv",
                                   number=10, #Number of pieces of your data
                                   repeats=3,
                                   summaryFunction = f1) #repeats=1 = "cv"
            #tuneGrid=enet.grid
)
beep()
plot(nb)
nb$results
