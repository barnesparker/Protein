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

# Logistic Regression
models <- protein.01 %>% 
  filter(!is.na(Response)) %>% 
  select(-c(Set, SiteNum)) %>%
  as.data.frame() %>% 
  drop_na() %>% 
  bestglm(IC = "AIC", family = binomial, method = "exhaustive")
models$BestModel
beep()


protein.m <- protein.01 %>% 
  filter(!is.na(Response)) %>% 
  select(-c(Set, SiteNum))
log_mod <- glm(Response ~ normalization+Iupred.score+Consensus, data = protein.m, family = binomial)
pred_data <- protein.01 %>% 
  filter(is.na(Response)) %>% 
  select(-c(Set, SiteNum))
pred.probs = predict.glm(log_mod, newdata = pred_data, type="response")
n_cutoff <- 100
cutoff <- seq(0.01,0.99,length=n_cutoff)
misclass <- rep(0,n_cutoff)

for(i in 1:n_cutoff){
  preds <- 1 * (pred.probs > cutoff[i])               ##### predict 
  conf.mat <- table(preds,pred_data)                  ##### get confusion matrix
  misclass[i] <- 1 - sum(diag(conf.mat))/sum(conf.mat)    #### get misclassification rate
}
best_cutoff <- cutoff[which.min(misclass)]
plot(cutoff,misclass,type="l",ylab="Misclassification Rate",xlab="Cutoff")
abline(v = best_cutoff, col = "red")

### Machine Learning

# Oblique Random Forest
.grid <- expand.grid(alpha=seq(0.4, 0.8, length=10),
                         lambda=seq(0, 0.02, length=10))
elnet <- train(form=imdb_score~.,
               data=(imdb.train %>% select(-Set, -movie_title)),
               method="glmnet",
               trControl=trainControl(method="repeatedcv",
                                      number=10, #Number of pieces of your data
                                      repeats=3), #repeats=1 = "cv"
               tuneGrid=enet.grid
)
plot(elnet)
elnet$bestTune