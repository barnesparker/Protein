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

### Logistic Regression

# Find the best model using an exhaustive search and Cross-validation 
# (Make sure Response is last variable)
# THIS TAKES A MINUTE OR SO
models <- protein.01 %>% 
  filter(!is.na(Response)) %>% 
  select(-c(Set, SiteNum)) %>%
  as.data.frame() %>% 
  bestglm(IC = "AIC", family = binomial, method = "exhaustive")
models$BestModel
beep()
# looks like only Consensus and normalization made the best model

# make a new data set of just significant variables found from bestglm
protein.m <- protein.01 %>% 
  select(Consensus, normalization, Response)

# Train the model
log_mod <- glm(Response ~., data = protein.m %>% filter(!is.na(Response)), family = binomial)

# Separate the data for prediction
pred_data <- protein.m %>% 
  filter(is.na(Response))

# make a vector of prediction probabilities
pred.probs = predict.glm(log_mod, newdata = pred_data, type="response")

### Now we will determine the best cutoff 

# how many cutoffs will there be?
n_cutoff <- 999
# make the 
cutoffs <- seq(0.001,0.99,length=n_cutoff)
# initizalize the vector for misclass rate
misclass <- rep(0,n_cutoff)

# find misclass rate for each cutoff
for(i in seq_along(cutoffs)){
  preds <- 1 * (pred.probs > cutoffs[i])               ##### predict 
  conf.mat <- table(preds, pred.probs)                  ##### get confusion matrix
  misclass[i] <- 1 - sum(diag(conf.mat))/sum(conf.mat)    #### get misclassification rate
}
# find best cutoff
best_cutoff <- cutoffs[which.min(misclass)]
# plot the misclass rates
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