### Helpful R functions ###

# scatterplot with smoother for binary data
scp_smooth <- function(data, x, y, ...) {
  data %>% ggplot(mapping=aes(x, y, col =...)) +
    geom_point() +
    geom_smooth()
}

mach_learn <- function(g, ) {
  g.grid <- expand.grid(n.trees = 100, interaction.depth=1, shrinkage=.1, n.minobsinnode = 10)
  gbm <- train(form=Response~.,
               data=(protein.train %>% select(-Set, -SiteNum)),
               method='gbm',
               trControl=trainControl(method="repeatedcv",
                                      number=10, #Number of pieces of your data
                                      repeats=3,
                                      summaryFunction = prSummary), #repeats=1 = "cv"
               tuneGrid=g.grid
  )
  beep()
}

