### Helpful R functions ###

# scatterplot with smoother for binary data
scp_smooth <- function(data, x, y, ...) {
  data %>% ggplot(mapping=aes(x, y, col =...)) +
    geom_point() +
    geom_smooth()
}

f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = TRUE)
  c(F1 = f1_val)
}

