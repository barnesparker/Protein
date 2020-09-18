### Helpful R functions ###

# scatterplot with smoother for binary data
scp_smooth <- function(data, x, y, ...) {
  data %>% ggplot(mapping=aes(x, y, col =...)) +
    geom_point() +
    geom_smooth()
}

# Linear Imputation -- DOESN'T WORK YET
lin_imp <- function(df, na_var, pred_vars) {
  # Make another tibble for model training
  pred_df <- df %>% select({{na_var}}, {{pred_vars}})
  
  # Train a linear model
  
  mod <- lm(pred_df)
  
  # impute na values with model
  df %>%  
    rowwise() %>%
    mutate("Hello" = replace_na({{na_var}}, predict(mod, newdata = tibble({{pred_vars}}))))
}

