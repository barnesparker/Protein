### Helpful R functions ###

# scatterplot with smoother for binary data
scp_smooth <- function(data, x, y, ...) {
  data %>% ggplot(mapping=aes(x, y, col =...)) +
    geom_point() +
    geom_smooth()
}