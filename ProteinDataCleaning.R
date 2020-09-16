### Cleaning the Protein Data Set

# libraries
library(tidyverse)
library(DataExplorer)
library(zoo)
library(gridExtra)
library(car)
### Load in data and create one tibble
protein.train <- read_csv("ProteinTrain.csv")
protein.test <- read_csv("ProteinTest.csv")
protein <- bind_rows(train=protein.train, test=protein.test, .id="Set")
# make Amino Acid a factor
protein <- protein %>% 
  mutate(Amino.Acid = as_factor(Amino.Acid))

# Summary
summary(protein)


# Missing values
plot_missing(protein)

# Replace missing "Consensus" values with the mean of it's SVM and ANN values
protein.c <- protein %>%
  mutate(PSSM = na.approx(PSSM, ANN)) %>% 
  rowwise() %>% 
  mutate(Consensus = mean(Consensus, mean(c(ANN, SVM))))

protein.d <- protein %>%
  mutate(PSSM = na.approx(PSSM, ANN), Consensus = na.approx(SVM, ANN)) 

# Correlation matrix
plot_correlation(protein.c, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))

# Jittered scatterplot of Consensus and Response, with color as Amino Acid
sp1 <- protein %>% ggplot(mapping=aes(x=SVM, y=Response, col = Amino.Acid)) +
  geom_point() +
  geom_smooth()

sp2 <- protein %>% ggplot(mapping=aes(x=normalization, y=Response, col = Amino.Acid)) +
  geom_point() +
  geom_smooth()

sp3 <- protein %>% ggplot(mapping=aes(x=Iupred.score, y=Response, col = Amino.Acid)) +
  geom_point() +
  geom_smooth()

grid.arrange(sp1, sp2, sp3, ncol = 2)
# table of Amino Acid and Response
addmargins(table(protein$Amino.Acid, protein$Response))

# Just drop ANN, PSSM, and Consensus because they are so heavily correlated to
# SVM and SVM has 0 missing 
#protein.c <- protein %>% select(-c(ANN, PSSM, Consensus))

### Visualization

plot_missing(protein.c)
# Histograms
protein.c %>% ggplot(mapping = aes(SVM)) + 
  geom_histogram(bins = 100)


protein.c %>% write_csv("ProteinCleaned.csv")




