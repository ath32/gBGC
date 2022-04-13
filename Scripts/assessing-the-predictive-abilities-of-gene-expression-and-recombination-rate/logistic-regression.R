#Import packages
library(dplyr)
library(caret)
library(regclass)

#Load data
source <- "./Project_Iso/Mouse/mouse_px_gc3.csv"
data <- read.csv(source, header = TRUE)

#Create dummies to allow logistic regression
dummies <- data
dummies$PxAbundance <- log(dummies$PxAbundance)
dummies <- dummies %>% filter(PxAbundance != '-Inf')

#Create and run model
model.1 <- glm(TAA~PxAbundance+GC, data = dummies, family = binomial)
print(summary(model.1))
cor.test(data$PxAbundance, data$GC, method='spearman')

#Calculate variance inflation factors
VIF(model.1)
