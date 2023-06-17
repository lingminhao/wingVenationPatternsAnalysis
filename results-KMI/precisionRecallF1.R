
# This is a R script to produce the precision and recall for each species in 
# 5-fold cross validation experiment

library(dplyr)
library(caret)

sp <- c("A.gressitti", "B.karnyi", "C.albiceps", "C.albiceps_mutant", "C.bezziana", "C.megacephala", "C.nigripes", 
        "C.rufifacies", "C.vicina", "L.alba", 'L.sericata', "S.aquila", "S.nudiseta", "S.princeps")

# For 80/20 train test split
gf_1.pro <- read.table(file = "./fold1/gf_1.pro", header = TRUE, sep = "") %>% filter(train == "n")
gf_2.pro <- read.table(file = "./fold2/gf_2.pro", header = TRUE, sep = "") %>% filter(train == "n")
gf_3.pro <- read.table(file = "./fold3/gf_3.pro", header = TRUE, sep = "") %>% filter(train == "n")
gf_4.pro <- read.table(file = "./fold4/gf_4.pro", header = TRUE, sep = "") %>% filter(train == "n")
gf_5.pro <- read.table(file = "./fold5/gf_5.pro", header = TRUE, sep = "") %>% filter(train == "n")
gf.pro <- rbind(gf_1.pro, gf_2.pro, gf_3.pro, gf_4.pro, gf_5.pro)

summaryMetric <- function(folds){
  cm <- table(factor(folds$predicted, levels = sp), factor(folds$observed, levels = sp))
  n = sum(cm) # number of instances
  nc = nrow(cm) # number of classes
  diag = diag(cm) # number of correctly classified instances per class 
  rowsums = apply(cm, 1, sum) # number of instances per class
  colsums = apply(cm, 2, sum) # number of predictions per class
  
  recall = diag / colsums 
  precision = diag / rowsums 
  f1 = 2 * precision * recall / (precision + recall) 
  
  return(data.frame(precision, recall, f1))
}
