library(parallel)
library(png)
library(caret)
library(IM)
library(MASS)
library(ggplot2)
library(dplyr)
library(UBL)
library(magick)

# GENERIC FUNCTION FOR THIS SCRIPT 

getKrawt <- function(img, n = 90, v = 0.5){
  # Input: img = a PNG image
  #        n = order of the Krawtchouk moment variant.
  # Output: Krawtchouk moment invariant matrix of order n for img
  img <- histeq(img)
  obj <- momentObj(I=img,type="krawt",order=c(n,n), v)
  Moments(obj)<- NULL
  Reconstruct(obj) <- dim(img)
  Invariant(obj) <- c(3, 90)
  return(obj)
} 

getK <- function(data){
  # Input: data = set of binary images in PNG format
  # Output: K = an invariant matrix with each row = vectorized Krawtchouk invariant matrix for each image in data.
  no.cores <- detectCores() - 2 
  clust <- makeCluster(no.cores) 
  clusterEvalQ(clust, {
    library(IM)
  })
  KMList <- parLapply(clust, data, getKrawt)
  KMList.vector.form <- lapply(KMList, function(k) as.vector(k@invariant[1:nrow(k@invariant)-1,1:ncol(k@invariant)-1]))
  K <- do.call(rbind, KMList.vector.form)
  return(K)
}

image_preprocess <- function(path){
  # Input: path = set of path of images 
  # Output: a dataframe with K obtained from getK associated with ID & species category for each row
  image_set <- lapply(path, readPNG) 
  sp <- as.factor(gsub('.[0-9]*.png','',path))
  return(data.frame(sp = sp, getK(image_set)))
}

## get the family of a species 
getFamily <- function(species){
  # Input: species = type of the species 
  # Output: the family of the species. 
  
  if(species %in% c("A.gressitti", "B.karnyi", "L.alba", "S.aquila", "S.princeps"))
  {"Sarcophagidae"} 
  else if (species %in% c("S.nudiseta"))
  {"Muscidae"} 
  else {"Calliphoridae"}
}

resizeImage <- function(path, width = "256"){
  # Input: path = set of path of (binarized) image.
  #        width = desired width after compressing the image. 
  # Output: compressed image with the desired width
  image <- image_read(paste0("./dataset_original/", path))
  image <- image_scale(image, width)
  image_write(image, paste0("./dataset_compressed/", path), format = "png")
}

# RESIZE AND COMPRESS THE IMAGE (only run once)
setwd("./dataset_original")
filename <- list.files()
setwd("..")
sapply(filename, resizeImage)

# LOAD THE COMPRESSED DATASET AND SPLIT INTO TRAINING AND TESTING DATASET 
# USING STRATIFIED SAMPLING (only run once)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd('dataset_compressed')
image_path <- list.files()

image_df <- image_preprocess(image_path)
setwd("..")

## First we save it into .rds file because it took time to preprocess the resized data 
saveRDS(image_df, file = "./image_df.rds")

## You can directly run this when you had the rds file
image_df <- readRDS("./image_df.rds")

## Get the IDs for k-folds
data <- image_df
set.seed(123)
folds <- createFolds(image_df$sp, k = 5, list = TRUE, returnTrain = FALSE)
      
# CHANGE THE FOLD NUMBER TO PERFORM DIFFERENT EXPERIMENT FOR CROSS-VALIDATION 
testingID <- folds[[1]]
trainingID <- setdiff(seq(nrow(image_df)), testingID)

training <- image_df[trainingID,]
testing <- image_df[testingID,]   # change this for female data or left wings data 

# TRAINING DATA PREPROCESSING 
## Scaling 
scaled_training <- scale(training[,-1])

## Q-mode PCA
pca.invar <- prcomp(scaled_training)
pca.invar.var <- pca.invar$sdev^2
pca.invar.var.per <- pca.invar.var / sum(pca.invar.var) * 100
pca.invar.scores <- data.frame(sp = training[,1], pca.invar$x[,1:which(cumsum(pca.invar.var.per) >= 92)[1]]) # Chosing the number of principal components
colnames(pca.invar.scores)[1] <- "sp"

## lda
lda.model <- lda(sp~., data=pca.invar.scores)

LD <- as.matrix(pca.invar.scores[,-1]) %*% lda.model$scaling

training_lda <- data.frame(sp = training[,1], LD)

## smote + enn
training_final <- SmoteClassif(sp~., training_lda)
training_final <- ENNClassif(sp~., training_final)[[1]]

# TESTING DATA PREPROCESSING 
matrix.mul <- scale(testing[,-1], center = attr(scaled_training, "scaled:center"), 
                    scale = attr(scaled_training, "scaled:scale")) %*% pca.invar$rotation[,1:which(cumsum(pca.invar.var.per) >= 92)[1]] %*% lda.model$scaling
testing_final <- data.frame(sp = testing[,1], matrix.mul)
rownames(testing_final) <- seq(nrow(testing_final))

# PREPARE THE DATA FOR GUIDE program
train <- data.frame(training_final, weight =  1)
test <- data.frame(testing_final, weight = 0)
trainANDtest <- rbind(train,test)
write.table(trainANDtest, "./wings.txt", sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE)
