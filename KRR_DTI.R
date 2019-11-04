########REQUIRED PACKAGES########
#################################
install.packages("Hmisc")
library(Hmisc)
install.packages("gplots")
library(gplots)
install.packages("ChemmineR")
source("https://bioconductor.org/biocLite.R")
biocLite("ChemmineR")
library("ChemmineR")
biocLite("ChemmineOB")
library("ChemmineOB")
source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")
library("Biostrings")
library("seqinr")
source("https://bioconductor.org/biocLite.R")
biocLite("GOSemSim")
library('GOSemSim')
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library('org.Hs.eg.db')
install.packages("png")
require("png")
source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")
require("clusterProfiler")
install.packages("kernlab")
library("kernlab")
install.packages("CVST")
library("CVST")
library("reshape")
library(tidyr)
source("https://bioconductor.org/biocLite.R")
biocLite("Rchemcpp")
install.packages("data.table")
library("data.table")
install.packages("gdata")
library("gdata")
install.packages("cvTools")
library("cvTools")



########DATASET LOADING AND PREPARATION########
###############################################
data <- read.table("DTI_100molecules.txt")
drug_id <- read.table("Drug_ids.txt",stringsAsFactors = FALSE)
row.names(data) <- drug_id$V1
protein_id <- read.table("Protein_ids.txt")
names(data) <- protein_id$V1



########BASIC STATISTICS AND PLOTS########
##########################################
summary(data)
describe(data)
data_matrix <- as.matrix(data)
heatmap(data_matrix)
hist(data_matrix)
heatmap.2(data_matrix,col=redgreen(75), scale="none", margins=c(5,10), density.info="none", trace="none")
heatmap.2(data_matrix,col=redgreen(1500), scale="none", margins=c(4,9), density.info="none", trace="none")
heatmap.2(data_matrix,col=redgreen(1500), scale="none", margins=c(2,2), density.info="none", trace="none",Rowv = FALSE,Colv = FALSE,dendrogram = "none", xlab="PROTEINS",ylab="DRUGS",labRow = FALSE,labCol = FALSE)



#######DRUG SMILES LOADING AND CONVERSION TO SDF########
########################################################
smiset <- read.SMIset("Drug_smiles.txt")
view(smiset[1:4]) 
cid(smiset[1:4]) 
drug_smiles_sdf <- smiles2sdf(smiset)
drug_smiles_sdf[[1]]



########GENERATING MOLECULAR FINGERPRINTS USING OpenBabel########
#################################################################
FP2 <- fingerprintOB(drug_smiles_sdf, "FP2")
FP3 <- fingerprintOB(drug_smiles_sdf, "FP3")
FP4 <- fingerprintOB(drug_smiles_sdf, "FP4")

#MACCS <- fingerprintOB(drug_smiles_sdf, "MACCS")
####Based on the used OS, R package, and ChemmineOB package version, MACCS may be not available



########BUILDING THREE FINGERPRINT-BASED TANIMOTO KERNELS########
#################################################################
FP2_Tanimoto_fpSim_Search <- data.frame(FirstDrug=as.numeric(),
                                        SecondDrug=as.numeric(), 
                                        Value=as.numeric(),stringsAsFactors = FALSE) 


for (i in 1:100) {
  for (j in 1:100)
  {
    temp <- fpSim(FP2[i], FP2[j], method="Tanimoto")
    FP2_Tanimoto_fpSim_Search <- rbind(FP2_Tanimoto_fpSim_Search, c(i,j,temp))
  }
}

names(FP2_Tanimoto_fpSim_Search)[1] <- "First.Drug"
names(FP2_Tanimoto_fpSim_Search)[2] <- "Second.Drug"
names(FP2_Tanimoto_fpSim_Search)[3] <- "TanimotoKernel"

FP2_Tanimoto_fpSim_Search$Second.Drug <- sub("^", "Second.Drug.",FP2_Tanimoto_fpSim_Search$Second.Drug)



FP3_Tanimoto_fpSim_Search <- data.frame(FirstDrug=as.numeric(),
                                        SecondDrug=as.numeric(), 
                                        Value=as.numeric(), 
                                        stringsAsFactors=FALSE) 


for (i in 1:100) {
  for (j in 1:100)
  {
      temp <- fpSim(FP3[i], FP3[j], method="Tanimoto")
      FP3_Tanimoto_fpSim_Search <- rbind(FP3_Tanimoto_fpSim_Search, c(i,j,temp))
  }
}

names(FP3_Tanimoto_fpSim_Search)[1] <- "First.Drug"
names(FP3_Tanimoto_fpSim_Search)[2] <- "Second.Drug"
names(FP3_Tanimoto_fpSim_Search)[3] <- "TanimotoKernel"



FP4_Tanimoto_fpSim_Search <- data.frame(FirstDrug=as.numeric(),
                                        SecondDrug=as.numeric(), 
                                        Value=as.numeric(), 
                                        stringsAsFactors=FALSE) 


for (i in 1:100) {
  for (j in 1:100)
  {
      temp <- fpSim(FP4[i], FP4[j], method="Tanimoto")
      FP4_Tanimoto_fpSim_Search <- rbind(FP4_Tanimoto_fpSim_Search, c(i,j,temp))
  }
}

names(FP4_Tanimoto_fpSim_Search)[1] <- "First.Drug"
names(FP4_Tanimoto_fpSim_Search)[2] <- "Second.Drug"
names(FP4_Tanimoto_fpSim_Search)[3] <- "TanimotoKernel"



########SUPPLEMENTARY ANALYSIS##############
########MOLECULAR STRUCTURE EXPLORATION#####
############################################
apset <- sdf2ap(drug_smiles_sdf)
data(apset)
cmp.search(apset, apset[1], type=3, cutoff = 0.3)
draw_sdf(drug_smiles_sdf[[1]])
plotStruc(drug_smiles_sdf[[1]])
groups(drug_smiles_sdf)
atomcount(drug_smiles_sdf[1:100])
atomcountMA(drug_smiles_sdf[1:100])



########MULTIPLE PAIRWISE ALIGNMENT - SMITH-WATERMAN#########
#############################################################
Protein_Seq <- read.fasta("Protein_sequneces.fa", seqtype="AA", as.string="TRUE")
Protein_Sequences_Only <- getSequence(Protein_Seq, as.string=TRUE)
data("BLOSUM62")
protein_dat <- read.table("Protein_sequneces.dat")

Smith_Waterman_Scores <- data.frame(Seq1=as.numeric(),
                                    Seq2=as.numeric(), 
                                    Score=as.numeric(), 
                                    stringsAsFactors=FALSE) 

for (i in 1:100) {
  for (j in 1:100)
  {
    
    t <- pairwiseAlignment(protein_dat[i,], protein_dat[j,],substitutionMatrix=BLOSUM62,type="local")
    Smith_Waterman_Scores <- rbind(Smith_Waterman_Scores, c(i,j,t@score))
    
  }
}

names(Smith_Waterman_Scores)[1] <- "First.Protein"
names(Smith_Waterman_Scores)[2] <- "Second.Protein"
names(Smith_Waterman_Scores)[3] <- "Score"



########NORMALIZATION OF SMITH-WATERMAN SCORES########
######################################################
dt <- data.table(Smith_Waterman_Scores)
dt.lookup <- dt[First.Protein == Second.Protein]
setkey(dt,"First.Protein" )
setkey(dt.lookup,"First.Protein" )
colnames(dt.lookup) <- c("First.Protein","Second.Protein","Score1")
dt <- dt[dt.lookup]
setkey(dt,"Second.Protein" )
setkey(dt.lookup,"Second.Protein")
colnames(dt.lookup) <- c("First.Protein","Second.Protein","Score2")
dt <- dt[dt.lookup][
  , Normalized :=  Score / (sqrt(Score1) * sqrt(Score2))][
    , .(First.Protein, Second.Protein, Normalized)]
dt <- dt[order(dt$First.Protein),]
Smith_Waterman_Scores <- as.data.frame(dt)



########UNIPROT TO ENTREZID CONVERSION / BIOLOGICAL PROCESS-BASED PROTEIN SIMILARITY########
#############Gene Ontology (GO) annotations-based similarities##############################
############################################################################################
protein_id_list <- readLines("Protein_ids.txt")
protein_id_list_entrez_ids <- bitr(protein_id_list, 'UNIPROT', 'ENTREZID', annoDb='org.Hs.eg.db')
mgeneSim_calculated_similarities <- mgeneSim(protein_id_list_entrez_ids$ENTREZID, organism="human", ont="BP")



############GENERIC STRING KERNEL#############
######OUTPUT FROM PYTHON SCRIPT###############
GS_Kernel <- read.table("GS_KERNEL_SIMILARITIES")



########COMPUTING KRONECKER PRODUCTS########
reshaped_FP2_Tanimoto_fpSim_Search <- reshape(FP2_Tanimoto_fpSim_Search, idvar="First.Drug", timevar="Second.Drug", direction="wide")
reshaped_Smith_Waterman_Scores <- reshape(Smith_Waterman_Scores, idvar="First.Protein", timevar="Second.Protein", direction="wide")

FP3_Tanimoto <- reshape(FP3_Tanimoto_fpSim_Search, idvar="First.Drug", timevar="Second.Drug", direction="wide")
FP4_Tanimoto <- reshape(FP4_Tanimoto_fpSim_Search, idvar="First.Drug", timevar="Second.Drug", direction="wide")


row.names(reshaped_Smith_Waterman_Scores) <- reshaped_Smith_Waterman_Scores$First.Protein
reshaped_Smith_Waterman_Scores <- reshaped_Smith_Waterman_Scores[,-1]

row.names(reshaped_FP2_Tanimoto_fpSim_Search) <- reshaped_FP2_Tanimoto_fpSim_Search$First.Drug
reshaped_FP2_Tanimoto_fpSim_Search <- reshaped_FP2_Tanimoto_fpSim_Search[,-1]

row.names(FP3_Tanimoto) <- FP3_Tanimoto$First.Drug
FP3_Tanimoto <- FP3_Tanimoto[,-1]

row.names(FP4_Tanimoto) <- FP4_Tanimoto$First.Drug
FP4_Tanimoto <- FP4_Tanimoto[,-1]

row.names(reshaped_FP2_Tanimoto_fpSim_Search) <- drug_id$V1
names(reshaped_FP2_Tanimoto_fpSim_Search) <- drug_id$V1

row.names(FP3_Tanimoto) <- drug_id$V1
names(FP3_Tanimoto) <- drug_id$V1

row.names(FP4_Tanimoto) <- drug_id$V1
names(FP4_Tanimoto) <- drug_id$V1

row.names(reshaped_Smith_Waterman_Scores) <- protein_id$V1
names(reshaped_Smith_Waterman_Scores) <- protein_id$V1

row.names(GS_Kernel) <- protein_id$V1
names(GS_Kernel) <- protein_id$V1

reshaped_Smith_Waterman_Scores_MATRIX <- as.matrix(reshaped_Smith_Waterman_Scores)
reshaped_FP2_Tanimoto_fpSim_Search_MATRIX <- as.matrix(reshaped_FP2_Tanimoto_fpSim_Search)
FP3_Tanimoto_MATRIX <- as.matrix(FP3_Tanimoto)
FP4_Tanimoto_MATRIX <- as.matrix(FP4_Tanimoto)
GS_Kernel_MATRIX <- as.matrix(GS_Kernel)

FP2_based_KroneckerProducts <- kronecker(reshaped_FP2_Tanimoto_fpSim_Search_MATRIX, reshaped_Smith_Waterman_Scores_MATRIX, make.dimnames=TRUE)
FP3_based_KroneckerProducts <- kronecker(FP3_Tanimoto_MATRIX, reshaped_Smith_Waterman_Scores_MATRIX, make.dimnames=TRUE)
FP4_based_KroneckerProducts <- kronecker(FP4_Tanimoto_MATRIX, reshaped_Smith_Waterman_Scores_MATRIX, make.dimnames=TRUE)

FP2_GSK_KroneckerProducts <- kronecker(reshaped_FP2_Tanimoto_fpSim_Search_MATRIX, GS_Kernel_MATRIX,make.dimnames = TRUE)
FP3_GSK_KroneckerProducts <- kronecker(FP3_Tanimoto_MATRIX, GS_Kernel_MATRIX,make.dimnames = TRUE)
FP4_GSK_KroneckerProducts <- kronecker(FP4_Tanimoto_MATRIX,GS_Kernel_MATRIX,make.dimnames = TRUE)



#########FP2/SW CROSS VALIDATION#########
#########################################
df_matrix <- FP2_based_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP2_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
#  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP2_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP2_complete_matrix <- t(FP2_complete_matrix)
FP2_sorted <- FP2_complete_matrix[order(rownames(FP2_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



#########FP3/SW CROSS VALIDATION#########
#########################################
df_matrix <- FP3_based_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP3_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
  #  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP3_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP3_complete_matrix <- t(FP3_complete_matrix)
FP3_sorted <- FP3_complete_matrix[order(rownames(FP3_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



#########FP4/SW CROSS VALIDATION#########
#########################################
df_matrix <- FP4_based_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP4_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
  #  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP4_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP4_complete_matrix <- t(FP4_complete_matrix)
FP4_sorted <- FP4_complete_matrix[order(rownames(FP4_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



#########FP2/GSK CROSS VALIDATION#########
##########################################
df_matrix <- FP2_GSK_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP2_GSK_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
  #  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP2_GSK_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP2_GSK_complete_matrix <- t(FP2_GSK_complete_matrix)
FP2_GSK_sorted <- FP2_GSK_complete_matrix[order(rownames(FP2_GSK_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



#########FP3/GSK CROSS VALIDATION#########
##########################################
df_matrix <- FP3_GSK_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP3_GSK_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
  #  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP3_GSK_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP3_GSK_complete_matrix <- t(FP3_GSK_complete_matrix)
FP3_GSK_sorted <- FP3_GSK_complete_matrix[order(rownames(FP3_GSK_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



#########FP4/GSK CROSS VALIDATION#########
#########################################
df_matrix <- FP4_GSK_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 5
folds <- cvFolds(nrow(df_matrix), K=k)
CV_FP4_GSK_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_one <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_two <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_three <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_four <- matrix(0, nrow=nrow(df_matrix),ncol=1)
set_five <- matrix(0, nrow=nrow(df_matrix),ncol=1)
for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
  #  predictions[folds$subsets[folds$which == i]] <- predict%*%alpha
  if(i==1) {set_one <- crossprod(alpha,t(predict))}
  if(i==2) {set_two <- crossprod(alpha,t(predict))}
  if(i==3) {set_three <- crossprod(alpha,t(predict))}
  if(i==4) {set_four <- crossprod(alpha,t(predict))}
  if(i==5) {set_five <- crossprod(alpha,t(predict))}
}
FP4_GSK_complete_matrix <- cbind(set_one,set_two,set_three,set_four,set_five)
FP4_GSK_complete_matrix <- t(FP4_GSK_complete_matrix)
FP4_GSK_sorted <- FP4_GSK_complete_matrix[order(rownames(FP4_GSK_complete_matrix)), ]
df_vector <- as.matrix(y)
df_vector <- na.omit(df_vector)
new_sorted <- df_vector[order(rownames(df_vector)),]
new_sorted <- as.matrix(new_sorted)



######FP3/GSK 50-FOLD-CROSS VALIDATION########
##############################################
ptm <- proc.time()
df_matrix <- FP3_GSK_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]

df_vector <- y
df_vector <- na.omit(df_vector)

k = 50
folds <- cvFolds(nrow(df_matrix), K=k)
FiftyFoldFP3GSK_predictions <- matrix(0, nrow=nrow(df_matrix),ncol=1)

for(i in 1:k){
  train <- df_matrix[folds$subsets[folds$which != i], folds$subsets[folds$which != i]]
  y_train <- df_vector[folds$subsets[folds$which!=i]]
  test <- df_matrix[folds$subsets[folds$which==i], folds$subsets[folds$which!=i]]
  y_test <- df_vector[folds$subsets[folds$which==i]]
  
  fold_size = nrow(train)
  I = diag(fold_size)
  lambda = 0.1
  alpha <- solve(train + lambda * I) %*% y_train
  
  predict <- as.matrix(test)
 FiftyFoldFP3GSK_predictions[folds$subsets[folds$which == i]] <- predict%*%alpha

}
proc.time() - ptm



#############INTEGRATING ALL KERNEL MODELS TO GET THE HIGHEST PRECISION POSSIBLE###########
#########FOLLOWING THE APPROACH PROPOSED BY KronRLS-MKL####################################
combined <- cbind(new_sorted, FP2_sorted,FP3_sorted,FP4_sorted,FP2_GSK_sorted,FP3_GSK_sorted,FP4_GSK_sorted)
m1 <- combined
t <- m1[,-1][cbind(1:nrow(m1), max.col(-abs(m1[,-1]-m1[,1])))]
tt <- cbind(combined[,1],t)
summary(tt)
rmse(tt[,1],tt[,2])
cor(tt)



#########PREDICTION OF NaN VALUES IN THE DTI DATASET##########
##############################################################
y <- unmatrix(as.matrix(data),byrow=FALSE)
y <- as.matrix(y)

df_matrix <- FP3_GSK_KroneckerProducts
df_vector <- y

df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]
df_vector <- na.omit(df_vector)

I <- diag(nrow(df_matrix))
lambda <- 0.1
alpha <- solve(lambda * I + df_matrix) %*% df_vector
K_test <- FP3_GSK_KroneckerProducts[, !(colnames(KFP_3GSK_roneckerProducts) %in% rownames(df_vector))]
K_test <- K_test[(rownames(K_test) %in% rownames(alpha)),]
PREDICTIONS_BASED_ON_FP3_GSK <- crossprod(alpha,K_test)



#########3-FOLD-CROSS VALIDATION OF FP2/SW#########
###################################################
df_vector <- y
df_vector <- na.omit(df_vector)
folding <- df_vector[1:2740,]
folding <- as.matrix(folding)
df_matrix <- FP2_based_KroneckerProducts
df_vector <- y
df_matrix <- df_matrix[,!is.nan(df_vector[,1])]
df_matrix <- df_matrix[!is.nan(df_vector[,1]),]
df_vector <- na.omit(df_vector)
K_training <- df_matrix[!(rownames(df_vector) %in% rownames(folding)), !(rownames(df_vector) %in% rownames(folding))]
y_fold_vector <- df_vector[!(rownames(df_vector) %in% rownames(folding)),]
y_fold_vector <- as.matrix(y_fold_vector)
I_fold <- diag(nrow(K_training))
alpha_fold_FP2 <- solve(lambda * I_fold + K_training) %*% y_fold_vector
fold_test <- df_matrix[rownames(df_matrix) %in% rownames(K_training), colnames(df_matrix) %in% rownames(folding)]
fold_predictions_FP2 <- crossprod(alpha_fold_FP3,fold_test)



#########PLUGGING IN PREDICTED VALUES#########
##############################################
a=t(data_matrix)
a[is.na(a)]=Predictions_FP2[1,]
plugged = t(a)
heatmap.2(plugged,col=redgreen(1500), scale="none", margins=c(4,9), density.info="none", trace="none",dendrogram = "none",Rowv = "none",Colv = "none")
heatmap.2(data_matrix,col=redgreen(1500), scale="none", margins=c(4,9), density.info="none", trace="none",dendrogram = "none",Rowv = "none",Colv = "none")