## Clear the workspace and import the data
rm(list = ls())
ptm <- proc.time()

#### load required libraries
library(PharmacoGx)
library(Biobase)
library(survcomp)
#### Input parameters
IC50_CutOff        <- 100     ## The cut off for the IC50 (this criteria does not affect the training much)
polDegree          <- 1      ## degree of polynomial for the gene expression data
TrainingFrac       <- 0.6     ## Fraction of data we want to use for training 
## (The rest will be used for cross validation)
RandCVResNum       <- 100     ## Number of random results to obtain RMS of 
## random results for CV data
CVItermax          <- 5000      ## Number of iterations for cross validation 
## among all the data
CV_frac             <- 0.2     ## fraction of data for cross validation
## load the data (Challenge data and L1000)
load("~/Desktop/PharmacoGx/data/Challenge.RData")
load("~/Desktop/PharmacoGx/data/L1000.RData")

################### Building the input matrix and output vector
###################
###################

## Define new matrices to store gene expression and gene names
CellName_Methylation <- colnames(Challenge@molecularProfiles$methylation.BS)
CellNum_Methylation <- length(colnames(Challenge@molecularProfiles$methylation.BS))

### Etracting methylation data for highly variance methylation
allMethyl_Num  <- length(exprs(
  Challenge@molecularProfiles$methylation.MS)[,1])
MethylationVar <- apply(exprs(Challenge@molecularProfiles$methylation.BS), 1, var)
Ind_TopMethylations <- which(MethylationVar >= sort(MethylationVar, decreasing=T)[1000], arr.ind=TRUE)

TopMethyl_Num  <- length(Ind_TopMethylations)


CellMethyl_Feature = list()
for(i in 1:CellNum_Methylation){
  CellMethyl_Feature[[CellName_Methylation[i]]] <- matrix(exprs(
    Challenge@molecularProfiles$methylation.MS)[Ind_TopMethylations,i], ncol = TopMethyl_Num)
}


## Define new matrices to store gene expression and gene names
CellName_RNA <- pData(Challenge@molecularProfiles$rna)[,24]
CellNum      <- length(CellName_RNA)

## Obtain average/minimum/maximum gene expression for 
## each gene between different cell lines

Ind_TopGenes <- which(fData(Challenge@molecularProfiles$rna)
                      [,4] %in% L1000[,3])
TopGene_Num  <- length(Ind_TopGenes)
allGene_Num  <- length(exprs(
  Challenge@molecularProfiles$rna)[,1])

## Gene Expression matrix for the cells and define a cell 
## feature based on Gene Expression data
CellExp_Feature = list()
for(i in 1:CellNum){
  CellExp_Feature[[CellName_RNA[i]]] <- matrix(exprs(
    Challenge@molecularProfiles$rna)[Ind_TopGenes,i], ncol = TopGene_Num)
}

### Highly significant features(***) among L1000 genes
SigFeat_GeneExp <- c(46,59,129,205,235,374,476,681,687,698,756,
                     764,777,799,878,1062)
Len_SigFeat     <- length(SigFeat_GeneExp)

### Siginificant features with higher pValue(**) among L1000 genes
SigFeat2_GeneExp <- c(9,15,18,49,50,67,83,94,96,118,123,157,161,
                      170,179,194,196,264,271,291,
                      302,345,411,446,450,457,464,477,490,500,535,
                      580,591,596,601,656,680,812,
                      813,953,958,960)
Len_SigFeat2     <- length(SigFeat2_GeneExp)

## New vector for monotherapy and drug combination information
DrugDrugCell_Syn  = list()
DrugCell_AUC      = list()
DrugCell_IC50     = list()

DrugCombGiven_Syn = c()
CellName_Comb     = c()
Drug_Name_Comb    = c()
DrugDrugCell_QC   = c()

for(i  in 1:length(Challenge@sensitivity$info[,1])){
  CellDrug_Name <- paste(Challenge@sensitivity$info[i,1], 
                         Challenge@sensitivity$info[i,2], sep="---")
  ## if it is one drug: CellName---DrugName   
  ## if it is drug combination CellName---Drug1///Drug2
  if(grepl("///",Challenge@sensitivity$info[i,2])){
    DrugDrugCell_Syn[[CellDrug_Name]] <- Challenge@sensitivity$profiles[i,4]
    
    if(!is.na(Challenge@sensitivity$profiles[i,4]) && 
       Challenge@sensitivity$info[i,9] == "1"){
      if(Challenge@sensitivity$profiles[i,4] > -200){
        DrugCombGiven_Syn <- c(DrugCombGiven_Syn, 
                               Challenge@sensitivity$profiles[i,4])
        CellName_Comb     <- c(CellName_Comb, Challenge@sensitivity$info[i,1])
        Drug_Name_Comb    <- c(Drug_Name_Comb, Challenge@sensitivity$info[i,2])        
      }
    }
    
  }else{
    DrugCell_AUC[[CellDrug_Name]]  <- Challenge@sensitivity$profiles[i,5]
    DrugCell_IC50[[CellDrug_Name]] <- Challenge@sensitivity$profiles[i,6]
  }
}

## building input matrix for regression


## input data include 2) Cell gene expression feature, 
## 3) Drug1-Celll AUC, 4) Drug1-Cell IC50,
## 5) Drug2-Cell AUC, 6) Drug2-Cell IC50
## First column of the matrix shoudl be all equal to 1
# InputMat_col1 <- matrix(rep(1:TestData_length,1), ncol = 1)

## output vecotr also is obtained here
InputMat_col1 = c()
InputMat_col2 = c()
InputMat_col3 = c()
InputMat_col4 = c()
InputMat_col5 = c()
OutputData    = c()

##### Building different columns of the input matrix
InputGeneExp     <- c(SigFeat_GeneExp, SigFeat2_GeneExp)
DrugComb_length  <- length(DrugCombGiven_Syn)
InputMethylation <- Ind_TopMethylations[1:5]

for(i in 1:DrugComb_length){
  DrugsNames <- strsplit(Drug_Name_Comb[i], "///")
  CellDrug1_Name <- paste(CellName_Comb[i], DrugsNames[[1]][1], sep="---")
  CellDrug2_Name <- paste(CellName_Comb[i], DrugsNames[[1]][2], sep="---")
  if(!is.null(CellExp_Feature[[CellName_Comb[i]]]) &&
     !is.null(CellMethyl_Feature[[CellName_Comb[i]]]) &&
     !is.na(DrugCell_AUC[[CellDrug1_Name]]) &&
     !is.na(DrugCell_IC50[[CellDrug1_Name]]) && 
     !is.na(DrugCell_AUC[[CellDrug2_Name]]) &&
     !is.na(DrugCell_IC50[[CellDrug2_Name]])){
    
    InputMat_col1 <- rbind(InputMat_col1, 1)
    InputMat_col2 <- rbind(InputMat_col2, 
                           CellExp_Feature[[CellName_Comb[i]]][1,InputGeneExp])
    
    if(DrugCell_IC50[[CellDrug1_Name]] > IC50_CutOff ){
      InputMat_col3 <- rbind(InputMat_col3, IC50_CutOff )
    }else{
      InputMat_col3 <- rbind(InputMat_col3, DrugCell_IC50[[CellDrug1_Name]])
    }
    
    if(DrugCell_IC50[[CellDrug2_Name]] > IC50_CutOff ){
      InputMat_col4 <- rbind(InputMat_col4, IC50_CutOff )
    }else{
      InputMat_col4 <- rbind(InputMat_col4, DrugCell_IC50[[CellDrug2_Name]])
    }
    
    InputMat_col5 <- rbind(InputMat_col5, 
                           CellMethyl_Feature[[CellName_Comb[i]]][1,InputMethylation])
    
    OutputData    <- rbind(OutputData, DrugCombGiven_Syn[i])
  }
}

## We try to normalize all the data to get a better convergence
Max_CellFeature <- max(abs(InputMat_col2))
MaxOutput       <- max(abs(OutputData))
OutputData      <- OutputData

### Make the complete polynomial of different degress  
### using matrix of gene expression data 
NcolExp     <- length(InputMat_col2[1,])
if(polDegree > 1){
  InputMatrix <- poly(matrix(InputMat_col2, ncol = NcolExp), degree = polDegree, raw = TRUE)
}else{
  InputMatrix <- matrix(InputMat_col2, ncol = NcolExp)
}

#### Input matrix for training and cross validation
NcolMethyl  <- length(InputMat_col5[1,])
InputMatrix <- cbind(InputMatrix, matrix(InputMat_col5, ncol = NcolMethyl),
                     matrix(InputMat_col3, ncol = 1), 
                     matrix(InputMat_col4, ncol = 1))

################### Training and post-processing
###################
###################

### Seperating training input matrix and outpur vector
CIndex_CV_Vector = c()
CIndex_Train_Vector = c()

for(CViter in 1:CVItermax){
  CvRand_Ind  <- sample(1:length(OutputData),floor((1-CV_frac)*length(OutputData)),replace=F)
  ### linear regression model for training
  lm.Train <- lm(OutputData[CvRand_Ind] ~ InputMatrix[CvRand_Ind,])
  
  #### Obtain the prediction values for the training data
  #### and compare them with the observed values
  TrainCoeff      <- matrix(coef(summary(lm.Train))[, 1], ncol = 1)
  InputMatrix_aux <- cbind(matrix(InputMat_col1, ncol = 1), InputMatrix)
  
  Train_Pred      <- InputMatrix_aux[CvRand_Ind,] %*% TrainCoeff
  Train_Obse      <- OutputData[CvRand_Ind]

  #### Obtain the prediction values for the cross validation data
  #### and compare them with the observed values
  CrossVal_Pred      <- InputMatrix_aux[-CvRand_Ind,] %*% TrainCoeff
  CrossVal_Obse      <- OutputData[-CvRand_Ind]

  ##### C Index calculartion
  CIndex_CrossVal <- (1 - concordance.index(CrossVal_Obse, 
                                            CrossVal_Pred, rep(1, length(CrossVal_Obse))
                                            , na.rm = TRUE)[[1]])
  CIndex_Train    <- (1 - concordance.index(Train_Obse, 
                                            Train_Pred, rep(1, length(Train_Obse))
                                            , na.rm = TRUE)[[1]])
  ### Cross validation among all the data
  CIndex_CV_Vector <- c(CIndex_CV_Vector, CIndex_CrossVal)
  CIndex_Train_Vector <-  c(CIndex_Train_Vector, CIndex_Train)
  
  
}
minRandCI <- min(CIndex_CV_Vector)
meanRandCI <- mean(CIndex_CV_Vector)
maxRandCI <- max(CIndex_CV_Vector)

### print the results
cat("Minimum CIndex for random cross validation:", minRandCI,  "\n")
cat("Average CIndex for random cross validation:", meanRandCI, "\n")
cat("Maximum CIndex for random cross validation:", maxRandCI,  "\n")

proc.time() - ptm