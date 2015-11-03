## Clear the workspace and import the data
rm(list = ls())
ptm <- proc.time()

#### load required libraries
library(PharmacoGx)
library(Biobase)
library(survcomp)
## load the data (Challenge data and L1000)
load("~/Desktop/PharmacoGx/data/Challenge.RData")
load("~/Desktop/PharmacoGx/data/L1000.RData")

#### Input parameters
IC50_CutOff        <- 100     ## The cut off for the IC50 (this criteria does not affect the training much)
polDegree          <- 2       ## degree of polynomial for the gene expression data
CVItermax          <- 1000    ## Number of iterations for cross validation (and training)
                              ## among all the data
CV_frac            <- 0.2     ## fraction of data for cross validation
FeatNum_Gene       <- 5      ## Number of features from gene expression data
FeatNum_Methyl     <- 5      ## Number of features from methylation data
TopMethyl_Num      <- 1000    
################### Building the input matrix and output vector
###################

## Define new matrices to store gene expression and gene names
CellName_Methylation <- colnames(Challenge@molecularProfiles$methylation.MA)
CellNum_Methylation <- ncol(Challenge@molecularProfiles$methylation.MA)

### Etracting methylation data for highly variance methylation
allMethyl_Num  <- nrow(exprs(Challenge@molecularProfiles$methylation.MA))
MethylationVar <- apply(exprs(Challenge@molecularProfiles$methylation.MA), 1, var)
Ind_TopMethylations <- order(MethylationVar, decreasing=T)[1:TopMethyl_Num]

CellMethyl_Feature = list()
for(i in 1:CellNum_Methylation){
  CellMethyl_Feature[[CellName_Methylation[i]]] <- t(exprs(
    Challenge@molecularProfiles$methylation.MA)[,i])
}

## Define new matrices to store gene expression and gene names
CellName_RNA <- pData(Challenge@molecularProfiles$rna)[,"cellid"] #### ID
CellNum      <- nrow(pData(Challenge@molecularProfiles$rna))

## Obtain average/minimum/maximum gene expression for 
## each gene between different cell lines

Ind_TopGenes <- which(fData(Challenge@molecularProfiles$rna)
                      [,"Symbol"] %in% L1000[,"pr_gene_symbol"])
TopGene_Num  <- length(Ind_TopGenes)
allGene_Num  <- nrow(exprs(Challenge@molecularProfiles$rna))

## Gene Expression matrix for the cells and define a cell 
## feature based on Gene Expression data
CellExp_Feature <- list()
for(i in 1:CellNum){
  CellExp_Feature[[CellName_RNA[i]]] <- t(exprs(Challenge@molecularProfiles$rna)[,i])
}

## New vector for monotherapy and drug combination information
DrugDrugCell_Syn  <- list()
DrugCell_AUC      <- list()
DrugCell_IC50     <- list()

DrugCombGiven_Syn <- c()
CellName_Comb     <- c()
Drug_Name_Comb    <- c()
DrugDrugCell_QC   <- c()

for(i  in 1:nrow(Challenge@sensitivity$info)){
  CellDrug_Name <- paste(Challenge@sensitivity$info[i,"cellid"], 
                         Challenge@sensitivity$info[i,"drugid"], sep="---")
  ## if it is one drug: CellName---DrugName   
  ## if it is drug combination CellName---Drug1///Drug2
  if(grepl("///",Challenge@sensitivity$info[i,"drugid"])){
    DrugDrugCell_Syn[[CellDrug_Name]] <- Challenge@sensitivity$profiles[i,"Synergy_score"]
    
    if(!is.na(Challenge@sensitivity$profiles[i,"Synergy_score"]) && 
       Challenge@sensitivity$info[i,"quality"] == "1"){
        DrugCombGiven_Syn <- c(DrugCombGiven_Syn, 
                               Challenge@sensitivity$profiles[i,"Synergy_score"])
        CellName_Comb     <- c(CellName_Comb, Challenge@sensitivity$info[i,"cellid"])
        Drug_Name_Comb    <- c(Drug_Name_Comb, Challenge@sensitivity$info[i,"drugid"])        
    }
    
  }else{
    DrugCell_AUC[[CellDrug_Name]]  <- Challenge@sensitivity$profiles[i,"auc_recomputed"]
    DrugCell_IC50[[CellDrug_Name]] <- Challenge@sensitivity$profiles[i,"ic50_recomputed"]
  }
}

## building input matrix for regression


## input data include 2) Cell gene expression feature, 
## 3) Drug1-Celll AUC, 4) Drug1-Cell IC50,
## 5) Drug2-Cell AUC, 6) Drug2-Cell IC50
## First column of the matrix shoudl be all equal to 1
# InputMat_col1 <- matrix(rep(1:TestData_length,1), ncol = 1)

## output vecotr also is obtained here
InputMat_col1 <- c()
InputMat_col2 <- c()
InputMat_col3 <- c()
InputMat_col4 <- c()
InputMat_col5 <- c()
OutputData    <- c()

##### Building different columns of the input matrix

##### If you want to put your features it is just required to change the 
##### indeces of the genes and methylation in InputGeneExp and 
#### InputMethylation, respectively
DrugComb_length  <- length(DrugCombGiven_Syn)
InputGeneExp     <- sample(1:allGene_Num,FeatNum_Gene,replace=F)
InputMethylation <- sample(1:allMethyl_Num,FeatNum_Methyl,replace=F)

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
      InputMat_col3 <- rbind(InputMat_col3, IC50_CutOff)
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
NcolExp     <- ncol(InputMat_col2)
NcolMethyl  <- ncol(InputMat_col5)
InputMatrix <- c()

# Adding gene expression features to the input matrix
if(polDegree > 1){
  InputMatrix <- cbind(InputMatrix, poly(matrix(InputMat_col2, ncol = NcolExp), 
                      degree = polDegree, raw = TRUE))
}else{
  InputMatrix <- cbind(InputMatrix, matrix(InputMat_col2, ncol = NcolExp))
}

# Adding methylation features to the input matrix
if(polDegree > 1){
  InputMatrix <- cbind(InputMatrix, poly(matrix(InputMat_col5, ncol = NcolMethyl), 
                                     degree = polDegree, raw = TRUE))
}else{
  InputMatrix <- cbind(InputMatrix, matrix(InputMat_col5, ncol = NcolMethyl))
}
#### Input matrix for training and cross validation
InputMatrix <- cbind(InputMatrix, matrix(InputMat_col3, ncol = 1), 
                     matrix(InputMat_col4, ncol = 1))

################### Training and post-processing
###################
###################

### Seperating training input matrix and outpur vector
CIndex_CV_Vector    <- c()
CIndex_Train_Vector <- c()

OutputTrain_Vec   <- c()
OutputCV_Vec      <- c()
PredictTrain_Vec  <- c()
PredictCV_Vec     <- c()

for(CViter in 1:CVItermax){
  CvRand_Ind  <- sample(1:length(OutputData),
                        floor((1-CV_frac)*length(OutputData)),replace=F)
  TrainOutput <- OutputData[CvRand_Ind]
  CrossValOutput <- OutputData[-CvRand_Ind]
  ### linear regression model for training
  lm.Train <- lm(TrainOutput ~ InputMatrix[CvRand_Ind,])
  
  #### Obtain the prediction values for the training data
  #### and compare them with the observed values
  TrainCoeff      <- matrix(coef(summary(lm.Train))[, 1], ncol = 1)
  InputMatrix_aux <- cbind(matrix(InputMat_col1, ncol = 1), InputMatrix)
  
  Train_Pred      <- matrix(fitted(lm.Train), ncol = 1)
  Train_Obse      <- TrainOutput
  
  #### Obtain the prediction values for the cross validation data
  #### and compare them with the observed values
  CrossVal_Pred      <- InputMatrix_aux[-CvRand_Ind,] %*% TrainCoeff
  CrossVal_Obse      <- CrossValOutput
  
  ##### C Index calculartion
  CIndex_CrossVal <- (1 - concordance.index(CrossVal_Obse, 
                                            CrossVal_Pred, rep(1, length(CrossVal_Obse))
                                            , na.rm = TRUE)[[1]])
  CIndex_Train    <- (1 - concordance.index(Train_Obse, 
                                            Train_Pred, rep(1, length(Train_Obse))
                                            , na.rm = TRUE)[[1]])
  ### Cross validation among all the data
  CIndex_CV_Vector    <-  c(CIndex_CV_Vector, CIndex_CrossVal)
  CIndex_Train_Vector <-  c(CIndex_Train_Vector, CIndex_Train)
  OutputTrain_Vec     <-  c(OutputTrain_Vec, Train_Obse)
  OutputCV_Vec        <-  c(OutputCV_Vec, CrossVal_Obse)
  PredictTrain_Vec    <-  c(PredictTrain_Vec, Train_Pred)
  PredictCV_Vec       <-  c(PredictCV_Vec, CrossVal_Pred)
}

## post-processing
minRandCI_CV   <- min(CIndex_CV_Vector)
meanRandCI_CV  <- mean(CIndex_CV_Vector)
maxRandCI_CV   <- max(CIndex_CV_Vector)

minRandCI_Train   <- min(CIndex_Train_Vector)
meanRandCI_Train  <- mean(CIndex_Train_Vector)
maxRandCI_Train   <- max(CIndex_Train_Vector)

### print CInd for training and cross validation data
cat("Minimum CIndex for random cross validation:", minRandCI_CV,  "\n")
cat("Average CIndex for random cross validation:", meanRandCI_CV, "\n")
cat("Maximum CIndex for random cross validation:", maxRandCI_CV,  "\n")

cat("Minimum CIndex for training data:", minRandCI_Train,  "\n")
cat("Average CIndex for training data:", meanRandCI_Train, "\n")
cat("Maximum CIndex for training data:", maxRandCI_Train,  "\n")

if(length(OutputTrain_Vec) < 1.2e5 && length(OutputCV_Vec) < 1.2e5){
  AllIteration_CI_Train <- (1 - concordance.index(OutputTrain_Vec,
                                                  PredictTrain_Vec, rep(1, length(OutputTrain_Vec))
                                                  , na.rm = TRUE)[[1]])
  
  AllIteration_CI_CV <- (1 - concordance.index(OutputCV_Vec, 
                                               PredictCV_Vec, rep(1, length(OutputCV_Vec))
                                               , na.rm = TRUE)[[1]])
  cat("CIndex for All Data cross validation:", AllIteration_CI_CV,  "\n")
  
  cat("CIndex for All Data training:", AllIteration_CI_Train,  "\n")
}

proc.time() - ptm