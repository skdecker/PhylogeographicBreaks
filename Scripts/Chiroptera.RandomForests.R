##Conduct Random Forests analysis on Chiroptera data to find predictors of phylogeographic breaks
#Load neccessary packages
library(randomForest)
library(caret)
library(mlbench)
library(Hmisc)
library(ggplot2)

#Set working directory
setwd("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/")

#Read in and check data
all<-read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/impute_30_withcrypticsIBD.csv", header=T) #this is the data with imputed missing values
colnames(all)

#fix classes for random forest
str(all) #check structure of data
all$Break=as.factor(all$Break)
all$PredCryp=as.factor(all$PredCryp)
all$IBD_Sig=as.factor(all$IBD_Sig)
all$IBD_pval=as.numeric(all$IBD_pval)
all$Nseq=as.numeric(all$Nseq)
all$Family=as.factor(all$Family)
all$DietBreadth_2=as.factor(all$DietBreadth_2)
all$TrophLev_2=as.factor(all$TrophLev_2)
all$ForagingStrat_2=as.factor(all$ForagingStrat_2)
all$UpperElev_m_2=as.numeric(all$UpperElev_m_2)
all$DisectMount_2=as.factor(all$DisectMount_2)
all$BiogeogRealm_2=as.factor(all$BiogeogRealm_2)
all$HabitatBreadth_2=as.factor(all$HabitatBreadth_2)

#remove variables we will not use (Species name, Gene)
nosp=all[,c(3:47)]

#------------------------------------------------------------------------------
#perform Recursive Feature Elimination (RFE) to identify the best set of variables for the model
set.seed(391)
x = as.data.frame(nosp[,c(1,5:45)]) #predictor variables
y = nosp[,(2)] #response variable
subsets = c(1:25) #try subsets of up to 25 variables
control = rfeControl(functions = rfFuncs, method = "repeatedCV", repeats = 5, verbose = FALSE) #set control for RFE
lmtry = rfe(x, y, sizes=subsets, rfeControl = control) #run RFE
lmtry #view results of RFE
to_keep = predictors(lmtry) #save predictor variables from RFE

#------------------------------------------------------------------------------
#Build 50 Random Forests classifiers to identify traits important in predicting phylogeographic breaks
randomnumbers = ceiling(runif(50, 1, 1000)) #set 50 random seeds
#loop through 50 random seed numbers
for (r in 1:length(randomnumbers)){
  
  #split into traning and testing datasets
  set.seed(randomnumbers[r]) #set seed based on random numbers generated above
  trainIndex<-createDataPartition(nosp$Break, p=.85, list=F, times=1) #partition data - 85% training, 15% testing
  training<-nosp[trainIndex,] #training data from indexes
  testing<-nosp[-trainIndex,] #testing data from indexes
  
  #set up training control using cross validation
  ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                       classProbs = TRUE,
                       summaryFunction = twoClassSummary)

  #Build classifier with variables selected by RFE, 2000 decision trees, and cross validation
  set.seed(randomnumbers[r]) #set seed based on random numbers generated above
  all_rf <- train(Break ~ Area + Max_Lat + PredCryp + AdultBrainMass_g_2  + Length_Lat 
                  + DisectMount_2 + WingLoad_Nm2_3 + TrophLev_2 
                  + Min_Lon + Wingspan_cm_3 + AdultMass_g_2 
                  + AETMean_mm_1 + BiogeogRealm_2 + HumPopDensChange_1 
                  + Family, data = training,
                  method = "rf",
                  ntree = 2000,
                  importance=TRUE,
                  na.action=na.omit,
                  metric = "ROC",
                  trControl = ctrl) 
  
  #Get final model from cross validation
  modeldf = as.data.frame(all_rf$finalModel$err.rate) #cross validation error rate
  modelmeans = colMeans(modeldf) #take mean of columns
  write.table(data.frame(randomnumbers[r], unname(modelmeans[1]), unname(modelmeans[2]), unname(modelmeans[3])), file="Chiroptera_trainingerror.csv", col.names=FALSE, sep=",", append=TRUE)
  
  #Use the final model to predict the test dataset
  pred<-predict(all_rf, newdata = testing)
  
  #Generate a confusion matrix from the predictions
  CM=confusionMatrix(data = pred, testing$Break , positive = "Y") #save confusion matrix
  CMdf=data.frame(t(as.data.frame(c(CM$overall, CM$byClass)))) #as a dataframe
  write.table(CMdf, file="Chiroptera_confusionmatrix.csv", sep=",", append=TRUE)
  
  #Get importance values from final model
  all_imp<-importance(all_rf$finalModel) #save variable importance
  write.table(as.data.frame(all_imp), file="Chiroptera_variableimportance.csv", sep=",", append=TRUE)
}

#------------------------------------------------------------------------------
#Build RF classifier with IBD result as response variable
#split into training and testing datasets
set.seed(391)
trainIndex<-createDataPartition(all$IBD_Sig, p=.85, list=F, times=1) #partition data - 85% training; 15% testing
training<-all[trainIndex,] #set training dataset
testing<-all[-trainIndex,] #set testing dataset

#set up training control using cross validation
ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)

#Build RF classifier using same predictor variables as phylogeographic breaks model
set.seed(391)
all_rf <- train(IBD_Sig ~ Area + Max_Lat + PredCryp + AdultBrainMass_g_2  + Length_Lat 
                + DisectMount_2 + WingLoad_Nm2_3 + TrophLev_2 
                + Min_Lon + Wingspan_cm_3 + AdultMass_g_2 
                + AETMean_mm_1 + BiogeogRealm_2 + HumPopDensChange_1 
                + Family, data = training,
                method = "rf",
                ntree = 2000,
                importance=TRUE,
                na.action=na.omit,
                metric = "ROC",
                trControl = ctrl)

#Get final model from cross validation
modeldf = as.data.frame(all_rf$finalModel$err.rate) #cross validation error rate
modelmeans = colMeans(modeldf) #means of columns
write.table(data.frame(unname(modelmeans[1]), unname(modelmeans[2]), unname(modelmeans[3])), file="Chir.IBDmodel.csv", col.names=FALSE, sep=",", append=TRUE)

#Use the final model to predict the test dataset
pred<-predict(all_rf, newdata = testing)

#Generate a confusion matrix from the predictions
CM=confusionMatrix(data = pred, testing$IBD_Sig , positive = "Y") #save test data confusion matrix
CMdf=data.frame(t(as.data.frame(c(CM$overall, CM$byClass)))) #as data frame
write.table(CMdf, file="Chir.IBDconfusionmat.csv", sep=",", append=TRUE)

#Get importance values from final model
all_imp<-importance(all_rf$finalModel) #save variable importance from final IBD RF model
write.table(as.data.frame(all_imp), file="Chir.IBDvarimp_reg.csv", sep=",", append=TRUE)