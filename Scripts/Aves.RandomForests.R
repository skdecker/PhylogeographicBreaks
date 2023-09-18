#RF with IBD as response variable - birds
library(randomForest)
library(caret)
library(mlbench)
library(Hmisc)
library(ggplot2)

setwd("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/")

#read data
birds.ibd<-read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/Aves_Traits.wIBD.csv", header=T)
colnames(birds.ibd)
attach(birds.ibd)

str(birds.ibd)

#fix classes for random forest
birds.ibd$Gene = as.factor(birds.ibd$Gene)
birds.ibd$Break = as.factor(birds.ibd$Break)
birds.ibd$IBD.P = as.numeric(birds.ibd$IBD.P)
birds.ibd$IBD.Sig = as.factor(birds.ibd$IBD.Sig)
birds.ibd$Family1 = as.factor(birds.ibd$Family1)
birds.ibd$Order1 = as.factor(birds.ibd$Order1)
birds.ibd$Habitat = as.factor(birds.ibd$Habitat)
birds.ibd$Habitat.Density = as.factor(birds.ibd$Habitat.Density)
birds.ibd$Migration = as.factor(birds.ibd$Migration)
birds.ibd$Trophic.Level = as.factor(birds.ibd$Trophic.Level)
birds.ibd$Trophic.Niche = as.factor(birds.ibd$Trophic.Niche)
birds.ibd$Primary.Lifestyle = as.factor(birds.ibd$Primary.Lifestyle)

#-------------------------------------------------------------------------------

#RFE
set.seed(391)
x = as.data.frame(birds.ibd[,c(8, 12:13, 20, 27, 29:37, 43:48, 52, 53)]) #removes highly correlated variables
y = birds.ibd[,(4)]
subsets = c(1:22)

control = rfeControl(functions = rfFuncs, method = "repeatedCV", repeats = 5, verbose = FALSE)

lmtry = rfe(x, y, sizes=subsets, rfeControl = control)
lmtry #Keep all 22
to_keep = predictors(lmtry)


#-------------------------------------------------------------------------------
#HEEEEEEERRRRRRRREEEEE *****************
#Random Forests with the 22 variables (after RFE and removal of (non-morpho) correlates > .7) 
#With Break as response variable and 50 independent model iterations
#split into training and testing datasets
randomnumbers #= ceiling(runif(50, 1, 1000)) #use the same random seed numbers as used for bat RF models
for (r in 1:length(randomnumbers)){
set.seed(randomnumbers[r])
trainIndex<-createDataPartition(birds.ibd$Break, p=.85, list=F, times=1)
training<-birds.ibd[trainIndex,]
testing<-birds.ibd[-trainIndex,]

#SUBSAMPLING INSIDE OF RESAMPLING (CROSS_VALIDATION):
ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
colnames(training)
#run rf with cross-validation + HabitatBreadth_2 + AdultLength_mm_2
set.seed(randomnumbers[r])
all_rf <- train(Break ~ Area + AbsMax_Lat + AbsMin_Lat + Order1 + Beak.Length_Nares + Beak.Length_Culmen
                + Beak.Width + Beak.Depth + Tarsus.Length + Wing.Length + Kipps.Distance
                + Secondary1 + Hand.Wing.Index + Tail.Length + Mass + Habitat + Habitat.Density
                + Migration + Trophic.Level + Trophic.Niche + Primary.Lifestyle + Centroid.Longitude
                + Range.Size, data = training,
                method = "rf",
                ntree = 2000,
                importance=TRUE,
                na.action=na.omit,
                metric = "ROC",
                trControl = ctrl)



#Get final model from cross validation
modeldf = as.data.frame(all_rf$finalModel$err.rate)
modelmeans = colMeans(modeldf)
write.table(data.frame(randomnumbers[r], unname(modelmeans[1]), unname(modelmeans[2]), unname(modelmeans[3])), file="Aves_TrainError_50iter.csv", col.names=FALSE, sep=",", append=TRUE)

#Use the final model to predict the test dataset
pred<-predict(all_rf, newdata = testing)

#Generate a confusion matrix from the predictions
CM=confusionMatrix(data = pred, testing$Break, positive = "Y")
CMdf=data.frame(t(as.data.frame(c(CM$overall, CM$byClass))))
write.table(CMdf, file="Aves_ConfMat_50iter.csv", sep=",", append=TRUE)

#Get importance values from final model
all_imp<-importance(all_rf$finalModel)
all_imp
write.table(as.data.frame(all_imp), file="Aves_VarImp_50iter.csv", sep=",", append=TRUE)
}

#-------------------------------------------------------------------------------

#Random Forests with the 22 variables (after RFE and removal of non-morpho correlates > .7) 
#With IBD as response variable
#split into traning and testing datasets
set.seed(391)
trainIndex<-createDataPartition(birds.ibd$IBD.Sig, p=.85, list=F, times=1)
training<-birds.ibd[trainIndex,]
testing<-birds.ibd[-trainIndex,]

#SUBSAMPLING INSIDE OF RESAMPLING (CROSS_VALIDATION):
ctrl <- trainControl(method = "repeatedcv", repeats = 5,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
colnames(training)
#run rf with cross-validation + HabitatBreadth_2 + AdultLength_mm_2
set.seed(391)
all_rf <- train(IBD.Sig ~ Area + AbsMax_Lat + AbsMin_Lat + Order1 + Beak.Length_Culmen + Beak.Length_Nares
                + Beak.Width + Beak.Depth + Tarsus.Length + Wing.Length + Kipps.Distance
                + Secondary1 + Hand.Wing.Index + Tail.Length + Mass + Habitat + Habitat.Density
                + Migration + Trophic.Level + Trophic.Niche + Primary.Lifestyle + Centroid.Longitude
                + Range.Size, data = training,
                method = "rf",
                ntree = 2000,
                importance=TRUE,
                na.action=na.omit,
                metric = "ROC",
                trControl = ctrl)



#Get final model from cross validation
modeldf = as.data.frame(all_rf$finalModel$err.rate)
modelmeans = colMeans(modeldf)
write.table(data.frame(unname(modelmeans[1]), unname(modelmeans[2]), unname(modelmeans[3])), file="Aves.IBDmodel.csv", col.names=FALSE, sep=",", append=TRUE)

#Use the final model to predict the test dataset
pred<-predict(all_rf, newdata = testing)

#Generate a confusion matrix from the predictions
CM=confusionMatrix(data = pred, testing$IBD.Sig , positive = "Y")
CMdf=data.frame(t(as.data.frame(c(CM$overall, CM$byClass))))
write.table(CMdf, file="Aves.IBDconfusionmat.csv", sep=",", append=TRUE)

#Get importance values from final model
all_imp<-importance(all_rf$finalModel)
all_imp
write.table(as.data.frame(all_imp), file="Aves.IBDvarimp.csv", sep=",", append=TRUE)
