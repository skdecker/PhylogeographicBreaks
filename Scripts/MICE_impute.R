#load packages
library(VIM)
library(mice)

#read in dataset 
data<-read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/impute_30_wcryp.csv", na.strings="NA")
data = data[1:126,]
summary(data)

#### EVALUATE MISSING DATA ####
aggr_plot <- aggr(data, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))

pMiss <- function(x){sum(is.na(x))/length(x)*100}
miss_percent<-apply(data,2,pMiss)
miss_percent = as.matrix(miss_percent)
write.csv(miss_percent,file="/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/chir_missingdata.csv", row.names = T)

#### DATA IMPUTATION ####

#skip columns <- data[,which(miss_percent <= 50)]
#skip columns 
#skip aggr_plot <- aggr(columns, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(columns), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
#skip write.csv(columns, file="/Users/decker.391/Downloads/trimmed.csv")
#subset dataset for imputation 
##new dataset contains all columns with missing values + plus genus and family columns (to use as predictors in the MICE imputation)
#sal_data<-data[c(1:4,52:57,64:67,72:73,83,88:89,91:93)]

#run initial mice with no iterations
colnames(data)
data=data[,c(1:23,25:52)] #remove references columns
colnames(data)
data.nodup=data[,c(1:45)] #remove duplicate columns
summary(data.nodup)

aggr_plot <- aggr(data.nodup, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, labels=names(data.nodup), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
colnames(data.nodup)

init_sal<-mice(data.nodup, maxit = 0)

#set mice accessors
meth<-init_sal$method
predM = init_sal$predictorMatrix

#change values for imputation method 
##by default, the method uses:
  #pmm, predictive mean matching (numeric data) 
  #logreg, logistic regression imputation (binary data, factor with 2 levels) 
  #polyreg, polytomous regression imputation for un- ordered categorical data (factor > 2 levels) 
  #polr, proportional odds model for (ordered, > 2 levels).

data.nodup$DietBreadth_2=as.factor(data.nodup$DietBreadth_2)
data.nodup$TrophLev_2=as.factor(data.nodup$TrophLev_2)
data.nodup$DisectMount_2=as.factor(data.nodup$DisectMount_2)
data.nodup$HabitatBreadth_2=as.factor(data.nodup$HabitatBreadth_2)
data.nodup$Species=as.factor(data.nodup$Species)
data.nodup$Gene=as.factor(data.nodup$Gene)
data.nodup$Family=as.factor(data.nodup$Family)
data.nodup$ForagingStrat_2=as.factor(data.nodup$ForagingStrat_2)
data.nodup$BiogeogRealm_2=as.factor(data.nodup$BiogeogRealm_2)
data.nodup$PredCryp=as.factor(data.nodup$PredCryp)

meth[c("DietBreadth_2")]="polyreg"
meth[c("TrophLev_2")]="polyreg"
meth[c("DisectMount_2")]="logreg"
meth[c("HabitatBreadth_2")]="polyreg"
meth[c("PredCryp")]="logreg"

#change values for prediction matrix (what variables are used as predictors for imputation)
##the columns set by predM=1 are used as a predictor for the target block (in the rows)
##columns with predM=0 are not used for prediction in imputation

predM[, c("Species")]=0 #NA
predM[, c("Gene")]=0 #NA
predM[, c("Break")]=0 #NA
predM[, c("Nseq")]=0 #NA
predM[, c("Area")]=0 #NA
predM[, c("Max_Lat")]=0 #NA
predM[, c("Min_Lat")]=0 #NA
predM[, c("Med_Lat")]=0 #NA
predM[, c("Abs_Max_Lat")]=0 #NA
predM[, c("Abs_Min_Lat")]=0 #NA
predM[, c("Length_Lat")]=0 #NA
predM[, c("Max_Lon")]=0 #NA
predM[, c("Min_Lon")]=0 #NA
predM[, c("Med_Lon")]=0 #NA
predM[, c("Length_Lon")]=0 #NA
predM[, c("Family")]=1 #NA
predM[, c("GRArea_km2_1")]=0 #pmm
predM[, c("GRMaxLat_1")]=0 #pmm
predM[, c("GRMinLat_1")]=0 #pmm
predM[, c("GRMidLat_1")]=0 #pmm
predM[, c("GRMaxLong_1")]=0 #pmm
predM[, c("GRMinLong_1")]=0 #pmm
predM[, c("GRMidLong_1")]=0 #pmm
predM[, c("HumPopDensMean_n.km2_1")]=0 #pmm
predM[, c("HumPopDensChange_1")]=0 #pmm
predM[, c("PrecipMean_mm_1")]=0 #pmm
predM[, c("TempMean_degC_1")]=0 #pmm
predM[, c("AETMean_mm_1")]=0 #pmm
predM[, c("PETMean_mm_1")]=0 #pmm
predM[, c("AdultMass_g_2")]=0 #pmm
predM[, c("AdultBrainMass_g_2")]=0 #pmm
predM[, c("AdultLength_mm_2")]=0 #pmm
predM[, c("AdultFA_mm_2")]=0 #pmm
predM[, c("LitterSize_2")]=0 #pmm
predM[, c("DietBreadth_2")]=0 #polyreg
predM[, c("TrophLev_2")]=0 #polyreg
predM[, c("ForagingStrat_2")]=0 #NA
predM[, c("UpperElev_m_2")]=0 #pmm
predM[, c("DisectMount_2")]=0 #logreg
predM[, c("BiogeogRealm_2")]=0 #NA
predM[, c("HabitatBreadth_2")]=0 #polyreg
predM[, c("Wingspan_cm_3")]=0 #pmm
predM[, c("AspectRatio_3")]=0 #pmm
predM[, c("WingLoad_Nm2_3")]=0 #pmm
predM[, c("PredCryp")]=0

#check collinearity
mice:::find.collinear(data.nodup, threshold = 0.99)

#do MICE imputation (with 15 iterations)
impute_30<-mice(data.nodup, method=meth, predictorMatrix=predM, maxit=30, print=T, remove.collinear=F)  

#check density plots for imputed traits
##density of the imputed data is showed in magenta while the density of the observed data is showed in blue
##the distributions should be similar
densityplot(impute_30)

#examine convergence of iterations for specific imputed variables
plot(impute_30, "PredCryp")


#join iterations and write imputed dataset to a csv file 
impute_30_complete<-complete(impute_30)
write.table(impute_30_complete, file="/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/impute_30_withcryptics.csv", col.names=T, row.names=F, sep=",")













