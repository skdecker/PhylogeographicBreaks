##Figure Making 
#Load Libraries
library("tools")
library("raster")
library("rgeos")
library("adegenet")
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("dplyr")
library("apex")
library("maptools")
library("ggmap")
library("stringr")
library("gridExtra")
library("svglite")
library("NatParksPalettes")

#------------------------------------------------------------------------------
#Figure 1
#Read in data with break coordinates from Monmonier's algorithm for birds
birdcoords1 = as.data.frame(read.csv('/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/NEWAvesbreaks.csv', header=FALSE, sep=","))
birdcoords2 = as.data.frame(read.csv('/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/NEWAvesbreaks2.csv', header=FALSE, sep=","))

#Read in trait data 
birdtraits = read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/Aves_Traits.csv", header=T)
birdspeciestomap = birdtraits$Species1[which(birdtraits$Break == "Y")] #Save species determined to have a break

#format species names in break dataset
for (i in 1:nrow(birdcoords1)){
  x=birdcoords1$V1[i]
  genus = str_split(x, "-", simplify=T)[,1]
  species = str_split(x, "-", simplify=T)[,2]
  birdcoords1$V4[i] = paste(genus, species, sep=" ")
}

#format species names in break dataset
for (i in 1:nrow(birdcoords2)){
  x=birdcoords2$V1[i]
  genus = str_split(x, "-", simplify=T)[,1]
  species = str_split(x, "-", simplify=T)[,2]
  birdcoords2$V4[i] = paste(genus, species, sep=" ")
}

#Read in data with break coordinates from Monmonier's algorithm for bats
batcoords1 = as.data.frame(read.csv('/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/NEWChirbreaks.csv', header=FALSE, sep=","))
batcoords2 = as.data.frame(read.csv('/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/NEWChirbreaks2.csv', header=FALSE, sep=","))

#Read in trait data
battraits = read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/impute_30_withcrypticsIBD.csv", header=T)
batspeciestomap = battraits$Species[which(battraits$Break == "Y")] #Save species determined to have a break

#Format species names in break dataset
for (i in 1:nrow(batcoords1)){
  x=batcoords1$V1[i]
  genus = str_split(x, "-", simplify=T)[,1]
  species = str_split(x, "-", simplify=T)[,2]
  batcoords1$V4[i] = paste(genus, species, sep=" ")
}

#Format species names in break dataset
for (i in 1:nrow(batcoords2)){
  x=batcoords2$V1[i]
  genus = str_split(x, "-", simplify=T)[,1]
  species = str_split(x, "-", simplify=T)[,2]
  batcoords2$V4[i] = paste(genus, species, sep=" ")
}

#plot breaks as lines on a world map
data("wrld_simpl")
plot(wrld_simpl, axes=TRUE, col="grey65", border=NA) 
for (x in birdspeciestomap){
  index = (which(x == birdcoords1$V4))
  lines(birdcoords1[index,c(3)], birdcoords1[index,c(2)], col="darkseagreen1", lwd=1.8)
} 
for (x in birdspeciestomap){
  index = (which(x == birdcoords2$V4))
  lines(birdcoords2[index,c(3)], birdcoords2[index,c(2)], col="darkseagreen1", lwd=1.8)
} 
for (x in batspeciestomap){
  index = (which(x == batcoords1$V4))
  lines(batcoords1[index,c(3)], batcoords1[index,c(2)], col="steelblue3", lwd=1.8)
}
for (x in batspeciestomap){
  index = (which(x == batcoords2$V4))
  lines(batcoords2[index,c(3)], batcoords2[index,c(2)], col="steelblue3", lwd=1.8)
}  

#------------------------------------------------------------------------------
#Figure 2 Variable Importance Plots
##bats
batvars = read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/1Chiroptera_VarImp_50iter.csv", header=TRUE)
batvars$Category = as.factor(batvars$Category) #save trait category as factor
batpos = batvars[which(batvars$avg.mda >= 0),] #save only traits with mda > 0

Cvarplot = ggplot(batpos, aes(x=reorder(LabelName, avg.mda), avg.mda, fill=Category, label=LabelName)) + 
  geom_col(color="black", linewidth=0.4) + 
  scale_fill_manual(values=c("Taxonomic"="#3F3F7B", "Geographic"="#278192", "Morphology"="#00B089", "Ecology"="#8FF7BD", "Environmental"="#34364F")) +
  geom_text(aes(y=avg.mda), angle=0, hjust=1.05, size=4)+
  coord_flip() +
  scale_y_reverse(limits=c(21,0)) +
  labs(x="", y="Mean Decrease in Accuracy") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill="none")

#save as svg
svglite(filename = "/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/Chir.VarImp1.svg", width=6, height=7)
Cvarplot
dev.off()

##birds
birdvars = read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/2Aves_VarImp_50iterWAvg.csv", header=TRUE) 
bird_pos = birdvars[which(birdvars$Avg.MDA > 0),] #save traits with mda > 0 
bird_pos$Category = as.factor(bird_pos$Category)  #save trait categories as factor

Avarplot = ggplot(bird_pos, aes(x=reorder(VarName, Avg.MDA), Avg.MDA, fill=Category, label=VarName)) + 
  geom_col(color="black", linewidth=0.4) + 
  scale_fill_manual(values=c("Taxonomy"="#3F3F7B", "Geographic"="#278192", "Morphology"="#00B089", "Ecology"="#8FF7BD", "Environmental"="#34364F")) +
  geom_text(aes(y=Avg.MDA), angle=0, hjust=1.05, size=4)+
  coord_flip() +
  scale_y_reverse(limits=c(16,0)) +
  labs(x="", y="Mean Decrease in Accuracy") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(fill="none")

#save as svg
svglite(filename = "/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/ModelResults/Aves.VarImp.svg", width=6, height=7)
Avarplot
dev.off()


#-----------------------------------------------------------------------------
#Figure 2 Box Plots
##Bats
#read in data
all<-read.csv("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/impute_30_withcrypticsIBD.csv", header=T)

##Area Sampled
areaplot = ggplot(all, aes(x=Break, y=Area)) + geom_boxplot(fill=c("#278192", "#0bced9")) + coord_flip() + labs(x="", y="Occurrence Area")+ scale_x_discrete(label = c("No Break", "Break"))+ theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##lengthlat
lenlatplot = ggplot(all, aes(x=Break, y=Length_Lat)) + geom_boxplot(fill=c("#278192", "#0bced9"))+ coord_flip() +labs(x="", y="Latitude Length")+ scale_x_discrete(label = c("No Break", "Break"))+ theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##maxLat
maxlatplot = ggplot(all, aes(x=Break, y=Max_Lat)) + geom_boxplot(fill=c("#278192", "#0bced9"))+ coord_flip() + labs(x="", y="Maximum Latitude")+ scale_x_discrete(label = c("No Break", "Break"))+ theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##brain Mass
brainmassplot = ggplot(all, aes(x=Break, y=AdultBrainMass_g_2)) + geom_boxplot(fill=c("#00B089", "#00d7a7"))+ coord_flip() + labs(x="", y="Brain Mass")+ scale_x_discrete(label = c("No Break", "Break"))+ theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##wingspan
wingspanplot = ggplot(all, aes(x=Break, y=Wingspan_cm_3)) + geom_boxplot(fill=c("#00B089", "#00d7a7"))+ coord_flip() + labs(x="", y="Wingspan") + scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

#plot all significant variables
setwd("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/Figures/")
svglite("4batsigvar.svg", width=5, height=8)
grid.arrange(areaplot, lenlatplot, maxlatplot, brainmassplot, wingspanplot, nrow=5)
dev.off()

##Birds
#read in data
birds = read.csv(file="/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Data/Aves_Traits.csv", header=TRUE)
birds$Break = as.factor(birds$Break)

##AbsmaxLat
absmaxlatplot = ggplot(birds, aes(x=Break, y=AbsMax_Lat)) + geom_boxplot(fill=c("#278192", "#0bced9")) + coord_flip() + labs(x="", y="Absolute Maximum Latitude")+ scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##MinLat
minlatplot = ggplot(birds, aes(x=Break, y=Min.Latitude)) + geom_boxplot(fill=c("#278192", "#0bced9")) + coord_flip() + labs(x="", y="Minimum Latitude")

##Beak_nares
beak_naresplot = ggplot(birds, aes(x=Break, y=Beak.Length_Nares)) + geom_boxplot(fill=c("#00B089", "#00d7a7")) + coord_flip() + labs(x="", y="Beak Length from Nares")+ scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##Beak_culmen
beak_culplot= ggplot(birds, aes(x=Break, y=Beak.Length_Culmen)) + geom_boxplot(fill=c("#00B089", "#00d7a7")) + coord_flip() + labs(x="", y="Beak Length from Culmen")+ scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##winglength
winglenplot= ggplot(birds, aes(x=Break, y=Wing.Length)) + geom_boxplot(fill=c("#00B089", "#00d7a7")) + coord_flip() + labs(x="", y="Wing Length")+ scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))

##tarsus length
tarsusplot= ggplot(birds, aes(x=Break, y=Tarsus.Length)) + geom_boxplot(fill=c("#00B089", "#00d7a7")) + coord_flip() + labs(x="", y="Tarsus Length")+ scale_x_discrete(label = c("No Break", "Break")) + theme(axis.text=element_text(size=12)) + theme(axis.title=element_text(size=14))


#plot all significant variables
setwd("/Users/decker.391/Documents/Research/Phylogatr/FinalVersions/Results/Figures/")
svglite("1birdsigvarboxes.svg", width=5, height=8)
grid.arrange(absmaxlatplot, beak_naresplot, beak_culplot, winglenplot, tarsusplot, nrow=5)
dev.off()
