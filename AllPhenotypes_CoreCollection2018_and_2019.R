###############################################################################################
#                 Script to do Bayz to analyse data faba bean core collection 2018 and 2019
###############################################################################################

# Load libraries
library(dplyr)
library(tidyr)
library(geosphere)
library(BayzR)
library(lme4)
library(ggplot2)
require(gridExtra)
library("lubridate")
library(rlist)
library("modules")
library(geosphere)

setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Phenotypic_analysis/20200311_AllPhenotypes_CoreCollection2018_and_2019")


##### functions

# Function to load data
function_to_load_data=function(filename){
  data=read.table(filename,sep="\t")
  colnames(data)=c("PLID","GPS","TRID","PlotSize","PlotSizewithpaths","GP_Origin","Genotype","PDID","PD_Description","Score",
                   "scorer","date","Comments")
  data[data == "\\N"] <- NA
  data[,13]=NULL
  return(as.data.frame(data))
}

# function to label trials NS or Se and add Year
makevariables=function(df){

  df=separate(data=df,col=GPS,into=c("lon","lat"),sep=";")
  df=separate(data=df,col=PlotSizewithpaths,into=c("length","width"),sep="mx")
  df=separate(data=df,col=width,into=c("width","leftover"),sep="m")
  df$leftover=NULL

  df$Location=rep(NA,nrow(df))
  df$Year=rep(NA,nrow(df))

  for (i in seq(1:nrow(df))){
    Nameofperson=as.character(df$scorer[i])

    if (startsWith(Nameofperson,"Jens")==T | startsWith(Nameofperson,"Andrea")==T | startsWith(Nameofperson,"Cathrine")==T){
      df$Location[i]="NS"
    }

    if (startsWith(Nameofperson,"Linda")==T | startsWith(Nameofperson,"Winnie")==T){
      df$Location[i]="Se"
    }

    date=df$date[i]
    df$Year[i]=substring(date,7,10)

    if (is.na(date)){
      df$Year[which(df$TRID==11)]='2018'
      df$Year[which(df$TRID==13)]='2018'
      df$Year[which(df$TRID==22)]='2019'
      df$Year[which(df$TRID==23)]='2019'
    }
  }

  return(df)

}


# function to predict neighbours of each plot
Neighboureffects=function(df){
  
  for (i in 1:nrow(df)){
    plotno=df$PLID[i]
    trno=df$TRID[i]
    Subset=df[which(df$TRID==trno),]
    
    # Take a plot and calculate its distance to all other plots in the same trial
    Subset_idx=which(Subset$PLID==plotno)
    lon1=as.numeric(as.character(Subset$lon[Subset_idx]))
    lat1=as.numeric(as.character(Subset$lat[Subset_idx]))
    if (length(lon1)>1 & length(unique(lon1)==1)){
      lon1=lon1[1] #if the same plot is registered twice 
      lat1=lon1[1]
    }
    NEI=cbind(rep(NA,nrow(Subset)),rep(NA,nrow(Subset)))
    colnames(NEI)=c("Neighbour","Distance")
    NEI=as.data.frame(NEI)
    for (j in seq(1:nrow(Subset))){
      neighbour=Subset$PLID[j]
      lon2=as.numeric(as.character(Subset$lon[j]))
      lat2=as.numeric(as.character(Subset$lat[j]))
      distance=distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
      NEI$Neighbour[j]=neighbour
      NEI$Distance[j]=distance
    }
    # Choose closest neighbours. As many as defined by parameters
    ClosestToMostDistantNei=order(NEI$Distance)
    PlotLength=as.numeric(Subset$length[1])
    PlotWidth=as.numeric(Subset$width[1])
    LargestDim=max(c(PlotLength,PlotWidth))
    
    CloseEnough=which(round(NEI$Distance,2)<=LargestDim & NEI$Distance!=0)
    CloseNeighbours=NEI[CloseEnough,]
    
    # find phenotypes of those neighbours
    NeiIdx=which(Subset$PLID %in% CloseNeighbours$Neighbour)
    NeiScores=df$Score[NeiIdx]
    
    df$AverageNeiScores[i]=mean(as.numeric(as.character(NeiScores)),na.rm=T)
  }
  return(df)
}

NeighboureffectsMoreDates=function(df){
  
  for (i in 1:nrow(df)){
    plotno=df$PLID[i]
    trno=df$TRID[i]
    date=df$date[i]
    Subset=df[which(df$TRID==trno & df$date==date),]
    
    # Take a plot and calculate its distance to all other plots in the same trial
    Subset_idx=which(Subset$PLID==plotno)
    lon1=as.numeric(as.character(Subset$lon[Subset_idx]))
    lat1=as.numeric(as.character(Subset$lat[Subset_idx]))
    NEI=cbind(rep(NA,nrow(Subset)),rep(NA,nrow(Subset)))
    colnames(NEI)=c("Neighbour","Distance")
    NEI=as.data.frame(NEI)
    for (j in seq(1:nrow(Subset))){
      neighbour=Subset$PLID[j]
      lon2=as.numeric(as.character(Subset$lon[j]))
      lat2=as.numeric(as.character(Subset$lat[j]))
      distance=distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
      NEI$Neighbour[j]=neighbour
      NEI$Distance[j]=distance
    }
    # Choose closest neighbours. As many as defined by parameters
    ClosestToMostDistantNei=order(NEI$Distance)
    PlotLength=as.numeric(Subset$length[1])
    PlotWidth=as.numeric(Subset$width[1])
    LargestDim=max(c(PlotLength,PlotWidth))
    
    CloseEnough=which(round(NEI$Distance,2)<=LargestDim & NEI$Distance!=0)
    CloseNeighbours=NEI[CloseEnough,]
    
    # find phenotypes of those neighbours
    NeiIdx=which(Subset$PLID %in% CloseNeighbours$Neighbour)
    NeiScores=df$Score[NeiIdx]
    
    df$AverageNeiScores[i]=mean(as.numeric(as.character(NeiScores)),na.rm=T)
  }
  return(df)
}



################## Height###################

# Load file
Heightdf=function_to_load_data("PlantHeightCore.txt")
dim(Heightdf)
head(Heightdf)
Heightdf=as.data.frame(Heightdf)

# Make variables ready
Heightdf_ready=makevariables(Heightdf)

Heightdf_ready=Neighboureffects(Heightdf_ready)
Heightdf_ready$AverageNeiScores

# if average neighbour effect is NA then set to average
NAforNEI=which(is.na(Heightdf_ready$AverageNeiScores))
PlotWithNoNeighbour=Heightdf_ready$PLID[NAforNEI]
TRID_here=Heightdf_ready$TRID[NAforNEI]
idx_in_this_tr=which(Heightdf_ready$TRID==TRID_here)
meanintr=mean(Heightdf_ready$AverageNeiScores[idx_in_this_tr],na.rm=T)
Heightdf_ready$AverageNeiScores[NAforNEI]=as.numeric(meanintr)


# Take a look at the data
Heightdf_ready$YearLoc=paste(Heightdf_ready$Location,Heightdf_ready$Year,sep="")
Trials=unique(Heightdf_ready$YearLoc)
length(Trials)

lines=unique(Heightdf_ready$GP_Origin)
length(lines) #124 lines in collection
counting_lines=count(Heightdf_ready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>9),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

Heightdf_filtered=Heightdf_ready[-which(Heightdf_ready$GP_Origin=="Kontu"),]
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="Lynx"),]
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="Taifun"),]
nrow(Heightdf_filtered)

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<9),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="x1"),]
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="x3"),]
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="x5"),]
Heightdf_filtered=Heightdf_filtered[-which(Heightdf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Heightdf_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 10 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
Heightdf_filtered$Score=as.numeric(as.character(Heightdf_filtered$Score))
nrow(Heightdf_filtered)
Heightdf_filtered=Heightdf_filtered[-which(is.na(Heightdf_filtered$Score)==T),]
nrow(Heightdf_filtered)



#histograms of observations from all trials
All=ggplot(data=Heightdf_filtered, aes(x=(Heightdf_filtered$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) +
  labs(x="Height (cm)", y="Observations") +
  geom_vline(xintercept = mean(Heightdf_filtered$Score,na.rm=T), col="green")

pdf("PlantHeight_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=Heightdf_filtered[which(Heightdf_filtered$TRID==11),]
nrow(NS18)

NS19=Heightdf_filtered[which(Heightdf_filtered$TRID==22),]
nrow(NS19)

Se19=Heightdf_filtered[which(Heightdf_filtered$TRID==23),]
nrow(Se19)


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + xlim(10,155) + ylim(0,60) +
  labs(x="Height (cm)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green", size=2) #looks normally distributed, mean 65 cm

p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + xlim(10,155) + ylim(0,60) +
  labs(x="Height (cm)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm

p3=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + xlim(10,155) + ylim(0,60) +
  labs(x="Height (cm)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) #looks close to normally distributed, mean 71 cm

pdf("PlantHeight.pdf", 15, 30)
grid.arrange(p1,p2,p3)
dev.off()


Heightdf_filtered$Score=as.numeric(as.character(Heightdf_filtered$Score))
Heightdf_filtered$lon=as.numeric(as.character(Heightdf_filtered$lon))
Heightdf_filtered$lat=as.numeric(as.character(Heightdf_filtered$lat))
Heightdf_filtered=droplevels(Heightdf_filtered)


# Designate block effect?
#NS18$block=NS18$lat
#dfnew=data.frame(x=as.numeric(NS18$lon), y=as.numeric(NS18$lat), row.names = seq(1:length(NS18$lat)))
#plot(dfnew, pch = 16, col = "red")
#library(gatepoints)
#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS18$block[as.numeric(selectedPoints)]=2

#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS18$block[as.numeric(selectedPoints)]=3

#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS18$block[as.numeric(selectedPoints)]=1

#NS19$block=NS19$lat
#dfnew=data.frame(x=as.numeric(NS19$lon), y=as.numeric(NS19$lat), row.names = seq(1:length(NS19$lat)))
#plot(dfnew, pch = 16, col = "red")
#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS19$block[as.numeric(selectedPoints)]=2

#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS19$block[as.numeric(selectedPoints)]=3

#selectedPoints <- fhs(dfnew, mark = TRUE)
#NS19$block[as.numeric(selectedPoints)]=1


# Fit phenotypical models in lme4
fit <- bayz((Score) ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) + freg(AverageNeiScores),
            data = Heightdf_filtered, chain=c(100000, 500, 10))



summary(fit) #line explains 75% when all together.there is a genotype-yearloc (gxe) effect of 8.0%
plot(fit)
BLUPS=fit$Estimates[7:123,1]
ind_names=rownames(fit$Estimates)[7:123]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Year%2018"),]
fixf2=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Year%2019"),]
fixf3=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%NS"),]
fixf4=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%Se"),]
fixf5=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit$Estimates[which(rownames(fit$Estimates)=="freg.AverageNeiScores"),]

fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8

fit$Parameters #to check sizes of parameters

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Heightdf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Heightdf_filtered$GP_Origin %in% colnames(GRM))==F)
Heightdf_filtered_again=Heightdf_filtered[-phenotypeswithnogeno,]

length(unique(Heightdf_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded <- bayz(Score ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores) +
              ranf(droplevels(GP_Origin), V=GRM_filtered),data = Heightdf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded)
fit_expanded$Parameters #to check sizes of parameters


rownames(fit_expanded$Estimates)
fixf1=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Year%2018"),]
fixf2=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Year%2019"),]
fixf3=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Location%NS"),]
fixf4=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Location%Se"),]
fixf5=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="freg.AverageNeiScores"),]
fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Heightdf_filtered$GP_Origin),as.numeric(as.character(Heightdf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_Height_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_Height_0518.csv",quote=F,row.names = F, col.names = F, sep=",")






################## Downy Mildew susceptibility ###################

# Load file
Downy=function_to_load_data("DownyCore.txt")
dim(Downy)
head(Downy)
Downydf=as.data.frame(Downy)

# Make variables ready
Downydf_ready=makevariables(Downydf)
Downydf_ready=Neighboureffects(Downydf_ready)

# if average neighbour effect is NA then set to average
NAforNEI=which(is.na(Downydf_ready$AverageNeiScores))
PlotWithNoNeighbour=Downydf_ready$PLID[NAforNEI]
TRID_here=Downydf_ready$TRID
idx_in_this_tr=which(Downydf_ready$TRID==TRID_here)
meanintr=mean(Downydf_ready$AverageNeiScores[idx_in_this_tr],na.rm=T)
Downydf_ready$AverageNeiScores[NAforNEI]=as.numeric(meanintr)

# Take a look at the data
Downydf_ready$YearLoc=paste(Downydf_ready$Location,Downydf_ready$Year,sep="")
Trials=unique(Downydf_ready$YearLoc)
length(Trials)

lines=unique(Downydf_ready$GP_Origin)
length(lines) #124 lines in collection
counting_lines=count(Downydf_ready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Trials))),]
MoreObservations

Donydf_filtered=Downydf_ready[-which(Downydf_ready$GP_Origin=="Kontu"),]
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="Lynx"),]
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="Taifun"),]
nrow(Donydf_filtered)

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Trials)),]
LessObservations
# Not core lines should be filtered out
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="x1"),]
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="x3"),]
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="x5"),]
Donydf_filtered=Donydf_filtered[-which(Donydf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Donydf_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Donydf_filtered$Score=as.numeric(as.character(Donydf_filtered$Score))
Donydf_filtered$Score=as.numeric(as.character(Donydf_filtered$Score))
nrow(Donydf_filtered)
errors=which(Donydf_filtered$Score<0 | Donydf_filtered$Score>9)
Donydf_filtered=Donydf_filtered[-errors,]
nrow(Donydf_filtered)
Donydf_filtered=Donydf_filtered[-which(is.na(Donydf_filtered$Score)==T),]

#histograms of observations from all trials
All=ggplot(data=Donydf_filtered, aes(x=(as.numeric(as.character(Donydf_filtered$Score))))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Susceptibility (character)", y="Observations") +
  geom_vline(xintercept = mean(as.numeric(as.character(Donydf_filtered$Score)),na.rm=T), col="green") +
  scale_x_continuous(breaks = seq(0, 9, 1))


pdf("DownyMildewSuscp_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS19=Donydf_filtered[which(Donydf_filtered$TRID==22),]
nrow(NS19)

Se19=Donydf_filtered[which(Donydf_filtered$TRID==23),]
nrow(Se19)


# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,100) +
  labs(x="Susceptibility (character)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,100) +
  labs(x="Susceptibility (character)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("DownyMildewSusceptibility.pdf", 15, 30)
grid.arrange(p2,p3)
dev.off()



# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each
Donydf_filtered=droplevels(Donydf_filtered)
# Fit phenotypical models in lme4
fit <- bayz((Score) ~ ranf(GP_Origin) + fixf(Location) + ran2f(GP_Origin,Location) +freg(AverageNeiScores),
            data = Donydf_filtered, chain=c(100000, 500, 10))

summary(fit) #line explains 75% when all together.there is a genotype-yearloc (gxe) effect of 8.0%

BLUPS=fit$Estimates[5:116,1]
ind_names=rownames(fit$Estimates)[5:116]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


fixf3=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%NS"),]
fixf4=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%Se"),]
fixf8=fit$Estimates[which(rownames(fit$Estimates)=="freg.AverageNeiScores"),]

fixf3
fixf4
fixf8

fit$Parameters #to check sizes of parameters

# Extend model to include GRM
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Donydf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Donydf_filtered$GP_Origin %in% colnames(GRM))==F)
Donydf_filtered_again=Donydf_filtered[-phenotypeswithnogeno,]

length(unique(Donydf_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

Donydf_filtered_again=droplevels(Donydf_filtered_again)

fit_expanded <- bayz(Score ~ ranf(GP_Origin) + fixf(Location) + ran2f(GP_Origin,Location) +freg(AverageNeiScores) +
                       ranf(droplevels(GP_Origin), V=GRM_filtered),data = Donydf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded)
fit_expanded$Parameters #to check sizes of parameters


rownames(fit_expanded$Estimates)
fixf4=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Location%Se"),]
fixf4
fixf6=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="freg.AverageNeiScores"),]
fixf6




# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Donydf_filtered$GP_Origin),as.numeric(as.character(Donydf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
phenotypefileheight_notcor=na.omit(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_DownyMildew_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_DownyMildew_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################## Earliness ###################

# Load file
Ear=function_to_load_data("EarlinessCore.txt")
dim(Ear)
head(Ear)
Eardf=as.data.frame(Ear)

# Make variables ready
Eardf_ready=makevariables(Eardf)
Eardf_ready=Neighboureffects(Eardf_ready)

# Take a look at the data
Eardf_ready$YearLoc=paste(Eardf_ready$Location,Eardf_ready$Year,sep="")
Trials=unique(Eardf_ready$YearLoc)
length(Trials)

lines=unique(Eardf_ready$GP_Origin)
length(lines) #128 lines in collection
counting_lines=count(Eardf_ready, vars =GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Trials))),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens

Eardf_filtered=Eardf_ready[-which(Eardf_ready$GP_Origin=="Kontu"),]
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="Lynx"),]
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="Taifun"),]
nrow(Eardf_filtered) #1207

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Trials)),]
LessObservations
# Core10, Core126, Core130, Core14, Core148, Core16, Core2, Core23, Core291, Core292, Core305, Core306, Core307, Core313, Core314, Core315, Core316 and Core319 only in 2019?
# Core52 only one obs. Sejet2018. Fits with data.
# Core5, 11 obs. Missing one from Sejet2018. Fits with data. One could suscept that Core52 is in fact core5, but we will not make any assumptions
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="Core52"),]


# Not core lines should be filtered out
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="x1"),]
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="x3"),]
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="x5"),]
Eardf_filtered=Eardf_filtered[-which(Eardf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Eardf_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Eardf_filtered$Score=as.numeric(as.character(Eardf_filtered$Score))
nrow(Eardf_filtered)
errors=which(Eardf_filtered$Score<0)
nodata=which(is.na(Eardf_filtered$Score)==T)
Eardf_filtered=Eardf_filtered[-nodata,]
nrow(Eardf_filtered) #1182 observations

#histograms of observations from all trials
All=ggplot(data=Eardf_filtered, aes(x=(as.integer(Eardf_filtered$Score)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Earliness of Flowering (Date in May)", y="Observations") +
  geom_vline(xintercept = mean(Eardf_filtered$Score,na.rm=T), col="green") +
  xlim(0,max(Eardf_filtered$Score)) +
  scale_x_continuous(breaks = seq(0, max(Eardf_filtered$Score), 4))

pdf("Earliness_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials

NS18=Eardf_filtered[which(Eardf_filtered$TRID==11),]
nrow(NS18)

Se18=Eardf_filtered[which(Eardf_filtered$TRID==13),]
nrow(Se18)

Se19=Eardf_filtered[which(Eardf_filtered$TRID==23),]
nrow(Se19)

NS19=Eardf_filtered[which(Eardf_filtered$TRID==22),]
nrow(NS19)



# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Eardf_filtered$Score)) + ylim(0,90) +
  labs(x="Earliness of Flowering (Date in May)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=Se18, aes(x=(Se18$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Eardf_filtered$Score)) + ylim(0,90) +
  labs(x="Earliness of Flowering (Date in May)", y="Observations") +
  ggtitle("Sejet 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se18$Score,na.rm=T), col="green",size=2)

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Eardf_filtered$Score)) + ylim(0,90) +
  labs(x="Earliness of Flowering (Date in May)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2)

p4=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Eardf_filtered$Score)) + ylim(0,90) +
  labs(x="Earliness of Flowering (Date in May)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2)

pdf("EarlinessOfFlowering.pdf", 15, 60)
grid.arrange(p1,p2,p3,p4,nrow=4)
dev.off()


Eardf_filtered$Score=as.numeric(as.character(Eardf_filtered$Score))
Eardf_filtered$lon=as.numeric(as.character(Eardf_filtered$lon))
Eardf_filtered$lat=as.numeric(as.character(Eardf_filtered$lat))
Eardf_filtered=droplevels(Eardf_filtered)
nrow(Eardf_filtered)
# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each

fit_earliness <- bayz((Score) ~ ranf(GP_Origin) + fixf(Location) +  fixf(Year) + fixf(YearLoc) + ran2f(GP_Origin,Location) +ran2f(GP_Origin,Year)+ran2f(GP_Origin,YearLoc)+freg(AverageNeiScores),
                  data = Eardf_filtered, chain=c(100000, 500, 10))

summary(fit_earliness) #line explains 58% when all together.there is a genotype-yearloc (gxe) effect of 19.2%


fixf3_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fix12_=fit_earliness$Estimates[which(rownames(fit_earliness$Estimates)=="freg.AverageNeiScores"),]
fix12_




BLUPS=fit_earliness$Estimates[7:123,1]
ind_names=rownames(fit_earliness$Estimates)[7:123]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_earliness$Parameters

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Eardf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Eardf_filtered$GP_Origin %in% colnames(GRM))==F)
Eardf_filtered_again=Eardf_filtered[-phenotypeswithnogeno,]

length(unique(Eardf_filtered_again$GP_Origin))==dim(GRM_filtered) # 99 individuals in both
Eardf_filtered_again=droplevels(Eardf_filtered_again)


fit_expanded_earliness <-bayz(Score ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores) +
                                ranf(droplevels(GP_Origin), V=GRM_filtered),data = Eardf_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_earliness)
fit_expanded_earliness$Parameters

fixf3_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fix12_=fit_expanded_earliness$Estimates[which(rownames(fit_expanded_earliness$Estimates)=="freg.AverageNeiScores"),]
fix12_



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Eardf_filtered$GP_Origin),as.numeric(as.character(Eardf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_Earliness_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_Earliness_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



################## End of flowering ###################

# Load file
End=function_to_load_data("Endflow.txt")
dim(End)
head(End)
Enddf=as.data.frame(End)

# Make variables ready
Enddf_ready=makevariables(Enddf)
Enddf_ready=Neighboureffects(Enddf_ready)


# Take a look at the data
Enddf_ready$YearLoc=paste(Enddf_ready$Location,Enddf_ready$Year,sep="")
Trials=unique(Enddf_ready$YearLoc)
length(Trials)

lines=unique(Enddf_ready$GP_Origin)
length(lines) #128 lines in collection
counting_lines=count(Enddf_ready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Trials))),]
MoreObservations
Enddf_filtered=Enddf_ready[-which(Enddf_ready$GP_Origin=="Kontu"),]
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="Lynx"),]
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="Taifun"),]
nrow(Enddf_filtered) #610

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Trials)),]
LessObservations
# Core10, Core126, Core130, Core14, Core148, Core16, Core2, Core23, Core291, Core292, Core305, Core306, Core307, Core313, Core314, Core315, Core316 and Core319 only in 2019?
# Core52 only one obs. Sejet2018. Fits with data.
# Core5, 11 obs. Missing one from Sejet2018. Fits with data. One could suscept that Core52 is in fact core5, but we will not make any assumptions


# Not core lines should be filtered out
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="x1"),]
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="x3"),]
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="x5"),]
Enddf_filtered=Enddf_filtered[-which(Enddf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Enddf_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens

###############Quality check of data##############
Enddf_filtered$Score=as.numeric(as.character(Enddf_filtered$Score))
nrow(Enddf_filtered) #598
errors=which(Enddf_filtered$Score<0)
nodata=which(is.na(Enddf_filtered$Score)==T)
Enddf_filtered=Enddf_filtered[-nodata,]
nrow(Enddf_filtered) #589 observations

#histograms of observations from all trials
All=ggplot(data=Enddf_filtered, aes(x=(as.integer(Enddf_filtered$Score)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="End of Flowering (Date in June)", y="Observations") +
  geom_vline(xintercept = mean(Enddf_filtered$Score,na.rm=T), col="green") +
  xlim(0,max(Enddf_filtered$Score)) +
  scale_x_continuous(breaks = seq(0, max(Enddf_filtered$Score), 4))

pdf("End_all.pdf", 15, 15)
grid.arrange(All)
dev.off()


#histograms of observations from each trial

Trials
NS19=Enddf_filtered[which(Enddf_filtered$TRID==22),]
nrow(NS19)

NS18=Enddf_filtered[which(Enddf_filtered$TRID==11),]
nrow(NS18)


# Make histogram plots
#histograms of observations from all trials, save these

p1=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Enddf_filtered$Score)) + ylim(0,80) +
  labs(x="End of Flowering (Date in June)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Enddf_filtered$Score)) + ylim(0,80) +
  labs(x="End of Flowering (Date in June)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2)

pdf("EndOfFlowering.pdf", 15, 30)
grid.arrange(p2,p1,nrow=2)
dev.off()


Enddf_filtered=droplevels(Enddf_filtered)

fit_endflow <- bayz((Score) ~ ranf(GP_Origin) +  fixf(Year)  +ran2f(GP_Origin,Year) +freg(AverageNeiScores),
                    data = Enddf_filtered, chain=c(100000, 500, 10))

summary(fit_endflow) #line explains 66%



BLUPS=fit_endflow$Estimates[5:121,1]
ind_names=rownames(fit_endflow$Estimates)[5:121]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_endflow$Parameters


fix1=fit_endflow$Estimates[which(rownames(fit_endflow$Estimates)=="fixf.Year%2019"),]
fix1
fix2=fit_endflow$Estimates[which(rownames(fit_endflow$Estimates)=="freg.AverageNeiScores"),]
fix2


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Enddf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Enddf_filtered$GP_Origin %in% colnames(GRM))==F)
Enddf_filtered_again=Enddf_filtered[-phenotypeswithnogeno,]

length(unique(Enddf_filtered_again$GP_Origin.x))==dim(GRM_filtered) # 98 individuals in both
Enddf_filtered_again=droplevels(Enddf_filtered_again)


fit_expanded_endflow <-bayz((Score) ~ ranf(GP_Origin) + fixf(Year)  +ran2f(GP_Origin,Year) +freg(AverageNeiScores) +
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = Enddf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_endflow)

fit_expanded_endflow$Parameters
fit_expanded_endflow$Estimates

fix1=fit_expanded_endflow$Estimates[which(rownames(fit_expanded_endflow$Estimates)=="fixf.Year%2019"),]
fix1
fix2=fit_expanded_endflow$Estimates[which(rownames(fit_expanded_endflow$Estimates)=="freg.AverageNeiScores"),]
fix2

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Enddf_filtered$GP_Origin),as.numeric(as.character(Enddf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_EndFlow_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_EndFlow_0518.csv",quote=F,row.names = F, col.names = F, sep=",")





################## Duration of flowering ###################

nrow(Enddf_filtered)
nrow(Eardf_filtered)

LenofFlow_df=merge(Enddf_filtered,Eardf_filtered,by="PLID")
nrow(LenofFlow_df)

LenofFlow_df$Score.x #end of flowering date in june
LenofFlow_df$Score.y # earliness of flowering date in may
LenofFlow_df$earlinessinjune=LenofFlow_df$Score.y-31

LenofFlow_df$ScoreLen=LenofFlow_df$Score.x-LenofFlow_df$earlinessinjune
LenofFlow_df=LenofFlow_df[-which(LenofFlow_df$ScoreLen==49),] #Sort out outlier

nrow(LenofFlow_df)
#histograms of observations from all trials
All=ggplot(data=LenofFlow_df, aes(x=(as.integer(LenofFlow_df$ScoreLen)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="End of Flowering (Date in June)", y="Observations") +
  geom_vline(xintercept = mean(LenofFlow_df$ScoreLen,na.rm=T), col="green") +
  xlim(0,max(LenofFlow_df$ScoreLen)) +
  scale_x_continuous(breaks = seq(0, max(LenofFlow_df$ScoreLen), 4))

pdf("Dur_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

Trials=unique(LenofFlow_df$YearLoc.x)
Trials

NS18=LenofFlow_df[which(LenofFlow_df$TRID.x==11),]
nrow(NS18)

NS19=LenofFlow_df[which(LenofFlow_df$TRID.x==22),]
nrow(NS19)




p1=ggplot(data=NS18, aes(x=(NS18$ScoreLen))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(NS18$ScoreLen)) + ylim(0,65) +
  labs(x="Duration of Flowering (Days)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$ScoreLen,na.rm=T), col="green",size=2)

p2=ggplot(data=NS19, aes(x=(NS19$ScoreLen))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(NS19$ScoreLen)) + ylim(0,80) +
  labs(x="Duration of Flowering (Days)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$ScoreLen,na.rm=T), col="green",size=2)


pdf("DurationOfFlowering_yeara.pdf", 15, 30)
grid.arrange(p1,p2)
dev.off()


Neighboureffects1=function(df){
  
  for (i in 1:nrow(df)){
    plotno=df$PLID[i]
    trno=df$TRID.x[i]
    Subset=df[which(df$TRID.x==trno),]
    
    #Take a plot and calculate its distance to all other plots in the same trial
    Subset_idx=which(Subset$PLID==plotno)
    lon1=as.numeric(as.character(Subset$lon.x[Subset_idx]))
    lat1=as.numeric(as.character(Subset$lat.x[Subset_idx]))
    NEI=cbind(rep(NA,nrow(Subset)),rep(NA,nrow(Subset)))
    colnames(NEI)=c("Neighbour","Distance")
    NEI=as.data.frame(NEI)
    for (j in seq(1:nrow(Subset))){
      neighbour=Subset$PLID[j]
      lon2=as.numeric(as.character(Subset$lon.x[j]))
      lat2=as.numeric(as.character(Subset$lat.x[j]))
      distance=distm(c(lon1, lat1), c(lon2, lat2), fun = distHaversine)
      NEI$Neighbour[j]=neighbour
      NEI$Distance[j]=distance
    }
    #Choose closest neighbours. As many as defined by parameters
    ClosestToMostDistantNei=order(NEI$Distance)
    PlotLength=as.numeric(Subset$length.x[1])
    PlotWidth=as.numeric(Subset$width.x[1])
    LargestDim=max(c(PlotLength,PlotWidth))
    
    CloseEnough=which(round(NEI$Distance,2)<=LargestDim & NEI$Distance!=0)
    CloseNeighbours=NEI[CloseEnough,]
    
    # find phenotypes of those neighbours
    NeiIdx=which(Subset$PLID %in% CloseNeighbours$Neighbour)
    NeiScores=df$Score.x[NeiIdx]
    
    df$AverageNeiScores[i]=mean(as.numeric(as.character(NeiScores)),na.rm=T)
  }
  return(df)
}

LenofFlow_df=Neighboureffects1(LenofFlow_df)

# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each
LenofFlow_df$ScoreLen=as.numeric(as.character(LenofFlow_df$ScoreLen))
LenofFlow_df$lon.x=as.numeric(as.character(LenofFlow_df$lon.x))
LenofFlow_df$lat.x=as.numeric(as.character(LenofFlow_df$lat.x))
LenofFlow_df=droplevels(LenofFlow_df)
nrow(LenofFlow_df)


fit_durflow <- bayz((ScoreLen) ~ ranf(GP_Origin.x) +  fixf(Year.x)  +ran2f(GP_Origin.x,Year.x) +freg(AverageNeiScores),
                      data = LenofFlow_df, chain=c(100000, 500, 10))

summary(fit_durflow) #line explains 66%


fix1=fit_durflow$Estimates[which(rownames(fit_durflow$Estimates)=="fixf.Year.x%2019"),]
fix1
fix2=fit_durflow$Estimates[which(rownames(fit_durflow$Estimates)=="freg.AverageNeiScores"),]
fix2

BLUPS=fit_durflow$Estimates[5:121,1]
ind_names=rownames(fit_durflow$Estimates)[5:121]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_durflow$Parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% LenofFlow_df$GP_Origin.x)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((LenofFlow_df$GP_Origin.x %in% colnames(GRM))==F)
LenofFlow_df_again=LenofFlow_df[-phenotypeswithnogeno,]

length(unique(LenofFlow_df_again$GP_Origin.x))==dim(GRM_filtered) # 98 individuals in both
LenofFlow_df_again=droplevels(LenofFlow_df_again)


fit_expanded_durflow <-bayz((ScoreLen) ~ ranf(GP_Origin.x) + fixf(Year.x)  +ran2f(GP_Origin.x,Year.x) +freg(AverageNeiScores) +
                          ranf(droplevels(GP_Origin.x), V=GRM_filtered),
                        data = LenofFlow_df_again, chain=c(200000, 1000, 100))

summary(fit_expanded_durflow)

fit_expanded_durflow$Parameters
fit_expanded_durflow$Estimates

fix1=fit_expanded_durflow$Estimates[which(rownames(fit_expanded_durflow$Estimates)=="fixf.Year.x%2019"),]
fix1
fix2=fit_expanded_durflow$Estimates[which(rownames(fit_expanded_durflow$Estimates)=="freg.AverageNeiScores"),]
fix2

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(LenofFlow_df$GP_Origin.x),as.numeric(as.character(LenofFlow_df$ScoreLen)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_FlowDur_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_FlowDur_0518.csv",quote=F,row.names = F, col.names = F, sep=",")


################## botrytis perc. ###################

# Load file
Botr_perc=function_to_load_data("BotrPercCore.txt")
dim(Botr_perc)
head(Botr_perc)
Botr_percdf=as.data.frame(Botr_perc)

# Make variables ready
Botr_percdf_ready=makevariables(Botr_percdf)
Botr_percdf_ready=NeighboureffectsMoreDates(Botr_percdf_ready)
Botr_percdf_ready$AverageNeiScores

# if average neighbour effect is NA then set to average
NAforNEI=which(is.na(Botr_percdf_ready$AverageNeiScores))
PlotWithNoNeighbour=Botr_percdf_ready$PLID[NAforNEI]
TRID_here=Botr_percdf_ready$TRID[NAforNEI]
idx_in_this_tr=which(Botr_percdf_ready$TRID==TRID_here)
meanintr=mean(Botr_percdf_ready$AverageNeiScores[idx_in_this_tr],na.rm=T)
Botr_percdf_ready$AverageNeiScores[NAforNEI]=as.numeric(meanintr)



# Take a look at the data
Botr_percdf_ready$YearLoc=paste(Botr_percdf_ready$Location,Botr_percdf_ready$Year,sep="")
Trials=unique(Botr_percdf_ready$YearLoc)
length(Trials) #only data from NS2019 but we have more dates

Dates=unique(Botr_percdf_ready$date)
length(Dates) #3 dates


lines=unique(Botr_percdf_ready$GP_Origin)
length(lines) #119 lines in collection
counting_lines=count(Botr_percdf_ready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Dates))),]
MoreObservations
Botr_Perc_filtered=Botr_percdf_ready[-which(Botr_percdf_ready$GP_Origin=="Kontu"),]
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="Lynx"),]
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="Taifun"),]
nrow(Botr_Perc_filtered) #1044

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Dates)),]
LessObservations
# Not core lines should be filtered out
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="x1"),]
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="x3"),]
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="x5"),]
Botr_Perc_filtered=Botr_Perc_filtered[-which(Botr_Perc_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Botr_Perc_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Botr_Perc_filtered$Score=as.numeric(as.character(Botr_Perc_filtered$Score))
nrow(Botr_Perc_filtered)
errors=which(Botr_Perc_filtered$Score<0 | Botr_Perc_filtered$Score>100)
nodata=which(is.na(Botr_Perc_filtered$Score)==T)
Botr_Perc_filtered=Botr_Perc_filtered[-nodata,]
nrow(Botr_Perc_filtered) #952 observations

# Remove observations where the affected leaf coverage gets smaller as the days progress. Then we assume it was probably not botrytis
unique_plotids=unique(Botr_Perc_filtered$PLID)
rm_error=list()

for (plid in unique_plotids){
  row_idx=which(Botr_Perc_filtered$PLID==plid)
  Scores=Botr_Perc_filtered$Score[row_idx]
  day=Botr_Perc_filtered$date[row_idx]
  dates=parse_date_time(day, orders = c("ymd", "dmy", "mdy"))
  minitable=cbind(row_idx,Scores,dates)
  minitable_sorted=minitable[order(dates),]

  if (length(day)>1){
    if (all(diff(minitable_sorted[,2])>=0)==F) {
      rm_error=list.append(rm_error, plid)
    }
  }
}
indexestoremove=which(Botr_Perc_filtered$PLID %in% rm_error)
Botr_Perc_filtered=Botr_Perc_filtered[-indexestoremove,]
nrow(Botr_Perc_filtered)

#histograms of observations from all trials
All=ggplot(data=Botr_Perc_filtered, aes(x=(as.integer(Botr_Perc_filtered$Score)))) +
  geom_histogram(binwidth=3,fill = "black",
                 alpha=.8) +
  labs(x="Botrytis susceptibility (% leaf coverage)", y="Observations") +
  geom_vline(xintercept = mean(Botr_Perc_filtered$Score,na.rm=T), col="green") +
  xlim(0,max(Botr_Perc_filtered$Score)) +
  scale_x_continuous(breaks = seq(0, max(Botr_Perc_filtered$Score), 4))

pdf("BotrytisPerc_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Dates

firstdate=Botr_Perc_filtered[which(Botr_Perc_filtered$date=="25/06/2019"),]
nrow(firstdate)

seconddate=Botr_Perc_filtered[which(Botr_Perc_filtered$date=="08/07/2019"),]
nrow(seconddate)

thirddate=Botr_Perc_filtered[which(Botr_Perc_filtered$date=="17/07/2019"),]
nrow(thirddate)


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=firstdate, aes(x=(firstdate$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Botr_Perc_filtered$Score)) + ylim(0,100) +
  labs(x="Botrytis susceptibility (% leaf coverage)", y="Observations") +
  ggtitle("Nordic Seed 2019, 25/06/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(firstdate$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=seconddate, aes(x=(seconddate$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Botr_Perc_filtered$Score)) + ylim(0,100) +
  labs(x="Botrytis susceptibility (% leaf coverage)", y="Observations") +
  ggtitle("Nordic Seed 2019, 08/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(seconddate$Score,na.rm=T), col="green",size=2)

p3=ggplot(data=thirddate, aes(x=(thirddate$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Botr_Perc_filtered$Score)) + ylim(0,100) +
  labs(x="Botrytis susceptibility (% leaf coverage)", y="Observations") +
  ggtitle("Nordic Seed 2019, 17/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(thirddate$Score,na.rm=T), col="green",size=2)

pdf("Botrytis_perc.pdf", 15, 30)
grid.arrange(p1,p2,p3,nrow=3)
dev.off()


Botr_Perc_filtered$Score=as.numeric(as.character(Botr_Perc_filtered$Score))
Botr_Perc_filtered$lon=as.numeric(as.character(Botr_Perc_filtered$lon))
Botr_Perc_filtered$lat=as.numeric(as.character(Botr_Perc_filtered$lat))
Botr_Perc_filtered=droplevels(Botr_Perc_filtered)

# Fit phenotypical models in bayz
fit_botrPerc <- bayz((Score) ~ ranf(GP_Origin) +fixf(date) +freg(AverageNeiScores),
            data = Botr_Perc_filtered, chain=c(100000, 500, 10))

summary(fit_botrPerc) #line explains 57%

plot(fit_botrPerc)
BLUPS=fit_botrPerc$Estimates[4:111,1]
ind_names=rownames(fit_botrPerc$Estimates)[4:111]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit_botrPerc$Estimates[which(rownames(fit_botrPerc$Estimates)=="fixf.date%17/07/2019"),]
fixf2=fit_botrPerc$Estimates[which(rownames(fit_botrPerc$Estimates)=="fixf.date%25/06/2019"),]
fixf3=fit_botrPerc$Estimates[which(rownames(fit_botrPerc$Estimates)=="freg.AverageNeiScores"),]

fixf1
fixf2
fixf3

fit_botrPerc$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Botr_Perc_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Botr_Perc_filtered$GP_Origin %in% colnames(GRM))==F)
Botr_Perc_filtered_again=Botr_Perc_filtered[-phenotypeswithnogeno,]

length(unique(Botr_Perc_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded_botrPerc <-bayz((Score) ~ ranf(droplevels(GP_Origin)) + fixf(date) + +freg(AverageNeiScores) +
                                ranf(droplevels(GP_Origin), V=GRM_filtered),
                              data = Botr_Perc_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_botrPerc)
fit_expanded_botrPerc$Parameters #to check sizes of parameters


fixf1=fit_expanded_botrPerc$Estimates[which(rownames(fit_expanded_botrPerc$Estimates)=="fixf.date%17/07/2019"),]
fixf2=fit_expanded_botrPerc$Estimates[which(rownames(fit_expanded_botrPerc$Estimates)=="fixf.date%25/06/2019"),]
fixf3=fit_botrPerc$Estimates[which(rownames(fit_botrPerc$Estimates)=="freg.AverageNeiScores"),]


fixf1
fixf2
fixf3

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Botr_Perc_filtered$GP_Origin),as.numeric(as.character(Botr_Perc_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrPerc_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrPerc_0518.csv",quote=F,row.names = F, col.names = F, sep=",")


### here 



##### Botrytis early date, 25/06/2019
firstdate=droplevels(firstdate)
fit_botrPerc_firstdate <- bayz((Score) ~ ranf(GP_Origin) + +freg(AverageNeiScores),
                     data = firstdate, chain=c(100000, 500, 10))

summary(fit_botrPerc_firstdate) #line explains 70%
fit_botrPerc_firstdate$Parameters
fit_botrPerc_firstdate$Estimates[which(rownames(fit_botrPerc_firstdate$Estimates)=="freg.AverageNeiScores"),]


BLUPS=fit_botrPerc_firstdate$Estimates[4:111,1]
ind_names=rownames(fit_botrPerc_firstdate$Estimates)[4:111]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% firstdate$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((firstdate$GP_Origin %in% colnames(GRM))==F)
Botr_Perc_filtered_again=firstdate[-phenotypeswithnogeno,]

length(unique(Botr_Perc_filtered_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

Botr_Perc_filtered_again=droplevels(Botr_Perc_filtered_again)
fit_expanded_botrPerc_firstdate <-bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores) +
                               ranf(droplevels(GP_Origin), V=GRM_filtered),
                             data = Botr_Perc_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_botrPerc_firstdate)
fit_expanded_botrPerc_firstdate$Parameters
fit_expanded_botrPerc_firstdate$Estimates[which(rownames(fit_expanded_botrPerc_firstdate$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(firstdate$GP_Origin),as.numeric(as.character(firstdate$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_botrPerc2506_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_botrPerc2506_0518.csv",quote=F,row.names = F, col.names = F, sep=",")









##### Botrytis early date, 08/07/2019
seconddate=droplevels(seconddate)
fit_botrPerc_seconddate <- bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores),
                               data = seconddate, chain=c(100000, 500, 10))

summary(fit_botrPerc_seconddate) #line explains 74%
fit_botrPerc_seconddate$Parameters
fit_botrPerc_seconddate$Estimates[which(rownames(fit_botrPerc_seconddate$Estimates)=="freg.AverageNeiScores"),]


BLUPS=fit_botrPerc_seconddate$Estimates[4:110,1]
ind_names=rownames(fit_botrPerc_seconddate$Estimates)[4:110]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seconddate$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seconddate$GP_Origin %in% colnames(GRM))==F)
Botr_Perc_filtered_again=seconddate[-phenotypeswithnogeno,]

length(unique(Botr_Perc_filtered_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

Botr_Perc_filtered_again=droplevels(Botr_Perc_filtered_again)
fit_expanded_botrPerc_seconddate <-bayz((Score) ~ ranf(GP_Origin) + freg(AverageNeiScores) +
                                         ranf(droplevels(GP_Origin), V=GRM_filtered),
                                       data = Botr_Perc_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_botrPerc_seconddate)
fit_expanded_botrPerc_seconddate$Parameters
fit_expanded_botrPerc_seconddate$Estimates[which(rownames(fit_expanded_botrPerc_seconddate$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seconddate$GP_Origin),as.numeric(as.character(seconddate$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_botrPerc0807_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_botrPerc0807_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




##### Botrytis early date, 17/07/2019
thirddate=droplevels(thirddate)
fit_botrPerc_thirddate <- bayz((Score) ~ ranf(GP_Origin) + freg(AverageNeiScores),
                                data = thirddate, chain=c(100000, 500, 10))

summary(fit_botrPerc_thirddate) #line explains 78%
fit_botrPerc_thirddate$Estimates[which(rownames(fit_botrPerc_thirddate$Estimates)=="freg.AverageNeiScores"),]

fit_botrPerc_thirddate$Parameters

BLUPS=fit_botrPerc_thirddate$Estimates[4:109,1]
ind_names=rownames(fit_botrPerc_thirddate$Estimates)[4:109]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% thirddate$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((thirddate$GP_Origin %in% colnames(GRM))==F)
Botr_Perc_filtered_again=thirddate[-phenotypeswithnogeno,]

length(unique(Botr_Perc_filtered_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

Botr_Perc_filtered_again=droplevels(Botr_Perc_filtered_again)
fit_expanded_botrPerc_thirddate <-bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores) +
                                          ranf(droplevels(GP_Origin), V=GRM_filtered),
                                        data = Botr_Perc_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_botrPerc_thirddate)
fit_expanded_botrPerc_thirddate$Estimates[which(rownames(fit_expanded_botrPerc_thirddate$Estimates)=="freg.AverageNeiScores"),]

fit_expanded_botrPerc_thirddate$Parameters

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(thirddate$GP_Origin),as.numeric(as.character(thirddate$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_botrPerc1707_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_botrPerc1707_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



######################################################################
#     Botrytis character upper leaves only. Latest date              #
######################################################################
# Load file
Botr_perc=function_to_load_data("BotrPercUpperSecCore.txt")
dim(Botr_perc)
head(Botr_perc)
Botr_percdf=as.data.frame(Botr_perc)


# Make variables ready
Botr_percdf_upperready=makevariables(Botr_percdf)
Botr_percdf_upperready=Neighboureffects(Botr_percdf_upperready)
Botr_percdf_upperready$AverageNeiScores

# Take a look at the data
Botr_percdf_upperready$YearLoc=paste(Botr_percdf_upperready$Location,Botr_percdf_upperready$Year,sep="")
Trials=unique(Botr_percdf_upperready$YearLoc)
length(Trials) #only data from NS2019

Dates=unique(Botr_percdf_upperready$date)
length(Dates) #only 1 date


lines=unique(Botr_percdf_upperready$GP_Origin)
length(lines) #119 lines in collection
counting_lines=count(Botr_percdf_upperready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Dates))),]
MoreObservations
Botr_Perc_filtered_upper=Botr_percdf_upperready[-which(Botr_percdf_upperready$GP_Origin=="Kontu"),]
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="Lynx"),]
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="Taifun"),]
nrow(Botr_Perc_filtered_upper) #348

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Dates)),]
LessObservations
# Not core lines should be filtered out
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="x1"),]
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="x3"),]
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="x5"),]
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-which(Botr_Perc_filtered_upper$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Botr_Perc_filtered_upper, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Botr_Perc_filtered_upper$Score=as.numeric(as.character(Botr_Perc_filtered_upper$Score))
nrow(Botr_Perc_filtered_upper)
errors=which(Botr_Perc_filtered_upper$Score<0 | Botr_Perc_filtered_upper$Score>100)
nodata=which(is.na(Botr_Perc_filtered_upper$Score)==T)
Botr_Perc_filtered_upper=Botr_Perc_filtered_upper[-nodata,]
nrow(Botr_Perc_filtered_upper) #952 observations


#histograms of observations from all trials
All= ggplot(data=Botr_Perc_filtered_upper, aes(x=(as.integer(Botr_Perc_filtered_upper$Score)))) +
  geom_histogram(binwidth=3,fill = "black",
                 alpha=.8) +
  labs(x="Botrytis susceptibility (% upper leaf coverage)", y="Observations") +
  geom_vline(xintercept = mean(Botr_Perc_filtered_upper$Score,na.rm=T), col="green") +
  xlim(0,max(Botr_Perc_filtered_upper$Score)) +
  scale_x_continuous(breaks = seq(0, max(Botr_Perc_filtered_upper$Score), 4))

pdf("BotrytisPerc_2307upperleaves.pdf", 15, 15)
grid.arrange(All)
dev.off()

# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each
Botr_Perc_filtered_upper$Score=as.numeric(as.character(Botr_Perc_filtered_upper$Score))
Botr_Perc_filtered_upper$lon=as.numeric(as.character(Botr_Perc_filtered_upper$lon))
Botr_Perc_filtered_upper$lat=as.numeric(as.character(Botr_Perc_filtered_upper$lat))
Botr_Perc_filtered_upper=droplevels(Botr_Perc_filtered_upper)

fit_botrPerc <- bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores),
                     data = Botr_Perc_filtered_upper, chain=c(100000, 500, 10))

fit_botrPerc$Estimates[which(rownames(fit_botrPerc$Estimates)=="freg.AverageNeiScores"),]

summary(fit_botrPerc) #line explains 57%

plot(fit_botrPerc)
BLUPS=fit_botrPerc$Estimates[4:104,1]
ind_names=rownames(fit_botrPerc$Estimates)[4:104]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_botrPerc$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Botr_Perc_filtered_upper$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Botr_Perc_filtered_upper$GP_Origin %in% colnames(GRM))==F)
Botr_Perc_filtered_upper_again=Botr_Perc_filtered_upper[-phenotypeswithnogeno,]

length(unique(Botr_Perc_filtered_upper_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded_botrPerc_upper <-bayz((Score) ~ ranf(droplevels(GP_Origin)) +freg(AverageNeiScores)+
                               ranf(droplevels(GP_Origin), V=GRM_filtered),
                             data = Botr_Perc_filtered_upper_again, chain=c(200000, 1000, 100))

summary(fit_expanded_botrPerc_upper)
fit_expanded_botrPerc_upper$Estimates[which(rownames(fit_expanded_botrPerc_upper$Estimates)=="freg.AverageNeiScores"),]

fit_expanded_botrPerc_upper$Parameters #to check sizes of parameters



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Botr_Perc_filtered_upper$GP_Origin),as.numeric(as.character(Botr_Perc_filtered_upper$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrPercUpperleavesOnly1707_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrPercUpperleavesOnly1707_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################## Botrytis character ###################

# Load file
BotrCha=function_to_load_data("BotrChaCore.txt")
dim(BotrCha)
head(BotrCha)
BotrChadf=as.data.frame(BotrCha)

# Make variables ready
BotrChadf_ready=makevariables(BotrChadf)
BotrChadf_ready=NeighboureffectsMoreDates(BotrChadf_ready)
BotrChadf_ready$AverageNeiScores

# Take a look at the data
BotrChadf_ready$YearLoc=paste(BotrChadf_ready$Location,BotrChadf_ready$Year,sep="")
Trials=unique(BotrChadf_ready$YearLoc)
length(Trials)

lines=unique(BotrChadf_ready$GP_Origin)
length(lines) #114 lines in collection
counting_lines=count(BotrChadf_ready, vars = GP_Origin)
Dates=unique(BotrChadf_ready$date)
length(Dates) #2 dates

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Dates))),]
MoreObservations
Botr_cha_filtered=BotrChadf_ready[-which(BotrChadf_ready$GP_Origin=="Kontu"),]
Botr_cha_filtered=Botr_cha_filtered[-which(Botr_cha_filtered$GP_Origin=="Lynx"),]
nrow(Botr_cha_filtered) #672

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Dates)),]
LessObservations

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Botr_cha_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Botr_cha_filtered$Score=as.numeric(as.character(Botr_cha_filtered$Score))
errors=which(Botr_cha_filtered$Score<0 | Botr_cha_filtered$Score>9)
nodata=which(is.na(Botr_cha_filtered$Score)==T)
Botr_cha_filtered=Botr_cha_filtered[-nodata,]
nrow(Botr_cha_filtered) #656 observations


# Remove observations where the affected leaf coverage gets smaller as the days progress. Then we assume it was probably not botrytis
unique_plotids=unique(Botr_cha_filtered$PLID)
rm_error=list()

for (plid in unique_plotids){
  row_idx=which(Botr_cha_filtered$PLID==plid)
  Scores=Botr_cha_filtered$Score[row_idx]
  day=Botr_cha_filtered$date[row_idx]
  dates=parse_date_time(day, orders = c("ymd", "dmy", "mdy"))
  minitable=cbind(row_idx,Scores,dates)
  minitable_sorted=minitable[order(dates),]

  if (length(day)>1){
    if (all(diff(minitable_sorted[,2])>=0)==F) {
      rm_error=list.append(rm_error, plid)
    }
  }
}
indexestoremove=which(Botr_cha_filtered$PLID %in% rm_error)
Botr_cha_filtered=Botr_cha_filtered[-indexestoremove,] #638 observations

#histograms of observations from all trials
All=ggplot(data=Botr_cha_filtered, aes(x=(as.integer(Botr_cha_filtered$Score)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  geom_vline(xintercept = mean(Botr_cha_filtered$Score,na.rm=T), col="green") +
  xlim(0,max(Botr_cha_filtered$Score)) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("BotrytisCha_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Dates

firstdate=Botr_cha_filtered[which(Botr_cha_filtered$date=="03/07/2019"),]
nrow(firstdate)

seconddate=Botr_cha_filtered[which(Botr_cha_filtered$date=="16/07/2019"),]
nrow(seconddate)


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=firstdate, aes(x=(firstdate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,150) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  ggtitle("Sejet 2019, 03/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(firstdate$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=seconddate, aes(x=(seconddate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,150) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  ggtitle("Sejet 2019, 16/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(seconddate$Score,na.rm=T), col="green",size=2)

pdf("Botrytis_character.pdf", 15, 30)
grid.arrange(p1,p2,nrow=2)
dev.off()


Botr_cha_filtered$Score=as.numeric(as.character(Botr_cha_filtered$Score))
Botr_cha_filtered$lon=as.numeric(as.character(Botr_cha_filtered$lon))
Botr_cha_filtered$lat=as.numeric(as.character(Botr_cha_filtered$lat))
Botr_cha_filtered=droplevels(Botr_cha_filtered)
nrow(Botr_cha_filtered)



# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each
# Fit phenotypical models in bayz
fit_botrCha <- bayz((Score) ~ ranf(GP_Origin) +fixf(date) +freg(AverageNeiScores),
                     data = Botr_cha_filtered, chain=c(100000, 500, 10))

summary(fit_botrCha) #line explains 57%

plot(fit_botrCha)
BLUPS=fit_botrCha$Estimates[4:113,1]
ind_names=rownames(fit_botrCha$Estimates)[4:113]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit_botrCha$Estimates[which(rownames(fit_botrCha$Estimates)=="fixf.date%03/07/2019"),]
fixf2=fit_botrCha$Estimates[which(rownames(fit_botrCha$Estimates)=="fixf.date%16/07/2019"),]
fixf3=fit_botrCha$Estimates[which(rownames(fit_botrCha$Estimates)=="freg.AverageNeiScores"),]

fixf1
fixf2
fixf3

fit_botrCha$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Botr_cha_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Botr_cha_filtered$GP_Origin %in% colnames(GRM))==F)
Botr_cha_filtered_again=Botr_cha_filtered[-phenotypeswithnogeno,]

length(unique(Botr_cha_filtered_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

fit_expanded_botrCha <-bayz((Score) ~ ranf(GP_Origin) + fixf(as.factor(date)) + freg(AverageNeiScores) +
                               ranf(droplevels(GP_Origin), V=GRM_filtered),
                             data = Botr_cha_filtered_again, chain=c(100000, 5000, 100))
summary(fit_expanded_botrCha)
fit_expanded_botrCha$Parameters #to check sizes of parameters


fixf1=fit_expanded_botrCha$Estimates[which(rownames(fit_expanded_botrCha$Estimates)=="fixf.as.factor(date)%03/07/2019"),]
fixf2=fit_expanded_botrCha$Estimates[which(rownames(fit_expanded_botrCha$Estimates)=="fixf.as.factor(date)%16/07/2019"),]
fixf3=fit_expanded_botrCha$Estimates[which(rownames(fit_expanded_botrCha$Estimates)=="freg.AverageNeiScores"),]

fixf1
fixf2
fixf3
# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Botr_cha_filtered$GP_Origin),as.numeric(as.character(Botr_cha_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrCha_all_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrCha_all_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



##### Botrytis early date, 03/07/2019
firstdate$Score=as.numeric(as.character(firstdate$Score))
firstdate$lon=as.numeric(as.character(firstdate$lon))
firstdate$lat=as.numeric(as.character(firstdate$lat))
firstdate=droplevels(firstdate)

fit_botrCha_firstdate <- bayz((Score) ~ ranf(GP_Origin)  +freg(AverageNeiScores),
                    data = firstdate, chain=c(100000, 500, 10))

summary(fit_botrCha_firstdate) #line explains 51%
fit_botrCha_firstdate$Estimates[which(rownames(fit_botrCha_firstdate$Estimates)=="freg.AverageNeiScores"),]


plot(fit_botrCha_firstdate)
BLUPS=fit_botrCha_firstdate$Estimates[4:113,1]
ind_names=rownames(fit_botrCha_firstdate$Estimates)[4:113]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_botrCha_firstdate$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% firstdate$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((firstdate$GP_Origin %in% colnames(GRM))==F)
Botr_cha_firstdate_again=firstdate[-phenotypeswithnogeno,]

length(unique(Botr_cha_firstdate_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

fit_expanded_botrCha_firstdate <-bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores) +
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = Botr_cha_firstdate_again, chain=c(100000, 5000, 100))
summary(fit_expanded_botrCha_firstdate)
fit_expanded_botrCha_firstdate$Parameters #to check sizes of parameters
fit_expanded_botrCha_firstdate$Estimates[which(rownames(fit_expanded_botrCha_firstdate$Estimates)=="freg.AverageNeiScores"),]



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(firstdate$GP_Origin),as.numeric(as.character(firstdate$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrCha_0307_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrCha_0307_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



##### Botrytis late date, 16/07/2019
seconddate$Score=as.numeric(as.character(seconddate$Score))
seconddate$lon=as.numeric(as.character(seconddate$lon))
seconddate$lat=as.numeric(as.character(seconddate$lat))
seconddate=droplevels(seconddate)

fit_botrCha_secdate <- bayz((Score) ~ ranf(GP_Origin)  +freg(AverageNeiScores),
                              data = seconddate, chain=c(100000, 500, 10))

summary(fit_botrCha_secdate) #line explains 53%
fit_botrCha_secdate$Estimates[which(rownames(fit_botrCha_secdate$Estimates)=="freg.AverageNeiScores"),]


plot(fit_botrCha_secdate)
BLUPS=fit_botrCha_secdate$Estimates[4:113,1]
ind_names=rownames(fit_botrCha_secdate$Estimates)[4:113]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_botrCha_secdate$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seconddate$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seconddate$GP_Origin %in% colnames(GRM))==F)
Botr_cha_secdate_again=seconddate[-phenotypeswithnogeno,]

length(unique(Botr_cha_secdate_again$GP_Origin))==dim(GRM_filtered) # 94 individuals in both

fit_expanded_botrCha_secdate <-bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores) +
                                        ranf(droplevels(GP_Origin), V=GRM_filtered),
                                      data = Botr_cha_secdate_again, chain=c(100000, 5000, 100))
summary(fit_expanded_botrCha_secdate)
fit_expanded_botrCha_secdate$Parameters #to check sizes of parameters
fit_expanded_botrCha_secdate$Estimates[which(rownames(fit_expanded_botrCha_secdate$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seconddate$GP_Origin),as.numeric(as.character(seconddate$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrCha_1607_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrCha_1607_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



###############################################
#     Botrytis perc converted to character
###############################################

# First load until Botr_Perc_filtered for the NS data and till Botr_Cha_filtered for the Sejet 2019 data
library(BotrytisConverter)
library(plyr)

Botr_Perc_filtered_1=Botr_Perc_filtered # run again if you need to
Botr_cha_filtered_1=Botr_cha_filtered # run again if you need to

for (i in seq(1:nrow(Botr_Perc_filtered_1))){
  convertedvalue=conversion_func_threshold(Botr_Perc_filtered_1$Score[i])
  Botr_Perc_filtered_1$ScoresConverted[i]=convertedvalue
}

Botr_Perc_filtered_1$Score=NULL
colnames(Botr_Perc_filtered_1)
Botr_Perc_filtered_1$Score=Botr_Perc_filtered_1$ScoresConverted
Botr_Perc_filtered_1$ScoresConverted=NULL
colnames(Botr_Perc_filtered_1)
colnames(Botr_cha_filtered_1)
AllBotrytisData=rbind.fill(Botr_Perc_filtered_1,Botr_cha_filtered_1)
nrow(AllBotrytisData)==nrow(Botr_Perc_filtered_1)+nrow(Botr_cha_filtered_1)

AllBotrytisData
AllBotrytisData=NeighboureffectsMoreDates(AllBotrytisData)

#histograms of observations from all trials
All=ggplot(data=AllBotrytisData, aes(x=(as.integer(AllBotrytisData$Score)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  geom_vline(xintercept = mean(AllBotrytisData$Score,na.rm=T), col="green") +
  xlim(0,max(AllBotrytisData$Score)) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("Botrytis_Cha_Converted_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Dates=unique(AllBotrytisData$date)
Dates

firstdate=AllBotrytisData[which(AllBotrytisData$date=="25/06/2019"),]
nrow(firstdate)

seconddate=AllBotrytisData[which(AllBotrytisData$date=="03/07/2019"),]
nrow(seconddate)

thirddate=AllBotrytisData[which(AllBotrytisData$date=="08/07/2019"),]
nrow(thirddate)

fourthdate=AllBotrytisData[which(AllBotrytisData$date=="16/07/2019"),]
nrow(fourthdate)

fifthdate=AllBotrytisData[which(AllBotrytisData$date=="17/07/2019"),]
nrow(fifthdate)


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=firstdate, aes(x=(firstdate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,200) +
  labs(x="Botrytis susceptibility (character, converted)", y="Observations") +
  ggtitle("NordicSeed 2019, 25/06/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(firstdate$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=seconddate, aes(x=(seconddate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,200) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  ggtitle("Sejet 2019, 03/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(seconddate$Score,na.rm=T), col="green",size=2)

p3=ggplot(data=thirddate, aes(x=(thirddate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,200) +
  labs(x="Botrytis susceptibility (character, converted)", y="Observations") +
  ggtitle("Nordic Seed 2019, 08/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(thirddate$Score,na.rm=T), col="green",size=2)

p4=ggplot(data=fourthdate, aes(x=(fourthdate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,200) +
  labs(x="Botrytis susceptibility (character)", y="Observations") +
  ggtitle("Sejet 2019, 16/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(fourthdate$Score,na.rm=T), col="green",size=2)

p5=ggplot(data=fifthdate, aes(x=(fifthdate$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,9) + ylim(0,200) +
  labs(x="Botrytis susceptibility (character, converted)", y="Observations") +
  ggtitle("Nordic Seed 2019, 17/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(fifthdate$Score,na.rm=T), col="green",size=2)

pdf("Botrytis_character.pdf", 15, 75)
grid.arrange(p1,p2,p3,p4,p5,nrow=5)
dev.off()



# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each
AllBotrytisData$Score=as.numeric(as.character(AllBotrytisData$Score))
AllBotrytisData$lon=as.numeric(as.character(AllBotrytisData$lon))
AllBotrytisData$lat=as.numeric(as.character(AllBotrytisData$lat))
AllBotrytisData=droplevels(AllBotrytisData)
nrow(AllBotrytisData)


fit_botr_ALL <- bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) + fixf(date) +fixf(Location) +freg(AverageNeiScores),
                    data = AllBotrytisData, chain=c(100000, 500, 10))

summary(fit_botr_ALL) #line explains 50%

fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.date%03/07/2019"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.date%08/07/2019"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.date%16/07/2019"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.date%17/07/2019"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.date%25/06/2019"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.Location%NS"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="fixf.Location%Se"),]
fit_botr_ALL$Estimates[which(rownames(fit_botr_ALL$Estimates)=="freg.AverageNeiScores"),]



BLUPS=fit_botr_ALL$Estimates[5:116,1]
ind_names=rownames(fit_botr_ALL$Estimates)[5:116]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% AllBotrytisData$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((AllBotrytisData$GP_Origin %in% colnames(GRM))==F)
AllBotrytisData_filteredagain=AllBotrytisData[-phenotypeswithnogeno,]

length(unique(AllBotrytisData_filteredagain$GP_Origin))==dim(GRM_filtered) # 94 individuals in both


AllBotrytisData_filteredagain=droplevels(AllBotrytisData_filteredagain)

fit_expanded_botrCha_conv <-bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) + fixf(date) +fixf(Location) +freg(AverageNeiScores) +
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = AllBotrytisData_filteredagain, chain=c(100000, 500, 100))
summary(fit_expanded_botrCha_conv)
fit_expanded_botrCha_conv$Parameters

fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.date%03/07/2019"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.date%08/07/2019"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.date%16/07/2019"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.date%17/07/2019"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.date%25/06/2019"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%NS"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%Se"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(AllBotrytisData$GP_Origin),as.numeric(as.character(AllBotrytisData$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrChaAllConverted_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrChaAllConverted_0518.csv",quote=F,row.names = F, col.names = F, sep=",")









################################
# Botrytis two first dates
################################

Dates=unique(AllBotrytisData$date)
Dates

earlydates_botr=rbind.fill(firstdate,seconddate)


earlydates_botr=droplevels(earlydates_botr)
nrow(earlydates_botr)


fit_botr_early <- bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) +fixf(Location) +freg(AverageNeiScores),
                     data = earlydates_botr, chain=c(100000, 500, 10))

summary(fit_botr_early) #line explains 50%

fit_botr_early$Estimates[which(rownames(fit_botr_early$Estimates)=="fixf.Location%NS"),]
fit_botr_early$Estimates[which(rownames(fit_botr_early$Estimates)=="fixf.Location%Se"),]
fit_botr_early$Estimates[which(rownames(fit_botr_early$Estimates)=="freg.AverageNeiScores"),]



BLUPS=fit_botr_early$Estimates[5:116,1]
ind_names=rownames(fit_botr_early$Estimates)[5:116]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_botr_early$Parameters
#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% earlydates_botr$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((earlydates_botr$GP_Origin %in% colnames(GRM))==F)
AllBotrytisData_filteredagain=earlydates_botr[-phenotypeswithnogeno,]

length(unique(AllBotrytisData_filteredagain$GP_Origin))==dim(GRM_filtered) # 94 individuals in both


AllBotrytisData_filteredagain=droplevels(AllBotrytisData_filteredagain)

fit_expanded_botrCha_conv <-bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) + fixf(date) +fixf(Location) +freg(AverageNeiScores) +
                                   ranf(droplevels(GP_Origin), V=GRM_filtered),
                                 data = AllBotrytisData_filteredagain, chain=c(100000, 500, 100))
summary(fit_expanded_botrCha_conv)
fit_expanded_botrCha_conv$Parameters

fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%NS"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%Se"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(earlydates_botr$GP_Origin),as.numeric(as.character(earlydates_botr$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrChaConverted_EarlyDates_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrChaAllConverted_EarlyDates_BLUPS_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################################
# Botrytis two latest dates
################################

Dates=unique(AllBotrytisData$date)
Dates

latedates_botr=rbind.fill(fourthdate,fifthdate)


latedates_botr=droplevels(latedates_botr)
nrow(latedates_botr)


fit_botr_late <- bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) +fixf(Location) +freg(AverageNeiScores),
                       data = latedates_botr, chain=c(100000, 500, 10))

summary(fit_botr_late) #line explains 50%

fit_botr_late$Estimates[which(rownames(fit_botr_late$Estimates)=="fixf.Location%NS"),]
fit_botr_late$Estimates[which(rownames(fit_botr_late$Estimates)=="fixf.Location%Se"),]
fit_botr_late$Estimates[which(rownames(fit_botr_late$Estimates)=="freg.AverageNeiScores"),]

BLUPS=fit_botr_late$Estimates[5:116,1]
ind_names=rownames(fit_botr_late$Estimates)[5:116]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_botr_late$Parameters
#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% latedates_botr$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((latedates_botr$GP_Origin %in% colnames(GRM))==F)
AllBotrytisData_filteredagain=latedates_botr[-phenotypeswithnogeno,]

length(unique(AllBotrytisData_filteredagain$GP_Origin))==dim(GRM_filtered) # 94 individuals in both


AllBotrytisData_filteredagain=droplevels(AllBotrytisData_filteredagain)

fit_expanded_botrCha_conv <-bayz((Score) ~ ranf(GP_Origin) + ran2f(GP_Origin,Location) + fixf(date) +fixf(Location) +freg(AverageNeiScores) +
                                   ranf(droplevels(GP_Origin), V=GRM_filtered),
                                 data = AllBotrytisData_filteredagain, chain=c(100000, 500, 100))

summary(fit_expanded_botrCha_conv)
fit_expanded_botrCha_conv$Parameters

fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%NS"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="fixf.Location%Se"),]
fit_expanded_botrCha_conv$Estimates[which(rownames(fit_expanded_botrCha_conv$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(latedates_botr$GP_Origin),as.numeric(as.character(latedates_botr$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_BotrChaConverted_LateDates_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_BotrChaAllConverted_LateDates_BLUPS_0518.csv",quote=F,row.names = F, col.names = F, sep=",")


################## Hilum color###################

# Load file
hilum=function_to_load_data("HilumColCore_updated.txt")
dim(hilum)
head(hilum)
hilumdf=as.data.frame(hilum)

# Make variables ready
hilumdf_ready=makevariables(hilumdf)
hilumdf_ready=Neighboureffects(hilumdf_ready)



# Take a look at the data
hilumdf_ready$YearLoc=paste(hilumdf_ready$Location,hilumdf_ready$Year,sep="")
Trials=unique(hilumdf_ready$YearLoc)
length(Trials)

lines=unique(hilumdf_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(hilumdf_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

Hilumdf_filtered=hilumdf_ready[-which(hilumdf_ready$GP_Origin=="Kontu"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Lynx"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Taifun"),]
nrow(Hilumdf_filtered) #1001

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations

#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="x1"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="x3"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="x5"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Hilumdf_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core113 four times in Se2019, as Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core119 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations
# In all trials there is a mix of 1 and 2s for this core line

#Core12 four times in Se2019, as Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core144 registered two a mix of 1 and 2
# In all trials there is a mix of 1 and 2s for this core line

#Core30 has 4 in 2018 NS

#Core305 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core342 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core347 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core36 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core407 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core46 four times in Se2019, s Linda registered one as a mix of 1 and 2. It then went into the database as two observations
# In all trials there is a mix of 1 and 2s for this core line

#Core68 four times in Se2019, as Linda registered one as a mix of 1 and 2. It then went into the database as two observations

#Core72 five times in Se2019, as Linda registered 2 plots as a mix of 1 and 2. It then went into the database as two observations

#Core78 four times in Se2019, as Linda registered one as a mix of 1 and 2. It then went into the database as two observations
#Keep them all as it seems like they really are a mix






###############Quality check of data##############
Hilumdf_filtered$Score=as.numeric(as.character(Hilumdf_filtered$Score))
nodata=which(is.na(Hilumdf_filtered$Score)==T)
Hilumdf_filtered=Hilumdf_filtered[-nodata,]
nrow(Hilumdf_filtered) #929


# Identify lines that show differences in this morphology

uniquelines=unique(Hilumdf_filtered$GP_Origin)
for (i in seq(1:length(uniquelines))){
  idxs=which(Hilumdf_filtered$GP_Origin==uniquelines[i])
  phenotypes=Hilumdf_filtered$Score[idxs]
  uniquephenotypes=length(unique(na.omit(phenotypes)))
  if (uniquephenotypes>1){
    print(uniquelines[i])
  }
}
# I think we should remove lines that are consistently different in all observations in one trial compared to ALL other observations.
# There might be a mix up for these
# these are core403, core402 and core 62
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Core403"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Core402"),]
Hilumdf_filtered=Hilumdf_filtered[-which(Hilumdf_filtered$GP_Origin=="Core62"),]


#histograms of observations from all trials
All=ggplot(data=Hilumdf_filtered, aes(x=(Hilumdf_filtered$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Hilum colour (1=light, 2=dark)", y="Observations") +
  geom_vline(xintercept = mean(Hilumdf_filtered$Score,na.rm=T), col="green")

pdf("Hilumcol_all.pdf", 15, 15)
grid.arrange(All)
dev.off()


#histograms of observations from each trial
Trials
NS18=Hilumdf_filtered[which(Hilumdf_filtered$TRID==11),]
nrow(NS18)

NS19=Hilumdf_filtered[which(Hilumdf_filtered$TRID==22),]
nrow(NS19)

Se19=Hilumdf_filtered[which(Hilumdf_filtered$TRID==23),]
nrow(Se19)

# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(1,2) + ylim(0,300) +
  labs(x="Hilum color (1: Light, 2: Dark)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(1,2) + ylim(0,300) +
  labs(x="Hilum color (1: Light, 2: Dark)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p4=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(1,2) + ylim(0,300) +
  labs(x="Hilum color (1: Light, 2: Dark)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))


pdf("HilumColour.pdf", 15, 45)
grid.arrange(p2,p3,p4)
dev.off()






Hilumdf_filtered$Score=as.numeric(as.character(Hilumdf_filtered$Score))
Hilumdf_filtered$lon=as.numeric(as.character(Hilumdf_filtered$lon))
Hilumdf_filtered$lat=as.numeric(as.character(Hilumdf_filtered$lat))
Hilumdf_filtered=droplevels(Hilumdf_filtered)
nrow(Hilumdf_filtered)

# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each

fit_hilumcol <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +fixf(YearLoc) +fixf(Location) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores),
                  data = Hilumdf_filtered, chain=c(100000, 500, 10))

summary(fit_hilumcol) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fixf3_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix11_=fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fit_hilumcol$Estimates[which(rownames(fit_hilumcol$Estimates)=="freg.AverageNeiScores"),]

fit_hilumcol$Parameters

BLUPS=fit_hilumcol$Estimates[7:120,1]
ind_names=rownames(fit_hilumcol$Estimates)[7:120]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Hilumdf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Hilumdf_filtered$GP_Origin %in% colnames(GRM))==F)
Hilumdf_filtered_again=Hilumdf_filtered[-phenotypeswithnogeno,]

length(unique(Hilumdf_filtered_again$GP_Origin))==dim(GRM_filtered) # 99 individuals in both
Hilumdf_filtered_again=droplevels(Hilumdf_filtered_again)


fit_expanded_hilum <- bayz(Score ~  fixf(Year) +fixf(Location) +fixf(YearLoc)  +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = Hilumdf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_hilum)
fit_expanded_hilum$Parameters

fixf3_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fit_expanded_hilum$Estimates[which(rownames(fit_expanded_hilum$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Hilumdf_filtered$GP_Origin),as.numeric(as.character(Hilumdf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_hilumcol_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_HilumCol_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################## Novitron damage ###################

# Load file
Novitron=function_to_load_data("NovitronDamageCore.txt")
dim(Novitron)
head(Novitron)
Novitrondf=as.data.frame(Novitron)

# Make variables ready
Novitrondf_ready=makevariables(Novitrondf)
Novitrondf_ready=Neighboureffects(Novitrondf_ready)



# Take a look at the data
Novitrondf_ready$YearLoc=paste(Novitrondf_ready$Location,Novitrondf_ready$Year,sep="")
Trials=unique(Novitrondf_ready$YearLoc)
length(Trials)

lines=unique(Novitrondf_ready$GP_Origin)
length(lines) #119 lines in collection
counting_lines=count(Novitrondf_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

Novitron_filtered=Novitrondf_ready[-which(Novitrondf_ready$GP_Origin=="Kontu"),]
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="Lynx"),]
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="Taifun"),]
nrow(Novitron_filtered) #348

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="x1"),]
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="x3"),]
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="x5"),]
Novitron_filtered=Novitron_filtered[-which(Novitron_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Novitron_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}



###############Quality check of data##############
Novitron_filtered$Score=as.numeric(as.character(Novitron_filtered$Score))
nrow(Novitron_filtered) #336
nodata=which(is.na(Novitron_filtered$Score)==T)
Novitron_filtered=Novitron_filtered[-nodata,]
nrow(Novitron_filtered) #1182 observations


#histograms of observations from each trial
Trials

NS19=Novitron_filtered[which(Novitron_filtered$TRID==22),]
nrow(NS19)

# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,5) + ylim(0,250) +
  labs(x="Novitron damage", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("NovitronDamage.pdf", 15, 15)
grid.arrange(p3)
dev.off()



Novitron_filtered$Score=as.numeric(as.character(Novitron_filtered$Score))
Novitron_filtered$lon=as.numeric(as.character(Novitron_filtered$lon))
Novitron_filtered$lat=as.numeric(as.character(Novitron_filtered$lat))
Novitron_filtered=droplevels(Novitron_filtered)
nrow(Novitron_filtered)
# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each

fit_NovDam <- bayz((Score) ~ ranf(GP_Origin) + freg(AverageNeiScores),
                      data = Novitron_filtered, chain=c(100000, 500, 10))

summary(fit_NovDam) #line explains 63% when all together.

fit_NovDam$Parameters

fit_NovDam$Estimates[which(rownames(fit_NovDam$Estimates)=="freg.AverageNeiScores"),]

BLUPS=fit_NovDam$Estimates[4:115,1]
ind_names=rownames(fit_NovDam$Estimates)[4:115]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Novitron_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Novitron_filtered$GP_Origin %in% colnames(GRM))==F)
Novitron_filtered_again=Novitron_filtered[-phenotypeswithnogeno,]

length(unique(Novitron_filtered_again$GP_Origin))==dim(GRM_filtered) # 99 individuals in both
Novitron_filtered_again=droplevels(Novitron_filtered_again)


fit_expanded_novitron <- bayz((as.factor(Score))~  ranf(GP_Origin) + freg(AverageNeiScores) +
                             ranf(droplevels(GP_Origin), V=GRM_filtered),
                           data = Novitron_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_novitron)
fit_expanded_novitron$Parameters
fit_expanded_novitron$Estimates[which(rownames(fit_expanded_novitron$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Novitron_filtered$GP_Origin),as.numeric(as.character(Novitron_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_Novitron_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_Novitron_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################## maturation date ###################

# Load file
Mat=function_to_load_data("MaturationCore.txt")
dim(Mat)
head(Mat)
Matdf=as.data.frame(Mat)

# Make variables ready
Matdf_ready=makevariables(Matdf)
Matdf_ready=Neighboureffects(Matdf_ready)



# Take a look at the data
Matdf_ready$YearLoc=paste(Matdf_ready$Location,Matdf_ready$Year,sep="")
Trials=unique(Matdf_ready$YearLoc)
length(Trials)

lines=unique(Matdf_ready$GP_Origin)
length(lines) #128 lines in collection
counting_lines=count(Matdf_ready, vars = GP_Origin)

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$n>(3*length(Trials))),]
MoreObservations
Matdf_filtered=Matdf_ready[-which(Matdf_ready$GP_Origin=="Kontu"),]
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="Lynx"),]
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="Taifun"),]
nrow(Matdf_filtered) #946

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<3*length(Trials)),]
LessObservations

# Not core lines should be filtered out
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="x1"),]
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="x3"),]
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="x5"),]
Matdf_filtered=Matdf_filtered[-which(Matdf_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Matdf_filtered, vars = GP_Origin)
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$n[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}


###############Quality check of data##############
Matdf_filtered$Score=as.numeric(as.character(Matdf_filtered$Score))
nrow(Matdf_filtered) #937
errors=which(Matdf_filtered$Score<0)
nodata=which(is.na(Matdf_filtered$Score)==T)
Matdf_filtered=Matdf_filtered[-nodata,]
nrow(Matdf_filtered) #907 observations

#histograms of observations from all trials
All= ggplot(data=Matdf_filtered, aes(x=(as.integer(Matdf_filtered$Score)))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="Maturation date (Date in July)", y="Observations") +
  geom_vline(xintercept = mean(Matdf_filtered$Score,na.rm=T), col="green") +
  xlim(0,max(Matdf_filtered$Score)) +
  scale_x_continuous(breaks = seq(0, max(Matdf_filtered$Score), 4))

pdf("Maturation_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials


NS19=Matdf_filtered[which(Matdf_filtered$TRID==22),]
nrow(NS19)

NS18=Matdf_filtered[which(Matdf_filtered$TRID==11),]
nrow(NS18)

Se19=Matdf_filtered[which(Matdf_filtered$TRID==23),]
nrow(Se19)



# Make histogram plots
#histograms of observations from all trials, save these

p1=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Matdf_filtered$Score)) + ylim(0,120) +
  labs(x="Maturation date (Date in July)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2)

p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Matdf_filtered$Score)) + ylim(0,120)+
  labs(x="Maturation date (Date in July)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2)

p3=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=2,fill = "black",
                 alpha=.8) + xlim(0,max(Matdf_filtered$Score)) + ylim(0,120) +
  labs(x="Maturation date (Date in July)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2)

pdf("Maturation.pdf", 15, 45)
grid.arrange(p1,p2,p3)
dev.off()

Matdf_filtered$Score=as.numeric(as.character(Matdf_filtered$Score))
Matdf_filtered$lon=as.numeric(as.character(Matdf_filtered$lon))
Matdf_filtered$lat=as.numeric(as.character(Matdf_filtered$lat))
Matdf_filtered=droplevels(Matdf_filtered)
nrow(Matdf_filtered)
# Fit phenotypical models in Bayz
# Ignore placement effect as it is a randomized blockdesign with 3 reps in each

fit_mat <- bayz((Score) ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores),
            data = Matdf_filtered, chain=c(100000, 500, 10))


summary(fit_mat) #line explains 69%


plot(fit_mat)
BLUPS=fit_mat$Estimates[7:123,1]
ind_names=rownames(fit_mat$Estimates)[7:123]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.Year%2018"),]
fixf2=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.Year%2019"),]
fixf3=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.Location%NS"),]
fixf4=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.Location%Se"),]
fixf5=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit_mat$Estimates[which(rownames(fit_mat$Estimates)=="freg.AverageNeiScores"),]
fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8

fit_mat$Parameters #to check sizes of parameters


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Matdf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Matdf_filtered$GP_Origin %in% colnames(GRM))==F)
Matdf_filtered_again=Matdf_filtered[-phenotypeswithnogeno,]

length(unique(Matdf_filtered_again$GP_Origin))==dim(GRM_filtered) # 99 individuals in both

Matdf_filtered_again=droplevels(Matdf_filtered_again)
nrow(Matdf_filtered_again)
fit_expanded_mat=bayz((Score) ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores) +
       + ranf(droplevels(GP_Origin), V=GRM_filtered),
     data = Matdf_filtered_again, chain=c(100000, 500, 10))


summary(fit_expanded_mat)
fit_expanded_mat$Parameters

fixf1=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.Year%2018"),]
fixf2=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.Year%2019"),]
fixf3=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.Location%NS"),]
fixf4=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.Location%Se"),]
fixf5=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit_expanded_mat$Estimates[which(rownames(fit_expanded_mat$Estimates)=="freg.AverageNeiScores"),]

fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8

fit_mat$Parameters #to check sizes of parameters

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Matdf_filtered$GP_Origin),as.numeric(as.character(Matdf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_Matu_NotCor_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_Matu_0518.csv",quote=F,row.names = F, col.names = F, sep=",")




################## TGW ###################

# Load file
TGW=function_to_load_data("TGW.txt")
dim(TGW)
head(TGW)
TGWdf=as.data.frame(TGW)

# Make variables ready
TGWf_ready=makevariables(TGWdf)
TGWf_ready=Neighboureffects(TGWf_ready)



# Take a look at the data
TGWf_ready$YearLoc=paste(TGWf_ready$Location,TGWf_ready$Year,sep="")
Trials=unique(TGWf_ready$YearLoc)
length(Trials)

lines=unique(TGWf_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(TGWf_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

TGW_filtered=TGWf_ready[-which(TGWf_ready$GP_Origin=="Kontu"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Lynx"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Taifun"),]
nrow(TGW_filtered) #946

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="x1"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="x3"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="x5"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(TGW_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
TGW_filtered$Score=as.numeric(as.character(TGW_filtered$Score))
nodata=which(is.na(TGW_filtered$Score)==T)
TGW_filtered=TGW_filtered[-nodata,]
nrow(TGW_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Core403"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Core402"),]
TGW_filtered=TGW_filtered[-which(TGW_filtered$GP_Origin=="Core62"),]
nrow(TGW_filtered) #885

#histograms of observations from all trials
All=ggplot(data=TGW_filtered, aes(x=(TGW_filtered$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) +
  labs(x="TGW (g)", y="Observations") +
  geom_vline(xintercept = mean(TGW_filtered$Score,na.rm=T), col="green")

pdf("TGW_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=TGW_filtered[which(TGW_filtered$TRID==11),]
nrow(NS18)

NS19=TGW_filtered[which(TGW_filtered$TRID==22),]
nrow(NS19)

Se19=TGW_filtered[which(TGW_filtered$TRID==23),]
nrow(Se19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=20,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,30) +
  labs(x="TGW (g)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=20,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,30) +
  labs(x="TGW (g)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p4=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=20,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,30) +
  labs(x="TGW (g)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("TGW.pdf", 15, 45)
grid.arrange(p2,p3,p4)
dev.off()


TGW_filtered$Score=as.numeric(as.character(TGW_filtered$Score))
TGW_filtered$lon=as.numeric(as.character(TGW_filtered$lon))
TGW_filtered$lat=as.numeric(as.character(TGW_filtered$lat))
TGW_filtered=droplevels(TGW_filtered)
nrow(TGW_filtered)


fit_TGW <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +fixf(YearLoc) +fixf(Location) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc)+freg(AverageNeiScores),
                     data = TGW_filtered, chain=c(100000, 500, 10))

summary(fit_TGW) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fixf3_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix11_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fix12_=fit_TGW$Estimates[which(rownames(fit_TGW$Estimates)=="freg.AverageNeiScores"),]
fix12_

fit_TGW$Parameters

BLUPS=fit_TGW$Estimates[7:120,1]
ind_names=rownames(fit_TGW$Estimates)[7:120]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% TGW_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((TGW_filtered$GP_Origin %in% colnames(GRM))==F)
TGW_filtered_again=TGW_filtered[-phenotypeswithnogeno,]
TGW_filtered_again=droplevels(TGW_filtered_again)

length(unique(TGW_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both


fit_expanded_TGW <- bayz(Score ~  fixf(Year) +fixf(Location) +fixf(YearLoc)  +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                             ranf(droplevels(GP_Origin), V=GRM_filtered),
                           data = TGW_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_TGW)
fit_expanded_TGW$Parameters

fixf3_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fit_expanded_TGW$Estimates[which(rownames(fit_expanded_TGW$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(TGW_filtered$GP_Origin),as.numeric(as.character(TGW_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_TGW_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_TGW_0518.csv",quote=F,row.names = F, col.names = F, sep=",")






################## Seed coat col ###################
# Load file
seedcoat=function_to_load_data("SeedCoatColCore.txt")
dim(seedcoat)
head(seedcoat)
seedcoatdf=as.data.frame(seedcoat)

# Make variables ready
seedcoat_ready=makevariables(seedcoatdf)
seedcoat_ready=Neighboureffects(seedcoat_ready)


# Take a look at the data
seedcoat_ready$YearLoc=paste(seedcoat_ready$Location,seedcoat_ready$Year,sep="")
Trials=unique(seedcoat_ready$YearLoc)
length(Trials)

lines=unique(seedcoat_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(seedcoat_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

seedcoat_filtered=seedcoat_ready[-which(seedcoat_ready$GP_Origin=="Kontu"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Lynx"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Taifun"),]
nrow(seedcoat_filtered) #960

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x1"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x3"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x5"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(seedcoat_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
seedcoat_filtered$Score=as.numeric(as.character(seedcoat_filtered$Score))
nodata=which(is.na(seedcoat_filtered$Score)==T)
seedcoat_filtered=seedcoat_filtered[-nodata,]
nrow(seedcoat_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core403"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core402"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core62"),]
nrow(seedcoat_filtered) #900

#histograms of observations from all trials
All=ggplot(data=seedcoat_filtered, aes(x=(seedcoat_filtered$Score))) +
  geom_histogram(binwidth=1, fill = "black",
                 alpha=.8) +
  labs(x="TGW (g)", y="Observations") +
  geom_vline(xintercept = mean(seedcoat_filtered$Score,na.rm=T), col="green")

pdf("SeedCoat_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=seedcoat_filtered[which(seedcoat_filtered$TRID==11),]
nrow(NS18)

NS19=seedcoat_filtered[which(seedcoat_filtered$TRID==22),]
nrow(NS19)

Se19=seedcoat_filtered[which(seedcoat_filtered$TRID==23),]
nrow(Se19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,300) +
  labs(x="Seed Coat Colour (1= White/gray, 2=Beige, 3=Green, 4=Violet, 5=Black)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,300) +
  labs(x="Seed Coat Colour (1= White/gray, 2=Beige, 3=Green, 4=Violet, 5=Black)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p4=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,300) +
  labs(x="Seed Coat Colour (1= White/gray, 2=Beige, 3=Green, 4=Violet, 5=Black)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("SeedCoatCol.pdf", 15, 45)
grid.arrange(p2,p3,p4)
dev.off()


seedcoat_filtered$Score=as.numeric(as.character(seedcoat_filtered$Score))
seedcoat_filtered$lon=as.numeric(as.character(seedcoat_filtered$lon))
seedcoat_filtered$lat=as.numeric(as.character(v$lat))
seedcoat_filtered=droplevels(seedcoat_filtered)
nrow(seedcoat_filtered)


fit_coatcol <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +fixf(YearLoc) +fixf(Location) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores),
                data = seedcoat_filtered, chain=c(100000, 500, 10))

summary(fit_coatcol) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fixf3_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix11_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_

fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="freg.AverageNeiScores"),]


fit_coatcol$Parameters

BLUPS=fit_coatcol$Estimates[7:120,1]
ind_names=rownames(fit_coatcol$Estimates)[7:120]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seedcoat_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seedcoat_filtered$GP_Origin %in% colnames(GRM))==F)
seedcoat_filtered_again=seedcoat_filtered[-phenotypeswithnogeno,]
seedcoat_filtered_again=droplevels(seedcoat_filtered_again)

length(unique(seedcoat_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both


fit_expanded_SeedCoat <- bayz(Score ~  fixf(Year) +fixf(Location) +fixf(YearLoc)  +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                           ranf(droplevels(GP_Origin), V=GRM_filtered),
                         data = seedcoat_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_SeedCoat)
fit_expanded_SeedCoat$Parameters

fixf3_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seedcoat_filtered$GP_Origin),as.numeric(as.character(seedcoat_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_SeedCoatCol_notcorrected_0516.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_SeedCoatCol_0516.csv",quote=F,row.names = F, col.names = F, sep=",")




################## Seed coat col black vs. not ###################
# Load file
seedcoat=function_to_load_data("SeedCoatColCore.txt")
dim(seedcoat)
head(seedcoat)
seedcoatdf=as.data.frame(seedcoat)

# Make variables ready
seedcoat_ready=makevariables(seedcoatdf)


# Take a look at the data
seedcoat_ready$YearLoc=paste(seedcoat_ready$Location,seedcoat_ready$Year,sep="")
Trials=unique(seedcoat_ready$YearLoc)
length(Trials)

lines=unique(seedcoat_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(seedcoat_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

seedcoat_filtered=seedcoat_ready[-which(seedcoat_ready$GP_Origin=="Kontu"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Lynx"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Taifun"),]
nrow(seedcoat_filtered) #960

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x1"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x3"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="x5"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(seedcoat_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
seedcoat_filtered$Score=as.numeric(as.character(seedcoat_filtered$Score))
nodata=which(is.na(seedcoat_filtered$Score)==T)
seedcoat_filtered=seedcoat_filtered[-nodata,]
nrow(seedcoat_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core403"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core402"),]
seedcoat_filtered=seedcoat_filtered[-which(seedcoat_filtered$GP_Origin=="Core62"),]
nrow(seedcoat_filtered) #900

#histograms of observations from each trial
Trials
NS18=seedcoat_filtered[which(seedcoat_filtered$TRID==11),]
nrow(NS18)

NS19=seedcoat_filtered[which(seedcoat_filtered$TRID==22),]
nrow(NS19)

Se19=seedcoat_filtered[which(seedcoat_filtered$TRID==23),]
nrow(Se19)

#Make all black 2 and all non-black 1
blackidx=which(seedcoat_filtered$Score==5)
nonblackidx=which(seedcoat_filtered$Score!=5)

seedcoat_filtered$Score[blackidx]=2
seedcoat_filtered$Score[nonblackidx]=1


seedcoat_filtered$Score=as.numeric(as.character(seedcoat_filtered$Score))
seedcoat_filtered$lon=as.numeric(as.character(seedcoat_filtered$lon))
seedcoat_filtered$lat=as.numeric(as.character(seedcoat_filtered$lat))
seedcoat_filtered=droplevels(seedcoat_filtered)
nrow(seedcoat_filtered)
seedcoat_filtered=Neighboureffects(seedcoat_filtered)


All=ggplot(data=seedcoat_filtered, aes(x=(seedcoat_filtered$Score))) +
  geom_histogram(binwidth=1, fill = "black",
                 alpha=.8) +
  labs(x="Seed coat (1= not black, 2=black)", y="Observations") +
  geom_vline(xintercept = mean(seedcoat_filtered$Score,na.rm=T), col="green")

pdf("SeedCoatBlack_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=seedcoat_filtered[which(seedcoat_filtered$TRID==11),]
nrow(NS18)

NS19=seedcoat_filtered[which(seedcoat_filtered$TRID==22),]
nrow(NS19)

Se19=seedcoat_filtered[which(seedcoat_filtered$TRID==23),]
nrow(Se19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(seedcoat_filtered$Score)) + ylim(0,320) +
  labs(x="Seed coat (1= not black, 2=black)", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(seedcoat_filtered$Score)) + ylim(0,320) +
  labs(x="Seed coat (1= not black, 2=black)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p4=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) + xlim(0,max(seedcoat_filtered$Score)) + ylim(0,320) +
  labs(x="Seed coat (1= not black, 2=black)", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("SeedCoatColBlack_.pdf", 15, 45)
grid.arrange(p2,p3,p4)
dev.off()



fit_coatcol <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +fixf(YearLoc) +fixf(Location) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores),
                    data = seedcoat_filtered, chain=c(100000, 500, 10))

summary(fit_coatcol) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fixf3_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix11_=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fixf8=fit_coatcol$Estimates[which(rownames(fit_coatcol$Estimates)=="freg.AverageNeiScores"),]


fit_coatcol$Parameters

BLUPS=fit_coatcol$Estimates[7:120,1]
ind_names=rownames(fit_coatcol$Estimates)[7:120]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seedcoat_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seedcoat_filtered$GP_Origin %in% colnames(GRM))==F)
seedcoat_filtered_again=seedcoat_filtered[-phenotypeswithnogeno,]
seedcoat_filtered_again=droplevels(seedcoat_filtered_again)

length(unique(seedcoat_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both


fit_expanded_SeedCoat <- bayz(Score ~  fixf(Year) +fixf(Location) +fixf(YearLoc)  +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +ran2f(GP_Origin,YearLoc) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                                ranf(droplevels(GP_Origin), V=GRM_filtered),
                              data = seedcoat_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_SeedCoat)
fit_expanded_SeedCoat$Parameters

fixf3_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Location%NS"),]
fixf3_
fixf4_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Location%Se"),]
fixf4_
fix7_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Year%2018"),]
fix7_
fix8_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.Year%2019"),]
fix8_
fix9_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%NS2019"),]
fix9_
fix10_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%Se2018"),]
fix10_
fix11_=fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="fixf.YearLoc%Se2019"),]
fix11_
fit_expanded_SeedCoat$Estimates[which(rownames(fit_expanded_SeedCoat$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seedcoat_filtered$GP_Origin),as.numeric(as.character(seedcoat_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_BlackSeeds_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_BlackSeeds_0518.csv",quote=F,row.names = F, col.names = F, sep=",")






################## Seed area ###################
# Load file
seedarea=function_to_load_data("SeedArea.txt")
dim(seedarea)
head(seedarea)
seedareadf=as.data.frame(seedarea)

# Make variables ready
seedarea_ready=makevariables(seedareadf)
seedarea_ready=Neighboureffects(seedarea_ready)


# Take a look at the data
seedarea_ready$YearLoc=paste(seedarea_ready$Location,seedarea_ready$Year,sep="")
Trials=unique(seedarea_ready$YearLoc)
length(Trials)

lines=unique(seedarea_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(seedarea_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

seedarea_filtered=seedarea_ready[-which(seedarea_ready$GP_Origin=="Kontu"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Lynx"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Taifun"),]
nrow(seedarea_filtered) #610

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="x1"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="x3"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="x5"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(seedarea_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
seedarea_filtered$Score=as.numeric(as.character(seedarea_filtered$Score))
nodata=which(is.na(seedarea_filtered$Score)==T)
seedarea_filtered=seedarea_filtered[-nodata,]
nrow(seedarea_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Core403"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Core402"),]
seedarea_filtered=seedarea_filtered[-which(seedarea_filtered$GP_Origin=="Core62"),]
nrow(seedarea_filtered) #567

#histograms of observations from all trials
All=ggplot(data=seedarea_filtered, aes(x=(seedarea_filtered$Score))) +
  geom_histogram(binwidth=10, fill = "black",
                 alpha=.8) +
  labs(x="Seed area", y="Observations") +
  geom_vline(xintercept = mean(seedarea_filtered$Score,na.rm=T), col="green")

pdf("Seed_Areaall.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=seedarea_filtered[which(seedarea_filtered$TRID==11),]
nrow(NS18)

NS19=seedarea_filtered[which(seedarea_filtered$TRID==22),]
nrow(NS19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=10,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,50) +
  labs(x="Seed area", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=10,fill = "black",
                 alpha=.8) + xlim(0,max(TGW_filtered$Score)) + ylim(0,50) +
  labs(x="Seed area", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("SeedArea.pdf", 15, 30)
grid.arrange(p2,p3)
dev.off()


seedarea_filtered$Score=as.numeric(as.character(seedarea_filtered$Score))
seedarea_filtered$lon=as.numeric(as.character(seedarea_filtered$lon))
seedarea_filtered$lat=as.numeric(as.character(seedarea_filtered$lat))
seedarea_filtered=droplevels(seedarea_filtered)
nrow(seedarea_filtered)


fit_Seedarea <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +ran2f(GP_Origin,Year)+freg(AverageNeiScores),
                    data = seedarea_filtered, chain=c(100000, 500, 10))

summary(fit_Seedarea) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fix8_=fit_Seedarea$Estimates[which(rownames(fit_Seedarea$Estimates)=="fixf.Year%2019"),]
fix8_
fit_Seedarea$Estimates[which(rownames(fit_Seedarea$Estimates)=="freg.AverageNeiScores"),]



fit_Seedarea$Parameters

BLUPS=fit_Seedarea$Estimates[5:118,1]
ind_names=rownames(fit_Seedarea$Estimates)[5:118]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seedarea_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seedarea_filtered$GP_Origin %in% colnames(GRM))==F)
seedarea_filtered_again=seedarea_filtered[-phenotypeswithnogeno,]
seedarea_filtered_again=droplevels(seedarea_filtered_again)

length(unique(seedarea_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both


fit_expanded_Seedarea <- bayz(Score ~  fixf(Year)  +ran2f(GP_Origin,Year) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                                ranf(droplevels(GP_Origin), V=GRM_filtered),
                              data = seedarea_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_Seedarea)
fit_expanded_Seedarea$Parameters


fix8_=fit_expanded_Seedarea$Estimates[which(rownames(fit_expanded_Seedarea$Estimates)=="fixf.Year%2019"),]
fix8_
fit_expanded_Seedarea$Estimates[which(rownames(fit_expanded_Seedarea$Estimates)=="freg.AverageNeiScores"),]



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seedarea_filtered$GP_Origin),as.numeric(as.character(seedarea_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_SeedArea_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_SeedArea_0518.csv",quote=F,row.names = F, col.names = F, sep=",")


################## Seed width ###################
# Load file
seedwidth=function_to_load_data("Seedwidthcore.txt")
dim(seedwidth)
head(seedwidth)
seedwidthdf=as.data.frame(seedwidth)

# Make variables ready
seedwidth_ready=makevariables(seedwidthdf)
seedwidth_ready=Neighboureffects(seedwidth_ready)


# Take a look at the data
seedwidth_ready$YearLoc=paste(seedwidth_ready$Location,seedwidth_ready$Year,sep="")
Trials=unique(seedwidth_ready$YearLoc)
length(Trials)

lines=unique(seedwidth_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(seedwidth_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

seedwidth_filtered=seedwidth_ready[-which(seedwidth_ready$GP_Origin=="Kontu"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Lynx"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Taifun"),]
nrow(seedwidth_filtered) #610

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="x1"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="x3"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="x5"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(seedwidth_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
seedwidth_filtered$Score=as.numeric(as.character(seedwidth_filtered$Score))
nodata=which(is.na(seedwidth_filtered$Score)==T)
seedwidth_filtered=seedwidth_filtered[-nodata,]
nrow(seedwidth_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Core403"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Core402"),]
seedwidth_filtered=seedwidth_filtered[-which(seedwidth_filtered$GP_Origin=="Core62"),]
nrow(seedwidth_filtered) #567

#histograms of observations from all trials
All=ggplot(data=seedwidth_filtered, aes(x=(seedwidth_filtered$Score))) +
  geom_histogram(binwidth=0.5, fill = "black",
                 alpha=.8) +
  labs(x="Seed width", y="Observations") +
  geom_vline(xintercept = mean(seedwidth_filtered$Score,na.rm=T), col="green")

pdf("Seed_widthall.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=seedwidth_filtered[which(seedwidth_filtered$TRID==11),]
nrow(NS18)

NS19=seedwidth_filtered[which(seedwidth_filtered$TRID==22),]
nrow(NS19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=0.5,fill = "black",
                 alpha=.8) + xlim(0,max(seedwidth_filtered$Score)) + ylim(0,75) +
  labs(x="Seed width", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=0.5,fill = "black",
                 alpha=.8) + xlim(0,max(seedwidth_filtered$Score)) + ylim(0,75) +
  labs(x="Seed width", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("Seedwidth.pdf", 15, 30)
grid.arrange(p2,p3)
dev.off()


seedwidth_filtered$Score=as.numeric(as.character(seedwidth_filtered$Score))
seedwidth_filtered$lon=as.numeric(as.character(seedwidth_filtered$lon))
seedwidth_filtered$lat=as.numeric(as.character(seedwidth_filtered$lat))
seedwidth_filtered=droplevels(seedwidth_filtered)
nrow(seedwidth_filtered)


fit_Seedwidth <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +ran2f(GP_Origin,Year)+freg(AverageNeiScores),
                     data = seedwidth_filtered, chain=c(100000, 500, 10))

summary(fit_Seedwidth) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fix8_=fit_Seedwidth$Estimates[which(rownames(fit_Seedwidth$Estimates)=="fixf.Year%2019"),]
fix8_
fit_Seedwidth$Estimates[which(rownames(fit_Seedwidth$Estimates)=="freg.AverageNeiScores"),]


fit_Seedwidth$Parameters

BLUPS=fit_Seedwidth$Estimates[5:118,1]
ind_names=rownames(fit_Seedwidth$Estimates)[5:118]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seedwidth_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seedwidth_filtered$GP_Origin %in% colnames(GRM))==F)
seedwidth_filtered_again=seedwidth_filtered[-phenotypeswithnogeno,]
seedwidth_filtered_again=droplevels(seedwidth_filtered_again)

length(unique(seedwidth_filtered_again$GP_Origin))==dim(GRM_filtered) # 98 individuals in both


fit_expanded_SeedWIDTH <- bayz(Score ~  fixf(Year)  +ran2f(GP_Origin,Year) +freg(AverageNeiScores)+ ranf(GP_Origin) +
                                ranf(droplevels(GP_Origin), V=GRM_filtered),
                              data = seedwidth_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_SeedWIDTH)
fit_expanded_SeedWIDTH$Parameters


fit_expanded_SeedWIDTH$Estimates[which(rownames(fit_expanded_SeedWIDTH$Estimates)=="fixf.Year%2019"),]
fit_expanded_SeedWIDTH$Estimates[which(rownames(fit_expanded_SeedWIDTH$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seedwidth_filtered$GP_Origin),as.numeric(as.character(seedwidth_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Seedwidth_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Seedwidth_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



################## Seed length ###################
# Load file
seedlength=function_to_load_data("SeedLengthcore.txt")
dim(seedlength)
head(seedlength)
seedlengthdf=as.data.frame(seedlength)

# Make variables ready
seedlength_ready=makevariables(seedlengthdf)
seedlength_ready=Neighboureffects(seedlength_ready)



# Take a look at the data
seedlength_ready$YearLoc=paste(seedlength_ready$Location,seedlength_ready$Year,sep="")
Trials=unique(seedlength_ready$YearLoc)
length(Trials)

lines=unique(seedlength_ready$GP_Origin)
length(lines) #123 lines in collection
counting_lines=count(seedlength_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

seedlength_filtered=seedlength_ready[-which(seedlength_ready$GP_Origin=="Kontu"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Lynx"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Taifun"),]
nrow(seedlength_filtered) #610

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="x1"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="x3"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="x5"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(seedlength_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
seedlength_filtered$Score=as.numeric(as.character(seedlength_filtered$Score))
nodata=which(is.na(seedlength_filtered$Score)==T)
seedlength_filtered=seedlength_filtered[-nodata,]
nrow(seedlength_filtered) 

# remove lines core403, core402 and core 62 as they consistently behaved different for all observations in one trial compared to ALL other observations, when looking at hilum color
# We concluded a mix up and we do not trust these 3 lines
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Core403"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Core402"),]
seedlength_filtered=seedlength_filtered[-which(seedlength_filtered$GP_Origin=="Core62"),]
nrow(seedlength_filtered) #567

#histograms of observations from all trials
All=ggplot(data=seedlength_filtered, aes(x=(seedlength_filtered$Score))) +
  geom_histogram(binwidth=0.5, fill = "black",
                 alpha=.8) +
  labs(x="Seed lengthth", y="Observations") +
  geom_vline(xintercept = mean(seedlength_filtered$Score,na.rm=T), col="green")

pdf("Seed_Lengthall.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=seedlength_filtered[which(seedlength_filtered$TRID==11),]
nrow(NS18)

NS19=seedlength_filtered[which(seedlength_filtered$TRID==22),]
nrow(NS19)



# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(binwidth=0.5,fill = "black",
                 alpha=.8) + xlim(0,max(seedlength_filtered$Score)) + ylim(0,40) +
  labs(x="Seed length", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=0.5,fill = "black",
                 alpha=.8) + xlim(0,max(seedlength_filtered$Score)) + ylim(0,40) +
  labs(x="Seed length", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("Seedlength.pdf", 15, 30)
grid.arrange(p2,p3)
dev.off()


seedlength_filtered$Score=as.numeric(as.character(seedlength_filtered$Score))
seedlength_filtered$lon=as.numeric(as.character(seedlength_filtered$lon))
seedlength_filtered$lat=as.numeric(as.character(seedlength_filtered$lat))
seedlength_filtered=droplevels(seedlength_filtered)
nrow(seedlength_filtered)


fit_Seedlength <- bayz((as.factor(Score)) ~ ranf(GP_Origin) + fixf(Year) +ran2f(GP_Origin,Year)+freg(AverageNeiScores),
                      data = seedlength_filtered, chain=c(100000, 500, 10))

summary(fit_Seedlength) #line explains 81% when all together.there is a genotype-yearloc (gxe) effect of 10.5%

fit_Seedlength$Estimates[which(rownames(fit_Seedlength$Estimates)=="fixf.Year%2019"),]
fit_Seedlength$Estimates[which(rownames(fit_Seedlength$Estimates)=="freg.AverageNeiScores"),]



fit_Seedlength$Parameters

BLUPS=fit_Seedlength$Estimates[5:118,1]
ind_names=rownames(fit_Seedlength$Estimates)[5:118]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS



#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% seedlength_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((seedlength_filtered$GP_Origin %in% colnames(GRM))==F)
seedlength_filtered_again=seedlength_filtered[-phenotypeswithnogeno,]
seedlength_filtered_again=droplevels(seedlength_filtered_again)

length(unique(seedlength_filtered_again$GP_Origin))==dim(GRM_filtered) # 98 individuals in both


fit_expanded_SeedLENGTH <- bayz(Score ~  fixf(Year)  +ran2f(GP_Origin,Year)+freg(AverageNeiScores) + ranf(GP_Origin) +
                                 ranf(droplevels(GP_Origin), V=GRM_filtered),
                               data = seedlength_filtered_again, chain=c(100000, 500, 100))

summary(fit_expanded_SeedLENGTH)
fit_expanded_SeedLENGTH$Parameters


fit_expanded_SeedLENGTH$Estimates[which(rownames(fit_expanded_SeedLENGTH$Estimates)=="fixf.Year%2019"),]
fit_expanded_SeedLENGTH$Estimates[which(rownames(fit_expanded_SeedLENGTH$Estimates)=="freg.AverageNeiScores"),]


# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(seedlength_filtered$GP_Origin),as.numeric(as.character(seedlength_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_SeedLength_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_SeedLength_0518.csv",quote=F,row.names = F, col.names = F, sep=",")





################## Internode length ###################

# Load file
Internodelength=function_to_load_data("InternodeLengthCore.txt")
dim(Internodelength)
head(Internodelength)
Internodelengthdf=as.data.frame(Internodelength)

# Make variables ready
Internodelength_ready=makevariables(Internodelengthdf)
Internodelength_ready=Neighboureffects(Internodelength_ready)


# Take a look at the data
Internodelength_ready$YearLoc=paste(Internodelength_ready$Location,Internodelength_ready$Year,sep="")
Trials=unique(Internodelength_ready$YearLoc)
length(Trials)

lines=unique(Internodelength_ready$GP_Origin)
length(lines) #119 lines in collection
counting_lines=count(Internodelength_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

iLength_filtered=Internodelength_ready[-which(Internodelength_ready$GP_Origin=="Kontu"),]
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="Lynx"),]
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="Taifun"),]
nrow(iLength_filtered) #610

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="x1"),]
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="x3"),]
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="x5"),]
iLength_filtered=iLength_filtered[-which(iLength_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(iLength_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
iLength_filtered$Score=as.numeric(as.character(iLength_filtered$Score))
nrow(iLength_filtered) #336
nodata=which(is.na(iLength_filtered$Score))
iLength_filtered=iLength_filtered[-nodata,]
nrow(iLength_filtered) #323

#histograms of observations from all trials
All=ggplot(data=iLength_filtered, aes(x=(iLength_filtered$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) +
  labs(x="Internode length (cm)", y="Observations") +
  geom_vline(xintercept = mean(iLength_filtered$Score,na.rm=T), col="green")

pdf("InternodeLength_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

NS19=iLength_filtered[which(iLength_filtered$TRID==22),]
nrow(NS19)

# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) +
  labs(x="Internode length (cm)", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("Internodelength_ns19.pdf", 30, 30)
grid.arrange(p3)
dev.off()


iLength_filtered$Score=as.numeric(as.character(iLength_filtered$Score))
iLength_filtered$lon=as.numeric(as.character(iLength_filtered$lon))
iLength_filtered$lat=as.numeric(as.character(iLength_filtered$lat))
nrow(iLength_filtered) #318
iLength_filtered=droplevels(iLength_filtered)
nrow(iLength_filtered)

# fit model
fit_iLength <- bayz(Score ~ ranf(GP_Origin)+freg(AverageNeiScores),data = iLength_filtered, chain=c(100000, 500, 10))
summary(fit_iLength) #line explains 82% when all together

fit_iLength$Parameters
fit_iLength$Estimates[which(rownames(fit_iLength$Estimates)=="freg.AverageNeiScores"),]



BLUPS=fit_iLength$Estimates[4:115,1]
ind_names=rownames(fit_iLength$Estimates)[4:115]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% iLength_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((iLength_filtered$GP_Origin %in% colnames(GRM))==F)
iLength_filtered_again=iLength_filtered[-phenotypeswithnogeno,]

iLength_filtered_again=droplevels(iLength_filtered_again)

length(unique(iLength_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both

fit_expanded_iLength<- bayz(Score ~  ranf(GP_Origin) +freg(AverageNeiScores)+
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = iLength_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded_iLength)
fit_expanded_iLength$Estimates[which(rownames(fit_expanded_iLength$Estimates)=="freg.AverageNeiScores"),]

fit_expanded_iLength$Parameters
# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(iLength_filtered$GP_Origin),as.numeric(as.character(iLength_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_internodelen_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_internodelen_0518.csv",quote=F,row.names = F, col.names = F, sep=",")







################## Ovules ###################
# Load file
Ovules=function_to_load_data("Ovulesprpod.txt")
dim(Ovules)
head(Ovules)
Ovdf=as.data.frame(Ovules)

# Make variables ready
Ovdf_ready=makevariables(Ovdf)
Ovdf_ready=Neighboureffects(Ovdf_ready)


# Take a look at the data
Ovdf_ready$YearLoc=paste(Ovdf_ready$Location,Ovdf_ready$Year,sep="")
Trials=unique(Ovdf_ready$YearLoc)
length(Trials)

lines=unique(Ovdf_ready$GP_Origin)
length(lines) #119 lines in collection
counting_lines=count(Ovdf_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>(length(Trials)*3)),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

OV_filtered=Ovdf_ready[-which(Ovdf_ready$GP_Origin=="Kontu"),]
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="Lynx"),]
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="Taifun"),]
nrow(OV_filtered) #348

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<(length(Trials)*3)),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="x1"),]
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="x3"),]
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="x5"),]
OV_filtered=OV_filtered[-which(OV_filtered$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(OV_filtered, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 7 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
OV_filtered$Score=as.numeric(as.character(OV_filtered$Score))
nodata=which(is.na(OV_filtered$Score)==T)
OV_filtered=OV_filtered[-nodata,]
OV_filtered=OV_filtered[-which(OV_filtered$Score>8),] #as these two are extreme outliers

nrow(OV_filtered) #323


NS19=OV_filtered[which(OV_filtered$TRID==22),]
nrow(NS19)

# Make histogram plots
# Make histogram plots
#histograms of observations from all trials, save these
p3=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=0.3,fill = "black",
                 alpha=.8) +
  labs(x="Ovules pr. pod", y="Observations") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) +
  scale_x_continuous(breaks = seq(0, 9, 1))

pdf("Ovules_ns19.pdf", 30, 30)
grid.arrange(p3)
dev.off()


OV_filtered$Score=as.numeric(as.character(OV_filtered$Score))
OV_filtered$lon=as.numeric(as.character(OV_filtered$lon))
OV_filtered$lat=as.numeric(as.character(OV_filtered$lat))
nrow(OV_filtered) #318
OV_filtered=droplevels(OV_filtered)
nrow(OV_filtered)

# fit model
fit_Ovules <- bayz(Score ~ ranf(GP_Origin)+freg(AverageNeiScores),data = OV_filtered, chain=c(100000, 500, 10))
summary(fit_Ovules) #

fit_Ovules$Parameters
fit_Ovules$Estimates[which(rownames(fit_Ovules$Estimates)=="freg.AverageNeiScores"),]


BLUPS=fit_Ovules$Estimates[4:115,1]
ind_names=rownames(fit_Ovules$Estimates)[4:115]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% OV_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((OV_filtered$GP_Origin %in% colnames(GRM))==F)
OV_filtered_again=OV_filtered[-phenotypeswithnogeno,]

OV_filtered_again=droplevels(OV_filtered_again)

length(unique(OV_filtered_again$GP_Origin))==dim(GRM_filtered) # 93 individuals in both

fit_expanded_OV<- bayz(Score ~  ranf(GP_Origin) +freg(AverageNeiScores)+
                              ranf(droplevels(GP_Origin), V=GRM_filtered),
                            data = OV_filtered_again, chain=c(10000, 500, 100))

summary(fit_expanded_OV)
fit_expanded_OV$Parameters
fit_expanded_OV$Estimates[which(rownames(fit_expanded_OV$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(OV_filtered$GP_Origin),as.numeric(as.character(OV_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_Ovulesprpod_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_Ovulesprpod_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



################## Branching ###################

# Load file
Branchdf=function_to_load_data("BranchingCore.txt")
dim(Branchdf)
head(Branchdf)
Branchdf=as.data.frame(Branchdf)

# Make variables ready
Branchdf_ready=makevariables(Branchdf)
Branchdf_ready=Neighboureffects(Branchdf_ready)


# Take a look at the data
Branchdf_ready$YearLoc=paste(Branchdf_ready$Location,Branchdf_ready$Year,sep="")
Trials=unique(Branchdf_ready$YearLoc)
length(Trials)

lines=unique(Branchdf_ready$GP_Origin)
length(lines) #124 lines in collection
counting_lines=count(Branchdf_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>9),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="Kontu"),]
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="Lynx"),]
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="Taifun"),]
nrow(Branchdf_ready)

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$freq<9),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="x1"),]
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="x3"),]
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="x5"),]
Branchdf_ready=Branchdf_ready[-which(Branchdf_ready$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Branchdf_ready, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 10 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
Branchdf_ready$Score=as.numeric(as.character(Branchdf_ready$Score))
nrow(Branchdf_ready)
Branchdf_filtered=Branchdf_ready[-which(is.na(Branchdf_ready$Score)==T),]
nrow(Branchdf_filtered)

#histograms of observations from all trials
All=ggplot(data=Branchdf_filtered, aes(x=(Branchdf_filtered$Score))) +
  geom_histogram(binwidth=0.2,fill = "black",
                 alpha=.8) +
  labs(x="Average number of fertile lateral shoots pr. plant", y="Observations") +
  geom_vline(xintercept = mean(Branchdf_filtered$Score,na.rm=T), col="green")

pdf("Branching_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials
NS18=Branchdf_filtered[which(Branchdf_filtered$TRID==11),]
nrow(NS18)

NS19=Branchdf_filtered[which(Branchdf_filtered$TRID==22),]
nrow(NS19)

Se19=Branchdf_filtered[which(Branchdf_filtered$TRID==23),]
nrow(Se19)


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + ylim(0,150) +
  labs(x="Average number of fertile lateral shoots pr. plant", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green", size=2) #looks normally distributed, mean 65 cm

p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + ylim(0,150) +
  labs(x="Average number of fertile lateral shoots pr. plant") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm

p3=ggplot(data=Se19, aes(x=(Se19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + ylim(0,150) +
  labs(x="Average number of fertile lateral shoots pr. plant", y="Observations") +
  ggtitle("Sejet 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(Se19$Score,na.rm=T), col="green",size=2) #looks close to normally distributed, mean 71 cm

pdf("Branching.pdf", 15, 30)
grid.arrange(p1,p2,p3)
dev.off()


Branchdf_filtered$Score=as.numeric(as.character(Branchdf_filtered$Score))
Branchdf_filtered$lon=as.numeric(as.character(Branchdf_filtered$lon))
Branchdf_filtered$lat=as.numeric(as.character(Branchdf_filtered$lat))
Branchdf_filtered=droplevels(Branchdf_filtered)

# Fit phenotypical models in lme4
fit <- bayz((Score) ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores),
            data = Branchdf_filtered, chain=c(100000, 500, 10))

summary(fit) #line explains 22% 
plot(fit)
BLUPS=fit$Estimates[7:123,1]
ind_names=rownames(fit$Estimates)[7:123]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Year%2018"),]
fixf2=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Year%2019"),]
fixf3=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%NS"),]
fixf4=fit$Estimates[which(rownames(fit$Estimates)=="fixf.Location%Se"),]
fixf5=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit$Estimates[which(rownames(fit$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit$Estimates[which(rownames(fit$Estimates)=="freg.AverageNeiScores"),]
fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8

fit$Parameters #to check sizes of parameters

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Branchdf_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Branchdf_filtered$GP_Origin %in% colnames(GRM))==F)
Branchdf_filtered_again=Branchdf_filtered[-phenotypeswithnogeno,]
Branchdf_filtered_again=droplevels(Branchdf_filtered_again)

length(unique(Branchdf_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded <- bayz(Score ~ ranf(GP_Origin) + fixf(Location) + fixf(Year) +fixf(YearLoc) +ran2f(GP_Origin,YearLoc) +ran2f(GP_Origin,Year) +ran2f(GP_Origin,Location) +freg(AverageNeiScores) +
                       ranf(droplevels(GP_Origin), V=GRM_filtered),data = Branchdf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded)
fit_expanded$Parameters #to check sizes of parameters


rownames(fit_expanded$Estimates)
fixf1=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Year%2018"),]
fixf2=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Year%2019"),]
fixf3=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Location%NS"),]
fixf4=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Location%Se"),]
fixf5=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%NS2018"),]
fixf6=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%NS2019"),]
fixf7=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.YearLoc%Se2019"),]
fixf8=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="freg.AverageNeiScores"),]
fixf1
fixf2
fixf3
fixf4
fixf5
fixf6
fixf7
fixf8



# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Branchdf_filtered$GP_Origin),as.numeric(as.character(Branchdf_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_branching_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_branching_0518.csv",quote=F,row.names = F, col.names = F, sep=",")



# only Nordic Seed Data as Sejet seem to score different
Branchdf_filtered_onlyNS=Branchdf_filtered[-which(Branchdf_filtered$TRID==23),]
nrow(Branchdf_filtered_onlyNS)

Branchdf_filtered_onlyNS=droplevels(Branchdf_filtered_onlyNS)

All=ggplot(data=Branchdf_filtered_onlyNS, aes(x=(Branchdf_filtered_onlyNS$Score))) +
  geom_histogram(binwidth=0.2,fill = "black",
                 alpha=.8) +
  labs(x="Average number of fertile lateral shoots pr. plant", y="Observations") +
  geom_vline(xintercept = mean(Branchdf_filtered_onlyNS$Score,na.rm=T), col="green")

pdf("Branching_NS.pdf", 15, 15)
grid.arrange(All)
dev.off()

p1=ggplot(data=NS18, aes(x=(NS18$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + ylim(0,200) +
  labs(x="Average number of fertile lateral shoots pr. plant", y="Observations") +
  ggtitle("Nordic Seed 2018") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS18$Score,na.rm=T), col="green", size=2) #looks normally distributed, mean 65 cm

p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(fill = "black",
                 alpha=.8) + ylim(0,150) +
  labs(x="Average number of fertile lateral shoots pr. plant") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm

pdf("BranchingNS.pdf", 15, 30)
grid.arrange(p1,p2)
dev.off()



# Fit phenotypical models in lme4
fit_NS <- bayz((Score) ~ ranf(GP_Origin) + fixf(Year)+ran2f(GP_Origin,Year) +freg(AverageNeiScores),
            data = Branchdf_filtered_onlyNS, chain=c(100000, 500, 10))


summary(fit_NS)
BLUPS=fit_NS$Estimates[5:121,1]
ind_names=rownames(fit_NS$Estimates)[5:121]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fixf1=fit_NS$Estimates[which(rownames(fit_NS$Estimates)=="fixf.Year%2018"),]
fixf2=fit_NS$Estimates[which(rownames(fit_NS$Estimates)=="fixf.Year%2019"),]
fixf1
fixf2
fit_NS$Estimates[which(rownames(fit_NS$Estimates)=="freg.AverageNeiScores"),]

fit_NS$Parameters #to check sizes of parameters

#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Branchdf_filtered_onlyNS$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Branchdf_filtered_onlyNS$GP_Origin %in% colnames(GRM))==F)
Branchdf_filtered_again=Branchdf_filtered_onlyNS[-phenotypeswithnogeno,]
Branchdf_filtered_again=droplevels(Branchdf_filtered_again)

length(unique(Branchdf_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded <- bayz(Score ~ ranf(GP_Origin) + fixf(Year)  +ran2f(GP_Origin,Year) +freg(AverageNeiScores) +
                       ranf(droplevels(GP_Origin), V=GRM_filtered),data = Branchdf_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded)
fit_expanded$Parameters #to check sizes of parameters
fixf2=fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.Year%2019"),]
fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="freg.AverageNeiScores"),]
fixf2

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Branchdf_filtered_onlyNS$GP_Origin),as.numeric(as.character(Branchdf_filtered_onlyNS$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_branchingOnlyNS_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_branchingOnlyNS_0518.csv",quote=F,row.names = F, col.names = F, sep=",")








################## Lodging ###################

# Load file
Lodging=function_to_load_data("LodgingCore.txt")
dim(Lodging)
head(Lodging)
Lodging=as.data.frame(Lodging)

# Make variables ready
Lodging_ready=makevariables(Lodging)
Lodging_ready=NeighboureffectsMoreDates(Lodging_ready)
Lodging_ready$AverageNeiScores

# Take a look at the data
Lodging_ready$YearLoc=paste(Lodging_ready$Location,Lodging_ready$Year,sep="")
Trials=unique(Lodging_ready$YearLoc)
length(Trials)
Dates=unique(Lodging_ready$date)
length(Dates) #2 dates

lines=unique(Lodging_ready$GP_Origin)
length(lines) #124 lines in collection
counting_lines=count(Lodging_ready, vars = "GP_Origin")

#Check lines with occurence over 9 (3 reps in 3 trials) times
MoreObservations=counting_lines[which(counting_lines$freq>9),]
MoreObservations
# core149 is 6 times in NS2018 can that be right? It is in datafile from Jens
# core30 is 4 times in NS2018 can that be right? It is in datafile from Jens
# Kontu, Lynx, Taifun are not actually core but added to help count field effects. We don't use that

Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="Kontu"),]
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="Lynx"),]
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="Taifun"),]
nrow(Lodging_ready)

#Check lines with occurence less than 9 (3 reps in 3 trials) times
LessObservations=counting_lines[which(counting_lines$n<9),]
LessObservations
#if present 6 times, no problem as NS2018 miss some cores
# Not core lines should be filtered out
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="x1"),]
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="x3"),]
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="x5"),]
Lodging_ready=Lodging_ready[-which(Lodging_ready$GP_Origin=="Fuego"),]

#Check if we have some lines that comes in a number of replicates not dividable by 3
counting_lines=count(Lodging_ready, vars = "GP_Origin")
for (i in seq(1:nrow(counting_lines))){
  dividableby3=counting_lines$freq[i]%%3
  if (dividableby3!=0){
    print(counting_lines[i,])
  }
}

#Core30 has 10 replicates, but we would expect that due to 4 replicates in NS2018


###############Quality check of data##############
Lodging_ready$Score=as.numeric(as.character(Lodging_ready$Score))
nrow(Lodging_ready)
Lod_filtered=Lodging_ready[-which(is.na(Lodging_ready$Score)==T),]
nrow(Lod_filtered)

#histograms of observations from all trials
All=ggplot(data=Lod_filtered, aes(x=(Lod_filtered$Score))) +
  geom_histogram(binwidth=1,fill = "black",
                 alpha=.8) +
  labs(x="0: no lodging 9: heavy lodging", y="Observations") +
  geom_vline(xintercept = mean(Lod_filtered$Score,na.rm=T), col="green")

pdf("Lodging_all.pdf", 15, 15)
grid.arrange(All)
dev.off()

#histograms of observations from each trial
Trials

NS19=Lod_filtered[which(Lod_filtered$TRID==22),]
nrow(NS19)

firstdate=Lod_filtered[which(Lod_filtered$date=="2019-07-17"),]
seconddate=Lod_filtered[which(Lod_filtered$date=="2019-08-15"),]

# Make histogram plots
#histograms of observations from all trials, save these
p2=ggplot(data=NS19, aes(x=(NS19$Score))) +
  geom_histogram(binwidth=1, fill = "black",
                 alpha=.8) + ylim(0,100) + xlim(0,9) +
  labs(x="0: no lodging 9: heavy lodging") +
  ggtitle("Nordic Seed 2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(NS19$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm

pdf("Lodging.pdf", 15, 15)
grid.arrange(p2)
dev.off()


# Make histogram plots
#histograms of observations from all trials, save these
p1=ggplot(data=firstdate, aes(x=(firstdate$Score))) +
  geom_histogram(binwidth=1, fill = "black",
                 alpha=.8) + ylim(0,40) + xlim(0,9) +
  labs(x="0: no lodging 9: heavy lodging") +
  ggtitle("Nordic Seed 17/07/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(firstdate$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm

p2=ggplot(data=seconddate, aes(x=(seconddate$Score))) +
  geom_histogram(binwidth=1, fill = "black",
                 alpha=.8) + ylim(0,40) + xlim(0,9) +
  labs(x="0: no lodging 9: heavy lodging") +
  ggtitle("Nordic Seed 15/08/2019") +
  theme(plot.title = element_text(size = 40, face = "bold"), axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20)) +
  geom_vline(xintercept = mean(seconddate$Score,na.rm=T), col="green",size=2) #looks very scewed towards high plants, mean 98 cm


pdf("Lodging_dates.pdf", 15, 30)
grid.arrange(p1,p2)
dev.off()

Lod_filtered$Score=as.numeric(as.character(Lod_filtered$Score))
Lod_filtered$lon=as.numeric(as.character(Lod_filtered$lon))
Lod_filtered$lat=as.numeric(as.character(Lod_filtered$lat))
Lod_filtered=droplevels(Lod_filtered)

# Fit phenotypical models in lme4
fit_lod <- bayz((Score) ~ ranf(GP_Origin) +freg(AverageNeiScores) +fixf(date),
            data = Lod_filtered, chain=c(100000, 500, 10))

summary(fit_lod) #line explains 22% 
BLUPS=fit_lod$Estimates[4:115,1]
ind_names=rownames(fit_lod$Estimates)[4:115]
BLUPswithnames=cbind(ind_names,BLUPS)
BLUPswithnames=as.data.frame(BLUPswithnames)
BLUPswithnames$BLUPS

fit_lod$Parameters #to check sizes of parameters
fit_lod$Estimates[which(rownames(fit_lod$Estimates)=="fixf.date%2019-08-15"),]
fit_lod$Estimates[which(rownames(fit_lod$Estimates)=="freg.AverageNeiScores"),]


#Hvis man har kinship matrixes (GRM) så kan få broadsense og narrow sense med at lave en ranf() med og uden en GRM, som gives med V=GRM.
GRM=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",header=T,row.names=1)
genotypeswithnopheno=which((colnames(GRM) %in% Lod_filtered$GP_Origin)==F)
GRM_filtered=GRM[-genotypeswithnopheno,-genotypeswithnopheno]
phenotypeswithnogeno=which((Lod_filtered$GP_Origin %in% colnames(GRM))==F)
Lod_filtered_again=Lod_filtered[-phenotypeswithnogeno,]
Lod_filtered_again=droplevels(Lod_filtered_again)

length(unique(Lod_filtered_again$GP_Origin))==dim(GRM_filtered) # 100 individuals in both

fit_expanded <- bayz(Score ~ ranf(GP_Origin) +freg(AverageNeiScores) + fixf(date) +
                       ranf(droplevels(GP_Origin), V=GRM_filtered),data = Lod_filtered_again, chain=c(200000, 1000, 100))

summary(fit_expanded)
fit_expanded$Parameters #to check sizes of parameters
fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="fixf.date%2019-08-15"),]
fit_expanded$Estimates[which(rownames(fit_expanded$Estimates)=="freg.AverageNeiScores"),]

# Extract only BLUPs for GWAS
BLUPswithnames
# delete every part of the name before "Core"
name=separate(BLUPswithnames,"ind_names",c("firsthalf","name"),sep="%")[,2]
BLUPswithnames$ind_names=name
BLUPswithnames

# Output Phenotype file for GWAS
phenotypefileheight_notcor=cbind(as.character(Lod_filtered$GP_Origin),as.numeric(as.character(Lod_filtered$Score)))
phenotypefileheight_notcor=as.data.frame(phenotypefileheight_notcor)
colnames(phenotypefileheight_notcor)=c("Genotype","phenotype")

# genotype average
GenotypeMeans2=aggregate(as.numeric(as.character(phenotypefileheight_notcor$phenotype)),list(phenotypefileheight_notcor$Genotype),mean)
#write.table(GenotypeMeans2,"Phenotypes_Core_lodging_notcorrected_0518.csv",quote=F,row.names = F, col.names = F, sep=",")
#write.table(BLUPswithnames,"BLUPs_Core_lodging_0518.csv",quote=F,row.names = F, col.names = F, sep=",")


