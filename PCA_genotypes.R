library(ggfortify)
library(tidyverse)
library(grid)
grid.newpage()




# Core

# load data
d1 <- read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/ForGRM.tfam", sep = " ", header=F)
d2 <- read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GenofileForGRMmatrix.csv", sep=",", header=T)

linenames=d1[,1]
length(linenames)==nrow(d2)
d3=cbind(linenames,d2)
colnames(d3)=c("Line",seq(1:(ncol(d2))))
d3[,1]=as.character(linenames)

df=as.data.frame(d3)
dim(df)

#====================================================================#
#     PCA on core  pop. genotypes                                    #
#====================================================================#

# we can not sort out missing genotypes as no SNPs would go into the model then
# Therefore lets set missing genotypes to the average, we do this by centering genotypes and then set them to 0.
dfscaled=data.matrix(df[,2:ncol(df)])
dfscaled_=scale(dfscaled,scale=F)
dfscaled_[is.na(dfscaled_)]= 0
merged=cbind.data.frame(linenames,dfscaled_)
variables <- merged[2:ncol(merged)]
merged_df=as.data.frame(merged)

#where did the lines come from
origin=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/db_NORFAB/ProFaba_db_test/ACTUAL_TABLES/GP_tableNew.txt",header=F,sep="\t")
head(origin)

merged_df$Origin=merged_df$linenames
merged_df$Origin=as.character(merged_df$Origin)
for (i in seq(1:nrow(merged_df))){
  line=merged_df[i,1]
  line=as.character(line)
  idx=which(as.character(origin$V2)==line)
  if (length(idx)>0){
    donor=as.character(origin$V4)[idx]
    merged_df$Origin[i]=donor
    
  }else{
    merged_df$Origin[i]="\\N"
  }
}
merged_df$Origin=as.factor(merged_df$Origin)
unique(merged_df$Origin)
       
merged_df$Origin[which(merged_df$Origin=="B.V.")]="Bert Vandenberg"
merged_df$Origin[which(merged_df$Origin=="ICaRDA")]="ICARDA"
unique(merged_df$Origin)

#remove zero variant columns
variances=apply(variables, 2, var)
length(variances)
variances=as.numeric(as.character(variances))
length(which(variances==0))
variables1=variables[,-which(variances==0)]
dim(variables1)

# Make PCA plot
rownames(merged_df)=as.character(linenames)
autoplot(prcomp(variables1), data = merged_df, label.size = 4.0) +theme_bw()
# save as 10 x 10 pdf
autoplot(prcomp(variables1), data = merged_df, shape = FALSE, label.size = 1.0)


# try colouring lines after who they came from

library(ggbiplot)
ggbiplot(prcomp(variables1), choices = 1:2 , var.axes =FALSE) +theme_bw()
ggbiplot(prcomp(variables1), choices = 3:4 , var.axes =FALSE) +theme_bw()

pca=prcomp(variables1)
colours=c("red","blue","yellow","orange","navy","forestgreen","black","lightblue","grey","dodgerblue3","gold4","darkslategray1","deeppink2","deeppink4","darkorchid1","darkkhaki","burlywood3","green")
colors <- colours[as.numeric(merged_df$Origin)]

s3d <-scatterplot3d(pca$x[, 1], pca$x[, 2],pca$x[, 3],xlab="Sepal.length",ylab="Sepal.width", zlab="Petal.length", pch = 16,color=colors)
legend("right", legend = unique(merged_df$linenames),
       col =  colors, pch = 16)
#save as 10 x 10

s3d <-scatterplot3d(pca$x[, 1], pca$x[, 2],pca$x[, 3],xlab="Sepal.length",ylab="Sepal.width", zlab="Petal.length", pch = 16)
#save as 10 x 10

rownames(merged_df)=as.numeric(linenames)

numericlinenames=str_split_fixed(linenames, "Core", 2)[,2]
rownames(merged_df)=as.numeric(numericlinenames)

merged_df_new=merged_df
merged_df_new[,1]=numericlinenames

autoplot(prcomp(variables1), data = merged_df_new, colour = 'Origin', shape = FALSE, label.size = 4.0)
#d.factanal <- factanal(variables1, factors = 3, scores = 'regression')
#autoplot(variables1, data = merged_df_new, colour = 'linenames')
         