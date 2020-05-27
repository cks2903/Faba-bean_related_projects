#Make GRM from genotype file
args=commandArgs(trailingOnly = TRUE)


geno <- read.table(args[1],sep=",",header=T)
dim(geno)

# VanRaden (2009) method 1
M=dim(geno)[2]
N=dim(geno)[1]

geno_plus1=geno+1

CALCULATE_P_MULTIPLIED_BY_Q=function(SNP_column){
  SNP_column=data.matrix(SNP_column)
  countsfirstallele=sum(SNP_column,na.rm=T)
  p_freq=countsfirstallele/(length(na.omit(SNP_column))*2)
  q_freq=1-p_freq
  value=p_freq*q_freq
  return(value)
}

# Apply function to all columns
# Sum across results
# 
results=mapply(CALCULATE_P_MULTIPLIED_BY_Q, geno_plus1)
SUM=sum(results)

# 
Z <- scale(geno_plus1,scale=FALSE) #Centered genotype matrix. Every column (SNP) get a mean of 0
Z[is.na(Z)] <- 0 #Set NA to zero in scaled matrix. That is setting them equal to avg.

#Z is the same as calculated in Vanraden method 1. You would get the same if you took 
#M-P wheere P is a matrix where each locus (i)
# has the same number all the way down and that number is 2(pi-0.5)
# I double checked this and it is indeed the same
# Now the mean effect of each SNP is 0

#Calculate GRM matrix
G=(Z%*%t(Z))/(2*SUM)


write.table(G,file="GRM_CoreCollection20200403_VanRaden1.csv",sep=",",col.names=F,row.names=F,quote=F)





# VanRaden (2009) method 1, but only SNPs with MAF>0.05

FINDLOWMAFS=function(SNP_column){
  SNP_column=data.matrix(SNP_column)
  countsfirstallele=sum(SNP_column,na.rm=T)
  p_freq=countsfirstallele/(length(na.omit(SNP_column))*2)
  
  if (p_freq<=0.05 | p_freq>=0.95){
    return(1)
  }
  else {
  	return(0)
  }
}

MAFFilter=mapply(FINDLOWMAFS, geno_plus1)
idxtoremove=which(MAFFilter==1)
length(idxtoremove) #10659 SNPs removed due to minor allele threshold


geno_filtered=geno_plus1[,-idxtoremove]
M_filtered=dim(geno_filtered)[2]
N=dim(geno_filtered)[1]
Z_filtered <- scale(geno_filtered[,-1],scale=FALSE) #Centered genotype matrix. Every column (SNP) get a mean of 0
Z_filtered[is.na(Z_filtered)] <- 0 #Set NA to zero in scaled matrix. That is setting them equal to avg.

results_filtered=mapply(CALCULATE_P_MULTIPLIED_BY_Q, geno_filtered)
SUM_filtered=sum(results_filtered)

#Calculate GRM matrix
G_MAFfilter=(Z_filtered%*%t(Z_filtered))/(2*SUM_filtered)
write.table(G_MAFfilter,file="GRM_CoreCollection20200403_VanRaden1MAFfilter.csv",sep=",",col.names=F,row.names=F,quote=F)

