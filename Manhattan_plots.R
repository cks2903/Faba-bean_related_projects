
####################################
# Manhattan plots on core lines    #
####################################

# Remember SNP position are arbitrary
# But SNPs in same gene are located next to each other

library(qqman)

setwd("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GWAS_Core/GWAS_Core_20200518_V4")

snpinfo=read.table("/Volumes/ST_MBG-PMg/Cathrine/NorFab/Genotypes_20200311/GWAS_Core/GWAS_Core_20200429_V3/20200423_Height_NotCor.map",sep="\t",head=F)
colnames(snpinfo)=c("CHR","SNP","BPCM","BP")
head(snpinfo)


# now load in GWAS result one trait at a time
files <- list.files(pattern = "\\.ps")

for (i in seq(1:length(files))){
  name=files[i]
  trait=strsplit(name, "\\.")[[1]][1]
  dateending=substr(trait, 8,9)
  trait_withoutdate=strsplit(trait, dateending)[[1]][2]
  
  results=read.table(name,sep="\t",head=F)
  colnames(results)=c("SNP","beta","SE","P")
  head(results)
  
  Merged=merge(results,snpinfo)
  head(Merged)
  Merged=as.data.frame(Merged)

  gwResults=cbind(as.character(Merged$SNP),Merged$CHR,Merged$BP,Merged$P)
  gwResults=as.data.frame(gwResults)
  colnames(gwResults)=c("SNP","CHR","BP","P")
  gwResults$CHR=as.numeric( gwResults$CHR)
  gwResults$BP=as.numeric(as.character(gwResults$BP))
  gwResults$P=as.numeric(as.character(gwResults$P))
  head(gwResults)
  
  filename=paste(trait_withoutdate,"_Manhattan",".pdf",sep="")
  pdf(filename, 150, 50)

  manhattan(gwResults, main = "Manhattan Plot", ylim = c(0, 10), cex = 5.0, 
            cex.axis = 0.9, genomewideline = -log10(0.05/nrow(gwResults)),suggestiveline=F)
            
  dev.off()

  

  filename2=paste(trait_withoutdate,"_qqplot",".pdf",sep="")
  pdf(filename2, 10, 10)
  qq(gwResults$P)
  dev.off()
}
