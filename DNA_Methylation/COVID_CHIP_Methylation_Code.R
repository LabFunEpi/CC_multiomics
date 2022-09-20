
#Install packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("minfi")

install.packages('CpGassoc')

install.packages("devtools")

library(devtools)
install_github("remap-cisreg/ReMapEnrich")

library(ReMapEnrich)

#Processing idat files
require(minfi)
targets<-read.metharray.sheet(baseDir)
RGset<-read.metharray.exp(base=baseDir,targets=targets)
pd<-pData(RGset)

#SWAN Normalization
MSet.swan<-preprocessSWAN(RGset) 
a<-getBeta(MSet.swan)
b<-detectionP(RGset,"m+u")
c<-getMeth(MSet.swan)
d<-getUnmeth(MSet.swan)

#Load sample pheno file containing sample names and group information
groups <- read.csv("GEO_metadata_Methylation.csv")

#Rename samples using sample names from pheno file
a2 <- a
colnames(a2) <- rownames(groups)

#Make sure your pheno file has a Condition column with the disease state (COVID_CHIP vs COVID_Only)
groups2 <- groups$Condition

#Differential Analysis
library(CpGassoc)

sss <- cpg.assoc(as.matrix(a2), indep = groups2)

#extracting results from sss and calling it sss2
sss2 <- sss$results 

#getting mean methylation for each row of Covid/Chip samples
meanCOV.CHIP <- rowMeans(a2[,group=="COVID_CHIP"]) 

#getting mean methylation for each row of Covid only samples
meanCOV.Only <- rowMeans(a2[,group=="COVID_Only"]) 

#calculating the difference between methylation of Covid/Chip and Covid only
sss2$DeltaBeta <- meanCOV.CHIP - meanCOV.Only 

#filtering for CpGs with p < 0.01 and deltaBeta > +/- 0.1
a2sig <- sss2[which(sss2$P.value < 0.01 & abs(sss2$DeltaBeta) > 0.1), ]

#Intersecting significant CpGs with ENCODE cCREs
#Create a bed file containing the significant CpGs and sort using bedtools
# sort -k 1,1 -k2,2n TET2mutCpGs.bed > TET2mutCpGs_sorted.bed

#Sort ENCODE cCRE bed file
# sort -k 1,1 -k2,2n encodeCRE.bed > encodeCRE_sorted.bed

#Intersect the CpGs with the cCREs
# bedtools intersect -a TET2mutCpGs_sorted.bed -b encodeCRE_sorted.bed -wa -wb > TET2mut_cCRE_intersect.bed

#Intersect significant CpGs with ReMAP myeloid library

library(ReMapEnrich)

query <- bedToGranges("TET2mutCpGs.bed")

remapCatalog <- bedToGranges("remap2022_nr_macs2_hg38_v1_0.bed.gz") 

enrichment.df <- enrichment(query, remapCatalog, fractionCatalog=1E-9)

enrichmentDotPlot(enrichment.df)




