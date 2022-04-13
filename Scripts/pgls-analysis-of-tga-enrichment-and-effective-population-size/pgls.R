library(ape)
library(geiger)
library(nlme)
library(phytools)
library(caper)

### PGLS ANALYSIS ###

setwd("./Project_Iso")
data<-read.table("./PGLS/taa_data.txt",header=F,sep="\t", stringsAsFactors = F)
tree<-read.nexus("./PGLS/species_tree.tre")
names(data) <- c("Species", "TAA_enrichment", "TGA_enrichment", "TAG_enrichment", 'GC')
rownames(data) <- data$Species
comp.data<-comparative.data(tree, data, names.col="Species", vcv.dim=3, warn.dropped=TRUE, vcv=TRUE)
model<-pgls(TGA_enrichment~GC, data=comp.data, lambda = "ML")
summary(model)
