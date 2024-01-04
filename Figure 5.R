
# Figure 5A- Clustered Heat Map
#Subsetting cell cycle genes
setwd("~/Dropbox (Sadler Lab)/Kirsten's shared folder/Bhavani/MANUSCRIPTS/2023_Uhrf1 and replication defects/DATA/Bioinformatics/CDK Inh_RNA Seq/Analysis CDK inhibitor vs DMSO treated_2023/Heat Map Fig 4A")
MainGDF <- read.table("dnmt1_DEG_lvr_human_dup removed.txt",
                      header=T, sep="\t")
head(MainGDF)
rownames(MainGDF) <- MainGDF$zfin_gene_id
UPSET <- read.table("UPSET genes from figure 2B.txt",
                header=T, sep="\t") 
head(UPSET)
subsetUPSET_dnmt1mut<- subset(MainGDF, MainGDF$zfin_gene_id %in% UPSET $gene.names)
write.csv(subsetUPSET_dnmt1mut, "subsetUPSET_dnmt1mut.csv")

#Heatmap
library ( DESeq2 )
library ( ggplot2 )
library(tidyverse)
library(RColorBrewer)
library(pheatmap)


install.packages("pheatmap")         
library("pheatmap")

cHeatmap <- read.csv("vsd(4).csv",
                       header=T, sep=",")
head(cHeatmap)

rownames(cHeatmap) <- cHeatmap$X

cHeatmap_matrix <- as.matrix(cHeatmap[2:7])
head(cHeatmap_matrix)

#scale the data to Z-scores (by row)
heat <- t(scale(t(cHeatmap_matrix)))
write.csv(heat, "heat.csv")

#Visualization
# set colour scheme 
#myCol <- colorRampPalette(c('dodgerblue', 'white', 'red'))(100)
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)

library("pheatmap")
cluster_cols=FALSE
pheatmap(heat) 

##Heatmap with Row Clusters
cluster_cols=FALSE
pheatmap(heat, cutree_rows = 10)

dev.off()
