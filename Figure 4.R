#Figure 4C- Clustered heat map
#Starting with normalized data 
#Subset the normalized data with all the DEGs and use this dataset for further analysis
# Install & load pheatmap package
library ( DESeq2 )
library ( ggplot2 )
library(tidyverse)
library(RColorBrewer)
library(pheatmap)


install.packages("pheatmap")         
library("pheatmap")

cHeatmap <- read.table("normalized counts_DEG_CDKInh_heatmap_data.txt",
                      header=T, sep="\t")
head(cHeatmap)

rownames(cHeatmap) <- geneExp$X

cHeatmap_matrix <- as.matrix(cHeatmap[2:19])
head(cHeatmap_matrix)

#scale the data to Z-scores (by row)
heat <- t(scale(t(cHeatmap_matrix)))
write.csv(heat, "heat.csv")

#Visualization
# set colour scheme 
#myCol <- colorRampPalette(c('dodgerblue', 'white', 'red'))(100)
myCol <- colorRampPalette(c('royalblue', 'white', 'red3'))(100)


library("pheatmap")
pheatmap(heat) 

##Heatmap with Row Clusters
cluster_cols=FALSE
pheatmap(heat, cutree_rows = 10)

dev.off()


#For getting gene information from clusters
# creating annotation from the dendrogram and export the data of heatmap clustering 
##########################################################################################
library("pheatmap")
my_heatmap_out <-pheatmap(heat, silent = TRUE) # with silent to do NOT plot the heatmap

#Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
my_heatmap_orderedGenes <- rownames(heat[my_heatmap_out$tree_row[["order"]],])

#Re-order original data (samples) to match ordering in heatmap (left-to-right)
my_heatmap_orderedSamples <- colnames(heat[,my_heatmap_out$tree_col[["order"]]])

# Plot Dendrogram (in case many row/genes it is not useful)
plot(my_heatmap_out$tree_row)
plot(my_heatmap_out$tree_col)

# If you want something like gene-to-cluster assignment, you can 'cut' 
# your row dendrogram into a pre-selected number of groups as follows
# 13 groups sorted by cluster numbers
my_heatmap_orderedGenes_cluster <- sort(cutree(my_heatmap_out$tree_row, k=10))
my_heatmap_orderedGenes_cluster <- cutree(my_heatmap_out$tree_row, k=10)

# count number of genes or extract gene list for each cluster
length(which(my_heatmap_orderedGenes_cluster=="2")) # cluster 1: 1840 genes
my_heatmap_orderedGenes_cluster[my_heatmap_orderedGenes_cluster==1]

# convert into data.frame for using in pheatmap
my_heatmap_orderedGenes_cluster_df <- data.frame(my_heatmap_orderedGenes_cluster)
colnames(my_heatmap_orderedGenes_cluster_df) <- c("cluster")

# add a word in front of cluster number so that colors in heatmap are not a progressive gradient 
my_heatmap_orderedGenes_cluster_df$cluster <- paste0("cluster", my_heatmap_orderedGenes_cluster_df$cluster)
##########################################################################################

# plot the heatmap with the annotation from the dendrogram clustering
##########################################################################################
library("pheatmap")
my_cluster_colour <- list(cluster = c(cluster1 = "#a6cee3", cluster2 = "#1f78b4", cluster3 = "#b2df8a",
                                      cluster4 = "#33a02c", cluster5 = "#fb9a99", cluster6 = "#e31a1c",
                                      cluster7 = "#fdbf6f", cluster8 = "#ff7f00", cluster9 = "#cab2d6",
                                      cluster10 = "#6a3d9a", cluster11 = "#FFDC15", cluster12 = "#b15928",
                                      cluster13 = "#808184"))

my_heatmap <- pheatmap(All_TPM_tissue_our_norm_filtered, 
                       annotation_row = my_heatmap_orderedGenes_cluster_df, 
                       annotation_colors = my_cluster_colour, 
                       cutree_rows = 13, show_rownames=F, cellwidth=15)

# creating function to save the object of the heatmap so don't need to run again for save it
save_pheatmap_png <- function(x, filename, width= 10, height= 7.5, units="in", res=300) {
  png(filename, width = width, height = height, units = units, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_eps <- function(x, filename, horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5) {
  postscript(filename, horizontal = horizontal, onefile = onefile, paper = paper, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(my_heatmap, "./Results/Heatmap_public_plusOur/plusJunior/RNA-seq_AllpublicTPM_OcBimac_ClustAll_plusOurAll_plusJunior_RNA-seq_slim15_Annot.png")
save_pheatmap_eps(my_heatmap, "./Results/Heatmap_public_plusOur/plusJunior/RNA-seq_AllpublicTPM_OcBimac_ClustAll_plusOurAll_plusJunior_RNA-seq_slim15_Annot.eps")
##########################################################################################

# add column with Gene_ID for exporting clusters
##########################################################################################
my_heatmap_orderedGenes_cluster_df_ID <- cbind(my_heatmap_orderedGenes_cluster_df,
                                               row.names(my_heatmap_orderedGenes_cluster_df) )
colnames(my_heatmap_orderedGenes_cluster_df_ID) <- c("cluster", "transcript_ID")
write.table(my_heatmap_orderedGenes_cluster_df_ID, file = "./Output_Data/Heatmap_orderedGenes_cluster_Gene_ID.txt",
            quote=F, sep='\t', row.names = F, col.names = T)

save(my_heatmap_orderedGenes_cluster_df_ID, file = "./R_Objects/Heatmap_orderedGenes_cluster_Gene_ID.RData")
##########################################################################################
#To get the clusters...do the same for every cluster c1 to c10
c2 <- my_heatmap_orderedGenes_cluster[my_heatmap_orderedGenes_cluster==2]
write.csv(c2,file = "c2.csv")
#To get the gene information from main list
> Maingenes_c2<- subset(MainGDF, MainGDF$ENSEMBLID %in% c2 $ENSEMBLID)
> write.csv(Maingenes_c2, "Maingenes_c2.csv")

#####################################################################################################
Figure 4D, E
#Subsetting from DEG lists for the clusters 
>setwd("~/Dropbox (Sadler Lab)/Kirsten's shared folder/Bhavani/MANUSCRIPTS/2023_Uhrf1 and replication defects/DATA/Bioinformatics/CDK Inh_RNA Seq/Analysis CDK inhibitor vs DMSO treated_2023/Heat Map Fig 4D")
>MainGDF <- read.table("DEG_CDKInh_both_na.txt",
                      header=T, sep="\t")
>head(MainGDF)
>rownames(MainGDF) <- geneExp$ENSEMBLID
>c10 <- read.csv("c10.csv",
                    header=T, sep="\t") 
>head(c10)
>Maingenes_c10<- subset(MainGDF, MainGDF$ENSEMBLID %in% c10 $ENSEMBLID)
>write.csv(Maingenes_c10, "Maingenes_c10.csv")


# REVIGO ON CLUSTERS
#Revigo - upload packages
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db

# Curate GO terms
hi272_sig <- subset(hi272, hi272$padj_hi272<0.05)
hi272_sig_UP <- subset(hi272_sig, hi272_sig$log2FC_hi272>0)
hi272_sig_DOWN <- subset(hi272_sig, hi272_sig$log2FC_hi272<0)

#Plotting poscor, DEGs in one, other or both
fig = list(UP=hi272_sig_UP$ensembl_id, DOWN=hi272_sig_DOWN$ensembl_id)

#Functions - just load as is
for_plotting_comparecluster = function(list,finalgo,no_top_IDs,text_wrapsize){
  finaldb = as.data.frame(finalgo)
  topids = vector()
  names = names(list)
  
  #extract top ids for each group
  for (i in 1:length(names)){
    name = names[i]
    top = finaldb[finaldb$Cluster==name,]
    top = top[order(top$p.adjust),]
    top = top$ID[1:no_top_IDs]
    topids = c(topids,top)
  }
  
  dropids = unique(finaldb$ID)
  dropids = dropids[!dropids %in% topids]
  #subset GO to only include top IDs
  plotdb = dropGO(finalgo,term=dropids)
  
  library(dplyr)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  
  plot_table = as.data.frame(plotdb)
  plot_table$GeneRatio = as.numeric(sub('/.*','',plot_table$GeneRatio))/as.numeric(sub('.*/','',plot_table$GeneRatio))
  plot_table$BgRatio = as.numeric(sub('/.*','',plot_table$BgRatio))/as.numeric(sub('.*/','',plot_table$BgRatio))
  plot_table$Description = str_wrap(plot_table$Description, width = text_wrapsize)
  plot_table$p.adjust = as.numeric(plot_table$p.adjust)
  plot_table$Description = factor(plot_table$Description,levels=rev(unique(plot_table$Description)))
  
  return(plot_table)
}

# create table for uploading on REVIGO
create_revigo_table = function(go, filename){
  db = as.data.frame(go)
  db = db[order(db$p.adjust),]
  keepgos = data.frame(ID=vector(),p.adjust=vector())
  #create list of unique GO IDs with lowest p-value of all groups
  for (i in 1:nrow(db)){
    if (!db$ID[i] %in% keepgos$ID){
      keepgos = rbind(keepgos,db[i,])
    }
  }
  write.table(keepgos[,c(2,6)],quote = F,sep = ' ',row.names = F,col.names=F,file = filename)
}

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv("Revigo_BP_Table.csv")
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$TermID[revi$Representative==-1]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}


fig_go = compareCluster(fig, fun = "enrichGO", keyType = "ENSEMBL", OrgDb = "org.Dr.eg.db" , ont = "BP",readable = TRUE) #change MF to CC or BP depending on what you want
create_revigo_table(fig_go, "fig_goidswithp.txt")

# Copy and paste this file (fig_goidswithp.txt) to REVIGO (http://revigo.irb.hr) (add Danio rerio species) and retrieve summary of GO terms as csv file and save it in the working directory as REVIGO_fig.csv
revi = read.delim("Revigo_BP_Table.csv")
revi=as.data.frame(revi)
db = as.data.frame(fig_go)

ids = db$ID[!db$ID %in% revi$TermID[revi$Dispensability<=0.5]]
dropped_go = dropGO(fig_go,term = ids)

scale = 1.5
wrapsize = 42

#Plot terms -17 terms
fig_plottable = for_plotting_comparecluster(fig,dropped_go,17,wrapsize) # change number depending on how many top terms you would like to be shown
write.table(fig_plottable, file = "GO_BP_10terms.txt", quote = F, sep = "\t ", row.names = F)

plot_fig = ggplot(fig_plottable, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="magenta", high="cyan", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c("UP", "DOWN"))+
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size = 22,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18),
        text = element_text(size=15),
        axis.title.x =element_text(size=22))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 18*scale),breaks=c(0.02,0.04,0.08,0.14))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.002,fontface='bold'),colour='white',show.legend = F)

plot_fig

postscript("hi272.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 11, height = 13)
plot_fig
dev.off()
#write.csv(fig_plottable, file="Terms_for_Upset.csv")

##################################################################################################




