
#Figure S3A- Volcano plot
res <- read.csv("CDKIControl_vs_DMSOControl_volcano.csv", row.names = 1)

res_plot <- res[order(res$padj), ]
res_plot$col  <- 'gray40'
res_plot$col[res_plot$log2FoldChange > 0 & res_plot$padj < 0.05] <- 'blue'
res_plot$col[res_plot$log2FoldChange < 0 & res_plot$padj < 0.05] <- 'orangered'
par(mar = c(4, 4.4, 5, 10), xpd = TRUE)
plot(res_plot$log2FoldChange, -log10(res_plot$padj), col = res_plot$col, pch = 19, cex =1, xlab = expression(log[2]("FC")), ylab = expression(-log[10]("FDR")), cex.lab=1.2, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, frame.plot = FALSE, box(bty="l"))
legend("topright", inset = c(-0.25, 0), legend = c("Not Significant", "Upregulated", "Downregulated"), col = c("gray", "orangered", "blue"), lty = 2:4, cex = 0.9, pch = 19, bty="n")

###################################################################################################
#Figure S3C- Crossplot

library(readr)
both_na <- read_csv("both_na.csv", col_types = cols(...1 = col_skip()))
View(both_na)
both_na <- na.omit(both_na)

library(ggplot2)

all_plot <-ggplot(both_na, 
                  aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) + 
                  geom_point(size=1.25, shape=16)  + theme_minimal() + coord_fixed() +  
                  geom_vline(xintercept = 0) + 
                  geom_hline(yintercept = 0) + 
                  scale_y_continuous(name="dmso", limits=c(-7, 10)) + 
                  scale_x_continuous(name="cdk", limits=c(-8, 10)) + 
                  scale_color_manual(values=c('orangered','purple', 'grey')) +
                  geom_smooth(method='lm', color="black", size=0.5, formula= y~x)
                  
all_plot

r <- cor(both_na$log2FoldChange.x, both_na$log2FoldChange.y, method = "pearson")
r
#r = 0.793

#FINDING R FOR ONLY DIFF EXP EITHER IN DMSO OR CDK TREATED
library(readr)
dmsoorcdk <- read_csv("dmsoorcdk.csv", col_types = cols(...1 = col_skip()))
View(dmsoorcdk)
dmsoorcdk <- na.omit(dmsoorcdk)

library(ggplot2)

all_plot <-ggplot(dmsoorcdk, 
                  aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) + 
  geom_point(size=1.25, shape=16)  + theme_minimal() + coord_fixed() +  
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous(name="dmso", limits=c(-7, 10)) + 
  scale_x_continuous(name="cdk", limits=c(-8, 10)) + 
  scale_color_manual(values=c('orangered','purple', 'grey')) +
  geom_smooth(method='lm', color="black", size=0.5, formula= y~x)

all_plot

r <- cor(dmsoorcdk$log2FoldChange.x, dmsoorcdk$log2FoldChange.y, method = "pearson")
r
#0.6665309

#Calculate and plot DEGs in both
Both <- subset(both_na, both_na$Legend == "both < 0.05")
all_plot <-ggplot(Both, 
                  aes(x=log2FoldChange.x, y=log2FoldChange.y, color=Legend)) + 
  geom_point(size=1.25, shape=16)  + theme_minimal() + coord_fixed() +  
  geom_vline(xintercept = 0) + 
  geom_hline(yintercept = 0) + 
  scale_y_continuous(name="dmso", limits=c(-7, 10)) + 
  scale_x_continuous(name="cdk", limits=c(-8, 10)) + 
  scale_color_manual(values=c('orangered')) +
  geom_smooth(method='lm', color="black", size=0.5, formula= y~x)

all_plot

r <- cor(Both$log2FoldChange.x, Both$log2FoldChange.y, method = "pearson")
r
#r = 0.947

#heatmap of overlapping DEGs

Both_heatmap <- Both[, -c(1,3,5,6,7,8,9,10,12,13,14,15,16)]
Both_heatmap <- na.omit(Both_heatmap)
#write.csv(Both_heatmap, file="Both_heatmap.csv")
#manually upload, converting log2fc columns to numeric

dat_matrix <- as.matrix(Both_heatmap[1: nrow(Both_heatmap), 2: ncol(Both_heatmap)])
rownames(dat_matrix) <- Both_heatmap$gene.names.x
colnames(dat_matrix) <- colnames(Both_heatmap)[2: ncol(Both_heatmap)]
head(dat_matrix) #looks good: gene.names are row names and log2FC are the row values

#visual
library(pheatmap)
pheatmap(dat_matrix, cluster_cols=F, fontsize_row=2, cellheight=1, cellwidth = 50, cutree_rows = 3)

####################################################################################################
#Figure S3D- TE Correlation plot
#library(readxl)
#DESeq2_TE_only_family <- read_excel("Input_Data/TE_analysis/DESeq2_TE_only_family.xlsx")
#View(DESeq2_TE_only_family)

# better to export from excel as comma separeted and import here as following
# Example:DESeq2_TE_only_family_clean.csv
DESeq2_TE_only_family <- read.csv("./Input_Data/SQuIRE/DESeq2_TE_only_family_clean.csv")
DESeq2_TE_only_locus <- read.csv("./Input_Data/SQuIRE/DESeq2_TE_only_locus.csv")
CDK_TE_only_family <- read.csv("Mutant_vs_Control_CDKI_family.csv")
DMSO_TE_only_family <- read.csv("Mutant_vs_Control_DMSO_family.csv")
# preparation of dataframe 
##################################################################
library("dplyr")
# Separate a column in multiple column base on special characters
library("tidyr")
DESeq2_TE_only_family_sep <- DESeq2_TE_only_family %>% 
  separate(TE_name, c("repName", "repFamily", "repClass"), sep = "([:])", extra = "merge", fill = "right", remove = FALSE)
View(DESeq2_TE_only_family_sep)

DESeq2_TE_only_locus_sep <- DESeq2_TE_only_locus %>% 
  separate(TE_ID, c("chr", "start", "end", "repName", "repFamily", "repClass", "milliDiv", "strand", "extra"), sep = "([|:,])", extra = "merge", fill = "right", remove = FALSE)
View(DESeq2_TE_only_locus_sep)

library("dplyr")
DESeq2_TE_only_family_sep_count <- DESeq2_TE_only_family_sep %>%
  group_by(repClass) %>% count()

DESeq2_TE_only_locus_sep_count <- DESeq2_TE_only_locus_sep %>%
  group_by(repClass) %>% count()

# Do the same for CDK_TE_only_family and DMSO_TE_only_family and save as CDK_TE_only_family_sep_count and DMSO_TE_only_family_sep_count

####################################################################
# Plotting - MA plot
#Merge both CDK_TE_only_family_sep_count and DMSO_TE_only_family_sep_count and save as "both_mergedTEfamily.csv"

postscript("./Results/MAplot/lvr_DMSO_PD_hi272_TE_correlation.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)
pdf("./Results/MAplot/lvr_DMSO_PD_hi272_TE_correlation.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 10, height = 7.5)

both_mergedTEfamily <- read.delim("both_mergedTEfamily.csv")

plot(both_mergedTEfamily$log2FC.x, both_mergedTEfamily$log2FC.y, pch=20, 
     xlab="Log2FC DMSO", ylab="log2FC PD", main="Liver PD - DMSO", col="grey")
with(subset(both_mergedTEfamily, repClass=="DNA"), 
     points(log2FC.x, log2FC.y, pch=20, col="orange"))
with(subset(both_mergedTEfamily, repClass=="LTR"), 
     points(log2FC.x, log2FC.y, pch=20, col="purple1"))
with(subset(both_mergedTEfamily, repClass=="LINE"), 
     points(log2FC.x, log2FC.y, pch=20, col="cadetblue1"))
with(subset(both_mergedTEfamily, repClass=="SINE"), 
     points(log2FC.x, log2FC.y, pch=20, col="dodgerblue1"))
abline(lm(both_mergedTEfamily$log2FC.y ~ both_mergedTEfamily$log2FC.x), lwd=2) 
dev.off()
r <- cor(both_mergedTEfamily$log2FC.y, Pboth_mergedTEfamily$log2FC.x, method = "pearson")
r


####################################################################################################
Figure S3F- REVIGO
#Upload & subset datasets
uhrf1 <- read.csv("DMSOMut12345_vs_DMSOCon1234_HUM.csv")
uhrf1_sig <- subset(uhrf1, uhrf1$padj<0.05)
uhrf1_sig_UP <-  subset(uhrf1_sig, uhrf1_sig$log2FoldChange>0)
uhrf1_sig_down <- subset(uhrf1_sig, uhrf1_sig$log2FoldChange<0)
  
cdkmut <- read.csv("Livhi272CDKMut_Livhi272CDKCON_HUM.csv")
cdkmut_sig <- subset(cdkmut, cdkmut$padj<0.05)
cdkmut_sig_UP <-  subset(cdkmut_sig, cdkmut_sig$log2FoldChange>0)
cdkmut_sig_down <-  subset(cdkmut_sig, cdkmut_sig$log2FoldChange<0)

#Revigo
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(ggplot2)
OrgDb = org.Dr.eg.db

# Curate GO terms
#Functions
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
#####

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv("Revigo_BP_Table.csv")
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$TermID[revi$Representative==-1]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}

# Perform GO on these clusters
fig = list(UP_uhrf1=uhrf1_sig_UP$X, UP_cdk1=cdkmut_sig_UP$zf_gene.names, DOWN_uhrf1=uhrf1_sig_down$X, DOWN_cdk1=cdkmut_sig_down$zf_gene.names)

fig_go = compareCluster(fig, fun = "enrichGO", keyType = "ENSEMBL", OrgDb = "org.Dr.eg.db" , ont = "BP",readable = TRUE)
create_revigo_table(fig_go, "fig_goidswithp.txt")

# Copy and paste this file to REVIGO (http://revigo.irb.hr) (add Danio rerio species) and retrieve summary of GO terms as csv file and save it in the working directory as REVIGO_fig.csv
revi = read.delim("Revigo_BP_Table.csv")
revi=as.data.frame(revi)
db = as.data.frame(fig_go)

ids = db$ID[!db$ID %in% revi$TermID[revi$Dispensability<=0.5]]
dropped_go = dropGO(fig_go,term = ids)

scale = 1.5
wrapsize = 42

#Plot terms - 8 terms
fig_plottable = for_plotting_comparecluster(fig,dropped_go,10,wrapsize) # change number depending on how many top terms you would like to be shown
write.table(fig_plottable, file = "GO_BP_10terms.txt", quote = F, sep = "\t ", row.names = F)

plot_fig = ggplot(fig_plottable, aes(x = Cluster, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  theme_bw(base_size = 14) +
  scale_color_continuous(name=str_wrap('Adjusted p-value', width = wrapsize),low="blue", high="red", guide=guide_colorbar(reverse=TRUE))+
  ylab(NULL) +
  xlab('')+
  scale_x_discrete(labels= c('UP uhrf1', "UP cdk", "DOWN uhrf1", "DOWN cdk"))+
  theme(legend.key.height = unit(30,'points'),
        axis.text.y = element_text(size = 22,angle=0, hjust=1),
        axis.text.x = element_text(angle = 45,hjust=1,size=18),
        text = element_text(size=15),
        axis.title.x =element_text(size=22))+
  scale_size(name='Ratio of Genes\nin Gene List',range=c(6*scale, 18*scale),breaks=c(0.02,0.04,0.08,0.14))+
  theme(strip.text.y = element_blank())+
  geom_text(aes(label=Count,size=0.002,fontface='bold'),colour='white',show.legend = F)

plot_fig

postscript("B.eps", horizontal = FALSE, onefile = FALSE, paper = "special", width = 11, height = 13)
plot_fig
dev.off()

###########################################################################################################################
