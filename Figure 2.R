#Figure 2A
#Upload the dataset
hi272 <- read.delim("hi272&dnmt1_lvr_human(2).txt")

#Revigo
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

drop_from_revigo = function(current_go,revigo_filename){
  revi = read.csv("Revigo_BP_Table.csv")
  db = as.data.frame(current_go)
  ids = db$ID[!db$ID %in% revi$TermID[revi$Representative==-1]]
  dropped_go = dropGO(current_go,term = ids)
  return(dropped_go)
}


fig_go = compareCluster(fig, fun = "enrichGO", keyType = "ENSEMBL", OrgDb = "org.Dr.eg.db" , ont = "BP",readable = TRUE) #change MF to CC or BP depending on what you want
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
.........................................................
#Figure 2B
#UPSET of genes on lists
#install.packages("genekitr")
library(genekitr)
library(ggplot2)

#Manually made a list of each term and it's genes. Saved file as .csv
library(readr)
Upset <- read_csv("UPSET.csv")
View(Upset)

set1 <- Upset$`mitotic cell cycle process`
set1 <- na.omit(set1)
length(set1)

set2 <- Upset$`DNA replication`
set2 <- na.omit(set2)
length(set2)

set3 <- Upset$`nuclear division`
set3 <- na.omit(set3)
length(set3)

set4 <- Upset$`regulation of cell cycle process`
set4 <- na.omit(set4)
length(set4)

set5 <- Upset$`cell division`
set5 <- na.omit(set5)
length(set5)


set7 <- Upset$negative.regulation.of.chromosome.organization
set7 <- na.omit(set7)
length(set7)

# five groups
la_gene_list <- list(
  mitotic.cell.cycle.process = set1, DNA.replication = set2, nuclear.division = set3, regulation.of.cell.cycle.process = set4,  cell.division = set5, negative.regulation.of.chromosome.organization = set6)


plotVenn(la_gene_list,
         use_venn = FALSE,
         main_text_size = 15,
         legend_text_size = 8,
         legend_position = 'left'
)

