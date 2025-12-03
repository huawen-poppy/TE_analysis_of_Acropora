library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(SeuratData)
library(SeuratDisk)
theme_set(theme_cowplot())

setwd("/Acropora_sc_analysis/seurat_analysis/data")


# this is for the new genome annotation version
Ahem<-readRDS('Ahem_after_heatmap_resolution_0.08.RDS')
Ahem<-RenameIdents(Ahem,"1"="Unknown","2"="Gastrodermal cell","3"="Calicoblast",
                   "4"="Unknown","5"="Neurons","6"="Nematocytes",
                   "7"="Immune cell","8"="Endosymbiotic cell","9"="Immune cell",
                   "10"="Gland cell")


DimPlot(Ahem,label=F)
#save as Ahem_cluster_with_labels

p1<-DimPlot(Ahem,reduction = "umap",group.by = "orig.ident")
p2<-DimPlot(Ahem,reduction = "umap",label = F,repel = TRUE)
p1+p2
#save as umap_group_cluster_labels

DimPlot(Ahem,reduction = "umap",split.by = "orig.ident")
#save as umap_group_cluster_labels2

DefaultAssay(Ahem)<-"RNA"
Ahem$Cell_type<-Idents(Ahem)
Ahem@meta.data$Cell_type<-as.character(Ahem@meta.data$Cell_type)


# calculate the cell munber in each celltype of each sample
Idents(Ahem)<-"orig.ident"
DimPlot(Ahem,reduction = "umap")

#Ahem$Cell_type<-as.factor(Ahem$Cell_type)
Idents(Ahem)<-"Cell_type"
table(Idents(Ahem)) #how many cells in each cluster
table(Ahem$orig.ident) #how many cells in each group
prop.table(table(Idents(Ahem)))  #proportion of cells in each group
table(Idents(Ahem),Ahem$orig.ident)   #how many cells in each group within cell types
prop.table(table(Idents(Ahem),Ahem$orig.ident),margin = 2) 
#the proportion of cells in each group within each celltype

# plot the cell proportion in each group of celltypes
Ahem@meta.data %>%
  group_by(orig.ident,Cell_type) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent,fill=Cell_type))+
  geom_col()+
  ggtitle("Percentage of cell per cell type")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("Group")+
  ylab("Percentage")
#save as percentage_cell_of_celltype_group


# barplot of each celltype between groups
Ahem@meta.data %>%
  group_by(orig.ident,Cell_type) %>%
  count() %>%
  group_by(Cell_type) %>%
  mutate(percent=n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=orig.ident,y=percent,fill=orig.ident))+
  facet_wrap(~Cell_type)+geom_bar(stat = "identity",width = 0.7)+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 10,colour = "black"))+
  scale_fill_discrete(name="Group")+ylab("Percentage")
# save a condition preference of each cluster in groups  
  
saveRDS(Ahem,"Ahem_after_DEG.RDS")



#####################################################
# below is to identify the celltype marker genes

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
ahem.markers <- FindAllMarkers(Ahem, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

ahem.markers<-ahem.markers[ahem.markers$p_val_adj<0.05,]

write.csv(ahem.markers,"./celltype_markers_all.csv",quote = F,row.names = F)
