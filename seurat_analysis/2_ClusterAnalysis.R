library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Acropora_sc_analysis/seurat_analysis/data")
#load the Seurat Object for each sample
filter_samples<-readRDS("filter_samples.RDS")

symbols<-read.delim("Ahem_tabulated_annots.tsv",sep = '\t')
symbols<-symbols[,c(1,4)]
head(symbols)

library(stringr)
symbols$Hit.description=str_extract(symbols$Hit.description,"(?<=GN=).*(?= PE)")

write.csv(symbols,'Acropora_genesymbol.csv',quote = F,row.names = F)

id2symbol<-read.csv("Acropora_genesymbol.csv",header = T)
samples_algae<-readRDS("filter_samples_algae.RDS")

colnames(id2symbol)<-c("gene","symbol")
#features<-SelectIntegrationFeatures(object.list = filter_samples)

# combine the data
combine.anchors<-FindIntegrationAnchors(object.list = filter_samples,anchor.features = 2000,dims = 1:30)
combine.integrated<-IntegrateData(anchorset = combine.anchors)

combine.anchors_algae<-FindIntegrationAnchors(object.list = samples_algae,anchor.features = 2000,dims = 1:30)
combine.integrated_algae<-IntegrateData(anchorset = combine.anchors_algae)

target_cell_name<-rownames(combine.integrated@meta.data)
target_cell_algae<-subset(combine.integrated_algae@meta.data,
                          rownames(combine.integrated_algae@meta.data) %in% target_cell_name)
colnames(target_cell_algae)[1:3]<-c("ident","UMI_count_algae","gene_count_algae")
target_cell_algae<-target_cell_algae[,1:3]

combine.integrated@meta.data<-merge(combine.integrated@meta.data,
                                    target_cell_algae,by="row.names",all.x = T)
combine.integrated@meta.data[is.na(combine.integrated@meta.data)]<-0
combine.integrated@meta.data<-combine.integrated@meta.data[,-5]
rownames(combine.integrated@meta.data)<-combine.integrated@meta.data$Row.names
combine.integrated@meta.data<-combine.integrated@meta.data[,-1]

combine.integrated@meta.data$row_name<-rownames(combine.integrated@meta.data)
combine.integrated@meta.data<-
  combine.integrated@meta.data[match(target_cell_name,combine.integrated@meta.data$row_name),]

DefaultAssay(combine.integrated)<-'integrated'

# summary the metadata in combine.integrated
tapply(combine.integrated@meta.data$UMI_count_algae,
       combine.integrated@meta.data$orig.ident,summary)

# run the standard workflow for visualization and clustering
combine.integrated <- ScaleData(combine.integrated, verbose = FALSE)
combine.integrated <- RunPCA(combine.integrated, npcs = 30, verbose = FALSE)
combine.integrated <- RunUMAP(combine.integrated, reduction = "pca", dims = 1:20)
combine.integrated <- RunTSNE(combine.integrated, reduction = "pca", dims = 1:20,check_duplicates=FALSE)
DimPlot(combine.integrated,reduction = "pca")

combine.integrated <- FindNeighbors(combine.integrated, reduction = "pca", dims = 1:20)

# check which resolution is best
test_cell<-FindClusters(object = combine.integrated,resolution = c(0.08,0.1,0.3,0.5,0.7,0.9,1.1,1.3),algorithm = 1,group.singletons = F)
Idents(test_cell)<-'integrated_snn_res.0.08'
DimPlot(test_cell,reduction = 'umap',label = T,label.size = 6)
library(clustree)
clustree(test_cell,prefix='integrated_snn_res.')

combine.integrated <- FindClusters(combine.integrated, resolution = 0.3,algorithm = 1,group.singletons = F)

saveRDS(combine.integrated,"combine.integrated_reso_0.5.RDS")
saveRDS(combine.integrated_algae,"combine.integrated_algae.RDS")

# Visualization
p1 <- DimPlot(combine.integrated, reduction = "umap",group.by = "orig.ident")
p2 <- DimPlot(combine.integrated, reduction = "umap", label = F, repel = TRUE)
p1 + p2
#save as umap_sample_clusters.png

# visualize two conditions side by side
DimPlot(combine.integrated, reduction = "umap", split.by = "orig.ident")
#save as umap_split_samples.png


Ahem<-combine.integrated
rm(combine.integrated)


saveRDS(Ahem,"Ahem_feature2000_algorithm1.RDS") #new is one used variable 2000 features to combine the data

Ahem@meta.data$algae_ratio=Ahem@meta.data$UMI_count_algae/(Ahem@meta.data$nCount_RNA+Ahem@meta.data$UMI_count_algae)

# check each cluster's algae ratio
aa=Ahem@meta.data

# boxplot the UMI_count_algae per samples in group
ggplot(aa,aes(orig.ident,UMI_count_algae,fill=factor(orig.ident)))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  labs(fill="Samples")

# boxplot the UMI count algae per group  
ggplot(aa,aes(orig.ident,UMI_count_algae))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("Group")

# barplot the algae UMI count ratio per cell in clusters
ggplot(aa,aes(x=seurat_clusters,y=algae_ratio,group=
                reorder(row.names(aa),as.numeric(seurat_clusters))))+
  geom_bar(stat = "identity",show.legend = F,position = position_dodge(preserve = "single"),aes(fill=seurat_clusters))+
  theme(axis.title.y = element_text("algae UMI ratio"),axis.title.x = element_text("seurat cluster cells"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylim(0,0.3)+
  xlab("Clusters")+
  ylab("Algae UMI count ratio")

# boxplot the algae UMI count ratio per cluster (remove the influences of the cell numbers)  
ggplot(aa,aes(x=seurat_clusters,y=algae_ratio,fill=factor(seurat_clusters)))+
  geom_boxplot()+
  xlab("Clusters")+
  ylab("Algae UMI count ratio")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# boxplot the algae UMI count per cluster
ggplot(aa, aes(x=seurat_clusters,y=UMI_count_algae,fill=factor(seurat_clusters)))+
  geom_boxplot(outlier.size = 0.4)+
  theme_bw()+
  xlab("Clusters")+
  ylab("Algae UMI count")+
  theme(legend.position = "none")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

# barplot the average algae UMI count/cell in each cluster
b1<-aggregate(x=aa$UMI_count_algae,by=list(aa$seurat_clusters),FUN=mean)
colnames(b1)<-c("cluster","mean_algae_UMI_count")
b1$cluster<-factor(b1$cluster,levels=c(0:15))
ggplot(b1,aes(x=cluster,y=mean_algae_UMI_count,group=reorder(row.names(b1),as.numeric(cluster))))+
  geom_bar(stat = "identity",show.legend = F,
           position = position_dodge(preserve = "single"),aes(fill=cluster))+
  theme(axis.title.y = element_text("Average algae UMI count/cell"),
        axis.title.x = element_text("Clusters"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  xlab("Clusters")+
  ylab("Average algae UMI count/cell")
ggsave("../plots/resolution_0.5/barplot_average_algae_UMI_count_per_cell.png")
#save as barplot_average_algae_UMI_count/cell


# identify markers for each cluster and plot heatmap 
Ahem@meta.data$seurat_clusters<-as.numeric(as.character(Ahem@meta.data$seurat_clusters))+1
Idents(Ahem)<-"seurat_clusters"
Ahem@active.ident<-factor(x=Ahem@active.ident,levels = 1:16)

# visualization
DimPlot(Ahem,reduction = "tsne",label = F) #save as tsne_cluster_with_labels.png
DimPlot(Ahem,reduction = 'umap',label = F) #save as umap_cluster_with_labels.png

DimPlot(Ahem, reduction = "tsne", split.by = "orig.ident")
#save as Ahem_tsne_samples_clusters.png

# find all the DE genes in each cluster and choose top 10 of each cluster
# choose the DE gene's top 35
DefaultAssay(Ahem)<-"RNA"
Ahem.markers <- FindAllMarkers(Ahem,only.pos = T,min.pct = 0.25,thresh.use=0.25)
Ahem.markers %>% group_by(cluster) %>% top_n(10,avg_log2FC) -> top10
Ahem@meta.data[,'seurat_clusters',drop=FALSE] %>% add_rownames() %>% group_by(seurat_clusters) %>% sample_n(35) -> index

id2symbol$gene<-gsub(pattern = "_",replacement = "-",x=id2symbol$gene)
order_info<-top10$gene


top10<-merge(top10,id2symbol,by.x="gene",by.y="gene",all.x=T)
for (i in 1:dim(top10)[1]){
  if (is.na(top10[i,8])){
    top10[i,8]=top10[i,1]
  }
}
top10_backup<-top10
top10<-top10[match(order_info,top10$gene),]

#all_genes<-rownames(Ahem)
Ahem<-NormalizeData(Ahem,verbose = F)
Ahem<-ScaleData(Ahem,features=all_genes)
DoHeatmap(Ahem,features = top10$gene,label = F)
DoHeatmap(Ahem,features = top10$gene,label = F,cells = index$rowname)

Ahem.markers<-Ahem.markers[Ahem.markers$p_val_adj<0.05,]
Ahem.markers<-as.data.frame(Ahem.markers)
Ahem.markers$index<-paste0(Ahem.markers$gene,"_",Ahem.markers$cluster)
order_info<-Ahem.markers$index
Ahem.markers<-merge(Ahem.markers,id2symbol,by.x="gene",by.y = "gene")
for (i in 1:dim(Ahem.markers)[1]){
  if (is.na(Ahem.markers[i,9])){
    Ahem.markers[i,9]=Ahem.markers[i,1]
  }
}
Ahem.markers<-Ahem.markers[match(order_info,Ahem.markers$index),]
DoHeatmap(Ahem,features=Ahem.markers$gene,label=F) #too ugly

write.csv(Ahem.markers,"resolution_0.5_cluster_markers_padj_val_0.05_integrated.csv")
saveRDS(Ahem,'Ahem_after_heatmap_resolution_0.5.RDS')
#write.csv(top10,"top10_markergene.csv",quote = F)


#'''
#Marker specificity
#We define that a gene is enriched in a cell if the expression of the gene 
#in the cell is higher than the average expression of the gene among all cells. 
#We further define the specificity of a marker gene in a cluster as the percentage 
#of cells , in which the marker gene is enriched as defined above within the cluster.
#'''
exp <- Ahem@assays$RNA@data
exp_spcificity <- exp >Matrix::rowMeans(exp)
temp_list <- list()
for (i in 1:16){
  sub_cluster <- Ahem.markers[Ahem.markers$cluster ==i,]
  marker_exp <- exp_spcificity[Ahem.markers[Ahem.markers$cluster ==i,]$gene,]
  spcificity_within_cluster <- Matrix::rowSums(marker_exp[,colnames(marker_exp) %in% rownames(Ahem@meta.data[Ahem@meta.data$seurat_clusters == i,])])/sum(Ahem@meta.data$seurat_clusters == i)
  spcificity_outside_cluster <- Matrix::rowSums(marker_exp[,!colnames(marker_exp) %in% rownames(Ahem@meta.data[!Ahem@meta.data$seurat_clusters == i,])])/sum(!Ahem@meta.data$seurat_clusters == i)
  temp <- data.frame(gene=rownames(marker_exp),spcificity_within_cluster=spcificity_within_cluster,spcificity_outside_cluster=spcificity_outside_cluster)
  sub_cluster <- left_join(sub_cluster, temp)
  temp_list[[i]]<-sub_cluster
}

Ahem.markers<- do.call(rbind,temp_list)

#Filter low specificity markers
Ahem.markers <- Ahem.markers[Ahem.markers$spcificity_within_cluster>0.5,]

#write.csv(Ahem.markers,"cluster_markers_specificity.RDS")
write.csv(Ahem.markers,"resolution_0.5_cluster_markers_specificity_with_filter_integrated.csv",row.names = T)


