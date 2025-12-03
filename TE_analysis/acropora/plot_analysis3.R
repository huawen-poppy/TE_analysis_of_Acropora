library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
te=readRDS('./pure_te.rds')
CC7=readRDS('./Ahem_after_annotation.RDS')

target=te
Idents(target)<-target$celltype
DefaultAssay(target)<-'RNA'
te_markers<-FindAllMarkers(target,only.pos = F)
te_markers<-te_markers[te_markers$p_val_adj<0.05,]

target2=CC7
Idents(target2)<-target2$Cell_type
DefaultAssay(target2)<-'RNA'
ahem_markers<-FindAllMarkers(target2,only.pos = F)
ahem_markers<-ahem_markers[ahem_markers$p_val_adj<0.05,]

for (celltype in unique(ahem_markers$cluster)){
  celltype=as.character(unique(ahem_markers$cluster))[7]
  a<-ahem_markers[ahem_markers$cluster==celltype,] %>%
    top_n(n=20,wt=avg_log2FC)
  target_gene<-a$gene
  b<-te_markers[te_markers$cluster==celltype,] %>%
    top_n(n=20,wt=avg_log2FC)
  target_te<-b$gene
  
  target_cell_ahem<-rownames(target2@meta.data[target2@meta.data$Cell_type==celltype,])
  target_cell_te<-rownames(target@meta.data[target@meta.data$celltype==celltype,])
  
  expr_matrix_te<-target@assays$RNA@counts[rownames(target@assays$RNA) %in% target_te,colnames(target@assays$RNA) %in% target_cell_te]
  expr_matrix_ahem<-target2@assays$RNA@counts[rownames(target2@assays$RNA) %in% target_gene,colnames(target2@assays$RNA) %in% target_cell_ahem]
  expr_matrix_te<-as.data.frame(as.matrix(expr_matrix_te))
  expr_matrix_ahem<-as.data.frame(as.matrix(expr_matrix_ahem))
  
  expr_matrix_ahem<-scale(t(expr_matrix_ahem))
  expr_matrix_te<-scale(t(expr_matrix_te))
  #expr_matrix_ahem$source<-'Gene'
  #expr_matrix_te$source<-'TE'
  colnames(expr_matrix_ahem)<-gsub('-1','',colnames(expr_matrix_ahem))
  #combined_matrix<-rbind(expr_matrix_ahem,expr_matrix_te)
  #combined_matrix<-as.data.frame(t(combined_matrix))
  #combined_matrix<-scale(combined_matrix)
  #cor_matrix<-cor(combined_matrix)
  cor_matrix<-cor(expr_matrix_ahem,expr_matrix_te)
  #cor_matrix[is.na(cor_matrix)] <- 0
  # Perform hierarchical clustering using "ward" method
  #hc <- hclust(as.dist(1 - cor_matrix), method = "ward.D2")
  
  library(pheatmap)
  library(circlize)
  cols<-colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  breaks <- seq(-1, 1, length.out = length(cols) + 1)
  pheatmap(cor_matrix, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
           colormap = cols)
  cellname=gsub(' ',"_",celltype)
  #ggsave(paste0(cellname,"_heatmap_top20_DEG_DETE.png"), plot = last_plot(), device = "png", width =16, height = 16, units = "in")
}
