library(Seurat)
library(data.table)
setwd('/TE_analysis/acropora/')
all_data=readRDS('./Ahem_after_annotation.RDS')
a1=readRDS('./Acropora1/A1_origin.rds')
a2=readRDS('./Acropora2/A2_origin.rds')
a3=readRDS('./Acropora3/A3_origin.rds')
a1@meta.data$orig.ident='a1'
a2@meta.data$orig.ident='a2'
a3@meta.data$orig.ident='a3'

# filter the data
al_names=colnames(all_data)[grepl("-1_1",colnames(all_data))]
a2_names=colnames(all_data)[grepl("-1_2",colnames(all_data))]
a3_names=colnames(all_data)[grepl("-1_3",colnames(all_data))]
a1_names=gsub('-1_1','',al_names)
a2_names=gsub('-1_2','',a2_names)
a3_names=gsub('-1_3','',a3_names)
a1=subset(a1,cells=a1_names)
a2=subset(a2,cells=a2_names)
a3=subset(a3,cells=a3_names)

selected_f<-rownames(a1)[Matrix::rowSums(a1)>3]
a1=subset(a1,features = selected_f)

selected_f<-rownames(a2)[Matrix::rowSums(a2)>3]
a2=subset(a2,features = selected_f)

selected_f<-rownames(a3)[Matrix::rowSums(a3)>3]
a3=subset(a3,features = selected_f)


te_gene_all=merge(a1,c(a2,a3))

# get all tes as gene set
te_genes=rownames(te_gene_all)[rownames(te_gene_all) %like% 'NODE-']
te_data=subset(te_gene_all,features = te_genes)

# split the dataset by sample and then integrate them
te_list=SplitObject(te_data,split.by = 'orig.ident')
te_list<-lapply(te_list, FUN=function(x) {
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x,selection.method = 'vst',nfeatures=2000)
})
features<-SelectIntegrationFeatures(te_list)
te_anchors<-FindIntegrationAnchors(te_list,anchor.features = features)
te_combine<-IntegrateData(te_anchors)
DefaultAssay(te_combine)<-'integrated'

te_combine=ScaleData(te_combine)
#te_combine=FindVariableFeatures(te_combine,selection.method = 'vst',nfeatures = 2000)
te_combine=RunPCA(te_combine,npcs=30)
te_combine=RunUMAP(te_combine,dims = 1:30,reduction = 'pca')
# cell type annotation information
label_info=all_data@meta.data[,c('orig.ident',"Cell_type")]
label_info$id=gsub('-1','',rownames(label_info))
unique(label_info$id==rownames(te_combine@meta.data))
te_combine@meta.data$celltype=label_info$Cell_type

DimPlot(te_combine,reduction = 'umap',group.by = 'celltype')
ggsave('./umap_pure_te_celltype.png')
saveRDS(te_combine,'pure_te.rds')


### get all the te groups 
### first get all the te annotation
te_names=read.table('/TE_analysis/acropora/Acropora1/transcr> ts.txt',header=F,comment.char="")
te_names=te_names[te_names$V1 %like% 'NODE',]
te_names=as.data.frame(te_names)
te_names$te=gsub('#.*','',te_names$te_names)
te_names$te=gsub('_','-',te_names$te)

te_matrix=as.matrix(te_data@assays$RNA@counts)
all_te_copy=rownames(te_matrix)
all_te_copy=as.data.frame(all_te_copy)
all_te_copy$id=1:nrow(all_te_copy)
all_te_copy=merge(all_te_copy,te_names,by.x='all_te_copy',by.y='te',all.x=T)
all_te_copy=all_te_copy[order(all_te_copy$id),]
all_te_copy$group=gsub('.*#','',all_te_copy$te_names)
all_te_copy$group=gsub('/.*','',all_te_copy$group)
all_te_copy$family=all_te_copy$te_names
all_te_copy$family=gsub('.*#','',all_te_copy$family)
all_te_copy$family=gsub('.*/','',all_te_copy$family)
write.csv(all_te_copy,'te_info_annotation.csv')

# genreate the te group matrix 
all_te_copy$group_final=all_te_copy$group
all_te_copy$group_final=gsub('tRNA|Satellite|Simple-repeat|snRNA|RC','Other',all_te_copy$group_final)
all_te_copy$group_final=gsub('\\?','',all_te_copy$group_final)
rownames(te_matrix)<-all_te_copy$group_final
te_matrix_group=rowsum(te_matrix,group = rownames(te_matrix))

te_group=CreateSeuratObject(te_matrix_group)
te_group$group=te_data$group
te_group$sample=te_data$sample
te_group$celltype=te_data$celltype

# since there are only 6 genes
# we cannot integrate the dataset using seurat algorithm
# then we just merge all the dataset as a whole
te_group<-NormalizeData(te_group)
te_group<-FindVariableFeatures(te_group,selection.method = 'vst',nfeatures = 6)
te_group<-ScaleData(te_group)
te_group=RunPCA(te_group,npcs=10)
te_group=RunUMAP(te_group,features = rownames(te_matrix_group),reduction = 'pca')
DimPlot(te_group,reduction = 'umap',group.by = 'celltype')
saveRDS(te_group,'pure_te_group.rds')

# generate the pure te family rds files
rownames(te_matrix)<-all_te_copy$family
te_matrix_family=rowsum(te_matrix,group = rownames(te_matrix))

te_family=CreateSeuratObject(te_matrix_family)
te_family$group=te_data$group
te_family$sample=te_data$sample
te_family$celltype=te_data$celltype
te_list=SplitObject(te_family,split.by = 'sample')
te_list<-lapply(te_list, FUN=function(x) {
  x<-NormalizeData(x)
  x<-FindVariableFeatures(x,selection.method = 'vst',nfeatures=94)
})
features<-SelectIntegrationFeatures(te_list)
te_anchors<-FindIntegrationAnchors(te_list,anchor.features = features)
te_combine<-IntegrateData(te_anchors)
DefaultAssay(te_combine)<-'integrated'

te_combine=ScaleData(te_combine)
#te_combine=FindVariableFeatures(te_combine,selection.method = 'vst',nfeatures = 2000)
te_combine=RunPCA(te_combine,npcs=20)
te_combine=RunUMAP(te_combine,dims = 1:10,reduction = 'pca')
DimPlot(te_combine,reduction = 'umap',group.by = 'celltype')
saveRDS(te_combine,'pure_te_family.rds')

