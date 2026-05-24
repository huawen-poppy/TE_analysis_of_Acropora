library(Seurat)
library(data.table)
setwd('/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/TE_analysis/acropora/')
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
te_names=read.table('/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/TE_analysis/acropora/Acropora1/transcripts.txt',header=F,comment.char="")
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

# origin analysis (direct to the genes)
all_data=NormalizeData(all_data)
all_data=ScaleData(all_data)

all_data=FindVariableFeatures(all_data,selection.method = 'vst',nfeatures = 2000)

CC7=readRDS('/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/seurat_analysis/CC7_raw_data_analysis_extra/data/CC7_after_DEG.RDS')
all_data@reductions$umap=CC7@reductions$umap
#all_data=RunPCA(all_data,features = VariableFeatures(all_data))
#all_data=RunUMAP(all_data,dims = 1:20)

DimPlot(all_data,reduction = 'umap',group.by = 'celltype')

# find the marker genes
Idents(all_data)=all_data$group
group_markers=FindAllMarkers(all_data,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)
group_markers=group_markers[group_markers$p_val_adj<0.05,]

group_markers[rownames(group_markers)[grepl('NODE-',rownames(group_markers))],]
# get the NODE-12214-Length-249#Unknown express higher in apo group

DotPlot(all_data,features = 'NODE-12214-Length-249#Unknown',group.by = 'group')
VlnPlot(all_data,features = 'NODE-12214-Length-249#Unknown',group.by = 'group')
all_data$celltype_stim=paste0(all_data$celltype,'_',all_data$group)

Idents(all_data)=all_data$celltype
celltype_markers=FindAllMarkers(all_data,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)
celltype_markers=celltype_markers[celltype_markers$p_val_adj<0.05,]

features=unique(celltype_markers[rownames(celltype_markers)[grepl('NODE-',rownames(celltype_markers))],]$gene)
DotPlot(all_data,features = features,group.by = 'celltype')+RotatedAxis()
FeaturePlot(all_data,features = features,reduction = 'umap')
VlnPlot(all_data,features = features,group.by = 'celltype')

all_data$celltype_stim=paste0(all_data$celltype,'_',all_data$group)
Idents(all_data)=all_data$celltype_stim

for (i in unique(all_data$celltype)) {
  a=FindMarkers(all_data,ident.1 = paste0(i,'_A'),ident.2 = paste0(i,'_S'),min.pct = 0.25,logfc.threshold = 0.25)
  a=a[a$p_val_adj<0.05,]
  assign(paste0(gsub(' ','_',i),'_stim'),a)
  write.csv(a,paste0('./seurat_analysis/DEGs_',gsub(' ','_',i),'.csv'),quote = F,row.names = F)
}

celltype_stim_markers=FindAllMarkers(all_data,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)
celltype_stim_markers
## check whether we need to subset only the TEs to analysis
count_df=all_data@assays$RNA@counts
aipgenes=rownames(count_df)[rownames(count_df) %like% 'AIPGENE']
te_genes=rownames(count_df)[!rownames(count_df) %like% 'AIPGENE']

#aipgenes_count=as.matrix(count_df[aipgenes,])
#tegenes_count=as.matrix(count_df[te_genes,])

boxplot(all_data@meta.data[,c(11,17)],horizontal = T,main='Boxplot for UMIs of TEs and aipgenes')

# sincs that the UMIs and amounts of the TEs are so different from the aipgenes
library(data.table)
te_genes=rownames(all_data)[!rownames(all_data) %like% 'AIPGENE']

te_data=subset(all_data,features = te_genes)
te_data=NormalizeData(te_data)
te_data=ScaleData(te_data)
Idents(te_data)=te_data$celltype
group_markers=FindAllMarkers(te_data,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)
group_markers=group_markers[group_markers$p_val_adj<0.05,]

group_markers[rownames(group_markers)[grepl('NODE-',rownames(group_markers))],]
# get the NODE-12214-Length-249#Unknown express higher in apo group


# make the gene name file
count_data=all_data@assays$RNA@counts
count_data=as.data.frame(as.matrix(count_data))

all_genes=rownames(count_data)
gene_category=c('AIPGENE','Unknown','DNA','LINE','LTR','RC',
                'Satellite','Simple-repeat','SINE','SINE?','snRNA','#tRNA')
for (type in gene_category) {
  assign(paste0('all_',type),all_genes[grepl(type,all_genes,fixed = T)])
}

all_SINE=all_SINE[!grepl('\\?',all_SINE)]

# calculate the sum count table for the category

for (type in gene_category) {
  target_count=count_data[get(paste0('all_',type)),]
  assign(paste0(type,'_SUM'),colSums(target_count))
}

count_data=count_data[all_AIPGENE,]
gene_category=gene_category[2:12]
for (type in gene_category) {
  count_data=rbind(count_data,get(paste0(type,'_SUM')))
  rownames(count_data)[dim(count_data)[1]]=paste0(type,'_SUM')
}

# generate the rds 
data_cate=CreateSeuratObject(counts = count_data)

# add the meta information
data_cate@meta.data$group=all_data@meta.data$group
data_cate@meta.data$sample=all_data@meta.data$sample
data_cate@meta.data$seurat_clusters=all_data@meta.data$seurat_clusters
data_cate@meta.data$celltype=all_data@meta.data$celltype

saveRDS(data_cate,'./seurat_analysis/te_groups.rds')
data_cate=readRDS('./seurat_analysis/te_groups.rds')
data_cate=NormalizeData(data_cate)
data_cate=ScaleData(data_cate)

# dotplot all the te groups expression
te_groups=rownames(data_cate)[!grepl('AIPGENE',rownames(data_cate))]
Idents(data_cate)=data_cate$celltype
DotPlot(data_cate,features=te_groups)+RotatedAxis()
VlnPlot(data_cate,features = te_groups)

Idents(data_cate)=data_cate$group
data_cate$celltype_group=paste0(data_cate$celltype,'_',data_cate$group)

Idents(data_cate)=as.factor(data_cate$celltype_group)

# find markers
group_markers=FindAllMarkers(data_cate,min.pct = 0.25,logfc.threshold = 0.25)

celltype_markers=FindAllMarkers(data_cate,min.pct = 0.25,logfc.threshold = 0.25)

write.csv(group_markers,'te_group_AS_markers.csv',row.names = T)
write.csv(celltype_markers,'te_group_celltype_markers.csv',row.names = T)

for (i in unique(data_cate$celltype)) {
  assign(i,FindMarkers(data_cate,ident.1 = paste0(i,'_A'),ident.2 = paste0(i,'_S'),min.pct = 0.25,logfc.threshold = 0.25))
}

