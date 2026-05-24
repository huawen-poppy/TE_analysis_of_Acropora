library(Seurat)
library(ggplot2)
library(ggpubr)
library(data.table)
cc7=readRDS('../merge_all_samples_subset.rds')
cc7
aip=readRDS('/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/seurat_analysis/CC7_raw_data_analysis_extra/data/CC7_after_DEG.RDS')
aip
cc7@reductions$umap=aip@reductions$umap
DimPlot(cc7,reduction = 'umap',group.by = 'celltype')
DimPlot(cc7,reduction = 'umap',group.by = 'group')

colnames(cc7@meta.data)
draw_data=cc7@meta.data[,c(4,9,17,18)]
ggplot(draw_data, aes(x=celltype, y=nCount_RNA_tes_all,fill=group)) + 
  geom_boxplot(outlier.shape = NA)+
  ylim(c(0,150))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab('Cell type')+
  ylab("Detected TEs' UMI amounts")+
  stat_compare_means(aes(group = group), label = "p.signif")
  
####### read the xenia data
xen=readRDS('/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/TE_analysis/xenia/seurat_analysis/xenia_integrated.rds')
DefaultAssay(xen)<-'RNA'
colnames(xen@meta.data)
count_data=as.matrix(xen@assays$RNA@counts)
tes_name=rownames(count_data)[rownames(count_data) %like% 'NODE-']
tes_count=count_data[tes_name,]
xen_te=CreateSeuratObject(tes_count)
xen_te@meta.data$celltype=xen@meta.data$celltype

xen_te$celltype=gsub('precursors','Progenitors',xen_te$celltype)
xen_te$celltype=gsub('epidermis','Epidermal cell',xen_te$celltype)
xen_te$celltype=gsub('gland','Gland cell',xen_te$celltype)
xen_te$celltype=gsub('gastrodermis','Gastrodermal cell',xen_te$celltype)
xen_te$celltype=gsub('neuron','Neurons',xen_te$celltype)
xen_te$celltype=gsub('muscle','Muscle cell',xen_te$celltype)
xen_te$celltype=gsub('cnidocyte','Nematocytes',xen_te$celltype)
xen_te$celltype=gsub('alga-hosting_cells','Endosymbiotic cell',xen_te$celltype)
xen_te$celltype=gsub('digestive_filaments','Gastrodermal cell',xen_te$celltype)

saveRDS(xen_te,'xenia_pure_te.rds')

draw_data2=xen_te@meta.data[,c(2,3,4)]
ggplot(draw_data2, aes(x=celltype, y=nFeature_RNA)) + 
  geom_boxplot()+
  #ylim(c(0,150))+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab('Cell type')+
  ylab("Detected TEs amounts")+
  stat_compare_means(aes(group=celltype), label = "p.signif", method="t.test", comparisons = combn(1:8, 2, FUN = list)) 


