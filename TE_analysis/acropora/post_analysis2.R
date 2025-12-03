library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
te=readRDS('./pure_te.rds')
te_group=readRDS('pure_te_group.rds')
te_family=readRDS('pure_te_family.rds')
CC7=readRDS('./Ahem_after_annotation.RDS')

DimPlot(CC7,reduction = 'umap',group.by = 'Cell_type')
DimPlot(te,reduction = 'umap',group.by = 'celltype')

# compute the signature score
# first is to calculate the differential expression TEs for each cell type
target=te
Idents(target)<-target$celltype
DefaultAssay(target)<-'RNA'
te_markers<-FindAllMarkers(target,only.pos = F)

te_markers<-te_markers[te_markers$p_val_adj<0.05,]
dim(te_markers)
write.csv(te_markers,'./data/pure_te_markers_celltype.csv',quote = F)

marker_anno=read.csv('./data/te_info_annotation.csv',header = T,row.names = 1)
te_markers=merge(te_markers,marker_anno,by.x = 'gene',by.y = 'all_te_copy',all.x = T)
te_markers$group=gsub('rRNA|RC|tRNA|Simple_repeat|RNA|snRNA|Satellite','Others',te_markers$group)

te_marker_category=te_markers[,c(7,10)]
plot_input=te_marker_category %>%
  group_by(cluster,group) %>%
  count()

ggplot(plot_input, aes(fill=group, y=n, x=cluster)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of differentially expressed features')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  RotatedAxis()
ggsave('./plots/pure_te_number_differential_express_features_celltype_te_category_barplot.png')

head(te_marker_category)

# signature score
Idents(te)<-te$celltype
for (i in unique(te_markers$cluster)) {
  test=AddModuleScore(te,features = list(te_markers[te_markers$cluster==i,]$gene),ctrl = 5,name = gsub(' ','_',i))
  #test@meta.data$celltype<-as.factor(test@meta.data$celltype)
  a=unique(test$celltype)
  test$celltype<-factor(test$celltype,levels = c(setdiff(a,i),i))
  ggplot(test@meta.data,aes(x=celltype,y=test@meta.data[,ncol(test@meta.data)]))+
    geom_boxplot()+#outlier.shape = NA)+
    ylab('Signature scores')+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
    RotatedAxis()+
    ggtitle(colnames(test@meta.data)[ncol(test@meta.data)])+
    stat_compare_means(aes(group=celltype), label = "p.signif", method="t.test", comparisons = combn(1:10,2,FUN = list)[c(9,17,24,30,35,39,42,44,45)])
  ggsave(paste0('./plots/signature_score_te_',gsub(" ",'_',i),'.png'))
}

## for heatmap different groups
suppressPackageStartupMessages({
  library(rlang)
})

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by)]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    group.use <- groups.use[, c(i, additional.group.by), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(foo=rep(x = levels(x = group.use[[i]]), times = lines.width))
      placeholder.groups[additional.group.by] = NA
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.levels <- levels(x = group.use[[i]])
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    #group.use = group.use[order(group.use[[i]]), , drop=F]
    group.use <- group.use[with(group.use, eval(parse(text=paste('order(', paste(c(i, additional.group.by), collapse=', '), ')', sep='')))), , drop=F]
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))){
        if (colname == group.by){
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))), "#FFFFFF",group.bar=T)
        } else {
          cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))),group.bar=T)
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = grid::gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        #temp <- as.data.frame(cols[[colname]][levels(group.use[[colname]])])
        #colnames(temp) <- 'color'
        #temp$x <- temp$y <- 1
        #temp[['name']] <- as.factor(rownames(temp))
        
        #temp <- ggplot(temp, aes(x=x, y=y, fill=name)) + geom_point(shape=21, size=5) + labs(fill=colname) + theme(legend.position = "bottom")
        #legend <- get_legend(temp)
        #multiplot(plot, legend, heights=3,1)
        
        if ((colname == group.by) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

te_markers<-read.csv('./data/pure_te_markers_celltype.csv',header = T,row.names = 1)
te_markers %>%
  group_by(cluster) %>%
  top_n(n=20,wt=avg_log2FC) -> top20

DoMultiBarHeatmap(te, assay = 'integrated', features = top20$gene, group.by='celltype')

plot_input=test@meta.data
ggplot(plot_input,aes(x=celltype,y=Cluster1))+
  geom_violin()


### plot distance
te_distance=read.csv('./data/same_strain_distance_df_for_plot_remove_duplicate_gene_distance_group_with_label.csv',header = T,row.names = 1)
te_distance=te_distance[,c(2,3,4,5,6,7)]
head(te_distance)
#colnames(te_distance)[1]<-'distance'
te_distance$type<-with(te_distance,ifelse(distance>2000,'Distal','Proximal'))

te_distance$cluster_te=paste0(te_distance$cluster,'_',te_distance$te_id)
#te_distance_input<-te_distance[,c(1,2,6,8)]
te_distance_new=te_distance %>%
  group_by(cluster_te) %>%
  slice_min(order_by = distance)

te_distance_new<-as.data.frame(te_distance_new)
#te_distance_new$cluster=gsub('_NODE-.*','',te_distance_new$cluster_te)
#te_distance_new$type<-with(te_distance_new,ifelse(distance>2000,'Distal','Proximal'))

te_distance_new<-te_distance_new[!duplicated(te_distance_new$cluster_te),]

#te_distance_new<-te_distance_new[,c(4,5)]
te_distance_new=as.data.frame(te_distance_new)
write.csv(te_distance_new,'./data/marker_te_cloest_genes_distance_2k_all_gene.csv') # all gene is that not only one cloest gene
plot_input=te_distance_new %>%
  dplyr::group_by(cluster,type) %>%
  dplyr::count()

library(plyr)
plot_input2<-ddply(plot_input,.(cluster),transform,percent=n/sum(n)*100)
# Format the labels and calculate their positions
plot_input2 = ddply(plot_input2, .(cluster), transform, pos = (cumsum(n) - 0.5 * n))
plot_input2$label = paste0(sprintf("%.2f", plot_input2$percent), "%")

ggplot(plot_input2, aes(fill=type, y=n, x=cluster)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of TEs')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  RotatedAxis()+
  ggtitle('Distance to protein-coding genes')+
  #scale_y_continuous(labels = scales::percent) +
  geom_text(aes(label = label), 
            position = position_stack(vjust = 0.5), size = 3.5)

ggsave('../plots/pure_te_number_TEs_celltype_distance_to_coding_genes_barplot_correct_celltype_group.png')

# add the gene open or not information to plot_input2
DefaultAssay(CC7)<-'RNA'

te_distance_new$gene_open='no'
for (i in 1:12537){
  target_gene=gsub("_",'-',te_distance_new$gene_id[i])
  target_celltype=te_distance$cluster[i]
  target_cells=rownames(CC7@meta.data[CC7@meta.data$Cell_type==target_celltype,])
  if (target_gene %in% rownames(CC7)){
    target_invest=CC7@assays$RNA@counts[target_gene,target_cells]
    if (sum(target_invest)>0) {
      te_distance_new$gene_open[i]='yes'
    }
  } 
}

table(te_distance_new$gene_open)

te_distance_new$celltype_gene_open=paste0(te_distance_new$cluster,'_',te_distance_new$gene_open)
te_distance_new<-data.frame(te_distance_new)
#te_distance_new_test<-te_distance_new[,c(10,7)]
te_distance_new$type<-with(te_distance_new,ifelse(distance>2000,'Distal','Proximal'))
a=te_distance_new %>%
  group_by(celltype_gene_open,type) %>%
  dplyr::summarise(distance_type = n())

a$celltype=gsub('_.*','',a$celltype_gene_open)
a$gene_open=gsub('.*_','',a$celltype_gene_open)

plot_data=a %>%
  group_by(celltype_gene_open) %>%
  mutate(percent=distance_type/sum(distance_type))

write.csv(plot_data,'./data/plot_data_for_distance_gene_open_2k-allgenes.csv')
write.csv(te_distance_new,'./data/te_gene_distance_open_info-2k-allgene.csv')

plot_data$gene_open<-gsub('no','TE+Gene-',plot_data$gene_open)
plot_data$gene_open<-gsub('yes','TE+Gene+',plot_data$gene_open)
library(webr)
for (celltype in unique(plot_data$celltype)){
  target=plot_data[plot_data$celltype==celltype,]
  #target=target[,c(2,6,3)]
  a=PieDonut(target, aes(gene_open,type,count=distance_type), 
             title = celltype,explode = 1,showPieName = F,
             pieLabelSize = 5,donutLabelSize = 5)
}

plot_data2=plot_data %>%
  group_by(celltype_gene_open) %>%
  dplyr::summarise(open_type=sum(distance_type))

plot_data2<-as.data.frame(plot_data2)

plot_data2$celltype=gsub('_.*','',plot_data2$celltype_gene_open)
plot_data2=plot_data2 %>%
  group_by(celltype) %>%
  mutate(open_percent=open_type/sum(open_type))
plot_data2=as.data.frame(plot_data2)

final=merge(plot_data,plot_data2,by='celltype_gene_open')

write.csv(final,'./data/data_for_dounut_plot.csv')

# plot the detected umis for genes and tes
origin_meta<-CC7@meta.data
te_meta<-te@meta.data
colnames(te_meta)[2]<-'nCount_TE'
colnames(te_meta)[3]<-'nFeature_TE'
te_meta$id<-rownames(te_meta)
origin_meta$id<-rownames(origin_meta)
origin_meta$id<-gsub('-1','',origin_meta$id)
plot_data<-merge(origin_meta,te_meta,by='id')
plot_data<-plot_data[,c(2,3,4,13,14,15)]
colnames(plot_data)[6]<-'celltype'

library(tidyr)
data_long<-gather(plot_data,nCount_type,nCount_value,nCount_RNA:nCount_TE,factor_key = T)

library(reshape2)
data_long<-melt(plot_data,id.vars = c('orig.ident.x','nFeature_RNA','nFeature_TE','celltype'))

ggplot(data_long, aes(fill=variable, y=value, x=orig.ident.x)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of detected features')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

aggregate(data_long$nFeature_RNA,by=list(data_long$orig.ident.x),FUN=sum)

data_long2<-melt(plot_data,id.vars = c('orig.ident.x','nCount_RNA','nCount_TE','celltype'))
ggplot(data_long2, aes(fill=variable, y=value, x=orig.ident.x)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of detected features')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

write.csv(plot_input2,'./data/plot_data_for_dounut_plot.csv')
