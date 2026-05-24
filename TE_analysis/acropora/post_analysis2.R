library(Seurat)
library(data.table)
library(ggplot2)
library(dplyr)
library(ggpubr)
te=readRDS('./pure_te.rds')
te_group=readRDS('pure_te_group.rds')
te_family=readRDS('pure_te_family.rds')
CC7=readRDS('./Ahem_after_annotation.RDS')
CC7@meta.data$Cell_type=factor(CC7@meta.data$Cell_type,levels=c('Unknown','Gastrodermal cell','Calicoblast','Neurons','Nematocytes','Immune cell','Endosymbiotic cell','Gland cell'))
p=DimPlot(CC7,reduction = 'umap',group.by = 'Cell_type',pt.size = 1.2)
ggsave('./fig1b_umap.png',plot = p,width = 8,height = 6,dpi=300)

Idents(CC7)=CC7$Cell_type
celltype_markers=FindAllMarkers(CC7,only.pos = F,min.pct = 0.25,logfc.threshold = 0.25)
celltype_markers=celltype_markers[celltype_markers$p_val_adj<0.05,]

features=unique(celltype_markers[rownames(celltype_markers)[grepl('NODE-',rownames(celltype_markers))],]$gene)
DotPlot(all_data,features = features,group.by = 'celltype')+RotatedAxis()

te@meta.data$celltype=factor(te@meta.data$celltype,levels=c('Unknown','Gastrodermal cell','Calicoblast','Neurons','Nematocytes','Immune cell','Endosymbiotic cell','Gland cell'))
p=DimPlot(te,reduction = 'umap',group.by = 'celltype',pt.size = 1.2)
ggsave('./fig3a_umap.png',plot = p,width = 8,height = 6,dpi=300)


# compute the signature score
# first is to calculate the differential expression TEs for each cell type
target=te
Idents(target)<-target$celltype
DefaultAssay(target)<-'RNA'
te_markers<-FindAllMarkers(target,only.pos = F)

te_markers<-te_markers[te_markers$p_val_adj<0.05,]
dim(te_markers)
write.csv(te_markers,'./pure_te_markers_celltype.csv',quote = F)

marker_anno=read.csv('./te_info_annotation.csv',header = T,row.names = 1)
te_markers=merge(te_markers,marker_anno,by.x = 'gene',by.y = 'all_te_copy',all.x = T)
te_markers$group=gsub('rRNA|RC|tRNA|Simple_repeat|RNA|snRNA|Satellite','Others',te_markers$group)

te_marker_category=te_markers[,c(7,10)]
plot_input=te_marker_category %>%
  dplyr::group_by(cluster,group) %>%
  count()

ggplot(plot_input, aes(fill=group, y=freq, x=cluster)) + 
  geom_bar(position="stack", stat="identity")+
  ylab('Number of differentially expressed features')+
  theme_bw()+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  RotatedAxis()
ggsave('./plots/pure_te_number_differential_express_features_celltype_te_category_barplot.png')

plot_input$group <- factor(
  plot_input$group,
  levels = c("DNA","LINE", "LTR","SINE", "Others","Unknown")
)

p <- ggplot(plot_input, aes(x = cluster, y = n, fill = group)) +
  geom_bar(
    stat = "identity",
    position = "stack",
    width = 0.75,
    color = "white",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Unknown" = "#CC79A7",
      "SINE"    = "#0072B2",
      "Others"  = "#56B4E9",
      "LTR"     = "#009E73",
      "LINE"    = "#E69F00",
      "DNA"     = "#D55E00"
    )
  ) +
  #scale_fill_brewer(palette = 'Dark2')+
  labs(
    x = NULL,
    y = "Number of differentially expressed features",
    fill = "Group"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    axis.title.y = element_text(size = 20, color = "black", face = "plain"),
    axis.text.x = element_text(
      size = 18,
      color = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 18, color = "black"),
    
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    legend.key.size = unit(0.5, "cm"),
    legend.position = 'top',
    plot.margin = margin(8, 8, 8, 8)
  )

p
ggsave(
  "fig2c_stacked_bar_DE_features.png",
  plot = p,
  width = 8,
  height = 7,
  dpi = 300
)

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
  ggsave(paste0('./signature_score_te_',gsub(" ",'_',i),'.png'))
}

for (i in unique(te_markers$cluster)) {
  
  # Add module score
  test <- AddModuleScore(
    te,
    features = list(te_markers[te_markers$cluster == i, ]$gene),
    ctrl = 5,
    name = paste0(gsub(" ", "_", i), "_score")
  )
  
  # Get the newly added score column
  score_col <- colnames(test@meta.data)[ncol(test@meta.data)]
  
  # Put the target cell type at the end of the x-axis
  celltype_levels <- unique(test$celltype)
  test$celltype <- factor(
    test$celltype,
    levels = c(setdiff(celltype_levels, i), i)
  )
  
  # Compare target cell type with all other cell types
  comparisons_use <- lapply(
    setdiff(levels(test$celltype), i),
    function(x) c(x, i)
  )
  
  p <- ggplot(
    test@meta.data,
    aes(x = celltype, y = .data[[score_col]])
  ) +
    geom_boxplot(
      width = 0.65,
      outlier.size = 0.4,
      outlier.alpha = 0.6,
      linewidth = 0.45,
      fill = "white",
      color = "black"
    ) +
    stat_compare_means(
      comparisons = comparisons_use,
      label = "p.signif",
      method = "t.test",
      size = 4.5,
      bracket.size = 0.35,
      tip.length = 0.01
    ) +
    labs(
      x = NULL,
      y = "Signature scores",
      title = gsub("_", " ", score_col)
    ) +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 18,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 14,
        color = "black",
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(
        size = 14,
        color = "black"
      ),
      axis.line = element_line(
        color = "black",
        linewidth = 0.6
      ),
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.5
      ),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8)
    )
  
  print(p)
  
  ggsave(
    filename = paste0("./signature_score_te_", gsub(" ", "_", i), ".png"),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

## drawing the median detected/expressed genes and te in the cell types
celltype_col <- "celltype"
DefaultAssay(te)='RNA'
CC7@meta.data$celltype=CC7@meta.data$Cell_type
# Calculate detected feature number per cell
gene_mat <- GetAssayData(
  CC7,
  assay = DefaultAssay(CC7),
  slot = "counts"
)

te_mat <- GetAssayData(
  te,
  assay = DefaultAssay(te),
  slot = "counts"
)

# 4. Calculate detected gene/TE number per cell
qc_df <- data.frame(
  cell = colnames(CC7),
  Cell_type = CC7@meta.data[[celltype_col]],
  detected_gene = Matrix::colSums(gene_mat > 0),
  detected_TE = Matrix::colSums(te_mat > 0)
)

# 5. Calculate median detected gene and TE number per cell type
plot_df <- qc_df %>%
  group_by(Cell_type) %>%
  summarise(
    median_detected_gene = median(detected_gene, na.rm = TRUE),
    median_detected_TE = median(detected_TE, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

# 6. Plot
library(ggrepel)


p <- ggplot(
  plot_df,
  aes(
    x = median_detected_TE,
    y = median_detected_gene,
    label = Cell_type
  )
) +
  geom_point(
    aes(color = Cell_type),
    size = 7,
    alpha = 0.9
  ) +
  geom_text_repel(
    size = 6,
    colour = "black",
    max.overlaps = Inf,
    box.padding = 0.6,
    point.padding = 0.8
  ) +
  labs(
    x = "Median number of detected TEs",
    y = "Median number of detected genes",
    color = "Cell type"
  ) +
  theme_classic(base_size = 20) +
  theme(
    text = element_text(colour = "black"),
    
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    
    axis.title.x = element_text(
      size = 22,
      colour = "black",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 22,
      colour = "black",
      margin = margin(r = 10)
    ),
    
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    
    plot.margin = margin(5,5,5,5)
  )+
  theme(legend.position = "none")

p

ggsave(
  "fig3b_median_detected_gene_vs_TE_colored_by_celltype.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

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

te_markers<-read.csv('./pure_te_markers_celltype.csv',header = T,row.names = 1)
te_markers %>%
  group_by(cluster) %>%
  top_n(n=20,wt=avg_log2FC) -> top20

DoMultiBarHeatmap(te, assay = 'integrated', features = top20$gene, group.by='celltype')

plot_input=test@meta.data
ggplot(plot_input,aes(x=celltype,y=Cluster1))+
  geom_violin()


### plot distance
te_distance=read.csv('./same_strain_distance_df_for_plot_remove_duplicate_gene_distance_group_with_label.csv',header = T,row.names = 1)
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

#te_distance_new<-te_distance_new[!duplicated(te_distance_new$cluster_te),]

#te_distance_new<-te_distance_new[,c(4,5)]
te_distance_new=as.data.frame(te_distance_new)
write.csv(te_distance_new,'./marker_te_cloest_genes_distance_2k_all_gene.csv') # all gene is that not only one cloest gene
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
ggsave('./pure_te_number_TEs_celltype_distance_to_coding_genes_barplot_correct_celltype_group.png')
library(ggplot2)
library(Seurat)

p <- ggplot(plot_input2, aes(fill = type, y = n, x = cluster)) + 
  geom_bar(
    position = "stack",
    stat = "identity",
    width = 0.75,
    colour = "black",
    linewidth = 0.25
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Number of closest genes to DETEs",
    fill = NULL
  ) +
  RotatedAxis() +
  #scale_fill_manual(
  #  values = c(
  #    "Distal" = "#CC79A7",
  #    "Proximal"    = "#0072B2",
  #  ))+
  theme_classic(base_size = 22) +
  theme(
    text = element_text(colour = "black"),
    
    plot.title = element_text(
      size = 24,
      face = "bold",
      colour = "black",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    
    axis.text.x = element_text(
      size = 20,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 18,
      colour = "black"
    ),
    
    axis.title.y = element_text(
      size = 24,
      colour = "black",
      #face = "bold",
      margin = margin(r = 10)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    legend.title = element_blank(),
    legend.text = element_text(
      size = 18,
      colour = "black"
    ),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "top",
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    plot.margin = margin(10, 15, 10, 10)
  )

p

ggsave(
  "./fig4b_pure_te_number_TEs_celltype_distance_to_coding_genes_barplot_correct_celltype_group.png",
  plot = p,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "white"
)

# add the gene open or not information to plot_input2
DefaultAssay(CC7)<-'RNA'

te_distance_new$gene_open='no'
for (i in 1:dim(te_distance_new)[1]){
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

write.csv(plot_data,'./plot_data_for_distance_gene_open_2k-allgenes.csv')
write.csv(te_distance_new,'./te_gene_distance_open_info-2k-allgene.csv')

plot_data$gene_open<-gsub('no','TE+Gene-',plot_data$gene_open)
plot_data$gene_open<-gsub('yes','TE+Gene+',plot_data$gene_open)
library(webr)
safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]", "_", x)
}
plot_data$gene_open=as.factor(plot_data$gene_open)
for (celltype in unique(plot_data$celltype)){
  target=plot_data[plot_data$celltype==celltype,]
  
  #target=target[,c(2,6,3)]
  #PieDonut(target, aes(gene_open,type,count=distance_type), 
  #           title = celltype,explode = 1,showPieName = F,
  #           pieLabelSize = 5,donutLabelSize = 5,ratioByGroup = F)
  p=PieDonut(
    target,
    aes(
      pies = gene_open,
      donuts = type,
      count = distance_type
    ),
    #title = celltype,
    showPieName = FALSE,

    # Important: remove or greatly reduce explode
    explode = 0.4,
    explodeDonut = TRUE,
    
    # Larger labels for publication
    pieLabelSize = 9,
    donutLabelSize = 8,
    titlesize = 8,showRatioThreshold = F,labelposition = 10
  )

  
  recorded_plot <- recordPlot()
  png(filename = paste0(
    "./donut_TE_gene_distance_",
    celltype,
    ".png"
  ), width = 10, height = 8, units = "in", res = 300)
  replayPlot(recorded_plot)
  dev.off()
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

write.csv(final,'./data_for_dounut_plot.csv')

# plot the detected umis for genes and tes
origin_meta<-CC7@meta.data
te_meta<-te@meta.data
colnames(te_meta)[2]<-'nCount_TE'
colnames(te_meta)[3]<-'nFeature_TE'
te_meta$id<-rownames(te_meta)
origin_meta$id<-rownames(origin_meta)
origin_meta$id<-gsub('-1','',origin_meta$id)
plot_data<-merge(origin_meta,te_meta,by='id')
plot_data<-plot_data[,c(2,3,4,14,15,16)]
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

# for the venn plot 
library(eulerr)
library(grid)

# Input exact region counts
venn_counts <- c(
  "All closest open genes" = 4459,
  "All DEGs" = 328,
  "TF-coding genes" = 152,
  "All closest open genes&All DEGs" = 205,
  "All closest open genes&TF-coding genes" = 18,
  "All DEGs&TF-coding genes" = 7,
  "All closest open genes&All DEGs&TF-coding genes" = 1
)

fit <- euler(venn_counts)

png(
  "fig4c_venn_no_set_labels.png",
  width = 10,
  height = 7,
  units = "in",
  res = 300,
  bg = "white"
)
plot(
  fit,
  fills = list(
    fill = c("#F8766D", "#7CAE00", "#619CFF"),
    alpha = 0.65
  ),
  edges = list(
    col = "white",
    lwd = 0.8
  ),
  labels = F,
  quantities = list(
    cex = 2.5,
    font = 1,
    col = "black"
  )
)

dev.off()

# plot sup fig2
consensus_fasta <- "TEconsensus.fa.classified"
repeatmasker_out <- "acropora_genome.fa.out"
genome_fai <- "acropora_genome.fa.fai"

library(tidyverse)
library(data.table)
library(patchwork)
library(scales)
options(scipen = 999)

genome_size <- fread(genome_fai, header = FALSE) %>%
  pull(V2) %>%
  sum()

genome_size

# -----------------------------
# Clean repeat class names
# -----------------------------
clean_repeat_class <- function(x) {
  x <- as.character(x)
  
  case_when(
    is.na(x) | x == "" | x %in% c("Unknown", "Unclassified", "NA") ~ "Unclassified",
    x %in% c("SINE?") ~ "SINE",
    x %in% c("RC") ~ "RC/Helitron",
    x %in% c("Simple_repeat", "Low_complexity") ~ "Simple/low-complexity repeats",
    x %in% c("RNA", "rRNA", "tRNA", "snRNA", "scRNA", "srpRNA") ~ "RNA repeats",
    TRUE ~ x
  )
}

# -----------------------------
# Read RepeatModeler consensus FASTA headers
# -----------------------------
headers <- readLines(consensus_fasta)
headers <- headers[grepl("^>", headers)]

consensus_df <- tibble(header = headers) %>%
  mutate(
    class_family = ifelse(
      grepl("#", header),
      sub(".*#([^[:space:]]+).*", "\\1", header),
      "Unclassified"
    ),
    TE_class_raw = sub("/.*", "", class_family),
    TE_class_clean = clean_repeat_class(TE_class_raw)
  )

consensus_summary_clean <- consensus_df %>%
  dplyr::group_by(TE_class_clean) %>%
  dplyr::summarise(
    n_consensus = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_consensus)) %>%
  mutate(
    label = comma(n_consensus)
  )

# Check total consensus number
sum(consensus_summary_clean$n_consensus)

consensus_summary_clean

# -----------------------------
# Read RepeatMasker .out file
# -----------------------------
rm <- fread(
  repeatmasker_out,
  skip = 3,
  fill = TRUE,
  header = FALSE,
  data.table = FALSE
)

rm <- rm[, 1:15]

colnames(rm) <- c(
  "SW_score", "perc_div", "perc_del", "perc_ins",
  "query_seq", "start", "end", "left",
  "strand", "repeat_name", "class_family",
  "repeat_start", "repeat_end", "repeat_left", "ID"
)

rm_clean <- rm %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    TE_class_raw = sub("/.*", "", class_family),
    TE_class_clean = clean_repeat_class(TE_class_raw)
  ) %>%
  filter(!is.na(start), !is.na(end), end >= start) %>%
  mutate(
    start0 = start - 1,
    end1 = end
  )

# -----------------------------
# Function to merge overlapping intervals
# -----------------------------
merge_intervals <- function(df) {
  df <- df %>%
    arrange(query_seq, start0, end1)
  
  result <- list()
  k <- 1
  
  for (chr in unique(df$query_seq)) {
    x <- df %>%
      filter(query_seq == chr) %>%
      arrange(start0, end1)
    
    if (nrow(x) == 0) next
    
    cur_start <- x$start0[1]
    cur_end <- x$end1[1]
    
    for (i in seq_len(nrow(x))) {
      s <- x$start0[i]
      e <- x$end1[i]
      
      if (s <= cur_end) {
        cur_end <- max(cur_end, e)
      } else {
        result[[k]] <- tibble(
          query_seq = chr,
          start0 = cur_start,
          end1 = cur_end
        )
        k <- k + 1
        cur_start <- s
        cur_end <- e
      }
    }
    
    result[[k]] <- tibble(
      query_seq = chr,
      start0 = cur_start,
      end1 = cur_end
    )
    k <- k + 1
  }
  
  bind_rows(result)
}

# -----------------------------
# Merge intervals separately for each cleaned repeat class
# -----------------------------
te_class_coverage_clean <- rm_clean %>%
  group_by(TE_class_clean) %>%
  group_modify(~ {
    merged <- merge_intervals(.x)
    tibble(
      merged_bp = sum(merged$end1 - merged$start0)
    )
  }) %>%
  ungroup() %>%
  mutate(
    genome_percent = merged_bp / genome_size * 100
  ) %>%
  arrange(desc(genome_percent)) %>%
  mutate(
    label = paste0(round(genome_percent, 2), "%")
  )

te_class_coverage_clean

# -----------------------------
# Total merged RepeatMasker coverage across all repeat classes
# -----------------------------
all_repeat_merged <- merge_intervals(rm_clean)

total_repeat_bp <- sum(all_repeat_merged$end1 - all_repeat_merged$start0)
total_repeat_percent <- total_repeat_bp / genome_size * 100

total_repeat_bp
total_repeat_percent

# -----------------------------
# Shared class order
# Prefer ordering by genome coverage, then add classes unique to Panel A
# -----------------------------
class_order <- te_class_coverage_clean %>%
  arrange(desc(genome_percent)) %>%
  pull(TE_class_clean)

extra_classes <- setdiff(
  consensus_summary_clean$TE_class_clean,
  class_order
)

class_order <- c(class_order, extra_classes)
consensus_summary_clean=consensus_summary_clean[order(consensus_summary_clean$n_consensus),]
consensus_summary_clean <- consensus_summary_clean %>%
  mutate(
    TE_class_clean = factor(TE_class_clean, levels = rev(class_order))
  )

te_class_coverage_clean=te_class_coverage_clean[order(te_class_coverage_clean$genome_percent),]
te_class_coverage_clean <- te_class_coverage_clean %>%
  mutate(
    TE_class_clean = factor(TE_class_clean, levels = rev(class_order))
  )

# -----------------------------
# Colors
# -----------------------------
te_colors <- c(
  "Other" = "#666666",
  "RNA repeats" = "#8C8C8C",
  "Simple/low-complexity repeats" = "#999999",
  "Satellite" = "#0072B2",
  "RC/Helitron" = "#D55E00",
  "SINE" = "#CC79A7",
  "LTR" = "#009E73",
  "LINE" = "#56B4E9",
  "DNA" = "#E69F00",
  "Unclassified" = "#BDBDBD"
)

# Add fallback colors for unexpected classes
all_classes <- union(
  as.character(consensus_summary_clean$TE_class_clean),
  as.character(te_class_coverage_clean$TE_class_clean)
)

missing_classes <- setdiff(all_classes, names(te_colors))

if (length(missing_classes) > 0) {
  extra_colors <- rep("#B3B3B3", length(missing_classes))
  names(extra_colors) <- missing_classes
  te_colors <- c(te_colors, extra_colors)
}


p1 <- ggplot(
  consensus_summary_clean,
  aes(x = TE_class_clean, y = n_consensus, fill = TE_class_clean)
) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(values = te_colors) +
  scale_y_continuous(
    labels = comma,
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "a",
    x = NULL,
    y = "Number of TE consensus sequences"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black", size = 13),
    axis.title = element_text(colour = "black", size = 15),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(face = "bold", size = 20)
  )

p1

p2 <- ggplot(
  te_class_coverage_clean,
  aes(x = TE_class_clean, y = genome_percent, fill = TE_class_clean)
) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = label),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(values = te_colors) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.20))
  ) +
  labs(
    title = "b",
    x = NULL,
    y = "Genome coverage by repeat class (%)"
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    axis.text = element_text(colour = "black", size = 13),
    axis.title = element_text(colour = "black", size = 15),
    axis.line = element_line(colour = "black"),
    plot.title = element_text(face = "bold", size = 20)
  )

p2

fig_s2 <- p1 + p2 + plot_layout(widths = c(1, 1))

fig_s2

ggsave(
  "Supplementary_Fig_S2_TE_annotation.png",
  fig_s2,
  width = 13,
  height = 6.5,
  dpi = 600
)
