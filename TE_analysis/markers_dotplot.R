library(Seurat)
adata=readRDS('Ahem_after_annotation.RDS')
DefaultAssay(adata)<-'integrated'
markers=FindAllMarkers(adata,only.pos = F)
markers=markers[markers$p_val_adj<0.05,]

library(dplyr)
top5_markers=markers %>%
  group_by(cluster) %>%
  slice_max(n=5,order_by=avg_log2FC)

anno=read.csv('Acropora_genesymbol.csv',header = T)
head(anno)
top5_markers$order=1:nrow(top5_markers)
top5_markers<-merge(top5_markers,anno,by.x = 'gene',by.y = 'Query',all.x = T)
unique(top5_markers$Hit.description)
top5_markers<-top5_markers[!duplicated(top5_markers),]

library(ggplot2)
top5_markers<-top5_markers[order(top5_markers$order),]

DotPlot(adata,features = top5_markers$gene)+
  RotatedAxis()+
  scale_x_discrete(labels=top5_markers$Hit.description)

te=readRDS('../../../AipSC_analysis_old_genome/TE_analysis/acropora/pure_te.rds')
Idents(te)<-te@meta.data$celltype
markers_te=FindAllMarkers(te,only.pos = T)
markers_te

library(Seurat)
library(dplyr)
library(ggplot2)

adata <- readRDS("Ahem_after_annotation.RDS")

# For marker detection, RNA assay is usually more appropriate than integrated assay
# If you really want to use integrated assay, change this back to "integrated"
DefaultAssay(adata) <- "integrated"

markers <- FindAllMarkers(
  adata,
  only.pos = TRUE
)

markers <- markers %>%
  filter(p_val_adj < 0.05)

top5_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(
    n = 5,
    order_by = avg_log2FC,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(order = row_number())

anno <- read.csv(
  "Acropora_genesymbol.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

# Keep one annotation per gene to avoid duplicated rows after joining
anno_unique <- anno %>%
  distinct(Query, .keep_all = TRUE)

top5_markers <- top5_markers %>%
  left_join(
    anno_unique,
    by = c("gene" = "Query")
  ) %>%
  arrange(order)

# Important: DotPlot cannot use duplicated genes
plot_markers <- top5_markers %>%
  filter(gene %in% rownames(adata)) %>%
  distinct(gene, .keep_all = TRUE)

# Create labels for x-axis
gene_labels <- plot_markers$Hit.description

# If annotation is missing, use gene ID instead
gene_labels[is.na(gene_labels) | gene_labels == ""] <- 
  plot_markers$gene[is.na(gene_labels) | gene_labels == ""]

names(gene_labels) <- plot_markers$gene


library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
adata=readRDS('Ahem_after_annotation.RDS')
DefaultAssay(adata)<-'integrated'
markers=FindAllMarkers(adata,only.pos = F)
markers=markers[markers$p_val_adj<0.05,]

library(dplyr)
top5_markers=markers %>%
  group_by(cluster) %>%
  slice_max(n=5,order_by=avg_log2FC)

anno=read.csv('Acropora_genesymbol.csv',header = T)
head(anno)
top5_markers$order=1:nrow(top5_markers)
top5_markers<-merge(top5_markers,anno,by.x = 'gene',by.y = 'Query',all.x = T)
unique(top5_markers$Hit.description)
top5_markers<-top5_markers[!duplicated(top5_markers),]

library(ggplot2)
top5_markers<-top5_markers[order(top5_markers$order),]
plot_markers <- top5_markers %>%
  filter(gene %in% rownames(adata)) %>%
  distinct(gene, .keep_all = TRUE)

# Create gene labels
gene_labels <- plot_markers$Hit.description

# If annotation is missing, use original gene ID
gene_labels[is.na(gene_labels) | gene_labels == ""] <- 
  plot_markers$gene[is.na(gene_labels) | gene_labels == ""]

# Optional: wrap long gene descriptions
gene_labels <- str_wrap(gene_labels, width = 18)

# Important: make labels named by gene ID
names(gene_labels) <- plot_markers$gene

p <- DotPlot(
  adata,
  features = plot_markers$gene,
  #dot.scale = 10,
  cols = c("lightgrey", "blue")
) +
  RotatedAxis() +
  scale_x_discrete(labels = gene_labels) +
  labs(
    x = "Marker genes",
    y = "Cell type"
  ) +
  theme_classic(base_size = 16) +
  theme(
    text = element_text(colour = "black"),
    
    axis.text.x = element_text(
      size = 16,
      colour = "black",
      angle = 60,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 22,
      colour = "black"
    ),
    
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    
    legend.title = element_text(
      size = 14,
      colour = "black",
    ),
    legend.text = element_text(
      size = 14,
      colour = "black"
    ),
    
    plot.margin = margin(5,5,5,5)
  )

p
ggsave(
  "marker_dotplot_publication.png",
  plot = p,
  width = 14,
  height = 8,
  dpi = 300
)