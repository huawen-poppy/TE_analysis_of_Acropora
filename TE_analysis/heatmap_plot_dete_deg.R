library(Seurat)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
rm(list=ls())

setwd('./Desktop/acropora_celltype/')
te  <- readRDS("./pure_te.rds")
CC7 <- readRDS("./Ahem_after_annotation.RDS")

# -----------------------------
# Set identities
# -----------------------------
target  <- te
Idents(target) <- target$celltype

target2 <- CC7
Idents(target2) <- target2$Cell_type

DefaultAssay(target) <- "RNA"
te_markers <- FindAllMarkers(target, only.pos = T)
te_markers <- te_markers[te_markers$p_val_adj < 0.05, ]

# RNA gene markers
DefaultAssay(target2) <- "RNA"
ahem_markers <- FindAllMarkers(target2, only.pos = T)
ahem_markers <- ahem_markers[ahem_markers$p_val_adj < 0.05, ]

te_feats = te_markers$gene
gene_feats = ahem_markers$gene
target  <- NormalizeData(target, verbose = FALSE)
target2 <- NormalizeData(target2, verbose = FALSE)

# Average expression per cell type (features x celltypes)
avg_te <- AverageExpression(target, assays = DefaultAssay(target), features = te_feats, verbose = FALSE)[[DefaultAssay(target)]]
avg_ge <- AverageExpression(target2, assays = "RNA", features = gene_feats, verbose = FALSE)[["RNA"]]

# Match cell types (columns)
common_ct <- intersect(colnames(avg_te), colnames(avg_ge))
avg_te <- avg_te[, common_ct, drop=FALSE]
avg_ge <- avg_ge[, common_ct, drop=FALSE]

# Scale per feature
avg_te_s <- t(scale(t(avg_te)))
avg_ge_s <- t(scale(t(avg_ge)))

# Correlation across cell types: (genes x TE)
cor_mat_ct <- cor(t(avg_ge_s), t(avg_te_s), method = "pearson")


# Row annotation: DEG -> cell type
row_map <- ahem_markers %>%
  dplyr::select(gene, cluster) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

row_anno <- data.frame(
  Celltype = as.factor(row_map$cluster[match(rownames(cor_mat_ct), row_map$gene)])
)
rownames(row_anno) <- rownames(cor_mat_ct)

# Column annotation: DETE -> cell type
col_map <- te_markers %>%
  dplyr::select(gene, cluster) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

col_anno <- data.frame(
  Celltype = as.factor(col_map$cluster[match(colnames(cor_mat_ct), col_map$gene)])
)
rownames(col_anno) <- colnames(cor_mat_ct)

library(RColorBrewer)

all_ct <- union(levels(row_anno$Celltype), levels(col_anno$Celltype))
pal <- colorRampPalette(brewer.pal(8, "Set1"))(length(all_ct))
ann_colors <- list(Celltype = setNames(pal, all_ct))

# Order rows and columns by cell type
row_order <- order(row_anno$Celltype)
col_order <- order(col_anno$Celltype)

cor_mat_ct <- cor_mat_ct[row_order, col_order]
row_anno   <- row_anno[row_order, , drop=FALSE]
col_anno   <- col_anno[col_order, , drop=FALSE]

pheatmap(
  cor_mat_ct,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE
)

UMAPPlot(te,group.by='celltype')
DimPlot(te,group.by = 'celltype')
---------------------------------------------
# use the same color for cell types as in the umap plot
library(ggplot2)

# Make the UMAP plot object (don’t just print it)
p_umap <- DimPlot(te, group.by = "celltype")
p_build <- ggplot_build(p_umap)

colors_hex <- unique(p_build$data[[1]]$colour)
groups <- unique(p_build$data[[1]]$group)

# Map the internal group IDs back to the actual cluster names (identities)
# The mapping from group number to identity is typically stored here:
group_to_ident <- unique(p_build$plot$data[, c("celltype")])
# Ensure the mapping is correctly ordered
group_to_ident <- levels(group_to_ident)

# Create a named character vector of colors
# This vector has cluster names as names and hex codes as values
cluster_colors <- colors_hex
for (i in 1:length(groups)) {
  names(cluster_colors)[i] <- group_to_ident[groups[i]]
}

# Print the named vector of colors
print(cluster_colors)
names(cluster_colors)=c('Gland cell','Calicoblast','Unknown','Endosymbiotic cell',
                        'Immune cell','Nematocytes','Gastrodermal cell','Neurons')


all_ct <- union(levels(row_anno$Celltype), levels(col_anno$Celltype))

# Keep only the relevant celltypes, preserve UMAP colors
ann_colors <- list(Celltype = cluster_colors)
saaave=pheatmap(
  cor_mat_ct,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = row_anno,
  annotation_col = col_anno,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE
)
saaave
ggsave("all_deg_dete_heatmap.png", plot = saaave,dpi=300, device = "png", width =7, height = 8, units = "in")
