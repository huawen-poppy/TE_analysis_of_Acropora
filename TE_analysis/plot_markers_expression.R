library(ComplexHeatmap)
setwd('/Users/zhonh0b/Desktop/acropora_celltype')
adata=readRDS('pure_te.rds')
marker=read.csv('te_markers_only_positive_celltype_RNA_scaled.csv',header=T,row.names = 1)
head(marker)

library(dplyr)
library(ggplot2)
marker=read.csv('celltype_markers_all.csv')
marker$up_down=ifelse(marker$avg_log2FC>0,'up_regulated','down_regulated')
df_summary=aggregate(up_down ~ cluster, data = marker, FUN = function(x) c(up = sum(x == "up_regulated"), down = sum(x == "down_regulated")))
df_summary=(do.call(cbind, df_summary))
df_summary=as.data.frame(df_summary)
df_summary
df_summary_long <- reshape(df_summary, varying = list(2:3), direction = "long", sep = "_")
df_summary_long
df_summary_long$time=ifelse(df_summary_long$time==1,'up_regulated','down_regulated')
colnames(df_summary_long)<-c('Celltype','Category','value','id')
df_summary_long$Category=as.factor(df_summary_long$Category)
df_summary_long=df_summary_long[,c(1,2,3)]
df_summary_long$Category<-as.factor(df_summary_long$Category)
df_summary_long$Category=gsub('down_regulated','Down_regulated',df_summary_long$Category)
df_summary_long$Category=gsub('up_regulated','Up_regulated',df_summary_long$Category)
df_summary_long$value=as.numeric(df_summary_long$value)
ggplot(df_summary_long, aes(x = Celltype, y = value, fill = Category))+ 
     geom_bar(stat="identity") +
     #scale_fill_manual(values = c("up" = "red", "down" = "blue")) +
     #xlab("Cell Type") +
     ylab("Number of DEGs") +
     theme(axis.text.x = element_text(size=18,angle = 45,vjust = 1,hjust = 1))+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=15))

ggsave('./for_paper/fig1d_deg_number.png',dpi = 300)

library(ggplot2)

df_summary_long$Category=gsub('_','-',df_summary_long$Category)
#df_summary_long$Category=factor(df_summary_long$Category,levels = c('Up-regulated','Down-regulated'))
ggplot(df_summary_long, aes(x = Celltype, y = value, fill = Category)) + 
  geom_bar(
    stat = "identity",
    width = 0.75,
    colour = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "Up-regulated" = "lightcoral",
      "Down-regulated" = "lightblue3"
    )
  ) +
  labs(
    x = NULL,
    y = "Number of DEGs",
    fill = NULL
  ) +
  theme_classic(base_size = 22) +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(
      size = 20,
      colour = "black",
      angle = 45,
      vjust = 1,
      hjust = 1
    ),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.title.y = element_text(size = 22, colour = "black", face = "bold"),
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.text = element_text(size = 18, colour = "black"),
    legend.position = "right"
  )

ggsave(
  "./for_paper/fig1d_DEG_barplot_publication.png",
  width = 10,
  height = 6,
  dpi = 600
)

ahem=readRDS('Ahem_after_annotation.RDS')
Idents(ahem)<-ahem$Cell_type
ahem$Cell_type<-factor(ahem$Cell_type,level=c('Unknown','Gastrodermal cell','Calicoblast','Neurons','Nematocytes','Immune cell','Endosymbiotic cell','Gland cell'))
ahem$orig.ident<-gsub('A','Sample',ahem$orig.ident)
ahem@meta.data %>%
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
  ylab("Percentage")+
  theme(axis.text.x = element_text(size=18))+
  theme(axis.text.y = element_text(size=18))+
  theme(axis.title.y = element_text(size=18))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size=16))+
  theme(legend.text = element_text(size=15))


p <- ahem@meta.data %>%
  group_by(orig.ident, Cell_type) %>%
  count() %>%
  group_by(orig.ident) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = orig.ident, y = percent, fill = Cell_type)) +
  geom_col(
    width = 0.75,
    colour = "black",
    linewidth = 0.25
  ) +
  #scale_fill_manual(values = celltype_colors) +
  labs(
    title = "Percentage of cells per cell type",
    x = NULL,
    y = "Percentage (%)",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 20) +
  theme(
    # Make all text black
    text = element_text(colour = "black"),
    
    # Title
    plot.title = element_text(
      size = 24,
      face = "bold",
      colour = "black",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    
    # Axis text
    axis.text.x = element_text(
      size = 18,
      colour = "black",
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(
      size = 18,
      colour = "black"
    ),
    
    # Axis titles
    axis.title.x = element_text(
      size = 20,
      colour = "black",
      face = "bold",
      margin = margin(t = 10)
    ),
    axis.title.y = element_text(
      size = 20,
      colour = "black",
      face = "bold",
      margin = margin(r = 10)
    ),
    
    # Axis lines and ticks
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    
    # Legend
    legend.title = element_text(
      size = 18,
      colour = "black",
      face = "bold"
    ),
    legend.text = element_text(
      size = 16,
      colour = "black"
    ),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "right",
    
    # Clean background
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    # Margins
    plot.margin = margin(10, 15, 10, 10)
  )

p
ggsave('./for_paper/fig1c_celltype_proportion.png',  plot = p,
       width = 8,
       height = 6,
       dpi = 600)
