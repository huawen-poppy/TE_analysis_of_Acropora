# TE_analysis_of_Acropora
The pipeline for TE analysis on the scRNA-seq data of Acropora hemprichii

Single cell RNA-seq data analysis
(1) cellranger analysis [preprocess the raw data]
(2) Seurat analysis [clustering, cell type annotation, DEG identification]
(3) scenic analysis [GRN identify for each cell type]

TE analysis
(1) Identify the TE [reference construction]
(2) TE expression quantification from the sc data
(3) downstream analysis on TE and gene expression 
