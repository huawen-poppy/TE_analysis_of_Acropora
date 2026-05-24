# TE_analysis_of_Acropora
The pipeline for TE analysis on the scRNA-seq data of Acropora hemprichii

# Below are the main analysis workflows
# Single cell RNA-seq data analysis
(1) cellranger analysis [preprocess the raw sequencing data to form the expression count tables]

(2) Seurat analysis [process the data to do clustering, cell type annotation, DEG identification]

(3) scenic analysis [GRN identify for each cell type]

# TE analysis
(1) Identify the TE [reference construction using the RepeatMasker, LongrepMasker]

(2) TE expression quantification from the sc data [use the pseudoaligment from kallisto]

(3) downstream analysis on TE and gene expression 
