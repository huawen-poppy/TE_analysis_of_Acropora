library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

rm(list=ls())
setwd("/Acropora_sc_analysis/seurat_analysis/data/input_data")

############## setup the Seurat Object ###############
# Load the dataset
samples<-c('A1','A2','A3') #can also use paste("A",1:3,sep="")
samples<-as.list(samples)
pth="/Acropora_sc_analysis/seurat_analysis/data/input_data/"

samples_algae<-c(paste0("A",1:3,"_algae"))
samples_algae<-as.list(samples_algae)

# the count matrix is stored in A1[['RNA']]@counts
# we can also check the information using the A1@meta.data
# orig.ident is the samples' name
# ncount_RNA is each cells UMI count
# nFeature_RNA is the gene count detected in each cell
for (i in 1:3){
  names(samples)[i]<-samples[[i]]
  variable=paste(samples[[i]],".data",sep="")
  read_path=paste(pth,samples[[i]],sep = "")
  variable<-Read10X(data.dir = read_path)
  #assign(samples[[i]],CreateSeuratObject(counts = variable,project = samples[[i]],min.cells = 3,min.features = 200))
  samples[[i]]<-CreateSeuratObject(counts = variable,project = samples[[i]],min.cells = 0,min.features = 0)
}

for (i in 1:3){
  names(samples_algae)[i]<-samples_algae[[i]]
  variable=paste(samples_algae[[i]],".data",sep="")
  read_path=paste(pth,samples_algae[[i]],sep = "")
  variable<-Read10X(data.dir = read_path)
  #assign(samples[[i]],CreateSeuratObject(counts = variable,project = samples[[i]],min.cells = 3,min.features = 200))
  samples_algae[[i]]<-CreateSeuratObject(counts = variable,project = samples_algae[[i]],
                                         min.cells = 0,min.features = 0)
}



# use PercentageFeatureSet() to calculate the motichondrial QC matrix
# use MT- starting genes as a set of mitochondrial genes
for (i in 1:3){
  #samples[[i]][["percent.mt"]]<-PercentageFeatureSet(samples[[i]],pattern="^MT-")
  samples_algae[[i]][["percent.mt"]]<-PercentageFeatureSet(samples_algae[[i]],pattern="^MT-")
}

for (i in 1:3){
  #samples[[i]]$percent.mt[is.na(samples[[i]]$percent.mt)]<-0
  samples_algae[[i]]$percent.mt[is.na(samples_algae[[i]]$percent.mt)]<-0
}

#check the total UMI count for each barcode distributed
#check the total UMI count among all cells for each gene distributed
for (i in 1:3) {
  print(i)
  print(summary(colSums(samples[[i]][['RNA']]@counts)))
  print(summary(rowSums(samples[[i]][['RNA']]@counts)))
  print(summary(colSums(samples_algae[[i]][['RNA']]@counts)))
  print(summary(rowSums(samples_algae[[i]][['RNA']]@counts)))
}


# find the infection point

# the number of the unique genes and total molecules are automatically calculated during CreateSeuratObject()
# we can find them sorted in the object meta data
# here show the QC matrics for the first 5 cells
for (i in 1:3){
  print(head(samples[[i]]@meta.data,5))
  print(head(samples_algae[[i]]@meta.data,5))
}

################# FILTER & Normalization & Scale THE DATA ##################
filter_samples<-c('A1','A2','A3')
filter_samples<-as.list(filter_samples)
filter_samples_algae<-c(paste0("A",1:3,"_algae"))
filter_samples_algae<-as.list(filter_samples_algae)
samples_algae_backup<-samples_algae

for (i in 1:3){
  samples[[i]]<-subset(samples[[i]],subset = nCount_RNA >= 200)
  samples_algae[[i]]<-subset(samples_algae[[i]],subset = nCount_RNA > 0)
  index<-samples[[i]]@meta.data$nCount_RNA
  uUMI_high <- unname(quantile(index,0.99))
  samples[[i]]<-subset(samples[[i]],subset = nCount_RNA < uUMI_high)
  #samples[[i]]<-subset(samples[[i]],subset = percent.mt<0.002)
  filter_samples[[i]]<-NormalizeData(samples[[i]])
  filter_samples_algae[[i]]<-NormalizeData(samples_algae[[i]])
  filter_samples[[i]]<-ScaleData(filter_samples[[i]])
  filter_samples[[i]]<-FindVariableFeatures(filter_samples[[i]],selection.method = 'vst',nfeatures = 2000,verbose = F)
  filter_samples_algae[[i]]<-FindVariableFeatures(filter_samples_algae[[i]],
                                                  selection.method = 'vst',nfeatures = 2000,verbose = F)
  print(dim(filter_samples[[i]]@meta.data))
  print(dim(filter_samples_algae[[i]]@meta.data))
}

saveRDS(filter_samples,"../filter_samples.RDS")
saveRDS(filter_samples_algae,"../filter_samples_algae.RDS")

#this one is the one without scaling in each sample
