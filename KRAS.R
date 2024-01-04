rm(list = ls())
library(Seurat)
library(clustree)
library(reshape2)
library(DoubletFinder)
library(SeuratDisk)
library(symphony)
library(data.table)
library(ggplot2)
library(ggthemes)
library(ggrastr)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
#set.seed(0815)
scRNA <- readRDS(file = "./data/draw_data.RDS")
#calculate the percentage of mitochondrion in every cell
scRNA <- PercentageFeatureSet(object = scRNA, pattern = "^MT-", col.name = "percent.mt")
#QC
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
#Normalization,SCale,FindVariable
scRNA  = SCTransform(scRNA, verbose = TRUE)
#PCA
scRNA = RunPCA(scRNA, verbose = FALSE)
#Plot PCA elbow
ElbowPlot(scRNA, ndims = 50)
#cluster
scRNA = RunUMAP(scRNA, dims = 1:30, verbose = FALSE)
scRNA = RunTSNE(scRNA, dims = 1:30, verbose = FALSE)
#best resolution
scRNA_test = FindNeighbors(scRNA, verbose = FALSE)
scRNA_test = FindClusters(scRNA_test, resolution = c(seq(0,1.6,0.1)), verbose = FALSE)
clustree(scRNA_test, prefix = "SCT_snn_res.")
#RESOLUTION=1.1
scRNA = FindNeighbors(scRNA, verbose = FALSE)
scRNA = FindClusters(scRNA, resolution = 1.1, verbose = FALSE)


#doubleFinder
#doubletFinder best PC
a=list()
b=list()
c=list()
for (i in ident_list[1:4]) {
  sc_one <- subset(scRNA,orig.ident == i)
  scRNA.list <- paramSweep_v3(sc_one, PCs = 1:30, sct = T)
  scRNA.list.stats <- summarizeSweep(scRNA.list, GT = FALSE)
  bcmvn_scRNA <- find.pK(scRNA.list.stats) 
  mpK<-as.numeric(as.vector(bcmvn_scRNA$pK[which.max(bcmvn_scRNA$BCmetric)]))
  a[i]=mpK
  #Homotypic Doublet Proportion Estimate -------------------------------------
  annotations <-sc_one@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  DoubletRate = ncol(sc_one)*8*1e-6 
  b[i]=DoubletRate
  nExp_poi <- round(DoubletRate*length(sc_one$seurat_clusters)) 
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  #identify doublet cell
  sc_one <- doubletFinder_v3(sc_one, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = T)
  sc_one <- doubletFinder_v3(sc_one, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
  c[i]=sc_one
  
}

#merge doubleFinder
c[[1]]@meta.data[,"DF_hi.lo"] <- 'Singlet'
c[[1]]@meta.data$DF_hi.lo[which(c[[1]]@meta.data$DF.classifications_0.25_0.08_1509 == "Doublet")] <- "Doublet-High Confidience"
c[[1]]@meta.data$DF_hi.lo[which(c[[1]]@meta.data$DF.classifications_0.25_0.08_1509 == "Doublet" & c[[1]]@meta.data$DF.classifications_0.25_0.08_1420 == "Singlet")] <- "Doublet-Low Confidience"
c[[2]]@meta.data[,"DF_hi.lo"] <- 'Singlet'
c[[2]]@meta.data$DF_hi.lo[which(c[[2]]@meta.data$DF.classifications_0.25_0.07_4522 == "Doublet" & c[[2]]@meta.data$DF.classifications_0.25_0.07_4259 == "Singlet")] <- "Doublet-Low Confidience"
c[[2]]@meta.data$DF_hi.lo[which(c[[2]]@meta.data$DF.classifications_0.25_0.07_4522 == "Doublet")] <- "Doublet-High Confidience"
c[[3]]@meta.data[,"DF_hi.lo"] <- 'Singlet'
c[[3]]@meta.data$DF_hi.lo[which(c[[3]]@meta.data$DF.classifications_0.25_0.07_643 == "Doublet" & c[[3]]@meta.data$DF.classifications_0.25_0.07_600 == "Singlet")] <- "Doublet-Low Confidience"
c[[3]]@meta.data$DF_hi.lo[which(c[[3]]@meta.data$DF.classifications_0.25_0.07_643 == "Doublet")] <- "Doublet-High Confidience"
c[[4]]@meta.data[,"DF_hi.lo"] <- 'Singlet'
c[[4]]@meta.data$DF_hi.lo[which(c[[4]]@meta.data$DF.classifications_0.25_0.09_1012 == "Doublet" & c[[4]]@meta.data$DF.classifications_0.25_0.09_951 == "Singlet")] <- "Doublet-Low Confidience"
c[[4]]@meta.data$DF_hi.lo[which(c[[4]]@meta.data$DF.classifications_0.25_0.09_1012 == "Doublet")] <- "Doublet-High Confidience"

merged <- merge(c[[1]],y=c(c[[2]],c[[3]],c[[4]]))
scRNA@meta.data$DF_hi.lo <- merged@meta.data$DF_hi.lo

## plot doubleFinder
DimPlot(scRNA, reduction = "umap", group.by ="DF_hi.lo",cols =c("black","red","gold"))+
  theme(legend.position = "bottom")+
  DimPlot(scRNA, reduction = "umap", group.by ="seurat_clusters")

#drop double
scRNA_Singlet <- subset(scRNA,DF_hi.lo == 'Singlet')

#cell type annotation
##symphone
ref_10x = readRDS("./pbmcs_10x_reference.rds")
query = mapQuery(
  exp_query = scRNA_Singlet@assays$SCT@counts,
  metadata_query = scRNA_Singlet@meta.data,
  ref_obj = ref_10x,
  vars = c("orig.ident"),
  do_normalize = F, 
  do_umap = T 
)
query = knnPredict(query, ref_10x, ref_10x$meta_data$cell_type, k = 5)
scRNA_Singlet$Pred.cluster <- query$meta_data$cell_type_pred_knn

##cell type
scRNA_Singlet1 <- scRNA_Singlet
Idents(scRNA_Singlet1) <- "seurat_clusters"
scRNA_Singlet1 <- RenameIdents(scRNA_Singlet1,"0"="NK","20"="Monocyte",
                               "1"="Monocyte", "21"="Macrophages",
                               "2"="T_CD8", "22"="Monocyte",
                               "3"= "Neutrophil", "23"="MK",
                               "4"= "Neutrophil", "24"="Neutrophil",
                               "5"= "T_CD4","25"="T",
                               "6"= "Macrophages", "26"="MK",
                               "7"= "T_CD8", "27"="T_CD4",
                               "8"= "Monocyte","28"="RBC",
                               "9"= "NK" ,"29"="B",
                               "10"="Neutrophil","30"="Monocyte",
                               "11"="Monocyte", "31"="MK",
                               "12"="Monocyte", "32"="RBC",
                               "13"= "T_CD8", "33"="Macrophages",
                               "14"= "NK", "34"="B",
                               "15"= "T_CD8","35"="Macrophages",
                               "16"= "Monocyte", "36"="DC",
                               "17"= "T", "37"="Neutrophil",
                               "18"= "Monocyte","38"="Monocyte",
                               "19"= "T_CD8","39"="Macrophages","40"="T"
)
scRNA_Singlet@meta.data$Cell_Type <- Idents(scRNA_Singlet1)
saveRDS(scRNA_Singlet,"./pbmc/scRNA_Singlet_mt20_celltype.RDS")
