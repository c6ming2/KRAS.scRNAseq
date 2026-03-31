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
set.seed(0206)
NSCLC_colors <- c('B'='#66C2A5','DC'='#FC8D62','RBC'='#FF6347','MP'='#4e2da6','T_CD8'='#968763',
                  'Neu'='#8DA0CB','MK'='#E78AC3','T_CD4'='#D1C656','T'='#f2ec72','NK'='#62AAEA',
                  'Monocyte'='#A6D854')
patient_colors <- c('P1-baseline'='#11e38c','P1-cycle1'='#a1def0','P1-cycle2'='#d1bfec',
                    'P1-cycle3'='#ed2bb1','P2-baseline'='#115d52','P2-cycle1'='#239eb3')
T_colors <- c('T_CD8'='#968763','T_CD4'='#D1C656','T cell'='#f2ec72')
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
scRNA.Singlet1 <- RenameIdents(scRNA.Singlet1,"0"="NK","1"="Monocyte",
                               "2"="T", "3"= "T","4"= "MP","5"= "Neu",
                               "6"= "NK", "7"= "T", "8"= "Monocyte","9"= "Neu",
                               "10"="Monocyte","11"="T", "31"="MP",
                               "14"= "MK",  "16"= "Monocyte",  "20"= "B", "19"= "Neu",
                               "21"= "Neu","22"= "MK", "24"= "DC")
scRNA.Singlet@meta.data$Cell_Type <- Idents(scRNA.Singlet1)
saveRDS(scRNA_Singlet,".scRNA_Singlet_mt20_celltype.RDS")

##cellchat
library(CellChat)
scRNA_Singlet$Cell_Type <- droplevels(scRNA_Singlet$Cell_Type)
Seuratobject=createCellChat(scRNA_Singlet@assays$SCT@data,meta = Seuratobject@meta.data,group.by = 'Cell_Type')
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
Seuratobject@DB <- CellChatDB
Seuratobject <- subsetData(Seuratobject)
Seuratobject <- identifyOverExpressedGenes(Seuratobject) # 相当于suerat中的FindAllMarkers
Seuratobject <- identifyOverExpressedInteractions(Seuratobject)
Seuratobject <- projectData(Seuratobject, PPI.human) 
Seuratobject <- computeCommunProb(Seuratobject)  
Seuratobject <- filterCommunication( Seuratobject, min.cells = 10)
Seuratobject <- computeCommunProbPathway(Seuratobject)
cellchat <- aggregateNet(cellchat)


##Tcell
TCELL <- subset(scRNA_Singlet,Cell_Type %in% c('T'))
TCELL  = SCTransform(TCELL, verbose = TRUE)
TCELL= RunPCA(TCELL, verbose = FALSE)
TCELL=RunHarmony(TCELL,'orig.ident',plot_convergence=T)
TCELL = RunUMAP(TCELL, dims = 1:30, reduction = 'harmony',verbose = FALSE)
TCELL_test = FindNeighbors(TCELL, reduction = 'harmony',verbose = FALSE)
TCELL_test = FindClusters(TCELL_test, resolution = c(seq(0,1.6,0.1)), verbose = FALSE)
clustree(TCELL_test, prefix = "SCT_snn_res.")
#RESOLUTION=0.2
TCELL = FindNeighbors(TCELL, reduction = 'harmony',verbose = FALSE)
TCELL = FindClusters(TCELL, resolution = 0.5, verbose = FALSE)
DotPlot(TCELL,features = c('CD8A','CD8B','CD4'))
VlnPlot(TCELL,features = c('CD8A','CD8B','CD4'))
TCELL$T.type <- plyr::mapvalues(TCELL$seurat_clusters,
                                c(0,1,2,3,4,5,6,7,8,9),
                                c('CD8T','CD8T','CD4T','CD4T','CD8T','Other','CD4T','CD8T','CD8T','CD4T'))
CD8T <- subset(TCELL,T.type=='CD8T')
CD8T  = SCTransform(CD8T, verbose = TRUE)
CD8T= RunPCA(CD8T, verbose = FALSE)
CD8T=RunHarmony(CD8T,'orig.ident',plot_convergence=T)
CD8T = RunUMAP(CD8T, dims = 1:30, reduction = 'harmony',verbose = FALSE)
#寻找最佳resolution
CD8T_test = FindNeighbors(CD8T, reduction = 'harmony',verbose = FALSE)
CD8T_test = FindClusters(CD8T_test, resolution = c(seq(0,1.6,0.1)), verbose = FALSE)
clustree(CD8T_test, prefix = "SCT_snn_res.")
#RESOLUTION=0.2
CD8T = FindNeighbors(CD8T, reduction = 'harmony',verbose = FALSE)
CD8T = FindClusters(CD8T, resolution = 0.3, verbose = FALSE)
saveRDS(CD8T,'CD8T.RDS')

###monocle
library(monocle)
data = as(as.matrix(CD8T@assays$RNA@counts), 'sparseMatrix')
pd = new('AnnotatedDataFrame', data = CD8T@meta.data)
fData = data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd = new('AnnotatedDataFrame', data = fData)
monocle_cds = newCellDataSet(data,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5,
                             expressionFamily = negbinomial.size())  
monocle_cds = estimateSizeFactors(monocle_cds)
monocle_cds = estimateDispersions(monocle_cds) 
disp_table = dispersionTable(monocle_cds)
unsup_clustering_genes = subset(disp_table, mean_expression >= 0.1 & dispersion_empirical>=1 * dispersion_fit)
monocle_cds = setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)
monocle_cds = reduceDimension(monocle_cds, max_components = 2,method = 'DDRTree')
monocle_cds = orderCells(monocle_cds)
saveRDS(monocle_cds,'CD8T_monocle.rds')

###slingshot############################
library(slingshot)
Idents(CD8T) <- "seurat_clusters"
monocyte_slim <- as.SingleCellExperiment(CD8T)
sim <- slingshot(monocyte_slim,clusterLabels = 'seurat_clusters',reducedDim='UMAP',start.clus = '1')
lin1=slingshot(sim,clusterLabels ='seurat_clusters',reducedDim = 'UMAP')

##SCENIC
data(list='motifAnnotations_hgnc_v9',package = 'RcisTarget')
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases_old",nCores=1)
exprMat  <-  as.matrix(CD8T@assays$RNA@data)
cellInfo <-  CD8T@meta.data[,c(10,2,3)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
genesKept <- geneFiltering(exprMat, scenicOptions,
                           minCountsPerGene = 3 * 0.01 * ncol(exprMat), 
                           minSamples = ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases_old",nCores=1)
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,coexMethod=c("top5perTarget"))
scenicOptions <- initializeScenic(org="hgnc", dbDir="cisTarget_databases_old",nCores=1)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") 
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="scenicOptions_GENE4.rds")

##hdWGCNA
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(UCell)
theme_set(theme_cowplot())
CD8T@active.assay <- "RNA"
CD8T <- ScaleData(CD8T, features=VariableFeatures(CD8T))
CD8T <- SetupForWGCNA(CD8T,gene_select = "fraction",fraction = 0.05,wgcna_name = "CLUSTER6")
CD8T <- MetacellsByGroups( CD8T = CD8T,group.by = c("seurat_clusters", "orig.ident"),reduction = 'harmony', # select the dimensionality reduction to perform KNN on
  k = 25, max_shared = 10,  ident.group = 'seurat_clusters')
CD8T <- NormalizeMetacells(CD8T)
CD8T <- SetDatExpr(CD8T,group_name = c(6), group.by='seurat_clusters',assay = 'RNA',  slot = 'data')
CD8T <- TestSoftPowers(CD8T,networkType = 'signed' )
CD8T <- ConstructNetwork(CD8T,overwrite_tom = T,tom_name = '6')
PlotDendrogram(CD8T, main='cluster 6 hdWGCNA Dendrogram')
TOM <- GetTOM(CD8T)
CD8T <- ModuleEigengenes(CD8T,group.by.vars="orig.ident")
hMEs <- GetMEs(CD8T)
MEs <- GetMEs(CD8T, harmonized=FALSE)
CD8T <- ModuleConnectivity( CD8T,group.by = 'seurat_clusters', group_name = '6')
CD8T <- ResetModuleNames(CD8T,new_name = "6-M")
modules <- GetModules(CD8T) %>% subset(module != 'grey')
hub_df <- GetHubGenes(CD8T, n_hubs = 10)
CD8T <- ModuleExprScore(CD8T, n_genes = 25,method='UCell')
saveRDS(CD8T,'hdwGCnA.rds')













