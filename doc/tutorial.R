## ----vignette-options, echo=FALSE, message=FALSE, warning=FALSE---------------
require(BiocStyle)

## ----load datasets, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE--------
library(NetID)
sce <- readRDS("blood_sce.rds")

## ----quick show example data, eval=TRUE, echo=TRUE, warning=FALSE-------------
TF <- read.csv("mouse-tfs.csv",header=TRUE,stringsAsFactors=FALSE)
sce

## ----Run NetID, eval=TRUE, echo=TRUE, warning=FALSE---------------------------
dyn.out <- RunNetID(sce, regulators = TF[,1], targets = TF[,1],netID_params = list(normalize=FALSE), dynamicInfer = FALSE)

## ----Build Seurat object, eval=TRUE, echo=TRUE, warning=FALSE-----------------
library(SummarizedExperiment)
library(Seurat)
Se <- CreateSeuratObject(counts = assays(sce)$spliced)

Se <- NormalizeData(Se, normalization.method = "LogNormalize", scale.factor = 10000)
Se <- FindVariableFeatures(Se, selection.method = "vst", nfeatures = 1000)

all.genes <- rownames(Se)
Se <- ScaleData(Se, features = all.genes)
Se <- RunPCA(Se, features = VariableFeatures(object = Se))
Se <- FindNeighbors(Se, dims = 1:10)
Se <- FindClusters(Se, resolution = 0.5)

## ----Run NetID with SNN, eval=TRUE, echo=TRUE, warning=FALSE------------------
#dyn.out <- RunNetID(Se, regulators = TF[,1], targets = TF[,1],netID_params = list(normalize=FALSE, SNN = TRUE), dynamicInfer = FALSE)

## ----Run model, eval=TRUE, echo=TRUE, warning=FALSE---------------------------
FateDynamic(sce,global_params = list(cluster_label = "celltype",min_counts = 10 ,nhvgs = 3000,velo = FALSE),palantir_params = list(start_cell = "TGATACGTGCTAGGAGCTT"))

## ----Run FateDynamic (SCENT), eval=TRUE, echo=TRUE, warning=FALSE-------------
#FateDynamic(sce,global_params = list(cluster_label = "celltype",min_counts = 10 ,nhvgs = 3000,velo = FALSE),palantir_params = list(method = "scent"))

## ----Run FateDynamic, eval=TRUE, echo=TRUE, warning=FALSE---------------------
dyn.out <- RunNetID(sce,regulators = TF[,1], targets = TF[,1],netID_params =
                    list(normalize=FALSE), dynamicInfer = TRUE,velo=FALSE)

## ----Run FateCausal, eval=TRUE, echo=TRUE, warning=FALSE----------------------
GRN <- FateCausal(dyn.out,L=30,aggre_method = "RobustRankAggreg")

## ----sessionInfo, echo=T------------------------------------------------------
sessionInfo()

