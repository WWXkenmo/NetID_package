## ----vignette-options, echo=FALSE, message=FALSE, warning=FALSE---------------
require(BiocStyle)

## ----load datasets, eval=TRUE, echo=TRUE, message=FALSE, warning=FALSE, message=FALSE----
library(NetID)
sce <- readRDS("blood_sce.rds")

## ----quick show example data, eval=TRUE, echo=TRUE, warning=FALSE, message=FALSE----
TF <- read.csv("mouse-tfs.csv",header=TRUE,stringsAsFactors=FALSE)
sce

## ----Validate cell type label, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
celltype = as.character(as.matrix(sce$celltype))
transform_to_underscore <- function(input_string) {
  transformed_string <- gsub(" ", "_", tolower(input_string))
  return(transformed_string)
}

# Apply the transformation function to the character sequence
transformed_sequence <- sapply(celltype, transform_to_underscore)
names(transformed_sequence) <- NULL
sce$celltype <- as.factor(transformed_sequence)

## ----Run NetID, eval=TRUE, echo=TRUE, warning=FALSE---------------------------
dyn.out <- RunNetID(sce, regulators = TF[,1], targets = TF[,1],
                    netID_params = list(normalize=FALSE), 
                    dynamicInfer = FALSE,no_cores = 4)

## ----Run NetID (set path), eval=FALSE, echo=TRUE, warning=FALSE,message=FALSE----
#  dyn.out <- RunNetID(sce, regulators = TF[,1], targets = TF[,1],
#                      netID_params = list(normalize=FALSE),
#                      dynamicInfer = FALSE, work_dir = "Set/the/Path/")

## ----Check NetID Parameters, eval=FALSE, echo=TRUE, warning=FALSE-------------
#  ?check_netID_params()

## ----Run NetID (geosketch results), eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
# install python "geosketch"
dyn.out <- RunNetID(sce, regulators = TF[,1], targets = TF[,1],
                    netID_params = list(sketch.method = "geosketch",
                                        normalize=FALSE), 
                    dynamicInfer = FALSE)

## ----return NetID outputs, eval=FALSE, echo=TRUE, warning=FALSE,message=FALSE----
#  names(dyn.out)

## ----Build Seurat object, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
library(SummarizedExperiment)
library(Seurat)
Se <- CreateSeuratObject(counts = assays(sce)$spliced)

Se <- NormalizeData(Se, normalization.method = "LogNormalize",
                    scale.factor = 10000)
Se <- FindVariableFeatures(Se, selection.method = "vst", 
                           nfeatures = 1000)

all.genes <- rownames(Se)
Se <- ScaleData(Se, features = all.genes)
Se <- RunPCA(Se, features = VariableFeatures(object = Se))
Se <- FindNeighbors(Se, dims = 1:10)
Se <- FindClusters(Se, resolution = 0.5)

Se@meta.data$celltype <- celltype
## or use 'seurat_clusters' as cluster_label in FateDynamic

## ----Run NetID with SNN, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
dyn.out <- RunNetID(Se, regulators = TF[,1], targets = TF[,1],
                    netID_params = list(normalize=FALSE, SNN = TRUE), 
                    dynamicInfer = FALSE)

## ----reset inout, eval=FALSE, echo=TRUE, warning=FALSE, message=FALSE---------
#  sce <- Se

## ----Run model, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE-------------
FateDynamic(sce,
            global_params = list(cluster_label = "celltype",
                                 min_counts = 10 ,
                                 nhvgs = 3000,
                                 velo = FALSE),
            palantir_params = list(start_cell = "TGATACGTGCTAGGAGCTT"))

## ----check FateDynamic parameters, eval=FALSE, echo=TRUE, warning=FALSE,message=FALSE----
#  ?check_global_params
#  ?check_palantir_params
#  ?check_cellrank_params

## ----Run FateDynamic (SCENT), eval=FALSE, echo=TRUE, warning=FALSE,message=FALSE----
#  FateDynamic(sce,
#              global_params = list(cluster_label = "celltype",
#                                   min_counts = 10 ,
#                                   nhvgs = 3000,
#                                   velo = FALSE),
#              palantir_params = list(method = "scent"))

## ----Run FateDynamic (set terminal state), eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
FateDynamic(sce,
            global_params = list(cluster_label = "celltype",
                                 min_counts = 10 ,
                                 nhvgs = 3000,
                                 velo = FALSE),
            palantir_params = 
              list(start_cell = "TGATACGTGCTAGGAGCTT",
                   terminal_state = c("erythroblasts_stage2",
                                      "neutrophils")))

## ----Run FateDynamic, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE-------
dyn.out <- RunNetID(sce,
                    regulators = TF[,1], 
                    targets = TF[,1],
                    netID_params =
                    list(normalize=FALSE,
                         sketch.method = "geosketch"), 
                    dynamicInfer = TRUE,
                    velo=FALSE)

## ----load fate probability and pseudotime, echo=FALSE, message=FALSE, warning=FALSE----
anndata <- reticulate::import("anndata", convert = FALSE)
ad = anndata$read_h5ad("./output/FateRes.h5ad")
fate_prob <- reticulate::py_to_r(ad$obs)
ID <- colnames(fate_prob)[-ncol(fate_prob)]
barcode <- rownames(fate_prob)
fate_prob <- as.matrix(fate_prob[, -ncol(fate_prob)])
rownames(fate_prob) <- barcode
colnames(fate_prob) <- ID

## ----show cell fate, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE--------
head(fate_prob)

## ----Plugin cell fate, eval=FALSE, echo=TRUE, warning=FALSE,message=FALSE-----
#  ## Run NetID
#  dyn.out <- RunNetID(sce,
#                      regulators = TF[,1],
#                      targets = TF[,1],
#                      netID_params =
#                      list(normalize=FALSE,
#                           sketch.method = "geosketch"),
#                      dynamicInfer = FALSE,
#                      velo=FALSE)
#  
#  ## Inject dynamic information
#  dyn.out$LineageClass <- LineageClassifier(fate_prob, maxState = 10, cut_off = 0)
#  dyn.out$fate_prob <- fate_prob # cell fate probability matrix

## ----plot cell fate probability, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
## load basis information
dyn.out$basis <- reducedDim(sce, "PCA")[,c(1,2)]

## For SCseq, e.g.
## dyn.out$basis <- sce@umap

## For Seurat, e.g.
## dyn.out$basis <- Se@reductions$umap_spliced@cell.embeddings

library(cowplot)
p1 = plotFateProb(dyn.out,basis=dyn.out$basis,basis_name = "PCA",
                  lineage = colnames(dyn.out$fate_prob)[1])
p2 = plotFateProb(dyn.out,basis=dyn.out$basis,basis_name = "PCA",
                  lineage = colnames(dyn.out$fate_prob)[2])
plot_grid(p1,p2)

## ----Run FateCausal, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE--------
GRN <- FateCausal(dyn.out,L=30)

## ----check network, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE---------
GRN$grn[[1]][1:10,1:10]

## ----statistical test of the module, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
result = ModuleFinder(GRN,"erythroblasts_stage2",c("Gata1"),
                      method = "spinglass")
result

## ----statistical test of the module2, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
result2 = ModuleFinder(GRN,"neutrophils",c("Gata1"),
                       method = "spinglass")
result2

## ----identify smallar module, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
result = ModuleFinder(GRN,"erythroblasts_stage2",c("Gata1"),
                      method = "spinglass",gamma = 0.4)
result

## ----visualize the module, eval=TRUE, echo=TRUE, warning=FALSE,message=FALSE----
graphPlot(GRN,result,"Gata1","erythroblasts_stage2")

## ----sessionInfo, echo=T------------------------------------------------------
sessionInfo()

