#' @title A end-to-end fate dynamic inference
#'
#' @description A end-to-end cell dynamic inference. Run this function before RunNetID, if set dynamicInfor in RunNetID equal to TRUE.
#'
#' @param sce
#' A SingleCellExperiment object need to contain the "count" assay. "spliced" and "unspliced" assays are required if run cellrank
#'
#' @param global_params
#' a list object could manually set the core parameters of FateDynamic
#' see the details of parameters through ?check_global_params
#'
#' @param palantir_params
#' a list object could manually set the core parameters of palantir
#' see the details of parameters through ?check_palantir_params
#'
#' @param cellrank_params
#' a list object could manually set the core parameters of cellrank
#' see the details of parameters through ?check_cellrank_params
#'
#' @param work_dir
#' user could specific the working directory through input the address, and the varID object (pruneKnn) would be saved in this directory.
#' Otherwise the local working directory would be used (getwd())
#' Default: NULL
#'
#'
#' @export
#'
FateDynamic <- function(sce,global_params = list(),palantir_params = list(),cellrank_params = list(),work_dir = NULL){
  ## Build Annodata object using SingleCellExperiments
  if(!is.null(work_dir)){setwd(work_dir)}
  env = environment()

  suppressPackageStartupMessages(require("scran"))
  suppressPackageStartupMessages(require("irlba"))
  suppressPackageStartupMessages(require("rsvd"))
  suppressPackageStartupMessages(require("reticulate"))
  suppressPackageStartupMessages(require("Seurat"))
  suppressPackageStartupMessages(require("Hmisc"))

  anndata <- reticulate::import('anndata', convert = FALSE)
  sc <- reticulate::import('scanpy', convert = FALSE)
  scv <- reticulate::import('scvelo', convert = FALSE)

  sce <- sce[,colSums(assays(sce)$spliced)>0]
  X <- assays(sce)$spliced
  obs <- as.data.frame(colData(sce))
  var <- data.frame(genes = rownames(sce))
  X <- X[!duplicated(var$genes),]
  var <- data.frame(genes = rownames(X))
  rownames(var) <- rownames(X)
  obsm <- NULL
  reductions <- reducedDimNames(sce)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(reducedDim(sce, name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(reducedDimNames(sce)))
  }

  layers <- list()
  for (layer in c("spliced","unspliced")) {
    mat <- assays(sce)[[layer]]
    mat <- mat[rownames(X),]
    layers[[layer]] <- Matrix::t(mat)
  }

  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )

  adata$write("adata.h5ad")
  ####
  # input global parameters
  global_params <- check_global_params(global_params)
  list2env(global_params,env)
  ## check parameters
  palantir_params <- check_palantir_params(palantir_params)
  list2env(palantir_params,env)
  cellrank_params <- check_cellrank_params(cellrank_params)
  list2env(cellrank_params,env)

  # input parameters
  writeLines(paste("Using palantir to perform cell fate analysis...",sep=""))
  n.cores = 8
  if(is.null(palantir_params[["start_cell"]])){
    if(method == "cytotrace"){
      writeLines("Using CyToTRACE to identify start cell...")
      res = CytoTRACE::CytoTRACE(as.matrix(X),ncores = n.cores)$CytoTRACE
      res = as.matrix(res)
    }
    if(method == "scent"){
      writeLines("Using SCENT to identify start cell...")
      object <- CreateSeuratObject(X)
      object <- NormalizeData(object)
      X_norm <- GetAssayData(object)
      rm(object);gc()
      if(length(intersect(list.files(getwd()),"PPI.rds")) == 1){
        writeLines("Find PPI object at local dictionary, Read PPI...")
        PPI <- readRDS("PPI.rds")
      }else{
        PPI <- getPPI_String(X_norm,species=species)
      }
      res = SCENT::CompCCAT(X_norm, PPI)
      saveRDS(PPI,file="PPI.rds")
      res = as.matrix(res)
      rownames(res) <- colnames(X_norm)
      rm(X_norm,PPI);gc()
    }
    if(method == "markergene"){
      res = as.matrix(X[root_gene,])
    }
    start_cell = rownames(res)[which.max(res)]
  }

  input = "adata.h5ad"
  fate_predict_py = capitalize(tolower(as.character(fate_predict)))
  plot = capitalize(tolower(as.character(plot)))
  #if(is.null(cluster_label)) cluster_label = ""
  script = system.file("/python/FateDynamic_py.py", package = "NetID")
  commandLines = paste('python ',script,' -i ',input,' -m ',min_counts,' -n ',nhvgs,' -pc ',npcs,' -k ', knn, ' -pl ',"True",' -dc ',ndcs,' -r ',start_cell,' -nw ',nwps,' -mo ',mode,' -p ',plot,' -f ',fate_predict_py, ' -c ', cluster_label, ' -w ',weight_connectivities,' -nc ', ncores, ' -t ',tolerance,sep="")
  #reticulate::py_run_file(system.file(commandLines, package = "NetID"))
  system(commandLines, wait=TRUE)

  if(velo){
    writeLines(paste("Using cellrank with ",mode," mode velocity to perform cell fate analysis...",sep=""))
    script = system.file("/python/FateDynamic_py.py", package = "NetID")
    commandLines = paste('python ',script,' -i ',input,' -m ',min_counts,' -n ',nhvgs,' -pc ',npcs,' -k ', knn, ' -pl ',"False",' -dc ',ndcs,' -r ',start_cell,' -nw ',nwps,' -mo ',mode,' -p ',plot,' -f ',fate_predict_py, ' -c ', cluster_label, ' -w ',weight_connectivities,' -nc ', ncores, ' -t ',tolerance,sep="")
    #reticulate::py_run_file(system.file(commandLines, package = "NetID"))
    system(commandLines, wait=TRUE)
  }
}

getPPI_String <- function (object = NULL, species = 9606, score_threshold = 600,
                           save = FALSE)
{
  suppressPackageStartupMessages(require("data.table"))
  suppressPackageStartupMessages(require("igraph"))

  linkFiles <- paste("https://stringdb-static.org/download/protein.links.v11.0/",
                     species, ".protein.links.v11.0.txt.gz", sep = "")
  if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(linkFiles)))) {
    if (!file.exists(basename(linkFiles)))
      download.file(linkFiles, destfile = basename(linkFiles))
    gf <- gzfile(basename(linkFiles), "rt")
  }
  PPI <- read.table(gf, header = T, sep = "")
  PPI[,1] <- as.factor(PPI[,1])
  PPI[,2] <- as.factor(PPI[,2])
  close(gf)
  infoFiles <- paste("https://stringdb-static.org/download/protein.info.v11.0/",
                     species, ".protein.info.v11.0.txt.gz", sep = "")
  if (!file.exists(sub(pattern = ".gz", replacement = "", x = basename(infoFiles)))) {
    if (!file.exists(basename(infoFiles)))
      download.file(infoFiles, destfile = basename(infoFiles))
    gf <- gzfile(basename(infoFiles), "rt")
  }
  Pinfo <- read.table(gf, header = T, sep = "\t", colClasses = c("character",
                                                                 "character", "NULL", "NULL"), quote = "", row.names = 1)
  close(gf)
  PPI <- subset(PPI, combined_score > score_threshold)
  ENSP1 <- levels(PPI[, 1])
  levels(PPI[, 1]) <- toupper(Pinfo[ENSP1, ])
  ENSP2 <- levels(PPI[, 2])
  levels(PPI[, 2]) <- toupper(Pinfo[ENSP2, ])
  if (!is.null(object)) {
    gene_data <- rownames(object)
    gene_data_upper <- toupper(gene_data)
    gene_data <- as.data.frame(unique(as.data.table(data.frame(gene_data,
                                                               gene_data_upper)), by = "gene_data_upper"))
    rownames(gene_data) <- gene_data[, 2]
    PPI <- PPI[which(is.element(PPI[, 1], gene_data[, 2])),
    ]
    PPI <- PPI[which(is.element(PPI[, 2], gene_data[, 2])),
    ]
    levels(PPI[, 1]) <- gene_data[levels(PPI[, 1]), 1]
    levels(PPI[, 2]) <- gene_data[levels(PPI[, 2]), 1]
  }
  nodes <- union(PPI[, 1], PPI[, 2])
  links <- PPI[, 1:2]
  net <- graph_from_data_frame(d = links, vertices = nodes,
                               directed = FALSE)
  net <- igraph::simplify(net)
  if (save) {
    saveRDS(as_adj(net), paste(species, "_ppi_matrix_STRING-11.0.Rda",
                               sep = ""))
  }
  file.remove(paste(species, ".protein.links.v11.0.txt.gz",
                    sep = ""))
  file.remove(paste(species, ".protein.info.v11.0.txt.gz",
                    sep = ""))
  return(as_adj(net))
}



