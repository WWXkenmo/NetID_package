#' @title An End-to-End Fate Dynamic Inference
#'
#' @description This function provides an end-to-end cell dynamic inference process. Run this function before RunNetID, if you set dynamicInfor in RunNetID equal to TRUE.
#'
#' @param sce
#' A SingleCellExperiment object that needs to contain the "count" assay. "spliced" and "unspliced" assays are required if running cellrank.
#'
#' @param global_params
#' A list object that can manually set the core parameters of FateDynamic. Refer to the details of parameters using ?check_global_params.
#'
#' @param palantir_params
#' A list object that can manually set the core parameters of palantir. Refer to the details of parameters using ?check_palantir_params.
#'
#' @param cellrank_params
#' A list object that can manually set the core parameters of cellrank. Refer to the details of parameters using ?check_cellrank_params.
#'
#' @param work_dir
#' Users can specify the working directory by inputting the address, and the varID object (pruneKnn) will be saved in this directory. Otherwise, the local working directory will be used (getwd()).
#' Default: NULL.
#'
#'
#' @export
#'
FateDynamic <- function(sce,global_params = list(),palantir_params = list(),cellrank_params = list(),work_dir = NULL){
  ## Build Annodata object using SingleCellExperiments
  if(!is.null(work_dir)){setwd(work_dir)}
  env = environment()

  anndata <- reticulate::import('anndata', convert = FALSE)
  sc <- reticulate::import('scanpy', convert = FALSE)
  scv <- reticulate::import('scvelo', convert = FALSE)

  sce <- sce[,colSums(SummarizedExperiment::assays(sce)$spliced)>0]
  X <- SummarizedExperiment::assays(sce)$spliced
  obs <- as.data.frame(SummarizedExperiment::colData(sce))
  var <- data.frame(genes = rownames(sce))
  X <- X[!duplicated(var$genes),]
  var <- data.frame(genes = rownames(X))
  rownames(var) <- rownames(X)
  obsm <- NULL
  reductions <- SingleCellExperiment::reducedDimNames(sce)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(SingleCellExperiment::reducedDim(sce, name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(SingleCellExperiment::reducedDimNames(sce)))
  }

  layers <- list()
  for (layer in c("spliced","unspliced")) {
    mat <- SummarizedExperiment::assays(sce)[[layer]]
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
      if (require(CytoTRACE)) {
      	res = CytoTRACE::CytoTRACE(as.matrix(X),ncores = n.cores)$CytoTRACE
      	res = as.matrix(res)
      }
    }
    if(method == "scent"){
      writeLines("Using SCENT to identify start cell...")
      object <- Seurat::CreateSeuratObject(X)
      object <- Seurat::NormalizeData(object)
      X_norm <- Seurat::GetAssayData(object)
      rm(object);gc()
      if(length(intersect(list.files(getwd()),"PPI.rds")) == 1){
        writeLines("Find PPI object at local dictionary, Read PPI...")
        PPI <- readRDS("PPI.rds")
      }else{
        PPI <- getPPI_String(X_norm,species=species)
      }
      if (require(SCENT)) {
      	res = SCENT::CompCCAT(X_norm, PPI)
      }
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
  ## terminal_state
  terminalstate = terminal_state[1]
  if(is.null(terminalstate)){
	terminalstate = "None"
  }else{
  for(i in 2:length(terminal_state)){
	terminalstate = paste(terminalstate," ",terminal_state[i],sep="")
  }
  }

  input = "adata.h5ad"
  fate_predict_py = Hmisc::capitalize(tolower(as.character(fate_predict)))
  plot = Hmisc::capitalize(tolower(as.character(plot)))
  #if(is.null(cluster_label)) cluster_label = ""
  script = system.file("/python/FateDynamic_py.py", package = "NetID")
  commandLines = paste('python ',script,' -i ',input,' -m ',min_counts,' -n ',nhvgs,' -pc ',npcs,' -k ', knn, ' -pl ',"True",' -dc ',ndcs,' -r ',start_cell,' -ts ',terminalstate,' -nw ',nwps,' -mo ',mode,' -p ',plot,' -f ',fate_predict_py, ' -c ', cluster_label, ' -w ',weight_connectivities,' -nc ', ncores, ' -t ',tolerance,sep="")
  #reticulate::py_run_file(system.file(commandLines, package = "NetID"))
  system(commandLines, wait=TRUE)

  if(velo){
    writeLines(paste("Using cellrank with ",mode," mode velocity to perform cell fate analysis...",sep=""))
    script = system.file("/python/FateDynamic_py.py", package = "NetID")
    commandLines = paste('python ',script,' -i ',input,' -m ',min_counts,' -n ',nhvgs,' -pc ',npcs,' -k ', knn, ' -pl ',"False",' -dc ',ndcs,' -r ',start_cell, ' -ts ',terminalstate, ' -nw ',nwps,' -mo ',mode,' -p ',plot,' -f ',fate_predict_py, ' -c ', cluster_label, ' -w ',weight_connectivities,' -nc ', ncores, ' -t ',tolerance,sep="")
    #reticulate::py_run_file(system.file(commandLines, package = "NetID"))
    system(commandLines, wait=TRUE)
  }
}

getPPI_String <- function (object = NULL, species = 9606, score_threshold = 600,
                           save = FALSE)
{
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
    gene_data <- as.data.frame(unique(data.table::as.data.table(data.frame(gene_data,
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
  net <- igraph::graph_from_data_frame(d = links, vertices = nodes,
                               directed = FALSE)
  net <- igraph::simplify(net)
  if (save) {
    saveRDS(igraph::as_adj(net), paste(species, "_ppi_matrix_STRING-11.0.Rda",
                               sep = ""))
  }
  file.remove(paste(species, ".protein.links.v11.0.txt.gz",
                    sep = ""))
  file.remove(paste(species, ".protein.info.v11.0.txt.gz",
                    sep = ""))
  return(igraph::as_adj(net))
}



