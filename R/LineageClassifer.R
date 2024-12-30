#' @title Classifies cells into lineages based on their fate probabilities using a clustering approach.
#'
#' @param fate_prob A matrix where rows are cells, columns are cell fates, and values are the probabilities of each cell being assigned to each fate.
#' @param cut_off Minimum threshold required for a cell to be assigned to a lineage. Default is 2.
#' @param maxState Maximum number of clusters to consider. Default is 5.
#' @param diffvar If TRUE, allows clusters to have different variances. Default is TRUE.
#' @param unique_assign If TRUE, assigns each cell to one specific lineage. If FALSE, allows overlap between lineages. Default is FALSE.
#'
#' @return A list where each element contains the cell IDs assigned to a specific lineage. Lineages with no cells are excluded.
#'
#' @details 
#' The function uses a Gaussian mixture model (via the `mclust` package) to group cells into clusters based on cell fate probabilities. 
#' It assigns cells to lineages based on the specified `cut_off` threshold and whether unique or shared assignments are allowed.
#'
#' @examples
#' # Example usage
#' fate_prob <- matrix(runif(100, 0.1, 0.9), nrow = 10, ncol = 10)
#' rownames(fate_prob) <- paste0("Cell", 1:10)
#' colnames(fate_prob) <- paste0("Fate", 1:10)
#' 
#' result <- LineageClassifier(fate_prob, cut_off = 2, maxState = 5, diffvar = TRUE, unique_assign = FALSE)
#'
#' @import mclust
#'
#' @export
#' 
LineageClassifer <- function(fate_prob,cut_off=2,maxState = 5, diffvar=TRUE, unique_assign = FALSE){
  sampleID <- rownames(fate_prob)
  cellfate <- colnames(fate_prob)
  #fate_prob <- fate_prob %*% diag(1/colMeans(fate_prob))
  #fate_prob <- diag(1/rowSums(fate_prob)) %*% fate_prob
  rownames(fate_prob) <- sampleID
  colnames(fate_prob) <- cellfate
  fateprob.v <- log2(1.000001+ fate_prob / (1.000001 - fate_prob))
  if(diffvar == TRUE){
    ## default assumes different variance for clusters
    mcl.o <- mclust::Mclust(fateprob.v, G = maxState)
  }
  else {
    mcl.o <- mclust::Mclust(fateprob.v, G = maxState, modelNames = c("E"))
  }
  mu.v <- mcl.o$param$mean
  for(i in 1:nrow(mu.v)){
    if(nrow(mu.v) == 2){
      mu.v[i,] <- mu.v[i,] / mcl.o$param$mean[-i,]
    }else{
      mu.v[i,] <- mu.v[i,] / colMeans(mcl.o$param$mean[-i,])
    }
  }
  colnames(mu.v) <- paste0("cluster",1:ncol(mu.v))
  class <- mcl.o$classification
  class <- paste0("cluster",class)
  print(mu.v)
  if(unique_assign){
    writeLines("Unique assign cell state into a specific lineage...")
    label <- NULL
    for(i in 1:ncol(mu.v)){
      if(max(mu.v[,i])>cut_off){
        label <- c(label, rownames(mu.v)[which.max(mu.v[,i])])
      }else{
        label <- c(label,"uncertain")
      }
    }
    lineage_list <- list()
    drop_fate <- NULL
    for(i in 1:nrow(mu.v)){
      lineage_list[[i]] <- colnames(mu.v)[which(label %in% c("uncertain",rownames(mu.v)[i]))]
      lineage_list[[i]] <- rownames(fate_prob)[which(class %in% lineage_list[[i]])]
      if(length(lineage_list[[i]])==0) drop_fate <- c(drop_fate,i)
    }
    names(lineage_list) <- rownames(mu.v)
    if(!is.null(drop_fate))lineage_list <- lineage_list[-drop_fate]
  }else{
    writeLines("Allow shared cell state between different lineage...")
    label_list <- list()
    for(i in 1:nrow(mu.v)){
      label_list[[i]] <- colnames(mu.v)[which(mu.v[i,]>cut_off)]
    }
    uncertain_cellstate <- colnames(mu.v)[colnames(mu.v) %in% unique(unlist(label_list)) == FALSE]
    for(i in 1:nrow(mu.v)) label_list[[i]] <- c(label_list[[i]],uncertain_cellstate)

    lineage_list <- list()
    drop_fate <- NULL
    for(i in 1:nrow(mu.v)){
      lineage_list[[i]] <- label_list[[i]]
      lineage_list[[i]] <- rownames(fate_prob)[which(class %in% lineage_list[[i]])]
      if(length(lineage_list[[i]])==0) drop_fate <- c(drop_fate,i)
    }
    names(lineage_list) <- rownames(mu.v)
    if(!is.null(drop_fate))lineage_list <- lineage_list[-drop_fate]
  }
  lineage_list
}