#' @title A check function for NetID
#'
#' @description a check function for FateDynamic function.
#' FateDynamic have provided global_params, user could specific the parameters the want, the parameters are listed as follow
#'
#' @param min_counts
#' Default: 20
#'
#' @param nhvgs
#' default: 2000
#'
#' @param npcs
#' default: 30
#'
#' @param knn
#' default: 30
#' 
#' @param velo
#' if run the velocity model (cellrank), default: TRUE
#'
#' @params fate_predict
#' default: TRUE
#' 
#' @params cluster_label
#' specific the attributes in SingleCellExperiment object that contains the cluster label information
#' 
#' @export
#'
check_global_params <- function(params){
  # min_counts
  # nhvgs
  # npcs
  # knn
  # velo
  # fate_predict
  # cluster_label
  params_space <- list(
    min_counts = 20,
	nhvgs = 2000,
	npcs = 30,
	knn = 30,
	velo = TRUE,
	fate_predict = TRUE,
	cluster_label = NULL
  )
  if(length(intersect(names(params),names(params_space)))!=length(names(params))){
    stop("Invalid global parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}