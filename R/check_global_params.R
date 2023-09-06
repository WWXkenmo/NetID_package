#' @title A Function to Validate NetID Parameters
#'
#' @description This function validates parameters for the FateDynamic function.
#' FateDynamic provides global_params; users can specify the desired parameters, which are listed as follows:
#'
#' @param min_counts
#' Default: 20.
#'
#' @param nhvgs
#' Default: 2000.
#'
#' @param npcs
#' Default: 30.
#'
#' @param knn
#' Default: 30.
#'
#' @param velo
#' If the velocity model (cellrank) is run, default: TRUE.
#'
#' @param fate_predict
#' Default: TRUE.
#'
#' @param cluster_label
#' Specify the attribute in the SingleCellExperiment object that contains cluster label information.
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
  if(length(intersect(names(params), names(params_space))) != length(names(params))){
    stop("Invalid global parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}
