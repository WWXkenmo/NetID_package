#' @title A Function to Validate NetID Parameters
#'
#' @description This function validates parameters for the RunNetID function.
#' RunNetID provides netID_params; users can specify the desired parameters, which are listed as follows:
#'
#' @param g
#' Manually set the genes users want to analyze in their dataset. Default: NULL.
#'
#' @param var
#' If using variable genes to calculate principal components, used by the geosketch method. Default: FALSE.
#'
#' @param do.prune
#' If prune KNN graph.
#'
#' @param SNN
#' Use the Shared Nearest-Neighbor graph. Default: FALSE.
#'
#' @param sampled_cells
#' The barcode or ID of sampled cells.
#'
#' @param sketch.method
#' Perform sketching sampling on single-cell datasets: "random", "geosketch", or "SeuratSketching".
#'
#' @param ndim
#' Dimensions of PCs, used by the geosketch method. Default: 30.
#'
#' @param n_cell
#' The number of sampled cells. Default: 500.
#'
#' @param Threshold_Num
#' The minimum number of nearest neighbors for each seed cell after assignments.
#'
#' @param normalize
#' If performing normalization to the count matrix. Default: FALSE.
#'
#' @param prior_net
#' If a binary matrix indicates the prior knowledge of gene regulation, where rows are regulators and columns are targets.
#'
#' @export
#'
check_netID_params <- function(params){
  # var: If using variable genes to calculate principal components, used by geosketch method, default: FALSE
  # sampled_cells: The barcode of sampled cells
  # sketch.method: Perform sketching sampling on single-cell datasets: "geosketch" or "SeuratSketching"
  # ndim: Dimensions of PCs, used by geosketch method, default: 30
  # n_cell: The number of sampled cells, default: 500
  # Threshold_Num: The minimum number of nearest neighbors for each seed cell after assignments
  
  params_space <- list(
    g = NULL,
    var = FALSE,
    do.prune = TRUE,
    SNN = FALSE,
    sampled_cells = NULL,
    sketch.method = "random",
    ndim = 30,
    n_cell = 500,
    Threshold_Num = 2,
    normalize = FALSE,
    prior_net = NULL
  )
  if(length(intersect(names(params), names(params_space))) != length(names(params))){
    stop("Invalid global parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}
