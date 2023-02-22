#' @title A check function for NetID
#'
#' @description a check function for RunNetID function.
#' RunNetID have provided netID_params, user could specific the parameters the want, the parameters are listed as follow
#'
#' @param var
#' if using variable gene to calculate principal components, used by geosketch method, default: FALSE
#'
#' @param sampled_cells
#' the barcode or ID of sampled cells
#'
#' @param sketch.method
#' perform sketching sampling on single cell datasets, "geosketch" or "SeuratSketching"
#'
#' @param ndim
#' dimensions of PCs, used by geosketch method, Default: 30
#' 
#' @param n_cell
#' the number of sampled cells, default: 500
#'
#' @param Threshold_Num
#' the minimum nn of each seed cells after assignments, default: 2
#'
#' @param normalize
#' if perform normalization to the count matrix, default: FALSE
#' 
#' @param prior_net
#' if a binary matrix indicate the prior knowledge of gene regulation, row is the regulator, column is the target.
#' 
#' @export
#'
check_netID_params <- function(params){
  # var: if using variable gene to calculate principal components, used by geosketch method, default: FALSE
  # sampled_cells: the barcode of sampled cells
  # sketch.method: perform sketching sampling on single cell datasets, "geosketch" or "SeuratSketching"
  # ndim: dimensions of PCs, used by geosketch method, default: 30
  # n_cell: the number of sampled cells, default: 500
  # Threshold_Num: the minimum nn of each seed cells after assignments
  
  params_space <- list(
    var = FALSE,
	sampled_cells = NULL,
	sketch.method = "geosketch",
	ndim = 30,
	n_cell = 500,
	Threshold_Num = 2,
	normalize = FALSE,
	prior_net = NULL
  )
  if(length(intersect(names(params),names(params_space)))!=length(names(params))){
    stop("Invalid global parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}