#' @title A check function for NetID
#'
#' @description a check function for FateDynamic function.
#' FateDynamic have provided cellrank_params, user could specific the parameters the want, the parameters are listed as follow
#'
#' @param mode
#' the velocity mode, default: stochastic
#'
#' @param weight_connectivities
#' default: 0.2
#'
#' @param ncores
#' default: 6
#'
#' @param plot
#' default: TRUE
#' 
#' @param tolerance
#' default: 1e-6
#'
#' @export
#'
check_cellrank_params <- function(params){
  # mode
  # plot
  # weight_connectivities
  # ncores
  # tolerance
  params_space <- list(
    mode = "stochastic",
	weight_connectivities = 0.2,
	ncores = 6,
	plot = TRUE,
	tolerance = 10^-6
  )
  if(length(intersect(names(params),names(params_space)))!=length(names(params))){
    stop("Invalid cellrank parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}