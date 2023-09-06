#' @title A function to check NetID
#'
#' @description This function checks the parameters for the FateDynamic function.
#' FateDynamic provides cellrank_params, allowing users to specify desired parameters. The available parameters are listed as follows:
#'
#' @param mode
#' The velocity mode. Default: stochastic.
#'
#' @param weight_connectivities
#' Default: 0.2.
#'
#' @param ncores
#' Default: 6.
#'
#' @param plot
#' Default: TRUE.
#'
#' @param tolerance
#' Default: 1e-6.
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
  if(length(intersect(names(params), names(params_space))) != length(names(params))){
    stop("Invalid cellrank parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}
