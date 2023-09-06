#' @title A Function to Validate NetID Parameters
#'
#' @description This function validates parameters for the FateDynamic function.
#' FateDynamic provides palantir_params; users can specify the desired parameters, which are listed as follows:
#'
#' @param ndcs
#' The number of diffusion components. Default: 10.
#'
#' @param start_cell
#' The barcode or ID of the root cell. Default: NULL.
#'
#' @param terminal_state
#' The terminal cell type labels.
#'
#' @param nwps
#' The number for waypoint sampling. Default: 500.
#'
#' @param method
#' Specify a method to find the root cell. Three options include: "scent", "cytotrace", and "markergene".
#'
#' @param root_gene
#' If choosing "markergene", please specify the ID of the gene.
#'
#' @param plot
#' If generating the plots in the working directory.
#'
#' @export
#'
check_palantir_params <- function(params){
  params_space <- list(
    ndcs = 10,
    start_cell = NULL,
    terminal_state = NULL,
    nwps = 500,
    method = "scent",
    species = 10090,
    root_gene = NULL,
    plot = TRUE
  )
  if(length(intersect(names(params), names(params_space))) != length(names(params))){
    stop("Invalid palantir parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}
