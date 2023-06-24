#' @title A check function for NetID
#'
#' @description a check function for FateDynamic function.
#' FateDynamic have provided palantir_params, user could specific the parameters the want, the parameters are listed as follow
#'
#' @param ndcs
#' the number of diffusion components, default: 10
#'
#' @param start_cell
#' the barcode or ID of root cell, default: NULL
#'
#' @param terminal_state
#' the terminal cell type labels
#'
#' @param nwps
#' number for waypoint sampling, default: 500
#'
#' @param method
#' specific a method to find the root cell, three options includes: "scent","cytotrace" and "markergene"
#' 
#' @param root_gene
#' if choose "markergene", please specific the ID of the gene.
#'
#' @param plot
#' if generate the plots in working directory.
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
  if(length(intersect(names(params),names(params_space)))!=length(names(params))){
    stop("Invalid palantir parameter list. Please check the input!")
  }
  params_space[names(params)] <- params
  params_space
}