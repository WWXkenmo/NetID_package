#' @title Plot Cell Labels and Fate Probability on Dimensional Reduction Space
#'
#' @description A simple function to plot cell fate probability and cell type labels.
#'
#' @param dyn.out
#' The output of RunNetID.
#'
#' @param basis
#' A 2D array with (cell*2) dimensions, where columns represent coordinates.
#'
#' @param celltype
#' Cell type labels.
#'
#' @param basis_name
#' The name of the dimensional reduction space.
#'
#' @param lineage
#' Lineage name.
#'
#' @examples
#' \dontrun{
#' ## Plot cell types
#' plotFateProb(dyn.out, basis = dyn.out$basis, basis_name = "PCA", celltype = celltype)
#' ## Plot fate probability
#' plotFateProb(dyn.out, basis = dyn.out$basis, basis_name = "PCA", lineage = colnames(dyn.out$fate_prob)[1])
#' }
#'
#' @export
#'
plotFateProb <- function(dyn.out,basis,celltype = NULL,basis_name = "Dim",lineage = NULL){
	# Create a data frame
	
	if(is.null(celltype)&is.null(lineage)){
		stop("need give color information. cell type label or lineage label")
	}
	
	if(!is.null(lineage)){
	fate_prob <- dyn.out$fate_prob[,lineage]
	data <- data.frame(FateProb = fate_prob, UMAP1 = basis[, 1], UMAP2 = basis[, 2])

	# Define a color palette (adjust as needed)
	color_palette <- viridis::viridis(256, option = "A")

	# Create the UMAP plot using ggplot2	
	p1 <- ggplot2::ggplot(data, ggplot2::aes(x = UMAP1, y = UMAP2, color = FateProb)) +
	  ggplot2::geom_point(size = 3, alpha = 0.7) +
	  viridis::scale_color_viridis(option = "A") +
	  ggplot2::labs(title = paste(lineage," fate probability",sep=""),
		   x = paste(basis_name,"-1",sep=""), y = paste(basis_name,"-2",sep=""), color = "FateProb") +
	  ggplot2::theme_minimal() +
	  ggplot2::theme(legend.position = "right")
	}
	if(!is.null(celltype)){
	data <- data.frame(celltype = celltype, UMAP1 = basis[, 1], UMAP2 = basis[, 2])

	# Define a color palette (adjust as needed)
	color_palette <- viridis::viridis(length(table(celltype)), option = "A")

	# Create the UMAP plot using ggplot2
	p1 <- ggplot2::ggplot(data, ggplot2::aes(x = UMAP1, y = UMAP2, color = celltype)) +
	  ggplot2::geom_point(size = 3, alpha = 0.7) +
	  ggplot2::scale_color_manual(values = color_palette) +
	  ggplot2::labs(title = "Cell Type",
		   x = paste(basis_name,"-1",sep=""), y = paste(basis_name,"-2",sep=""), color = "Cell Type Label") +
	  ggplot2::theme_minimal() +
	  ggplot2::theme(legend.position = "right")
	}
	p1
}
