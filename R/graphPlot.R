#' @title Plot Module Network
#'
#' @description Visualize the identified module.
#'
#' @param fate_grn
#' A list object containing the lineage-specific gene regulatory network (output of FateCausal).
#'
#' @param testRes
#' Module test result (output of the ModuleFinder).
#'
#' @param TF
#' The name of the transcription factor.
#'
#' @param lineage
#' The name of the lineage.
#'
#' @param sign
#' If visualizing the activation or inhibition regulation.
#'
#' @param cutoff
#' The minimum threshold of the network edges.
#'
#' @examples
#' \dontrun{
#' graphPlot(grn_list, result, "Klf1", "Ery")
#' }
#'
#' @export
#'
graphPlot <- function(fate_grn,testRes,TF,lineage,sign = TRUE,cutoff = 0){
# Extract corresponding module
# induce graph from the fate grn

	if(sign){
		lineage_grn = fate_grn$sign_grn[[lineage]]
	}else{
		lineage_grn = fate_grn$grn[[lineage]]
	}
	
	target = rownames(lineage_grn)
	regulator = colnames(lineage_grn)
	if(length(target) != length(regulator)){
	  writeLines("targets and regulators not the same length, create a square adjacency matrix")
	  lineage_grn = convertToSquareMatrix(lineage_grn)
	  target = rownames(lineage_grn)
	  regulator = colnames(lineage_grn)
	}
	if(sum(target == regulator) != length(regulator)){
	  writeLines("targets and regulators not the same")
	  lineage_grn = lineage_grn[target, target]
	}
	## extract module
	moduleName = testRes$ModuleList[[TF]]
	## induce graph
	adj = lineage_grn[moduleName,moduleName]
	dat = Matrix::summary(Matrix::Matrix(t(adj),sparse=TRUE))
	dat = as.data.frame(dat)
	dat$i = colnames(adj)[dat$i]
	dat$j = rownames(adj)[dat$j]
	dat$reg_type = rep(NA,nrow(dat))
	dat = dat[abs(dat$x)>cutoff,]
	dat$reg_type[dat$x > 0] = "activation"
	dat$reg_type[dat$x < 0] = "inhibition"
	
	colnames(dat) <- c("from","to","x","reg_type")
	dat <- dat[,-3]
	graph <- igraph::graph.data.frame(dat, directed = TRUE, vertices = data.frame(gene = unique(c(dat$from,dat$to)), gene_type = "Target"))
	
	# Customize appearance
	igraph::V(graph)$color <- ifelse(igraph::V(graph)$gene_type == "TF", "orange", "lightblue")
	igraph::V(graph)$size <- ifelse(igraph::V(graph)$gene_type == "TF", 30, 20)
	igraph::V(graph)$label <- igraph::V(graph)$name
	# Customize edge attributes and appearance
	igraph::E(graph)$color <- ifelse(igraph::E(graph)$reg_type == "activation", "green", "red")
	igraph::E(graph)$width <- ifelse(igraph::E(graph)$reg_type == "activation", 2, 1)

	# Choose a force-directed layout
	layout <- igraph::layout_with_fr(graph)

	# Plot and save the graph
	igraph::V(graph)$label.cex <- 0.8  # Reduce node ID size
	igraph::V(graph)$label.color <- "black"  # Change node ID color to black
	igraph::V(graph)$color[which(igraph::V(graph)$name == TF)] <- "red"
	plot(graph, layout = layout, edge.arrow.size = 0.5, main = paste(lineage," specific ", TF," regulatory network",sep=""))
	if(sign){
	legend("bottomright", legend = c("Activation", "Inhibition"), col = c("green", "red"), lwd = 1, x.intersp = 1)
	}
}

convertToSquareMatrix <- function(adj_matrix) {
    # Get the unique labels from rows and columns
    labels <- union(rownames(adj_matrix), colnames(adj_matrix))
    
    # Create an empty square matrix with dimensions equal to the number of unique labels
    n <- length(labels)
    square_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(labels, labels))
    
    # Copy values from the labeled adjacency matrix to the square matrix
    square_matrix[rownames(adj_matrix), colnames(adj_matrix)] <- adj_matrix
    
    return(square_matrix)
}

