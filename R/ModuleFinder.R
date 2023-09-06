#' @title Module Identification
#'
#' @description Identify modules from the gene regulatory network.
#'
#' @param fate_grn
#' A list object containing the lineage-specific gene regulatory network (output of FateCausal).
#'
#' @param lineage
#' Specify the name of the lineage.
#'
#' @param root_tf
#' The transcription factors for which the module needs to be identified.
#'
#' @param method
#' Module identification method: "rw" (random walk) or "spinglass" (spinglass).
#' Default: spinglass.
#'
#' @param steps
#' Tunable parameter for random walk. Description: the number of steps to take.
#' Default: 30.
#'
#' @param gamma
#' Tunable parameter for spinglass. Description: Real constant representing the gamma argument of the algorithm. This balances the importance of present and non-present edges in a community. A community is a set of vertices having many edges within the community and few edges outside. The default value of 1.0 makes existing and non-existing links equally important. Smaller values emphasize existing links, while greater values emphasize missing links.
#' Default: 0.5.
#'
#' @param sign
#' Whether to consider activation or inhibition regulation.
#'
#' @return A list containing the following objects:
#' @param ModuleList
#' A list of identified modules.
#'
#' @param ModuleTestResult
#' A data.frame containing the statistical test results of the modules.
#'
#' @examples
#' \dontrun{
#' result = ModuleFinder(fate_grn, "Ery", c("Gata1"), method = "spinglass")
#' }
#'
#' @export
#'
ModuleFinder <- function(fate_grn,lineage,root_tf,method = "spinglass",steps = 30,gamma = 0.5,sign = TRUE,output_type = "list"){
	if(sign){
	writeLines("Using signed network, the regulation type is determined by granger regression model")
	## using signed network
	GRN_weight = fate_grn$sign_grn
	}else{
	writeLines("Using unsigned network")
	GRN_weight = fate_grn$grn
	}
	
	idTypes <- c("rw","spinglass")
    msg <- paste0("should be one of ", paste(idTypes, collapse = ", "),
        ".")
    if (!method %in% idTypes) {
        stop("'method' ", msg)
    }
	
	idTypes2 <- c("graph","list")
    msg <- paste0("should be one of ", paste(idTypes2, collapse = ", "),
        ".")
    if (!output_type %in% idTypes2) {
        stop("'output_type' ", msg)
    }
	## extract the network
	grn = GRN_weight[[lineage]] # row is the target, column is the regulator
	## convert to sign network (activation, inhibition and no regulation)
	lineage_grn = abs(grn)
	
	## check if row and target gene is the same
	
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
	
	g = igraph::graph_from_adjacency_matrix(t(lineage_grn), mode = "directed", weighted = T)
	components <- igraph::decompose.graph(g)
	g <- components[[which.max(sapply(components, igraph::vcount))]]
	## filter out TFs
	root_tf = intersect(igraph::V(g)$name,root_tf)
	
	## perform random walk
	ntop <- length(root_tf)
	motif_list <- list()
	modN.v <- sizeN.v <- vector();
	for(start in root_tf){
		
		if(method == "rw"){
		community <- igraph::random_walk(
			  g,
			  start,
			  steps,
			  mode = "out",
		)}
		if(method == "spinglass"){
		community <- igraph::spinglass.community(g,spins=25,start.temp=1,stop.temp=0.1,cool.fact=0.99,update.rule=c("config"),gamma=gamma,vertex=start);
		}
		
		## extract the nodes within the module
		if(method == "spinglass"){community <- community$community}
		if(method == "rw"){community <- igraph::as_ids(community)}
		
		if(length(community)!=0){
		if(output_type == "list"){
		motif_list[[start]] = c(community)
		}
		community <- igraph::induced.subgraph(g,community);
		if(output_type == "graph"){
		motif_list[[start]] = community
		}
		modN.v = c(modN.v, mean(abs(igraph::get.edge.attribute(community,name="weight"))))
		if(output_type == "graph"){
		sizeN.v <- c(sizeN.v, length(igraph::V(community)$name))
		}
		if(output_type == "list"){
		sizeN.v <- c(sizeN.v, length(community))
		}
		}else{
		motif_list[[start]] = numeric(0)
		modN.v = c(modN.v, 0)
		sizeN.v <- c(sizeN.v,0)
		}
	}
	names(modN.v) <- root_tf;
	writeLines("Modularity values=");
	print(modN.v);
	
	## filter out empty module
	indicator = sizeN.v>0 & modN.v > 0
	modN.v = modN.v[indicator]
	motif_list = motif_list[indicator]
	root_tf = root_tf[indicator]
	sizeN.v = sizeN.v[indicator]
	ntop <- length(root_tf)
	
	if(sum(sizeN.v)==0){
	 stop("Apologies, it appears that the centralize module couldn't be found for the inputted TFs.")
	}else{
	## Perform statistical test to each module
	nMC = 1000
	writeLines("Starting Monte Carlo Runs");
	modNmc.m <- matrix(nrow=ntop,ncol=nMC);
	tmpEW.v <- igraph::get.edge.attribute(g,name="weight")
	for(m in 1:ntop){
	  nN <- sizeN.v[m];
	  if( (nN> 1)){
	  for(run in 1:nMC){
	   tmpEW.v.perm <- tmpEW.v[sample(1:length(tmpEW.v),length(tmpEW.v),replace=FALSE)];   
	   g.perm <- igraph::set.edge.attribute(g,"weight",value=tmpEW.v.perm)
	   if(output_type=="graph"){
	   subgrW.o <- igraph::induced.subgraph(g.perm,V(motif_list[[m]])$name);
	   }
	   if(output_type=="list"){
	   subgrW.o <- igraph::induced.subgraph(g.perm,motif_list[[m]]);
	   }
	   modNmc.m[m,run] <- mean(igraph::get.edge.attribute(subgrW.o,name="weight"));
	  }
	  }
	  writeLines(paste("Done for seed/module ",root_tf[m],sep=""));
	}
	
	modNpv.v <- rep(1,ntop);
	for(v in 1:ntop){
	  if( (sizeN.v[v] > 1) && (sizeN.v[v]< 1000)){
		modNpv.v[v] <- length(which(modNmc.m[v,] > modN.v[v]))/nMC;
	  }
	}
	names(modNpv.v) <- root_tf;
	
	ModuleList <- NULL
	for(i in 1:length(motif_list)){
		ModuleList[[i]] = c(igraph::V(g)$name[motif_list[[i]]])
	}
	names(ModuleList) <- names(motif_list)
	
	### Summarize the results and output
	tab <- data.frame(TF = root_tf, Module_size = sizeN.v, ModularityScore = modN.v, P_value = modNpv.v)
	
	res <- list()
	res$ModuleList = ModuleList
	res$ModuleTestResult = tab
	
	res	
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
