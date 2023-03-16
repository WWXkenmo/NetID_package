#' @title Run NetID model
#'
#' @description This function is for generating the cell fate specific weighted gene regulatory network (GRN). Using the GRN skeleton and fate probability matrix, FateCausal would learn the granger causal using fate probability to order gene expression.
#'
#' @param dyn.out
#' The output of NetID.
#'
#' @param velo_infor
#' if use velocity information to infer gene regulatory network
#' Default: FALSE
#'
#' @param L
#' The number of lags to utilize
#' Default: 30
#'
#' @param alpha
#' use ridger (0) or lasso regression (1)
#' Default: 0
#'
#' @param lambda
#' regularization parameters
#' Default: 100
#'
#' @param cutoff
#' when use velocity aid GRN construction, FateCausal would convert the scaled velocity value > |cutoff| or < -|cutoff| into 1 and 0 to denote induction or repression. And used lagged gene expression to predict the induction or repression state.
#'
#' @param weight
#' the final GRN is defined as (1-w)*Granger_coefficient + w*GENIE3_coefficient, if aggre_method = "manual"
#' Default: 0.2
#'
#' @param redirected
#' if re-estimate the causal direction through granger causal test.
#' Default: FALSE
#'
#' @param fate_method
#' user need to specific which fate probability they prefer to use, "cellrank" or "palantir"
#' Default: "palantir"
#'
#' @param aggre_method
#' use which aggregation method, "RobustRankAggreg", "SUMMA" or "manual". If use "manual", need to specific weight parameter.
#'
#' @param restart
#' if re-estimate the GRN. Default: FALSE
#'
#' @param n.cores
#' The number of cores. Default: 8
#'
#' @param work_dir
#' user could specific the working directory through input the address, and the inferred GRN would be saved in this directory.
#' Otherwise the local working directory would be used (getwd())
#' Default: NULL
#'
#' @return a list contains following objects
#' @param GRN
#' a list of cell-fate specific GRN
#'
#' @param coef
#' a list of granger causal correlation matrix
#'
#' @examples
#' \dontrun{
#' GRN.filter <- FateCausal(dyn.out,velo_infor = FALSE,L=80,alpha = 0,weight=0.8,lambda = 100,aggre_method = "RobustRankAggreg",restart = TRUE)
#' }
#'
#' @export
#'
FateCausal <- function(dyn.out,velo_infor=FALSE,L=30,alpha = 0,lambda=100,cutoff=0.25,weight=0.2,redirected=FALSE,fate_method = "palantir",aggre_method = "RobustRankAggreg",restart=FALSE,n.cores = 8,work_dir = NULL){
  if(!is.null(work_dir)){setwd(work_dir)}

  suppressPackageStartupMessages(require("Matrix"))
  suppressPackageStartupMessages(require("lmtest"))
  suppressPackageStartupMessages(require("glmnet"))
  suppressPackageStartupMessages(require("summa"))
  suppressPackageStartupMessages(require("RobustRankAggreg"))

  if(velo_infor) velo_m <- dyn.out$velo_m
  #velo_m <- velo_m[apply(velo_m,1,function(x){sum(is.nan(x))})==0,]
  fate_prob <- as.matrix(dyn.out$fate_prob)
  if(velo_infor) fate_prob_velo <- as.matrix(dyn.out$fate_prob_velo)
  pseudotime <- dyn.out$pseudotime
  metaCells <- which(rownames(fate_prob) %in% colnames(dyn.out$skeleton$metaExp))
  metaExp <- dyn.out$skeleton$metaExp

  NN <- t(dyn.out$varID_res$NN)
  qval <- t(dyn.out$varID_res$pvM)

  ## Causal Inference
  ## Generate new directed skeleton through granger causal model
  if(redirected){
    dyn.out$skeleton$skeleton <- t(GrangerCausalSkeleton(t(dyn.out$skeleton$skeleton),metaExp,fate_prob,K=L,LineageClass,cutoff = cutoff))
  }

  if(velo_infor){
    genes <- intersect(rownames(velo_m),rownames(metaExp))
    writeLines(paste("Filter ",length(genes)," targets and ",nrow(metaExp)," regulators",sep=""))
    dyn.out$skeleton$skeleton <- t(t(dyn.out$skeleton$skeleton)[genes,])
    dyn.out$skeleton$GENIE3_net <- t(t(dyn.out$skeleton$GENIE3_net)[genes,])
  }
  fate_prob_aggre <- NULL
  for(i in 1:nrow(NN)){
    neighbor <- NN[i,]
    #neighbor <- c(neighbor[1],neighbor[2:length(neighbor)][qval[i,]>0.05])
    prob <- fate_prob[neighbor,]
    if(length(neighbor)!=1){
      prob <- colMeans(as.matrix(prob))
    }
    fate_prob_aggre <- rbind(fate_prob_aggre,prob)
  }
  rownames(fate_prob_aggre) <- rownames(fate_prob)

  if(velo_infor){
    fate_prob_velo_aggre <- NULL
    for(i in 1:nrow(NN)){
      neighbor <- NN[i,]
      #neighbor <- c(neighbor[1],neighbor[2:length(neighbor)][qval[i,]>0.05])
      prob <- fate_prob_velo[neighbor,]
      if(length(neighbor)!=1){
        prob <- colMeans(as.matrix(prob))
      }
      fate_prob_velo_aggre <- rbind(fate_prob_velo_aggre,prob)
    }
    rownames(fate_prob_velo_aggre) <- rownames(fate_prob_velo)
  }

  if(fate_method == "palantir"){
    LineageClass <- dyn.out$LineageClass
    fate_prob <- fate_prob_aggre
  }
  if(fate_method == "cellrank"){
    LineageClass <- dyn.out$LineageClass_velo
    fate_prob <- fate_prob_velo_aggre
  }

  if(fate_method %in% c("palantir","cellrank") == FALSE){
    stop("Invalid fate_method!")
  }


  GEP <- metaExp
  if(length(intersect(list.files(getwd()),"grn_list.rds")) == 1 && !restart){
    writeLines("Find grn_list object at local dictionary, Read grn_list...")
    grn_list <- readRDS("grn_list.rds")
  }else{
    grn_list <- list()
    fateid <- NULL
    for(i in 1:ncol(fate_prob)){
      if(length(LineageClass[[i]])>0){
        grn <- t(dyn.out$skeleton$skeleton)
        metaCells_fate = intersect(LineageClass[[i]],rownames(fate_prob)[metaCells])
        metaCells_fate = which(rownames(fate_prob) %in% metaCells_fate)
        GEP_sub <- GEP[,rownames(fate_prob)[metaCells_fate][order(fate_prob[metaCells_fate,i],decreasing=FALSE)]]
        cat(paste("Inferring GRN on Lineage ",i,",ncells = ",ncol(GEP_sub)," ,readout = Gene Expression..\n",sep=""))

        for(n in 1:nrow(GEP)){
          GEP_sub[n,] <- scale(GEP_sub[n,])
        }
        GEP_sub[is.nan(GEP_sub)] <- 0

        grn_list[[i]] <- do.call(rbind,parallel::mclapply(rownames(grn),function(j){
          pa <- c(j,colnames(grn)[grn[j,]!=0])
          if(length(pa) == 1){grn_vec <- grn[j,]}else{
            grn_vec = iLasso_raw(GEP_sub,V=NULL,L=L,grn_vec = grn[j,],alpha = alpha, lambda = lambda,pa,method="ridger")$grn_vec
            return(
              grn_vec
            )
          }},mc.cores=n.cores,mc.preschedule=TRUE))
        rownames(grn_list[[i]]) <- rownames(grn)
        fateid <- c(fateid, colnames(fate_prob)[i])
        saveRDS(grn_list,file="grn_list.rds")
      }
    }

    names(grn_list) <- fateid
  }

  ### build velo net flow
  if(velo_infor){
    if(length(intersect(list.files(getwd()),"grn_velo_list.rds")) == 1 && !restart){
      writeLines("Find grn_velo_list object at local dictionary, Read grn_velo_list...")
      grn_velo_list <- readRDS("grn_velo_list.rds")
    }else{
      genes <- intersect(rownames(metaExp),rownames(velo_m))
      V <- velo_m[genes,]
      grn_velo_list <- list()
      for(i in 1:ncol(fate_prob)){
        grn <- t(dyn.out$skeleton$skeleton)
        metaCells_fate = intersect(LineageClass[[i]],rownames(fate_prob)[metaCells])
        metaCells_fate = which(rownames(fate_prob) %in% metaCells_fate)
        GEP_sub <- GEP[,rownames(fate_prob)[metaCells_fate][order(fate_prob[metaCells_fate,i],decreasing=FALSE)]]
        cat(paste("Inferring GRN on Lineage ",i,",ncells = ",ncol(GEP_sub)," ,readout = RNA velocity...\n",sep=""))

        for(n in 1:nrow(GEP)){
          GEP_sub[n,] <- scale(GEP_sub[n,])
        }
        GEP_sub[is.nan(GEP_sub)] <- 0
        V_sub <- V[,colnames(GEP_sub)]
        V_sub[is.nan(V_sub)] <- 0
        grn_velo_list[[i]] <- do.call(rbind,parallel::mclapply(rownames(grn),function(j){
          pa <- c(j,colnames(grn)[grn[j,]!=0])
          if(length(pa) == 1){grn_vec <- grn[j,]}else{
            grn_vec = iLasso_raw(GEP_sub,V=V_sub,L=L,grn_vec = grn[j,],lambda = lambda,pa,method="ridger",velo_net=TRUE,cutoff=cutoff)$grn_vec
            return(
              grn_vec
            )
          }},mc.cores=n.cores,mc.preschedule=TRUE))
        rownames(grn_velo_list[[i]]) <- rownames(grn)
      }
      names(grn_velo_list) <- colnames(fate_prob)
      saveRDS(grn_velo_list,file="grn_velo_list.rds")
    }
  }

  ### Aggregate network inferred from gene expression and RNA velocity
  ## estimate weights
  extract <- function(x,M){x <- as.matrix(x);v <- M[x[1],x[2]];v}
  if(velo_infor){
    scaffoldNet <- as(as.matrix(t(dyn.out$skeleton$skeleton)[rownames(grn_list[[1]]),colnames(grn_list[[1]])]), "dgCMatrix")
    scaffoldNet <- Matrix::summary(scaffoldNet)
    prediction_m <- NULL
    for(i in 1:ncol(fate_prob)){
      prediction_m <- cbind(prediction_m,cbind(apply(scaffoldNet,MARGIN=1,extract,M=grn_list[[i]]),apply(scaffoldNet,MARGIN=1,extract,M=grn_velo_list[[i]])))
    }
    weight_vec_velo <- summa(prediction_m,"rank")@estimated_performance
    #weight_vec <- rep(1,length(weight_vec))
  }

  ### ordering GENIE3 network
  cat("Integrating Network Skeleton And Granger Causal Coefficience...\n")
  GENIE3_grn <- t(dyn.out$skeleton$GENIE3_net) * t(dyn.out$skeleton$skeleton)
  GENIE3_grn <- GENIE3_grn[rownames(grn_list[[1]]),colnames(grn_list[[1]])]
  #GENIE3_grn <- GENIE3_grn[genes,]
  GENIE3_score <- as.data.frame(Matrix::summary(as(GENIE3_grn, "sparseMatrix")))
  GENIE3_score$x <- rankScore(GENIE3_score$x)
  GENIE3_base <- sparseMatrix(
    i = GENIE3_score$i,
    j = GENIE3_score$j,
    x = GENIE3_score$x,dims = dim(GENIE3_grn)
  )
  rownames(GENIE3_base) <- rownames(GENIE3_grn)
  colnames(GENIE3_base) <- colnames(GENIE3_grn)

  ## ordering network inferred from
  grn_list_final <- list()
  for(i in 1:ncol(fate_prob)){
    ### generate the rank score
    grn <- grn_list[[i]]
    grn_score <- as.data.frame(Matrix::summary(as(grn, "sparseMatrix")))
    grn_score$x <- rankScore(grn_score$x)
    m1 <- sparseMatrix(
      i = grn_score$i,
      j = grn_score$j,
      x = grn_score$x,dims = dim(grn)
    )
    rownames(m1) <- rownames(grn)
    colnames(m1) <- colnames(grn)

    ### generate the rank score
    if(velo_infor){
      grn <- grn_velo_list[[i]]
      grn_score <- as.data.frame(Matrix::summary(as(grn, "sparseMatrix")))
      grn_score$x <- rankScore(grn_score$x)
      m2 <- sparseMatrix(
        i = grn_score$i,
        j = grn_score$j,
        x = grn_score$x,dims = dim(grn)
      )
      rownames(m2) <- rownames(grn)
      colnames(m2) <- colnames(grn)
    }


    if(velo_infor){
      weights_normalize <- c(weight_vec_velo[2*i-1],weight_vec_velo[2*i])
      weights_normalize <- weights_normalize / sum(weights_normalize)
      granger_net <- (weights_normalize[1]*m1 + weights_normalize[2]*m2)
    }else{
      granger_net <- m1
    }

    ### aggregate through SUMMA
    if(aggre_method == "SUMMA"){
      scaffoldNet <- as(as.matrix(t(dyn.out$skeleton$skeleton)[rownames(grn_list[[1]]),colnames(grn_list[[1]])]), "dgCMatrix")
      scaffoldNet <- Matrix::summary(scaffoldNet)
      prediction_m <- cbind(apply(scaffoldNet,MARGIN=1,extract,M=granger_net),apply(scaffoldNet,MARGIN=1,extract,M=GENIE3_base))
      weight_vec <- summa(prediction_m,"rank")@estimated_performance
      weight_vec <- weight_vec / sum(weight_vec)
      rm(prediction_m);gc()
      grn_list_final[[i]] <- granger_net*weight_vec[1]+GENIE3_base*weight_vec[2]
    }
    if(aggre_method == "RobustRankAggreg"){
      scaffoldNet <- as(as.matrix(t(dyn.out$skeleton$skeleton)[rownames(grn_list[[1]]),colnames(grn_list[[1]])]), "dgCMatrix")
      scaffoldNet <- Matrix::summary(scaffoldNet)
      prediction_m <- cbind(apply(scaffoldNet,MARGIN=1,extract,M=granger_net),apply(scaffoldNet,MARGIN=1,extract,M=GENIE3_base))
      linkList <- list(c(1:nrow(prediction_m))[order(prediction_m[,1],decreasing=TRUE)],
                       c(1:nrow(prediction_m))[order(prediction_m[,2],decreasing=TRUE)])
      linkList[[1]] <- as.character(linkList[[1]][which(prediction_m[linkList[[1]],1]!=0)])
      linkList[[2]] <- as.character(linkList[[2]][which(prediction_m[linkList[[2]],2]!=0)])

      r = rankMatrix(linkList)
      score <- aggregateRanks(rmat = r)
      scaffoldNet$x[as.integer(score$Name)] <- rankScore(1 / score$Score)
      aggreGRN <- sparseMatrix(
        i = scaffoldNet$i,
        j = scaffoldNet$j,
        x = scaffoldNet$x,dims = dim(GENIE3_base)
      )
      rownames(aggreGRN) <- rownames(GENIE3_base)
      colnames(aggreGRN) <- colnames(GENIE3_base)
      grn_list_final[[i]] <- aggreGRN
    }
    if(aggre_method == "manual"){
      grn_list_final[[i]] <- granger_net*(1-weight)+GENIE3_base*weight
    }
  }
  names(grn_list_final) <- names(grn_list)

  res <- NULL
  res$GRN <- grn_list_final
  res$coef <- grn_list
  res
}

countNonZero <- function(vec){
  vec <- vec - min(vec)
  sum(vec!=0)
}
rankScore <- function(score){
   score <- abs(score)
   score <- (1/(rank(-score))^2)*sign(score)
   score
}
scaleVar <- function(x){
   scale(x)
}

iLasso_raw <- function(X,V,L,pa,grn_vec,ntree=1000,alpha = 0, lambda = 100,method="XGBoost",velo_net = FALSE,cutoff=0.1){
  # Load time series
  # The column of X need to be ordered 1:t

  ttime = (L+1):ncol(X)

  if(!velo_net){
  tval = X[pa[1],ttime]
  }else{
	tval_pos <- tval_neg <- rep(0,length(ttime))
	tval_pos[scaleVar(V[pa[1],ttime])>cutoff] <- 1
	tval_neg[-scaleVar(V[pa[1],ttime])>cutoff] <- 1
  }
  regulators = pa[-1]
  numregs = length(regulators)

  idx = NULL
  for(i in regulators){
    idx = c(idx, rep(i,L))
  }
  idx = paste0(idx,rep(paste0("_",1:L,"L",sep=""),length(regulators)),sep="")
  #idx = paste0(rep(idx,length(ttime)),"_",rep(ttime,each=length(idx)),sep="")

  dat <- do.call(rbind,lapply(ttime,function(tstamp){
    X_sub = X[regulators,(tstamp-L):(tstamp-1)]
    X_sub = as.numeric(t(X_sub))
    #writeLines(paste("Done"," ", i,sep=""))
    return(
        X_sub
      )
  }))

  colnames(dat) = idx

  ### running lasso regression model
  if(!velo_net){count <- countNonZero(tval)}else{count <- max(sum(tval_pos),sum(tval_neg))}

  if(count>=3){
  #cvfit = cv.glmnet(dat, tval, family = "gaussian",alpha = alpha,type.measure = "mse", nfolds = 20)
  if(method == "ridger"){
  if(!velo_net){
  fit = glmnet(dat, tval, family = "gaussian",alpha = alpha,lamba = lambda)
  coef = coef(fit, s = lambda)[-1]
  }else{
    if(sum(tval_pos)>3) fit_pos = glmnet(dat, tval_pos, family = "binomial",alpha = alpha,lamba = lambda)
	if(sum(tval_neg)>3) fit_neg = glmnet(dat, tval_neg, family = "binomial",alpha = alpha,lamba = lambda)
	if(sum(tval_pos)>3){coef_pos = coef(fit_pos, s = lambda)[-1]}else{coef_pos = rep(0,L*numregs)}
	if(sum(tval_neg)>3){coef_neg = coef(fit_neg, s = lambda)[-1]}else{coef_neg = rep(0,L*numregs)}
	coef = (coef_pos + coef_neg)/2
  }
  }
  ### Using XGBoost regression

  }else{
  coef <- rep(0,L*numregs)
  }
  #
  coef.m = rowMeans(t(matrix(coef,nrow = L,ncol = numregs)))
  grn_vec[pa[-1]] = coef.m

  res <- list()
  res$coef.m <- coef.m
  res$coef <- coef
  res$grn_vec <- grn_vec
  #res$dat <- dat
  #res$tval <- tval
  res
}

GrangerCausalSkeleton <- function(skeleton,metaExp,fate_prob,K,LineageClass,cutoff=0.05,unique=FALSE,n.cores=5){
  ####
  # Skeleton: row is the targets, column is the regulators

  require("lmtest")
  require("Matrix")

  ### Build Undirected Skeleton
  #skeleton <- skeleton[,rownames(skeleton)] + t(skeleton[,rownames(skeleton)])
  skeleton <- t(skeleton)
  skeleton[skeleton!=0] <- 1
  Genes <- rownames(skeleton)

  sn <- function(x){names(x) <- x; x}
  phaseInfor <- do.call(c,parallel::mclapply(sn(Genes),function(gn){
  x<- metaExp[gn,]
   return(
     list(gep = x)
   )
  },mc.cores = 3,mc.preschedule = T))
  names(phaseInfor) <- Genes
  phaseInfor <- na.omit(phaseInfor)

  g_pair <- summary(as(as.matrix(skeleton),"dgCMatrix"))
  g_pair <- as.data.frame(g_pair)
  g_pair$i_name <- rownames(skeleton)[g_pair$i]
  g_pair$j_name <- colnames(skeleton)[g_pair$j]
  g_pair <- cbind(g_pair$i_name,g_pair$j_name)
  colnames(g_pair) <- c("Gene1","Gene2")
  g_pair <- as.data.frame(g_pair)
  genes <- unique(c(g_pair[,1],g_pair[,2]))

  p_tab_list <- skeleton_list <- list()

  for(ft in 1:ncol(fate_prob)){
  t = fate_prob[intersect(colnames(metaExp),LineageClass[[ft]]),ft]
  sn <- function(x){ names(x) <- x;x }
  idx = names(t)
  # To accelarte the computation, we extract the top quantile

  p_tab <- do.call(rbind,parallel::mclapply(sn(rownames(g_pair)),function(i){
	  pair <- as.matrix(g_pair[i,])
      GeneA <- pair[1]
      GeneB <- pair[2]
      ## Test gene A -> gene B
      #writeLines("Test causal from gene A to gene B")

	  a = phaseInfor[[GeneA]]
	  b = phaseInfor[[GeneB]]
	  a <- a[idx]
      b <- b[idx]
      a <- a[order(t,decreasing=FALSE)]
      b <- b[order(t,decreasing=FALSE)]
      if(sum(a!=0)<=5 || sum(b!=0) <=5){
	    p_a2b = p_b2a = 1
	  }else{
      p_a2b <- grangertest(a,b,K)$`Pr(>F)`[2]

      ## Test gene B -> gene A
      #writeLines("Test causal from gene B to gene A")
      p_b2a <- grangertest(b,a,K)$`Pr(>F)`[2]
      }


    dir = NULL
    if(unique){
    if(p_a2b<=p_b2a) dir = 1
    if(p_a2b>=p_b2a) dir = -1
    #if(min(c(p_a2b,p_b2a))>0.05) dir = 0
    pval = min(c(p_a2b,p_b2a))
    }else{
    if(p_a2b<=cutoff) dir = c(dir,1)
    if(p_a2b>cutoff) dir = c(dir,0)
    if(p_b2a<=cutoff) dir = c(dir,-1)
    if(p_b2a>cutoff) dir = c(dir,0)
    pval = c(p_a2b,p_b2a)
	pval = p_a2b
    }
    if(unique){
    res <- c(pair,dir,pval)
    }else{
    res <- rbind(c(pair,dir[1],pval[1]))
    }
    #writeLines(paste("Done"," ", i,sep=""))
    return(
        res
      )
  },mc.cores=n.cores,mc.preschedule=TRUE))

  p_tab <- as.data.frame(p_tab)
  colnames(p_tab) <- c("Gene1","Gene2","direction","p-value")
  p_tab$Gene1 <- as.character(as.matrix(p_tab$Gene1))
  p_tab$Gene2 <- as.character(as.matrix(p_tab$Gene2))
  p_tab$direction <- as.character(as.matrix(p_tab$direction))
  p_tab$`p-value` <- as.character(as.matrix(p_tab$`p-value`))
  p_tab_list[[ft]] <- p_tab

  ## Encode skeleton structure
  skeleton_res <- net_encode(genes,p_tab,sig=cutoff)
  skeleton_list[[ft]] <- skeleton_res[colnames(skeleton),rownames(skeleton)]
  }

  skeleton_list
}

net_encode <- function(genes, p_tab,sig){
  ## the row of returned network is target
  ## and the column of returned network is regulator
  NN_mat = matrix(0, nrow = length(genes), ncol = length(genes))
  rownames(NN_mat) <- genes
  colnames(NN_mat) <- genes
  p_tab$`p-value` <- p.adjust(p_tab$`p-value`,method="BH")
  p_tab <- p_tab[p_tab$`p-value`<sig,]
  for(i in 1:nrow(p_tab)){
    pair <- as.matrix(p_tab[i,c(1,2)])
    gene1 <- pair[1]
    gene2 <- pair[2]
    direction <- p_tab[i,"direction"]
    ##if direction = "-1", then gene2 -> gene1, gene2 is regulator
    ##if direction = "1", then gene1 -> gene2, gene1 is regulator
    if(direction == "-1"){
      NN_mat[gene1,gene2] <- -log10(as.numeric(p_tab[i,4]))
    }
    if(direction == "1"){
      NN_mat[gene2,gene1] <- -log10(as.numeric(p_tab[i,4]))
    }
  }
  NN_mat
}
