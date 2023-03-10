% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_palantir_params.R
\name{check_palantir_params}
\alias{check_palantir_params}
\title{A check function for NetID}
\usage{
check_palantir_params(params)
}
\arguments{
\item{ndcs}{the number of diffusion components, default: 10}

\item{start_cell}{the barcode or ID of root cell, default: NULL}

\item{nwps}{number for waypoint sampling, default: 500}

\item{method}{specific a method to find the root cell, three options includes: "scent","cytotrace" and "markergene"}

\item{root_gene}{if choose "markergene", please specific the ID of the gene.}

\item{plot}{if generate the plots in working directory.}
}
\description{
a check function for FateDynamic function.
FateDynamic have provided palantir_params, user could specific the parameters the want, the parameters are listed as follow
}
