#' @rawNamespace import(Matrix, except = image)
#' @import extraDistr
#' @import Rcpp
#' @import graphics 
#' @importFrom splines ns 
#' @importFrom mclust adjustedRandIndex classError surfacePlot
#' @importFrom aricode NMI
#' @importFrom progressr handlers handler_progress with_progress progressor
#' @importFrom rlang new_environment duplicate
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession sequential
#' @importFrom grDevices gray rainbow heat.colors
#' @importFrom methods as
#' @importFrom utils menu packageVersion
#' @importFrom stringdist amatch
#' @importFrom igraph graph_from_adjacency_matrix as_edgelist
#' @importFrom scales alpha
NULL

