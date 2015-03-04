#' @title bicon
#' @docType package
#' @useDynLib bicon
#' @import Matrix
#' @import ggplot2
#' @import deldir
#' @import scales
#' @importFrom akima interp
#' @importFrom gpclib area.poly
#' @importFrom Rcpp sourceCpp
#' @import network
#' @name bicon
#' @description This is a supplement package to the paper 'Multivariate Spatial Covariance Models: A Conditional Approach' by Cressie and Zammit-Mangion (2015, submitted). This package can be used for the construction of covariance matrices for bivariate processes modelled using the conditional approach. The covariance functions are assumed to be Matern (3/2) covariance functions whilst functionality is provided for constructing interaction matrices using the bisquare function as the interaction function. The accompanying vignettes may be used to replicate the studies found in the paper.
NULL
