#'  @docType class
#'  @title Parent class for basis functions
#' @description A basis function contains three slots, the first is a list of parameters (\code{pars}), the second is the number (\code{n}) of basis functions and the third is a list of functions (\code{fn}). The pars field usually contains parameters which appear in fn.
#' @keywords Basis functions
#' @rdname Basisclass
setClass("Basis", representation(pars="list",n = "numeric",fn="list"))

#'  @docType class
#'  @title Finite element basis 
#' @description FEBasis inherits from the virtual class Basis. A finite element basis is initialised using the function \code{initFEbasis}.
#' 
#' @keywords Basis functions
#' @rdname FEBasisclass
setClass("FEBasis",contains="Basis")

