#' @title Initialise a finite element basis
#' 
#' @description This function initialises an object of class \code{FEBasis} which defines a set of elemental `tent' basis functions over a pre-specified triangulation in 2-D
#'
#' @param p \code{n} \eqn{\times} 2 matrix of vertex locations.
#' @param t \code{m} \eqn{\times} 3 matrix of triangulations. Each row identifies which rows of \code{p} make up each triangle.
#' @param K \code{n} \eqn{\times} \code{n} connectivity matrix.
#' @return Object of class \code{FEBasis} 
#' @keywords finite elements, basis functions
#' @export
#' @examples
#' ## See vignette min_max_T
initFEbasis = function(p,t,K) {
    fn <- pars <- list()
    pars$p <- p
    pars$t <- t
    pars$K <- K
    df <- data.frame(x = pars$p[,1],
                     y = pars$p[,2],
                     n = 1:nrow(p))
    pars$vars <- df
    # Do tessellation
    Voronoi <- deldir(pars$p[,1],
                      pars$p[,2],
                      plotit='F',
                      sort=F,
                      rw=c(min(pars$p[,1])-0.00001,
                           max(pars$p[,1])+0.00001,
                           min(pars$p[,2])-0.00001,
                           max(pars$p[,2])+.00001))
    pars$pol <- PolygonfromVoronoi(Voronoi,pars$p)
    
    pars$vars$area_tess = rep(0,nrow(p))
    for (i in 1:nrow(p)) {
        pars$vars$area_tess[i] <- area.poly(pars$pol[[i]])
    }
    this_basis <- new("FEBasis", pars=pars, n=nrow(p), fn=fn)
    return(this_basis)
}



PolygonfromVoronoi <- function(Voronoi,p) {
    n = dim(p)[1]
    polygons <- vector("list",n)
    for (i in 1:n)  {
        X <- subset(Voronoi$dirsgs,(ind1 ==i | ind2 == i))
        X <- matrix(c(X$x1,X$x2,X$y1,X$y2,X$bp1,X$bp2),ncol=3)
        X <- unique(X)
        if(sum(X[,3])>0) {
            X <- rbind(X,c(p[i,],0))
        }
        edges <- X[,(1:2)]
        
        edges <- edges[chull(edges), ]
        polygons[[i]] <- as(edges,"gpc.poly")
    }
    return(polygons)
    
}


#' @title Assignment methods
#' @name [<--methods
#' @docType methods
#' @rdname assignment
#' @description Methods for \code{"[<-"}, i.e., extraction or subsetting of elements in the data frame of the mesh object
#' @examples 
#' ## See vignette min_max_T
NULL




#' @export
setMethod(f="[", signature="Basis", definition=function(x,i,j) {return(x@pars$vars[i][,])})

#' @rdname assignment
#' @aliases [<-,Basis,ANY,ANY-method
#' @export
setMethod(f="[<-", signature="Basis", definition=function(x,i,j,value) {
    x@pars$vars[i] <- value
    return(x)})

setGeneric("getDf", function(.Object) standardGeneric("getDf"))

setMethod("getDf",signature(.Object="FEBasis"),function(.Object) {
    return(.Object@pars$vars)
})

