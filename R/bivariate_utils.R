###-----------------------------
### Bisquare functions
###-----------------------------

#' @title Bisquare function
#' @name bisquare_fns
#' @aliases bisquare_1d
#' @aliases bisquare_2d
#' @aliases bisquare_B
#' 
#' @description The bisquare function
#' @param h displacement (1d)
#' @param h1 first component of displacement vector (2d)
#' @param h2 second component of displacement vector (2d)
#' @param delta shift parameter(s)
#' @param r aperture parameter
#' @param A gain parameter (amplitude)
#' @param areas area associated with each column in B
#' @param n1 number of rows
#' @param n2 number of columns
#' @details The bisquare function (shifted by \eqn{\Delta}) is given by \deqn{b(s,v) \equiv \left\{\begin{array}{ll} A\{1 - (|v- s - \Delta|/r)^2\}^2, &| v -s  - \Delta| \le r \\ 0, & \textrm{otherwise}. \end{array} \right.}{b(s,v) =  A{1 - (|v- s - d|/r)^2}^2  if | v -s  - d| <= r,  0 otherwise}
#' The function \code{bisquare_1d} accepts \code{h} in any shape or size, while for \code{bisquare_2d}, \code{h1} and \code{h2} need to be vectors of the same size. The parameter \code{delta} needs to be equal to one or two in length, depending on the spatial dimension.
#' 
#' The function \code{bisquare_B} is used to construct the matrix \eqn{B} given the bisquare parameters. It is meant to be used in problems of 2 dimensions.
#' @export
#' @examples 
#' h <- seq(-10,10,0.1)
#' y <- bisquare_1d(h=h,delta=3,r=4)
#' plot(h,y)
bisquare_1d <- function(h,delta = 0,r = 1,A=1) { 
  y <- abs(h - delta)
  A*(1-(y/r)^2)^2 * (y < r)
}


#' @export
#' @rdname bisquare_fns
bisquare_2d <- function(h1,h2,delta = c(0,0),r = 1, A=1) {
    bisquare_call(h1, h2, delta, r, A)
}

#' @export
#' @rdname bisquare_fns
bisquare_B <- function(h1,
                       h2,
                       delta = c(0,0),
                       r = 1,
                       A = 1,
                       areas = rep(1,length(h1)),
                       n1 = 10L,
                       n2 = 10L) {
  this_z <- areas*bisquare_2d(h1,h2,delta = delta,r = r,A = A)
  matrix(this_z,n2,n1,byrow=T) #faster than doing transpose
}

### Bisquare computation in R (slow)
bisquare_notC <- function(h,delta,r,A) { 
    y <- t(t(h)-delta)
    y <- sqrt(y[,1]^2 + y[,2]^2)
    A*(1-(y/r)^2)^2 * (y < r)
}


###--------------------
### Matrix construction
###--------------------

#' @title Construct full (bivariate) covariance/precision matrix 
#' @name makeCov
#' @aliases makeSY
#' @aliases makeQY
#' 
#' @description Construct the covariance or precision matrix for the bivariate model constructed using the conditional approach.
#' @param r vector of distances
#' @param var variance of C
#' @param var1 variance of C_{11}
#' @param var2 variance of C_{2|1}
#' @param kappa scale of C
#' @param kappa1 scale of C_{11}
#' @param kappa2 scale of C_{2|1}
#' @param B interaction matrix
#' @return Covariance (or precision) matrix
#' @details  Both C_{11} and C_{2|1} are Matern covariance functions with smoothness parameter equal to 3/2. Covariance matrices are computed from Matern covariance functions using the vector of distances \code{r}, so that Sigma[1,1] = cov(Y1(s),Y1(s+r[1])), Sigma[1 + n,1] = cov(Y2(s),Y1(s+r[1])) and so on. Currently the grids on which Y1 and Y2 are evaluated need to be identical. 
#' 
#' The matrix \eqn{B} is the interaction matrix. The full covariance matrix returned is \deqn{\Sigma =   \left( \begin{array}{cc}\Sigma_{11} & \Sigma_{11}B' \\ B \Sigma_{11} & \Sigma_{2\mid 1}+B\Sigma_{11}B'  \end{array}\right).}{Sigma = [Sigma_{11} & Sigma_{11}B' ; B Sigma_{11} & Sigma_{2|1} + B Sigma_{11}B'].}
#' @export
#' @examples
#' s <- 0 : 99
#' D <- as.matrix(dist(s))
#' r <- as.vector(D)
#'
#' ## Assume the interaction matrix is the identity
#' B <- diag(100)
#' Sigma <- makeSY(r=r,var1=1,var2=1,kappa1=0.5,kappa2=0.1,B=B)
#' image(Sigma)
makeSY <- function(r, var1,var2,kappa1,kappa2,B) {
    S11 <- makeS(r,var1,kappa1)
    S2_1 <-  makeS(r,var2,kappa2)
    BS_21 <- B %*% S11
    rbind(cbind(S11,t(BS_21)),cbind(BS_21, tcrossprod(BS_21,B) + S2_1))
}

#' @rdname makeCov
#' @export
makeQY <- function(r,var1,var2,kappa1,kappa2,B) {
    Q11 <- makeQ(r,var1,kappa1)
    Q2_1 <-  makeQ(r,var2,kappa2)
    B <- B
    BQ_21 <- crossprod(B, Q2_1)
    rBind(cBind(crossprod(chol(Q2_1) %*% B) + Q11,-BQ_21),cBind(-t(BQ_21), Q2_1))
}

#' @rdname makeCov
#' @export
makeQ <- function(r,var,kappa) {
    chol2inv(chol(Matern32(r,var,kappa)))
}

#' @rdname makeCov
#' @export
makeS <- function(r,var,kappa) {
    Matern32(r,var,kappa)
}
###-----------------------------
### Matern functions
###-----------------------------


Matern32 <- function(r,var,kappa) { 
  X <- covMat3_call(r,var,kappa)
  matrix(X,nrow=sqrt(length(r)))
}



###--------------------
### Matrix tools
###--------------------

#' @title Create a sparse identity matrix
#'
#' @description Creates a sparse identity matrix of size n x n
#' @param n size of matrix
#' @export
#' @examples 
#' require(Matrix)
#' Q <- Imat(4)
Imat <- function (n) 
{
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}


#' @title Find the log determinant
#'
#' @description Find the log determinant of a matrix Q from its Cholesky factor L.
#' @param L the Cholesky factor of Q
#' @export
#' @examples 
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' logdet(chol(Q))
logdet <- function(L)  {
    ## Find the log-determinant of Q from its Cholesky L
    diagL <- diag(L)
    return(2*sum(log(diagL)))
}



### ------------------
### Distribution tuning
### ------------------

#' @title Quantile difference when tuning Gamma distribution hyper-parameters
#'
#' @description See details below.
#' 
#' @param pars the shape and rate parameters (in that order) of the Gamma distribution
#' @param p5 the 5-th percentile of the desired error distribution
#' @param p95 the 95-th percentile of the desired error distribution
#' @return A measure of discrepancy between the Gamma distribution with parameters \code{par} and the desired 5-th 95-th percentiles of the error signal.
#' @details Assume that in our model there is an error term (for example, observation error) \eqn{e \sim \mathcal{N}(0,\sigma^2)} and we are required to estimate the precision, \eqn{\sigma^{-2}}. Assume further that we have sufficient prior knowledge to know that \eqn{F_e(0.05)} = \code{p5} and \eqn{F_e(0.95)} = \code{p95}, where \code{F_e} is the cumulative distribution function of \eqn{e} and \code{p5} and \code{p95} are known. This function helps configure a prior distribution over \eqn{\sigma^{-2}} which reflects this prior belief. Specifically, let \eqn{\Lambda_{\alpha,\beta}} be the cumulative distribution function of the Gamma distribution with parameters \eqn{\alpha} and \eqn{\beta}. Then this function returns the discrepancy 
#' \deqn{\frac{\Lambda_{\alpha,\beta}(0.05) - 1/(p95)^2}{1/(p95)^2}  + \frac{\Lambda_{\alpha,\beta}(0.95) - 1/(p5)^2}{1/(p5)^2}}
#' 
#' For example, we can see how representative a Gamma distribution over the precision with shape and rate parameters in \code{pars} is of our prior belief that \code{p5} and \code{p95} are 1 and 5, respectively, by calling \code{Findalphabeta_gamma(pars,1,5)}. This function can hance be passed to an optimisation routine to find better values for \code{pars}.
#' @keywords Gamma distribution, prior elicitation
#' @export
#' @examples
#'
#' require(actuar)
#' # Find the Gamma distribution over the precision corresponding to the 
#' # prior belief of the error (1/sqrt(precision)) lying between p5=2, p95=5
#' initpars <- c(5,0.1)
#' hyp_pars <- optim(par=initpars,Findalphabeta_gamma, p5=1, p95=5)
#'
#' # Now simulate from a Gamma with these parameters and verify quantiles
#' X <- rgamma(shape = hyp_pars$par[1], rate = hyp_pars$par[2],n=10000)
#' print( quantile(1/sqrt(X),c(0.05,0.95)))
Findalphabeta_gamma <- function(pars,p5,p95) {
    if(any(pars<0)) {
        return(Inf)
    } else {
        return( sum((qgamma(c(0.05,0.95),shape=pars[1],rate=pars[2]) - c(1/p95^2,1/p5^2))^2/c(1/p95^2,1/p5^2)))
    }
}