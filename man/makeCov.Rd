% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bivariate_utils.R
\name{makeCov}
\alias{makeCov}
\alias{makeQ}
\alias{makeQY}
\alias{makeS}
\alias{makeSY}
\title{Construct full (bivariate) covariance/precision matrix}
\usage{
makeSY(r, var1, var2, kappa1, kappa2, B, nu1 = 3/2, nu2 = 3/2)

makeQY(r, var1, var2, kappa1, kappa2, B, nu1 = 3/2, nu2 = 3/2)

makeQ(r, var, kappa, nu)

makeS(r, var, kappa, nu)
}
\arguments{
\item{r}{vector of distances}

\item{var1}{variance of C_{11}}

\item{var2}{variance of C_{2|1}}

\item{kappa1}{scale of C_{11}}

\item{kappa2}{scale of C_{2|1}}

\item{B}{interaction matrix}

\item{var}{variance of C}

\item{kappa}{scale of C}
}
\value{
Covariance (or precision) matrix
}
\description{
Construct the covariance or precision matrix for the bivariate model constructed using the conditional approach.
}
\details{
Both C_{11} and C_{2|1} are Matern covariance functions with smoothness parameter equal to 3/2. Covariance matrices are computed from Matern covariance functions using the vector of distances \code{r}, so that Sigma[1,1] = cov(Y1(s),Y1(s+r[1])), Sigma[1 + n,1] = cov(Y2(s),Y1(s+r[1])) and so on. Currently the grids on which Y1 and Y2 are evaluated need to be identical.

The matrix \eqn{B} is the interaction matrix. The full covariance matrix returned is \deqn{\Sigma =   \left( \begin{array}{cc}\Sigma_{11} & \Sigma_{11}B' \\ B \Sigma_{11} & \Sigma_{2\mid 1}+B\Sigma_{11}B'  \end{array}\right).}{Sigma = [Sigma_{11} & Sigma_{11}B' ; B Sigma_{11} & Sigma_{2|1} + B Sigma_{11}B'].}
}
\examples{
s <- 0 : 99
D <- as.matrix(dist(s))
r <- as.vector(D)

## Assume the interaction matrix is the identity
B <- diag(100)
Sigma <- makeSY(r=r,var1=1,var2=1,kappa1=0.5,kappa2=0.1,B=B)
image(Sigma)
}

