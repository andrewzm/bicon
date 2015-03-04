bicon
=====

This page hosts a supplement package to the paper 'Multivariate Spatial
    Covariance Models: A Conditional Approach' by Cressie and Zammit-Mangion
    (2015, submitted). This package can be used for the construction of
    covariance matrices for bivariate processes modelled using the
    conditional approach. The covariance functions are assumed to be Matern
    (3/2) covariance functions whilst functionality is provided for
    constructing interaction matrices using the bisquare function as the
    interaction function. The accompanying vignettes may be used to replicate
    the studies found in the paper.

To download the vignettes without reproducing them on your machine, please view them in the `vignettes` folder or by directly click on the following links:

[Vignette 1 (Section 3.2)](https://github.com/andrewzm/bicon/blob/master/vignettes/bivariate_sim.pdf?raw=true)

[Vignette 2 (Section 5)](https://github.com/andrewzm/bicon/blob/master/vignettes/min_max_T.pdf?raw=true)
    
If you wish to reproduce the results, you will need to install the package and its dependencies.

To install the package please install `devtools` and then type

    library(devtools)
    install_github("andrewzm/bicon",build_vignettes=T,dependencies=T))
    
If all dependencies are installed on your machine you should be able to run and compile the vignettes. To view the vignettes please type

    library(bicon)
    vignette()
    
and select the vignettes under `bicon`.

References
-----

Cressie, N., \& Zammit-Mangion, A. (2015). Multivariate spatial covariance models: A conditional approach.
Submitted.