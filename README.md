bicon
=====

This page hosts a supplement package to the paper 'Multivariate Spatial
    Covariance Models: A Conditional Approach' by Cressie and Zammit-Mangion
    (2015, in preparation). This package can be used for the construction of
    covariance matrices for bivariate processes modelled using the
    conditional approach. The covariance functions are assumed to be Matern
    (3/2) covariance functions whilst functionality is provided for
    constructing interaction matrices using the bisquare function as the
    interaction function. The accompanying vignettes may be used to replicate
    the studies found in the paper.
    
To install the package please install `devtools` and then type

    library(devtools)
    install_github("andrewzm/bicon",build_vignettes=T,dependencies=T))
    
If all dependencies are installed on your machine you should be able to run and compile the vignettes. To view the vignettes please type

    library(bicon)
    vignette()
    
and select the vignettes under `bicon`.

If you are unable to install bicon, you may still download the vignettes directly from this github page, from the `vignettes` folder.