#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  
#include <Rcpp.h>
using namespace Rcpp;

double bisquare_scalar(double h1, double h2, double delta1, double delta2, double r, double A) {
    double d1 = h1 - delta1;
    double d2 = h2 - delta2;
    double d = sqrt(pow(d1,2) + pow(d2,2));
    return (A * pow(1-pow(d/r,2),2)*(d < r));
}


// [[Rcpp::export]]
NumericVector bisquare_call(NumericVector h1, NumericVector h2, NumericVector delta, double r, double A) {
    int n = h1.size();
    NumericVector z(n);
    
    for (int i=0; i<n; i++) {
        z[i] = bisquare_scalar(h1[i],h2[i],delta[0],delta[1],r,A);
    }
    return z;
}
