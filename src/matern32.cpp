#include <R.h> 		 /* required  */
#include <Rmath.h>  	 /* for distribution functions etc. */  
#include <Rcpp.h>
using namespace Rcpp;

double covMat3_scalar(double d, double var, double kappa) {
	return (var*(1+kappa*d)*exp(-d*kappa));
}

void covMat3(double *out, double *h, int *n, double *var, double *kappa) {		
	for (int i=0; i<*n; ++i) {		
		out[i] = covMat3_scalar(h[i],*var,*kappa);
	}
}

// [[Rcpp::export]]
NumericVector covMat3_call(NumericVector h, double var, double kappa) {
  int n = h.size();
  NumericVector out(n);
  
  for (int i=0; i<n; ++i) {    
    out[i] = covMat3_scalar(h[i],var,kappa);
  }
  return(out);
}