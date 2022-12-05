#include <Rcpp.h>
using namespace Rcpp;

//' @title Gibbs sampling of bivariate normal variables using Rcpp
//' @description Gibbs sampling of bivariate normal variables using Rcpp
//' @param x0 Initial sampling value of x. 
//' @param y0 Initial sampling value of y. 
//' @param N Number of samples.
//' @return Sample results for binary variables.
// [[Rcpp::export]]
NumericMatrix RandomC(double x0,double y0, int N) {
  NumericMatrix out(N,2);
  out(0,0) = x0;
  out(0,1) = y0;
  for(int i=1;i < N;i++){
    NumericVector temp = rnorm(1,0.9*out(i-1,1),sqrt(0.19));
    out(i,0) = temp(0);
    NumericVector temp2 = rnorm(1,0.9*out(i,0),sqrt(0.19));
    out(i,1) = temp2(0);
  }
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//