#include <Rcpp.h>
using namespace Rcpp;

//' @title Gauss-Seidel iterative method using Rcpp
//' @description Gauss-Seidel iterative method using Rcpp
//' @param A Matrix of coefficients.
//' @param b The right-hand side of a linear system. 
//' @param x0 Iterate the initial value. 
//' @param N Number of iterations.
//' @return Results of iterative solutions for linear systems.
//' @examples
//' \dontrun{
//' set.seed(1)
//' n <- 1e3
//' A <- sample(0:1,n^2,replace = TRUE,prob=c(0.95,0.05))
//' dim(A) <- c(n,n)
//' A <- -A + diag(rowSums(A)+1)
//' M_x <- rnorm(n)
//' b <- A %*% M_x
//' x0 <- array(rep(0,n))
//' x <- Gauss_Seidel(A,b,x0,1e3)
//' }
//' @export
// [[Rcpp::export]]
NumericVector Gauss_Seidel(NumericMatrix A,NumericVector b,NumericVector x0,int N=10) {
  int n = A.ncol();
  NumericVector x_old = x0 - 10;
  NumericVector x_new = x0;
  for(int i=1;i<N;i++){
    x_old = x_new;
    for(int j=0; j<n; j++)
    {
      double sum=0;
      for(int k=0; k<n; k++)
      {
        if(j>k)
        {
          sum+=(-1)*A(j,k)*x_new(k);
        }
        else if(j<k)
        {
          sum+=(-1)*A(j,k)*x_old(k);
        }
      }
      x_new[j]=(sum+b(j))/A(j,j);
    }
  }
  return x_new;
}