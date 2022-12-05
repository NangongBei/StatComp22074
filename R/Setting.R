#' @title A illustration dataset
#' @name data
#' @description A dataset used to illustrate the performance of \code{vaccR} and \code{vaccC}.
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' tm <- microbenchmark::microbenchmark(
#'   vR = vaccR(age,female,ily),
#'   vC = vaccC(age,female,ily)
#' )
#' print(summary(tm)[,c(1,3,5,6)])
#' }
NULL

#' @title Linear system with n=100
#' @name matrix_data_1e2
#' @description A dataset of linear system with n=100.
#' 
NULL

#' @title Linear system with n=500
#' @name matrix_data_500
#' @description A dataset of linear system with n=500.
#' 
NULL

#' @title Linear system with n=1000
#' @name matrix_data_1e3
#' @description A dataset of linear system with n=1000.
#' 
NULL

#' @title HW1_work function
#' @name benchmarks
#' @description Use R package \code{microbenchmark} to compare the performance of C functions (\code{gibbsR} and \code{vaccR}) and Cpp functions (\code{gibbsC} and \code{vaccC}).
#' @examples
#' \dontrun{
#' x <- rdistPareto(5,2,2,FALSE)
#' print(x)
#' }
#' @import microbenchmark knitr rmarkdown boot bootstrap DAAG lmtest mediation pander xtable
#' @importFrom Rcpp evalCpp
#' @importFrom stats rexp runif rgamma rnorm lm var
#' @importFrom graphics hist lines abline legend
#' @useDynLib StatComp22074
NULL

#' Return random numbers of Pareto(a,b) distributes
#'
#' @param n random sample size
#' @param a parameters a
#' @param b parameters b
#' @param picture plot histogram
#'
#' @return random numbers of Pareto(a,b) distributes 
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rdistPareto(5,2,2,FALSE)
#' print(x)
#' }
#' 
rdistPareto = function(n = 1e3,a = 2,b = 2,picture = FALSE){
  u <- runif(n) # Generate uniformly distributed random numbers
  x <- b/(u)^(1/a) # calculate F^(-1)(U')
  if(picture){
    # Draw a histogram
    hist(x, prob = TRUE, breaks = 'scott',main = expression(f(x)==a*b^a/x^(a+1),x>=b))
    # Sampling and plotting probability density curves
    sampling <- seq(b,max(x),.1)
    lines(sampling,a*b^a/sampling^(a+1))
  }
  return(x)
}

#' Return random numbers of Beta(alpha,beta) distributes
#'
#' @param n random sample size
#' @param alpha parameters alpha
#' @param beta parameters beta
#' @param picture plot histogram
#'
#' @return random numbers of Beta(a,b) distributes
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rdistBETA(5,2,2,FALSE)
#' print(x)
#' }
#' 
rdistBETA = function(n = 1e3,alpha = 3,beta = 2,picture = FALSE){
  sample_total <- m <- 0 # draw the mth sample and Output sampling total
  y <- numeric(n) # Record the generated random numbers
  while(m < n)
  {
    u <- runif(1) # Generate random numbers from U(0,1)
    sample_total <- sample_total + 1
    x <- runif(1) # Generate random numbers from g(x)
    if(u <= x^(alpha-1)*(1-x)^(beta-1)) # Determine whether u <= rho(x)
    {
      # accept x and record it
      m <- m + 1
      y[m] <- x
    }
  }
  if(picture){
    # Draw a histogram
    hist(y, prob = TRUE, breaks = 'scott',main = expression(f(x)==x^(alpha-1)*(1-x)^(beta-1)/B(alpha,beta)))
    # Sampling and plotting probability density curves
    sampling <- seq(0,1,.01)
    lines(sampling,sampling^(alpha-1)*(1-sampling)^(beta-1)/beta(alpha,beta))
  }
  return(y)
}

#' Generate random numbers of mixture Gamma_Exp distributes
#'
#' @param n random sample size
#' @param r parameters r
#' @param beta parameters beta
#' @param picture plot histogram
#'
#' @return random numbers of mixture Gamma_Exp distributes
#' @export
#'
#' @examples
#' \dontrun{
#' x <- rdistGammaExp(5,4,2,FALSE)
#' print(x)
#' }
#' 
rdistGammaExp = function(n = 1e3,r = 4,beta = 2,picture = FALSE){
  lambda_sampling <- rgamma(n,r,beta) # Generate gamma distributed random numbers
  y <- rexp(lambda_sampling) # the length of x and lambda is n = 1000
  if(picture){
    # Draw a histogram to show the random sample
    hist(y, prob = TRUE, breaks = 'scott',main = expression(y %~% Exp(lambda %~% Gamma(r,beta))))
    # Sampling and plotting probability density curves
    sampling <- seq(0,10,.01)
    lines(sampling,r*beta^r/(beta+sampling)^(r+1))
  }
  return(y)
}

