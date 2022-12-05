#' Precondition conjugate gradient method with AMG as the preconditioner
#'
#' @param A Matrix of coefficients. It is required to be in matrix format.
#' @param b The right-hand side of a linear system. It is required to be in matrix format.
#' @param x Iterate the initial value. The default is an all-0 vector with the same dimension as the matrix A. It is required to be in array format.
#' @param tol Iterator convergence precision. The default value is \code{1e-5}.
#' @param max_iter The maximum number of iteration steps. If this deployment is exceeded and the accuracy is not reached, the iteration stops. The default value is \code{1e3}.
#' @param Ah Coefficient matrix \code{A} Matrix of different scales on each layer of grid output in \code{Multigrid_Setup()}.
#' @param Qh A projection array between different grids output in \code{Multigrid_Setup()}.
#' @param presmooth A pre-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v1 The number of iteration steps of the pre-smooth iterator. The default value is 2.
#' @param postsmooth A post-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v2 The number of iteration steps of the post-smooth iterator. The default value is 2.
#' @param solver The solver for the coarsest layer. The default is \code{'exact'}.
#' 
#'
#' @return Returns a list. The list includes the result vector \code{x}, number of iterations \code{iter_time}, the corresponding \code{error} and the flag \code{tag} of convergence or not.
#' @export
#'
#' 
PCG_MG <- function(A,b,x,tol=1e-5,max_iter=1e3,Ah,Qh,
                   presmooth='Gauss-Seidel',v1=2,
                   postsmooth='Gauss-Seidel',v2=2,
                   solver='exact'){
  if(length(x) <= 1){
    x = array(rep(0,dim(A)[2]))
  }
  tag <- 1
  r_new = b - A %*% x
  z = Multigrid_iter(Ah,r_new,Qh,x0=x,presmooth,v1,postsmooth,v2,solver)
  z_new = z$x
  p = z_new
  k = 0
  while(1){
    r_old = r_new
    z_old = z_new
    alpha = t(r_old) %*% z_old / (p %*% A %*% p)
    x = x + alpha[1,1] * p
    r_new = r_old - alpha[1,1] * A %*% p
    error = norm(r_new)/norm(b)
    if(error<tol){
      break
    }
    z = Multigrid_iter(Ah,r_new,Qh,x0=x,presmooth,v1,postsmooth,v2,solver)
    z_new = z$x
    beta = t(r_new) %*% z_new / t(r_old) %*% z_old
    p = z_old + beta[1,1] * p
    k = k + 1
    if(sum((alpha[1,1] * p)^2)/sum(x^2)<tol | k >= max_iter){
      tag <- 0
      break
    }
  }
  return(list(x=x,error=error,iter_time=k,tag=tag))
}

#' Conjugate gradient method
#'
#' @param A Matrix of coefficients. It is required to be in matrix format.
#' @param b The right-hand side of a linear system. It is required to be in matrix format.
#' @param x Iterate the initial value. The default is an all-0 vector with the same dimension as the matrix A. It is required to be in array format.
#' @param tol Iterator convergence precision. The default value is \code{1e-5}.
#' @param max_iter The maximum number of iteration steps. If this deployment is exceeded and the accuracy is not reached, the iteration stops. The default value is \code{1e3}.
#' 
#' @return Returns a list. The list includes the result vector \code{x}, number of iterations \code{iter_time}, the corresponding \code{error} and the flag \code{tag} of convergence or not.
#' @export
#'
CG <- function(A,b,x=0,tol=1e-5,max_iter=1e3){
  if(length(x) <= 1){
    x = rep(0,dim(A)[2])
  }
  tag <- 1
  r_new = b - A %*% x
  p = r_new
  k = 0
  while(1){
    r_old = r_new
    alpha = (t(r_old) %*% r_old) / (t(p) %*% A %*% p)
    x = x + alpha[1,1] * p
    r_new = r_old - alpha[1,1] * A %*% p
    error = ((t(r_new)%*%r_new)/(t(b)%*%b))[1,1]
    if(error<tol){
      break
    }
    beta = (t(r_new) %*% r_new) / (t(r_old) %*% r_old)
    p = r_old + beta[1,1] * p
    k = k + 1
    if(sum((alpha[1,1] * p)^2)/sum(x^2)<tol | k >= max_iter){
      tag <- 0
      break
    }
  }
  return(list(x=x,error=error,iter_time=k,tag=tag))
}
