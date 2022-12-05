#' Smooth function
#' @description Multiple grid each step of the iterator selection, you can set the iteration step number and precision information.
#'
#' @param Matrix_A Matrix of coefficients. It is required to be in matrix format.
#' @param Matrix_b The right-hand side of a linear system. It is required to be in matrix format.
#' @param x0 Iterate the initial value. The default is an all-0 vector with the same dimension as the matrix A. It is required to be in array format.
#' @param N The number of iteration steps of the smoother. If the precision \code{tol} is set, the iteration step number parameter is meaningless. The default value is 100.
#' @param tol Iterator convergence precision, if not set, iterate by the number of iteration steps N. The default value is -1, that is, this parameter is not set.
#' @param max_iter Maximum number of iteration steps. The default value is \code{1e3}.
#' @param SmoothMethod Method of smoothing, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#'
#' @return Returns a list. The list includes the result vector \code{x}, number of iterations \code{iter_time}, the corresponding \code{error}, precision of convergence \code{tol} and the flag \code{tag} of convergence or not.
#' @export
#'
#' @examples
#' \dontrun{
#' data(matrix_data_1e2)
#' attach(data)
#' x <- Smoother(A,b,x0=0,N=10,SmoothMethod='Gauss-Seidel')
#' }
#' 
Smoother <- function(Matrix_A, Matrix_b,x0=0,N=100,tol=-1,max_iter=1e3,SmoothMethod='Gauss-Seidel'){
  if(!is.array(Matrix_A)){
    return(warning('A is not matrix.'))
  }
  else if(!is.array(Matrix_b)){
    return(warning('b is not array.'))
  }
  else if(N %% 1 != 0){
    return(warning('The number of iterations should be an integer.'))
  }
  else if(dim(Matrix_A)[1]!=dim(Matrix_b)[1]){
    return(warning('Dimension of A and b is not match.'))
  }
  else if(dim(Matrix_A)[1]!=dim(Matrix_A)[2]){
    return(warning('A is not a square matrix.'))
  }
  else if(length(x0) <= 1){
    x0 = array(rep(0,dim(Matrix_A)[2]))
  }
  else if(!is.array(x0)){
    return(warning('x0 is not array.'))
  }
  else if(dim(Matrix_A)[2]!=dim(x0)){
    return(warning('Dimension of A and x0 is not match.'))
  }
  tag <- 1
  
  m <- dim(Matrix_b)[1];n <- dim(x0)
  x_old <- x0 - 10
  x_new <- x0
  if(SmoothMethod == 'Gauss-Seidel'){
    if(tol<0){
      for(k in 1:N){
        x_old <- x_new
        x_new[1] = (Matrix_b[1,1] - sum(Matrix_A[1,2:n]*x_old[2:n]))/Matrix_A[1,1]
        for(i in 2:(n-1)){
          x_new[i] = (Matrix_b[i,1] - sum(Matrix_A[i,1:(i-1)]*x_new[1:(i-1)]) - sum(Matrix_A[i,(i+1):n]*x_old[(i+1):n]))/Matrix_A[i,i]
        }
        x_new[n] = (Matrix_b[n,1] - sum(Matrix_A[n,1:(n-1)]*x_new[1:(n-1)]))/Matrix_A[n,n]
      }
    }
    else{
      k <- 0
      while(norm(Matrix_A%*%x_new-Matrix_b)/norm(Matrix_b)>=tol){
        x_old <- x_new
        x_new[1] = (Matrix_b[1,1] - sum(Matrix_A[1,2:n]*x_old[2:n]))/Matrix_A[1,1]
        for(i in 2:(n-1)){
          x_new[i] = (Matrix_b[i,1] - sum(Matrix_A[i,1:(i-1)]*x_new[1:(i-1)]) - sum(Matrix_A[i,(i+1):n]*x_old[(i+1):n]))/Matrix_A[i,i]
        }
        x_new[n] = (Matrix_b[n,1] - sum(Matrix_A[n,1:(n-1)]*x_new[1:(n-1)]))/Matrix_A[n,n]
        k <- k+1
        if(norm(matrix(x_old - x_new))/norm(matrix(x_old))< tol){
          break
        }
        if(k>=max_iter){
          tag <- 0
          break
        }
      }
    }
  }
  else if(SmoothMethod == 'Jacobi'){
    if(tol<0){
      for(k in 1:N){
        x_old <- x_new
        x_new = array((Matrix_b[,1] - Matrix_A %*% x_old + diag(diag(Matrix_A))%*%x_old)/diag(Matrix_A))
      }
    }
    else{
      k <- 0
      while(norm(Matrix_A%*%x_new-Matrix_b)/norm(Matrix_b)>=tol){
        x_old <- x_new
        x_new = array((Matrix_b[,1] - Matrix_A %*% x_old + diag(diag(Matrix_A))%*%x_old)/diag(Matrix_A))
        k <- k + 1
        if(norm(matrix(x_old - x_new))/norm(matrix(x_old))< tol){
          break
        }
        if(k>=max_iter){
          tag <- 0
          break
        }    
      }
    }
  }
  else if(SmoothMethod == 'exact'){
    k = -1
    tag = -1
    x_new = array(solve(Matrix_A)%*%Matrix_b)
    x_old = x_new
  }
  else{
    return(warning('No this smoothing method.'))
  }
  return(list(x=x_new,iter_time=k,tag=tag,tol = norm(matrix(x_old - x_new))/norm(matrix(x_old)),error=norm(Matrix_A%*%x_new-Matrix_b)/norm(Matrix_b)))
}
