#' Multiple grid building functions
#' 
#' @description Create different levels of multiple grids.
#'
#' @param A Matrix of coefficients. It is required to be in matrix format.
#' @param layer The maximum number of layers. If the input is less than 1, the changes between different layers are designed according to the sparsity of matrix A by default. The default value is \code{0}.
#' @param CoarsenType The coarsening method. The default is '\code{Projection}'.
#' @param InterpType Interpolation method. The default is '\code{Projection}'
#' @param pri Whether to print the details of each layer. The default value is \code{FALSE}.
#'
#' @return Returns a list. The list includes the coefficient matrix \code{Ah}, the projection matrix \code{Qh}, the number of non-zero elements \code{nnz}, the \code{dense} for each layer, and the number of \code{layer}s.
#' @export
#'
#' @examples
#' \dontrun{
#' data(matrix_data_1e2)
#' attach(data)
#' before <- Multigrid_Setup(A,layer = 10,CoarsenType='Projection',
#' InterpType='Projection',pri = TRUE)
#' }
#' 
Multigrid_Setup <- function(A,layer=0,CoarsenType='Projection',InterpType='Projection',pri=FALSE){
  Ah <- list(A)
  Qh <- c()
  dimension_A <- dim(A)[1]
  nnz <- sum(A!=0)
  dense <- c(nnz/(dimension_A^2))
  q <- max(floor((1 - dense)*5),1)
  i <- 1
  # Set up A
  while(q > 1){
    Q <- restriction(Ah[[i]],n=q,CoarsenType=CoarsenType)
    Qh <- c(Qh,list(Q))
    Ah <- c(Ah,list(Q%*%Ah[[i]]%*%t(Q)))
    dimension_A <- c(dimension_A,dim(Ah[[i+1]])[1])
    nnz <- c(nnz,sum(Ah[[i+1]]!=0))
    dense <- c(dense,nnz[i+1]/(dimension_A[i+1]^2))
    i <- i + 1
    if(i == layer){
      break
    }
    q <- max(floor((1 - dense[i])*10),1)
  }
  if(pri){
    pridata <- data.frame(dim_of_A=dimension_A,nnz=nnz,dense=dense,layer=1:i-1)
    print(pridata)
  }
  return(list(Ah=Ah,Qh=Qh,nnz=nnz,dense=dense,layer=i))
}

#' A single cycle of the multigrid method
#' @description A single step loop of the multiple grid method.
#'
#' @param Ah Coefficient matrix \code{A} Matrix of different scales on each layer of grid output in \code{Multigrid_Setup()}.
#' @param b The right-hand side of a linear system. It is required to be in matrix format.
#' @param Qh A projection array between different grids output in \code{Multigrid_Setup()}.
#' @param x0 Iterate the initial value. It is required to be in array format. The default is an all-0 vector with the same dimension as the matrix A.
#' @param presmooth A pre-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v1 The number of iteration steps of the pre-smooth iterator. The default value is 2.
#' @param postsmooth A post-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v2 The number of iteration steps of the post-smooth iterator. The default value is 2.
#' @param solver The solver for the coarsest layer. The default is \code{'exact'}.
#'
#' @return Returns a list. The list includes the result vector \code{x} after a cycle and the corresponding \code{error}.
#' @export
#'
#' 
Multigrid_iter <- function(Ah,b,Qh,x0=0,presmooth='Gauss-Seidel',v1=2,
                           postsmooth='Gauss-Seidel',v2=2,
                           solver='exact'){
  fh <- list(b)
  vh <- c()
  layer <- length(Qh)
  # Coarsening
  for(i in 1:layer){
    x <- Smoother(Ah[[i]],fh[[i]],x0,N=v1,SmoothMethod=presmooth)
    vh <- c(vh,list(x$x))
    f <- Qh[[i]]%*%(fh[[i]] - Ah[[i]] %*% vh[[i]])
    fh <- c(fh,list(f))
    x0 <- if(length(x0)>1) array(Qh[[i]]%*%x0) else x0
  }
  # solver
  layer_vh <- Smoother(Ah[[layer+1]],fh[[layer+1]],SmoothMethod=solver)
  vh <- c(vh,list(layer_vh$x))
  vsh <- vh
  # interpolation
  for(i in layer:1){
    f <- vh[[i]] + array(t(Qh[[i]])%*%vsh[[i+1]])
    x <- Smoother(Ah[[i]],fh[[i]],f,N=v2,SmoothMethod=postsmooth)
    vsh[[i]] <- x$x
  }
  return(list(x=vsh[[1]],error=norm(Ah[[1]]%*%vsh[[1]]-b)/norm(b)))
}

#' Multiple grid method solver
#' @description Solve according to the constructed grid.
#' 
#' @param Ah Coefficient matrix \code{A} Matrix of different scales on each layer of grid output in \code{Multigrid_Setup()}.
#' @param b The right-hand side of a linear system. It is required to be in matrix format.
#' @param Qh A projection array between different grids output in \code{Multigrid_Setup()}.
#' @param x0 Iterate the initial value. It is required to be in array format. The default is an all-0 vector with the same dimension as the matrix A.
#' @param presmooth A pre-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v1 The number of iteration steps of the pre-smooth iterator. The default value is 2.
#' @param postsmooth A post-smooth iterator, optionally \code{c('Gauss-Seidel', 'Jacobi', 'exact')}. The default value is \code{'Gauss-Seidel'}.
#' @param v2 The number of iteration steps of the post-smooth iterator. The default value is 2.
#' @param solver The solver for the coarsest layer. The default is \code{'exact'}.
#' @param tol Precision of iteration results. The default value is \code{1e-5}.
#' @param max_iter The maximum number of iteration steps. If this deployment is exceeded and the accuracy is not reached, the iteration stops. The default value is \code{1e3}.
#' @param Precond Whether to use multiple grids as the preconditioner, if \code{TRUE} use the preconditioned conjugate gradient method. The default value is \code{FALSE}.
#'
#' @return Returns a list. The list includes the result vector \code{x}, number of iterations \code{iter_time}, the corresponding \code{error} and the flag \code{tag} of convergence or not.
#' @export
#'
#' @examples
#' \dontrun{
#' data(matrix_data_1e2)
#' attach(data)
#' before <- Multigrid_Setup(A,layer = 10,CoarsenType='Projection',InterpType='Projection',pri = TRUE)
#' x <- Multigrid_Solve(Ah=before$Ah,b,Qh=before$Qh,x0=0,presmooth='Gauss-Seidel',v1=5,
#' postsmooth='Gauss-Seidel',v2=5,solver='exact',tol=1e-6,max_iter=1e3,Precond=FALSE)
#' }
#' 
Multigrid_Solve <- function(Ah,b,Qh,x0=0,presmooth='Gauss-Seidel',v1=2,
                            postsmooth='Gauss-Seidel',v2=2,
                            solver='exact',tol=1e-5,max_iter=1e3,Precond=FALSE){
  error <- 1e3
  iter_time <- 0
  tag <- 1
  if(Precond){
    result <- PCG_MG(Ah[[1]],b,x0,tol,max_iter,Ah,Qh,presmooth,v1,postsmooth,v2,solver)
    tag <- result$tag
    iter_time <- result$iter_time
    x0 <- result$x
    error <- result$error
  }
  else{
    while(error>tol){
      result <- Multigrid_iter(Ah,b,Qh,x0,presmooth,v1,postsmooth,v2,solver)
      x0 <- result$x
      error <- result$error
      iter_time <- iter_time + 1
      if((norm(matrix(result$x - x0))/norm(matrix(result$x))<tol) | (iter_time >= max_iter)){
        tag <- 0
        break
      }
    }
  }
  return(list(x=x0,iter_time=iter_time,error=error,tag=tag))
}

