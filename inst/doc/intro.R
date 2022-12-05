## -----------------------------------------------------------------------------
Multigrid_Setup <- function(A,layer=0,CoarsenType='Projection',InterpType='Projection',pri=FALSE){
  Ah <- list(A)
  Qh <- c()
  dimension_A <- dim(A)[1]
  nnz <- sum(A!=0) # count non-zero number
  dense <- c(nnz/(dimension_A^2)) # calculate dense of matrix
  q <- max(floor((1 - dense)*5),1) # Adjust coarsening degree according to dense
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
  # Print specific details of each layer change
  if(pri){
    pridata <- data.frame(dim_of_A=dimension_A,nnz=nnz,dense=dense,layer=1:i-1)
    print(pridata)
  }
  return(list(Ah=Ah,Qh=Qh,nnz=nnz,dense=dense,layer=i))
}

## -----------------------------------------------------------------------------
# One-cycle
Multigrid_iter <- function(Ah,b,Qh,x0=0,presmooth='Gauss-Seidel',v1=2,
                           postsmooth='Gauss-Seidel',v2=2,solver='exact'){
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

## -----------------------------------------------------------------------------
Multigrid_Solve <- function(Ah,b,Qh,x0=0,presmooth='Gauss-Seidel',v1=2,
                            postsmooth='Gauss-Seidel',v2=2,solver='exact',
                            tol=1e-5,max_iter=1e3,Precond=FALSE){
  error <- 1e3
  iter_time <- 0
  tag <- 1
  # precondition method
  if(Precond){
    # use PCG method
    result <- PCG_MG(Ah[[1]],b,x0,tol,max_iter,Ah,Qh,presmooth,v1,postsmooth,v2,solver)
    tag <- result$tag
    iter_time <- result$iter_time
    x0 <- result$x
    error <- result$error
  }
  # not precondition method
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
  # show the result
  return(list(x=x0,iter_time=iter_time,error=error,tag=tag))
}

## -----------------------------------------------------------------------------
Smoother <- function(Matrix_A, Matrix_b,x0=0,N=100,tol=-1,max_iter=1e3,SmoothMethod='Gauss-Seidel'){
  # Check whether the input meets the requirements
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
  tag <- 1 # tag whether converge
  m <- dim(Matrix_b)[1];n <- dim(x0)
  x_old <- x0 - 10
  x_new <- x0
  # Gauss-Seidel
  if(SmoothMethod == 'Gauss-Seidel'){
    if(tol<0){
      for(k in 1:N){
        x_old <- x_new
        x_new[1] = (Matrix_b[1,1] - sum(Matrix_A[1,2:n]*x_old[2:n]))/Matrix_A[1,1]
        for(i in 2:(n-1)){
          x_new[i] = (Matrix_b[i,1] - sum(Matrix_A[i,1:(i-1)]*x_new[1:(i-1)]) 
                      - sum(Matrix_A[i,(i+1):n]*x_old[(i+1):n]))/Matrix_A[i,i]
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
          x_new[i] = (Matrix_b[i,1] - sum(Matrix_A[i,1:(i-1)]*x_new[1:(i-1)]) 
                      - sum(Matrix_A[i,(i+1):n]*x_old[(i+1):n]))/Matrix_A[i,i]
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
  # Jacobi 
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
  # exact solve
  else if(SmoothMethod == 'exact'){
    k = -1
    tag = -1
    x_new = array(solve(Matrix_A)%*%Matrix_b)
    x_old = x_new
  }
  else{
    return(warning('No this smoothing method.'))
  }
  return(list(x=x_new,iter_time=k,tag=tag,tol = norm(matrix(x_old - x_new))/norm(matrix(x_old)),
              error=norm(Matrix_A%*%x_new-Matrix_b)/norm(Matrix_b)))
}

## -----------------------------------------------------------------------------
restriction <- function(Matrix_A,n=2,CoarsenType='Projection'){
  dimC <- dim(Matrix_A)
  if(dimC[1]%%n==0){
    Q <- matrix(rep(0,dimC[1]/n*dimC[1]),dimC[1]/n,dimC[1])
    for(i in 1:(dimC[1]/n)){
      Q[i,(n*(i-1)+1):(n*i)] <- 1/sqrt(n)
    }
  }
  else{
    q <- dimC[1] %/% n + 1
    r <- dimC[1] %% n
    Q <- matrix(rep(0,q*dimC[1]),q,dimC[1])
    for(i in 1:(q-1)){
      Q[i,(n*(i-1)+1):(n*i)] <- 1/sqrt(n)
    }
    Q[q,(n*i+1):dimC[1]] <- 1/sqrt(r)
  }
  return(Q)
}

## -----------------------------------------------------------------------------
PCG_MG <- function(A,b,x,tol=1e-5,max_iter=1e3,Ah,Qh,
                   presmooth='Gauss-Seidel',v1=2,
                   postsmooth='Gauss-Seidel',v2=2,solver='exact'){
  if(length(x) <= 1){
    x = array(rep(0,dim(A)[2]))
  }
  tag <- 1
  r_new = b - A %*% x
  # preconditioner
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
    # end repeat
    if(error<tol){
      break
    }
    # preconditioner
    z = Multigrid_iter(Ah,r_new,Qh,x0=x,presmooth,v1,postsmooth,v2,solver)
    z_new = z$x
    beta = t(r_new) %*% z_new / t(r_old) %*% z_old
    p = z_old + beta[1,1] * p
    k = k + 1
    # exceed max iteration time
    if(sum((alpha[1,1] * p)^2)/sum(x^2)<tol | k >= max_iter){
      tag <- 0
      break
    }
  }
  return(list(x=x,error=error,iter_time=k,tag=tag))
}

## ----eval=TRUE----------------------------------------------------------------
library(Rcpp)
library(StatComp22074)

data("matrix_data_1e3")
attach(data)
set.seed(1)
n <- dim(A)[1]
x0 <- array(rep(0,n))

before <- Multigrid_Setup(A,layer = 10,CoarsenType='Projection',
                          InterpType='Projection',pri = TRUE)
start_time1 <- Sys.time()
x_notPre <- Multigrid_Solve(Ah=before$Ah,b,Qh=before$Qh,x0=x0,
                            presmooth='Gauss-Seidel',v1=5,postsmooth='Gauss-Seidel',v2=5,
                            solver='exact',tol=1e-6,max_iter=1e3,Precond=FALSE)
end_time1 <- Sys.time()
start_time2 <- Sys.time()
x_Pre <- Multigrid_Solve(Ah=before$Ah,b,Qh=before$Qh,x0=x0,
                         presmooth='Gauss-Seidel',v1=5,postsmooth='Gauss-Seidel',v2=5,
                         solver='exact',tol=1e-6,max_iter=1e3,Precond=TRUE)
end_time2 <- Sys.time()

N <- 5e2
start_time3 <- Sys.time()
x_R <- Smoother(A,b,x0,tol=1e-10,max_iter=N,SmoothMethod = 'Gauss-Seidel')
end_time3 <- Sys.time()

start_time4 <- Sys.time()
x_Cpp <- Gauss_Seidel(A,b,x0,N)
end_time4 <- Sys.time()

data_Count <- data.frame(method=c('AMG','Precondition AMG','Gauss Seidel_R','Gauss Seidel_Rcpp'),
                         time=c(end_time1-start_time1,end_time2-start_time2,end_time3-start_time3,
                                end_time4-start_time4),
                         iter_time = c(x_notPre$iter_time,x_Pre$iter_time,x_R$iter_time,N),
                         error=c(x_notPre$error,x_Pre$error,x_R$error,norm(b-A%*%x_Cpp)/norm(b)))

knitr::kable(data_Count)

