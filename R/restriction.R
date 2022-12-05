#' Grid coarsening function
#'
#' @param Matrix_A The coefficient matrix to be coarsened. It is required to be in matrix format.
#' @param n Degree of coarsening, that is, the number of nodes merged each time. The default is 2.
#' @param CoarsenType The coarsening method. The default is '\code{Projection}'.
#'
#' @return The projection matrix between the two layers.
#' @export
#'
#' @examples
#' \dontrun{
#' data(matrix_data_1e2)
#' attach(data)
#' A <- data$A
#' b <- data$b
#' Q <- restriction(A,n=2,CoarsenType='Projection')
#' }
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

