\name{StatComp22074-package}
\alias{StatComp22074-package}
\alias{StatComp22074}
\docType{package}
\title{
  Short Introduction of StatComp22074
}
\description{
  This package is a multigrid method based linear system solver and a presentation of previously completed jobs.
}
\details{
StatComp22074 is a simple R package, which mainly realizes the solution of the linear system Ax=b based on the Algebraic Multigrid (AMG) method. At the same time, it reproduces several iterative algorithms including Gauss-Seidel iteration method, Jacobi iteration method, Preconditioned Conjugate Gradient (PCG) method and so on. At the same time, it realizes the Coarsening algorithm based on projection, and realizes the comparison between algorithms.
}
\author{
22074

Maintainer: 22074 <1105084640@qq.com>
}
\references{
  This optional section can contain literature or other references for
  background information.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
  data(matrix_data_1e2)
  attach(data)
  before <- Multigrid_Setup(A,layer = 10,
  CoarsenType='Projection',
  InterpType='Projection',pri = TRUE)
  x <- Multigrid_Solve(Ah=before$Ah,b,Qh=before$Qh,
  x0=0,presmooth='Gauss-Seidel',v1=5,
  postsmooth='Gauss-Seidel',v2=5,solver='exact',
  tol=1e-6,max_iter=1e3,Precond=FALSE)
  }
}
