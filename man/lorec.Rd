\name{lorec}
\alias{lorec}
\title{LOw Rank and sparsE Covariance estimation}
\description{
  Estimate covariance matrices that contain low rank and sparse components
}
\usage{
lorec(Sig, L, S, lambda, delta,  thr=1.0e-4, maxit=1e4)
}
\arguments{
  \item{Sig}{Covariance matrix:p by p matrix (symmetric)}
  \item{L}{Initializing the low rank component. Default diag(Sig).}
  \item{S}{Initializing the sparse component. Default diag(Sig).}
  \item{lambda}{(Non-negative) regularization parameter (scalar) for the
    low rank component, \code{L}, via  the nuclear
    norm penalty.}
  \item{delta}{(Non-negative) regularization parameter (scalar) for
    the sparse component, \code{S}, via the L1
    penalty.  The diagonal of   \code{S} is not penalized.}
  \item{thr}{Threshold for convergence. 
    Iterations stop when the relative Frobenius norm change of
    consecutive updates is less than the threshold. Default is 1e-4.}
  \item{maxit}{Maximum number of iterations. Default 10,000.}
}

\details{
  Estimate a low rank plus sparse covariance matrix using a composite
  penalty, nuclear norm plus  L1 norm (lasso).  This
  covariance structure can  be verified in many classical
  models, such as factor and random effect models.  The algorithm is
  based on  Nesterov's method, suitable for large-scale problems with low memory
  costs.  It achieves the optimal global convergence rate of smooth
  problems under the black-box model.
}
\value{
  A list with components
  \item{L}{Estimated low rank component}
  \item{S}{Estimated sparse component}
}


\references{
Xi Luo (2011).  High Dimensional Low Rank and Sparse
Covariance Matrix Estimation via Convex Minimization. Technical Report,
Department of Biostatistics and Center of Statistical Sciences, Brown
University. arXiv: 1111.1133.
}
\examples{
set.seed(100)
Sig <- matrix(0.8, 4,4) + diag(0.2,4)
x<-matrix(rnorm(50*20),ncol=4)%*%chol(Sig)
s<- var(x)
a<-lorec(s, lambda=.01, delta=0.01)
}
\keyword{multivariate}
\keyword{models}


