\name{calc_IMSE}
\alias{calc_IMSE}


\title{
  Integrated Mean-Square Error Given a New Inducing Point
}
\description{
  Calculates the Integrated Mean-Square Error (IMSE) given a set of data points,
  inducing point design, and new proposed inducing point location.
}
\usage{
calc_IMSE(xm1, Xm = NULL, X, theta = NULL, g = 1e-4,
          integral_bounds = NULL, epsK = sqrt(.Machine$double.eps),
          epsQ = 1e-5, mult = NULL)
}
\arguments{
  \item{xm1}{a vector containg the location of a proposed inducing point}
  \item{Xm}{
        optional design \code{matrix} of existing inducing points; \code{ncol(Xm) = length(xm1)}}
  \item{X}{
        the design \code{matrix} of input locations;
        \code{ncol(X) = length(xm1)}
  }
  \item{theta}{
        the lengthscale parameter (positive number) in a Gaussian
        correlation function;  a (default) \code{NULL} value sets
        the lengthscale at the square of the 10th percentile of
      pairwise distances between input locations \code{X}
      (similar to \code{\link{darg}} in \code{\link{laGP}} package)
  }
  \item{g}{
    the nugget parameter (positive number) in the covariance
  }
  \item{integral_bounds}{
    a 2 by d \code{matrix} containing the domain bounds for the data; first row contains minimum values for each dimension, second row contains maximum values; if \code{integral_bounds} is \code{NULL}, defaults to range of the input locations \code{X}
  }
  \item{epsK}{
      a small positive number added to the diagonal of the correlation matrix, of inducing points, K,
      for numerically stability for inversion
    }
  \item{epsQ}{
        a small positive number added to the diagonal of the Q \code{matrix} (see Cole (2021))
        for numerically stability for inversion
  }
  \item{mult}{
      an optional vector of length \code{nrow(X)} that contains the number of replicates for each design location in \code{X}
  }
}

\details{
  The function calculates the integrated mean-square error over the provided domain
  (\code{integral_bounds}). The IMSE is calculated in closed-form.
}
\value{
  the integrated mean-sqaure error
}
\references{
  D.A. Cole, R.B. Christianson, and R.B. Gramacy (2021).
  \emph{Locally Induced Gaussian Processes for Large-Scale Simulation Experiments}
  Statistics and Computing, 31(3), 1-21; preprint on arXiv:2008.12857;
  \url{https://arxiv.org/abs/2008.12857}

}
\author{
  D. Austin Cole \email{austin.cole8@vt.edu}
}

\examples{

## Build a set of input locations and existing inducing points
X = matrix(runif(100), ncol=2)
Xm = matrix(runif(10), ncol=2)

integral_bounds <- rbind(c(0,0), c(1,1))
xm1_new <- c(.4, .2)

## Calculate the integrated mean-square error
calc_IMSE(xm1=xm1_new, Xm=Xm, X=X,
          integral_bounds=integral_bounds)

## without an existing inducing point design
calc_IMSE(xm1=xm1_new, Xm=NULL, X=X,
          integral_bounds=integral_bounds)
}
