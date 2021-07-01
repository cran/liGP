\name{calc_wIMSE}
\alias{calc_wIMSE}


\title{
Weighted Integrated Mean-Square Error Given a New Inducing Point
}
\description{
  Calculates the Weighted Integrated Mean-Square Error (wIMSE) given a prediction location,
  local neighborhood, design of inducing points, and new proposed inducing point location.
}
\usage{
calc_wIMSE(xm1, Xm = NULL, Xn, theta = NULL, g = 1e-4,
           w_mean, w_var = NULL, integral_bounds = NULL,
           epsK = sqrt(.Machine$double.eps),
           epsQ = 1e-5, mult = NULL)
}
\arguments{
  \item{xm1}{a vector containg the location of a proposed inducing point}
  \item{Xm}{
        a design \code{matrix} of existing inducing points; \code{ncol(Xm) = length(xm1)}}
  \item{Xn}{
        a \code{matrix} of the local neighborhood;
        \code{ncol(Xn) = length(xm1)}}
  \item{theta}{
        the lengthscale parameter (positive number) in a Gaussian
        correlation function;  a (default) \code{NULL} value sets
        the lengthscale at the square of the 10th percentile of
      pairwise distances between neighborhood points
      (similar to \code{\link{darg}}
      in \code{\link{laGP}} package)
  }
  \item{g}{
    the nugget parameter (positive number) in the covariance
  }
  \item{w_mean}{
    a vector of the mean (center) of the Gaussian weight; \code{length(w_mean)} should equal \code{ncol(Xn)}
  }
  \item{w_var}{
    a positive number or vector of positive numbers (length equal to \code{ncol(Xn)}) denoting the variance(s) in the Gaussian weight. If \code{NULL} (default), \code{theta} is used.
  }
  \item{integral_bounds}{
    a 2 by d \code{matrix} containing the domain bounds for the data; first row contains minimum values for each dimension, second row contains maximum values; if \code{integral_bounds}=\code{NULL}, defaults to range of the local neighborhood \code{Xn}
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
      a vector of length \code{nrow(X)} that contains the number of replicates for each design location in \code{X}
  }
}

\details{
  The function calculates the integrated mean-square error with a Gaussian weight with mean \code{w_mean}
  (i.e. predictive location) and variance \code{w_var}. By using a Gaussian weight along with a Gaussian kernel for the GP, the wIMSE is calculated in closed-form.
}
\value{
  the weighted integrated mean-sqaure error
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
## Build a "local neighborhood" and existing inducing point design
X_center <- c(.5, .5)
Xn <- matrix(runif(100), ncol=2)
Xm <- matrix(runif(10), ncol=2)

integral_bounds <- rbind(c(0,0), c(1,1))
xm1_new <- c(.4, .2)

## Calculate the weighted integrated mean-square error
calc_wIMSE(xm1=xm1_new,  Xm=Xm, Xn=Xn, w_mean=X_center,
           integral_bounds=integral_bounds)

## Define weight's variance
calc_wIMSE(xm1=xm1_new, Xm=Xm, Xn=Xn, w_mean=X_center,
           w_var=c(.1, .2), integral_bounds=integral_bounds)

## Without an exisiting inducing point design
calc_wIMSE(xm1=xm1_new, Xm=NULL, Xn=Xn, w_mean=X_center,
           integral_bounds=integral_bounds)
}
