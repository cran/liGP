\name{optIP.wIMSE}
\alias{optIP.wIMSE}


\title{
Sequential Selection of an Inducing Point Design by Optimizing Weighted Integrates Mean-Sqaure Error
}
\description{
  Optimizes the weighted integrated mean-square error (wIMSE) surface to sequentially select
  inducing points for a given predictive location and local neighborhood.
}
\usage{
optIP.wIMSE(Xn, M, theta = NULL, g = 1e-4, w_mean, w_var = NULL,
            ip_bounds = NULL, integral_bounds = NULL, num_multistart = 15,
            fix_xm1 = TRUE, epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
            mult = NULL, verbose = TRUE)
}

\arguments{
  \item{Xn}{
        a \code{matrix} of the local neighborhood; \code{nrow(Xn)=N}
  }
  \item{M}{
        the positive integer number of inducing points placed for each local neighborhood; \code{M} should be less than \code{N}
  }

  \item{theta}{
        the lengthscale parameter (positive number) in a Gaussian
        correlation function; a (default) \code{NULL} value sets
        the lengthscale at the square of the 10th percentile of
      pairwise distances between neighborhood points
      (see \code{\link{darg}} in \code{\link{laGP}} package)
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
  \item{ip_bounds}{
    a 2 by d \code{matrix} containing the range of possible values for inducing points; first row contains minimum values for each dimension, second row contains maximum values; if \code{ip_bounds} is \code{NULL}, defaults to range of the local neighborhood \code{Xn}
  }
  \item{integral_bounds}{
    a 2 by d \code{matrix} containing the domain bounds for the data; first row contains minimum values for each dimension, second row contains maximum values; if \code{integral_bounds=NULL}, defaults to range of the local neighborhood \code{Xn}
  }
  \item{num_multistart}{
      a positive integer indicating the number of starting locations used to optimize wIMSE
      for each inducing point
    }
  \item{fix_xm1}{
    an indicator of whether or not the first inducing point should be fixed at \code{w_mean} (\code{TRUE}, default) or optimized (\code{FALSE})
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
  \item{verbose}{
    if \code{TRUE}, outputs the progress of the number of inducing points optimally placed
  }
}

\details{
  The function sequentially places \code{M} inducing points around the local neighborhood (\code{Xn}). Inducing points are placed based on the minimum in the wIMSE surface integrating over \code{integral_bounds}. Hyperparameters \code{theta,g} are fixed.
}
\value{
  The output is a \code{list} with the following components:

\item{Xm}{a \code{matrix} of the locally optimal inducing point locations}
\item{wimse}{a vector of the wIMSE progress at each inducing point selection step}
\item{time}{a scalar giving the passage of wall-clock time elapsed
    for (substantive parts of) the calculation}
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
## a "computer experiment"

## Simple 2-d Herbie's Tooth function used in Cole et al (2020);
## thanks to Lee, Gramacy, Taddy, and others who have used it before

## Build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)

\donttest{
library(laGP)

## Build a local neighborhood
Xstar <- matrix(c(.4, -1.1), nrow=1)
n <- 100; m <- 10
Xstar_to_X_dists <- distance(Xstar, X)
quant_dists <- quantile(Xstar_to_X_dists, n/nrow(X))
Xn <- X[Xstar_to_X_dists < quant_dists,]
Yn <- Y[Xstar_to_X_dists < quant_dists]

theta <- darg(NULL, Xn)$start
integral_bounds <- rbind(c(-2,-2), c(2,2))

## Optimize inducing point locations
Xm.wimse1 <- optIP.wIMSE(Xn, M=m, Xn=, theta=theta, w_mean=Xstar,
                        integral_bounds=integral_bounds)

## Use a different variance for weight
Xm.wimse2 <- optIP.wIMSE(Xn, M=m, Xn=, theta=theta, w_mean=Xstar,
                        w_var = c(theta/2, 3*theta),
                        integral_bounds=integral_bounds)

## Plot inducing point design and neighborhood
ranges <- apply(Xn, 2, range)
plot(Xn, pch = 16, cex=.5, col='grey',
     xlim=ranges[,1]+c(-.02, .02), ylim=ranges[,2]+c(-.02, .15),
     xlab='x1', ylab = 'x2',
     main='ALC-optimal Inducing Point Design')
points(Xstar[1], Xstar[2], pch=16)
points(Xm.wimse1$Xm, pch=2, col=3, lwd=2)
points(Xm.wimse2$Xm, pch=3, col=4, lwd=2)
legend('topleft', legend=c('Xstar','Neighborhood','Xm with w_var=theta',
                           'Xm with nonisotropic weight'),
       col=c(1,'grey',3,4), pch=c(16,16,2,3), lty=NA, lwd=2, ncol=2)

## View wIMSE progress
plot(1:m, log(Xm.wimse1$wimse), type='l', xlab='inducing point number',
     ylab='log wIMSE',main='wIMSE optimization progress')
}
}
