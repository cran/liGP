\name{optIP.ALC}
\alias{optIP.ALC}


\title{
Sequential Selection of an Inducing Point Design by Optimizing Active Learning Cohn
}
\description{
  Optimizes the ALC surface to sequentially select inducing points for a given predictive location
  and local neighborhood. ALC can be based solely on the predictive location or an additional set of
  reference locations.
}
\usage{
optIP.ALC(Xc, Xref = NULL, M, Xn, Yn, theta = NULL, g = 1e-4,
          ip_bounds = NULL, num_thread = 1, num_multistart = 15,
          epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01,
          rep_list = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xc}{a vector containing the predictive location used as the center of the design/neighborhood}
  \item{Xref}{a \code{matrix} containing other reference locations used in the predictive variance
  sum that is minimized}
  \item{M}{
        the positive integer number of inducing points placed for each local neighborhood; \code{M} should be less than \code{N}
  }
  \item{Xn}{
        a \code{matrix} of the local neighborhood of \code{N} nearest neighbors to \code{Xc}
  }
  \item{Yn}{
        a vector of the corresponding responses to \code{Xn}; \code{length(Yn)=nrow(Xn)}
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
  \item{ip_bounds}{
    a 2 by d \code{matrix} containing the range of possible values for inducing points; first row contains minimum values for each dimension, second row contains maximum values; if \code{ip_bounds} is \code{NULL}, defaults to range of the local neighborhood \code{Xn}
  }
  \item{num_thread}{
      a scalar positive integer indicating the number of threads to use for parallel processing
      for the multistart search: \code{num_thread<=num_multistart}
    }
  \item{num_multistart}{
      a positive integer indicating the number of starting locations used to optimize the ALC
      and find the global minimum
    }
  \item{epsK}{
      a small positive number added to the diagonal of the correlation matrix, of inducing points, K,
      for numerically stability for inversion
    }
  \item{epsQ}{
        a small positive number added to the diagonal of the Q \code{matrix} (see Cole (2021))
        for numerically stability for inversion
  }
  \item{tol}{
    a positive number to serve as the tolerance level
    for covergence of ALC when optimizing the location
    of the next inducing point
  }
  \item{rep_list}{
  a list from \code{find_reps} in the \code{hetGP} package that details the replicates in Xn and their associated Yn
  }
  \item{verbose}{
    if \code{TRUE}, outputs the progress of the number of inducing points optimally placed
  }
}

\details{
  The function sequentially places \code{M} inducing points around the local neighborhood (\code{Xn}) of \code{Xc}. The first inducing point is placed at \code{Xc}. The remaining points and placed based on the minimum in the ALC surface using \code{rbind(Xc, Xref)} as a reference set for the predictive variance. Hyperparameters \code{theta,g} are fixed.
}
\value{
  The output is a \code{list} with the following components.

\item{Xm}{a \code{matrix} of the locally optimal inducing point locations}
\item{alc}{a vector of the ALC progress at each inducing point selection step}
\item{time}{a scalar giving the passage of wall-clock time elapsed
    for (substantive parts of) the calculation}
}
\references{
  D.A. Cole, R.B. Christianson, and R.B. Gramacy (2021).
  \emph{Locally Induced Gaussian Processes for Large-Scale Simulation Experiments}
  Statistics and Computing, 31(3), 1-21; preprint on arXiv:2008.12857;
  \url{https://arxiv.org/abs/2008.12857}

  S. Seo, M. Wallat, T. Graepel, and K. Obermayer (2000).
  \emph{Gaussian Process Regression: Active Data Selection and Test Point Rejection}
  In Mustererkennung 2000, 27-34. New York, NY: Springer-Verlag.

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

## Optimize inducing point locations
Xm.alc <- optIP.ALC(Xstar, Xref=NULL, M=m, Xn=Xn, Yn=Yn, theta=theta)

## Define reference locations for ALC sum
Xref <- as.matrix(expand.grid(Xstar[1]+c(-.05, 0, .05), Xstar[2]+c(-.05, 0, .05)))
Xm.alc_with_Xref <- optIP.ALC(Xstar, Xref=Xref, M=m, Xn=Xn, Yn=Yn, theta=theta)

## Plot inducing point design and neighborhood
ranges <- apply(Xn, 2, range)
plot(Xn, pch = 16, cex=.5, col='grey',
     xlim=ranges[,1]+c(-.02, .02), ylim=ranges[,2]+c(-.02, .15),
     xlab='x1', ylab = 'x2',
     main='ALC-optimal Inducing Point Design')
points(Xstar[1], Xstar[2], pch=16)
points(Xm.alc$Xm, pch=2, col=3, lwd=2)
points(Xm.alc_with_Xref$Xm, pch=3, col=4, lwd=2)
legend('topleft', col=c(1,'grey',3,4), pch=c(16,16,2,3), lty=NA, lwd=2, ncol=2,
       legend=c('Xstar','Neighborhood','Xm based on Xstar','Xm based on Xref'))

## View  ALC progress
plot(1:m, Xm.alc$alc, type='l', xlab='Inducing point number',
     ylab='ALC',main='ALC optimization progress')
}
}

