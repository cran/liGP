\name{giGP}
\alias{giGP}


\title{
Global Inducing Point Approximate GP Regression For Many Predictive Locations
}
\description{
  Facilitates Gaussian process inference and prediction at a set of predictive locations
  through the implementation of an inducing point design. Optimizes hyperparameters and
  returns the moments of the predictive equations.
}
\usage{
giGP(XX, X = NULL, Y = NULL, Xm, g = 1e-6, theta = NULL, nu = NULL,
     epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01, reps = FALSE)
}
\arguments{
  \item{XX}{
      a \code{matrix} of out-of-sample
      predictive locations with \code{ncol(XX)=ncol(X)}
  }
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of all input locations. If using a list for \code{reps}, this entry is not used.
}
\item{Y}{
  a vector of all responses/dependent values with \code{length(Y)=nrow(X)}. If using a list for \code{reps}, this entry is not used.
}
  \item{Xm}{
    a \code{matrix} containing the inducing points design with \code{ncol(Xm)=ncol(X)}.
  }
 \item{g}{
        an initial setting or fixed value for the nugget parameter. In order to optimize the nugget, a list can be provided that includes:
\itemize{
  \item \code{start} -- starting value to initialize the nugget
  \item \code{min} -- minimum value in the allowable range for the nugget
  \item \code{max} -- maximum value in the allowable range for the nugget
  \item \code{ab} -- shape and rate parameters specifying a Gamma prior for the nugget
}
If \code{ab} is not provided, a prior is not placed with the likelihood for optimization. If \code{min} and \code{max} aren't provided, the nugget is not optimized. If a single positive scalar is provided, the nugget is fixed for all predictions. If \code{NULL}, an initial setting is based on \code{garg} in the \code{laGP} package.
  }
  \item{theta}{
        an initial setting or fixed value for the lengthscale parameter. A (default) \code{NULL} value generates an initial setting based on \code{darg} in the \code{laGP} package. Similarly, a list can be provided that includes:
\itemize{
  \item \code{start} -- starting value to initialize the lengthscale
  \item \code{min} -- minimum value in the allowable range for the lengthscale
  \item \code{max} -- maximum value in the allowable range for the lengthscale
  \item \code{ab} -- shape and rate parameters specifying a Gamma prior for the lengthscale
}
If \code{ab} is not provided, a prior is not placed with the likelihood for optimization. If \code{min} and \code{max} aren't provided, the lengthscale is not optimized. If a single positive scalar is provided, the lengthscale is fixed for all predictions.
  }
  \item{nu}{
    a positive number used to set the scale parameter;
    default (\code{NULL}) calculates the maximum likelihood estimator
  }
  \item{epsK}{
      a small positive number added to the diagonal of the correlation \code{matrix} of inducing points
      for numerically stability for inversion. The value is automatically increased if needed.
    }
  \item{epsQ}{
        a small positive number added to the diagonal
        of the Q \code{matrix} (see Cole (2021))
        for numerically stability for inversion. The value is automatically increased if needed.
  }
  \item{tol}{
    a positive number to serve as the tolerance level
    for covergence of the log-likelihood when optimizing
    the hyperparameter(s) theta, g
  }
    \item{reps}{
   a notification of replicate design locations in the data set.
   If \code{TRUE}, the unique design locations are used for the calculations along with the average response for each unique design location. Alternatively, \code{reps} can be a list from \code{find_reps} in the \code{hetGP} package. In this case, \code{X} and \code{Y} are not used.
}
}
\details{
  The function uses the likelihood and predictive equations derived in
  Snelson and Ghahramani (2006) to fit a induced Gaussian Process for predictions.
  All the data \code{\{X,Y\}} and inducing points \code{Xm} are used for each prediction.

}
\value{
  The output is a \code{list} with the following components.

  \item{mean}{a vector of predictive means of length \code{nrow(XX)}}
  \item{var}{a vector of predictive variances of length
    \code{nrow(XX)}}
  \item{nu}{a vector of values of the scale parameter of length
    \code{nrow(XX)}}
  \item{g}{a full version of the \code{g} argument}
  \item{theta}{a full version of the \code{theta} argument}
  \item{mle}{if \code{g} and/or \code{theta} is optimized, a
    \code{matrix} containing the values found for these parameters
    and the number of required iterations,
    for each predictive location in \code{XX}}
  \item{eps}{a vector of the jitter values used on the correlation matrix and Q matrix}
  \item{time}{a scalar giving the passage of wall-clock time elapsed
    for (substantive parts of) the calculation}
}
\references{

  D.A. Cole, R.B. Christianson, and R.B. Gramacy (2021).
  \emph{Locally Induced Gaussian Processes for Large-Scale Simulation Experiments}
  Statistics and Computing, 31(3), 1-21; preprint on arXiv:2008.12857;
  \url{https://arxiv.org/abs/2008.12857}

  E. Snelson Z. and Ghahramani. (2006).
  \emph{Sparse Gaussian Processes using Pseudo-inputs}
  Advances in Neural Information Processing Systems 18 , 1257-1264.
}
\author{
    D. Austin Cole \email{austin.cole8@vt.edu}
}

\examples{
## "1D Toy Problem"
## Test function from Forrester et al (2008);
library(hetGP); library(lhs)
X <- matrix(seq(0, 1, length=1000))
Y <- f1d(X)
XX <- matrix(seq(0, 1, length=100))
YY <- f1d(XX)
Xm <- randomLHS(10,1)

out <- giGP(XX=XX, X=X, Y=Y, Xm=Xm, theta=.1)
par(mfrow=c(1,2))
plot(X, Y, type='l', lwd=4, ylim=c(-8, 16))
lines(XX, out$mean, lwd=3, lty=2, col=2)
points(Xm, rep(-8, 10), lwd=2, pch=3, col=3)
legend('topleft', legend=c('Test Function', 'Predicted mean', 'Inducing Points'),
       lty=c(1, 2, NA), col=1:3, pch=c(NA, NA, 3), lwd=2)

plot(XX, YY - out$mean, ylab='Error', type = 'l')

##---------------------------------------------------------------------------##
## a "computer experiment"

## Simple 2-d Herbie's Tooth function used in Cole et al (2020);
## thanks to Lee, Gramacy, Taddy, and others who have used it before
library(lhs)
## Build up a design with N=~40K locations
x <- seq(-2, 2, by=0.05)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)

## Build a inducing point template centered at origin
Xm <- 4*randomLHS(30, 2) - 2

## Predictive grid with N'=400 locations,
xx <- seq(-1.975, 1.975, length=20)
XX <- as.matrix(expand.grid(xx, xx))
YY <- herbtooth(XX)

\donttest{
## Get the predictive equations, first with fixed lengthscale and nugget
out <- giGP(XX=XX, X=X, Y=Y, Xm=Xm, theta=.1)
## RMSE
sqrt(mean((out$mean - YY)^2))


## Refine with optimizing the lengthscale
theta_list <- list(start = .1, min = .05, max = 5)
out2 <- giGP(XX=XX, X=X, Y=Y, Xm=Xm, theta=theta_list)
## RMSE
sqrt(mean((out2$mean - YY)^2))

## Visualize the results
orig_par <- par()
par(mfrow=c(1,3))
image(xx, xx, matrix(out2$mean, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Predicted mean")
image(xx, xx, matrix(out2$mean-YY, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Bias")
image(xx, xx, sqrt(out2$nu) *matrix(sqrt(out$var), nrow=length(xx)),
      col=heat.colors(128), xlab="x1", ylab="x2", main="Stand. Dev.")
points(Xm, pch=3, col=3, lwd=2)
par(orig_par)
}
}
