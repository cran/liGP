\name{loiGP}
\alias{loiGP}


\title{
Locally Optimized Inducing Point Approximate GP Regression For Many Predictive Locations
}
\description{
  Facilitates localized Gaussian process inference and prediction at a large
  set of predictive locations, by opimizing a local set of inducing points
  for each predictive location's local neighborhood and then calling \code{\link{giGP}}.
}
\usage{
loiGP(XX, X = NULL, Y = NULL, M, N, g = 1e-6, theta = NULL, nu = NULL,
      method = c('wimse','alc'), integral_bounds = NULL, num_thread = 1,
      epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01, reps = FALSE)
}
\arguments{
   \item{XX}{
    a \code{matrix} of out-of-sample
    predictive locations with \code{ncol(XX) = ncol(X)}; \code{loiGP} calls \code{giGP} for
    each row of \code{XX}, independently
  }
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of all input locations. If \code{reps} is a list, this entry is not used.
}
\item{Y}{
  a vector of all responses/dependent values with \code{length(Y)=nrow(X)}. If \code{reps} is a list, this entry is not used.
}
  \item{M}{
        the positive integer number of inducing points placed for each local neighborhood; \code{M} should be less than \code{N}
  }
  \item{N}{
        the positive integer number of Nearest Neighbor (NN) locations used to build a local neighborhood
  }
 \item{g}{
        an initial setting or fixed value for the nugget parameter. In order to optimize g, a list can be provided that includes:
\itemize{
  \item \code{start} -- starting value to initialize the nugget
  \item \code{min} -- minimum value in the allowable range for the nugget
  \item \code{max} -- maximum value in the allowable range for the nugget
  \item \code{ab} -- shape and rate parameters specifying a Gamma prior for the nugget
}
If \code{ab} is not provided, a prior is not placed with the likelihood for optimization. If \code{min} and \code{max} aren't provided, the nugget is not optimized. A  \code{NULL} value generates an initial setting based on \code{garg} in the \code{laGP} package. If a single positive scalar is provided, the nugget is fixed for all predictions. Alternatively, a vector of nuggets whose length equals \code{nrow(XX)} can be provided to fix distinct nuggets for each prediction.
  }
  \item{theta}{
        an initial setting or fixed value for the lengthscale parameter. A (default) \code{NULL} value generates an initial setting based on \code{darg} in the \code{laGP} package. Similarly, a list can be provided that includes:
\itemize{
  \item \code{start} -- starting value to initialize the lengthscale
  \item \code{min} -- minimum value in the allowable range for the lengthscale
  \item \code{max} -- maximum value in the allowable range for the lengthscale
  \item \code{ab} -- shape and rate parameters specifying a Gamma prior for the lengthscale
}
If \code{ab} is not provided, a prior is not placed with the likelihood for optimization. If \code{min} and \code{max} aren't provided, the lengthscale is not optimized. If a single positive scalar is provided, the lengthscale is fixed for all predictions. Alternatively, a vector of lengthscales whose length equals \code{nrow(XX)} can be provided to fix distinct lengthscales for each prediction.
  }
  \item{nu}{
    a positive number used to set the scale parameter;
    default (\code{NULL}) calculates the maximum likelihood estimator
  }
  \item{method}{
    specifies the method by which the inducing point template is built. In brief,
    wIMSE (\code{"wimse"}, default) minimizes the weighted integrated mean-sqaure error,
    and ALC (\code{"alc"}) minimizes predictive variance at the preditive location.
  }
  \item{integral_bounds}{
    a 2 by d \code{matrix} of the domain bounds of the data (used in the calculation of \code{wimse});
    the first row contains
    minimum values, the second row the maximum values; only relevant when \code{method="wimse"}; if not
    provided, defaults to the range of each column of \code{X}
  }
  \item{num_thread}{
      a scalar positive integer indicating the number of threads to use for parallel processing
    }
  \item{epsK}{
      a small positive number added to the diagonal of the correlation matrix, of inducing points, K,
      for numerically stability for inversion. It is automatically increased if neccessary for each prediction.
    }
  \item{epsQ}{
        a small positive number added to the diagonal of the Q \code{matrix} (see Cole (2021))
        for numerically stability for inversion. It is automatically increased if neccessary for each prediction.
  }
  \item{tol}{
    a positive number to serve as the tolerance level
    for covergence of the log-likelihood when optimizing
    the hyperparameter(s) theta and/or g
  }
  \item{reps}{
   a notification of replicate design locations in the data set. If \code{TRUE}, the unique design locations are used for the calculations along with the average response for each unique design location. Alternatively, \code{reps} can be a list from \code{find_reps} in the \code{hetGP} package. In this case, \code{X} and \code{Y} are not used.
}
}
\details{
  This function builds a unique inducing point design to accompany the local neighborhood for each preditive location in \code{XX}. It then invokes \code{\link{giGP}} for each row of \code{XX} with \code{X=Xn, Y=Yn} from the corresponding local neighborhood and locally optimial inducing point design. For further information, see \code{\link{giGP}}.
}
\value{
  The output is a \code{list} with the following components:

 \item{mean}{a vector of predictive means of length \code{nrow(XX)}}
 \item{var}{a vector of predictive variances of length
    \code{nrow(XX)}}
 \item{nu}{a vector of values of the scale parameter of length
    \code{nrow(XX)}}
 \item{g}{a full version of the \code{g} argument}
 \item{theta}{a full version of the \code{theta} argument}
\item{Xm}{a \code{list} of inducing point designs; each entry in the list is a
  \code{matrix} containing  \code{M} locally optimized inducing points; \code{length(Xm)=nrow(XX)}}
\item{eps}{a matrix of \code{epsK} and \code{epsQ} (jitter) values used for each prediction, \code{nrow(eps)=nrow(XX)}}
 \item{mle}{if \code{g} and/or \code{theta} is optimized, a
    \code{matrix} containing the values found for these parameters
    and the number of required iterations,
    for each predictive location in \code{XX}}
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
library(hetGP); library(lhs)
X <- matrix(seq(0, 1, length=1000))
Y <- f1d(X)
XX <- matrix(seq(.01, .99, length=50))
YY <- f1d(XX)


n <- 80
m <- 10
int_bounds <- matrix(c(0,1))
\donttest{
out <- loiGP(XX=XX, X=X, Y=Y, M=m, N=n, method='wimse',
             integral_bounds=int_bounds)

## Plot predicted mean and error
orig_par <- par()
par(mfrow=c(1,2))
plot(X, Y, type='l', lwd=4, ylim=c(-8, 16))
lines(XX, out$mean, lwd=3, lty=2, col=2)
legend('topleft', legend=c('Test Function', 'Predicted mean'),
       lty=1:2, col=1:2, lwd=2)

plot(XX, YY - out$mean, xlab='X', ylab='Error', type = 'l')
par(orig_par)
}
}
