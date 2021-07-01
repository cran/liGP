\name{liGP_gauss_measure}
\alias{liGP_gauss_measure}


\title{
Localized Inducing Point Approximate GP Regression For a Gaussian Measure
}
\description{
  Facilitates locally induced Gaussian process inference and prediction across a Gaussian measure by:
  building one local neighborhood around the measure, shifting an inducing
  point template, optimizing hyperparameters, calculating GP mean predictions, and estimating the
  integral through a discrete set or quadrature.
}
\usage{
liGP_gauss_measure(xstar, X, Y, Xm.t, N, gauss_sd, measure_bounds = c(-Inf, Inf),
      g = 1e-6, epsi = NULL, epsK = 1e-6, epsQ = 1e-5, seq_length = 20)
}

\arguments{
  \item{xstar}{
      a one-row \code{matrix} of the mean of the Gaussian measure.
  }
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of all input locations.
}
\item{Y}{
  a vector of all responses/dependent values with \code{length(Y)=nrow(X)}.
}
  \item{Xm.t}{
    a \code{matrix} containing the \code{M} inducing points template with \code{ncol(Xmt) = ncol(X)}.
  }
    \item{N}{
          the positive integer number of nearest neighbor (NN) locations used to build a local neighborhood; \code{N} should be greater than \code{M}
  }
  \item{gauss_sd}{
    a vector of standard deviations for the Gaussian measure with with \code{length(gauss_sd)=nrow(X)}. Note: at this time, the Gaussian measure must only have one nonzero standard deviation (i.e. the Gaussian measure is a slice).
  }
    \item{measure_bounds}{
    a vector of the bounds of the Gaussian measure for the single dimension with a nonzero standard deviation. This is only used if \code{epsi} is \code{NULL}.
  }
  \item{g}{
        a fixed value for the nugget parameter.
  }
  \item{epsi}{
      an optional vector of Gaussian noise drawn from gauss_sd used with \code{xstar} to generate a set of predictive locations for estimating the integral. If not provided, the \code{integrate} function is called to perform to estimate the integral.
  }

  \item{epsK}{
      a small positive number added to the diagonal of the correlation \code{matrix} of inducing points
      for numerically stability for inversion
    }
  \item{epsQ}{
        a small positive number added to the diagonal
        of the Q \code{matrix} (see Cole (2021)) of inducing points
        for numerically stability for inversion
  }
  \item{seq_length}{
    a positive integer used to build sequences of this length in the nondegenerate dimension for the purpose of building a local neighbhorhood. This sequence is not used in prediction.
  }
}
\details{
  This function is built to deal with the special class of problems where liGP is used to predict and integrate over a degenerate Gaussian measure where only one dimension has a nonzero standard deviation.


}
\value{
  the pointwise estimate for the mean prediction over the Gaussian measure
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


\seealso{
  \code{\link[laGP]{darg}},
  \code{\link[stats]{integrate}}
}
\examples{

##---------------------------------------------------------------------------##
## a "computer experiment"

## Simple 2-d Herbie's Tooth function used in Cole et al (2020);
## thanks to Lee, Gramacy, Taddy, and others who have used it before

## Build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)

## Build a inducing point template centered at origin
X_m <- matrix(runif(20), ncol=2)
Xmt <- scale_ipTemplate(X, N=100, space_fill_design=X_m, method="qnorm")$Xm.t


## predictive center
xx <- matrix(c(.5, .5), ncol=2)

## Standard deviation of gaussian measure with random draws
gauss_sd <- c(0, .1)
epsi <- rnorm(30, sd = gauss_sd[2])

## Get the predictive equations, first with fixed lengthscale and nugget
out <- liGP_gauss_measure(xx, X=X, Y=Y, Xm.t=Xmt, N=100,
                          gauss_sd=gauss_sd, epsi=epsi)


## Refine with using integrate function
out2 <- liGP_gauss_measure(xx, X=X, Y=Y, Xm.t=Xmt, N=100, gauss_sd=gauss_sd)


}
