\name{liGP.forloop}
\alias{liGP.forloop}


\title{
Localized Inducing Point Approximate GP Regression For Many Predictive Locations
}
\description{
  Facilitates locally induced Gaussian process inference and prediction at a large
  set of predictive locations by: building local neighborhoods, shifting an inducing
  point template, optimizing hyperparameters, and calculating GP predictive equations.
}
\usage{
liGP.forloop(XX, X = NULL, Y = NULL, Xm.t, N, g = 1e-6, theta = NULL,
     nu = NULL, epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
     tol = .01, reps = FALSE, Xni.return = FALSE)
}

\arguments{
  \item{XX}{
      a \code{matrix} of out-of-sample
      predictive locations with \code{ncol(XX) = ncol(X)}
  }
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of all input locations. If \code{reps} is a list, this entry is not used.
}
\item{Y}{
  a vector of all responses/dependent values with \code{length(Y)=nrow(X)}. If \code{reps} is a list, this entry is not used.
}
  \item{Xm.t}{
    a \code{matrix} containing the \code{M} inducing points template with \code{ncol(Xm.t) = ncol(X)}. See 'Note' for more.
  }
    \item{N}{
          the positive integer number of nearest neighbor (NN) locations used to build a local neighborhood; \code{N} should be greater than \code{M}. See 'Note' for more.
  }
  \item{g}{
        an initial setting or fixed value for the nugget parameter. In order to optimize g, a list can be provided that includes:
\itemize{
  \item \code{start} -- starting value to initialize the nugget
  \item \code{min} -- minimum value in the allowable range for the nugget
  \item \code{max} -- maximum value in the allowable range for the nugget
  \item \code{ab} -- shape and rate parameters specifying a Gamma prior for the nugget
}
If \code{ab} is not provided, a prior is not placed with the likelihood for optimization. If \code{min} and \code{max} aren't provided, the nugget is not optimized.  A \code{NULL} value generates a list based on \code{garg} in the \pkg{laGP} package. If a single positive scalar is provided, the nugget is fixed for all predictions. Alternatively, a vector of nuggets whose length equals \code{nrow(XX)} can be provided to fix distinct nuggets for each prediction.
  }
  \item{theta}{
        an initial setting or fixed value for the lengthscale parameter. A (default) \code{NULL} value generates an initial setting based on \code{darg} in the \pkg{laGP} package. Similarly, a list can be provided that includes:
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
  \item{epsK}{
      a small positive number added to the diagonal of the correlation \code{matrix} of inducing points
      for numerically stability for inversion. It is automatically increased if neccessary for each prediction.
    }
  \item{epsQ}{
        a small positive number added to the diagonal
        of the Q \code{matrix} (see Cole (2021)) of inducing points
        for numerically stability for inversion. It is automatically increased if neccessary for each prediction.
  }
  \item{tol}{
    a positive number to serve as the tolerance level
    for covergence of the log-likelihood when optimizing
    the hyperparameter(s) theta, g
  }
  \item{reps}{
   a notification of replicate design locations in the data set. If \code{TRUE}, the unique design locations are used for the calculations along with the average response for each unique design location. Alternatively, \code{reps} can be a list from \code{find_reps} in the \pkg{hetGP} package. In this case, \code{X} and \code{Y} are not used.
}
  \item{Xni.return}{
    A scalar logical indicating whether or not a vector of indices into \code{X} (or \code{X0} if a reps list is supplied), specifying the chosen sub-design, should be returned on output.
  }
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
  \item{Xm.t}{the input for \code{Xm.t}}
  \item{eps}{a matrix of \code{epsK} and \code{epsQ} (jitter) values used for each prediction, \code{nrow(eps)=nrow(XX)}}
  \item{mle}{if \code{g} and/or \code{theta} is optimized, a
    \code{matrix} containing the values found for these parameters
    and the number of required iterations,
    for each predictive location in \code{XX}}
    \item{Xni}{
    when Xni.return = TRUE, this field contains a vector of indices of length \code{N} into \code{X} (or \code{X0}) indicating the sub-design (neighborhood) chosen. If \code{nrow(XX)>1}, a matrix is returned with each row matched with the corresponding row of \code{XX}
    }
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

\note{
When selecting the neighborhood size (N) and number of inducing points in
\code{Xm.t}, there is no general rule that works for all problems. However,
for lower dimensions (dim<9) the following values seem to perform well:
N = 100 + 10*dim, M = 10*dim
}
\seealso{
  \code{\link[laGP]{darg}}, \code{\link[laGP]{garg}},
  \code{\link[hetGP]{find_reps}},
  \code{\link[parallel]{makeCluster}}, \code{\link[parallel]{clusterApply}}
}
\examples{
## See LIGP examples
}
