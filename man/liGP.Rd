\name{liGP}
\alias{liGP}


\title{
Localized Inducing Point Approximate GP Regression For Many Predictive Locations
}
\description{
  Facilitates locally induced Gaussian process inference and prediction at a large
  set of predictive locations by: building local neighborhoods, shifting an inducing
  point template, optimizing hyperparameters, and calculating GP predictive equations.
}
\usage{
liGP(XX, X = NULL, Y = NULL, Xm.t, N, g = 1e-6, theta = NULL,
     nu = NULL, num_thread = 1, epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
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
  \item{num_thread}{
      a scalar positive integer indicating the number of threads to use for parallel processing
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
\details{
  When \code{num_threads > 1}, the predictions are performed in parallel using \code{foreach} with clusters created by \pkg{parallel}.

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
## "1D Toy Problem"
## Test function from Forrester et al (2008);
library(hetGP); library(lhs)
X <- matrix(seq(0, 1, length=1000))
Y <- f1d(X)
XX <- matrix(seq(0, 1, length=100))
YY <- f1d(XX)

## Create inducing point template
lhs_design <- randomLHS(9,1)
n <- 80
Xmt <- scale_ipTemplate(X, n, space_fill_design=lhs_design, method='qnorm')$Xm.t

out <- liGP(XX=XX, X=X, Y=Y, Xm=Xmt, N=n, theta=.1)

## Plot predicted mean and error
orig_par <- par()
par(mfrow=c(1,2))
plot(X, Y, type='l', lwd=4, ylim=c(-8, 16),
     main='LIGP fit to Test Function')
lines(XX, out$mean, lwd=3, lty=2, col=2)
legend('topleft', legend=c('Test Function', 'Predicted mean'),
       lty=1:2, col=1:2, lwd=2)

plot(XX, YY - out$mean, xlab='X', ylab='Error', type = 'l',
     main='Predicted Error')
par(orig_par)

##
## Generate new data from function with same mean and non-constant noise
fY <- function(x) { f1d(x) + rnorm(length(x), sd=(1.1 + sin(2*pi*x))) }
Y2 <- fY(X)

## Estimate lengthscale and nugget in predictions
library(laGP)
theta_prior <- darg(NULL, X)
g_prior <- garg(list(mle=TRUE), Y2)
out2 <- liGP(XX=XX, X=X, Y=Y2, Xm=Xmt, N=n, theta=theta_prior,
             g=g_prior, epsK=1e-5)

## Plot predictived mean and confidence intervals
plot(X, Y2, col='grey', cex=.5, pch=16,
     main='LIGP fit to heteroskedastic data', ylab='Y')
lines(X, Y, lwd=2)
lines(XX, out2$mean, lwd=2, lty=2, col=2)
lines(XX, out2$mean + 1.96*sqrt(out2$nu*out2$var), lwd=1, lty=4, col=2)
lines(XX, out2$mean - 1.96*sqrt(out2$nu*out2$var), lwd=1, lty=4, col=2)
legend('topleft', legend=c('Noisy data','Function mean', 'Predicted mean',
                           'Predicted 95 percent confidence interval'),
       lwd=2, lty=c(NA,1,2,3), pch=c(16,NA,NA,NA), col=c('grey',1,2,2))


## View mean and variance errors
par(mfrow=c(1,2))
plot(XX, YY - out2$mean, xlab='X', ylab='Mean Error', type = 'l')
plot(XX, (1.1 + sin(2*pi*XX))^2 - (out2$nu*out2$var),
     xlab='X', ylab='Variance Error', type = 'l')
par(orig_par)

##
## Generate new data with replicates
mults <- sample(2:10, nrow(X), replace=TRUE)
X.reps <- X[rep(1:nrow(X), mults),]
Y.reps <- fY(X.reps)
g_prior <- garg(list(mle=TRUE), Y.reps)

## Generate rep list from hetGP
rep_list <- find_reps(X.reps, Y.reps)

out3 <- liGP(XX=XX, Xm=Xmt, N=n, theta=theta_prior,
             g=g_prior, epsK=1e-5, reps = rep_list)


## Plot predictived mean and confidence intervals
plot(X.reps, Y.reps, col='grey', cex=.5, pch=16,
     main='LIGP fit to data with replicates', xlab='X', ylab='Y')
lines(X, Y, lwd=2)
lines(XX, out3$mean, lwd=2, lty=2, col=2)
lines(XX, out3$mean + 1.96*sqrt(out3$nu*out3$var), lwd=1, lty=4, col=2)
lines(XX, out3$mean - 1.96*sqrt(out3$nu*out3$var), lwd=1, lty=4, col=2)
legend('topleft', legend=c('Noisy data','Function mean', 'Predicted mean',
                           'Predicted 95 percent confidence interval'),
       lwd=2, lty=c(NA,1,2,3), pch=c(16,NA,NA,NA), col=c('grey',1,2,2))

\donttest{
##---------------------------------------------------------------------------##
## a "computer experiment"

## Simple 2-d Herbie's Tooth function used in Cole et al (2021);
## thanks to Lee, Gramacy, Taddy, and others who have used it before

## Build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)

## Build a inducing point template centered at origin
Xm <- matrix(runif(20), ncol=2)
Xmt <- scale_ipTemplate(X=X, N=100, method="chr", space_fill_design = Xm)$Xm.t

## Predictive grid with N'=400 locations
xx <- seq(-1.975, 1.975, length=20)
XX <- as.matrix(expand.grid(xx, xx))
YY <- herbtooth(XX)

## Get the predictive equations, first with fixed lengthscale and nugget
out <- liGP(XX=XX, X=X, Y=Y, Xm.t=Xmt, N=100, Xni.return=TRUE)
## RMSE
sqrt(mean((out$mean - YY)^2))

## View one local neighborhood
xylim <- apply(X[out$Xni[33,],], 2, range)
plot(X[,1], X[,2], pch=16, col='grey', cex=.5,
     xlim=xylim[,1] + c(-.05, .05), ylim=xylim[,2] + c(-.05, .05),
     xlab='X1', ylab='X2')
points(X[out$Xni[33,],1], X[out$Xni[33,],2], pch=16)
points(XX[33,1], XX[33,2], col=3, pch=17, cex=1.5)
points(sweep(Xmt, 2, XX[33,,drop=FALSE], '+'), pch=18, col=2)
legend('topleft', legend=c('Predictive location', 'Data not in neighborhoood',
                          'Neighborhood', 'Inducing points'),
       pch=c(17, 16, 16, 18), col=c(3, 'grey',1,2), cex=1.3)



##
## Refine with optimizing the lengthscale
theta_list <- darg(NULL, X)
out2 <- liGP(XX=XX, X=X, Y=Y, Xm.t=Xmt, N=100, theta=theta_prior)

## RMSE
sqrt(mean((out2$mean - YY)^2))

## Visualize the results
par(mfrow=c(1,3))
image(xx, xx, matrix(out2$mean, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Predictive Mean")
image(xx, xx, matrix(out2$mean-YY, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Bias")
image(xx, xx, matrix(sqrt(out2$nu*out$var), nrow=length(xx)),
      col=heat.colors(128), xlab="x1", ylab="x2", main="Stand. Dev.")
par(orig_par)

##
## Predictions from noisy training data with replicates
Xreps <- X[rep(1:nrow(X), 5),]
Ynoisy <- herbtooth(Xreps) + rnorm(nrow(Xreps), sd=.02)

library(hetGP)
reps_list <- find_reps(Xreps, Ynoisy)

## Priors for theta and g
theta_prior <- darg(NULL, Xreps)
g_prior <- garg(list(mle=TRUE), Ynoisy)

## Predictions with estimated nugget
out_noisydata <- liGP(XX, Xm.t = Xmt, N = 100, g=g_prior, theta=theta_prior,
                      reps = reps_list)
## RMSE
estimated_noise <- sqrt(out_noisydata$mle[,1]*out_noisydata$nu)
sqrt(mean((estimated_noise - .02)^2))


## Visualize the results
par(mfrow=c(1,3))
image(xx, xx, matrix(out_noisydata$mean, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Predictive Mean")
image(xx, xx, matrix(out_noisydata$mean-YY, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="Bias")
image(xx, xx, matrix(sqrt(out_noisydata$nu*out_noisydata$var) - .02, nrow=length(xx)),
      col=heat.colors(128), xlab="x1", ylab="x2", main="Stand. Dev. Error")
par(orig_par)
}

}
