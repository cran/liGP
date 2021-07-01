\name{qnormscale}
\alias{qnormscale}


\title{
Scaling of Inducing Point Design based on Inverse Gaussian CDF
}
\description{
  Scales a set of proposed inducing point locations in [0,1]^d to center
  around a reference location, returning the scaled design
}
\usage{
qnormscale(X, mean, sd)
}

\arguments{
  \item{X}{a \code{matrix} or containing a proposed inducing point design in [0,1]^d
}
\item{mean}{
  a vector representing the reference location to act as the center of the scaling;
  \code{length(mean) = ncol(X)}
}
  \item{sd}{
        a scalar or vector determining the standard deviation for each dimension of
        the Gaussian CDF
  }
}
\details{

  This function scales a set of proposed inducing points in [0,1]^d to be
  centered and concentrated around a reference location. The proposed inducing
  points are interpreted as quantiles of one-dimensional
  Gaussian distributions centered at the reference location with the standard deviation
  provided by the user. For each dimension \code{\link[stats]{qnorm}} is invoked to rescale the
  inducing points.

}
\value{
  a \code{matrix} of the scaled set of inducing points
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
  \code{\link[stats]{qnorm}}
}
\examples{
## Generate data and define xmean
X <- matrix(runif(30), ncol=2)
xmean <- c(0.3, 0.4) # doesn't need to be in [0,1]^2

## Scale centered at xmean with different standard deviations
X_scaled1 <- qnormscale(X, mean=xmean, sd=.1)
X_scaled2 <- qnormscale(X, mean=xmean, sd=c(.05,.15))

## View scaled X
plot(X, xlab='X1', ylab='X2')
points(xmean[1], xmean[2], pch=16)
points(X_scaled1, pch=2, col=3, lwd=2)
points(X_scaled2, pch=3, col=4, lwd=2)
legend('topright',legend = c('Original X','xmean', 'Xscaled1', 'Xscaled2'),
       pch = c(1,16,2,3), col= c(1,1,3,4), lwd=2, lty=NA)

}

