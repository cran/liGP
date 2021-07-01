\name{scale_ipTemplate}
\alias{scale_ipTemplate}


\title{
Inducing points design scaling for a local neighborhood template
}

\description{
Scales a design of inducing points around the center of the design matrix and its local neighborhood. The output is an inducing points design centered at the origin that can be used as a template for predictions anywhere in the design space (with a local neighborhood of the same size).  Method include scaling by a circumscribed hyperrectangle (\code{chr}) and an inverse Gaussian CDF (\code{qnorm}).}

\usage{
scale_ipTemplate(X, N, space_fill_design, method = c('qnorm','chr'))
}
\arguments{
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations
}
  \item{N}{
        the positive integer number of Nearest Neighbor (NN) locations used to build a local neighborhood
}
\item{space_fill_design}{a \code{matrix} in [0,1]^d with \code{M} rows and number of columns = \code{ncol(X)} that is scaled and centered to create the inducing points template
}
\item{method}{
  the method by which the inducing point template is scaled. In brief,
  cHR (\code{"chr"}) scales \code{space_fill_design} to circumscribe the neighborhood and
  qNorm (\code{"qnorm"}) scales \code{space_fill_design} by the inverse Gaussian CDF.
}

}
\details{
  This function calls separate subroutines for certain methods. When \code{method=qnorm}, \code{\link{qnormscale}} is called. The mean of the Gaussian distribution is the
  median of the design matrix. The standard deviation of the Gaussian distribution is one-third of the
  maximum distance from the median of the design matrix to the neighborhood points for each dimension.

  For each inducing point design, the origin (i.e. predictive location) is appended to the scaled inducing point design. Thus, the resulting design contains \code{M+1} inducing points.
}
\value{
  The output is a \code{list} with the following components:

\item{Xm.t}{a \code{matrix} of \code{M+1} inducing points centered at the origin}
\item{Xn}{a \code{matrix} of the local neighborhood at the center of the design}
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


\seealso{
  \code{\link[liGP]{qnormscale}}
}

\examples{
## Build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))

## Create inducing point template, first with
## circumscribed hyperrectangle (cHR) scaling
M = 10
Xm <- matrix(runif(2*M), ncol=2)
chr_temp <- scale_ipTemplate(X, N=100, space_fill_design=Xm, method="chr")
Xm.t_chr <- chr_temp$Xm.t
Xn <- chr_temp$Xn

## Now create template with Inverse Gaussian CDF scaling
qnorm_temp <- scale_ipTemplate(X, N=100, space_fill_design=Xm, method="qnorm")
Xm.t_qnorm <- qnorm_temp$Xm.t

\donttest{
## View different scaled template designs
X_center <- apply(X, 2, median)
ranges <- apply(Xn, 2, range)
plot(Xn[,1], Xn[,2], pch=16, cex=.5, col='grey',
     xlim=ranges[,1], ylim=ranges[,2]+c(0,.1),
     xlab = 'x1', ylab = 'x2', main='Scaled Inducing Point templates')
points(X_center[1],X_center[2], pch=16)
points(Xm.t_chr, col=3, pch=2, lwd=2)
points(Xm.t_qnorm, col=4, pch=3, lwd=2)
legend('topleft', pch=c(16,16,2,3), lty=NA, lwd=2, col=c(1,'grey',3,4), ncol=2,
       legend=c('Xcenter', 'Neighborhood','cHR IP template','qNorm IP template'))
}
}



