\name{scale_gauss_measure_ipTemplate}
\alias{scale_gauss_measure_ipTemplate}


\title{
Inducing points design scaling for a Gaussian measure local neighborhood template
}

\description{
Scales a design of inducing points around a Gaussian measure whose mean is the center of the design matrix and its local neighborhood. The output is an inducing points design centered at the origin that can be used as a template for predictions anywhere in the design space (with a local neighborhood of the same size).  Method include scaling by a circumscribed hyperrectangle (\code{chr}) and an inverse Gaussian CDF (\code{qnorm}).}

\usage{
scale_gauss_measure_ipTemplate(X, N, gauss_sd, space_fill_design,
                               method = c('qnorm','chr'), seq_length=20)
}
\arguments{
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations
}
  \item{N}{
        the positive integer number of Nearest Neighbor (NN) locations used to
        build a local neighborhood
}
\item{gauss_sd}{
      a vector of standard deviations for the Gaussian measure with with
      \code{length(gauss_sd)=nrow(X)}. Note: at this time, the Gaussian measure
      must only have one nonzero standard deviation (i.e. the Gaussian measure
      is a slice).
}
\item{space_fill_design}{a \code{matrix} in [0,1]^d with \code{M} rows and
number of columns = \code{ncol(X)} that is scaled and centered to create the
inducing points template
}
\item{method}{
  the method by which the inducing point template is scaled. In brief,
  cHR (\code{"chr"}) scales \code{space_fill_design} to circumscribe the
  neighborhood and qNorm (\code{"qnorm"}) scales \code{space_fill_design} by
  the inverse Gaussian CDF.
}
\item{seq_length}{
  an integer that defines the sequence length used to represent the gaussian
  measure when building the neighbhorhood.
}
}
\details{
  This function is built to deal with the special class of problems where liGP
  is used to predict and integrate over a degenerate Gaussian measure where only
  one dimension has a nonzero standard deviation. Separate subroutines are
  called for different methods. When \code{method=qnorm},
  \code{\link{qnormscale}} is called. The mean of the Gaussian distribution is
  the median of the design matrix. The standard deviation of the Gaussian
  distribution is one-third of the maximum distance from the median of the
  design matrix to the neighborhood points for each dimension.

  For each inducing point design, the origin (i.e. predictive location) is
  appended to the scaled inducing point design. Thus, the resulting design
  contains \code{M+1} inducing points.
}
\value{
  The output is a \code{list} with the following components.

\item{Xm.t}{a \code{matrix} of \code{M+1} inducing points centered at the origin}
\item{Xn}{a \code{matrix} of the local neighborhood at the center of the design}
\item{Xc}{a \code{matrix} of the center of the design used to build the local
neighborhood and inducing point template}
\item{gauss_sd}{the \code{gauss_sd} used to generate the local neighborhood}
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
X_center <- apply(X, 2, median)
## Create inducing point template, first with
## circumscribed hyperrectangle (cHR) scaling
M = 10
Xm <- matrix(runif(2*M), ncol=2)

## Create template with Inverse Gaussian CDF scaling
qnorm_temp <- scale_ipTemplate(X, N=100, space_fill_design=Xm, method="qnorm")
Xm.t_qnorm <- qnorm_temp$Xm.t
Xn <- qnorm_temp$Xn


## Create template with Inverse Gaussian CDF scaling
gauss_sd <-  c(0, .05)
qnorm_temp_gauss <- scale_gauss_measure_ipTemplate(X, N=100, gauss_sd=gauss_sd,
                                             space_fill_design=Xm,
                                             method="qnorm")
Xm.t_qnorm_gauss <- qnorm_temp_gauss$Xm.t
Xn_gauss <- qnorm_temp_gauss$Xn


\donttest{
## View different scaled template designs
ylim <- range(Xn_gauss[,2]) + c(-.03, .05)
plot(Xn, pch=16, cex=.5, col='grey',
     xlab = 'x1', ylab = 'x2', ylim = ylim,
     main='Locally optimized IP template based on Gaussian measure')
points(Xn_gauss, cex=.7)
points(X_center[1], X_center[2], pch=16, cex=1.5)
points(Xm.t_qnorm, pch=2, lwd=2, col=3)
points(Xm.t_qnorm_gauss, pch=6, lwd=2, col=2)
legend('topleft', pch = c(16, 1, 2, 3), col = c('grey', 1, 3, 2),
       legend=c('Local neighborhood (qNorm)',
                'Local neighborhood (Gauss measure)',
                'qnorm ip design',
                'Gaussian measure ip design'))
}
}



