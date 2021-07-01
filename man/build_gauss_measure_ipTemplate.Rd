\name{build_gauss_measure_ipTemplate}
\alias{build_gauss_measure_ipTemplate}


\title{
Inducing point template design for a Gaussian measure built through sequential optimization
}

\description{
Constructs a design of inducing points around a Gaussian measure whose mean is the center of the design matrix and its local neighborhood. The output is an inducing point design centered at the origin that can be used as a template for predictions anywhere in the design space (with a local neighborhood of the same size). The inducing points are sequentially selected by optimizing "\code{wimse}", weighted Integrated Mean Squared Error. }

\usage{
build_gauss_measure_ipTemplate(X = NULL, Y = NULL, M, N, gauss_sd,
                               theta = NULL, g = 1e-4, seq_length=20,
                               ip_bounds = NULL, integral_bounds = NULL,
                               num_multistart = 20,
                               epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                               reps = FALSE, verbose = TRUE)
}
\arguments{
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations. If using a list for \code{reps}, this entry is not used
}
\item{Y}{
  a vector of responses/dependent values with \code{length(Y)=nrow(X)}. If using a list for \code{reps}, this entry is not used
}
\item{M}{
    a positive integer number of inducing points; \code{M} should be less than \code{N}
}
\item{N}{
        the positive integer number of Nearest Neighbor (NN) locations used to build a local neighborhood
}
\item{gauss_sd}{
      a vector of standard deviations for the Gaussian measure with \code{length(gauss_sd)=nrow(X)}. Note: at this time, the Gaussian measure must only have one nonzero standard deviation (i.e. the Gaussian measure is a slice)
}
\item{theta}{
      the lengthscale parameter (positive number) in a Gaussian
      correlation function; a (default) \code{NULL} value sets
      the lengthscale at the square of the 10th percentile of
    pairwise distances between neighborhood points
    (similar to \code{\link{darg}}
    in \code{\link{laGP}} package)
}
\item{g}{
  the nugget parameter (positive number) in a covariance
}
\item{seq_length}{
   a positive integer used to build sequences of this length in the nondegenerate dimensions for the purpose of building a local neighbhorhood.
}
\item{ip_bounds}{a 2 by d \code{matrix} of the bounds used in the optimization of inducing points;
the first row contains minimum values, the second row the maximum values;
if not provided, the bounds of the center's local neighborhood are used}
\item{integral_bounds}{
    a 2 by d \code{matrix} of the bounds used in the calculation of \code{wimse}; the first row contains minimum values, the second row the maximum values; only relevant when \code{method="wimse"}; if not provided, defaults to the range of each column of \code{X}
}
\item{num_multistart}{
    a scalar positive integer indicating the number of multistart points
    used to optimize each inducing point
  }
  \item{epsK}{
    a small positive number added to the diagonal of the correlation \code{matrix} of inducing points
    for numerically stability for inversion
  }
\item{epsQ}{
      a small positive number added to the diagonal of the Q \code{matrix} (see Cole (2021))
    for numerically stability for inversion
  }
\item{reps}{
   a notification of replicate design locations in the data set. If \code{TRUE}, the unique design locations are used for the calculations along with the average response for each unique design location. Alternatively, \code{reps} can be a list from \code{find_reps} in the \code{hetGP} package. In this case, \code{X} and \code{Y} are not used.
}
\item{verbose}{
    when \code{TRUE}, prints the current number of inducing
    points selected during the sequential optimization process
}
}
\details{
    This function is built to deal with the special class of problems where liGP is used to predict and integrate over a degenerate Gaussian measure where only one dimension has a nonzero standard deviation. To build the wIMSE inducing point design, the function
  \code{\link{optIP.wIMSE}} is called with the reference point being the median of the design matrix.

  For each inducing point design, the first inducing point is placed at the predictive location (i.e. the
  origin).
}
\value{
  The output is a \code{list} with the following components.

\item{Xm.t}{a \code{matrix} of \code{M} inducing points centered at the origin}
\item{Xn}{a \code{matrix} of the local neighborhood at the center of the design}
\item{Xc}{a \code{matrix} of the center of the design used to build the local neighborhood and inducing point template}
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
  \code{\link[liGP]{optIP.wIMSE}}}


\examples{
## "2D Toy Problem"
## Herbie's Tooth function used in Cole et al (2020);
## thanks to Lee, Gramacy, Taddy, and others who have used it before

## build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)
X_center <- apply(X, 2, median)
gauss_sd <- c(0, .05)
\donttest{
## build a inducing point template, first with original weighted Integrated Mean-Square Error
int_bounds <- rbind(c(-2,-2), c(2,2))
wimse.out <- build_ipTemplate(X, Y, N=100, M=10, method='wimse',
                              integral_bounds=int_bounds)
Xm.t_wimse <- wimse.out$Xm.t
Xn <- wimse.out$Xn


wimse_gauss.out <- build_gauss_measure_ipTemplate(X, Y, N=100, M=10,
                                                  gauss_sd = gauss_sd,
                                                  integral_bounds=int_bounds)
Xm.t_wimse_gauss <- wimse_gauss.out$Xm.t
Xn_gauss <- wimse_gauss.out$Xn

## plot locally optimized inducing point templates
ylim <- range(Xn_gauss[,2]) + c(-.03, .05)
plot(Xn, pch=16, cex=.5, col='grey',
     xlab = 'x1', ylab = 'x2', ylim = ylim,
     main='Locally optimized IP template based on Gaussian measure')
points(Xn_gauss, cex=.7)
points(X_center[1], X_center[2], pch=16, cex=1.5)
points(Xm.t_wimse, pch=2, lwd=2, col=3)
points(Xm.t_wimse_gauss, pch=6, lwd=2, col=2)
legend('topleft', pch = c(16, 1, 2, 3), col = c('grey', 1, 3, 2),
       legend=c('Local neighborhood (wIMSE)',
                'Local neighborhood (Gauss measure)',
                'wIMSE ip design',
                'Gaussian measure ip design'))
}
}



