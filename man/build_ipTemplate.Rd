\name{build_ipTemplate}
\alias{build_ipTemplate}


\title{
Inducing point template design built through sequential optimization
}

\description{
Constructs a design of inducing points around the center of the design matrix and its local neighborhood. The output is an inducing point design centered at the origin that can be used as a template for predictions anywhere in the design space (with a local neighborhood of the same size). Different criteria are available to optimize the inducing points. The methods "\code{wimse}" and "\code{alc}" use weighted Integrated Mean Squared Error and Active Learning Cohn respectively to sequentially select inducing points. }

\usage{
build_ipTemplate(X = NULL, Y = NULL, M, N, theta = NULL, g = 1e-4,
                 method = c('wimse','alc'), ip_bounds = NULL,
                 integral_bounds = NULL, num_thread = 1, num_multistart = 20,
                 w_var = NULL, epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                 reps = FALSE, verbose = TRUE)
}
\arguments{
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations. If using a list for \code{reps}, this entry is not used.
}
\item{Y}{
  a vector of responses/dependent values with \code{length(Y)=nrow(X)}. If using a list for \code{reps}, this entry is not used.
}
\item{M}{
    a positive integer number of inducing points; \code{M} should be less than \code{N}
}
\item{N}{
        a positive integer number of Nearest Neighbor (NN) locations used to build a local neighborhood
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
\item{method}{
  specifies the method by which the inducing point template is built. In brief,
  wIMSE (\code{"wimse"}, default) minimizes the weighted integrated predictive variance
  and ALC (\code{"alc"}) minimizes predictive variance
}
\item{ip_bounds}{a 2 by d \code{matrix} of the bounds used in the optimization of inducing points;
the first row contains minimum values, the second row the maximum values;
if not provided, the bounds of the center's local neighborhood are used}
\item{integral_bounds}{
    a 2 by d \code{matrix} of the bounds used in the calculation of \code{wimse}; the first row contains minimum values, the second row the maximum values; only relevant when \code{method="wimse"}; if not provided, defaults to the range of each column of \code{X}
}
\item{num_thread}{
    a scalar positive integer indicating the number of GPUs
    available for calculating ALC; only relevant when \code{method="alc"}
  }
\item{num_multistart}{
    a scalar positive integer indicating the number of multistart points
    used to optimize each inducing point with wIMSE or ALC
  }
  \item{w_var}{
    a scalar positive number used as the variance for the Gaussian weight in wIMSE. If \code{NULL}, \code{theta} is used.
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
  This function calls separate subroutines for certain methods. For \code{method=wimse}, the function
  \code{\link{optIP.wIMSE}} is called with the reference point being the median of the design matrix. If
  \code{method=alc}, \code{\link{optIP.ALC}} is called with the predictive variance being minimized at the
  median of the design matrix. For any inducing point design, the first inducing point is placed at the predictive location (i.e. the
  origin).
}
\value{
  The output is a \code{list} with the following components.

\item{Xm.t}{a \code{matrix} of \code{M} inducing points centered at the origin}
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
  \code{\link[liGP]{optIP.wIMSE}},\code{\link[liGP]{optIP.ALC}}
}

\examples{

## "1D Toy Problem"
## Test function from Forrester et al (2008);
library(hetGP)
X <- as.matrix(seq(0, 1, length=1000))
Y <- f1d(X)
int_bounds <- matrix(c(0, 1))

## Center of design space used to build inducing point templates
X_center <- median(X)

\donttest{
## Optimize inducing points with weighted Integrated Mean-Square Error
wimse.out <- build_ipTemplate(X, Y, N=100, M=10, method="wimse", integral_bounds=int_bounds)
Xm.t_wimse <- wimse.out$Xm.t

## now optimize inducing points using Active Learning Cohn
alc.out <- build_ipTemplate(X, Y, N=100, M=10, method="alc", integral_bounds=int_bounds)
Xm.t_alc <- alc.out$Xm.t
Xn <- alc.out$Xn ## X_center neighborhood

## plot locally optimized inducing point templates
plot(X, Y, pch=16, cex=.5, col='grey')
points(Xn, f1d(Xn), col=2)
points(Xm.t_wimse + X_center, rep(-4, 10), pch=2, col=3)
points(Xm.t_alc + X_center, rep(-5, 10), pch=3, col=4)
legend('topleft', pch = c(16, 16, 2, 3), col = c('grey', 2, 3, 4),
       legend=c('Data', 'Local neighborhood', 'wIMSE inducing point design',
                'ALC inducing point design'))
}

## "2D Toy Problem"
## Herbie's Tooth function used in Cole et al (2020);
## thanks to Lee, Gramacy, Taddy, and others who have used it before

## build up a design with N=~40K locations
x <- seq(-2, 2, by=0.02)
X <- as.matrix(expand.grid(x, x))
Y <- herbtooth(X)
X_center <- apply(X, 2, median)
\donttest{
## build a inducing point template, first with weighted Integrated Mean-Square Error
int_bounds <- rbind(c(-2,-2), c(2,2))
wimse.out <- build_ipTemplate(X, Y, N=100, M=10, method="wimse", integral_bounds=int_bounds)
Xm.t_wimse <- wimse.out$Xm.t

## now optimize inducing points using Active Learning Cohn
alc.out <- build_ipTemplate(X, Y, N=100, M=10, method="alc", integral_bounds=int_bounds)
Xm.t_alc <- alc.out$Xm.t
Xn <- alc.out$Xn

## plot locally optimized inducing point templates
plot(Xn, pch=16, cex=.5, col='grey',
     xlab = 'x1', ylab = 'x2', main='Locally optimized IP templates')
points(X_center[1], X_center[2], pch=16, cex=1.5)
points(Xm.t_wimse, pch=2, lwd=2, col=3)
points(Xm.t_alc, pch =3, lwd=2, col=4)
legend('topleft', pch = c(16, 2, 3), col = c('grey', 3, 4),
       legend=c('Local neighborhood', 'wIMSE inducing point design',
                'ALC inducing point design'))
}
}



