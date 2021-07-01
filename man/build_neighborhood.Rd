\name{build_neighborhood}
\alias{build_neighborhood}


\title{
Nearest Neighbor (NN) data subset given a center
}

\description{
Constructs a neighborhood of points that are a subset of the data for a given center (i.e. predictive) location. }

\usage{
build_neighborhood(N, xx = NULL, X = NULL, Y = NULL, reps_list = NULL)
}
\arguments{
\item{N}{
        the positive integer number of Nearest Neighbor (NN) locations used to build a local neighborhood
}
\item{xx}{a row matrix of the location of the neighborhood's center. If \code{NULL}, the
center (median) of the data is used}
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations. If \code{reps_list} is supplied, this entry is not used.
}
\item{Y}{
  a vector of responses/dependent values with \code{length(Y)=nrow(X)}. If \code{reps_list} is supplied, this entry is not used.
}

\item{reps_list}{
   a list from \code{find_reps} in the \code{hetGP} package, that includes the entries \code{X0} and \code{Z0}. In this case, \code{X} and \code{Y} are not used.
}

}
\details{
  This function builds a local neighborhood around the center \code{xx}. If \code{X} is supplied, the \code{N} NN points are found and chosen. If \code{reps_list} is supplied, \code{N} unique data locations \code{X0} are supplied, along with their averaged responses (\code{Z0}) and original responses (\code{Zlist}).
}
\value{
  The output is a \code{list} with the following components:

\item{xx}{a row \code{matrix} of the neighborhood's center}
If \code{reps_list=NULL},
\item{Xn}{a \code{matrix} of the local neighborhood's design points.}
\item{Yn}{a \code{vector} of the local neighborhood's responses. Only provided when \code{Y} is provided.}
If \code{reps_list} is provided,
\item{Xn0}{a \code{matrix} of the local neighborhood's unique design locations.}
\item{Yn0}{a \code{vector} of averages observations at \code{Xn0}.}
\item{mult}{a vector of the number of replicates at \code{Xn0}.}
\item{Yn_list}{a \code{list} where each element corresponds to observations at a design in \code{Xn0}.}
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
## "2D Toy Problem"
## Herbie's Tooth function used in Cole et al (2021);
## thanks to Lee, Gramacy, Taddy, and others who have used it before
library(hetGP)

## Build data with replicates
x <- seq(-2, 2, by=0.05)
X <- as.matrix(expand.grid(x, x))
X <- rbind(X, X)
Y <- herbtooth(X) + rnorm(nrow(X), sd = .02)
reps_list <- find_reps(X, Y)
xx <- matrix(c(-0.12, 1.53), nrow=1)

## Build neighborhoods
neighborhood1 <- build_neighborhood(N=100, xx=xx, X=X, Y=Y)
neighborhood2 <- build_neighborhood(N=100, xx=xx, reps_list=reps_list)

## Compare neighborhood sizes
Xn0_range <- apply(neighborhood2$Xn0, 2, range)
plot(X, xlim = Xn0_range[,1] + c(-.15, .15), ylim = Xn0_range[,2] + c(-.1, .25),
     pch=3)
points(neighborhood2$Xn0, pch=16, col='grey')
points(neighborhood1$Xn, col=2, lwd=2)
points(xx, pch=17, col=3, cex=1.5)
legend('topleft', ncol=2, pch=c(3, 17, 16, 1), col=c(1, 3, 'grey', 2),
       legend=c('Design locations', 'Neighborhood center',
                'Xn based on unique locations', 'Xn ignoring unique locations'))
}



