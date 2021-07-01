\name{herbtooth}
\alias{herbtooth}


\title{
Herbie's Tooth function
}
\description{
  Two-dimensional function whose surface resembles a molar, with multiple local minima/maxima
}

\arguments{
  \item{X}{a \code{matrix} or \code{data.frame} containing
    the full (large) design matrix of input locations in [-2,2]^2
}
}
\details{

  A non-stationary function with many local minima/maxima that is difficult to model with a global
  model.
}
\value{
    a vector of evaluations of the Herbie's tooth function, length \code{= nrow(X)}

}
\references{

  H.K.H. Lee, R.B. Gramacy, C. Linkletter, and G. Gray. 2011.
  \emph{Optimization Subject to Hidden Constraints via Statistical Emulation}
  Pacific Journal of Optimization 7 (3): 467-78.

}
\author{
  D. Austin Cole \email{austin.cole8@vt.edu}
}



\examples{
X <- matrix(runif(200, min = -2, max = 2), ncol=2)
herbtooth(X)
}

