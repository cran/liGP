\name{borehole}
\alias{borehole}


\title{
Borehole equation data generator
}
\description{
  Generates a vector of outputs from the borehole function, a common testing function for computer experiments. The function models water flow through a borehole.
}
\usage{
borehole(X)
}

\arguments{
  \item{X}{a \code{matrix} containing
    the full (large) design matrix of input locations in [0,1]^8
}
}
\details{

  For more details, see Worley, B. A. (1987). Deterministic uncertainty analysis (No. CONF-871101-30). Oak Ridge National Lab., TN (USA)

}
\value{
  a vector of evaluations of the borehole function, length \code{= nrow(X)}

}
\references{

  Harper, W. V., & Gupta, S. K. (1983).
  \emph{Sensitivity/uncertainty analysis of a borehole scenario comparing Latin Hypercube Sampling and deterministic sensitivity approaches}
  (No. BMI/ONWI-516). Battelle Memorial Inst., Columbus, OH (USA). Office of Nuclear Waste Isolation.

  Worley, B. A. (1987).
  \emph{Deterministic uncertainty analysis}
  (No. CONF-871101-30). Oak Ridge National Lab., TN (USA).

}


\examples{
X <- matrix(runif(800), ncol=8)
borehole(X)
}
