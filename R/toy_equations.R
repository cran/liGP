###################################################################
#-----Using inducing points in local approximate GP framework-----#
###################################################################

## Herbie's tooth function
## X: a matrix of design locations
herbtooth <- function(X)
{
  g <- function(z)
    return(exp(-(z - 1)^2) + exp(-0.8*(z + 1)^2) - 0.05*sin(8*(z + 0.1)))
  return(-apply(apply(X, 2, g), 1, prod))
}



## Used to generate data in [0,1]^8
borehole <- function(X)
{
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # X = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################
  ## Assume inputs are on the unit cube
  ## rw in [.05,.15]
  ## r in [100, 50000]
  ## Tu in [63070,115600]
  ## Hu in [990, 1110]
  ## Tl in [63.1, 116]
  ## Hl in [700, 820]
  ## L in [1120, 1680]
  ## Kw in [9855, 12045]
  rw <- X[,1] * (.15 - .05) + .05
  r  <- X[,2] * (50000 - 100) + 100
  Tu <- X[,3] * (115600 - 63070) + 63070
  Hu <- X[,4] * (1110 - 990) + 990
  Tl <- X[,5] * (116 - 63.1) + 63.1
  Hl <- X[,6] * (820 - 700) + 700
  L  <- X[,7] * (1680 - 1120) + 1120
  Kw <- X[,8] * (12045 - 9855) + 9855

  frac1 <- 2 * pi * Tu * (Hu-Hl)

  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)

  y <- frac1 / frac2

  return(y)

}
