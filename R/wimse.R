## wimse_development
## Functions for integrated weighted variance
eps <- sqrt(.Machine$double.eps)

## erf:
##
## the error function, integrating the normal distribution
erf <- function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}

## optIP.wIMSE:
##
## function that optimizes M inducing points using the weighted IMSE
## centered at w_mean
optIP.wIMSE <- function(Xn, M, theta = NULL, g = 1e-4, w_mean, w_var = NULL,
                        ip_bounds = NULL, integral_bounds = NULL,
                        num_multistart = 15, fix_xm1 = TRUE,
                        epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                        mult = NULL, verbose = TRUE){
  w_mean <- as.vector(w_mean)
  #### Sanity checks ####
  if(length(w_mean)!=ncol(Xn))
    stop('length(w_mean) and ncol(Xn) don\'t match')
  if(M > nrow(Xn))
    warning('Number of inducing points (M) > Neighborhood size (N)')
  if (theta <= 0)
    stop('theta should be a positive number')
  if (g < 0)
    stop('g must be positive')
  if(!is.null(ip_bounds)){
    if(ncol(ip_bounds)!=ncol(Xn))
      stop('ip_bounds and Xn do not have the same number of columns')
    if(nrow(ip_bounds) !=2)
      stop('ip_bounds should be a matrix with two rows')
    if(sum(ip_bounds[1,] < ip_bounds[2,]) < ncol(Xn))
      stop('At least one dimensions bounds in ip_bounds is incorrect')
  }
  if(ncol(integral_bounds)!= ncol(Xn))
    stop('intgral_bounds and Xn do not have the same number of columns')
  if(num_multistart <0 | num_multistart %% 1 !=0)
    stop('num_multistart is not a positive integer')
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')
  if(!is.null(w_var))
    if((length(w_var) != ncol(Xn)) & (length(w_var) != 1))
      stop('length(w_var) doesn\'t match ncol(Xn)')

  ## Fill in missing values with defaults
  if (is.null(integral_bounds)) integral_bounds <- apply(Xn, 2, range)
  if (is.null(ip_bounds)) ip_bounds <- apply(Xn, 2, range)
  if (is.null(theta)) theta <- quantile(dist(Xn), .1)^2
  m0 <- ifelse (fix_xm1, 2, 1)
  ## For timing
  t1 <- proc.time()[3]

  ## Initial integrals and wIMSE
  wimse_track <- vector(length=M)
  n_opt_its <- matrix(nrow=num_multistart, ncol=M+1-m0)

  if(fix_xm1){
    ## Sets first inducing point at w_mean
    Xm <- matrix(w_mean, nrow=1)
    W <- weightW.ligp(Xm, theta, w_mean, w_var, grad=FALSE,
                      integral_bounds=integral_bounds)
    wimse_track[1] <- grad.wIMSE_xm1(xm1=w_mean, Xm, Xn, theta, g,
                                     w_mean, w_var,
                                     integral_bounds=integral_bounds,
                                     epsQ=epsQ, epsK=epsK, mult=mult)$wimse
    KQ <- calcKmQm(Xm, Xn, theta, g, epsQ=epsQ, epsK=epsK, mults=mult)
  } else {Xm <- KQ <- W <- NULL}

  ## Sequentially selects inducing points
  for (m in m0:M){
    xm1.opt <- selectIP.wIMSE(Xm, Xn, theta, g, w_mean, w_var, W, KQ,
                              ip_bounds, integral_bounds, num_multistart)
    n_opt_its[,m-1] <- xm1.opt$its
    wimse_track[m] <- xm1.opt$wIMSE.xm1

    ## Updates W and KQ with newly selected inducing point
    W <- weightW.ligp(Xm, theta, w_mean, w_var, W, xm1.opt$xm1,
                      grad=FALSE, integral_bounds=integral_bounds)
    if(m == 1){
      KQ <- calcKmQm(xm1.opt$xm1, Xn, theta, g, epsQ=epsQ, epsK=epsK, mults=mult)
    } else {
      KQ <- updateKmQm(xm1.opt$xm1, Xm, Xn, theta, g, KQ=KQ)
    }

    Xm <- rbind(Xm, xm1.opt$xm1)
    if (verbose) print(paste("Number of selected inducing points:", m))
  }

  ## For timing
  t2 <- proc.time()[3]

  return(list(Xm=Xm, wimse=exp(wimse_track), time=t2-t1))
}

## selectIP.wIMSE:
##
## function that selects an inducing point by optimizing weighted IMSE
## centerd at w_mean through the use of a multi-start scheme and optim
selectIP.wIMSE <- function(Xm, Xn, theta, g, w_mean, w_var = NULL,
                           W, KQ, ip_bounds, integral_bounds,
                           num_multistart, epsK = sqrt(.Machine$double.eps),
                           epsQ = 1e-5){
  ## Set of multistart points in ip_bounds
  xm.poss <- matrix(runif(num_multistart * ncol(Xn)), ncol=ncol(Xn))
  xm.poss <- t((ip_bounds[2,] - ip_bounds[1,]) * t(xm.poss) + ip_bounds[1,])
  poss.mat <- matrix(nrow=nrow(xm.poss), ncol=ncol(Xn) + 1)
  its <- vector(length=nrow(xm.poss))

  for (j in 1:nrow(xm.poss)){
    out <- opt.wIMSE_xm1(xm.poss[j,,drop=FALSE], Xm, Xn, theta, g, w_mean,
                         w_var, W=W, epsQ=epsQ, epsK=epsK, lower=ip_bounds[1,],
                         upper=ip_bounds[2,], integral_bounds=integral_bounds,
                         KQ=KQ)
    poss.mat[j,] <- c(out$opt$par, out$opt$value)
    its[j] <- out$its
  }
  xm1Index <- which.min(poss.mat[,ncol(poss.mat)])
  xm1.new <- poss.mat[xm1Index,-ncol(poss.mat),drop=FALSE]
  wIMSE.new <- poss.mat[xm1Index,ncol(poss.mat)]

  return(list(xm1=xm1.new, wIMSE.xm1=wIMSE.new, its=its))
}



## wimseByGrid:
##
## function that estimates the weighted IMSE centered at w_mean
## through a grid (Xgrid) of reference locations
wimseByGrid <- function(xm1, Xgrid, Xm, Xn, theta, g = 1e-4, w_mean,
                        w_var = NULL, epsQ = 1e-5,
                        epsK = sqrt(.Machine$double.eps), mult = NULL){
  w_mean <- matrix(w_mean, nrow=1)
  ###### Sanity Checks #####
  if(ncol(xm1)!=ncol(Xm))
    stop('xm1 and Xm do not have the same number of columns')
  if(ncol(Xn)!=ncol(Xm))
    stop('Xn and Xm do not have the same number of columns')
  if(length(w_mean)!=ncol(Xm))
    stop('length(w_mean) doesn\'t match ncol(Xm)')
  if(ncol(Xgrid)!=ncol(Xm))
    stop('Xgrid and Xm do not have the same number of columns')
  if(!is.null(w_var)){
    w_var <- as.vector(w_var)
    if(length(w_var)!=ncol(Xm))
      stop('length(w_var) doesn\'t match ncol(Xm)')
  } else w_var <- theta

  Xm1 <- rbind(Xm, xm1)
  kxgrid <- cov_gen(Xgrid, w_mean, w_var)
  KQ <- calcKmQm(Xm1, Xn, theta, g, epsK=epsK, epsQ=epsQ, mults=mult)
  K_XXm <- cov_gen(Xgrid, Xm1, theta)
  QmiCXX <- solve(KQ$Qm, t(K_XXm))
  Lam.new <- rep(g, nrow(Xgrid))
  KmiCXX <- solve(KQ$Km, t(K_XXm))
  sd2 <- 1 + Lam.new -colSums(t(K_XXm) * (KmiCXX - QmiCXX))

  return(mean(sd2*(kxgrid)))
}


## opt.wIMSE_xm1:
##
## function that optimizes the next inducing point location
## using wIMSE and its gradient
opt.wIMSE_xm1 <- function(xm1, Xm, Xn, theta, g, w_mean, w_var = NULL,
                          W = NULL, lower, upper, maxit = 100,
                          integral_bounds = NULL, epsQ = 1e-5,
                          epsK = sqrt(.Machine$double.eps),
                          LOG = TRUE, KQ = NULL){
  ## objective (and derivative saved)
  ###### Sanity checks ######
  if(is.null(dim(xm1))) xm1 <- matrix(xm1, nrow=1)
  if (!is.null(Xm))
    if(ncol(xm1)!=ncol(Xm))
      stop('xm1 and Xm do not have the same number of columns')
  if(ncol(Xn)!=ncol(xm1))
    stop('Xn and xm1 do not have the same number of columns')
  if(length(w_mean)!=ncol(Xn))
    stop('length(w_mean) and ncol(Xn) do not match')
  if(!is.null(w_var))
    if(length(w_var)!=ncol(Xn))
      stop('length(w_var) and ncol(Xn) do not match')
  if(!is.null(W) & !is.null(Xm))
    if(ncol(W)!=nrow(Xm))
      stop('W supplied doesn\'t match number of entries in Xm')

  deriv <- NULL; its <- 0
  f <- function(xm1, Xm, Xn, theta, nug, w_mean, w_var, W, epsQ, epsK,
                integral_bounds, LOG, KQ){

    out <- grad.wIMSE_xm1(xm1, Xm, Xn, theta, g=nug, w_mean, w_var,
                          W, epsQ=epsQ, epsK=epsK,
                          integral_bounds=integral_bounds, KQ=KQ)

    its <<- its + 1
    if (LOG){
      deriv <<-list(x=log(out$wimse), df=out$dwimse/out$wimse)
      return(log(out$wimse))
    } else {
      deriv <<-list(x=out$wimse, df=out$dwimse)
      return(out$wimse)
    }
  }

  ## derivative read from global variable
  df <- function(xm1, Xm, Xn, theta, nug, w_mean, w_var, W, epsQ, epsK,
                 integral_bounds, LOG, KQ) {
    if(any(xm1 != deriv$wimse))
      stop("xs don't match for successive f and df calls")
    return(deriv$df)
  }
  control <- list(maxit=maxit, pgtol= 1e-6) ## set up controls

  opt <- optim(xm1, f, df, lower=lower, upper=upper, method="L-BFGS-B",
               control=control, Xm=Xm, Xn=Xn, theta=theta, nug=g,
               w_mean=w_mean, w_var=w_var, W=W, epsQ=epsQ, epsK=epsK,
               integral_bounds=integral_bounds, LOG=LOG, KQ=KQ)

  return(list(opt=opt, its=its))
}

## weightW.ligp:
##
## function that calculates (or updates) the W matrix used
## for closed form wIMSE calculation
weightW.ligp <- function(Xm, theta, w_mean, w_var = NULL, W = NULL,
                         xm1 = NULL, grad = TRUE, integral_bounds){
  w_mean <- as.vector(w_mean)

  if(!is.null(xm1)) {
    xm1 <- matrix(xm1, nrow=1)
    ndim <- ncol(xm1)
  } else ndim <- ncol(Xm)

  ###### Sanity Checks #####
  if(length(w_mean)!=ndim)
    stop('length(w_mean) doesn\'t match the number of inducing point dimensions')
  if(ncol(integral_bounds)!=ndim)
    stop('ncol(integral_bounds) doesn\'t match the number of ',
         'inducing point dimensions')
  if (!is.null(xm1) & !is.null(Xm))
    if(ncol(xm1)!=ncol(Xm))
      stop('xm1 and Xm do not have the same number of columns')

  ## Deals with user input variance for weight, including zero variance
  if(!is.null(w_var)) {
    w_var <- as.vector(w_var)
    if(length(w_var)!= ndim)
      stop('length(w_var) doesn\'t match the number of inducing point dimensions')
    nonzero_wvar <- w_var !=0
    if (sum(nonzero_wvar) < ndim){
      w_var <- w_var[nonzero_wvar]
      w_mean <- w_mean[nonzero_wvar]
      Xm0 <- Xm[, !nonzero_wvar, drop=FALSE]
      Xm <- Xm[, nonzero_wvar, drop=FALSE]
      zero_var_int_bounds <- integral_bounds[, !nonzero_wvar, drop=FALSE]
      volume_constant <- prod(1/(zero_var_int_bounds[2,] -
                                   zero_var_int_bounds[1,]))
      integral_bounds <- integral_bounds[, nonzero_wvar, drop=FALSE]
      if(!is.null(xm1)) {
        xm10 <- xm1[, !nonzero_wvar, drop=FALSE]
        xm1 <- xm1[, nonzero_wvar, drop=FALSE]
      } else xm10 <- NULL
    }
  } else nonzero_wvar <- ndim
  lowBounds <- integral_bounds[1,]; upBounds <- integral_bounds[2,]
    ## updates W matrix if old matrix is supplied
    if(!is.null(W) & !is.null(xm1)){
      wwij <- weightWij.ligp(xm1, Xm, theta, w_mean, w_var=w_var,
                             integral_bounds=integral_bounds)
      wwii <- weightWij.ligp(xm1, theta=theta, w_mean=w_mean,
                             w_var=w_var, integral_bounds=integral_bounds)
      if(sum(nonzero_wvar) < ndim){
        wwij$w <- wwij$w * Wij.ligp(xm10, Xm0, theta, zero_var_int_bounds)$w
        wwii$w <- wwii$w * Wij.ligp(xm10, NULL, theta, zero_var_int_bounds)$w
      }

      W1 <- cbind(W, wwij$w); W2 <- c(wwij$w, wwii$w)
      W <- rbind(W1, W2)

      ## dW: each column corresponds to a gradient of last column of W
      ##  wrt a different dimension
      if(grad) dW <- rbind(wwij$dw, wwii$dw)

      ## Performs all calculations of W and gradients
    } else {
      if(!is.null(xm1)) Xm <- rbind(Xm, matrix(xm1, nrow=1))
      if(grad) dW <- matrix(1, ncol=ncol(Xm), nrow=nrow(Xm))
      W <- matrix(1, ncol=nrow(Xm), nrow=nrow(Xm))

      ## Computations if given w_var
      if(!is.null(w_var)){
        th_plus_2wvar <- 2*w_var + theta
        theta_times_wvar <- w_var*theta
        denom_frac <- sqrt((2*w_var + theta)/(w_var*theta))
        erf_denom <- theta_times_wvar*sqrt((2*w_var + theta)/(w_var*theta))
        num_lower <- lowBounds * th_plus_2wvar
        num_upper <- upBounds * th_plus_2wvar
        wmean2_wvar <- w_mean^2/w_var
        exp_denom <- theta_times_wvar*th_plus_2wvar
      }

      for(i in 1:nrow(W)){
        if(is.null(w_var)){
          sum3 <- sweep(Xm, 2, w_mean + Xm[i,], '+')
          erfs <- erf((sum3 - 3*lowBounds)/sqrt(3*theta)) -
            erf((sum3 - 3*upBounds)/sqrt(3*theta))
          if(grad){
            derf <-  2/sqrt(3*pi*theta) * (
              exp(-(sum3[ncol(W),]-3*lowBounds)^2/(3*theta)) -
                exp(-(sum3[ncol(W),]-3*upBounds)^2/(3*theta)))
          }

        } else {
          weighted_sum3_a <- sweep(Xm, 2, w_var,'*')
          weighted_sum3 <- sweep(weighted_sum3_a, 2, theta*w_mean +
                                   w_var*Xm[i,], '+')
          erfs <- erf(sweep(sweep(weighted_sum3, 2, num_lower),
                            2, erf_denom,"/")) -
            erf(sweep(sweep(weighted_sum3, 2, num_upper), 2, erf_denom,"/"))
          if(grad){
            derf <- 2/(sqrt(pi)*theta*denom_frac)*(
              exp(-(weighted_sum3[ncol(W),] - num_lower)^2/exp_denom) -
                exp(-(weighted_sum3[ncol(W),] - num_upper)^2/exp_denom))
          }

        }
        dw <- w <- matrix(nrow=ncol(W), ncol=ncol(Xm))

        for(j in 1:ncol(W)){
          if(is.null(w_var)){
            ## Calculates w(x_ik,x_jk) for each dimension
            w[j,] <- sqrt(pi *theta/12) *
              exp((2/(3*theta))*(w_mean*(Xm[i,] + Xm[j,] - w_mean) - Xm[i,]^2
                                 - Xm[j,]^2 + Xm[i,]*Xm[j,]))
            ## Calculates dw(i,j)/dxj,k (last entry in Xm)
            if(grad){
              if (j==nrow(W)){
                dw[j,] <- w[j,] * (
                  2/(3*theta) * (w_mean - 2*Xm[j,] + Xm[i,]) * erfs[j,] + derf)
                if(j==i) dw <- 2 * dw
              }
            }

          } else {
            ## Calculates w(x_ik,x_jk) for each dimension
            w[j,] <- sqrt(pi)/(2*denom_frac) *
              exp((weighted_sum3[j,]^2)/exp_denom -
                    wmean2_wvar - (Xm[i,]^2 + Xm[j,]^2)/theta)
            if(grad){
              ## Calculates for dw(i,j)/dxj,k (last entry in Xm)
              if (j == nrow(W)){
                dw[j,] <- w[j,] * (
                  -2/(theta*th_plus_2wvar) * (
                    (w_var + theta)*Xm[j,]- w_var*Xm[i,] -
                      theta*w_mean) * erfs[j,] + derf)
                if(j == i) dw <- 2 * dw
              }
            }
          }

        }
        w <- w * erfs

        ## stores product of w(x_ik,x_jk) across k
        W[i,] <- apply(w, 1, prod)

        ## Gradients for pt x_i for each dimension wrt x_j
        if(grad){
          for (k in 1:ncol(Xm))
            dW[i, k] <- prod(w[ncol(W),-k]) * dw[ncol(W),k]
        }
        ##############
      }
      if (sum(nonzero_wvar) < ndim){
        W0 <- W.ligp(Xm0, theta, xm1=xm10, integral_bounds=zero_var_int_bounds)
        W <- W * W0
      }
    }

    if(grad){
      dW_full <- dW
      if(exists('nonzero_wvar'))
        if (sum(nonzero_wvar) < ndim){
          dW_full <- matrix(0, ncol=ndim, nrow=nrow(W))
          dW_full[,nonzero_wvar] <- dW
        }
    }

    ## dW: gradient of W wrt last row of Xm (or xm1 if supplied)
    ## each column is the last(only nonzero) row/column of dW_d matrix
    ifelse(grad, return(list(W=W, dW=dW_full)), return(W))
  }

## weightWij.ligp:
##
## function that calculates individual components of the W matrix used
## for closed form wIMSE calculation
weightWij.ligp <- function(x1, X2 = NULL, theta, w_mean, w_var = NULL,
                           integral_bounds){
  lowBounds <- integral_bounds[1,]
  upBounds <- integral_bounds[2,]

  x1 <- matrix(x1, nrow=1)
  w_mean <- matrix(w_mean, nrow=1)
  doubleint <- is.null(X2)
  if(is.null(X2)){X2 <- x1} else {X2 <- as.matrix(X2)}

  ###### Sanity Checks #####
  if(ncol(x1)!=ncol(X2))
    stop('x1 and X2 do not have the same number of columns')
  if(length(w_mean)!=ncol(x1))
    stop('length(w_mean) doesn\'t match ncol(x1)')
  if(!is.null(w_var))
    if(length(w_var)!=ncol(x1))
      stop('length(w_var) doesn\'t match ncol(x1)')

  Dw <- dw <- w <- matrix(1, nrow=nrow(X2), ncol=ncol(X2))

  ## Calculations given w_var
  if(!is.null(w_var)){
    th_plus_2wvar <- 2*w_var + theta
    theta_times_wvar <- w_var*theta
    denom_frac <- sqrt((2*w_var + theta)/(w_var*theta))
    erf_denom <- theta_times_wvar*sqrt((2*w_var + theta)/(w_var*theta))
    num_lower <- lowBounds * th_plus_2wvar
    num_upper <- upBounds * th_plus_2wvar
    wmean2_wvar <- w_mean^2/w_var
    exp_denom <- theta_times_wvar*th_plus_2wvar
  }
  ## Separate wwij for each dimension
  for(j in 1:nrow(X2)){

    if(is.null(w_var)){
      ## Store in matrix w and them multiply
      w[j,] <- sqrt(pi *theta/12) *
        exp((2/(3*theta))*(w_mean*(x1[1,] + X2[j,] - w_mean) - x1[1,]^2
                           - X2[j,]^2 + x1[1,]*X2[j,]))
      sum3 <- w_mean + x1[1,] + X2[j,]
      erfs <- (erf((sum3 - 3*lowBounds)/sqrt(3*theta)) -
                 erf((sum3 - 3*upBounds)/sqrt(3*theta)))
      ## Calculation for dw(i,j)/dxi,k
      derf <-  2/sqrt(3*pi*theta) * (exp(-(sum3-3*lowBounds)^2/(3*theta)) -
                                       exp(-(sum3-3*upBounds)^2/(3*theta)))
      dw[j,] <- w[j,] * (2/(3*theta) *
                           (w_mean - 2*x1[1,] + X2[j,]) * erfs + derf)
    } else {
      weighted_sum3 <- theta*w_mean + w_var*(x1[1,] + X2[j,])
      erfs <- erf((weighted_sum3 - num_lower)/erf_denom) -
        erf((weighted_sum3 - num_upper)/erf_denom)
      w[j,] <- sqrt(pi)/(2*denom_frac) *
        exp((weighted_sum3^2)/exp_denom - wmean2_wvar -
              (x1[1,]^2 + X2[j,]^2)/theta)
      ## Calculataion for dw(i,j)/dxi,k
      derf <- 2/(sqrt(pi)*theta*denom_frac)*
        (exp(-(weighted_sum3 - num_lower)^2/exp_denom) -
           exp(-(weighted_sum3 - num_upper)^2/exp_denom))
      dw[j,] <- w[j,] * (-2/(theta*th_plus_2wvar) *
                           ((w_var + theta)*x1[1,] - w_var*X2[j,] -
                              theta*w_mean) * erfs + derf)
    }

    if(doubleint) dw[j,] <- 2 * dw[j,]
    w[j,] <- w[j,] * erfs
  }
  wj <- apply(w, 1, prod)

  dim_vec <- 1:ncol(X2)
  for(k in dim_vec){
    for (i in dim_vec[-k]){
      Dw[,k] <- Dw[,k] * w[,i]
    }
  }

  Dw <- Dw * dw
  return(list(w=wj, dw=Dw))
}


## calc_wIMSE:
##
## function that calculates the weighted IMSE given a new inducing
## point location
calc_wIMSE <- function(xm1, Xm = NULL, Xn, theta = NULL, g = 1e-4,
                       w_mean, w_var = NULL, integral_bounds = NULL,
                       epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                       mult = NULL){
  if (is.null(dim(xm1))) xm1 <- matrix(xm1, nrow=1)
  if(is.null(integral_bounds)) integral_bounds <- apply(Xn, 2, range)
  if(is.null(theta)) theta <- quantile(dist(Xn), .1)^2
  w_mean  <- as.vector(w_mean)

  ###### Sanity Checks #####
  if(!is.null(Xm))
    if(ncol(xm1)!=ncol(Xm))
      stop('xm1 and Xm do not have the same number of columns')
  if(ncol(Xn)!=ncol(xm1))
    stop('Xn and xm1 do not have the same number of columns')
  if(length(w_mean)!=ncol(xm1))
    stop('length(w_mean) doesn\'t match ncol(xm1)')
  if(!is.null(w_var)){
    w_var  <- as.vector(w_var)
    if(length(w_var)!=ncol(xm1))
      stop('length(w_var) doesn\'t match ncol(xm1)')
  }

  Xm1 <- rbind(Xm, xm1)
  KQ <- calcKmQm(Xm1, Xn, theta, g, epsK=epsK, epsQ=epsQ, mults=mult)
  W <- weightW.ligp(Xm1, theta, w_mean, w_var,
                    grad=FALSE, integral_bounds=integral_bounds)
  if(is.null(w_var)) w_var <- rep(theta, ncol(xm1))
  wimse <- prod((1 + g) * sqrt(w_var*pi)/2)

  for (k in 1:ncol(Xm1)){
    wimse <- wimse * (erf((w_mean[k] - integral_bounds[1,k])/sqrt(w_var[k])) -
                        erf((w_mean[k] - integral_bounds[2,k])/sqrt(w_var[k])))
  }
  wimse <- wimse - sum(diag(KQ$Kmi %*% W - solve(KQ$Qm, W)))

  return(wimse)
}


## W.ligp:
##
## function that calculates (or updates) the W matrix used
## for closed form IMSE calculation
W.ligp <- function(Xm, theta, W = NULL, xm1 = NULL, integral_bounds = NULL){
  if(is.null(integral_bounds)){
    lowBounds <- rep(0,ncol(Xm)); upBounds <- rep(1,ncol(Xm))
  } else lowBounds <- integral_bounds[1,]; upBounds <- integral_bounds[2,]

  ## updates W matrix if old matrix is supplied
  if (!is.null(W) & !is.null(xm1)){
    wij <- Wij.ligp(xm1, Xm, theta, integral_bounds=integral_bounds)
    wii <- Wij.ligp(xm1, theta=theta, integral_bounds=integral_bounds)
    W1 <- cbind(W, wij$w); W2 <- c(wij$w, wii$w)
    W <- rbind(W1, W2)
  } else {
    if(!is.null(xm1)) Xm <- rbind(Xm, matrix(xm1, nrow=1))
    W <- matrix(1, ncol=nrow(Xm), nrow=nrow(Xm))
    for(i in 1:nrow(W)){
      for(j in 1:ncol(W)){
        ## Calculates w(x_ik,x_jk) for each dimension
        w <- sqrt(pi *theta/8) *
          exp((1/(2*theta))*(- Xm[i,]^2 - Xm[j,]^2 + 2*Xm[i,]*Xm[j,]))
        erfs <- erf((Xm[j,] + Xm[i,] - 2*lowBounds)/sqrt(2*theta)) -
          erf((Xm[j,] + Xm[i,] - 2*upBounds)/sqrt(2*theta))

        w <- w * erfs
        ## stores product of w(x_ik,x_jk) across k
        W[i,j] <- prod(w)
      }
    }
  }

  return(W)
}


## Wij.ligp:
##
## function that calculates individual components of the W matrix used
## for closed form IMSE calculation
Wij.ligp <- function(x1, X2 = NULL, theta, integral_bounds){
  lowBounds <- integral_bounds[1,]
  upBounds <- integral_bounds[2,]

  x1 <- matrix(x1, nrow=1);
  doubleint <- is.null(X2)
  if(is.null(X2)){X2 <- x1} else {X2 <- as.matrix(X2)}

  w <- matrix(1, nrow=nrow(X2), ncol=ncol(X2))
  ## Separate wij for each dimension
  for(j in 1:nrow(X2)){
    ## Store in matrix and then multiply
    w[j,] <- sqrt(pi *theta/8) *
      exp((1/(2*theta))*(- x1[1,]^2 - X2[j,]^2 + 2*x1[1,]*X2[j,]))
    erfs <- (erf((x1[1,] + X2[j,] - 2*lowBounds)/sqrt(2*theta)) -
               erf((x1[1,] + X2[j,] - 2*upBounds)/sqrt(2*theta)))
    w[j,] <- w[j,] * erfs
  }
  wj <- apply(w, 1, prod)

  return(list(w=wj))
}


### calc_IMSE:
##
## function that calculates the IMSE given a new inducing
## point location
calc_IMSE <- function(xm1, Xm = NULL, X, theta = NULL, g = 1e-4,
                      integral_bounds = NULL, epsK = sqrt(.Machine$double.eps),
                      epsQ = 1e-5, mult = NULL){
  if (is.null(dim(xm1))) xm1 <- matrix(xm1, nrow=1)
  if(is.null(integral_bounds))
    integral_bounds <- apply(X, 2, range)
    if (is.null(Xm)) {Xm1 <- matrix(xm1, nrow=1)
    } else Xm1 <- rbind(Xm, xm1)
    if(is.null(theta)) theta <- quantile(dist(X), .1)^2

    KQ <- calcKmQm(Xm=Xm1, Xn=X, theta=theta, g=g, epsK=epsK, epsQ=epsQ,
                   mults=mult)
    W <- W.ligp(Xm, theta, xm1=xm1, integral_bounds=integral_bounds)

    imse <- (1 + g - sum(diag(KQ$Kmi %*% W - solve(KQ$Qm, W))))
    return(imse)
}



## grad.wIMSE_xm1:
##
## function that calculates the gradient of the weighted IMSE
## with respect to an additional inducing point xm1
grad.wIMSE_xm1 <- function(xm1, Xm = NULL, Xn, theta, g, w_mean,
                           w_var = NULL, W = NULL, integral_bounds = NULL,
                           KQ = NULL, epsK = sqrt(.Machine$double.eps),
                           epsQ = 1e-5, grad = TRUE, mult = NULL){

  w_mean  <- as.vector(w_mean)
  xm1 <- matrix(xm1, nrow=1)

  if(is.null(integral_bounds))
    integral_bounds <- rbind(rep(0,ncol(Xn)), rep(1,ncol(Xn)))

  lowBounds <- integral_bounds[1,]; upBounds <- integral_bounds[2,]

  ###### Sanity Checks #####
  if(!is.null(Xm))
    if(ncol(xm1)!=ncol(Xm))
      stop('xm1 and Xm do not have the same number of columns')
  if(ncol(Xn)!=ncol(xm1))
    stop('Xn and xm1 do not have the same number of columns')
  if(length(w_mean)!=ncol(xm1))
    stop('length(w_mean) doesn\'t match ncol(xm1)')
  if(!is.null(w_var))
    if(length(w_var)!=ncol(xm1))
      stop('length(w_var) doesn\'t match ncol(xm1)')
  if(!is.null(W) & !is.null(Xm))
    if(ncol(W)!=nrow(Xm))
      stop('W supplied doesn\'t match number of entries in Xm')

  Xm1 <- rbind(Xm, xm1)
  numDimen <- ncol(Xm1)

  if(is.null(KQ)){
    KQ <- calcKmQm(Xm1, Xn, theta, g, epsK=epsK, epsQ=epsQ,
                   mults=mult)
  } else KQ <- updateKmQm(xm1, Xm, Xn, theta, g, KQ=KQ)

  Qmi <- solve(KQ$Qm)
  KmiQmi <- KQ$Kmi - Qmi
  wWij <- weightW.ligp(Xm, theta, w_mean, w_var, W, xm1,
                       integral_bounds=integral_bounds, grad=grad)

  ifelse(grad, W <- wWij$W, W <- wWij)

  ## Calculate integrated weighted variance
  if(is.null(w_var)) w_var <- rep(theta, numDimen)

  nonzero_wvar <- w_var !=0
  if (sum(nonzero_wvar) > 0){
    wimse <- prod((1 + g) * sqrt(w_var[nonzero_wvar]*pi) /2)
  } else wimse <- 1


  for (k in 1:numDimen){
    if (w_var[k] != 0){
      wimse <- wimse * (erf((w_mean[k] - lowBounds[k])/sqrt(w_var[k])) -
                          erf((w_mean[k] - upBounds[k])/sqrt(w_var[k])))
    }
  }
  # return(c(wimse, sum(diag(KQ$Kmi %*% W - solve(KQ$Qm, W)))))
  wimse <- wimse - sum(diag(KQ$Kmi %*% W - solve(KQ$Qm, W)))
  if (!grad) return(wimse)

  ## Calculate gradients of wimse
  dwimse <- vector(length=numDimen)

  for (k in 1:numDimen){
    dKmi <- gradKmi_xmi(nrow(Xm1), k, Xm1, theta, KQ$Kmi)
    dQm <- gradQm_xmi(nrow(Xm1), k, Xm1, Xn, theta, KQ, dKmi)
    dwimse[k] <-  - sum((dKmi$dKmi + Qmi %*% dQm %*% Qmi) * W)
    if (nonzero_wvar[k]){
      dW <- matrix(0, nrow=nrow(Xm1), ncol=nrow(Xm1))
      dW[,nrow(Xm1)] <- wWij$dW[,k]
      dW[nrow(Xm1),] <- wWij$dW[,k]
      dwimse[k] <-  dwimse[k]  - sum(KmiQmi * dW)
    }
  }

  return(list(wimse=wimse, dwimse=dwimse))
}
