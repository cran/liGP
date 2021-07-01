###################################################################
#-----Using inducing points in local approximate GP framework-----#
###################################################################
eps <- sqrt(.Machine$double.eps)


## calcALC.giGP:
##
## function that calculates ALC based on a set of proposed inducing
## points xm1. ALC is calculated based on a series of reference
## locations, Xref.
calcALC.giGP <- function (xm1, Xm = NULL, Xref, Xn, Yn, theta, g,
                          KQ = NULL, K_xrefm = NULL, nu = NULL,
                          rep_list = NULL){
  if(is.null(dim(xm1))) xm1 <- matrix(xm1, nrow=1)
  ## predicted variance for Xref (design reference set)
  if (is.null(K_xrefm)) K_xrefm <- hetGP::cov_gen(Xref, Xm, theta=theta)

  ## Loops & calculates ALC for multiple points in Xref
  ALC <- vector(length=nrow(xm1))
  for (i in 1:nrow(xm1)){
    Xm1 <-rbind(Xm, xm1[i,])
    xm1i <- xm1[i,,drop=FALSE]
    ## Cross-covariance b/w design reference set & amended inducing pts
    K_xrefm1 <- hetGP::cov_gen(X1=Xref, X2=xm1i, theta=theta)

    if (is.null(Xm)) { # Checks if there are no previous inducing pts
      KQ <- calcKmQm(Xm=Xm1, Xn=Xn, theta, g, mults=rep_list$mult)
      K_xref <- K_xrefm1
    } else {
      if (!is.null(KQ)) {KQ <- updateKmQm(xm1i, Xm, Xn, theta, g, KQ)
      } else KQ <- calcKmQm(Xm=Xm1, Xn, theta, g, mults=rep_list$mult)
      K_xref <- cbind(K_xrefm, K_xrefm1) # Amends K_xrefm to K_Xref,m+1
    }

    if (is.null(nu))
      nu <- calcNu(Yn, KQ$Qm, KQ$OL, KQ$OLiKnm, rep_list=rep_list)

    QmiKxref <- try(solve(KQ$Qm,t(K_xref)), silent=TRUE)
    if(class(QmiKxref)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    KmiKxref <- try(solve(KQ$Km,t(K_xref)), silent=TRUE)
    if(class(KmiKxref)[1] == 'try-error')
      stop('K matrix is numerically instable. Consider increasing epsK.')

    alc <- nu * (1 + g - sum(t(K_xref) * (KmiKxref - QmiKxref))/nrow(Xref))
    ALC[i] <- alc
  }

  return(list(ALC=ALC, nu=nu))
}


## calcALC:
##
## function that calculates the standard definition of ALC based
## on a series of reference locations, Xref. This does not include
## inducing points.
calcALC <- function(Xref, X, Y, theta, g = 1e-4, mult = NULL){
  N <- ifelse(is.null(mult), nrow(X), sum(mult))

  Kx <- hetGP::cov_gen(X, theta=theta)
  if(is.null(mult)){
    Sigma <- ifelse (N == 1, Kx + g, Kx + diag(rep(g, N)))
  } else{
    Sigma <- Kx + diag(rep(g, nrow(X))/mult)
  }

  SigmaI <- solve(Sigma)
  K.xref.x <- hetGP::cov_gen(Xref, X, theta=theta)

  nu_hat <- t(Y) %*% SigmaI %*% Y / N
  alc.stat <- nu_hat * (1 + g - K.xref.x %*% SigmaI %*% t(K.xref.x))


  return(alc.stat)
}

## gradALC_xm1:
##
## function that calculates the gradient of ALC with respect to an
## additional inducing point xm1. ALC is based on a series of reference
## locations, Xref.
gradALC_xm1 <- function(xm1, Xm, Xref, Xn, Yn, theta, g = 1e-4, nu = NULL,
                        KQ = NULL, K_xrefm = NULL, rep_list = NULL){

  if(is.null(dim(xm1))) {
    if (length(xm1) != ncol(Xm)) {
      xm1 <- matrix(xm1, ncol=ncol(Xm))
    } else xm1 <- matrix(xm1, nrow=1)
  }
  Xm1 <-rbind(Xm, xm1)
  m1 <- nrow(Xm1); N <- nrow(Xn)

  if (is.null(K_xrefm)) K_xrefm <- hetGP::cov_gen(X1=Xref, X2=Xm,
                                                  theta=theta)

  ## Calculates ALC for Xref
  ## Cross-covariance b/w design reference set & amended inducing pts
  K_xref.m1 <- hetGP::cov_gen(X1=Xref, X2=xm1, theta=theta)

  if (is.null(Xm)) { # Checks if there are no previous inducing pts
    KQ <- calcKmQm(Xm1, Xn=Xn, theta=theta, g=g, mults=rep_list$mult)
    K_xref <- K_xref.m1
  } else {
    if(!is.null(KQ)) {KQ <- updateKmQm(xm1, Xm, Xn, theta=theta,
                                       g=g, KQ=KQ)
    } else  KQ <- calcKmQm(Xm1, Xn=Xn, theta=theta, g=g, mults=rep_list$mult)
    K_xref <- cbind(K_xrefm, K_xref.m1) # Amends K_xrefm to K_Xref,m+1
  }
  Qm1iKxref <- try(solve(KQ$Qm, t(K_xref)), silent=TRUE)
  if (class(Qm1iKxref)[1] == 'try-error')
    stop('Q matrix is numerically instable. Consider increasing epsQ.')
  Km1iKxref <- try(solve(KQ$Km, t(K_xref)), silent=TRUE)
  if (class(Km1iKxref)[1] == 'try-error')
    stop('K matrix is numerically instable. Consider increasing epsK.')

  if (is.null(nu))
    nu <- calcNu(Yn, KQ$Qm, KQ$OL, KQ$OLiKnm, rep_list=rep_list)

  alc <- nu*(1 + g - sum(t(K_xref) * (Km1iKxref- Qm1iKxref))/nrow(Xref)) #


  ## Calculates gradient wrt xm1 for each dimension
  ## All values have M+1 inducing points (being differentiated wrt last xm1)
  d.alc <- matrix(data=0, nrow=nrow(xm1), ncol=ncol(xm1))
  K_mi <- solve(KQ$Km) #jitter already in KQ$Km
  for (i in (nrow(Xm) + 1):nrow(Xm1)){
    for (di in 1:ncol(xm1)) {
      dKnm <- gradKmn_xmi(i, di, Xm=Xm1, Xn=Xn, theta=theta) #vector (length=n)
      dK_xref.m1 <- gradKmn_xmi(i, di, Xm=Xm1, Xn=Xref, theta=theta) #vector (length=nrow(Xref))
      part1 <- dK_xref.m1 * matrix(Km1iKxref[i,] - Qm1iKxref[i,], nrow=1) #
      dK <- gradKmi_xmi(i, di, Xm1, theta, Kmi=K_mi)
      dKm1i <- dK$dKmi; dKm1 <- dK$dKm #(M+1 by M+1 matrix; M+1 length vector)
      dQm1 <- gradQm_xmi(i, di, Xm=Xm1, Xn=Xn, theta=theta, KQ, dK)

      part2 <- sum(Qm1iKxref * (dQm1 %*% Qm1iKxref)) +
        sum(t(K_xref) * (dK$dKmi %*% t(K_xref)))
      d.alc[i-nrow(Xm),di] <- drop(2*sum(part1) + part2)
    }
  }

  d.alc <- nu*d.alc

  return(list(alc=drop(alc), dalcs=-d.alc))
}

## obj.ALC.giGP:
##
## objective function used to optimize ALC for an additional inducing point
obj.ALC.giGP <- function(xm1, Xm = NULL, Xref, Xn, Yn, theta, g,  KQ = NULL, K_xrefm = NULL,
                         rep_list = NULL){
  return(log(drop(calcALC.giGP(xm1, Xm, Xref, Xn, Yn, theta, g,  KQ, K_xrefm,
                               rep_list=rep_list))))
}

## optALC_xm1:
##
## function that optimizes the location of the m+1st inducing point
## based on ALC (calculated using a series of reference locations Xref)
optALC_xm1 <- function(xm1, Xm = NULL, Xref, Xn, Yn, theta, g, nu = NULL,
                       ip_bounds, maxit = 100, KQ = NULL, K_xrefm = NULL,
                       tol = .01, rep_list = NULL){
  ## objective (and derivative saved)

  deriv <- NULL
  f <- function(xm1, Xm, Xref, Xn, Yn, theta, nug, nu, KQ,
                K_xrefm, rep_list){
    out <- gradALC_xm1(xm1, Xm, Xref, Xn, Yn, theta, g=nug, nu,
                       KQ, K_xrefm, rep_list=rep_list)

    deriv <<-list(x=xm1, df=out$dalcs/out$alc)
    return(log(out$alc))

  }

  ## derivative read from global variable
  df <- function(xm1, Xm, Xref, Xn, Yn, theta, nug, nu, KQ, K_xrefm, rep_list) {
    if(any(xm1 != deriv$x)) stop("xs don't match for successive f and df calls")
    return(deriv$df)
  }


  control <- list(maxit=maxit, pgtol=tol)
  xm1 <- as.vector(xm1)
  opt <- optim(xm1, f, df, lower=ip_bounds[1,], upper=ip_bounds[2,],
               method="L-BFGS-B", control=control, Xm=Xm, Xref=Xref, Xn=Xn,
               Yn=Yn, theta=theta, nug=g, nu=nu, KQ=KQ, K_xrefm=K_xrefm,
               rep_list=rep_list)

  return(opt)
}

## optIP.ALC:
##
## function that optimizes the locations for M inducing points centered
## around Xc and Xn using ALC (calculated using a reference set Xref)
optIP.ALC <- function(Xc, Xref = NULL, M, Xn, Yn, theta = NULL, g = 1e-4,
                      ip_bounds = NULL, num_thread = 1, num_multistart = 15,
                      epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01,
                      rep_list = NULL, verbose = TRUE){
  t1 <- proc.time()[3]
  ## Sanity checks
  if (is.null(dim(Xc))) Xc <- matrix(Xc, nrow=1)
  if(M > nrow(Xn))
    warning('Number of inducing points (M) > Neighborhood size (N)')
  if (is.null(ip_bounds)) ip_bounds <- apply(Xn, 2, range)
  if(is.null(theta)) theta <- quantile(dist(Xn), .1)^2

  dim <- ncol(Xn)
  Xref <- rbind(Xc, Xref)
  Xm <- Xc

  alc_track <- vector(length=M)
  alc.init <- calcALC.giGP(xm1=Xc, Xm=NULL, Xref=Xref, Xn, Yn,
                                  theta, g, nu=NULL, rep_list=rep_list)
  alc_track[1] <- log(alc.init$ALC)

  KQ <- calcKmQm(Xm, Xn, theta, g, epsK=epsK, epsQ=epsQ, mults=rep_list$mult)
  K_xrefm <- hetGP::cov_gen(X1=Xref, X2=Xm, theta=theta)

  for(i in 2:M){
    ## Optimization of ALC surface to find next inducing point location
    multiStart_bnds <- boundBox(Xc, Xm, ip_bounds[1,], ip_bounds[2,], perc=.1)

    ## Points within certain % of window from Xc
    xm1s <- matrix(runif(num_multistart*dim), ncol=dim)
    xm1.start <- t((multiStart_bnds$windU - multiStart_bnds$windL) * t(xm1s) +
                     multiStart_bnds$windL)
    xm1s <- matrix(runif(num_multistart*dim), ncol=dim)
    ## Points within box drawn by existing Xm
    xm2.start <- t((multiStart_bnds$upper - multiStart_bnds$lower) * t(xm1s) +
                     multiStart_bnds$lower)
    xm1 <- rbind(xm1.start, xm2.start)

    if (num_thread > 1){
      cl <- makeCluster(num_thread)
      registerDoParallel(cl)
      poss.x <- foreach(j=1:(nrow(xm1)), .errorhandling="pass",
                        .combine='rbind', .multicombine=TRUE,
                        .packages=c('laGP', 'hetGP')) %dopar% {
                          out <- try(optALC_xm1(xm1[j,,drop=FALSE]+1e-7, Xm,
                                                Xref, Xn, Yn, theta, g,
                                                ip_bounds=ip_bounds,
                                                maxit=100, KQ=KQ,
                                                K_xrefm=K_xrefm, tol=tol,
                                                rep_list=rep_list),
                                     silent=TRUE)
                          c(out$par, out$value)
                        }
      stopCluster(cl)
    } else {
      poss.x <- matrix(nrow=nrow(xm1), ncol=dim + 1)
      for (j in 1:nrow(xm1)){
        out <- try(optALC_xm1(xm1[j,,drop=FALSE]+1e-7, Xm, Xref, Xn, Yn,
                              theta, g, ip_bounds=ip_bounds, maxit=100, KQ=KQ,
                              K_xrefm=K_xrefm, tol=tol, rep_list=rep_list),
                   silent=TRUE)
        poss.x[j,] <- c(out$par, out$value)
      }
    }

    w.min <- which.min(poss.x[,ncol(poss.x)])
    xm_new <- poss.x[w.min,-ncol(poss.x),drop=FALSE]

    KQ <- updateKmQm(xm_new, Xm, Xn, theta, g, KQ=KQ)
    k_xrefm1 <- hetGP::cov_gen(Xref, xm_new, theta=theta)
    K_xrefm <- cbind(K_xrefm, k_xrefm1)
    Xm <- rbind(Xm, xm_new)
    if (verbose == TRUE) print(paste("Number of selected inducing points:",i))

    alc_track[i] <- poss.x[w.min,ncol(poss.x)]
  }

  t2 <- proc.time()[3]
  return(list(Xm=Xm, alc=exp(alc_track), time=t2-t1))
}
