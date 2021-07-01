###################################################################
#-----Using inducing points in local approximate GP framework-----#
###################################################################
eps <- sqrt(.Machine$double.eps)
partial_cov_gen <- utils::getFromNamespace("partial_cov_gen", "hetGP")

## optNegLogLik_theta:
##
## function that optimzes the negative log-likelihood
## with respect to theta (isotropic lengthscale)
optNegLogLik_theta <- function(theta, Xm, X, Y, g, epsQ = 1e-5,
                               epsK = sqrt(.Machine$double.eps), tol = .01,
                               rep_list = NULL){
  # browser()
  f <- function(theta, Xm, X, Y, g, epsQ, epsK, rep_list, theta_prior){

    nll <- negLogLik_theta(theta, Xm, X, Y, g,  epsQ = epsQ,
                           epsK = epsK, rep_list = rep_list,
                           theta_prior = theta_prior)
    its <<- its + 1

    #traced_points <<- rbind(traced_points, c(theta, nll))
    return(nll)
  }
  its <- 0

  out <- try(optimize(f, lower = theta$min, upper = theta$max,
                  tol = tol, Xm = Xm, X = X, Y = Y, g = g, epsQ = epsQ,
                  epsK = epsK, rep_list = rep_list, theta_prior = theta$ab),
             silent=TRUE)

  increase_epsK <- increase_epsQ <- 0
  while (class(out)[1] =='try-error' &
         (epsK*(10^increase_epsK) < 1e-3 | epsQ*(10^increase_epsQ) < 1e-3)) {
    if (epsQ < 1e-3){
      increase_epsQ <- increase_epsQ + 1
      out <- try(optimize(f, lower = theta$min, upper = theta$max,
                          tol = tol, Xm = Xm, X = X, Y = Y, g = g,
                          epsQ = epsQ*(10^increase_epsQ),
                          epsK = epsK, rep_list = rep_list,
                          theta_prior = theta$ab), silent = TRUE)
    } else {
      increase_epsQ <- 0
      increase_epsK <- increase_epsK + 1
      out <- try(optimize(f, lower = theta$min, upper = theta$max,
                          tol = tol, Xm = Xm, X = X, Y = Y, g = g, epsQ = epsQ,
                          epsK = epsK*(10^increase_epsK), rep_list = rep_list,
                          theta_prior = theta$ab), silent = TRUE)

    }

  }

  if (class(out)[1] == 'try-error'){
    stop('Covariance matrix approximation is unstable.',
         ' Consider decreasing n or inducing point design')
  }


  theta = out$minimum; nll = out$objective

  return(list(theta = theta, nll = nll, its = its,
              epsK=epsK*(10^(increase_epsK)),
              epsQ = epsQ*(10^(increase_epsQ))))
}

## optNegLogLik_g
##
## function that optimzes the negative log-likelihood
## with respect to g (nugget)
optNegLogLik_g <- function(nug, Xm, X, Y, theta,
                           epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                           tol = .01, rep_list = NULL){

  deriv <- NULL; its <- 0
  f <- function(nug, Xm, X, Y, theta, epsK, epsQ, rep_list, g_prior){
    nll <- negLogLik_g(nug, Xm, X, Y, theta, epsK = epsK, epsQ = epsQ,
                            rep_list = rep_list, g_prior = g_prior)
    its <<- its + 1

    return(nll)

  }

  out <- try(optimize(f, lower=nug$min, upper=nug$max, tol=tol,
                  Xm=Xm, X=X, Y=Y, theta=theta, epsQ=epsQ,
                  epsK=epsK, rep_list=rep_list, g_prior=nug$ab),
             silent=TRUE)

  ## Increase jitter so that matrices can be inverted and nugget is optimized
  increase_epsK <- increase_epsQ <- 0
  while (class(out)[1]=='try-error' &
         (epsK*(10^increase_epsK) < 1e-3 | epsQ*(10^increase_epsQ) < 1e-3)) {
    if (epsQ < 1e-3){
      increase_epsQ <- increase_epsQ + 1
      out <- try(optimize(f, lower=nug$min, upper=nug$max, tol=tol,
                          Xm=Xm, X=X, Y=Y, theta=theta,
                          epsQ=epsQ*(10^increase_epsQ),
                          epsK=epsK, rep_list=rep_list,
                          g_prior=nug$ab), silent=TRUE)

    } else {
      increase_epsQ <- 0
      increase_epsK <- increase_epsK + 1
      out <- try(optimize(f, lower=nug$min, upper=nug$max, tol=tol,
                          Xm=Xm, X=X, Y=Y, theta=theta, epsQ=epsQ,
                          epsK=epsK*(10^increase_epsK), rep_list=rep_list,
                          g_prior=nug$ab), silent=TRUE)
    }

  }

  if (class(out)[1] == 'try-error')
    stop('Covariance matrix approximation is unstable.',
         ' Consider decreasing n or inducing point design')

  g <- out$minimum; nll <- out$objective

  return(list(g=g, nll=nll, its=its,
              epsK=epsK*(10^(increase_epsK)),
              epsQ=epsQ*(10^(increase_epsQ))))
}

## optNegLogLik_gtheta:
##
## function that optimzes the negative log-likelihood
## with respect to the vector of values g (nugget)
## and theta (lengthscale)
optNegLogLik_g_and_theta <- function(nug, theta, Xm, X, Y,
                                     epsK = sqrt(.Machine$double.eps),
                                     epsQ = 1e-5, tol = .01, rep_list = NULL){

  deriv <- NULL; its <- 0
  f <- function(gtheta_val, Xm, X, Y, g_prior, theta_prior, epsK = epsK,
                epsQ = epsQ, rep_list = rep_list){

    out <- grad_NegLogLik_gtheta(gtheta_val, Xm, X, Y, g_prior=g_prior,
                                  theta_prior=theta_prior,
                                  epsQ=epsQ, epsK=epsK,
                                  rep_list=rep_list)
    its <<- its + 1
    deriv <<-list(x=out$nll, df=out$grad_nll)

    return(out$nll)

  }

  ## derivative read from global variable
  df <- function(gtheta_val, Xm, X, Y, g_prior, theta_prior, epsK = epsK,
                 epsQ = epsQ, rep_list = rep_list) {
    return(deriv$df)
  }

  out <- try(optim(c(nug$start, theta$start), f, df, method='L-BFGS-B',
                   lower=c(nug$min, theta$min),
                   upper=c(nug$max, theta$max),
                   control=list(maxit=100, factr=tol*1e13),
                   Xm=Xm, X=X, Y=Y, g_prior=nug$ab,
                   theta_prior=theta$ab, epsK=epsK, epsQ=epsQ,
                   rep_list=rep_list), silent=TRUE)
  ## Increase jitter so that matrices can be inverted and nugget & theta are
  ## optimized
  increase_epsK <- increase_epsQ <- 1
  while (class(out)[1] =='try-error' &
         (epsK*(10^increase_epsK) < 1e-3 | epsQ*(10^increase_epsQ) < 1e-3)) {
    if (epsQ < 1e-3){
      out <- try(optim(c(nug$start, theta$start), f, df, method='L-BFGS-B',
                       lower=c(nug$min, theta$min),
                       upper=c(nug$max, theta$max),
                       control= list(maxit=100, factr=tol*1e13),
                       Xm=Xm, X=X, Y=Y, g_prior=nug$ab,
                       theta_prior=theta$ab, epsK=epsK,
                       epsQ=epsQ*(10^increase_epsQ),
                       rep_list=rep_list), silent=TRUE)
      increase_epsQ <- increase_epsQ + 1
    } else {
      increase_epsQ <- 1
      out <- try(optim(c(nug$start, theta$start), f, df, method='L-BFGS-B',
                       lower=c(nug$min, theta$min),
                       upper=c(nug$max, theta$max),
                       control= list(maxit=100, factr=tol*1e13),
                       Xm=Xm, X=X, Y=Y, g_prior=nug$ab,
                       theta_prior=theta$ab, epsK=epsK*(10^increase_epsK),
                       epsQ=epsQ, rep_list=rep_list), silent=TRUE)
      increase_epsK <- increase_epsK + 1
    }

  }

  if (class(out)[1] == 'try-error')
    stop('Covariance matrix approximation is unstable.',
         ' Consider decreasing n or inducing point design')
  g <- out$par[1]
  theta <- out$par[2]
  nll <- out$value
  return(list(g=g, theta=theta, nll=nll, its=its,
              epsK=epsK*(10^(increase_epsK-1)),
              epsQ=epsQ*(10^(increase_epsQ-1))))
}



## negLogLik_theta:
##
## function that calculates the negative log-likelihood
## in relation to theta (isotropic lengthscale)
negLogLik_theta <- function(theta, Xm, X, Y, g = sqrt(.Machine$double.eps),
                            epsQ = eps, epsK = sqrt(.Machine$double.eps),
                            theta_prior = NULL, rep_list = NULL){

  if (is.null(dim(X))) X <- as.matrix(X, ncol=length(X))
  if (is.null(dim(Xm))) Xm <- matrix(Xm, ncol=ncol(X))
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)

  KQ <- calcKmQm(Xm=Xm, Xn=X, theta=theta, g, epsQ=epsQ,
                 epsK=epsK, inv=FALSE, mults=rep_list$mult)

  if (is.null(rep_list)){
    N <- nrow(X)
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), KmnOLiY),
                      silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    YSiY <- t(Y) %*% (Y / KQ$OL) -
      (t(Y) %*% KQ$OLiKnm) %*% (QmiKmnOLiY)
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(log(KQ$OL))

  } else {
    N <- sum(rep_list$mult)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')

    YSiY <- t(rep_list$Z) %*% (rep_list$Z /rep(KQ$OL,rep_list$mult)) -
      t(KmnOLiY) %*% QmiKmnOLiY
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(rep(log(KQ$OL), rep_list$mult))
  }

  nlogL <- drop(N/2*(log(2*pi) + log(YSiY) + logdet/N + 1))

  if(!is.null(theta_prior)){
    prior_density <- dgamma(theta, shape=theta_prior[1],
                          rate=theta_prior[2], log=TRUE)
    nlogL <- nlogL  - prior_density
  }

  return(nlogL)
}



## negLogLik_g:
##
## function that calculates the negative log-likelihood
## in relation to g (nugget)
negLogLik_g <- function(g, Xm, X, Y, theta,
                        epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                        g_prior = NULL, rep_list = NULL){

  if (is.null(dim(X))) X <- as.matrix(X, ncol=length(X))
  if (is.null(dim(Xm))) Xm <- matrix(Xm, ncol=ncol(X))
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  N <- nrow(X)

  KQ <- calcKmQm(Xm=Xm, Xn=X, theta=theta, g, epsQ=epsQ,
                 epsK=epsK, inv=FALSE, mults=rep_list$mult)

  if (is.null(rep_list)){
    N <- nrow(X)
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), KmnOLiY),
                      silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    YSiY <- t(Y) %*% (Y / KQ$OL) -
      (t(Y) %*% KQ$OLiKnm) %*% (QmiKmnOLiY)
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(log(KQ$OL))

  } else {
    N <- sum(rep_list$mult)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')

    YSiY <- t(rep_list$Z) %*% (rep_list$Z /rep(KQ$OL,rep_list$mult)) -
      t(KmnOLiY) %*% QmiKmnOLiY
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(rep(log(KQ$OL), rep_list$mult))

  }

  nlogL <- drop(N/2*(log(2*pi) + log(YSiY) + logdet/N + 1))

  if(!is.null(g_prior)){
    prior_density <- dgamma(g, shape=g_prior[1],
                            rate=g_prior[2], log=TRUE)
    nlogL <- nlogL - prior_density
  }

  return(nlogL)
}


## gradKmn_theta:
##
## function that calculates the gradient of the covariance
## matrix K(Xm, X) with respect to a dimension of theta
gradKmn_theta <- function(dimen, Xm, X = NULL, Kmn, theta){
  if (length(unique(theta)) != 1) {
    theta_vec <- vector(length=length(theta))
    theta_vec[dimen] <- theta[dimen]
  } else theta_vec <- theta
  grad.Kmn.th <- partial_cov_gen(Xm, theta_vec, type="Gaussian",
                                 arg="theta_k", X2=X) * Kmn
  return(grad.Kmn.th)
}


## grad_NegLogLik_theta:
##
## function that calculates the gradient of the negative
## log likelihood with respect to theta (isotropic lengthscale)
grad_NegLogLik_theta <- function(theta, Xm, X, Y, g,
                                 epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                                 theta_prior = NULL, rep_list = NULL){
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  dtheta <- length(theta); N <- nrow(X)

  KQ <- calcKmQm(Xm, Xn=X, theta=theta, g=g, epsK=epsK, epsQ=epsQ,
                  mults=rep_list$mult)

  if (is.null(rep_list)){
    N <- nrow(X)
    OLiY <- Y / KQ$OL
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    YSiY <- t(Y) %*% OLiY - t(KmnOLiY) %*% QmiKmnOLiY

  } else {
    N <- length(rep_list$Z)
    AY <- rep_list$mult * Y
    OLiY <- AY/KQ$OL
    OLi0 <- 1/rep(KQ$OL,rep_list$mult)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), KmnOLiY),
                      silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    OLiY0 <- OLi0 * rep_list$Z

    YSiY <- t(rep_list$Z) %*% OLiY0 - t(KmnOLiY) %*% QmiKmnOLiY
  }

  dLL <- vector()
  for (k in 1:dtheta){
    dKm <- gradKmn_theta(k, Xm, NULL, KQ$Km, theta) #M by M matrix
    dKmn <- gradKmn_theta(k, Xm, X, KQ$Kmn, theta) #M by N matrix
    dOL <- -2 * colSums(dKmn * KQ$Kmi_Kmn) +
      colSums(KQ$Kmi_Kmn * (dKm %*% KQ$Kmi_Kmn)) #vector

    if(is.null(rep_list)){
      dKmn_OLiKnm <- dKmn %*% KQ$OLiKnm
      dQm <- dKm + dKmn_OLiKnm + t(dKmn_OLiKnm) -
        t(KQ$OLiKnm) %*% (dOL * KQ$OLiKnm)
    } else {
      dKmn_AOLiKnm <- dKmn %*% (AOLiKnm)
      dQm <- dKm + dKmn_AOLiKnm + t(dKmn_AOLiKnm) -
        t(KQ$OLiKnm) %*% (dOL * AOLiKnm)
    }

    Qmi_dQm <- try(solve(KQ$Qm, dQm), silent=TRUE)
    if (class(Qmi_dQm)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    dOL_OLiY <- dOL * OLiY

    d_OLi <- -t(dOL_OLiY) %*% KQ$OLiKnm %*% QmiKmnOLiY
    d_Kmn <- t(OLiY) %*% t(dKmn) %*% QmiKmnOLiY
    d_Qmi <-  -t(KmnOLiY) %*% Qmi_dQm %*% QmiKmnOLiY
    Kmi_dKm <- solve(KQ$Km, dKm)

    if(is.null(rep_list)){
      dYSiY <- - t(OLiY) %*% dOL_OLiY - (2*d_OLi + 2*d_Kmn + d_Qmi)
      dLL[k] <- .5 * (N*dYSiY/YSiY + sum(diag(Qmi_dQm))  -
                            sum(diag(Kmi_dKm)) + sum(dOL/KQ$OL))
      logdet <- determinant(KQ$Qm)$modulus[1] -
        determinant(KQ$Km)$modulus[1] + sum(log(KQ$OL))
    } else {
      dOL0 <- rep(dOL, rep_list$mult)
      dYSiY <- - t(OLiY0) %*% (dOL0 * OLiY0) - (2*d_OLi + 2*d_Kmn + d_Qmi)
      dLL[k] <- .5 * (N*dYSiY/YSiY + sum(diag(Qmi_dQm)) -
                          sum(diag(Kmi_dKm)) + sum(dOL0*OLi0))
    }

   }

  if (!is.null(theta_prior)){
    grad_prior_density <- (theta_prior[1]-1)/theta - theta_prior[2]
    dLL <- dLL - grad_prior_density
  }

  return(dLL)
}


## grad_NegLogLik_g:
##
## function that calculates the gradient of the negative
## log likelihood with respect to g (nugget)
grad_NegLogLik_g <- function(g, Xm, X, Y, theta,
                             epsK = sqrt(.Machine$double.eps),
                             epsQ = 1e-5, g_prior = NULL, rep_list = NULL){

  if (is.null(dim(X))) X <- as.matrix(X, ncol=length(X))
  if (is.null(dim(Xm))) Xm <- matrix(Xm, ncol=ncol(X))
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)

  KQ <- calcKmQm(Xm=Xm, Xn=X, theta=theta, g,
                 epsQ=epsQ, epsK=epsK, inv=FALSE,
                 mults=rep_list$mult)

  if (is.null(rep_list)){
    N <- nrow(X)
    OLiY <- Y / KQ$OL
    dlogdet_OL <- sum(1/KQ$OL)
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    dQm <- - t(KQ$OLiKnm) %*% KQ$OLiKnm

    YSiY <- t(Y) %*% OLiY - t(KmnOLiY) %*% QmiKmnOLiY
    dYSiY <- - t(OLiY) %*% (OLiY - 2 * KQ$OLiKnm %*% QmiKmnOLiY) -
      t(QmiKmnOLiY) %*% (t(KQ$OLiKnm) %*% KQ$OLiKnm) %*% QmiKmnOLiY
    dlogdet <- sum(diag(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), dQm))) +
      dlogdet_OL

  } else {
    N <- length(rep_list$Z)
    AY <- rep_list$mult * Y
    OLiY <- AY/KQ$OL
    OLi0 <- 1/rep(KQ$OL, rep_list$mult)
    dlogdet_OL <- sum(OLi0)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    OLiY0 <- OLi0 * rep_list$Z
    dQm <- - t(KQ$OLiKnm) %*% (rep_list$mult * KQ$OLiKnm)

    YSiY <- t(rep_list$Z) %*% OLiY0 - t(KmnOLiY) %*% QmiKmnOLiY
    dYSiY <- - t(OLiY0) %*% OLiY0 +
      t(OLiY) %*% (2 * KQ$OLiKnm %*% QmiKmnOLiY) -
      t(QmiKmnOLiY) %*% (t(KQ$OLiKnm) %*% AOLiKnm) %*% QmiKmnOLiY
    dlogdet <- sum(diag(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), dQm))) +
      dlogdet_OL
  }

  dnll_g <- drop(.5*(N*dYSiY/YSiY + dlogdet))

  if(!is.null(g_prior)){
    grad_prior_density <- (g_prior[1]-1)/g - g_prior[2]
    dnll_g <- dnll_g - grad_prior_density
  }

  return(dnll_g)
}


## grad_NegLogLik_gtheta:
##
## function that calculates the gradient of the negative
## log likelihood with respect to the vector gtheta, where
## g is the nugget and theta is the isotropic lengthscale
grad_NegLogLik_gtheta <- function(gtheta, Xm, X, Y, epsQ = 1e-5,
                                  epsK = sqrt(.Machine$double.eps),
                                  rep_list = NULL, g_prior = NULL,
                                  theta_prior = NULL){

  g <- gtheta[1]; theta <- gtheta[2]

  if (!is.null(g_prior)){
    g_prior_density <- dgamma(g, shape=g_prior[1], rate=g_prior[2], log=TRUE)
    grad_g_prior <- (g_prior[1]-1)/g - g_prior[2]
  } else {g_prior_density <- grad_g_prior <- 0}
  if (!is.null(theta_prior)){
    theta_prior_density <- dgamma(theta, shape=theta_prior[1],
                                  rate=theta_prior[2], log=TRUE)
    grad_theta_prior <- (theta_prior[1]-1)/theta - theta_prior[2]
  } else {theta_prior_density <- grad_theta_prior <- 0}

  if (is.null(dim(X))) X <- as.matrix(X, ncol=length(X))
  if (is.null(dim(Xm))) Xm <- matrix(Xm, ncol=ncol(X))
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)

  KQ <- calcKmQm(Xm=Xm, Xn=X, theta=theta, g, epsQ=epsQ, epsK=epsK,
                  mults=rep_list$mult, inv=FALSE)

  ## Initial calcs and gradient wrt g
  if (is.null(rep_list)){
    N <- nrow(X)
    OLiY <- Y / KQ$OL
    dlogdet_OL_g <- sum(1/KQ$OL)
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), KmnOLiY),
                      silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    dQm_g <- - t(KQ$OLiKnm) %*% KQ$OLiKnm

    YSiY <- t(Y) %*% OLiY - t(KmnOLiY) %*% QmiKmnOLiY
    dYSiY_g <- - t(OLiY) %*% (OLiY - 2 * KQ$OLiKnm %*% QmiKmnOLiY) -
      t(QmiKmnOLiY) %*% (t(KQ$OLiKnm) %*% KQ$OLiKnm) %*% QmiKmnOLiY
    dlogdet_g <- sum(diag(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), dQm_g))) +
      dlogdet_OL_g

  } else {
    N <- length(rep_list$Z)
    AY <- rep_list$mult * Y
    OLiY <- AY/KQ$OL
    OLi0 <- 1/rep(KQ$OL,rep_list$mult)
    dlogdet_OL_g <- sum(OLi0)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), KmnOLiY),
                      silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    OLiY0 <- OLi0 * rep_list$Z
    dQm_g <- - t(KQ$OLiKnm) %*% (rep_list$mult * KQ$OLiKnm)

    YSiY <- t(rep_list$Z) %*% OLiY0 - t(KmnOLiY) %*% QmiKmnOLiY
    dYSiY_g <- - t(OLiY0) %*% OLiY0 +
      t(OLiY) %*% (2 * KQ$OLiKnm %*% QmiKmnOLiY) -
      t(QmiKmnOLiY) %*% (t(KQ$OLiKnm) %*% AOLiKnm) %*% QmiKmnOLiY
    dlogdet_g <- sum(diag(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), dQm_g))) +
      dlogdet_OL_g
  }

  dnll_g <- drop(N/2*(dYSiY_g/YSiY + dlogdet_g/N)) - grad_g_prior


  ## Gradient wrt theta
  dKm_th <- gradKmn_theta(1, Xm, NULL, KQ$Km, theta) #M by M matrix
  dK_mn_th <- gradKmn_theta(1, Xm, X, KQ$Kmn, theta) #M by N matrix
  dOL_th <- -2 * colSums(dK_mn_th * KQ$Kmi_Kmn) +
    colSums(KQ$Kmi_Kmn * (dKm_th %*% KQ$Kmi_Kmn)) #vector

  if(is.null(rep_list)){
    dKmn_OLiKnm_th <- dK_mn_th %*% KQ$OLiKnm
    dQm_th <- dKm_th + dKmn_OLiKnm_th + t(dKmn_OLiKnm_th) -
      t(KQ$OLiKnm) %*% (dOL_th * KQ$OLiKnm)
  } else {
    dKmn_AOLiKnm_th <- dK_mn_th %*% (AOLiKnm)
    dQm_th <- dKm_th + dKmn_AOLiKnm_th + t(dKmn_AOLiKnm_th) -
      t(KQ$OLiKnm) %*% (dOL_th * AOLiKnm)
  }

  Qmi_dQm_th <- try(
    solve(KQ$Qm + diag(rep(epsQ, nrow(Xm))), dQm_th),
                    silent=TRUE)
  if (class(Qmi_dQm_th)[1] == 'try-error')
    stop('Q matrix is numerically instable. Consider increasing epsQ.')

  dOL_OLiY_th <- dOL_th * OLiY
  d_OLi_th <- -t(dOL_OLiY_th) %*% KQ$OLiKnm %*% QmiKmnOLiY
  d_Kmn_th <- t(OLiY) %*% t(dK_mn_th) %*% QmiKmnOLiY
  d_Qmi_th <-  -t(KmnOLiY) %*% Qmi_dQm_th %*% QmiKmnOLiY

  Kmi_dKm_th <- solve(KQ$Km, dKm_th)

  if(is.null(rep_list)){
    dYSiY_th <- - t(OLiY) %*% dOL_OLiY_th - (2*d_OLi_th + 2*d_Kmn_th + d_Qmi_th)
    dnll_th <- .5 * (N*dYSiY_th/YSiY + sum(diag(Qmi_dQm_th))  -
                       sum(diag(Kmi_dKm_th)) +
                       sum(dOL_th/KQ$OL)) - grad_theta_prior
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(log(KQ$OL))
  } else {
    dOL0_th <- rep(dOL_th, rep_list$mult)
    dYSiY_th <- - t(OLiY0) %*% (dOL0_th * OLiY0) -
      (2*d_OLi_th + 2*d_Kmn_th + d_Qmi_th)
    dnll_th <- .5 * (N*dYSiY_th/YSiY + sum(diag(Qmi_dQm_th)) -
                       sum(diag(Kmi_dKm_th)) +
                       sum(dOL0_th*OLi0)) - grad_theta_prior
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(rep(log(KQ$OL), rep_list$mult))
  }

  ## Negative log-likelihood
  nll <- drop(N/2*(log(2*pi) + log(YSiY) + logdet/N + 1)) -
    g_prior_density - theta_prior_density

  return( list(nll=nll, grad_nll=c(dnll_g, dnll_th)) )
}




## caclNu:
##
## function that calculates the Maximum Likelihood Estimator
# for nu (scale hyperparameter)
calcNu <- function(Y, Qm, OL, OLiKnm, rep_list = NULL){
  KmnOLiY <- t(OLiKnm) %*% Y
  QmiKOLy <- solve(Qm , KmnOLiY)
  if(!is.null(rep_list)){
    nu <- (sum(rep_list$Z^2 / rep(OL, rep_list$mult)) -
                    t(KmnOLiY) %*% QmiKOLy)/length(rep_list$Z)
  } else {
    nu <- (sum(Y^2 / OL) - t(KmnOLiY) %*% QmiKOLy)/ length(Y)
  }
  return(drop(nu))
}

## calcKmQm:
##
## Calculates the matrices Km (covariance matrix of inducing points)
## and Qm (for prediction) using epsK, epsQ for numerical stability
calcKmQm <- function(Xm, Xn, theta, g = sqrt(.Machine$double.eps),
                     epsQ = 1e-5, epsK = sqrt(.Machine$double.eps),
                     inv = TRUE, mults = NULL){

  eps <- sqrt(.Machine$double.eps)
  M <- nrow(Xm)
  if(M > 1){
    K_m <- cov_gen(Xm, theta=theta) + diag(rep(epsK, M))
  } else {K_m <- as.matrix(1)}
  K_mn <- cov_gen(Xm, Xn, theta=theta)

  if(inv) {Kmi <- solve(K_m)}else Kmi <- NULL
  Kmi_Kmn <- solve(K_m, K_mn)

  Om <- 1 - colSums(K_mn * Kmi_Kmn)
  Om[Om < 0] <- 0

  L  <- g
  OL <- L + Om
  OLiKnm <- 1/OL * t(K_mn)
  if (!is.null(mults)){
    Q_m <- K_m + K_mn %*% (mults * OLiKnm)
  } else {
    Q_m <- K_m + K_mn %*% OLiKnm
  }

  return(list(Kmn=K_mn, Km=K_m, Kmi=Kmi, Qm=Q_m, OLiKnm=OLiKnm,
              L=L, Om=Om, OL=OL, Xm=Xm, epsQ=epsQ, epsK=epsK,
              Kmi_Kmn=Kmi_Kmn, mults=mults))
}


## updateKmQm:
##
## function that using matrix partitioning to efficiently
## update Km & Qm to K_{m+1} & Q_{m+1}
updateKmQm <- function(xm1, Xm, Xn, theta, g = sqrt(.Machine$double.eps), KQ){
  eps <- sqrt(.Machine$double.eps)
  if(is.null(dim(xm1))) xm1 <- matrix(xm1, nrow=1)
  K_m.m1 <- cov_gen(Xm, xm1, theta=theta)
  gam <- KQ$Kmi %*% K_m.m1
  rho <- drop(1 + KQ$epsK - t(K_m.m1) %*% gam)
  Km1 <- rbind(cbind(KQ$Km, K_m.m1),cbind(t(K_m.m1),1 + KQ$epsK))
  Km1i <- rbind(cbind(KQ$Kmi + gam %*% t(gam) / rho, -gam/rho),
                cbind(-t(gam)/rho, 1/rho))
  KmnG <- t(KQ$Kmn) %*% gam
  Kmnew.n <- cov_gen(xm1, Xn, theta=theta)
  Km1n <- rbind(KQ$Kmn, Kmnew.n)

  Om1 <- KQ$Om - (drop(KmnG^2) - 2*colSums(t(KmnG) * Kmnew.n) +
                    drop(Kmnew.n^2))/rho
  OL <- Om1 + KQ$L
  OLi <- 1/OL
  if (is.null(KQ$mults)){
    OLiKnm1 <- t(Km1n) * OLi
  } else {
    OLiKnm1 <- t(Km1n) * (KQ$mults * OLi)
  }
  Q_ms <- KQ$Km + KQ$Kmn %*% OLiKnm1[,1:nrow(Xm)]

  LKnmnew <- OLiKnm1[,nrow(Xm)+1,drop=FALSE]
  psi <- 1 + Kmnew.n %*% LKnmnew
  psi <- psi + KQ$epsQ
  gam <- K_m.m1 + KQ$Kmn %*% LKnmnew
  Qm1 <- rbind(cbind(Q_ms, gam), cbind(t(gam), psi))

  return(list(Kmn=Km1n, Km=Km1, Kmi=Km1i, Qm=Qm1, OLiKnm=OLiKnm1,
              epsQ=KQ$epsQ, epsK=KQ$epsK, L=KQ$L, Om=Om1, OL=OL,
              Xm=rbind(Xm, xm1), mults=KQ$mults))
}


## gradLogLik_xm:
##
## function that calculates the gradient of the log-likelihood
## with respect to one or multiple inducing points being changed
gradLogLik_xm <- function (Xm, X, Y, theta, g, Xm.fixed = NULL){
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  if(is.null(dim(X))) X <- matrix(X, ncol=length(X))
  Xm <- matrix(Xm, ncol=ncol(X))
  if (!is.null(Xm.fixed)){
    Xm.fixed <- matrix(Xm.fixed, ncol=ncol(X))
    Xm.change <- Xm
    Xm <- rbind(Xm.change, Xm.fixed)
  }
  eps <- sqrt(.Machine$double.eps); N <- nrow(X)
  #If gp fit is isotropic, makes theta into a vector of length d
  if(ncol(X) != length(theta)) theta <- rep(theta, ncol(X))

  KQ <- calcKmQm(Xm, X, theta, g)
  K_mi <- solve(KQ$Km + diag(eps, nrow(Xm)))
  KmiKmn <-try(solve(KQ$Km, KQ$Kmn), silent=TRUE)
  if (class(KmiKmn)[1] == 'try-error')
      stop('K matrix is numerically instable. Consider increasing epsK.')

  if (is.null(Xm.fixed)){
    dLogL_mat <- matrix(nrow=nrow(Xm), ncol=ncol(Xm))
    iseq <- 1:nrow(Xm)
  } else {
    dLogL_mat <- matrix(nrow=nrow(Xm.change), ncol=ncol(Xm))
    iseq <- 1:nrow(Xm.change)
  }

  OLiKnm <- t(KQ$Kmn) / KQ$OL
  Liy <- Y / KQ$OL

  OLiKnmy <- OLiKnm %*% Y
  QmiOLiKnmy <- try(solve(KQ$Qm, OLiKnmy), silent=TRUE)
  if(class(QmiOLiKnmy)[1] == 'try-error')
    stop('Q matrix is numerically instable. Consider increasing epsQ.')

  Siy <- Liy - OLiKnm %*% QmiOLiKnmy
  CCSy <- KmiKmn %*% Siy
  ySiy <- t(Y) %*% Siy

  ## Creates a matrix of partial derivatives with
  ## same dimensions as Xm
  for (i in iseq){
    for (di in 1:ncol(Xm)){
      if (nrow(Xm) > 1) {
        dC <- gradKmi_xmi(i, di, Xm, theta, K_mi)
        dK_mi <- dC$dKmi
        dK_m <- matrix(0,nrow=nrow(Xm), ncol=nrow(Xm))
        dK_m[i,] <- dK_m[,i] <- dC$dKm
      } else {dK_mi <- 0}
      dK_mn <- matrix(gradKmn_xmi(i, di, Xm, X, theta), nrow=1) # 1 by N matrix of nonzero row of dK_mn
      dO_m <- -2 * as.vector(dK_mn) * KmiKmn[i,] - colSums(KQ$Kmn * (dK_mi %*% KQ$Kmn)) # vector

      part1a <- matrix(rep(0, nrow(Xm)))
      part1a[i,1] <- dK_mn %*% Siy
      part1 <- t(part1a) %*% CCSy
      dySiy <- -( t(Siy) %*% (dO_m * Siy) + 2 * part1 - t(CCSy) %*% dK_m %*% CCSy )
      part2 <- matrix(0, nrow=nrow(Xm), ncol=nrow(Xm))
      part2[i, ] <-  dK_mn %*% OLiKnm
      dQm <- dK_m + part2 + t(part2) - t(OLiKnm) %*% (dO_m * OLiKnm)

      Qmi.dQm <- try(solve(KQ$Qm, dQm), silent=TRUE)
      if (class(Qmi.dQm)[1] == 'try-error')
        stop('Q matrix is numerically instable. Consider increasing epsQ.')

      dl1 <- sum(diag(Qmi.dQm))
      dl2 <- sum(KQ$Kmi * dK_m)
      dl3 <- sum(dO_m/KQ$OL)
      dLogL_mat[i, di] <- .5 * (N*dySiy/ySiy + dl1 - dl2 + dl3)
    }

  }
  return(dLogL_mat)
}


## negLogLik:
##
## function that calculates the negative log-likelihood given
## an inducing point design. Uses the Woodbury formulas.
negLogLik <- function(Xm, X, Y, theta, g=1e-4, rep_list = NULL, epsQ = 1e-5){
  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  if(is.null(dim(X))) X <- as.matrix(X, ncol=length(X))
  Xm <- matrix(Xm, ncol=ncol(X))
  N <- nrow(X)

  KQ <-calcKmQm(Xm, X, theta, g, mults=rep_list$mult)

  if(is.null(rep_list)){
    OLiy <- Y / KQ$OL
    KmnOLiy <- KQ$Kmn %*% OLiy

    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(log(KQ$OL))

    QmiKmnOLiy <- try(solve(KQ$Qm, KmnOLiy), silent=TRUE)
    if (class(QmiKmnOLiy)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')
    nlogL <-as.numeric(N/2*(log(2*pi)+log(t(Y) %*% OLiy -
                             t(KmnOLiy) %*% QmiKmnOLiy) + logdet/N + 1))
  } else {
    N <- sum(rep_list$mult)
    AOLiKnm <- rep_list$mult * KQ$OLiKnm
    KmnOLiY <- t(AOLiKnm) %*% Y
    QmiKmnOLiY <- try(solve(KQ$Qm + diag(rep(epsQ, nrow(Xm)))
                            , KmnOLiY), silent=TRUE)
    if (class(QmiKmnOLiY)[1] == 'try-error')
      stop('Q matrix is numerically instable. Consider increasing epsQ.')

    YSiY <- t(rep_list$Z) %*% (rep_list$Z /rep(KQ$OL,rep_list$mult)) -
      t(KmnOLiY) %*% QmiKmnOLiY
    logdet <- determinant(KQ$Qm)$modulus[1] -
      determinant(KQ$Km)$modulus[1] + sum(rep(log(KQ$OL), rep_list$mult))

    nlogL <- drop(N/2*(log(2*pi) + log(YSiY) + logdet/N + 1))
  }

  return(nlogL)
}



## boundBox:
##
## function used to create a bounding box for multi-start points
## used for sequential optimization of inducing points
boundBox <- function(Xref, Xm, lower, upper, perc = NULL){

  ## Calculates % of design space size in each direction (+/-)
  window <- perc * max(upper-lower)
  windL <- Xref - window; windU <- Xref + window

  diff <- abs(sweep(Xm, 2, Xref)) ## Distances for each dim of Xm from Xref
  max.diffs <- apply(diff, 2, max) ## Max differences
  lowlim <- uplim <- vector()
  ## Ensures boundaries don't extend beyond design space
  for (j in 1:ncol(Xm)){
    lowlim[j] <- max(lower[j], Xref[,j] - max.diffs[j])
    uplim[j] <- min(upper[j], Xref[,j] + max.diffs[j])
  }
  return(list(lower=lowlim, upper=uplim,
              windL=as.vector(windL), windU=as.vector(windU)))
}

## gradOmega_xmi:
##
## function that calculates the gradient of the Omega vector
## with respect to a dimension of the last inducing point
gradOmega_xmi <-function(i, Kmn, Kmi, dKnm, dKmi){
  # part1 <- Kmi[,i,drop=FALSE] %*% t(dKnm) #M+1 by n matrix
  # part2 <- dKmi %*% Kmn #M+1 by n matrix
  part3 <- (2*Kmi[,i,drop=FALSE]) %*% t(dKnm) + dKmi %*% Kmn #M+1 by n matrix
  d.om <- -colSums(Kmn * part3) #Equivalent to diag(Knm %*% part3)

  return(d.om)
}


## gradQm_xmi:
##
## function that calculates the gradient of the Q matrix
## with respect to a dimension of the last inducing point
gradQm_xmi <- function(i, di, Xm, Xn, theta, KQ, dC){
  dK_mn <- gradKmn_xmi(i, di, Xm, Xn, theta) #vector of length n

  Kmn.OLi.dKnm <- t(KQ$OLiKnm) %*% matrix(dK_mn) #col vector for col i in M+1 by M+1 Matrix
  K_mi <- KQ$Kmi #already included jitter + diag(eps, nrow(Xm)))
  dOmega <- gradOmega_xmi(i, KQ$Kmn, K_mi, dK_mn, dC$dKmi) #vector of length n

  part2 <- -t(KQ$OLiKnm) %*% (dOmega * KQ$OLiKnm)
  part2[i,] <- part2[i,] + as.vector(Kmn.OLi.dKnm + dC$dKm)
  part2[,i] <- part2[,i] + as.vector(Kmn.OLi.dKnm + dC$dKm)
  return(part2)
}


## gradKmn_theta:
##
## function that calculates the gradient of the covariance
## matrix K(Xm, X) with respect to theta in a particular dimension
gradKmn_theta <- function(dimen, Xm, X = NULL, Kmn, theta){
  if (length(unique(theta)) != 1) {
    theta_vec <- vector(length=length(theta))
    theta_vec[dimen] <- theta[dimen]
  } else theta_vec <- theta
  grad.Kmn.th <- partial_cov_gen(Xm, theta_vec, type="Gaussian",
                                 arg="theta_k", X2=X) * Kmn
  return(grad.Kmn.th)
}


## gradKern_xmi:
##
## function that calculates the gradient of the Gaussian kernel
## with respect to a dimension of Xmi
gradKern_xmi <- function(dimen, X, Xmi, theta){
  # if (length(theta)!=length(Xmi)) theta <- rep(theta, length(Xmi))
  final <- -2 * (Xmi[dimen]-X[,dimen])/theta *
    exp(-rowSums((X/sqrt(theta) - rep(Xmi/sqrt(theta),
                                      rep.int(nrow(X), length(Xmi))))^2))
  return(final)
}


##  Derivative of K_M wrt Xm[m,]
## gradKmi_xmi:
##
## function that calculates and returns the gradient of the inverse
## of Km with respect to Xm[m,] in a particular dimension
gradKmi_xmi <- function (m, dimen, Xm, theta, Kmi){
  if(is.null(dim(Xm))) Xm <- matrix(Xm, ncol=length(Xm))

  dKm_xm <- as.matrix(gradKern_xmi(dimen, Xm, Xm[m,], theta)) #vector of derivatives of K_m wrt inducing point m

  dKm.Kmi <- rbind(dKm_xm[-m,,drop=FALSE] %*% Kmi[m,],t(dKm_xm) %*% Kmi)

  return(list(dKmi=-Kmi %*% dKm.Kmi, dKm=dKm_xm))
}


## gradKmn_xmi:
##
## function that calculates and returns the gradient of the covariance
## matrix K(Xm, Xn) with respect to Xm[m,] in a particular dimension
gradKmn_xmi <- function(m, dimen, Xm, Xn, theta){
  if(is.null(dim(Xm))) Xm <- matrix(Xm, ncol=length(Xm))

  K_nmi <- vector(length=nrow(Xn))
  dK_m.n <- gradKern_xmi(dimen, Xn, Xm[m,], theta)

  return(dK_m.n)
}
