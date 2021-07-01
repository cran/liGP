eps <- sqrt(.Machine$double.eps)
## loadX:
##
## bulk copy a big X over to the C-side so that
## multiple closest calls on it can be quick
#' @useDynLib liGP loadX_R
loadX <- function(X, start=50)
{
  out <- .C("loadX_R",
            m=as.integer(ncol(X)),
            n=as.integer(nrow(X)),
            X=as.double(t(X)))
}




#' @useDynLib liGP unloadX_R
unloadX <- function()
{
  .C("unloadX_R")
}

## closestIndices:
##
## R interface to closest_indices_R, which uses sorting or quick select
## in order to facilitate laGP calculations on local subsets; this
## interface is primarily provided for debugging purposes

closest <- function(Xref, n=50)
{
  out <- .C("closest_R",
            m=as.integer(ncol(Xref)),
            start=as.integer(n),
            Xref=as.double(t(Xref)),
            nref=as.integer(nrow(Xref)),
            oD=integer(n))

  return(out$oD)
}


## liGP:
##
## Produces independent predictions for XX in parallel using local induced GPs
## with a single inducing point design template
liGP <- function(XX, X = NULL, Y = NULL, Xm.t, N, g = 1e-6,
                 theta = NULL, nu = NULL, num_thread = 1,
                 epsK = sqrt(.Machine$double.eps),
                 epsQ = 1e-5, tol = .01, reps = FALSE, Xni.return = FALSE){
  ## Collects data from X,Y or reps_list
  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if (reps){
    Xorig <- X; Y_orig <- Y
    reps <- find_reps(X, Y)
    X <- reps$X0
    Y <- reps$Z0
  } else reps <- FALSE

  ###### Sanity checks ######
  if(ncol(X)!=ncol(XX)) ## Dimension check
    stop('X and XX have a different number of columns')
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if (min(Xm.t) >0 || max(Xm.t) < 0) ## Attempt to verify input is a template
    warning('Xm.t doesn\'t seem to be a template centered at the origin')
  if(ncol(Xm.t)!=ncol(X)) ## Dimension check
    stop('X and Xm.t have a different number of columns')
  if(nrow(Xm.t) > N) ## Number of IP doesn't exceed neighborhood size
    stop('Number of inducing points (nrow(Xm.t)) exceeds neighborhood size (N)')
  if(N > nrow(X)) ## Neighborhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  if (!is.null(nu))
    if(nu <=0)
      stop('nu should be positive')
  if(num_thread %% 1 != 0 | num_thread < 1)
    stop('num_thread is not a positive integer')
  if (num_thread > nrow(XX)){
    warning(paste("num_thread exceeds number of predictive locations. num_thread set to",nrow(XX)))
    num_thread <- nrow(XX)
  }
  available_cores <- detectCores()
  if (num_thread > available_cores){
    warning(paste("num_thread exceeds number of available cores. num_thread set to",available_cores))
    num_thread <- available_cores
  }
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')
  if (tol <= 0)
    stop('tol should be a positive number')

  ## Initializes theta, g with same methods as laGP
  if (is.null(theta)) {
    Xc <- matrix(apply(X, 2, median), nrow=1)
    Xc_to_X_dists <- distance(Xc, X)
    quant_dists <- quantile(Xc_to_X_dists, N/nrow(X))
    Xn <- X[Xc_to_X_dists < quant_dists,,drop=FALSE]
    theta <- darg(NULL, Xn)
  }
  if (is.null(g))
      g <- ifelse(is.list(reps),
                  garg(list(mle=TRUE), reps$Z),
                  garg(list(mle=TRUE), Y))

  if(is.list(reps)) {X <- Y <- NULL}

  ## For timing
  t1 <- proc.time()[3]

  ## Makes predictions for each row of XX
  if (num_thread == 1){
    calc_mle <- ifelse(!is.list(g) & !is.list(theta), FALSE, TRUE)

    ligp_preds <- liGP.forloop(XX, X=X, Y=Y, Xm.t=Xm.t, N=N,
                       theta=theta, g=g, nu=nu, epsK=epsK,
                       epsQ=epsQ, tol=tol, reps=reps, Xni.return=Xni.return)
    out <- list(mean=ligp_preds$mean, var=ligp_preds$var, nu=ligp_preds$nu,
                g=g, theta=theta, Xm.t=Xm.t, eps=ligp_preds$eps)

    if(calc_mle)  out$mle <- ligp_preds$mle
    if(Xni.return)  out$Xni <- ligp_preds$Xni

  } else {
    numPred <- nrow(XX)
    list_to_hold_results <- list(list(), list(), list(), list(), list(),
                                list(), list())

    if (Xni.return) list_to_hold_results[[8]] <- list()

    if(!is.list(g) & !is.list(theta)){
      calc_mle <- FALSE
    } else {
      calc_mle <- TRUE
      if(is.list(g) & is.list(theta)){
        pred.mle <- matrix(nrow=numPred, ncol=3)
      } else pred.mle <- matrix(nrow=numPred, ncol=2)
      list_to_hold_results[[length(list_to_hold_results)+1]] <- list()
    }

    ## Splits XX into num_thread groups for efficient use of multiple GPUs
    groups <- split(sample(1:numPred), sample(1:numPred) %% num_thread)

    cl <- makeCluster(num_thread)
    registerDoParallel(cl)
    results <- foreach(i=1:num_thread, .errorhandling="pass", .combine='rbind',
                       .multicombine=TRUE, .packages=c('liGP','laGP', 'hetGP'),
                       .init=list_to_hold_results) %dopar% {
                         liGP.forloop(XX=XX[unlist(groups[i]),,drop=FALSE],
                                      X=X, Y=Y, Xm.t=Xm.t, N=N,
                                      theta=theta, g=g, nu=nu,
                                      epsQ=epsQ, epsK=epsK, tol=tol,
                                      reps=reps, Xni.return=Xni.return)
                       }
    stopCluster(cl)

    ## Collects results from different clusters
    pred.mean <- pred.var <- pred.nu <- rep(0, nrow(XX))
    if(Xni.return) Xni <- matrix(nrow=nrow(XX), ncol=N)
    eps.mat <- matrix(nrow=nrow(XX), ncol=2)
    colnames(eps.mat) <- c('epsK', 'epsQ')

    for (i in 1:num_thread){
      group_rows <- as.vector(unlist(groups[i]))
      pred.mean[group_rows] <- unlist(results[i+1,1])
      pred.var[group_rows] <- unlist(results[i+1,2])
      pred.nu[group_rows] <- unlist(results[i+1,3])
      eps.mat[group_rows,] <- as.matrix(results[i+1,7][[1]])
      if(calc_mle){
        gtheta_mat <- as.matrix(results[i+1,8][[1]])
        pred.mle[group_rows,] <- gtheta_mat
        if(Xni.return){
          Xni_mat <- as.matrix(results[i+1,9][[1]])
          Xni[group_rows,] <- Xni_mat
        }
      } else {
        if(Xni.return){
          Xni_mat <- as.matrix(results[i+1,8][[1]])
          Xni[group_rows,] <- Xni_mat
        }
      }
    }

    out <- list(mean=pred.mean, var=pred.var, nu=pred.nu,
                g=g, theta=theta, Xm.t=Xm.t, eps=eps.mat)

    ## Only returns mle if theta and/or g is optimized
    if(calc_mle){
      cols <- c('g', 'theta', 'its')
      col_select <- c(is.list(g), is.list(theta), TRUE)
      colnames(pred.mle) <- cols[col_select]
      out$mle <- as.data.frame(pred.mle)
    }
    if(Xni.return) out$Xni <- Xni
  }

  t2 <- proc.time()[3]

  out$time <- t2 - t1
  return(out)
}


## liGP.forloop: returns mean and variance predictions for Xref through a loop
##
## Produces independent predictions for XX, including local induced GPs
## with a single inducing point design template
liGP.forloop <- function(XX, X = NULL, Y = NULL, Xm.t, N, g = 1e-6, theta = NULL,
                 nu = NULL, epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                 tol = .01, reps = FALSE, Xni.return = FALSE){
  ## Collects data from X,Y or reps_list
  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if (reps){
    Xorig <- X; Y_orig <- Y
    reps <- find_reps(X, Y)
    X <- reps$X0
    Y <- reps$Z0
  } else reps <- FALSE

  ## Initialize g & theta or perform input checks
  if(!is.list(theta)){
    if(is.null(theta)){
      Xc <- matrix(apply(X, 2, median), nrow=1)
      Xc_to_X_dists <- distance(Xc, X)
      quant_dists <- quantile(Xc_to_X_dists, N/nrow(X))
      Xn <- X[Xc_to_X_dists < quant_dists,,drop=FALSE]
      theta <- unname(quantile(dist(Xn),.1)^2)
    }

    if(length(theta) == 1) {
      theta0 <- rep(theta, nrow(XX))
    } else if(length(theta) == nrow(XX)){
      theta0 <- theta
    } else{
      stop('An incorrect number of theta values was supplied.')
    }
  }
  if (is.null(g))
    g <- ifelse(is.list(reps),
                garg(list(mle=TRUE), reps$Z),
                garg(list(mle=TRUE), Y))
  if(!is.list(g)){
    if(length(g) == 1) {
      g0 <- rep(g, nrow(XX))
    } else if(length(g) == nrow(XX)){
      g0 <- g
    } else{
      stop('An incorrect number of g values was supplied.')
    }
  }

  if (is.null(dim(XX))) XX <- matrix(XX, nrow=1)
  if (!is.vector(Y)) Y <- as.vector(Y)
  #### Sanity checks ####
  if(ncol(X)!=ncol(XX)) ## Dimension check
    stop('X and XX have a different number of columns')
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if(ncol(Xm.t)!=ncol(X))
    stop('X and Xm.t have a different number of columns')
  if(nrow(Xm.t) > N) ## Number of IP doesn't exceed neighborhood size
    stop('Number of inducing points (nrow(Xm.t)) exceeds neighborhood size (N)')
  if(N > nrow(X)) ## Neighborhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  if (!is.null(nu))
    if(nu <=0)
      stop('nu should be positive')
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')
  if (tol <= 0)
    stop('tol should be a positive number')

  ##Storage of results
  nrow_xx <- nrow(XX)
  mean_vec <- var_vec <- nu_vec <- rep(0, nrow_xx)
  if(is.list(g)){
    if(is.list(theta)){
      mle <- data.frame(g=rep(0, nrow_xx), theta=rep(0, nrow_xx),
                       its=rep(0, nrow_xx))
    } else mle <- data.frame(g=rep(0, nrow_xx), its=rep(0, nrow_xx))
  } else if(is.list(theta)) mle <- data.frame(theta=rep(0, nrow_xx),
                                              its=rep(0, nrow_xx))
  if (Xni.return) Xni <- matrix(nrow=nrow(XX), ncol=N)

  ## Calls R-side function to load X into memory
  if (N < nrow(X)) loadX(X)

  ## For timing
  t1 <- proc.time()[3]

  for (i in 1:nrow(XX)){

    ## Builds local neighborhood for XX[i,]
    if (N < nrow(X)){
      cI <- closest(XX[i,,drop=FALSE], n=N)
      Xn <- X[cI,,drop=FALSE]; Yn <- Y[cI]
      if(is.list(reps)){
        reps_n_list <- list(mult=reps$mult[cI],
                         Z=matrix(c(unlist(reps$Zlist[cI]))))
      }
      if (Xni.return)  Xni[i,] <- (1:nrow(X))[cI]
    } else {
      Xn <- X; Yn <- Y
      if(is.list(reps))
        reps_n_list <- list(mult=reps$mult, Z=matrix(reps$Z))
      if (Xni.return)  Xni[i,] <- 1:nrow(X)
    }
    if (!is.list(reps)) reps_n_list <- NULL

    if(is.null(dim(Yn))) Yn <- matrix(Yn, ncol=1)

    ## Centers inducing point design at XX[i,]
    Xm <- Xm.t
    for(j in 1:ncol(Xm.t))
      Xm[,j] <- Xm.t[,j] + XX[i,j]

    ## Optimizes hyperparameter(s) with log-likelihood if
    ## optimization bounds provided
    if(is.list(theta)){
      if(is.list(g)){
        out <- optNegLogLik_g_and_theta(g, theta, Xm, Xn, Yn,
                                        epsQ=epsQ, epsK=epsK,
                                        tol=tol, rep_list=reps_n_list)
        mle[i,] <- c(out$g, out$theta, out$its)

      } else{
        out <- optNegLogLik_theta(theta, Xm, Xn, Yn, g0[i],
                                  epsQ=epsQ, epsK=epsK,
                                  tol=tol, rep_list=reps_n_list)
        mle[i,] <- c(out$theta, out$its)
        out$g <- g0[i]
      }
    } else {
      if(is.list(g)){
        out <- optNegLogLik_g(g, Xm, Xn, Yn, theta=theta0[i],
                              epsQ=epsQ, epsK=epsK, tol=tol,
                              rep_list=reps_n_list)
        mle[i,] <- c(out$g, out$its)
        out$theta <- theta0[i]
      } else {
        out <- list(g=g0[i], theta=theta0[i], its=0)
      }
    }

    ## Builds matrices for prediction
    KQ <- calcKmQm(Xm, Xn, out$theta, g=out$g, epsK=epsK, epsQ=epsQ,
                    inv=FALSE, mults=reps_n_list$mult)
    K_mXX <- cov_gen(Xm, XX[i,,drop=FALSE], theta=out$theta)
    QmiKXX <- try(solve(KQ$Qm, K_mXX), silent=TRUE)
    KmiKXX <- try(solve(KQ$Km, K_mXX), silent=TRUE)

    if(is.null(reps_n_list)){
      KmnOLiY <- t(KQ$OLiKnm) %*% Yn
    } else {
      KmnOLiY <- t(KQ$OLiKnm) %*% (reps_n_list$mult * Yn)
    }
    QmiKOLy <- try(solve(KQ$Qm, KmnOLiY), silent=TRUE)

    ## Increases epsQ or epsK if needed for matrix inversions
    eps_new <- matrix(c(epsK, epsQ), nrow=nrow(XX), ncol=2, byrow=TRUE)
    if (class(QmiKXX)[1]== "try-error" | class(KmiKXX)[1]== "try-error" |
        class(QmiKOLy)[1]== "try-error" ){
      ## Takes epsK that was adjusted in hyperparameter optimization
      if(exists('out$epsK')){
        eps_new[i,] <- c(out$epsK, epsQ)
        KQ <- calcKmQm(Xm, Xn, out$theta, g=out$g, epsK=out$epsK, epsQ=out$epsQ,
                       inv=FALSE, mults=reps_n_list$mult)
        K_mXX <- cov_gen(Xm, XX[i,,drop=FALSE], theta=out$theta)
        QmiKXX <- try(solve(KQ$Qm + diag(out$epsQ, nrow(Xm)), K_mXX),
                      silent=TRUE)
        KmiKXX <- try(solve(KQ$Km, K_mXX), silent=TRUE)

        if(is.null(reps_n_list)){
          KmnOLiY <- t(KQ$OLiKnm) %*% Yn
        } else {
          KmnOLiY <- t(KQ$OLiKnm) %*% (reps_n_list$mult * Yn)
        }
        QmiKOLy <- try(solve(KQ$Qm + diag(out$epsQ, nrow(Xm)), KmnOLiY),
                       silent=TRUE)
      } else {
        ## Incrementally increases jitter
        increase_epsK <- increase_epsQ <- 0
        matrix_inverse_calcs <- function(epsK, epsQ){
          KQ <- calcKmQm(Xm, Xn, out$theta, g=out$g, epsK=epsK, epsQ=epsQ,
                         inv=FALSE, mults=reps_n_list$mult)
          K_mXX <- cov_gen(Xm, XX[i,,drop=FALSE], theta=out$theta)
          QmiKXX <- try(solve(KQ$Qm + diag(epsQ, nrow(Xm)), K_mXX),
                         silent=TRUE)
          KmiKXX <- try(solve(KQ$Km, K_mXX), silent=TRUE)

          if(is.null(reps_n_list)){
            KmnOLiY <- t(KQ$OLiKnm) %*% Yn
          } else {
            KmnOLiY <- t(KQ$OLiKnm) %*% (reps_n_list$mult * Yn)
          }
          QmiKOLy <- try(solve(KQ$Qm + diag(epsQ, nrow(Xm)), KmnOLiY),
                          silent=TRUE)
          return(list(QmiKXX = QmiKXX, KmiKXX = KmiKXX, QmiKOLy = QmiKOLy))
        }

        while ((class(QmiKXX)[1]=="try-error" | class(KmiKXX)[1]=="try-error" |
                class(QmiKOLy)[1]=="try-error")
               & ((epsK*(10^increase_epsK) < 1e-3) |
                  (epsQ*(10^increase_epsQ) < 1e-3))) {
          if (epsQ*(10^increase_epsQ) < 1e-3){
            increase_epsQ <- increase_epsQ + 1
            new_matrices <- matrix_inverse_calcs(epsK=epsK,
                                                 epsQ=epsQ*(10^increase_epsQ))
            QmiKXX <- new_matrices$QmiKXX; KmiKXX <- new_matrices$KmiKXX
            QmiKOLy <- new_matrices$QmiKOLy
          } else {
            increase_epsK <- increase_epsK + 1
            new_matrices <- matrix_inverse_calcs(epsK*(10^increase_epsK),
                                                 epsQ=epsQ)
            QmiKXX <- new_matrices$QmiKXX; KmiKXX <- new_matrices$KmiKXX
            QmiKOLy <- new_matrices$QmiKOLy
          }

        }

        if (increase_epsK > 1){
          eps_new[i,] <- c(epsK*(10^increase_epsK), epsQ)
        } else {
          eps_new[i,] <- c(epsK, epsQ*(10^increase_epsQ))
        }
      }

    }

    if (class(QmiKXX)[1]=="try-error" | class(KmiKXX)[1]=="try-error" |
        class(QmiKOLy)[1]=="try-error" )
      stop('Covariance matrix approximation is unstable.',
           ' Consider decreasing N or the inducing point design size')

    ## Predicted mean & variance for XX[i,]
    ## Calculates MLE value of nu if not fixed
    if (is.null(reps_n_list)){
      mean_XX <- t(QmiKXX) %*% t(KQ$OLiKnm) %*% Yn
      if (is.null(nu)){
        nu_vec[i] <- (sum(Yn^2 / KQ$OL) - t(KmnOLiY) %*% QmiKOLy)/length(Yn)
      } else nu_vec[i] <- nu
    } else {

      mean_XX <- t(QmiKXX) %*% t(KQ$OLiKnm) %*% (reps_n_list$mult * Yn)
      if (is.null(nu)){
        nu_vec[i] <- (sum(reps_n_list$Z^2 / rep(KQ$OL, reps_n_list$mult)) -
                        t(KmnOLiY) %*% QmiKOLy)/length(reps_n_list$Z)
      } else nu_vec[i] <- nu
    }

    var_vec[i] <- 1 + out$g - colSums(K_mXX * (KmiKXX - QmiKXX))
    mean_vec[i] <- mean_XX

  }

  ## R-side function to remove X from memory
  if (N < nrow(X)) unloadX()

  t2 <- proc.time()[3]
  out_list <- list(mean=mean_vec, var=var_vec, nu=nu_vec, g=g, theta=theta,
                   time=t2-t1, eps=eps_new)

  if (is.list(theta) | is.list(g)) out_list$mle <- mle
  if (Xni.return) out_list$Xni <- Xni

  return(out_list)
}


## loiGP:
##
## Generates local optimal inducing point designs for each XX[i,] based on the
## neighborhood. Then using inducing point design and local neighborhood to
## generate prediction. Can be run in parallel.
loiGP <- function(XX, X = NULL, Y = NULL, M, N, g = 1e-6, theta = NULL,
                  nu = NULL, method = c('wimse','alc'),
                  integral_bounds = NULL, num_thread = 1,
                  epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01,
                  reps = FALSE){
  ## Collects data from X,Y or reps_list
  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if (reps){
    Xorig <- X; Y_orig <- Y
    reps <- find_reps(X, Y)
    X <- reps$X0
    Y <- reps$Z0
  } else reps=FALSE

  ###### Sanity checks ######
  if(ncol(X)!=ncol(XX)) ## Dimension check
    stop('X and XX have a different number of columns')
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if(N > nrow(X)) ## Neighorhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  if(M > N) ## Number of inducing points should be smaller than neighborhood size
    warning('M is greater than N')
  if (!is.null(nu))
    if(nu <=0)
      stop('nu should be positive')
  if(!method %in% c('wimse', 'alc'))
    stop('A valid method was not given. Choices include: wimse and alc')
  if(method =='wimse'){
    Xrange <- apply(X, 2, range)
    if (!is.null(integral_bounds))
      if(sum(Xrange[1,] < integral_bounds[1,]) > 0 || sum(Xrange[2,] > integral_bounds[2,]) > 0)
        stop('X outside integration bounds')
    if(is.null(integral_bounds))
      integral_bounds <- Xrange
  }

  ## Checks that number of threads is valid
  if(num_thread %% 1 != 0 | num_thread < 1)
    stop('num_thread is not a positive integer')
  if (num_thread > nrow(XX)){
    warning(paste("num_thread exceeds number of predictive locations. num_thread set to",
                  nrow(XX)))
    num_thread <- nrow(XX)
  }
  available_cores <- detectCores()
  if (num_thread > available_cores){
    warning(paste("num_thread exceeds number of available cores. num_thread set to",
                  available_cores))
    num_thread <- available_cores
  }
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')
  if (tol <= 0)
    stop('tol should be a positive number')

  ## For timing
  t1 <- proc.time()[3]

  ## Builds unique inducing point designs and makes predictions for each row of XX
  i <- NULL
  if (method == 'alc'){
    cl <- makeCluster(num_thread)
    registerDoParallel(cl)
    results <- foreach(i=1:nrow(XX), .errorhandling="pass", .combine='rbind',
                       .multicombine=TRUE, .packages=c('liGP','laGP', 'hetGP'),
                       .init=list(list(),list(),list(),list(),list(),list(),list())) %dopar% {
                         XXi <- XX[i,,drop=FALSE]
                         ## Builds neighborhood for XXi
                         XXi_to_X_dists <- distance(XXi, X)
                         dist_quant <- quantile(XXi_to_X_dists, N/nrow(X))
                         Xn <- X[XXi_to_X_dists < dist_quant,,drop=FALSE]
                         Yn <- Y[XXi_to_X_dists < dist_quant]
                         if (reps){
                           reps_n_list <- list(mult=reps$mult[XXi_to_X_dists < dist_quant],
                                            Z=matrix(c(unlist(reps$Zlist[XXi_to_X_dists < dist_quant]))))
                         } else reps_n_list <- NULL

                         ## Initializes theta based on local neighborhood if not provided
                         if (is.null(theta)) theta <- quantile(dist(Xn), .1)^2
                         if (is.null(g)) g <- laGP::garg(list(mle=TRUE), Yn)

                         ## Sets initial values for inducing point optimization within neighborhood area
                         ip_bounds <- apply(Xn, 2, range)

                         ## Builds locally optimal inducing point design with ALC
                         optXm <- optIP.ALC(XX[i,,drop=FALSE], Xref=NULL, M=M,
                                            Xn=Xn, Yn=Yn, theta=theta, g=g,
                                            ip_bounds=ip_bounds, verbose=FALSE)

                         ## Optimizes hyperparameters and calculates XXi prediction
                         preds <- giGP(XX[i,,drop=FALSE], Xm=optXm$Xm, Xn, Yn,
                                       theta=theta, g=g, nu=nu, epsQ=epsQ,
                                       epsK=epsK, tol=tol, reps=reps_n_list)

                         list(mean=preds$mean, var=preds$var, theta=theta,
                              nu=preds$nu, Xm=optXm$Xm, eps = preds$eps,
                              mle=preds$mle)
                       }
    stopCluster(cl)
  } else if(method == "wimse"){
    cl <- makeCluster(num_thread)
    registerDoParallel(cl)
    results <- foreach(i=1:(nrow(XX)), .errorhandling="pass", .combine='rbind',
                       .multicombine=TRUE, .packages=c('liGP','laGP', 'hetGP'),
                       .init=list(list(),list(),list(),list(),list(), list(), list())) %dopar% {

                         XXi <- XX[i,,drop=FALSE]
                         ## Builds neighborhood for XXi
                         XXi_to_X_dists <- distance(XXi, X)
                         dist_quant <- quantile(XXi_to_X_dists, N/nrow(X))
                         Xn <- X[XXi_to_X_dists < dist_quant,,drop=FALSE]
                         Yn <- Y[XXi_to_X_dists < dist_quant]
                         if (reps){
                           reps_n_list <- list(mult=reps$mult[XXi_to_X_dists < dist_quant],
                                            Z=matrix(c(unlist(reps$Zlist[XXi_to_X_dists < dist_quant]))))
                         } else reps_n_list <- FALSE

                         ## Initializes theta based on local neighborhood
                         ## if not provided
                         if (is.null(theta)) theta <- quantile(dist(Xn), .1)^2
                         if (is.null(g)) g <- garg(list(mle=TRUE), Yn)

                         ## Sets initial values for inducing point optimization
                         ## within neighborhood area
                         ip_bounds <- apply(Xn, 2, range)
                         if(ncol(Xn)==1) ip_bounds <- matrix(ip_bounds)
                         g_start <- ifelse(is.list(g), g$start, g)
                         theta_start <- ifelse(is.list(theta), theta$start, theta)
                         ## Builds locally optimal inducing point design with weighted IMSE
                         optXm <- optIP.wIMSE(XX[i,,drop=FALSE], M=M, Xn=Xn,
                                              theta=theta_start, g=g_start,
                                              ip_bounds=ip_bounds,
                                              integral_bounds=integral_bounds,
                                              epsK=epsK, epsQ=epsQ,
                                              verbose=FALSE)

                         ## Optimizes hyperparameters and calculates XXi prediction
                         preds <- giGP(XX[i,,drop=FALSE], Xm=optXm$Xm, Xn, Yn,
                                            theta=theta, g=g, nu=nu, epsQ=epsQ,
                                            epsK=epsK, tol=tol, reps=reps_n_list)

                         list(mean=preds$mean, var=preds$var, theta = theta,
                              nu=preds$nu, Xm=optXm$Xm, eps = preds$eps,
                              mle=preds$mle)
                       }
    stopCluster(cl)
  }

  ## Collects results from clusters
  pred.mean <- unlist(results[2:(nrow(XX)+1)])
  pred.var <- unlist(results[(nrow(XX)+2):(2*(nrow(XX)+1))])
  pred.theta <- unlist(results[(2*(nrow(XX)+1)+2):(3*(nrow(XX)+1))])
  pred.nu <- unlist(results[(3*(nrow(XX)+1)+2):(4*(nrow(XX)+1))])
  pred.Xm<- results[(4*(nrow(XX)+1)+2):(5*(nrow(XX)+1))]
  pred.eps<- matrix(unlist(results[(5*(nrow(XX)+1)+2):(6*(nrow(XX)+1))]),
                    byrow=TRUE, ncol=2)
  colnames(pred.eps) <- c('epsK', 'epsQ')
  ## For timing
  t2 <- proc.time()[3]

  if(is.list(theta) & is.list(g)){
    pred.mle <- matrix(unlist(results[(6*(nrow(XX)+1)+2):(7*(nrow(XX)+1))]),
                       ncol=3, byrow=TRUE)
  } else if(is.list(theta) | is.list(g)){
    pred.mle <- matrix(unlist(results[(6*(nrow(XX)+1)+2):(7*(nrow(XX)+1))]),
                       ncol=2, byrow=TRUE)
  } else {
    return(list(mean=pred.mean, var=pred.var, nu=pred.nu, g=g,
                theta=theta, Xm=pred.Xm, eps=pred.eps, time=t2-t1))
  }

  return(list(mean=pred.mean, var=pred.var, nu=pred.nu, g=g, theta=pred.theta,
              Xm=pred.Xm, eps=pred.eps, mle=as.data.frame(pred.mle), time=t2-t1))
}


## giGP:
##
## Generates predictions for XX based on a single (global) inducing point design
giGP <- function(XX, X = NULL, Y = NULL, Xm, g = 1e-6, theta = NULL, nu = NULL,
                 epsK = sqrt(.Machine$double.eps), epsQ = 1e-5, tol = .01,
                 reps = FALSE){

  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if(reps){
    Xorig <- X; Y_orig <- Y
    reps <- find_reps(X, Y)
    X <- reps$X0
    Y <- reps$Z0
  } else reps <- NULL


  ###### Sanity checks ######
  if(ncol(X)!=ncol(XX)) ## Dimension check
    stop('X and XX have a different number of columns')
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if(ncol(Xm)!=ncol(X)) ## Dimension check
    stop('Xm and X have a different number of columns')
  if (!is.null(nu))
    if(nu <=0)
      stop('nu should be positive')
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')
  if (tol <= 0)
    stop('tol should be a positive number')

  ## For timing
  t1 <- proc.time()[3]

  if(is.null(dim(Y))) Y <- matrix(Y, ncol=1)
  if(is.null(theta)) theta <- quantile(dist(X), .1)^2
  if(is.null(g)) g <- garg(list(mle=TRUE), Y)

  ## Optimizes hyperparameter(s) with log-likelihood if
  ## optimization bounds provided
  if(!is.list(theta)){
    if(!is.list(g)){
      out <- list(theta=theta, g=g, its=0)
    } else {
      out <- optNegLogLik_g(g, Xm, X, Y, theta, epsQ=epsQ,
                            epsK=epsK, tol=tol, rep_list=reps)
      out$theta <- theta
      mle <- list(g=out$g, its=out$its)
    }
  } else {
    if(!is.list(g)){
      out <- optNegLogLik_theta(theta, Xm, X, Y, g, epsQ=epsQ,
                                epsK=epsK, tol=tol, rep_list=reps)
      out$g <- g
      mle <- list(theta=out$theta, its=out$its)
    } else {
      out <- optNegLogLik_g_and_theta(g, theta, Xm, X, Y, epsQ=epsQ,
                                             epsK=epsK, tol=tol, rep_list=reps)
      mle <- list(g=out$g, theta=out$theta, its=out$its)
    }
  }

  ## Builds matrices for prediction
  KQ <- calcKmQm(Xm, X, theta=out$theta, g=out$g,
                 epsQ=epsQ, epsK=epsK, inv=FALSE, mults=reps$mult)

  K_mXX <- cov_gen(Xm, XX, theta=out$theta)
  QmiKXX <- try(solve(KQ$Qm, K_mXX), silent=TRUE)
  KmiKXX <- try(solve(KQ$Km, K_mXX), silent=TRUE)


  if(is.null(reps)){
    KmnOLiY <- t(KQ$OLiKnm) %*% Y
  } else {
    KmnOLiY <- t(KQ$OLiKnm) %*% (reps$mult * Y)
  }
  QmiKOLy <- try(solve(KQ$Qm, KmnOLiY), silent=TRUE)

  if (class(QmiKXX)[1]=="try-error" | class(KmiKXX)[1]=="try-error" |
      class(QmiKOLy)[1]=="try-error" ){
    if (exists('out$epsK')){
      eps_new <- list(K=out$epsK, Q=out$epsK)
      KQ <- calcKmQm(Xm, X, out$theta, g=out$g, epsK=out$epsK, epsQ=out$epsQ,
                     inv=FALSE, mults=reps$mult)
      K_mXX <- cov_gen(Xm, XX, theta=out$theta)
      QmiKXX <- try(solve(KQ$Qm + diag(out$epsQ, nrow(Xm)), K_mXX), silent=TRUE)
      KmiKXX <- try(solve(KQ$Km, K_mXX), silent=TRUE)

      KmnOLiY <- ifelse(is.null(reps), t(KQ$OLiKnm) %*% Y,
                        t(KQ$OLiKnm) %*% (reps$mult * Y))
      QmiKOLy <- try(solve(KQ$Qm + diag(out$epsQ, nrow(Xm)), KmnOLiY),
                     silent=TRUE)
    } else {
      increase_epsK <- increase_epsQ <- 0
      matrix_inverse_calcs <- function(epsK, epsQ){
        KQ <- calcKmQm(Xm, X, out$theta, g=out$g, epsK=epsK, epsQ=epsQ,
                       inv=FALSE, mults=reps$mult)
        K_mXX <- cov_gen(Xm, XX, theta=out$theta)
        QmiKXX <<- try(solve(KQ$Qm + diag(epsQ, nrow(Xm)), K_mXX), silent=TRUE)
        KmiKXX <<- try(solve(KQ$Km, K_mXX), silent=TRUE)

        if(is.null(reps)){
          KmnOLiY <<- t(KQ$OLiKnm) %*% Y
        } else {
          KmnOLiY <<- t(KQ$OLiKnm) %*% (reps$mult * Y)
        }
        QmiKOLy <<- try(solve(KQ$Qm + diag(epsQ, nrow(Xm)), KmnOLiY), silent=TRUE)
      }

      while ((class(QmiKXX)[1]=="try-error" | class(KmiKXX)[1]=="try-error" |
              class(QmiKOLy)[1]=="try-error")
             & ((epsK*(10^increase_epsK) < 1e-3) |
                (epsQ*(10^increase_epsQ) < 1e-3))) {
        if (epsQ*(10^increase_epsQ) < 1e-3){
          increase_epsQ <- increase_epsQ + 1
          try(matrix_inverse_calcs(epsK=epsK, epsQ=epsQ*(10^increase_epsQ)),
              silent=TRUE)
        } else {
          increase_epsK <- increase_epsK + 1
          try(matrix_inverse_calcs(epsK=epsK*(10^increase_epsK), epsQ=epsQ),
              silent=TRUE)
        }

      }
      if (increase_epsK > 1){
        eps_new <- list(K=epsK*(10^increase_epsK), Q=epsQ)
      } else {
        eps_new <- list(K=epsK, Q=epsQ*(10^increase_epsQ))
      }

    }
    if (class(QmiKXX)[1]=="try-error" | class(KmiKXX)[1]=="try-error" |
        class(QmiKOLy)[1]=="try-error" ){
      stop('Covariance matrix approximation is unstable.',
           ' Consider decreasing n or inducing point design')
    }

  }
  if (!exists('eps_new')) eps_new <- list(K=epsK, Q=epsQ)

  ## Calculates predictive mean and variance for each XX[i,]
  var_vec <- 1 + out$g - colSums(K_mXX * (KmiKXX - QmiKXX))

  if (is.null(reps)){
    mean_vec <- t(QmiKXX) %*% t(KQ$OLiKnm) %*% Y
    if(is.null(nu))
      nu <- (sum(Y^2 / KQ$OL) - t(KmnOLiY) %*% QmiKOLy)/length(Y)
  } else {
    mean_vec <- t(QmiKXX) %*% t(KQ$OLiKnm) %*% (reps$mult * Y)
    if(is.null(nu))
      nu <- (sum(reps$Z^2 / rep(KQ$OL, reps$mult)) -
                      t(KmnOLiY) %*% QmiKOLy)/length(reps$Z)
  }


  ## For timing
  t2 <- proc.time()[3]

  if(!is.list(g) & !is.list(theta)){
    return(list(mean=mean_vec, var=var_vec, nu=drop(nu), g=g, theta=theta,
                eps=eps_new, time=t2 - t1))
  } else {
    return(list(mean=mean_vec, var=var_vec, nu=drop(nu), g=g, theta=theta,
                mle=mle, eps=eps_new, time=t2 - t1))
  }

}


## wgigp:
##
## Returns prediction multiplied by a Gaussian weight
wgiGP <- function(xstar_i, xstar, gauss_sd, Xn, Yn, Xm, theta, g,
                  epsK = 1e-6, epsQ = 1e-5){

  nonzero_dim <- which(gauss_sd != 0)
  XX <- t(replicate(length(xstar_i), as.vector(xstar)))
  XX[,nonzero_dim] <- xstar_i

  gigp_model <- giGP(XX, Xn, Yn, Xm, theta=theta, g=g, epsK=epsK, epsQ=epsQ)
  pred_xstar <- gigp_model$mean

  weight <- dnorm(xstar_i, xstar[,nonzero_dim], gauss_sd[nonzero_dim])
  return(weight*pred_xstar)
}

## ligp_gauss_measure:
##
## Estimates the average response over a Gaussian measure centered at xstar.
## Uses either a vector of Gaussian noise to generate predictive locations
## whose predictions are averaged or calls "integrate" to use quadrature to
## estimate.
liGP_gauss_measure <- function(xstar, X, Y, Xm.t, N, gauss_sd,
                               measure_bounds = c(-Inf, Inf),
                               g = 1e-6, epsi = NULL,
                               epsK = 1e-6, epsQ = 1e-5, seq_length = 20){
  if(is.null(dim(xstar))) xstar <- matrix(xstar, nrow=1)
  nonzero_dim <- which(gauss_sd != 0)

  # Construct reference set for Gaussian measure
  ndim <- ncol(X)
  dfs <- list()
  for (i in 1:ndim){
    if (i == nonzero_dim) {
      dfs[[i]] <- seq(xstar[,i]-2*gauss_sd[i], xstar[,i]+2*gauss_sd[i],
                      length=seq_length)
    } else{
      dfs[[i]] <- xstar[,i]
    }
  }
  xstar_measure <- as.matrix(expand.grid(dfs[1:ndim]))

  # Build xstar neighborhood
  xx_dists <- distance(xstar_measure, X)
  min_dists <- apply(xx_dists, 2, min)
  quant <- quantile(min_dists, N/nrow(X))
  closest_indices <- min_dists < quant
  Xn <- X[closest_indices,]
  Yn <- Y[closest_indices]

  neighborhood_box <- apply(Xn, 2, range)
  Xn_theta <- darg(NULL, Xn)
  Xn_theta$max <- 20*Xn_theta$max
  Xm <- sweep(Xm.t, 2, as.vector(xstar), '+')

  if (is.null(epsi)){
    theta_out <- optNegLogLik_theta(Xn_theta, Xm, Xn, Yn, g,
                                           epsK=epsK, epsQ=epsQ)

    measure_est <- integrate(wgiGP, lower=measure_bounds[1],
                             upper=measure_bounds[2],
                             xstar=xstar, gauss_sd=gauss_sd,
                             Xn=Xn, Yn=Yn, Xm=Xm,
                             theta=theta_out$theta, g=g,
                             epsK=epsK, epsQ=epsQ)
    estimate <- measure_est$value
  } else {
    XX <- t(replicate(length(epsi), as.vector(xstar)))
    XX[,nonzero_dim] <- XX[,nonzero_dim] + epsi

    preds <- giGP(XX, Xn, Yn, Xm, theta=Xn_theta, g=g, epsK=epsK, epsQ=epsQ)

    estimate <- mean(preds$mean)
  }

  return(estimate)
}

## build_neighborhood:
##
## Builds local neighborhood based on closest nearest neighbor design points.
## If a reps list is provided, the neighborhood is built with the unique
## design locations.
build_neighborhood <- function(N, xx = NULL, X = NULL, Y = NULL,
                               reps_list = NULL){
  if (is.null(X)){
      if (is.null(reps_list$X0)){
        stop('X or reps_list$X0 must be supplied')
      } else {
      if(is.null(xx)) xx <- matrix(apply(reps_list$X0, 2, median), nrow=1)

      xx_to_X_dists <- distance(xx, reps_list$X0)
      quant_dists <- quantile(xx_to_X_dists, N/nrow(reps_list$X0))
      Xn <- reps_list$X0[xx_to_X_dists < quant_dists,,drop=FALSE]
      neighborhood <- list(xx=xx, Xn0=Xn)
      if (!is.null(reps_list$Z0))
        neighborhood$Yn0 <- reps_list$Z0[xx_to_X_dists < quant_dists]
      if (!is.null(reps_list$mult))
        neighborhood$mult <- reps_list$mult[xx_to_X_dists < quant_dists]
      if(!is.null(reps_list$Zlist)){
        neighborhood$Yn_list <- list()
        selected_rows <- which(xx_to_X_dists < quant_dists)
        for (i in 1:N)
          neighborhood$Yn_list[[i]] <- reps_list$Zlist[[selected_rows[i]]]
      }
    }
  }
  else {
    if (is.null(xx)) xx <- matrix(apply(X, 2, median), nrow=1)
    xx_to_X_dists <- distance(xx, X)
    quant_dists <- quantile(xx_to_X_dists, N/nrow(X))
    Xn <- X[xx_to_X_dists < quant_dists,,drop=FALSE]
    neighborhood <- list(xx=xx, Xn=Xn)
    if (!is.null(Y)) neighborhood$Yn <- Y[xx_to_X_dists < quant_dists]
  }

  return(neighborhood)
}
