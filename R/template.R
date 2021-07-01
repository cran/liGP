eps <- sqrt(.Machine$double.eps)

## build_ipTemplate:
##
## Creates an inducing points design optimized with ALC or wIMSE,
## then is returned centered at the origin
build_ipTemplate <- function(X = NULL, Y = NULL, M, N, theta = NULL, g = 1e-4,
                             method = c('wimse','alc'),
                             ip_bounds = NULL, integral_bounds = NULL,
                             num_thread = 1, num_multistart = 20, w_var = NULL,
                             epsK = sqrt(.Machine$double.eps), epsQ = 1e-5,
                             reps = FALSE, verbose = TRUE){
  ## Collects data from X,Y or reps_list
  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    reps_list <- reps
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if (reps){
    Xorig <- X; Y_orig <- Y
    reps_list <- find_reps(X, Y)
    X <- reps_list$X0
    Y <- reps_list$Z0
  } else reps_list <- FALSE

  ###### Sanity checks ######
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if (M > N)
    warning('Number of inducing points (M) > Neighborhood size (N)')
  if(N > nrow(X)) ## Neighborhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  if(!is.null(theta))
    if (theta <= 0)
      stop('theta should be a positive number')
  if (g < 0)
    stop('g must be positive')
  if(!method %in% c('wimse', 'alc'))
    stop('A valid method was not given. Choices include: wimse and alc')
  if(method=='wimse'){
    Xrange <- apply(X, 2, range)
    if (!is.null(integral_bounds))
      if(sum(Xrange[1,] < integral_bounds[1,]) > 0 ||
         sum(Xrange[2,] > integral_bounds[2,]) > 0)
        stop('X outside integration bounds')
    if(is.null(integral_bounds))
      integral_bounds <- Xrange
  }
  if(!is.null(ip_bounds)){
    if(nrow(ip_bounds) !=2)
      stop('ip_bounds should be a matrix with two rows')
    if(sum(ip_bounds[1,] < ip_bounds[2,]) < ncol(X))
      stop('At least one dimensions bounds in ip_bounds is incorrect')
  }
  if(num_multistart < 0 | num_multistart %% 1 !=0)
    stop('num_multistart is not a positive integer')
  ## Checks that number of threads is valid
  if(num_thread %% 1 != 0 | num_thread < 1)
    stop('num_thread is not a positive integer')
  if (num_thread > num_multistart){
    warning(paste("num_thread > num_multistart. num_thread set to",num_multistart))
    num_thread <- num_multistart
  }
  available_cores <- detectCores()
  if (num_thread > available_cores){
    warning(paste("num_thread exceeds number of available cores.",
                  "num_thread set to", available_cores))
    num_thread <- available_cores
  }
  if (epsK <= 0)
    stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')


  ## For timing
  t1 <- proc.time()[3]

  ## Builds neighborhood at center of design
  Xc <- matrix(apply(X, 2, median), nrow=1)

  if(is.list(reps_list)){
    reps_n_list <- build_neighborhood(N, Xc, reps_list = reps_list)
    Xn <- reps_n_list$Xn
    Yn <- reps_n_list$Yn
  } else {
    rep_n_list <- NULL
    neighborhood <- build_neighborhood(N, Xc, X, Y)
    Xn <- neighborhood$Xn; Yn <- neighborhood$Yn
  }

  neighborhood_box <- apply(Xn, 2, range)

  if(is.null(ip_bounds))
    ip_bounds <- neighborhood_box
  if (is.null(theta)) theta <- quantile(dist(Xn),.1)^2


  if (method == 'alc'){
    low.bound <- neighborhood_box[1,]; upp.bound <- neighborhood_box[2,]
    ## Builds inducing point design by optimizing ALC
    p <- optIP.ALC(Xc=Xc, Xref=NULL, M=M, Xn=Xn,
                   Yn=Yn, theta=theta, g=g,
                   ip_bounds=ip_bounds, num_thread=num_thread,
                   num_multistart=num_multistart,
                   verbose=verbose, epsQ=epsQ, epsK=epsK,
                   rep_list=rep_n_list)
    Xm.t <- sweep(p$Xm, 2, Xc)

  } else {
    ## Builds inducing point design by optimizing weighted IMSE
    p <- optIP.wIMSE(Xn=Xn, M=M, theta=theta,
                     g=g, w_mean=Xc, ip_bounds=ip_bounds,
                     integral_bounds=integral_bounds, w_var=w_var,
                     num_multistart=num_multistart, verbose=verbose,
                     epsQ=epsQ, epsK=epsK, mult=rep_n_list$mult)
    Xm.t <- sweep(p$Xm, 2, Xc)
  }
  ## For timing
  t2 <- proc.time()[3]

  return(list(Xm.t=Xm.t, Xn=Xn, Xc=Xc, time=t2-t1))
}


## scale_ipTemplate:
##
## Scales a inducing points design in [0,1]^d to fill local neighborhood.
## Returns template design centered at the origin.
scale_ipTemplate <- function(X, N, space_fill_design,
                             method = c('qnorm', 'chr')){
  t1 <- proc.time()[3]
  ###### Sanity checks ######
  if(N > nrow(X)) ## Neighorhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  if (ncol(space_fill_design) != ncol(X))
    stop('A space filling design was supplied with an incorrect ',
         'number of columns.')
  if (nrow(space_fill_design) > N)
    warning('Size of space_filling_design > Neighborhood size (N)')
  if (!method %in% c('qnorm','chr'))
    stop('A valid method was not given. Choices include: qnorm, chr')

  Xc <- matrix(apply(X, 2, median), nrow=1)
  Xn <- build_neighborhood(N, Xc, X)$Xn
  neighborhood_box <- apply(Xn, 2, range)

  if(method == 'qnorm'){
    ## Scales design by inverse normal CDF
    dist_from_Xc <- sweep(neighborhood_box, 2, Xc)
    qnorm_sd <- apply(abs(dist_from_Xc), 2, max)/3

    Xm.qnorm <- qnormscale(space_fill_design, rep(0, ncol(X)), qnorm_sd)
    Xm.t <- rbind(rep(0,ncol(Xn)), Xm.qnorm)
  } else {
    ## Scales design to a circumscribed hyperrectangle
    Xm.t <- sweep(space_fill_design - .5, 2,
                  neighborhood_box[2,] - neighborhood_box[1,], '*')
    Xm.t <- rbind(rep(0, ncol(Xn)), Xm.t)
  }

  ## For timing
  t2 <- proc.time()[3]

  return(list(Xm.t=Xm.t, Xn=Xn, time=t2-t1))
}


## qnormscale:
##
## The function that scales X to center around mean
## with standard deviations determined by sd, which can be vectorized
qnormscale <- function(X, mean, sd){
  m <- ncol(X)
  ## Sanity checks
  if(length(mean) == 1) {mean <- rep(mean, m)
  } else if(length(mean) != m) stop("X and mean dimension mismatch")
  if(length(sd) == 1) {sd <- rep(sd, m)
  } else if(length(sd) != m) stop("X and sd dimension mismatch")

  ## Scale each dimension independently
  for(j in 1:ncol(X))
    X[,j] <- qnorm(X[,j], mean=mean[j], sd=sd[j])

  ## Return scaled matrix
  return(X)
}

## build_gauss_measure_ipTemplate:
##
## Creates an inducing points design based on a local neighborhood
## for a Gaussian measure slice. Inducing points are optimized with
## wIMSE, and then returned centered at the origin
build_gauss_measure_ipTemplate <- function(X = NULL, Y = NULL, M, N, gauss_sd,
                                           theta = NULL, g = 1e-4,
                                           seq_length = 20, ip_bounds = NULL,
                                           integral_bounds = NULL,
                                           num_multistart = 20,
                                           epsK = sqrt(.Machine$double.eps),
                                           epsQ = 1e-5, reps = FALSE,
                                           verbose = TRUE){
  if(is.list(reps)){
    if(is.null(reps$X0)) stop('reps doesn\'t include \'X0\' in list')
    if(is.null(reps$Z0)) stop('reps doesn\'t include \'Z0\' in list')
    if(is.null(reps$mult)) stop('reps doesn\'t include \'mult\' in list')
    if(is.null(reps$Z)) stop('reps doesn\'t include \'Z\' in list')
    if(is.null(reps$Zlist)) stop('reps doesn\'t include \'Zlist\' in list')
    reps_list <- reps
    X <- reps$X0
    Y <- reps$Z0
  } else if (is.null(X) | is.null(Y)){
    stop('X and Y are required')
  } else if (reps){
    Xorig <- X; Y_orig <- Y
    reps_list <- find_reps(X, Y)
    X <- reps_list$X0
    Y <- reps_list$Z0
  } else reps_list <- FALSE

  ###### Sanity checks ######
  if(nrow(X)!=length(Y)) ## Number of entries check
    stop('Number of entries in Y doesn\'t match nrows of X')
  if (M > N)
    warning('Number of inducing points (M) > Neighborhood size (N)')
  if(N > nrow(X)) ## Neighorhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  nonzero_dim <- which(gauss_sd!=0)
  if(length(nonzero_dim) > 1)
    stop('The Gaussian measure can only have a non-zero gauss_sd ',
         'in one dimension.')
  if(length(gauss_sd) != ncol(X))
    stop('The number of entries and gauss_sd and ncol(X) do not match')
  if(!is.null(theta))
    if (theta <= 0)
      stop('theta should be a positive number')
  if (g < 0)
    stop('g must be positive')

  Xrange <- apply(X, 2, range)
  if (!is.null(integral_bounds))
    if(sum(Xrange[1,] < integral_bounds[1,]) > 0 ||
       sum(Xrange[2,] > integral_bounds[2,]) > 0)
      stop('X outside integration bounds')
  if(is.null(integral_bounds))
    integral_bounds <- Xrange
  if(!is.null(ip_bounds)){
    if(nrow(ip_bounds) !=2)
      stop('ip_bounds should be a matrix with two rows')
    if(sum(ip_bounds[1,] < ip_bounds[2,]) < ncol(X))
      stop('At least one dimensions bounds in ip_bounds is incorrect')
  }
  if(num_multistart < 0 | num_multistart %% 1 !=0)
    if (epsK <= 0)
      stop('epsK should be a positive number')
  if (epsQ <= 0)
    stop('epsQ should be a positive number')

  ## For timing
  t1 <- proc.time()[3]

  ##-----------------------------------------
  ## Builds neighborhood at center of design
  Xc <- matrix(apply(X, 2, median), nrow=1)

  # Construct reference set for Gaussian measure
  ndim <- ncol(X)
  dfs <- list()
  for (i in 1:ndim){
    if (i == nonzero_dim) {
      dfs[[i]] <- seq(Xc[,i] - 2*gauss_sd[i], Xc[,i] + 2*gauss_sd[i],
                      length=seq_length)
    } else{
      dfs[[i]] <- Xc[,i]
    }
  }
  Xc_measure <- as.matrix(expand.grid(dfs[1:ndim]))

  # Build Xc neighborhood
  if(N == nrow(X)){
    Xn <- X
  } else{
    xx_dists <- distance(Xc_measure, X)
    min_dists <- apply(xx_dists, 2, min)
    quant <- quantile(min_dists, N/nrow(X))
    closest_indices <- min_dists < quant
    Xn <- X[closest_indices,]
  }

  neighborhood_box <- apply(Xn, 2, range)
  Xnc_theta <- darg(NULL, Xn)$start

  ## Change gauss_sd to allow some weight in dimensions where it's zero
  nonzero2zero.ratio <- (neighborhood_box[2, -nonzero_dim] -
                           neighborhood_box[1, -nonzero_dim])/
    (neighborhood_box[2, nonzero_dim] - neighborhood_box[1, nonzero_dim])
  gauss_sd[-nonzero_dim] <- nonzero2zero.ratio*gauss_sd[nonzero_dim]

  if(is.list(reps_list)){
    rep_n_list <- list(mult=reps_list$mult[closest_indices],
                       Z=matrix(c(unlist(reps_list$Zlist[closest_indices]))))
  } else rep_n_list <- NULL
  if(is.null(ip_bounds))
    ip_bounds <- neighborhood_box
  ##------------------------------------------------------------
  ## Builds inducing point design by optimizing weighted IMSE
  Xm.wimse <- try(optIP.wIMSE(Xn=Xn, M=M, theta=Xnc_theta, g=g,
                         w_mean=Xc, w_var=gauss_sd^2,
                         ip_bounds=ip_bounds,
                         integral_bounds=integral_bounds,
                         num_multistart=num_multistart, verbose=verbose,
                         epsQ=epsQ, epsK=epsK, mult=rep_n_list$mult)$Xm,
             silent=TRUE)
  increase_epsK <- increase_epsQ <- 1
  while (class(Xm.wimse)=='try-error' & (epsK < 1e-3 & epsQ < 1e-3)) {
    if (epsQ < 1e-3){
      Xm.wimse <- try(optIP.wIMSE(Xn=Xn, M=M, theta=Xnc_theta, g=g,
                             w_mean=Xc, w_var=gauss_sd^2,
                             ip_bounds=ip_bounds,
                             integral_bounds=integral_bounds,
                             num_multistart=num_multistart, verbose=verbose,
                             epsQ=epsQ*(10^increase_epsQ),
                             epsK=epsK, mult=rep_n_list$mult)$Xm, silent=TRUE)
      increase_epsQ <- increase_epsQ + 1
    } else {
      increase_epsQ <- 1
      Xm.wimse <- try(optIP.wIMSE(Xn=Xn, M=M, theta=Xnc_theta, g=g,
                             w_mean=Xc, w_var=gauss_sd^2,
                             ip_bounds=ip_bounds,
                             integral_bounds=integral_bounds,
                             num_multistart=num_multistart, verbose=verbose,
                             epsQ=epsQ, epsK=epsK*(10^increase_epsK),
                             mult=rep_n_list$mult)$Xm, silent=TRUE)
      increase_epsK <- increase_epsK + 1
    }

  }

  Xm.t <- sweep(Xm.wimse, 2, Xc)

  ## For timing
  t2 <- proc.time()[3]

  return(list(Xm.t=Xm.t, Xn=Xn, Xc=Xc, gauss_sd=gauss_sd, time=t2-t1))
}

## scale_gauss_measure_ipTemplate:
##
## Scales a inducing points design in [0,1]^d to fill local neighborhood.
## Returns template design centered at the origin.
scale_gauss_measure_ipTemplate <- function(X, N, gauss_sd,
                                           space_fill_design,
                                           method = c('qnorm','chr'),
                                           seq_length=20){
  t1 <- proc.time()[3]
  ###### Sanity checks ######
  if(N > nrow(X)) ## Neighorhood cannot be bigger than data size
    stop('N is greater than the number of rows in X')
  nonzero_dim <- which(gauss_sd!=0)
  if(length(nonzero_dim) > 1)
    stop('The Gaussian measure can only have a non-zero gauss_sd ',
         'in one dimension.')
  if(length(gauss_sd) != ncol(X))
    stop('The number of entries and gauss_sd and ncol(X) do not match')
  if (ncol(space_fill_design) != ncol(X))
    stop('A space filling design was supplied with an incorrect ',
         'number of columns.')
  if (nrow(space_fill_design) > N)
    warning('Size of space_filling_design > Neighborhood size (N)')
  if (!method %in% c('qnorm','chr'))
    stop('A valid method was not given. Choices include: qnorm, chr')

  ##-----------------------------------------
  ## Builds neighborhood at center of design
  Xc <- matrix(apply(X, 2, median), nrow=1)

  # Construct reference set for Gaussian measure
  ndim <- ncol(X)
  dfs <- list()
  for (i in 1:ndim){
    if (i == nonzero_dim) {
      dfs[[i]] <- seq(Xc[,i] - 2*gauss_sd[i], Xc[,i] + 2*gauss_sd[i],
                      length=seq_length)
    } else{
      dfs[[i]] <- Xc[,i]
    }
  }
  Xc_measure <- as.matrix(expand.grid(dfs[1:ndim]))

  # Build Xc neighborhood
  if(N == nrow(X)){
    Xn <- X
  } else{
    xx_dists <- distance(Xc_measure, X)
    min_dists <- apply(xx_dists, 2, min)
    quant <- quantile(min_dists, N/nrow(X))
    closest_indices <- min_dists < quant
    Xn <- X[closest_indices,]
  }
  neighborhood_box <- apply(Xn, 2, range)
  Xnc_theta <- darg(NULL, Xn)$start

  ## Change gauss_sd to allow some weight in dimensions where it's zero
  nonzero2zero.ratio <- (neighborhood_box[2,-nonzero_dim] -
                           neighborhood_box[1,-nonzero_dim])/
    (neighborhood_box[2,nonzero_dim] - neighborhood_box[1,nonzero_dim])
  gauss_sd[-nonzero_dim] <- nonzero2zero.ratio*gauss_sd[nonzero_dim]


  if(method == 'qnorm'){
    ## Scales design by inverse normal CDF
    dist_from_Xc <- sweep(neighborhood_box, 2, Xc)
    qnorm_sd <- apply(abs(dist_from_Xc), 2, max)/3

    Xm.qnorm <- qnormscale(space_fill_design, rep(0, ncol(X)), qnorm_sd)
    Xm.t <- rbind(rep(0,ncol(Xn)), Xm.qnorm)
  } else {
    ## Scales design to a circumscribed hyperrectangle
    Xm.t <- sweep(space_fill_design - .5, 2,
                  neighborhood_box[2,] - neighborhood_box[1,], '*')
    Xm.t <- rbind(rep(0, ncol(Xn)), Xm.t)
  }

  ## For timing
  t2 <- proc.time()[3]

  return(list(Xm.t=Xm.t, Xn=Xn, Xc=Xc, gauss_sd=gauss_sd, time=t2 - t1))
}

