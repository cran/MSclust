
#####################
## Update of Gamma ##
#####################

Gamma.update <- function(X, z = NULL, w = NULL, mu, Lambda, Gamma0, method = "BFGS", iter.max=5000, reltol=1e-15, trace=0){
  
  # w is a nxd matrix of weights
  
  if(is.vector(X))
    X <- matrix(X, ncol = 1)
  if(is.data.frame(X))
    X <- data.matrix(X)
  if(!is.matrix(X))
    X <- matrix(X, nrow = length(X), ncol = 1, dimnames = list(names(X), deparse(substitute(X))))
  if(!is.numeric(X))
    stop("numeric matrix/vector expected for X")
  if(any(is.na(X)))
    stop("No NAs allowed")
  if(is.null(Gamma0))
    stop("An initial value for (the eigenvectors matrix) Gamma0 is needed")
  
  n <- nrow(X)
  d <- ncol(X)
  
  if(is.null(w)){
    w <- array(1,c(n,d),dimnames=list(1:n,paste("dim.",1:d,sep="")))
  }
  if(is.null(z)){
    z <- array(1,c(n),dimnames=list(1:n))
  }
  
  # -------------- #
  # Initialization #
  # -------------- #
  
  # Gamma
  
  # starting-world
  
  # if(is.null(Gamma0)){
  #   #S <- crossprod(sqrt(weights/sum(weights))*(X-matrix(rep(mu,n),n,d,byrow=TRUE)))
  #   S     <- cov(X)
  #   temp  <- eigen(S)
  #   Gamma <- temp$vectors
  # }else Gamma <- Gamma0
  
  # working-world
  
  tempfac <- QPl(Gamma0)
  P       <- tempfac$P
  initial.values <- tempfac$l
  
  # objective function
  
  f <- function(par, P, X, z, w, mu, Lambda, n){
    
    # ----- #
    # Gamma #
    # ----- #
    
    Gamma <- PlQ(P = P, l = par)
    
    # ---------------------------------- #
    # Objective function to be minimized #
    # ---------------------------------- #
    
    f.obj <- 0
    for(i in 1:n){
      f.obj <- f.obj + psych::tr( Gamma %*% diag( w[i,]/diag(Lambda) ) %*% t(Gamma) %*% base::tcrossprod( sqrt(z[i])*(X[i,] - mu) ) )
    }
    
    return(f.obj)
    
  }
  
  res <- optim(par = initial.values, fn = f, P = P, X = X, z = z, w = w, mu = mu, Lambda = Lambda, n = n, method = method, control = list(maxit = iter.max, reltol = reltol, trace = trace))
  
  f.obj <- res$value
  est   <- res$par
  
  Gamma <- PlQ(P = P, l = est)
  
  return(Gamma)
  
}

###########################################
## Generate random pxp covariance matrix ##
###########################################

covariance.matrix <- function(d){

  if(d>1) {
    sigma <- matrix(rnorm(d*d), ncol=d)
    sigma <- crossprod(sigma)+ diag(rep(0.1, d))
  }
  else sigma <- matrix(rchisq(1,1), ncol=d)

  return(sigma)

}

######################################################################################################################################################################################
################################################# MSCN mixture #######################################################################################################################
######################################################################################################################################################################################

#######################
## Density of a MSCN ##
#######################

dmscn <- function(x, mu = NULL, L=NULL, G=NULL, Sigma=NULL, alpha = NULL, eta = NULL){
  Lambda=L
  Gamma=G
  d = ifelse(is.null(x=ncol(x)),1,ncol(x))
  # Lambda is a matrix with diagonal elements
  if(is.null(mu)) mu=rep(0,d)
  if(is.null(alpha)) alpha = rep(0.99,d)
  if(is.null(eta)) eta = rep(1.01,d)
  
  if(min(alpha)<0 | max(alpha)>1)
    stop("alpha must be a vector with elements in (0,1)")
  if(min(eta)<=1)
    stop("eta must be a vector with elements greater than 1")

  if(is.vector(mu))
    d <- length(mu)
  else
    d <- 1

  if(!is.null(Sigma)){
    temp   <- eigen(Sigma)
    Lambda <- diag(temp$values)
    Gamma  <- temp$vectors
  }
  if(is.null(Sigma)&is.null(G)){Sigma=diag(d)
  temp   <- eigen(Sigma)
  Lambda <- diag(temp$values)
  Gamma  <- temp$vectors
  }
  
  if(is.matrix(x))
    q <- nrow(x)
  if(is.vector(x) & d > 1)
    q <- 1
  if(is.vector(x) & d == 1)
    q <- length(x)

  dens <- array(0,c(q,d),dimnames=list(1:q,paste("X.",1:d,sep="")))

  for(i in 1:q){
    if(is.matrix(x))
      xnew <- t(Gamma) %*% (x[i,]-mu)
    if(is.vector(x) & d > 1)
      xnew <- t(Gamma) %*% (x-mu)
    if(is.vector(x) & d == 1)
      xnew <- x[i] - mu

    for(j in 1:d){
      dens[i,j] <- alpha[j]*dnorm(xnew[j],mean=0,sqrt(Lambda[j,j])) + (1-alpha[j])*dnorm(xnew[j],mean=0,sqrt(eta[j]*Lambda[j,j]))
    }
  }

  PDF <- apply(dens,1,prod)
  PDF <- (PDF<10^(-323))*10^(-323)+(PDF>10^(-323))*PDF

  return(PDF)

}

###################################
## Random Generation from a MSCN ##
###################################

rmscn <- function(n, d = 2, mu = rep(0,d), L=NULL, G=NULL, Sigma=diag(d), alpha = rep(0.99,d), eta = rep(1.01,d)){
  Lambda=L
  Gamma=G
  # Lambda is a vector with diagonal elements

  if(min(alpha)<0 | max(alpha)>1)
    stop("alpha must be a vector with elements in (0,1)")
  if(min(eta)<=1)
    stop("eta must be a vector with elements greater than 1")

  if(is.vector(mu))
    d <- length(mu)
  else
    d <- 1

  if(!is.null(Sigma)){
    temp   <- eigen(Sigma)
    Lambda <- diag(temp$values)
    Gamma  <- temp$vectors
  }

  X <- good <- array(0,c(n,d),dimnames=list(1:n,paste("X.",1:d,sep="")))
  for(j in 1:d)
    good[,j] <- rbinom(n=n,size=1,prob=alpha[j])
  Delta <- SigmaCont <- array(0,c(n,d,d),dimnames=list(1:n,paste("X.",1:d,sep=""),paste("X.",1:d,sep="")))
  for(i in 1:n){
    Delta[i,,]     <- diag((good[i,]+(1-good[i,])/eta)^(-1))
    SigmaCont[i,,] <- Gamma %*% Delta[i,,] %*% Lambda %*% t(Gamma)
    X[i,]          <- mnormt::rmnorm(n = 1, mean = mu, varcov = SigmaCont[i,,])
  }

  return(X)

}

##########################
## ML for MSCN mixtures ##
##########################

mscn <- function(
  X,                 # data matrix
  k,                 # number of mixture components
  ini = "km", # or manual; in that case, a partition (init.clus) must be provided
  sz = NULL, #IF INITIALIZATION IS MANNUAL, this matrix contains the starting value for z
  al= c(0.5,0.99),   # minimum and maximum proportion of good points in each group for the contaminated normal distribution
  eta.min = 1.01, # minimum value for inflation value for the variance covariance matrix for the "bad" points
  m = "BFGS", #  method for the optimization of the eigenvector matrix, see optim for other options
  stop = c(10^-5,200),     # the 2 stopping rules tollerance and maximum number of iterations
  VB = FALSE       # If TRUE, tracing information on the progress of the optimization is produced; see optim() for details
  #and plot of the log-likelihood versus iterations
)
{md=m
  start.z=sz
  initialization=ini
  method =md
  alpha.min =al[1]
  alpha.max = al[2]
  tol=stop[1]
  iter.max=stop[2]
  trace = 0
  plotll = FALSE
  if (VB==TRUE){
    trace = 1
    plotll = TRUE
    }
    
  
  
  
    if(is.data.frame(X))
    X <- data.matrix(X)
  if(!is.matrix(X))
    X <- matrix(X, nrow = length(X), ncol = 1, dimnames = list(names(X), deparse(substitute(X))))
  if(!is.numeric(X))
    stop("numeric matrix/vector expected for X")
  if(any(is.na(X)))
    stop("No NAs allowed")
  if(is.null(k))
    stop("The number of groups k is NULL")

  n <- nrow(X) # number of units
  d <- ncol(X) # number of variables
  
  Xnew <- array(0,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  
  # --------------------- #
  # Definition of objects #
  # --------------------- #

  # E-step quantities

  z          <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  v          <- array(1/n,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  zvgood     <- array(1,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  zvbad      <- array(0,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  correction <- array(1,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))

  # M-step Parameters

  prior  <- numeric(k)
  mu     <- array(0,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  mu.mat <- array(0,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Sigma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Lambda <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Gamma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  alpha  <- array(alpha.max,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  eta    <- array(eta.min,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))

  # n/(n+1)
  # (n+1)/n
  
  # Distribution

  densX <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # fX
  dens  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # weights*fX

  # ------------------------- #
  # posteriors initialization #
  # ------------------------- #

  if(initialization=="random.soft"){

    z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k)
    z  <- z/rowSums(z)             # soft posterior probabilities (n x k)

  }

  if(initialization=="random.hard"){

    z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)

  }

  if(initialization=="manual"){ # z.start can be both soft and hard initialization

    z  <- start.z # soft or hard classification

  }

  if(initialization=="mclust"){

    if(d==1)
      mclustfit <- mclust::Mclust(data=X, G=k, modelNames="V")
    if(d>1)
      mclustfit <- mclust::Mclust(data=X, G=k, modelNames="VVV")

    z <- mclustfit$z

  }

  if(initialization=="km"){

    if(k>1){
      km       <- kmeans(X,k)
      groups   <- km$cluster
      z        <- mclust::unmap(classification=groups)
    }
    if(k==1){
      groups   <- rep(1,n)
      z        <- matrix(groups,nrow=n,ncol=k)
    }

  }

  if(initialization=="pam"){

    lk     <- cluster::pam(x=X,k=k)
    groups <- lk$clustering
    z      <- mclust::unmap(classification=groups)

  }

  for(j in 1:k){
    for(h in 1:d){
      zvgood[,h,j] <- z[,j]*v[,h,j]
      zvbad[,h,j]  <- z[,j]*(1-v[,h,j])
      correction[,h,j] <- v[,h,j] + (1-v[,h,j])/eta[h,j]
    }
  }

  # initialization of Gamma

  for(j in 1:k){
    xbar <- colSums(z[,j]/sum(z[,j])*X)
    S    <- crossprod(sqrt(z[,j]/sum(z[,j]))*(X-matrix(rep(xbar,n),n,d,byrow=TRUE)))
    #S <- as.matrix(forceSymmetric(x=S))
    temp  <- eigen(S)
    Gamma[,,j] <- temp$vectors
  }

  # EM algorithm --------------------------------------------------

  # Preliminary definition of convergence criterions

  check     <- 0
  iteration <- 1
  loglik    <- -Inf

  while(check<1){

    # ++++++ #
    # M-Step #
    # ++++++ #

    # Parameters

    # pi

    prior <- colMeans(z)

    # alpha

    size <- colSums(z)
    for(j in 1:k){
      for(h in 1:d){
        alpha[h,j] <- max(alpha.min,sum(zvgood[,h,j])/size[j])
      }
    }

    # mu
    
    mu     <- array(0,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
    for(j in 1:k){
      den <- colSums(z[,j]*correction[,,j])
      for(i in 1:n){
        Delta.norm <- diag(z[i,j]*correction[i,,j]/den)
        mu[,j] <- mu[,j] + Gamma[,,j] %*% Delta.norm %*% t(Gamma[,,j]) %*% X[i,]
      }
      mu.mat[,,j] <- matrix(mu[,j],nrow=n,ncol=d,byrow=TRUE) # matrix version of mu_j
    }
    
    # mu <- rep(0,d)
    # den <- colSums(w)
    # for(i in 1:n){
    #   Delta.norm <- diag(w[i,]/den)
    #   mu <- mu + Gamma %*% Delta.norm %*% t(Gamma) %*% X[i,]
    # }
    # 
    # for(j in 1:k){
    #   for(h in 1:d){
    #     mu[h,j]  <- sum(z[,j]*correction[,h,j]/sum(z[,j]*correction[,h,j])*X[,h])
    #   }
    # }
    
    # Lambda & theta
    
    #ptm <- proc.time()
    for(j in 1:k){
      Xnew[,,j] <- t(t(Gamma[,,j]) %*% t(X-mu.mat[,,j]))
    }
    
    # Lambda (scalar update)
    
    for(j in 1:k){
      for(h in 1:d){
        Lambda[h,h,j] <- sum(z[,j]*correction[,h,j]*Xnew[,h,j]^2)/size[j]
      }
    }
    
    # eta
    
    for(j in 1:k){
      for(h in 1:d){
        eta.new <- sum(zvbad[,h,j]/sum(zvbad[,h,j])*Xnew[,h,j]^2)/Lambda[h,h,j]
        eta[h,j] <- max(eta.min,eta.new)
      }
    }
    
    # Check on alpha and eta
    
    for(j in 1:k){
      for(h in 1:d){
        if(eta[h,j] == eta.min)
          alpha[h,j] <- alpha.max
      }
    }
    
    # Gamma

    for(j in 1:k){
      Gamma[,,j] <- Gamma.update(X = X, z = z[,j], w = correction[,,j], mu=mu[,j], Lambda=Lambda[,,j], Gamma0 = Gamma[,,j], method = method, iter.max = iter.max, reltol = tol, trace = trace)
    }
  

    # +++++++++++++++++++ #
    # Density Computation #
    # +++++++++++++++++++ #

    for(j in 1:k){

      densX[,j] <- dmscn(x = X, mu = mu[,j], L = Lambda[,,j], G=Gamma[,,j], Sigma=NULL, alpha = alpha[,j], eta = eta[,j])
      dens[,j]  <- prior[j]*densX[,j]

    }

    # ------ #
    # E-step #
    # ------ #

    z <- dens/matrix(rep(rowSums(dens),k),ncol=k)
    
    for(j in 1:k){
      Xnew[,,j] <- t(t(Gamma[,,j]) %*% t(X-mu.mat[,,j]))
    }
    
    for(i in 1:n){
      for(j in 1:k){
        for(h in 1:d){
          den <- alpha[h,j]*dnorm(Xnew[i,h,j],mean=0,sqrt(Lambda[h,h,j])) + (1-alpha[h,j])*dnorm(Xnew[i,h,j],mean=0,sqrt(eta[h,j]*Lambda[h,h,j]))
          v[i,h,j] <- alpha[h,j]*dnorm(Xnew[i,h,j],mean=0,sqrt(Lambda[h,h,j]))/den
        }
      }
    }
    v[is.nan(v)] <- 0 # to avoid 0/0 in v

    for(j in 1:k){
      for(h in 1:d){
        zvgood[,h,j] <- z[,j]*v[,h,j]
        zvbad[,h,j]  <- z[,j]*(1-v[,h,j])
        correction[,h,j] <- v[,h,j]+(1-v[,h,j])/eta[h,j]
      }
    }

    # ------------------------------------- #
    # Global - Observed-data log-likelihood #
    # ------------------------------------- #

    llvalues <- sum(log(rowSums(dens)))
    loglik   <- c(loglik,llvalues)

    iteration <- iteration + 1
    cat("*")

    if(iteration == iter.max |  (loglik[iteration]-loglik[iteration-1])<tol)
      check <- 1

  }

  # cat("\n")

  # Sigma

  for(j in 1:k){
    Sigma[,,j] <- Gamma[,,j] %*% Lambda[,,j] %*% t(Gamma[,,j])
  }

  # -------------------- #
  # Final Log-likelihood #
  # -------------------- #

  final.loglik <- loglik[iteration]

  # ---------------------------------------------------------- #
  # Plot to check monotonicity of the penalized log-likelihood #
  # ---------------------------------------------------------- #

  if(plotll){
    oldpar <- par(no.readonly = TRUE)   
    on.exit(par(oldpar))          
    par(mai=c(0.84,0.8,0.012,0.004))
    par(las = 3)
    par(cex.axis=0.7)
    par(cex.lab=1.2)
    plot(1:(iteration-1), loglik[2:iteration], type="l",axes = FALSE, xlab="iterations", ylab="Objective function", lwd=2)
    axis(1, at = 1:(iteration-1), labels = 1:(iteration-1))
    axis(2)
    box(col = "black")

  }

  # The EM-algorithm is finished #

  # --------------------- #
  # Classification Matrix #
  # --------------------- #

  cluster <- apply(z,1,which.max)
  detect1 <- array(1/n,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  detect2 <- array(1/n,c(n,d),dimnames=list(1:n,paste("dim.",1:d,sep="")))
  for(h in 1:d){
    for(i in 1:n){
      detect2[i,h] <- ifelse(v[i,h,cluster[i]] > 1/2,"*","bad")
    }
    for(j in 1:k){
      detect1[,h,j] <- ifelse(v[,h,j] > 1/2,"*","bad")
    }
  }

  # -------------------- #
  # Number of parameters #
  # -------------------- #

  npar <- (k-1) + d*k + d*k + d*(d-1)/2*k + d*k + d*k

  # -------------------- #
  # Information criteria #
  # -------------------- #

  AIC <- - 2*final.loglik + 2*npar
  BIC <- - 2*final.loglik + npar*log(n)
  
  KIC <- - 2*final.loglik + 3*(npar+1)
  KICc <- - 2*final.loglik + 2*(npar + 1)*n/(n-npar -2) - n*digamma((n-npar)/2) + n*log(n/2)
  
  AIC3 <- - 2*final.loglik + 3*npar
  CAIC <- - 2*final.loglik + npar*(1+log(n))
  AICc <- - 2*final.loglik + 2*npar*n/(n - npar -1)
  
  #if (n-npar-1>0) {
  #  AICc  <- - 2*final.loglik + 2*npar - (2*npar*(npar+1))/(n-npar-1)
  #  AICu  <- AICc - n*log(n/(n-npar-1))
  #}
  #else {
  #  AICc  <- AICu <- -Inf
  #}
  # ICL (teigen style)
  
  ent <- apply(z,1,max)
  ICL <- BIC - sum(ent*log(ent))
  
  AWE <- - 2*(final.loglik + sum(ent*log(ent))) + 2*npar*(3/2 + log(n))
  CLC <- - 2*final.loglik + 2*sum(ent*log(ent))
  
#  
  
  ##z.const <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
  #hard.z  <- (matrix(rep(apply(z[bestD,],1,max),k),no.trim,k,byrow=F)==z[bestD,])*1

  ##ECM <- sum(hard.z*log(z.const))
  #ECM <- sum(hard.z*log(z[bestD,]))
  #ICL <- BIC + 2*ECM

  result <- list(
    X              = X,
    n              = n,
    d              = d,
    k              = k,
    cluster        = cluster,
    #detect1        = detect1,
    detect        = detect2,
    npar           = npar,
    pi          = prior,
    mu             = mu,
    Lambda         = Lambda,
    Gamma          = Gamma,
    Sigma          = Sigma,
    alpha          = alpha,
    eta            = eta,
    z              = z,
    v              = v,
    weight     = correction,##v+(1-v)/eta
    iter.stop      = iteration,
    loglik         = final.loglik,
    llv= loglik,
    AIC            = AIC,
    BIC            = BIC,
    ICL            = ICL,
    KIC            = KIC,
    KICc           = KICc,
    AWE            = AWE,
    AIC3           = AIC3,
    CAIC           = CAIC,
    AICc           = AICc,
    #AICu           = AICu,
    CLC            = CLC,
    call           = match.call()
  )
  class(result) <- "MSclust" 

  return(result)

}

#####################################################################################################################################################################################
################################################# MSt mixture #######################################################################################################################
#####################################################################################################################################################################################

dt <- function(x, mu = 0, scale = 1, theta = 100, formula = "direct"){
  
  # scale corresponds to sd of dnorm()
  
  if(scale < 0)
    stop("scale must be greater than, or equal to, 0")
  if(theta < 0)
    stop("theta must be greater than, or equal to, 0")
  
  delta <- (x-mu)^2/scale^2
  
  d <- 1
  
  if(formula=="direct"){
    
    # Equation (5) in Peel & McLachlan (2000), Statistics & Computing
    
    num <- gamma((theta + 1)/2)/scale
    den <- sqrt(pi*theta)*gamma(theta/2)*(1+delta/theta)^((theta + 1)/2)
    
    PDF <- num/den
    
  }
  
  if(formula=="indirect"){
    
    intf <- function(w,df,x,mu,scale){
      dnorm(x,mean=mu,sd=scale/sqrt(w))*dgamma(w,shape=df/2,rate=df/2)
    }
    
    PDF <- sapply(1:length(x), function(i) stats::integrate(intf, lower=0, upper=Inf, df = theta, x = x[i], mu = mu, scale = scale)$value)
    
  }
  
  return(PDF)
  
}

dmst <- function(x, mu = NULL, L = NULL, G = NULL, Sigma = NULL, theta = NULL, formula = "direct"){
  Lambda=L
  Gamma=G
   d = ifelse(is.null(x=ncol(x)),1,ncol(x))
  # Lambda is a matrix with diagonal elements
  if(is.null(mu)) mu=rep(0,d)
  if(is.null(theta)) theta = rep(100,d)
  # Lambda is a vector with diagonal elements
  
  if(min(theta) < 1)
    stop("each element in theta must be greater than, or equal to, 1")
  

  
  if(!is.null(Sigma)){
    temp   <- eigen(Sigma)
    Lambda <- diag(temp$values)
    Gamma  <- temp$vectors
  }
   if(is.null(Sigma)&is.null(G)){Sigma=diag(d)
   temp   <- eigen(Sigma)
   Lambda <- diag(temp$values)
   Gamma  <- temp$vectors}
   
  if(is.matrix(x))
    q <- nrow(x)
  if(is.vector(x) & d > 1)
    q <- 1
  if(is.vector(x) & d == 1)
    q <- length(x)
  
  dens <- array(0,c(q,d),dimnames=list(1:q,paste("X.",1:d,sep="")))
  
  for(i in 1:q){
    
    if(is.matrix(x))
      xnew <- t(Gamma) %*% (x[i,]-mu)
    if(is.vector(x) & d > 1)
      xnew <- t(Gamma) %*% (x-mu)
    if(is.vector(x) & d == 1)
      xnew <- x[i] - mu
    
    for(h in 1:d)
      dens[i,h] <- dt(x = xnew[h], mu = 0, scale = sqrt(Lambda[h,h]), theta = theta[h], formula = formula)
    
  }
  
  PDF <- apply(dens,1,prod)
  #PDF <- (PDF<10^(-322))*10^(-322)+(PDF>10^(-322))*PDF
  
  return(PDF)
  
}

rmst <- function(n, d = 2, mu = rep(0,d), L = NULL, G = NULL, Sigma = diag(d), theta = rep(100,d), n.dens = "dmnorm"){
  Lambda=L
  Gamma=G
  norm.dens=n.dens
  # Lambda is a vector with diagonal elements
  
  if(min(theta) < 1)
    stop("each element in theta must be greater than, or equal to, 1")
  
  if(is.vector(mu))
    d <- length(mu)
  if(!is.vector(mu))
    d <- 1
  
  if(!is.null(Sigma)){
    temp   <- eigen(Sigma)
    Lambda <- diag(temp$values)
    Gamma  <- temp$vectors
  }
  
  X <- w <- array(0,c(n,d),dimnames=list(1:n,paste("X.",1:d,sep="")))
  for(h in 1:d)
    w[,h] <- rgamma(n = n, shape = theta[h]/2, rate = theta[h]/2)
  
  (1 + stats::rexp(n = n, theta[h]))
  
  Delta <- SigmaMod <- array(0,c(n,d,d),dimnames=list(1:n,paste("X.",1:d,sep=""),paste("X.",1:d,sep="")))
  for(i in 1:n){
    Delta[i,,] <- diag(1/w[i,])
    SigmaMod[i,,] <- Gamma %*% Delta[i,,] %*% Lambda %*% t(Gamma)
    if(norm.dens == "dmnorm")
      X[i,] <- mnormt::rmnorm(n = 1, mean = mu, varcov = SigmaMod[i,,])
    if(norm.dens == "dmvnorm")
      X[i,] <- mvtnorm::rmvnorm(n = 1, mean = mu, sigma = SigmaMod[i,,])
  }
  
  return(X)
  
}

df.Mstep <- function(z, w, size, theta.min = 1, theta0 = NULL, theta.old = NULL, formula = "direct", method = "BFGS", iter.max = 5000, reltol = 1e-15, trace = 0){
  
  # w is a vector of weights
  
  if(!is.vector(w))
    w <- as.vector(w)
  if(!is.vector(z))
    z <- as.vector(z)
  if(!is.numeric(w))
    stop("numeric vector expected for w")
  if(!is.numeric(z))
    stop("numeric vector expected for z")
  if(any(is.na(w)) | any(is.na(z)))
    stop("No NAs allowed for both z and w")
  
  n <- length(z)
  
  # -------------- #
  # Initialization #
  # -------------- #
  
  # theta
  
  # starting-world
  
  if(is.null(theta0)){
    theta <- 100
  }else{
    if(theta0 < 1)
      stop("theta0 must be greater than, or equal to, theta.min")
    theta <- theta0
  }
  
  # working-world
  
  initial.value <- log(theta-theta.min)
  
  # objective function
  
  f <- function(par, z, w, theta.old, theta.min, size){
    
    # ----- #
    # theta #
    # ----- #
    
    theta <- exp(par) + theta.min
    
    # -------------- #
    # log-likelihood #
    # -------------- #
    
    # formula (5.47) in McLachlan & Krishnan (2007)
    
    f.obj <- -size*log(gamma(theta/2)) + theta*size/2*log(theta/2) + theta*size/2*( sum(z*(log(w) - w))/size + digamma((theta.old + 1)/2) - log((theta.old + 1)/2) ) 
    
    # f.obj <- 0
    # for(i in 1:n){
    #   f.obj <- f.obj + z[i]*( -log(gamma(theta/2)) + theta/2*log(theta/2) + theta/2*( sum(log(w) - w) + digamma((theta.old + 1)/2) - log((theta.old + 1)/2) ) )
    # }
    
    # f.obj <- -n*log(gamma(theta/2)) + n*theta/2*log(theta/2) + n*theta/2*(sum(log(w) - w)/n + digamma((theta.old + 1)/2) - log((theta.old + 1)/2))
    
    return(f.obj)
    
  }
  
  res <- optim(par = initial.value, fn = f, z = z, w = w, theta.old = theta.old, theta.min = theta.min, size = size, method = method, control = list(maxit=iter.max,reltol=reltol,trace=trace,fnscale=-1))
  
  f.obj <- res$value
  est   <- res$par
  
  theta <- exp(est) + theta.min
  
  return(theta)
  
}

#########################
## ML for MSt mixtures ##
#########################

mst <- function(
  X,                 # data matrix
  k,                 # number of mixture components
  ini = "km", # or manual; in that case, a partition (init.clus) must be provided
  sz = NULL,
  df.min = 1,      # minimum proportion of good points in each group for the contaminated normal distribution
  dfU = "num", # criterion to update the degrees of freedom
  frm = "dir",
  m = "BFGS",
  stop = c(10^-5,200),     # the 2 stopping rules tollerance and maximum number of iterations
  VB = FALSE # If TRUE, tracing information on the progress of the optimization is produced; see optim() for details
  #and plot of the log-likelihood versus iterations
  )
{
  frml=frm
  md=m
  if(dfU == "num"){dfU = "numeric"}
  df.update=dfU
  if(frml == "dir"){frml = "direct"}
  formula=frml
  start.z=sz
  initialization=ini
  method =md
  tol=stop[1]
  iter.max=stop[2]
  trace = 0
  plotll = FALSE
  if (VB==TRUE){
    trace = 1
    plotll = TRUE
  }
  
  
  
  if(is.data.frame(X))
    X <- data.matrix(X)
  if(!is.matrix(X))
    X <- matrix(X, nrow = length(X), ncol = 1, dimnames = list(names(X), deparse(substitute(X))))
  if(!is.numeric(X))
    stop("numeric matrix/vector expected for X")
  if(any(is.na(X)))
    stop("No NAs allowed")
  if(is.null(k))
    stop("The number of groups k is NULL")
  
  n <- nrow(X) # number of units
  d <- ncol(X) # number of variables
  
  Xnew <- array(0,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  
  # --------------------- #
  # Definition of objects #
  # --------------------- #
  
  # E-step quantities
  
  z     <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep="")))
  v     <- array(1,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  delta <- array(NA,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))

  # M-step Parameters
  
  prior  <- numeric(k)
  mu     <- array(0,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  mu.mat <- array(0,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Sigma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Lambda <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  Gamma  <- array(0,c(d,d,k),dimnames=list(paste("X.",1:d,sep=""),paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  df     <- array(100,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
  
  # Distribution
  
  densX <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # fX
  dens  <- array(0,c(n,k),dimnames=list(1:n,paste("comp.",1:k,sep=""))) # weights*fX
  
  # ------------------------- #
  # posteriors initialization #
  # ------------------------- #
  
  if(initialization=="random.soft"){
    
    z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k)
    z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
    
  }
  
  if(initialization=="random.hard"){
    
    z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)
    
  }
  
  if(initialization=="manual"){ # z.start can be both soft and hard initialization
    
    z  <- start.z # soft or hard classification
    
  }
  
  if(initialization=="mclust"){
    
    if(d==1)
      mclustfit <- mclust::Mclust(data=X, G=k, modelNames="V")
    if(d>1)
      mclustfit <- mclust::Mclust(data=X, G=k, modelNames="VVV")
    
    z <- mclustfit$z
    
  }
  
  if(initialization=="km"){
    
    if(k>1){
      km       <- kmeans(X,k)
      groups   <- km$cluster
      z        <- mclust::unmap(classification=groups)
    }
    if(k==1){
      groups   <- rep(1,n)
      z        <- matrix(groups,nrow=n,ncol=k)
    }
    
  }
  
  if(initialization=="pam"){
    
    lk     <- cluster::pam(x=X,k=k)
    groups <- lk$clustering
    z      <- mclust::unmap(classification=groups)
    
  }
  
  # initialization of Gamma
  
  for(j in 1:k){
    xbar <- colSums(z[,j]/sum(z[,j])*X)
    S    <- crossprod(sqrt(z[,j]/sum(z[,j]))*(X-matrix(rep(xbar,n),n,d,byrow=TRUE)))
    #S <- as.matrix(forceSymmetric(x=S))
    temp  <- eigen(S)
    Gamma[,,j] <- temp$vectors
  }
  
  # EM algorithm --------------------------------------------------
  
  # Preliminary definition of convergence criterions
  
  check     <- 0
  iteration <- 1
  loglik    <- -Inf
  
  while(check<1){
    
    # ++++++ #
    # M-Step #
    # ++++++ #
    
    # Parameters
    
    # pi
    
    prior <- colMeans(z)
    size  <- colSums(z)
    
    # mu
    
    mu <- array(0,c(d,k),dimnames=list(paste("X.",1:d,sep=""),paste("comp.",1:k,sep="")))
    for(j in 1:k){
      den <- colSums(z[,j]*v[,,j])
      for(i in 1:n){
        Delta.norm <- diag(z[i,j]*v[i,,j]/den)
        mu[,j] <- mu[,j] + Gamma[,,j] %*% Delta.norm %*% t(Gamma[,,j]) %*% X[i,]
      }
      mu.mat[,,j] <- matrix(mu[,j],nrow=n,ncol=d,byrow=TRUE) # matrix version of mu_j
    }
    
    # mu <- rep(0,d)
    # den <- colSums(w)
    # for(i in 1:n){
    #   Delta.norm <- diag(w[i,]/den)
    #   mu <- mu + Gamma %*% Delta.norm %*% t(Gamma) %*% X[i,]
    # }
    # 
    # for(j in 1:k){
    #   for(h in 1:d){
    #     mu[h,j]  <- sum(z[,j]*correction[,h,j]/sum(z[,j]*correction[,h,j])*X[,h])
    #   }
    # }
    
    # Lambda & df
    
    for(j in 1:k){
      Xnew[,,j] <- t(t(Gamma[,,j]) %*% t(X-mu.mat[,,j]))
    }
    
    # Lambda (scalar update)
    
    for(j in 1:k){
      for(h in 1:d){
        Lambda[h,h,j] <- sum(z[,j]*v[,h,j]*Xnew[,h,j]^2)/size[j]
      }
    }
    
    # df
    
    for(j in 1:k){
      for(h in 1:d){
        
        if(df.update == "numeric"){
          df[h,j] <- df.Mstep(z = z[,j], w = v[,h,j], size = size[j], theta.min = df.min, theta0 = df[h,j], theta.old = df[h,j], formula = formula, method = method, iter.max = iter.max, reltol = tol, trace = trace)
        }
        
        if(df.update == "approx"){ 
          kk      <- -1 - sum(z[,j]*(log(v[,h,j]) - v[,h,j]))/size[j] - digamma((df[h,j] + 1)/2) + log((df[h,j] + 1)/2)
          df[h,j] <- (-exp(kk) + 2*exp(kk)*(exp(digamma(df[h,j]/2)) - (df[h,j]/2 - 1/2)))/(1-exp(kk))
        }
        
      }
    }
    
    # for(j in 1:k){
    #   for(h in 1:d){
    #     eta.new <- sum(zvbad[,h,j]/sum(zvbad[,h,j])*Xnew[,h,j]^2)/Lambda[h,h,j]
    #     eta[h,j] <- max(eta.min,eta.new)
    #   }
    # }
    
    # Gamma
    
    for(j in 1:k){
      Gamma[,,j] <- Gamma.update(X = X, z = z[,j], w = v[,,j], mu = mu[,j], Lambda = Lambda[,,j], Gamma0 = Gamma[,,j], method = method, iter.max = iter.max, reltol = tol, trace = trace)
    }
    #print(Gamma)
    
    # +++++++++++++++++++ #
    # Density Computation #
    # +++++++++++++++++++ #
    
    for(j in 1:k){
      
      densX[,j] <- dmst(x = X, mu = mu[,j], L = Lambda[,,j], G = Gamma[,,j], Sigma = NULL, theta = df[,j], formula = formula)
      dens[,j]  <- prior[j]*densX[,j]
      
    }
    
    # ------ #
    # E-step #
    # ------ #
    
    z <- dens/matrix(rep(rowSums(dens),k),ncol=k)
    
    for(j in 1:k){
      Xnew[,,j] <- t(t(Gamma[,,j]) %*% t(X-mu.mat[,,j]))
    }
    
    for(i in 1:n){
      for(j in 1:k){
        delta[i,,j] <- Xnew[i,,j]^2/diag(Lambda[,,j])
        for(h in 1:d){
          
          num      <- df[h,j] + 1
          den      <- df[h,j] + delta[i,h,j]
          v[i,h,j] <- num/den
          
        }
      }
    }
    #v[is.nan(v)] <- 0 # to avoid 0/0 in v
    
    # ------------------------------------- #
    # Global - Observed-data log-likelihood #
    # ------------------------------------- #
    
    llvalues <- sum(log(rowSums(dens)))
    #print(llvalues)
    loglik   <- c(loglik,llvalues)
    
    iteration <- iteration + 1
    #cat("*")
    
    if(iteration == iter.max | (loglik[iteration]-loglik[iteration-1])<tol)
      check <- 1
    
  }
  
  # cat("\n")
  
  # Sigma
  
  for(j in 1:k){
    Sigma[,,j] <- Gamma[,,j] %*% Lambda[,,j] %*% t(Gamma[,,j])
  }
  
  # -------------------- #
  # Final Log-likelihood #
  # -------------------- #
  
  final.loglik <- loglik[iteration]
  
  # ---------------------------------------------------------- #
  # Plot to check monotonicity of the penalized log-likelihood #
  # ---------------------------------------------------------- #
  
  if(plotll){
    oldpar <- par(no.readonly = TRUE)    # code line i
    on.exit(par(oldpar))            # code line i + 1
    par(mai=c(0.84,0.8,0.012,0.004))
    par(las = 3)
    par(cex.axis=0.7)
    par(cex.lab=1.2)
    plot(1:(iteration-1), loglik[2:iteration], type="l",axes = FALSE, xlab="iterations", ylab="Objective function", lwd=2)
    axis(1, at = 1:(iteration-1), labels = 1:(iteration-1))
    axis(2)
    box(col = "black")
    
  }
  
  # The EM-algorithm is finished #
  
  # --------------------- #
  # Classification Matrix # DA AGGIORNARE
  # --------------------- #
  
  cluster <- apply(z,1,which.max)
  detect1 <- array(1/n,c(n,d,k),dimnames=list(1:n,paste("dim.",1:d,sep=""),paste("comp.",1:k,sep="")))
  detect2 <- array(1/n,c(n,d),dimnames=list(1:n,paste("dim.",1:d,sep="")))
  for(h in 1:d){
    for(i in 1:n){
      detect2[i,h] <- ifelse(v[i,h,cluster[i]] > 1/2,"*","bad")
    }
    for(j in 1:k){
      detect1[,h,j] <- ifelse(v[,h,j] > 1/2,"*","bad")
    }
  }
  
  # -------------------- #
  # Number of parameters #
  # -------------------- #
  
  npar <- (k-1) + d*k + d*k + d*(d-1)/2*k + d*k
  
  # -------------------- #
  # Information criteria #
  # -------------------- #
  
  AIC <- - 2*final.loglik + 2*npar
  BIC <- - 2*final.loglik + npar*log(n)
  
  KIC <- - 2*final.loglik + 3*(npar+1)
  KICc <- - 2*final.loglik + 2*(npar + 1)*n/(n-npar -2) - n*digamma((n-npar)/2) + n*log(n/2)
  
  AIC3 <- - 2*final.loglik + 3*npar
  CAIC <- - 2*final.loglik + npar*(1+log(n))
  AICc <- - 2*final.loglik + 2*npar*n/(n - npar -1)
  
  if (n-npar-1<0) {
    AICc <- - 2*final.loglik + 2*npar*n/(n - npar -1)
  # because we use AICc when npar is larger n
    
  #  AICu  <- AICc - n*log(n/(n-npar-1))
  }
  else {
    AICc   <- Inf
  }
  
  # ICL (teigen style)
  
  ent <- apply(z,1,max)
  ICL <- BIC - sum(ent*log(ent))
  
  AWE <- - 2*(final.loglik + sum(ent*log(ent))) + 2*npar*(3/2 + log(n))
  CLC <- - 2*final.loglik + 2*sum(ent*log(ent))
  
  ##z.const <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
  #hard.z  <- (matrix(rep(apply(z[bestD,],1,max),k),no.trim,k,byrow=F)==z[bestD,])*1
  
  ##ECM <- sum(hard.z*log(z.const))
  #ECM <- sum(hard.z*log(z[bestD,]))
  #ICL <- BIC + 2*ECM
  
  result <- list(
    X              = X,
    n              = n,
    d              = d,
    k              = k,
    cluster        = cluster,
    #detect1        = detect1,
    detect        = detect2,##if the point is bad or not per each principal component given the cluster membership
    npar           = npar,
    pi          = prior,
    mu             = mu,
    Lambda         = Lambda,
    Gamma          = Gamma,
    Sigma          = Sigma,
    df             = df,
    z              = z,
    v              = v,
    iter.stop      = iteration,
    loglik         = final.loglik,
    AIC            = AIC,
    BIC            = BIC,
    ICL            = ICL,
    KIC            = KIC,
    KICc           = KICc,
    AWE            = AWE,
    AIC3           = AIC3,
    CAIC           = CAIC,
    AICc           = AICc,
    #AICu           = AICu,
    CLC            = CLC,
    call           = match.call()
  )
  class(result) <- "MSclust"
  
  return(result)
  
}