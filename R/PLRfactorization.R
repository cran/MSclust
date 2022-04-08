#-----------------------------#
# Strictly elimination matrix #
#-----------------------------#

Lbar <- function(n){
  
  Lb      <- matrix(0,(n*(n-1)/2),n^2)
  indcol  <- sapply(0:(n-2), function(i) 
    sapply((i*n+(i+2)):(n+i*n), function(j) j )
  )
  
  c1      <- 1:(n*(n-1)/2)
  c2      <- unlist(indcol)
  ind     <- cbind(c1,c2)
  Lb[ind] <- 1
  
  return(Lb)
  
}

#-------------------------------------#
# Stricly half-vectorization operator #
#-------------------------------------#

bvech <- function(A){   # A: square marrix
  
  n   <- ncol(A)
  res <- Lbar(n)%*%c(A)
  
  return(res)
  
}

#--------------#
# From Q to Pl #
#--------------#

QPl <- function(Q   # Orthogonal matrix
                ){
  
  luQ    <- Matrix::lu(Q)
  eluQ   <- Matrix::expand(luQ)
  l      <- bvech(as.matrix(eluQ$L))
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
               P=P,  # Permutation matrix
               l=l   # Vector of n(n-1)/2 elements
              )
            )
  
}


#--------------#
# From Q to PL #
#--------------#

QPL <- function(Q   # Orthogonal matrix
                ){
  
  luQ    <- lu(Q)
  eluQ   <- expand(luQ)
  L      <- eluQ$L
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
              P=P,  # Permutation matrix
              L=L   # Unit lower triangular matrix
              )
         )
  
}

#--------------#
# From Pl to Q #
#--------------#

PlQ <- function(P,  # Permutation matrix
                l   # Vector of n(n-1)/2 elements
){
  
  n  <- (1+sqrt(1+8*length(l)))/2
  L  <- diag(n) + matrix(t(Lbar(n))%*%l,n,n) 
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}

#--------------#
# From PL to Q #
#--------------#

PLQ <- function(P,  # Permutation matrix
                L   # Unit lower triangular matrix
                ){
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}


