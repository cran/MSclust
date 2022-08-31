plot.MSclust <- function(x, ...) {
  
  pch_good       <- 16
  pch_bad        <- 8
  
  
  
  model    <- x$model
  clusters <- as.character(x$cluster)
  outliers <- x$detect
  dat      <- x$X
  d        <- ncol(dat)
  d        <- ifelse(is.null(d), 1, d)
  
  if (d > 1) {
    if (d < 10) {
      print(ggparcoord(data = cbind(dat, clusters=clusters), mapping=ggplot2::aes(color=as.factor(clusters)), columns = 1:d) +
              theme_bw() + ggtitle('Parallel Coordinate Plot')+scale_color_discrete("Clusters",labels=levels(clusters)))
      
      if (d > 2) {
        print(pairs(dat, col = clusters, 
                    main = 'Cluster Memberships'))
      } else {
        print(plot(dat, col = clusters, 
                   main = 'Cluster Memberships'))
      }
    }
    else {
      print(ggparcoord(data = cbind(clusters=clusters, dat), mapping=ggplot2::aes(color=as.factor(clusters)), columns = 2:11) +
              theme_bw() + ggtitle('Parallel Coordinate Plot - First 10 varaibles')+scale_color_discrete("Clusters",labels=levels(clusters)))
    }
  }
  
  #++++ Log-likelihood over iterations ++++#
  # plot(x$loglik, type = 'b', pch = 16, xlab = 'Iteration', ylab = 'Log-Likelihood')
}



summary.MSclust <- function(object, ...) {
  
  cat('\nIterations:', object$iter.stop)
  
  cat("\n\nClustering table:")
  print(table(object$cluster))
  
  
  cat('\nMixing proportions:\n')
  print(object$pi)
  
  cat('\nComponent means:\n')
  print(object$mu)
  
  cat('\nComponent variances:\n')
  print(object$Sigma)
  
  cat('\nInformation Criteria:\n')
  print(data.frame(
    AIC  = object$AIC,
    BIC  = object$BIC,
    KIC  = object$KIC,
    KICc = object$KICc,
    AIC3 = object$AIC3,
    CAIC = object$CAIC,
    AICc = object$AICc,
    ICL  = object$ICL,
    AWE  = object$AWE,
    CLC  = object$CLC
  ))
}
