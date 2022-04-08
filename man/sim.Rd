
\docType{data}
\name{sim}
\alias{sim}
\title{A Mixture of Two Normal Distributions with outliers}
\format{
A matrix with 600 observations (rows) and 2 variables (columns).
  The first 180 rows belong to cluster 1, and the last 420 rows belong to cluster 2.
\describe{
  \item{V1}{variable 1.}
  \item{V2}{variable 2.}
}
}
\source{
Punzo, A. & Tortora, C. (2021). \emph{Multiple scaled contaminated normal distribution and its application in clustering.} \emph{Statistical Modelling}, \bold{21}(4): 332--358.
}
\usage{
data(sim)
}
\description{
A simulated mixture of two normal distributions with mean (0,0) and (2,4), respectively, the fist covariance matrix has diagonlas equal to 2 and covariance 1.5, the second 0.5 and 0.1 respectively. Moreover, 18 observations have been transformed in outliers. For details, see Punzo, A. & Tortora, C. (2021). \emph{Multiple scaled contaminated normal distribution and its application in clustering.} \emph{Statistical Modelling}, \bold{21}(4): 332--358.
}
\keyword{datasets}
