\name{rmst}
\alias{rmst}
\alias{dmst}
\title{Multiple Scaled Student-t Distribution 
}
\description{Probability density function and pseudo-random number generation for the multiple scaled Student-t distribution.

}
\usage{
  dmst(x, mu = NULL, L = NULL, G = NULL, Sigma = NULL, theta = NULL, formula = "direct")
rmst(n,d=2,mu=rep(0,d),L=NULL,G=NULL,Sigma=diag(d),theta=rep(100,d),n.dens="dmnorm")
}


\arguments{
	\item{x}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.}
  \item{n}{The number of observations to be generated.}
  \item{d}{A number specifing the dimension.}
  \item{mu}{Either a vector of length \code{d}, representing the mean value, or (except for \code{rmst}) a matrix whose rows represent different mean vectors; if it is a matrix, its dimensions must match those of \code{x}.}
  \item{L}{Lambda diagonal \code{d}-dimensional matrix of the eigenvalues of \code{Sigma}.}
  \item{G}{Gamma orthogonal \code{d}-dimensional matrix whose columns are the normalized eigenvectors of \code{Sigma}.}
  \item{Sigma}{A symmetric positive-definite \code{d}-dimensional matrix representing the scale matrix of the distribution; a vector of length 1 is also allowed (in this case, \code{d = 1} is set). Identity matrix by default.}
  \item{theta}{Vector of dimesion \code{d} containing the degrees of freedom.}  
  \item{n.dens}{\code{"dmnorm"} or \code{"dmvnorm"}, dependening on the way the density of the normal distribution is computed.}
  \item{formula}{\code{"direct"} or \code{"indirect"}; if \code{"direct"}, then Equation (5) in Peel & McLachlan (2000), Statistics & Computing is used.}
}
\value{
  \item{dmscn}{returns a vector of density values.}
  \item{rmscn}{returns a matrix of \code{n} rows of observations.}
}

\references{
Punzo, A., & Tortora, C. (2021). \emph{Multiple scaled contaminated normal distribution and its application in clustering.} \emph{Statistical Modelling}, \bold{21}(4): 332--358.

Forbes, F. & Wraith, D. (2014). \emph{A new family of multivariate heavy-tailed distributions with variable marginal amounts of tailweight: application to robust clustering.} \emph{Statistics and Computing}, \bold{24}(6), 971--984. %% for mst
}
\author{Cristina Tortora and Antonio Punzo}

\examples{
x <- matrix(c(0,0),1,2)
theta <- c(5,20)
density <- dmst(x = x, theta = theta)
density

n <- 100
random <- rmst(n = n, theta = theta)
plot(random)
}

