\name{rmscn}
\alias{rmscn}
\alias{dmscn}
\title{Multiple Scaled Contaminated Normal Distribution 
}
\description{Probability density function and pseudo random number generation for the multiple scaled contaminated normal distribution.

}
\usage{
  dmscn(x,  mu = NULL, L = NULL, G = NULL, Sigma = NULL, alpha = NULL, eta = NULL)
rmscn(n,d=2,mu=rep(0,d),L=NULL,G=NULL,Sigma=diag(d),alpha=rep(0.99,d),eta=rep(1.01,d))
}


\arguments{
  \item{x}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.}
 \item{n}{The number of random vectors to be generated.}
    \item{d}{A number specifing the dimenstion.}
  \item{mu}{Either a vector of length \code{d}, representing the mean value, or (except for rmscn) a matrix whose rows represent different mean vectors; if it is a matrix, its dimensions must match those of \code{x}.}
 
  \item{L}{Lambda diagonal \code{d}-dimensional matrix of the eigenvalues of Sigma.}
  \item{G}{Gamma orthogonal \code{d}-dimensional matrix whose columns are the normalized eigenvectors of Sigma.}
  \item{Sigma}{A symmetric positive-definite \code{d}-dimensional matrix representing the scale matrix of the distribution; a vector of length 1 is also allowed (in this case, \code{d = 1} is set). Identity matrix by default.}
  \item{alpha}{\code{d}-dimensional vector containing the proportion of good observations; it must be a number between 0 and 1.}
  \item{eta}{\code{d}-dimensional vector containing  the degree of contamination; it should be a number greater than 1.}   
  
}
\value{
  \item{dmscn}{returns a vector of density values.}
  \item{rmscn}{returns a matrix of \code{n} rows of observations.}
}

\references{
Punzo, A. & Tortora, C. (2021). \emph{Multiple scaled contaminated normal distribution and its application in clustering.} \emph{Statistical Modelling}, \bold{21}(4): 332--358.
}
\author{Cristina Tortora and Antonio Punzo}

\examples{
x <- matrix(c(0,0),1,2)
alpha <- c(0.8,0.6)
eta <- c(2,4)
density <- dmscn(x = x, alpha = alpha, eta = eta)
density

n <- 100
random <- rmscn(n = n, alpha = alpha, eta = eta)
plot(random)
}

