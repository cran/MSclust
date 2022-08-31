\name{mscn}
\alias{mscn}
\title{Mixtures of Multiple Scaled Contaminated Normal Distributions.}
\description{Fits a mixture of multiple scaled contaminated normal distributions to the given data.}
\usage{
mscn(X,k,ini="km",sz=NULL,al=c(0.5,0.99),eta.min=1.01,m="BFGS",stop=c(10^-5,200),VB=FALSE)
}


\arguments{
  \item{X}{A matrix or data frame such that rows correspond to observations and columns correspond to variables.}
  \item{k}{The number of clusters.}
  \item{ini}{Using kmeans by default or \code{"pam"} for partition around medoids, \code{"mclust"} for Gaussian mixture models, \code{"random.soft"} or \code{"random.hard"} for random or manual; if \code{"manual"}, a partition (\code{sz}) must be provided.}
  \item{sz}{If initialization is \code{"manual"}, this matrix contains the starting values for \code{z}.}
  \item{al}{2-dimensional vector containing minimum and maximum proportion of good points in each group for the contaminated normal distribution.}
  \item{eta.min}{Minimum value for inflation parameter for the covariance matrix for the bad points.}
  \item{m}{Method for the optimization of the eigenvector matrix, see \code{\link{optim}} for other options.}
  \item{stop}{2-dimensional vector with the Aitken criterion stopping rule and maximum number of iterations.}  
  \item{VB}{If \code{TRUE}, tracing information on the progress of the optimization is produced; see \code{\link{optim}} for details and plotting of the log-likelihood versus iterations.}  
  
}
\value{
  \item{X}{Data used for clustering.}
  \item{n}{The number of observations in the data.}
  \item{d}{The number of features in the data.}
  \item{k}{Value corresponding to the number of components.}
  \item{cluster}{Vector of group membership as determined by the model.}
  \item{detect}{Detect if the point is bad or not per each principal component given the cluster membership.}
  \item{npar}{The number of parameters.}
  \item{mu}{Either a vector of length \code{d}, representing the mean value, or (except for \code{rmscn}) a matrix whose rows represent different mean vectors; if it is a matrix, its dimensions must match those of \code{x}.}
  \item{Lambda}{Orthogonal matrix whose columns are the normalized eigenvectors of Sigma.}
  \item{Gamma}{Diagonal matrix of the eigenvalues of Sigma.}
  \item{Sigma}{A symmetric positive-definite matrix representing the scale matrix of the distribution.}
  \item{alpha}{Proportion of good observations.}
  \item{eta}{Degree of contamination.}
  \item{z}{The component membership of each observations.}
  \item{v}{The indicator if an observation is good or bad with respect to each dimension; 1 is good, and 0 means bad.}
  \item{weight}{The matrix of the expected value of the characteristic weights; corespond to the value of \code{v+(1-v)/eta}.}
  \item{iter.stop}{The number of iterations until convergence for the model.}
  \item{loglik}{The log-likelihood corresponding to the model.}
  \item{AIC}{The Akaike's Information Criterion of the model.}
  \item{BIC}{The Bayesian Information Criterion of the model.}
  \item{ICL}{The Integrated Completed Likelihood of the model.}
  \item{KIC}{The Kullback Information Criterion of the model.}
  \item{KICc}{The Bias correction of the Kullback Information Criterion of the model.}
  \item{AWE}{The Approximate Weight of Evidence of the model.}
  \item{AIC3}{Another version of Akaike's Information Criterion of the model.}
  \item{CAIC}{The Consistent Akaike's Information Criterion of the model.}
  \item{AICc}{The AIC version which is used when sample size \code{n} is small relative to \code{d}.}
  \item{CLC}{The Classification Likelihood Criterion of the model.}
}

\references{
  Punzo, A. & Tortora, C. (2021). \emph{Multiple scaled contaminated normal distribution and its application in clustering.} \emph{Statistical Modelling}, \bold{21}(4): 332--358.
	}
\author{Cristina Tortora and Antonio Punzo}


\examples{
## Not run:
\dontrun{
data(sim)
result <- mscn(X = sim, k = 2)
plot(result)
summary(result)}
## End(Not run) 
}

