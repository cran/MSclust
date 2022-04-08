
\name{summary.MSclust}
\alias{summary.MSclust}
\title{Summary for MSclust}
\usage{
\method{summary}{MSclust}(object, ...)
}
\arguments{
\item{object}{A \code{MSclust} object.}

\item{...}{Arguments to be passed to methods, such as graphical parameters.}
}
\value{
No return value, called to summarize the fitted model's results
}
\description{
Summarizes main information regarding a \code{MSclust} object.
}
\details{
Information includes clustering table, total outliers, outliers per cluster, mixing proportions,
  component means and variances.
}
\examples{
## Not run: 
\dontrun{
data(sim)
result <- mscn(X = sim, k = 2, initialization = "kmeans", method = "BFGS")
summary(result)}
## End(Not run) 
}
