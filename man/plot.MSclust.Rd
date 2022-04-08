\name{plot.MSclust}
\alias{plot.MSclust}
\title{MSclust Plotting}
\usage{
\method{plot}{MSclust}(x, ...)
}
\arguments{
\item{x}{A \code{MSclust} object.}

\item{...}{Arguments to be passed to methods, such as graphical parameters.}
}
\value{
No return value, called to visualize the fitted model's results
}
\description{
MSclust Plotting
}
\examples{

## Not run: 
\dontrun{
data(sim)
result <- mscn(X = sim, k = 2, initialization = "kmeans", method = "BFGS")
plot(result)}
## End(Not run) 

}
