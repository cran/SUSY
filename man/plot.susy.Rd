\name{plot.susy}
\alias{plot.susy}
\title{ susy plot method }
\description{
  Generate plot(s) for a \code{susy} object.
}
\usage{
  \method{plot}{susy}(x, type=c(4, 5), \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{type}{ Numeric, specifies the types of plot, defaults to \code{c(4, 5)}. \enumerate{
    \item{\emph{GMcrosscorrs}}
    \item{\emph{synchrony by segments}}
    \item{\emph{GM-Z}}
    \item{\emph{time series plot}}
    \item{\emph{Z not abs}}
  }}
  \item{\dots}{ Ignored. }
}
\details{
  Method can generate multiple types of plots by providing numeric vector to \code{type} argument. Note it will generate plots for each pair (cross computation) in \code{x}, so the final number of plots is \code{length(x) * length(type)}.
}
\value{
  Returns \code{NULL} invisibly. Generate plot(s) as a side effect.
}
\seealso{\code{\link{susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n, 300, 330),
  var2 = runif(n, 300, 330),
  var3 = runif(n, 300, 330)
)
res = susy(data, segment=30L, Hz=15L, permutation=TRUE)
plot(res, type=c(3,5))
}
\keyword{ data }
