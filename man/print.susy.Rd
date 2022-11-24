\name{print.susy}
\alias{print.susy}
\title{ susy print method }
\description{
  Prints information about an \code{susy} object.
}
\usage{
  \method{print}{susy}(x, corr.no.abs=TRUE, legacy=FALSE, \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{corr.no.abs}{ Logical, defaults to \code{TRUE} display correlation without the absolute value. }
  \item{legacy}{ Logical, defaults to \code{FALSE}, when \code{TRUE} print will produce an output that matches the output of legacy SUSY implementation. }
  \item{\dots}{ Extra arguments passed to \code{print.data.frame} method. }
}
\value{
  Returns \code{x} invisibly. Display output to console as a side effect.
}
\seealso{\code{\link{susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n, 300, 330),
  var2 = runif(n, 300, 330)
)
res = susy(data, segment=30L, Hz=15L)
res
print(res, corr.no.abs=FALSE)
print(res, digits=4)
print(res, legacy=TRUE)
}
\keyword{ data }
