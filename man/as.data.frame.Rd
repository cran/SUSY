\name{as.data.frame.susy}
\alias{as.data.frame.susy}
\title{ susy to data.frame conversion method }
\description{
  Turns \code{susy} class object into a \code{data.frame}.
}
\usage{
  \method{as.data.frame}{susy}(x, row.names=NULL, optional=FALSE, corr.no.abs=TRUE, \dots)
}
\arguments{
  \item{x}{ A \code{susy} object. }
  \item{row.names}{ Ignored, only for consistency to generic \code{as.data.frame} method. }
  \item{optional}{ Ignored, only for consistency to generic \code{as.data.frame} method. }
  \item{corr.no.abs}{ Logical, defaults to \code{TRUE} display correlation without the absolute value. }
  \item{\dots}{ Ignored. }
}
\value{
  Returns \code{data.frame}.
}
\seealso{\code{\link{susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n, 300, 330),
  var2 = runif(n, 300, 330)
)
res = susy(data, segment=30L, Hz=15L)
as.data.frame(res)
}
\keyword{ data }
