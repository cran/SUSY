\name{susy}
\alias{susy}
\title{ Surrogate Synchrony }
\description{
  Cross-correlations of two time series are computed up to a specific lag in seconds \code{maxlag}. Cross-correlation is done within segment of the time series. The size of segments \code{segment} can be chosen in seconds. Aggregation is then performed by transforming correlations to Fisher's Z, computing mean Z in each segment, then across all segments of the time series. Segment shuffling is used to create surrogate time series, on which the same computations are run. This provides effect sizes \code{ES}. SUSY provides these different synchrony measures for each twin time series: \code{mean Z} and \code{ES of mean Z}; \code{mean absolute_Z} and \code{ES of mean absolute_Z}.
}
\usage{
susy(x, segment, Hz, maxlag=3L, permutation=FALSE,
  restrict.surrogates=FALSE, surrogates.total=500)
}
\arguments{
  \item{x}{ A data.frame of numeric columns. }
  \item{segment}{ Integer, size in seconds. Must not be smaller than \code{2 * maxlag}, must not be larger than half the the time series (\code{nrow(x)/2}). }
  \item{Hz}{ Integer, frames per second (sampling rate). }
  \item{maxlag}{ Integer, maximum lag for \code{\link[stats]{ccf}} in seconds. Default \code{3} seconds. }
  \item{permutation}{ Logical, default \code{FALSE} requires \code{x} to have even number of columns which are processed in pairs (1-2, 3-4, etc.). When \code{permutation} is \code{TRUE} then function computes all pairs combinations between columns provided in \code{x} (\code{n*(n-1)/2} pairs). }
  \item{restrict.surrogates}{ Logical, default \code{FALSE}. Restrict the number of surrogates or not. }
  \item{surrogates.total}{ Numeric, the number of generated surrogates, default \code{500}. Ignored when \code{restrict.surrogates} is \code{FALSE} (default). }
}
\details{
  Segments are non-overlapping, and the number of segments that fit into the time series may have a remainder (usually a few seconds at the end of the time series), which is not considered.
}
\value{
  Object of class \code{susy} is returned. Each cross correlation pair is an element in resulting object.
}
\seealso{\code{\link{plot.susy}}, \code{\link{as.data.frame.susy}}, \code{\link{print.susy}}}
\examples{
n = 1000
data = data.frame(
  var1 = runif(n, 300, 330),
  var2 = runif(n, 300, 330),
  var3 = runif(n, 300, 330)
)

## use only first two columns
res = susy(data[, 1:2], segment=30L, Hz=15L)
length(res)
names(res)

## use all columns and permutation
res = susy(data, segment=30L, Hz=15L, permutation=TRUE)
length(res)
names(res)

## print susy
res
print(res, legacy=TRUE)

## plot susy
plot(res)
plot(res, type=1:2)
}
\keyword{ data }
