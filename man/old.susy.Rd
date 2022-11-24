\name{old.susy}
\alias{old.susy}
\title{ Surrogate Synchrony - legacy}
\description{
  An internal function, a legacy implementation to compute Surrogate Synchrony, please use \code{\link{susy}} function instead.\cr
  Cross-correlations are computed up to a specific lag in seconds \code{Maxlag}, then aggregated within a chosen \code{Segment} of e.g. 30s. Aggregation is performed by transforming correlations to Fisher's Z, then computing mean Z in each segment, then across all segments of the time series. Segment shuffling is used to create surrogate time series, on which the same computations are run. This provides effect sizes \code{ES}. \code{File}: The pairs of time series are in the columns of the file, variable names are in the header line. If \code{Automatic} is clicked, the synchrony is computed of all adjacent pairs of columns in the file. If \code{Automatic} is unclicked, you may choose the two columns to be analyzed for synchrony, and two plots are additionally prepared. SUSY provides two different synchrony measures of each twin time series: mean Z and ES of mean Z; mean absolute_Z and ES of mean absolute_Z.
}
\usage{
old.susy(
  filename="MEA-topstage-C3-Beethoven.txt", separator=" ", missings=".", datahead=TRUE,
  spalte1=2, spalte2=5, epoche=30, fps=15, maxlag=3, automatic=FALSE,
  plotid=3, plotid2=5, pseudosVereinfachen=FALSE, anzahlPseudosGesamt=500,
  zufallsdaten=FALSE, anzahlzufallsdaten=10000, korrOhneBetragAnzeigen=TRUE,
  permutation=FALSE
)
}
\arguments{
  \item{filename}{ Character scalar, name of a source file. Defaults to \code{"MEA-topstage-C3-Beethoven.txt"}. }
  \item{separator}{ Character scalar, source file data separator. Defaults to \code{" "}. }
  \item{missings}{ Character scalar, source file data missings marker. Defalts to \code{"."}. }
  \item{datahead}{ Logical scalar, tells if source data has a header row. Defaults to \code{TRUE}. }
  \item{spalte1}{ Numeric, first column number of data to use. Defaults to \code{2}. }
  \item{spalte2}{ Numeric, second column number of data to use. Defaults to \code{5}. }
  \item{epoche}{ Epochengroesse (SEGMENT SIZE) in s (Epoche darf nicht kleiner sein als Zeitfenster/fps). Epoche*fps darf nicht groesser sein als die Anzahl Daten (timepoints) geteilt durch 2. Beispiel 1: timepoints: 1000, fps = 10, also darf Epoche nicht groesser sein als 50s (=1000/(10*2)). Beispiel 2: maxlag : 90, fps = 10, also darf Epoche nicht KLEINER sein als 9s (=90/10). }
  \item{fps}{ Numeric, frames per second, sampling rate. }
  \item{maxlag}{ Numeric, window size for CCF, maximum lag in seconds, only whole number. }
  \item{automatic}{ Logical, if \code{TRUE} then all pairs of columns are computed, therefore number of columns must be even and also a plot is not generated. }
  \item{plotid}{ Numeric, specifies the type of the first plot, defaults to \code{3}. \code{1}: \emph{GMcrosscorrs}, \code{2}: \emph{synchrony by epoch}, \code{3}: \emph{GM-Z}, \code{4}: \emph{time series plot}, \code{5}: \emph{z not abs}. }
  \item{plotid2}{ Numeric, specifies the type of the first plot, defaults to \code{5}. }
  \item{pseudosVereinfachen}{ Logical, default \code{FALSE}. Simplify calculation of the pseudo data. }
  \item{anzahlPseudosGesamt}{ Numeric, default \code{500}. Ignored when \code{pseudosVereinfachen} is \code{FALSE} (default) }
  \item{zufallsdaten}{ Logical, if \code{TRUE} then use random data, normally distributed numbers in the range [0,50]. }
  \item{anzahlzufallsdaten}{ Numeric, number of rows in random data. }
  \item{korrOhneBetragAnzeigen}{ Logical, defaults to \code{TRUE} display correlation without the amount. }
  \item{permutation}{ Logical, if \code{TRUE} then function computes all pairs combinations between columns provided in \code{x} (\code{n*(n-1)/2} pairs) }
}
\note{
  Function does not handle missing values.
}
\value{
  Returns \code{NULL} invisibly. A side effects include:
  \itemize{
    \item{ output to console, its form depends on \code{korrOhneBetragAnzeigen} argument. }
    \item{ plots, according to \code{plotid1} and \code{plotid2} arguments. }
    \item{ file \code{"CrossCorrelations.txt"} written to current working directory, if \code{plotid1} or \code{plotid2} specifies plot type \code{1}. }
  }
}
\seealso{\code{\link{susy}}}
\keyword{ internal }
