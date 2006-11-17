\name{readDiag}
\alias{readDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read Argos DIAG format}
\description{
  Create a dataframe from Argos diag files. 
}
\usage{
readDiag(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{one or more names of DIAG files }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame with 8 columns. 
  \item{lon1, lat1 }{first pair of coordinates}
  \item{lon2, lat2}{second pair of coordinates}
  \item{gmt}{DateTimes as POSIXct}
  \item{id}{Platform Transmitting Terminal (PTT) ID}
  \item{lq}{Argos location quality class}
  \item{iq}{some other thing}

}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 \section{Warning }{There are other data in these files not returned by this function} 
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
\keyword{IO}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
