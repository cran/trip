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

\value{
A data frame with 8 columns. 
  \item{lon1, lat1 }{first pair of coordinates}
  \item{lon2, lat2}{second pair of coordinates}
  \item{gmt}{DateTimes as POSIXct}
  \item{id}{Platform Transmitting Terminal (PTT) ID}
  \item{lq}{Argos location quality class}
  \item{iq}{some other thing}

}

\keyword{IO}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
