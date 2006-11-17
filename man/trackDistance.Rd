\name{trackDistance}
\alias{trackDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Determine distance along a track}
\description{
  Calculate the distance between subsequent 2-D coordinates.
}
\usage{
trackDistance(track, longlat = FALSE, push = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{track}{matrix of 2-columns, with x/y coordinates }
  \item{longlat}{logical - are the data longlat?}
  \item{push}{an offset for measuring distance between this and the next (push) 
 location}
}
\details{
  This uses the (fast) \code{spDistsN1}. 
}
\value{
Vector of distances between coordinates. 
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ \code{\link[sp]{spDistsN1}} for the underlying code }
\examples{

}
\keyword{manip}% at least one, from doc/KEYWORDS

