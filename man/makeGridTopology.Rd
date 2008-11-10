\name{makeGridTopology}
\alias{makeGridTopology}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a GridTopology from a Spatial object}
\description{
  Sensible defaults are assumed, to match the extents of data to a manageable
  grid. Approximations for kilometres in longlat can be made using \code{cellsize}. 
}
\usage{
makeGridTopology(obj, cells.dim = c(100, 100), xlim = NULL, ylim = NULL, buffer = 0.1, cellsize = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{ ~~Describe \code{obj} here~~ }
  \item{cells.dim}{ ~~Describe \code{cells.dim} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{ylim}{ ~~Describe \code{ylim} here~~ }
  \item{buffer}{ ~~Describe \code{buffer} here~~ }
  \item{cellsize}{ ~~Describe \code{cellsize} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

}
\keyword{manip}% at least one, from doc/KEYWORDS

