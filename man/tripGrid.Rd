\name{tripGrid}
\alias{tripGrid}
\alias{interpequal}
\alias{countPoints}
\alias{kdePoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a grid of time spent}
\description{
  Create a grid of time spent from an object of class \code{trip} by approximating
  the time between locations for separate trip events. 
}
\usage{
tripGrid(x, grid = NULL, method = "count", dur = NULL,  ...)
interpequal(x, dur = NULL, quiet = FALSE)
countPoints(x, dur = 1, grid = NULL) 
kdePoints(x, h = NULL, grid =NULL, resetTime = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class trip}
  \item{grid}{GridTopology - will be generated automatically if NULL}
  \item{method}{name of method for quantifying time spent, see Details }
  \item{\dots}{other arguments passed to \code{interpequal} 
  or \code{kdePoints}}
  
  
  %% interpequal
  \item{dur}{The \"dur\"ation of time used to interpolate between available 
  locations (see Details)}
  \item{quiet}{logical - report on difference between time summed and time in trip?}
  
  %% count/kdePoints

  \item{h}{numeric vector of two elements specifying bandwidth for kernel density}
  \item{resetTime}{logical - reset the values of the kde grid to match the sum
  of the total time?}
  }
\details{
The intention is for \code{tripGrid} to be used, rather than the underlying functions
for interpolation or gridding. 

Trip locations are first interpolated, based on an equal-time spacing between
records. These interpolated points are then "binned" to a grid of cells.  The
time spacing is specified by the "dur"ation argument to \code{interpequal} in
seconds (i.e. \code{dur = 3600} is used for 1 hour). Shorter time periods will
require longer computation with a closer approximation to the total time spent in
the gridded result. 

Currently there are methods "count" and "kde" for quantifying time spent, corresponding
to the functions "countPoints" and "kdePoints". "kde" uses
kernel density to smooth the locations, "count" simply counts the points falling
in a grid cell. 
}
\value{
 \code{tripGrid} returns an object of class \code{SpatialGridDataFrame}, with one
 column "z" containing the time spent in each cell (usually assumed to be in hours).
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

