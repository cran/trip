\name{image2Grid}
\alias{image2Grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create a Grid from an image}
\description{
  A simple function to convert an image (xyz list) into a 
  SpatialPointsDataFrame. 
}
\usage{
image2Grid(x, p4 = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an list with components "x" and "y" (cor and "z" like that used by \code{image}
  and \code{contour}}
  \item{p4}{a PROJ.4 string specifying the coordinate system used}
}
\details{
 The coordinate vectors "x" and "y" may have the same length as the dimensions of
 matrix "z"
}
\value{
 A \code{SpatialGridDataFrame} with one column, "z". 
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
     # taken from image example
     x <- 10*(1:nrow(volcano))
     y <- 10*(1:ncol(volcano))
     im <- list(x = x, y = y, z = volcano)
     g <- image2Grid(im)
     image(im, col = terrain.colors(100), axes = FALSE)
     contour(g, add = TRUE)
     
     ## just in case there is any doubt
     image(g, xlim = c(100, 140), ylim = c(100, 140))
     abline(h = im$y+5, v = im$x+5)

}
\keyword{manip}% at least one, from doc/KEYWORDS
\keyword{list}
