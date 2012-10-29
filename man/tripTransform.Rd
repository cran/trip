\name{tripTransform}
\alias{tripTransform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reproject trip objects. 
}
\description{
  Projection transformation based on CRS strings in rgdal
}
\usage{
tripTransform(x, crs, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 trip object with projection metadata
}
  \item{crs}{
 CRS object, or PROJ.4 string accepted by \code{\link[sp]{CRS}}
}
  \item{\dots}{
 Further arguments to \code{\link[rgdal]{spTransform}}
}

}
\value{
trip object, with spatial coordinates reprojected
}

\author{
Michael D. Sumner
}
\note{
this is basically a cheat as I don't want to Depends on rgdal for trip, but
cannot figure out the proper way to define things such that Suggests is
enough for spTransform methods on trips
}

\seealso{
\code{\link[rgdal]{spTransform}}
}
\examples{
d <- data.frame(x = 1:10, y = rnorm(10), tms = Sys.time() + 1:10, id = gl(2, 5))
coordinates(d) <- ~x+y

tr <- trip(d, c("tms", "id"))
proj4string(tr) <- CRS("+proj=laea +lon_0=146")

tripTransform(tr, "+proj=longlat")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }

