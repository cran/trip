
\name{as.trip.SpatialLinesDataFrame}
\alias{as.trip.SpatialLinesDataFrame}
\alias{coerce,trip,SpatialLinesDataFrame-method}

\title{Coercion between trip objects and sp line objects}

\description{coercion between classes}

 \section{Methods}{
   \describe{
     \item{coerce}{\code{signature(from = "trip", to = "SpatialLinesDataFrame")}}
     }
   }

\usage{
as.trip.SpatialLinesDataFrame(from)
%%as(from, "SpatialLinesDataFrame")
}

\arguments{
  \item{from}{trip object}

}

\author{Michael D. Sumner}
\examples{
d <- data.frame(x = 1:10, y = rnorm(10), tms = Sys.time() + 1:10, id = gl(2, 5))

coordinates(d) <- ~x+y

tr <- trip(d, c("tms", "id"))

as.trip.SpatialLinesDataFrame(tr)

as(tr, "SpatialLinesDataFrame")

}
\keyword{spatial}
\keyword{manip}