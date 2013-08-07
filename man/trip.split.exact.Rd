\name{trip.split.exact}
\alias{as.ltraj.trip}
\alias{as.SpatialLinesDataFrame.trip}
\alias{as.trip.SpatialLinesDataFrame}
\alias{trip-deprecated}
\alias{trip.split.exact}
\alias{tripTransform}
\title{Deprecated functions in trip}
\usage{
  as.SpatialLinesDataFrame.trip(from)

  trip.split.exact(x, dates)

  as.ltraj.trip(xy)

  as.trip.SpatialLinesDataFrame(from)

  tripTransform(x, crs, ...)
}
\arguments{
  \item{x}{see \code{\link{cut.trip}}}

  \item{dates}{see \code{\link{cut.trip}}}

  \item{from}{trip object}

  \item{crs}{CRS object, or PROJ.4 string accepted by
  \code{\link[sp]{CRS}}}

  \item{\dots}{Further arguments to
  \code{\link[rgdal]{spTransform}}}

  \item{xy}{\code{trip} object}
}
\description{
  These functions will be declared defunct in a future
  release.
}
\seealso{
  \code{\link{cut.trip}}, \code{\link{as.Other}}
}

