\name{trackAngle}
\alias{trackAngle}
\alias{trackAngle.default}
\alias{trackAngle.trip}
\title{Determine distances or angles along a track}
\usage{
  trackAngle(x)
}
\arguments{
  \item{x}{trip object, or matrix of 2-columns, with x/y
  coordinates}
}
\value{
  Vector of angles (degrees) between coordinates.
}
\description{
  Calculate the angles between subsequent 2-D coordinates
  using Great Circle distance (spherical) methods.
}
\details{
  If \code{x} is a trip object, the return result has an
  extra element for the start and end point of each
  individual trip, with value NA.

  This is an optimized hybrid of "raster::bearing" and
  \code{\link[maptools]{gzAzimuth}}.
}
\author{
  Robert Hijmans, Roger Bivand
}

