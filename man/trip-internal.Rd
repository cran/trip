\docType{data}
\name{trip-internal}
\alias{.abdali}
\alias{.distances}
\alias{.g2ow}
\alias{.gcdist.c}
\alias{.intpFun}
\alias{.oc.col}
\alias{.single.trip.split}
\alias{.tripRbind}
\alias{.validTORdata}
\alias{trip-internal}
\title{Internal trip Functions}
\format{chr [1:256] "#FFFFFF" "#A00AFF" "#A013FF" "#A01DFF" "#A027FF" ...}
\usage{
  .validTORdata(object)

  .intpFun(x)

  .distances(x)

  .gcdist.c(lon1, lat1, lon2, lat2)

  .abdali(p1, p2)

  .tripRbind(obj, x)

  .single.trip.split(tr1, boundary.dates)

  .oc.col

  .g2ow(x)
}
\description{
  Internal trip functions
}
\details{
  These are not to be called by the user (or in some cases
  are just waiting for proper documentation to be written).
}
\keyword{datasets}
\keyword{internal}

