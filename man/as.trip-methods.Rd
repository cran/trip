\docType{methods}
\name{as.trip-methods}
\alias{as.trip}
\alias{as.trip-methods}
\alias{as.trip,ltraj-method}
\alias{coerce,trip,ltraj-method}
\alias{ltraj2trip}
\title{Coercion from other classes to \code{trip} objects}
\arguments{
  \item{x,ltr}{ltraj object}

  \item{list()}{Arguments passed to other methods. Ignored
  for \code{ltraj} method.}
}
\description{
  Coercing objects to \code{trip} class
}
\section{Methods}{
  \describe{

  \item{coerce}{\code{signature(from="ltraj", to="trip")}}

  \item{as.trip}{\code{signature(x="ltraj")}}

  }
}
\examples{
## Continuing the example from '?trip-methods:
utils::example("trip-methods", package="trip",
               ask=FALSE, echo=FALSE)

if (require(adehabitatLT)) {
    ##l <- as.ltraj.trip(tr)
    ##ltraj2trip(l)
    ##as.trip(l)
}
}

