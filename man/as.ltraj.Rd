
\name{as.ltraj.trip}
\alias{as.ltraj.trip}
\alias{ltraj2trip}
\alias{coerce,trip,ltraj-method}
\alias{coerce,ltraj,trip-method}
\title{Coercion between trip objects and ltraj objects}

\description{coercion between classes}

 \section{Methods}{
   \describe{
     \item{coerce}{\code{signature(from = "trip", to = "ltraj")}}
     \item{coerce}{\code{signature(from = "ltraj", to = "trip")}}
     }
   }

   \usage{
as.ltraj.trip(xy, typeII = TRUE, slsp = "remove")
ltraj2trip(ltr)
}

\arguments{
  \item{xy}{trip object}
 % \item{date}{see \code{\link[adehabitat]{as.ltraj}}}
 % \item{id}{see \code{\link[adehabitat]{as.ltraj}}}
 % \item{burst}{see \code{\link[adehabitat]{as.ltraj}}}
  \item{typeII}{see \code{\link[adehabitat]{as.ltraj}}}
  \item{slsp}{details for the \code{\link[adehabitat]{ltraj}} turning angles}
  \item{ltr}{ltraj object}
}

\author{Michael D. Sumner}
\examples{
d <- data.frame(x = 1:10, y = rnorm(10), tms = Sys.time() + 1:10, id = gl(2, 5))

coordinates(d) <- ~x+y

tr <- trip(d, c("tms", "id"))
require(adehabitat)
l <- as.ltraj.trip(tr)

ltraj2trip(l)


}
\keyword{spatial}
\keyword{manip}