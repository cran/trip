
\name{as.ppp.trip}
\alias{as.ppp.trip}
\alias{as.psp.trip}
\alias{as.ppp}
\alias{as.psp}
\alias{coerce,trip,psp-method}
\alias{coerce,trip,ppp-method}
\title{Coercion between trip objects and spatstat objects}

\description{coercion between classes}

 \section{Methods}{
   \describe{
     \item{coerce}{\code{signature(from = "trip", to = "ppp")}}
     \item{coerce}{\code{signature(from = "trip", to = "psp")}}
     }
   }

\usage{
\method{as.ppp}{trip}(X, ..., fatal)
\method{as.psp}{trip}(x, ..., from, to)
}

\arguments{
  \item{X,x}{trip object}
  \item{...}{Ignored}
  \item{fatal}{Logical value, see Details of as.ppp}
  \item{from,to}{See as.psp}
}

\author{Michael D. Sumner}
\examples{
d <- data.frame(x = 1:10, y = rnorm(10), tms = Sys.time() + 1:10, id = gl(2, 5))

coordinates(d) <- ~x+y

tr <- trip(d, c("tms", "id"))

require(spatstat)

as.ppp(tr)

as.psp(tr)

}
\keyword{spatial}
\keyword{manip}