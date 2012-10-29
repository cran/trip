\name{trip}
\alias{trip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to handle animal track data, organized as "trip"s}
\description{
  Extend the basic functionality of a \code{Spatial} data frame by specifying the
  data columns that define the "TimeOrdered" quality of the records.
}
\usage{
trip(obj, TORnames)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A \code{SpatialPointsDataFrame}, containing at least two columns
    with the DateTime and ID data as per \code{TORnames}}
  \item{TORnames}{Either an object of \code{TimeOrderedRecords}, or a 2-element
  character vector specifying the DateTime and ID column of \code{obj}}
}

\value{
 A trip object, with the usual slots of a \code{SpatialPointsDataFrame} and the
 added \code{TimeOrderedRecords}. For the most part this can be treated as a
 \code{data.frame} with \code{Spatial} coordinates.
}

\author{Michael D. Sumner}

\seealso{ \code{\link{speedfilter}}, and  \code{\link{tripGrid}} for simple(istic)
speed filtering and spatial time spent gridding.
}
\examples{
d <- data.frame(x = 1:10, y = rnorm(10), tms = Sys.time() + 1:10, id = gl(2, 5))
coordinates(d) <- ~x+y
tr <- trip(d, c("tms", "id"))

\dontrun{

## a simple example with the common fixes required for basic track data

dat <- read.csv("trackfile.csv")
names(dat)  ## e.g. [1] "long" "lat"  "seal"  "date"    "local"     "lq"
library(sp)
coordinates(dat) <- c("long", "lat")

## date/times may be in a particular time zone, please check
dat$gmt <- as.POSIXct(strptime(paste(dat$date, dat$local),
                      "%d-%b-%y %H:%M:%S"), tz = "GMT")


## if there are problems in the data, this will error
tr <- trip(dat, c("gmt", "seal"))


## the following code tries to fix common problems

## remove completely-duplicated rows
dat <- dat[!duplicated(dat), ]

## order the rows by seal, then by time
dat <- dat[order(dat$seal, dat$gmt), ]

## fudge duplicated times
dat$gmt <- adjust.duplicateTimes(dat$gmt, dat$seal)


## finally, convert to Spatial and create trip object
coordinates(dat) <- c("long", "lat")
tr <- trip(dat, c("gmt", "seal"))



}

%% data examples with adding columns, subsetting etc.
}
\keyword{manip}% at least one, from doc/KEYWORDS

