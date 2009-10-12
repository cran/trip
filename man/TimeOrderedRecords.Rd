\name{TimeOrderedRecords}
\alias{TimeOrderedRecords}
\alias{getTORnames}
\alias{getTimeID}
\title{Functions to specify and obtain DateTime and ID data from within (Spatial) data frames.}
\description{
  Specifies the DateTime and ID data within a data frame for objects of class \code{trip}.  Functions
  also for obtaining the names of the columns used, and the data itself. 
}
\usage{
TimeOrderedRecords(x)
getTORnames(obj)
getTimeID(obj)
}
\arguments{
  \item{x}{Character vector of 2-elements, specifying the data columns of DateTimes and IDs.}
  \item{obj}{\code{trip} object.}
}
\details{
  These simple functions are for creating classes with TimeOrdered data. The main use is for
  SpatialPointsDataFrames which are used with TimeOrderedRecords for the class \code{trip}.
}
\value{
  \code{TimeOrderedRecords} returns an object with a 2-element character vector, specifying the columns names. 
  \code{getTORnames} obtains the column names from an object extending the class, and \code{getTimeID} returns the
  data as a data frame from an object extending the class. 
}

\seealso{\code{\link{trip}}, for the use of this class with \code{\link[sp]{SpatialPointsDataFrame}}}
\examples{

tor <- TimeOrderedRecords(c("time", "id"))
getTORnames(tor)

}
\keyword{manip }% at least one, from doc/KEYWORDS

