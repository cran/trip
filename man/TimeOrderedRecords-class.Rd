\name{TimeOrderedRecords-class}
\docType{class}
\alias{TimeOrderedRecords-class}
\alias{trip,ANY,TimeOrderedRecords-method}
\alias{trip,trip,TimeOrderedRecords-method}
\title{Class "TimeOrderedRecords"}
\description{A simple class to act as a place-holder for DateTime and ID records in spatial data}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("TimeOrderedRecords", TOR.columns)}.
  \code{TOR.columns } are a 2-element character vector specifying the DateTime and ID columns in an
  object of class \code{trip}.
}
\section{Slots}{
	 \describe{
    \item{\code{TOR.columns}:}{2-element vector of class \code{"character"}}
  }
}
\section{Methods}{
  \describe{
    \item{trip}{\code{signature(obj = "ANY", TORnames = "TimeOrderedRecords")}: create a \code{trip} 
      object from a data frame }
    \item{trip}{\code{signature(obj = "trip", TORnames = "TimeOrderedRecords")}: create a \code{trip}
      object from an existing \code{trip} object }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{Michael D. Sumner}
\note{Future versions may change significantly, this class is very basic and could probably 
be implemented in a better way. Specifying TOR columns by formula would be a useful addition. }


\seealso{
	See Also \code{\link{trip}} for creating trip objects, and 
	 \code{\link{trip-class}} for the class. 
}
\examples{
tor <- new("TimeOrderedRecords", TOR.columns = c("datetime", "ID"))
tor <- TimeOrderedRecords(c("datetime", "ID"))
}
\keyword{classes}
