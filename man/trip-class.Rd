\name{trip-class}
\docType{class}
\alias{trip-class}
%\alias{[,trip-method}
\alias{[,trip,ANY,ANY,ANY-method}
\alias{lines,trip-method}
\alias{names.trip}
\alias{names<-.trip}
\alias{$<-,trip,character}
\alias{$<-,trip,character-method}
\alias{[[<-,trip,ANY,missing-method}
\alias{plot,trip,missing-method}
\alias{points,trip-method}
\alias{recenter,trip-method}
\alias{show,trip-method}
%\alias{spTransform,trip,CRS-method}
\alias{summary,trip-method}
\alias{text,trip-method}
\alias{subset,trip-method}
\alias{trip,trip,ANY-method}
\alias{trip,SpatialPointsDataFrame,ANY-method}
\alias{trip,ANY,ANY}
%\alias{trip,ANY,ANY,ANY}
\alias{trip,trip-method}
%\alias{$.trip,$<-.trip,[[.trip,[[<-.trip}
\alias{as.data.frame.trip}
\title{Class "trip" ~~~ }
\description{An extension of \code{"SpatialPointsDataFrame"} by including \code{"TimeOrderedRecords"}.
The records within the data frame are explicitly ordered by DateTime data within IDs.}
\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{trip(obj = "SpatialPointsDataFrame", TORnames = "TimeOrderedRecords"}.
  The object contains all the slots present within a \code{SpatialPointsDataFrame}, particularly \code{data}
  which contains columns of at least those specified by \code{TOR.columns}
}
\section{Slots}{
	 \describe{
    \item{\code{TOR.columns}:}{Object of class \code{"character"} specifiying the DateTime and ID columns (in that order) 
    in \code{data} }
    \item{\code{data}:}{Object of class \code{"data.frame"} the native data object for a Spatial data frame}
    
    }
    Also,  other slots usual to a \code{SpatialPointsDataFrame}
}
\section{Extends}{
Class \code{"TimeOrderedRecords"}, directly.
Class \code{"SpatialPointsDataFrame"}, directly.
Class \code{"SpatialPoints"}, by class \code{"SpatialPointsDataFrame"}.
Class \code{"Spatial"}, by class \code{"SpatialPointsDataFrame"}.
}
\section{Methods}{ Most of the methods available are by virtue of the sp package.
 Some, such as \code{split.data.frame} have been added to SPDF so that trip has the
 same functionality. 
  \describe{
     \item{trip}{\code{signature("ANY")}: try to return a trip object}
      \item{trip}{\code{signature(obj = "SpatialPointsDataFrame", TORnames = "ANY"}}{The usual construction}

    \item{[}{\code{signature(x = "trip")}: subset rows or columns as per 
    \code{SpatialPointsDataFrame}. 
  In the case that \code{TimeOrderedRecords} columns are dropped, the object 
  reverts to the straight Spatial version. }
    \item{lines}{\code{signature(x = "trip")}: add lines to a plot with separate 
  colours for each trip}
   %% \item{names}{\code{signature(x = "trip")}: return the column names of the data}
    \item{plot}{\code{signature(x = "trip", y = "missing")}: plot as 
  \code{SpatialPoints} }
    \item{points}{\code{signature(x = "trip")}: add points to a plot using the 
  Spatial coordinates }
    \item{recenter}{\code{signature(obj = "trip")}: perform coordinate recentering, 
  from the [-180,180] convention to [0, 360 }
    \item{show}{\code{signature(object = "trip")}: print a short summary of the trip data}
   % \item{spTransform}{\code{signature(x = "trip", "CRS")}: transform coordinates via rgal PROJ.4 }
    \item{summary}{\code{signature(object = "trip")}: print a summary as per 
  \code{SpatialPointsDataFrame} including a summary of the trip data}
    \item{text}{\code{signature(x = "trip")}: add text to a plot using Spatial coordinates}
    \item{trip}{\code{signature(obj = "trip", TORnames = "TimeOrderedRecords")}: (re)-create a trip object using \code{TimeOrderedRecords}}
    \item{trip}{\code{signature(obj = "trip", TORnames = "ANY")}: (re)-create a trip object by some other means }
    
%   \item{subset}{\code{signature(x = "trip", subset, select, drop = FALSE, ...)}: subset a trip in the expected manner. 
   \item{subset}{\code{signature(x = "trip",  ...)}: subset a trip in the expected manner. 
    %\item{split}{\code{signature(x = "trip")}: split a trip like a data frame}
   }
}}

\author{Michael D. Sumner }


\section{Warning }{There are some kludges to allow \code{trip} to do things, such as replace \code{POSIXt}
column data using \code{"$<-.trip"} and \code{"[[<-.trip"} which should not be
necessary once sp implements the new \code{data.frame} class of R >=2.4.0.} 

\seealso{
 \code{\link{trip}} for examples of directly using the class. 
}

\keyword{classes}
