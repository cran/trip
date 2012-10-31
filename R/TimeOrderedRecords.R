
setClass("TimeOrderedRecords", representation(TOR.columns = "character"))

setValidity("TimeOrderedRecords", function(object) {

	if (!is.character(object@TOR.columns) | !is.vector(object@TOR.columns))
		stop("TimeOrderedRecords data names must be character vector")
	## also support length == 1?
	if (length(object@TOR.columns) > 2) stop("TimeOrderedRecords data names must be of length 2")
	TRUE
})

TimeOrderedRecords <- function(x) {
	new("TimeOrderedRecords", TOR.columns = x)
}

getTORnames <- function(obj) obj@TOR.columns

getTimeID <- function(obj) as.data.frame(obj)[, getTORnames(obj)]

