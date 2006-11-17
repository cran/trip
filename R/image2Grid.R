image2Grid <- function(x, p4 = NA) {

	if (!all(names(x) %in% c("x", "y", "z"))) stop("image must have components x, y, and z")
	
	cells.dim <- dim(x$z)
	xx <- x$x
	yy <- x$y
	lx <- length(xx)
	ly <- length(yy)
	if (all(c(lx, ly) == (cells.dim + 1))) {
	print("corners")
		xx <- xx[-1] - diff(xx[1:2])/2
		yy <- yy[-1] - diff(yy[1:2])/2
	}
	SpatialGridDataFrame(GridTopology(c(xx[1], yy[1]), c(diff(xx[1:2]), diff(yy[1:2])),
		cells.dim), data.frame(z = as.vector(x$z[,ncol(x$z):1])), proj4string = CRS(p4))
}