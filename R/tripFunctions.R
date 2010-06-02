
 # readArgos

	 # adjust.duplicateTimes

	 # argos.sigma

	 # sepIdGaps

 # tripGrid
 # countPoints
 # kdePoints

 # interpequal

 # speedfilter

 # makeGridTopology

 # trackDistance


## need to clean up the "internal" functions, and ensure the arguments are passed in correctly
##  - and figure out which arguments are really useful anyway

forceCompliance <- function(x, tor) {
	isSpatial <- is(x, "SpatialPointsDataFrame")
	if (isSpatial) {

		crd.nrs <- x@coords.nrs
		x <- as.data.frame(x)
	}



	levs <- unique(x[[tor[2]]])
	tooshort <- tapply(x[[1]], x[[tor[2]]], function(x) length(x) < 3)
	x <- x[x[[tor[2]]] %in% levs[!tooshort], ]




	x <- x[!duplicated(x), ]
	x <- x[order(x[[tor[2]]], x[[tor[1]]]), ]
	x[[tor[1]]] <- adjust.duplicateTimes(x[[tor[1]]], x[[tor[2]]])
	if (isSpatial) {
		coordinates(x) <- crd.nrs
		x <- trip(x, tor)
	}
	x
}

interpequal <- function(x, dur = NULL, quiet = FALSE) {

  #equalTime<- function (coords, time, id = factor(rep(1, length(x))), dur = NULL)
#	{

if (!is(x, "trip")) stop("only trip objects supported")
if (is.null(dur)) stop("equal time duration must be specified \"dur = ?\"")
	## x must be a single trip
	#intp <- equalTime(coordinates(x), time = x[[x@TOR.columns[1]]],
	#					id = x[[x@TOR.columns[2]]],
	#dur = dur)

	tor <- x@TOR.columns
	time <- x[[tor[1]]]
	id <- factor(x[[tor[2]]])

	coords <- coordinates(x)
	x <- coords[,1]
	y <- coords[,2]

	levs <- levels(id)
	newPts <- NULL
	#if (is.null(dur))
	#   dur <- as.numeric(min(unlist(tapply(as.integer(time),
	#            id, diff))))
	    intpFun <- function(x) {
	        len <- round(x[3] + 1)
	        new <- seq(x[1], x[2], length = len)
	        if (len > 1)
	            new[-len]
	        else new
	    }
	    for (sub in levs) {
	        ind <- id == sub
	        xx <- x[ind]
	        yy <- y[ind]
	        tms <- time[ind]
	        dt <- diff(as.numeric(tms))
	        dtn <- dt/dur
	        ax <- cbind(xx, c(xx[-1], xx[length(xx)]), c(dtn, 0))
	        ay <- cbind(yy, c(yy[-1], yy[length(yy)]), c(dtn, 0))
	        intime <- as.numeric(tms) - min(as.numeric(tms))
	        at <- cbind(intime, c(intime[-1], intime[length(intime)]),
	            c(dtn, 0))
	        nx <- unlist(apply(ax, 1, intpFun))
	        ny <- unlist(apply(ay, 1, intpFun))
	        nt <- unlist(apply(at, 1, intpFun)) + min(tms)
	        ni <- factor(rep(sub, length = length(nt)))
	        newPts <- rbind(newPts, data.frame(x = nx, y = ny, time = nt,
	            id = ni))
	    }
	    origTotal <- sum(tapply(time, id, function(x) diff(range(as.numeric(x)))))
	    newTotal <- nrow(newPts) * dur
	    uType = "hours"
	    hTotal <- sum(tapply(time, id, function(x) difftime(range(x)[2],
	        range(x)[1], units = uType)))
	    if (!quiet) {
	    cat("lost seconds = ", as.integer(origTotal - newTotal),
	        " out of a total ", hTotal, " ", uType, "\n")
	    }
#}

	coordinates(newPts) <- c("x", "y")
	names(newPts) <- tor
	newPts
}

tripGrid.interp <- function(x, grid = NULL, method = "count", dur = NULL, ...) {

	method <- paste(method, "Points", sep = "")
	if (!exists(method)) stop("no such method: ", method)
	cat("Using method ", method, "\n\n")

	if (is.null(grid)) grid <- makeGridTopology(x)
	res <- SpatialGridDataFrame(grid, data.frame(z = rep(0, prod(grid@cells.dim))), CRS(proj4string(x)))

	tor <- x@TOR.columns
       	trip.list <- split(x[, tor], x[[tor[2]]])

	cnt <- 0
	for (this in trip.list) {
		this <- interpequal(this, dur = dur, quiet = TRUE)
		cnt <- cnt + nrow(this)
		res$z <- res$z + do.call(method, list(x = trip(this, tor), grid = grid, ...))$z
	}

	#origTotal <- sum(tapply(x[[x@TOR.columns[1]]], x[[x@TOR.columns[2]]], function(x) diff(range(as.numeric(x)))))
	#newTotal <- sum(res$z) * 3600
	#uType = "hours"
	#hTotal <- sum(tapply(x[[x@TOR.columns[1]]], x[[x@TOR.columns[2]]], function(x) difftime(range(x)[2],
	#	        range(x)[1], units = uType)))

	#cat("lost seconds = ", as.integer(origTotal - newTotal),
	#	       " out of a total ", hTotal, " ", uType, "\n")


	if (method == "countPoints") res$z <- res$z * dur
	#if (hours) res$z <- res$z/3600
	res
}

## TODO:
##a version of tripGrid that takes Lines, so as.SpatialLinesDataFrame.trip, and then to grid

## allow sigmas argument for density version, for each line segment


## replaces tripGrid, old version is now called tripGrid.interp

tripGrid <-
function (x, grid = NULL, method = "pixellate", ...)
{
    if (method %in% c("kde", "count"))
        warning("kde and count methods no longer supported from trip_1.1-6 and will be ignored, see ?tripGrid.interp for legacy function")
    if (!is.null(list(...)$dur))
        stop("dur(ation) not necessary for this function from trip_1.1-6 and will be ignored - time sum is now exact\n see ?tripGrid.interp 
for legacy function")
    require(spatstat)
  g2ow <- function(x) {
        mn <- x@cellcentre.offset - x@cellsize/2
        mx <- mn + x@cells.dim * x@cellsize
        owin(c(mn[1], mx[1]), c(mn[2], mx[2]), mask = matrix(TRUE, x@cells.dim[2], x@cells.dim[1]),
	xy =list(x = seq(mn[1], mx[1], length = x@cells.dim[1]), y = seq(mn[2], mx[2], length = x@cells.dim[2])) )

#owin(c(mn[2], mx[2]), c(mn[1], mx[1]), mask = matrix(TRUE, x@cells.dim[2], x@cells.dim[1]),
#	xy =list(x = seq(mn[1], mx[1], length = x@cells.dim[1]), y = seq(mn[2], mx[2], length = x@cells.dim[2])) )
    }
    if (is.null(grid))
        grid <- makeGridTopology(x)
    res <- as.image.SpatialGridDataFrame(SpatialGridDataFrame(grid,
        data.frame(z = rep(0, prod(grid@cells.dim)))))
    tor <- x@TOR.columns
    trip.list <- split.data.frame(x[, tor], x[[tor[2]]])
    ow <- g2ow(grid)
    sm <- 0
    zero.lengths <- FALSE
    sz <- 0
    for (this in trip.list) {
        xs <- coordinates(this)[, 1]
        ys <- coordinates(this)[, 2]
        dt <- diff(unclass(this[[tor[1]]]))
        sm <- sm + sum(dt)
        x.psp <- psp(xs[-length(xs)], ys[-length(ys)], xs[-1],
            ys[-1], window = ow)
        lngths <- lengths.psp(x.psp)


        if (any(!lngths > 0)) {
            zero.lengths <- TRUE
            zeros <- which(!lngths > 0)
            cc <- coordinates(this)[zeros, , drop = FALSE]
            x.ppp <- ppp(cc[, 1], cc[, 2], window = ow)
            if (method == "pixellate") {
                v <- pixellate(x.ppp, W = ow, weights = dt[zeros])$v
            }
            if (method == "density") {
                v <- density(x.ppp, ...)$v
            }

            res$z <- res$z + t(v)
            sz <- sz + sum(dt[zeros])
        }
        weights <- dt/ifelse(lngths > 0, lngths, .Machine$double.eps)
        if (method == "pixellate") {
            v <- pixellate(x.psp, W = ow, weights = weights)$v
        }
        if (method == "density") {
            v <- density(x.psp, ...)$v
        }
        res$z <- res$z + t(v)
    }
    if (zero.lengths) {
        warning("zero length lines present, time durations summed into cells assuming point-presence of degenerate line segment")
        cat("\n")
        cat(paste("Total time of trips:", sm, "\n"))
        cat(paste("Total time without zero length lines:", sm -
            sz, "\n"))
    }
    image2Grid(res, p4 = proj4string(x))
}




## tripGrid <- function(x, grid = NULL, method = "pixellate",...) {


##     ## deal with legacy
##     if (method %in% c("kde", "count")) warning("kde and count methods no longer supported from trip_1.1-6 and will be ignored, see ?tripGrid.interp for legacy function")

##     if (!is.null(list(...)$dur)) stop("dur(ation) not necessary for this function from trip_1.1-6 and will be ignored - time sum is now exact\n see ?tripGrid.interp for legacy function")
##     require(spatstat)
##     ## SpatialGridDataFrame to owin
##     g2ow <- function(x) {
##         mn <- x@cellcentre.offset - x@cellsize/2
##         mx <- mn + x@cells.dim * x@cellsize
##         owin(c(mn[1], mx[1]), c(mn[2], mx[2]), mask = matrix(TRUE, x@cells.dim[1], x@cells.dim[2]))
##     }
##     if (is.null(grid)) grid <- makeGridTopology(x)
##     res <- as.image.SpatialGridDataFrame(SpatialGridDataFrame(grid, data.frame(z = rep(0, prod(grid@cells.dim)))))
##     tor <- x@TOR.columns

##     trip.list <- split.data.frame(x[, tor], x[[tor[2]]])
##     ow <- g2ow(grid)
##     sm <- 0
##     zero.lengths <-  FALSE
##     sz <- 0
##     for (this in trip.list) {
##         xs <- coordinates(this)[,1]
##         ys <- coordinates(this)[,2]
##         dt <- diff(unclass(this[[tor[1]]]))
##         sm <- sm + sum(dt)

##         x.psp <- psp(xs[-length(xs)], ys[-length(ys)], xs[-1], ys [-1], window = ow)
##         lngths <- lengths.psp(x.psp)

##            if (any(!lngths > 0)) {
##                zero.lengths <- TRUE
##                zeros <- which(!lngths > 0)
##                cc <- coordinates(this)[zeros,,drop = FALSE]
##                #idx <- getGridIndex(cc, grid)
##                x.ppp <- ppp(cc[,1], cc[,2], window = ow)
##                if (method == "pixellate") {
##                    v <- pixellate(x.ppp, W = ow, weights = dt[zeros])$v
##                }
##                if (method == "density") {
##                    v <- density(x.ppp, ...)$v
##                }
##                res$z <- res$z + t(v)
##                sz <- sz + sum(dt[zeros])
##            }


##         ## line lengths may be zero
##         weights <- dt/ifelse(lngths > 0, lngths, .Machine$double.eps)
##         if (method == "pixellate") {
##             v <- pixellate(x.psp, W = ow, weights = weights)$v
##         }
##         if (method == "density") {
##             v <- density(x.psp,  ...)$v
##         }
##         res$z <- res$z + t(v)

##         #    x.psp <- psp(xs[i-1], ys[i-1], xs[i], ys[i], window = ow)
##         #    v <- dt[i-1] * pixellate(x.psp)$v / lengths.psp(x.psp)
##         #    res$z <- res$z + v
##         #}
##     }

##    # if (zero.lengths) {
##     #       if (zero.lengths) {
##      #          cat("\n")
## #
##  #           warning(paste("discrepancy in time sum is due to zero length line segments in trip", if (length(trip.list) > 1) "s", sep = ""))
##   #         }

##     if (zero.lengths) {
##         warning("zero length lines present, time durations summed into cells assuming point-presence of degenerate line segment")
##         cat("\n")
##         cat(paste("Total time of trips:", sm, "\n"))
##         cat(paste("Total time without zero length lines:", sm - sz, "\n"))

##     }
##     image2Grid(res, p4 = proj4string(x))
## }


kdePoints <- function (x, h = NULL, grid =NULL, resetTime = TRUE, ...)
{
    coords <- coordinates(x)
    xx <- coords[ , 1]
    yy <- coords[ , 2]
    time <- x[[getTORnames(x)[1]]]
    id <- x[[getTORnames(x)[2]]]
    timesum <- sum(tapply(time, id, function(x) diff(range(unclass(x)))))

    ## must acknowledge MASS for this
    if (missing(h)) {
        bandwidth.nrd <- function(x) {
            r <- quantile(x, c(0.25, 0.75))
            h <- (r[2] - r[1])/1.34
            4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
        }
        h <- c(bandwidth.nrd(xx), bandwidth.nrd(yy))/10
    }
    if (is.null(grid))  grid <- makeGridTopology(coords, ...)
    ## use bbox here
    dimXY <- grid@cells.dim
    nx <- nrow(x)
    gcs <- coordinatevalues(grid)
    gx <- gcs$s1 + grid@cellsize[1]
    gy <- gcs$s2 + grid@cellsize[2]
    ax <- outer(gx, xx, "-")/h[1]
    ay <- outer(gy, yy, "-")/h[2]
    z <- (matrix(dnorm(ax), dimXY[1], nx) %*% t(matrix(dnorm(ay),
        dimXY[2], nx)))/(nx * h[1] * h[2])
    if (resetTime) z <- (z * timesum/sum(z)) / 3600
    SpatialGridDataFrame(grid, data.frame(z = as.vector(z)), CRS(proj4string(x)))
}

countPoints <-
function (x, dur = 1, grid = NULL)
{
    coords <- coordinates(x)
    xx <- coords[ , 1]
    yy <- coords[ , 2]

    if (is.null(grid))  grid <- makeGridTopology(coords)

    #grd <- mkGrid(x, y, ...)
    #orig <- c(grd$x[1], grd$y[1])
    orig <- grid@cellcentre.offset - grid@cellsize/2

    #scl <- c(diff(grd$x)[1], diff(grd$y)[1])
    scl <- grid@cellsize
    xdim <- grid@cells.dim[1]
    ydim <- grid@cells.dim[2]
    xmin <- orig[1]
    xmax <- orig[1] + (xdim + 1) * scl[1]
    ymin <- orig[2]
    ymax <- orig[2] + (ydim + 1) * scl[2]
    xlim <- c(xmin, xmax)
    ylim <- c(ymin, ymax)

    if (xlim[1] < xmin || xlim[2] > (xmax) || ylim[1] <
        ymin || ylim[2] > (ymax))
        stop("Data are out of bounds")
    cps <- ceiling(cbind((xx - orig[1])/scl[1], (yy - orig[2])/scl[2]))
    tps <- tabulate((cps[, 1] - 1) * ydim + cps[, 2], xdim *
        ydim)
    mps <- matrix(tps, ydim, xdim)
    z <- t(mps)
    SpatialGridDataFrame(grid, data.frame(z = as.vector(z[,ncol(z):1])), CRS(proj4string(x)))
}


makeGridTopology <-
function (obj, cells.dim =c(100, 100),
    xlim =NULL, ylim = NULL, buffer = 0, cellsize = NULL, adjust2longlat = FALSE)
{
   if ((is.null(xlim) | is.null(ylim)) & missing(obj)) stop("require at least a Spatial object, matrix object, or xlim and ylim")
   if (!missing(obj)) bb <- bbox(obj)
   if (!is.null(xlim) & !is.null(ylim)) buffer <- 0
   if (is.null(xlim)) xlim <- bb[1,]
   if (is.null(ylim)) ylim <- bb[2,]
## PROBLEMS
   ## determination is boundary based, but grid is cell based
   ## break down into simpler functions, then recombine including longlat adjust
   ## gridFromNothing - world1 ?

   ## gridFromLimits
   ## gridFromLimits/dims
   ## gridFromLimits/cellsize
   ##
   ## gridFromDims?
   ## gridFromCellsize?
   ## gridFromDims/Cellsize?

   #proj <- NA
   #if (!missing(obj)) proj <- is.projected(obj)
   #if (is.na(proj)) {
  # 	warning("coordinate system unknown, assuming longlat")
  # 	proj <- FALSE
  # }
   if (is.null(cellsize) & adjust2longlat) warning("cellsize not provided with adjust2longlat, ignoring")
   if (!is.null(cellsize)) {
       if (!length(cellsize) == 2) stop("cellsize must be of length 2")
           if (adjust2longlat) {
               cellsize <- c(cellsize[1]/(cos((pi/180) * mean(ylim)) * 1.852 * 60), cellsize[2]/(1.852 * 60))
               if (any(!cellsize > 0)) stop("longlat adjustment resulted in invalid cellsize. Does it really make sense for these latitude limits? \n", paste(format(ylim), collapse = ","))
           }
		   xvalues <- seq(xlim[1], xlim[2] + cellsize[1], by = cellsize[1])
		   yvalues <- seq(ylim[1], ylim[2] + cellsize[2], by = cellsize[2])
		   xlim <- range(xvalues)
    		   ylim <- range(yvalues)
    		   cells.dim <- c(length(xvalues), length(yvalues))
	} else   cellsize <- c(diff(xlim), diff(ylim))/(cells.dim - 1)

    if (buffer > 0) {
        addXY <- ceiling(cellsize * buffer)
        xlim <- xlim + c(-addXY[1], addXY[1])
        ylim <- ylim + c(-addXY[2], addXY[2])
        cellsize <- c(diff(xlim), diff(ylim))/(cells.dim - 1)
    }
  new("GridTopology", cellcentre.offset = c(min(xlim), min(ylim)),
         cellsize = cellsize, cells.dim = as.integer(cells.dim))
}


## gridFromCellsize <- function(cell.size = c(1, 1)) {
##     GridTopology(c(180, -90) + cell.size/2, cell.size
## }
## gridFromDims <- function(cells.dim = c(360, 180)) {
##    GridTopology(c(0.5, -89.5), c(1, 1), cells.dim))
## }
## gridFromNothing <- function() {
##     GridTopology(c(0.5, -89.5), c(1, 1), c(360, 180))
## }
## ## mk.lims <- function(xlim, ylim) {

## ## }



## default 100x100
#gt <- mkGT(tr)

## obj and xlim
#gt <- mkGT(tr, xlim = c(100, 200))

## obj and ylim
#gt <- mkGT(tr, ylim = c(-80, 20))

## obj and xlim and ylim
#gt <- mkGT(tr, xlim = c(100, 200), ylim = c(-80, 20))

## xlim and ylim
#gt <- mkGT(xlim = c(100, 200), ylim = c(-80, 20))

## object and cellsize
#gt <- mkGT(tr, cellsize = c(50, 50))

## xlim and ylim and cellsize
#gt <- mkGT(xlim = c(100, 200), ylim = c(-80, 20), cellsize = c(50, 50))



#trg <- tripGrid(tr, dur = 36000, grid = gt)





#trackDistance <-
#function (track, longlat = FALSE)
#{
#    if (!is.matrix(track))
#        stop("track must be two-column matrix")
#    if (ncol(track) != 2)
#        stop("track must be two-column matrix")
#    n1 <- nrow(track) - 1
#    if (n1 < 2)
#        stop("less than two points")
#    res <- numeric(n1)
#    for (i in seq(along = res))
#      res[i] <- spDistsN1(track[i,,drop = FALSE],
#                          track[(i + 1),,drop = FALSE ], longlat = longlat)
#    res
#}


trackDistance <-
function (track, longlat = FALSE, push = 1)
{

	#print(longlat)
	#track <- coordinates(spData)
    if (!is.matrix(track))
        stop("track must be two-column matrix")
    if (ncol(track) != 2)
        stop("track must be two-column matrix")
    n1 <- nrow(track) - 1
    if (n1 < 1)
        stop("less than two points")
    res <- numeric(n1 - push + 1)
    for (i in seq(along = res))
      res[i] <- spDistsN1(track[i,,drop = FALSE],
                          track[(i + push),,drop = FALSE ], longlat = longlat)
    res
}


speedfilter <- function (x, max.speed = NULL, test = FALSE)
{
    if (!is(x, "trip"))
        stop("only trip objects supported")
    projected <- is.projected(x)
    if (is.na(projected)) {
        projected <- FALSE
        warning("coordinate system is NA, assuming longlat . . .")
    }
    FUN <- function(x, aadBUG = FALSE) {
        sqrt(sum((x)^2, na.rm = FALSE)/(if (aadBUG) 1 else length(x)))
    }
    if (is.null(max.speed)) {
        print("no max.speed given, nothing to do here")
        return(x)
    }
    longlat <- !projected
    coords <- coordinates(x)
    time = x[[x@TOR.columns[1]]]
    id = factor(x[[x@TOR.columns[2]]])
    x <- coords[, 1]
    y <- coords[, 2]
    aadBUG = FALSE
    pprm <- 3
    grps <- levels(id)
    if (length(x) != length(y))
        stop("x and y vectors must be of same\nlength")
    if (length(x) != length(time))
        stop("Length of times not equal to number of points")
    okFULL <- NULL
if (test) res <- list(speed = numeric(0), rms = numeric(0))
    for (sub in grps) {
        ind <- id == sub
        xy <- matrix(c(x[ind], y[ind]), ncol = 2)
        tms <- time[ind]
        npts <- nrow(xy)

        if (pprm%%2 == 0 || pprm < 3)
            stop("Points per running mean should be odd and greater than 3, pprm = 3")
        RMS <- rep(max.speed + 1, npts)
        offset <- pprm - 1
        ok <- rep(TRUE, npts)
        if (npts < (pprm + 1)) {
              warning("Not enough points to filter ID: \"", sub, "\"\n continuing . . . \n")
              okFULL <- c(okFULL, ok)

              next;
            }
        index <- 1:npts
iter <- 1
        while (any(RMS > max.speed, na.rm = TRUE)) {
            n <- length(which(ok))
            speed1 <- trackDistance(xy[ok, ], longlat = longlat)/(diff(unclass(tms[ok]))/3600)
            speed2 <- trackDistance(xy[ok, ], longlat = longlat, push = 2)/((unclass(tms[ok][-c(1, 2)]) -
                unclass(tms[ok][-c(n - 1, n)]))/3600)
            thisIndex <- index[ok]
            npts <- length(speed1)
            if (npts < pprm) {
                next
            }
            sub1 <- rep(1:2, npts - offset) + rep(1:(npts - offset),
                each = 2)
            sub2 <- rep(c(0, 2), npts - offset) + rep(1:(npts -
                offset), each = 2)
            rmsRows <- cbind(matrix(speed1[sub1], ncol = offset,
                byrow = TRUE), matrix(speed2[sub2], ncol = offset,
                byrow = TRUE))
            RMS <- c(rep(0, offset), apply(rmsRows, 1, FUN, aadBUG = aadBUG))

	if (test & iter == 1) {
                res$speed <- c(res$speed, 0, speed1)
		    res$rms <- c(res$rms, 0, RMS)
			break
	}
iter <- iter + 1
            bad <- RMS > max.speed
            bad[is.na(bad)] <- FALSE
            segs <- cumsum(c(0, abs(diff(bad))))
            segs[RMS <= max.speed] <- NA
            peaks <- tapply(RMS, segs, which.max)
            for (i in levels(as.factor(segs))) {
                RMS[segs == i & !is.na(segs)][peaks[[i]]] <- NA
            }
            RMS[1] <- 0
            RMS[length(RMS)] <- 0
            ok[thisIndex][is.na(RMS)] <- FALSE
        }
        okFULL <- c(okFULL, ok)
    }
	if (test) return(res)
    filt <- okFULL
    filt
}




adjust.duplicateTimes <- function (time, id)
{
    dups <- unlist(tapply(time, id, duplicated), use.names = FALSE)
    if (any(dups)) {
        time[dups] <- time[dups] + 1
        time <- Recall(time, id)
    }
    time
}


argos.sigma <- function(x, sigma = c(100, 80, 50, 20, 10, 4,  2),
	                   adjust = 111.12) {
	sigma <- sigma / adjust
	names(sigma) <- levels(x)
	sigma[x]
}

readArgos <-
  ## add "correct.all" argument - just return data frame if it fails, with
  ## suggestions of how to sort/fix it
function (x, correct.all = TRUE, dtFormat = "%Y-%m-%d %H:%M:%S",
    tz = "GMT", duplicateTimes.eps = 1e-2, p4 = "+proj=longlat +ellps=WGS84", verbose = FALSE)
{
    dout <- NULL
    for (con in x) {
        old.opt <- options(warn = -1)
        dlines <- strsplit(readLines(con), "\\s+", perl = TRUE)
        options(old.opt)
        loclines <- sapply(dlines, length) == 12
        if (any(loclines)) {
            dfm <- matrix(unlist(dlines[sapply(dlines, length) ==
                12]), ncol = 12, byrow = TRUE)

            if (dfm[1,7] == "LC") {
            	cat("file ", con, " appears to be a diag file, skipping. Use readDiag to obtain a dataframe. \n\n")
            	next
            }
            df <- vector("list", 12)
            names(df) <- c("prognum", "ptt", "nlines", "nsensor",
                "satname", "class", "date", "time", "latitude",
                "longitude", "altitude", "transfreq")
            for (i in c(1:4, 9:12)) df[[i]] <- as.numeric(dfm[, i])
            for (i in 5:6) df[[i]] <- factor(dfm[, i])
            for (i in 7:8) df[[i]] <- dfm[, i]
            df <- as.data.frame(df)
            df$gmt <- as.POSIXct(strptime(paste(df$date, df$time),
                dtFormat), tz)
            dout <- rbind(dout, df)
        }
        else {
            cat("Problem with file: ", con, " skipping\n")
        }
    }
    if (is.null(dout))
        stop("No data to return: check the files")
    if (correct.all) {
      ## should add a reporting mechanism for these as well
      ##  and return a data.frame if any of the tests fail
      ## sort them
        dout <- dout[order(dout$ptt, dout$gmt), ]
      ## remove duplicate rows
        dout <- dout[!duplicated(dout), ]
        ## adjust duplicate times (now that they are sorted properly)

        dt.by.id <- unlist(tapply(dout$gmt, dout$ptt, function(x) c(-1, diff(x))))

        dup.by.eps <- which(abs(dt.by.id) < duplicateTimes.eps)

        if (length(dup.by.eps) >= 1) {
          if (verbose) {
              cat("Adjusting duplicate times\n.....\n")
          for (i in  dup.by.eps) {ind <- i + (-2:1);print(cbind(dout[ind,c("ptt", "gmt", "class")], row.number = ind))}
      }
          dout$gmt <- adjust.duplicateTimes(dout$gmt, dout$ptt)
          if (verbose) {
              cat("\n  Adjusted records now: \n\n")
              for (i in  dup.by.eps) {ind <- i + (-2:1);print(cbind(dout[ind,c("ptt", "gmt", "class")], row.number = ind))}
          }

        }
        if(any(dout$longitude > 180)) {
        	cat("\nLongitudes contain values greater than 180, assuming proj.4 +over\n\n")
        	p4 <- "+proj=longlat +ellps=WGS84 +over"
        	}
        dout$class <- ordered(dout$class, levels = c("Z", "B", "A", "0", "1", "2", "3"))

        coordinates(dout) <- c("longitude", "latitude")
        proj4string(dout) <- CRS(p4)
        #tor <- TimeOrderedRecords(c("gmt", "ptt"))
        test <- try(dout <- trip(dout, c("gmt", "ptt")))
        if (!is(test, "trip")) {cat("\n\n\n Data not validated: returning object of class ", class(dout), "\n");return(dout)}

        ## for now, only return spdftor if correct.all is TRUE
        cat("\n\n\n Data fully validated: returning object of class ", class(dout), "\n")
        return(dout)
        }
    cat("\n\n\n Data not validated: returning object of class ", class(dout), "\n")
    dout
}



sepIdGaps <- function(id, gapdata, minGap = 3600 * 24 * 7) {

	toSep <- tapply(gapdata, id, function(x) which(diff(unclass(x) ) > minGap))

	tripID <- split(as.character(id), id)


	for (i in 1:length(tripID)) {
		this <- toSep[[i]]
		thisID <- tripID[[i]][1]
			if (length(this) > 0) {

				for (n in 1:length(this)) {
					tripID[[i]][(this[n]+1):length(tripID[[i]])] <- paste(thisID, n + 1, sep = "_")
				}

		}

	}
	unsplit(tripID, id)
}
