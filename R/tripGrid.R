
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
            ## trim psp objects (0-lines give NaNs)

            zero.lengths <- TRUE
            zeros <- which(!lngths > 0)
            cc <- coordinates(this)[zeros, , drop = FALSE]
            op <- options(warn = -1)
            x.ppp <- ppp(cc[, 1], cc[, 2], window = ow)
            options(op)
            if (method == "pixellate") {
                v <- pixellate(x.ppp, W = ow, weights = dt[zeros])$v
            }
            if (method == "density") {
                v <- density(x.ppp, ...)$v


            }

            res$z <- res$z + t(v)
            sz <- sz + sum(dt[zeros])
        }
        x.psp <- x.psp[lngths > 0]
        weights <- dt/ifelse(lngths > 0, lngths, .Machine$double.eps)
        weights <- weights[lngths > 0]
        if (method == "pixellate") {
            v <- pixellate(x.psp, W = ow, weights = weights)$v
        }
        if (method == "density") {
            #v <- density(x.psp, ...)$v
                for (li in 1:x.psp$n) {
                 dens <- density(x.psp[li], ...)$v
                 if (li == 1) {
                     v <- dens
                 } else {
                     v <- v + dens * dt[li]
                 }

             }
            }
        res$z <- res$z + t(v)


    }
    if (zero.lengths) {
        warning("zero length lines present, time durations binned into cells assuming point-presence of degenerate line segment")
        cat("\n")
        if (method == "pixellate") {
        cat(paste("Total time of trips:", sm, "\n"))
        cat(paste("Total time without zero length lines:", sm -
            sz, "\n"))
    }
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
