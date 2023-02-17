###############################################################################
# R file containing the functions necessary for working with Statec data
# Most of them are designed to handle time series well
#
# Author: Andrei V. Kostyrka
# Entity: Statec
# Dates: 2022-01 -- 2023-02
###############################################################################

.ef <- function(e) return(NULL) # For try-catching and debugging
.efv <- function(e) {print(e); return(NULL)} # For try-catching and debugging

#' Compute lagged series padded with NA
#'
#' @param x A numeric vector, matrix or data frame
#' @param lag An integer: lag order
#'
#' @return An object of the same dimensions as the input.
#' @export
#'
#' @examples
#' cbind(1:10, myLag(1:10))
myLag <- function(x, lag = 1L) {
  if (lag == 0L) return(x)
  if (length(dim(x)) == 2) {
    ret <- apply(x, 2, myLag, lag = lag)
    if (is.data.frame(x)) ret <- as.data.frame(ret)
    if ("mts" %in% class(x)) ret <- stats::ts(ret, start = stats::start(x), freq = stats::frequency(x))
    return(ret)
  } else if (length(dim(x)) == 0) {
    xo <- c(rep(NA, lag), x[1:(length(x)-lag)])
    if (!is.null(names(x))) names(xo)[1] <- names(x)[1]
    if ("ts" %in% class(x)) xo <- stats::ts(xo, start = stats::start(x), freq = stats::frequency(x))
    return(xo)
  } else stop("The input should be a matrix, data frame, or a vector.")
}

#' Compute differenced series padded with NA
#'
#' @param x A numeric vector, matrix or data frame
#' @param lag An integer: lag order for differencing (e.g. 12 for seasonal differences of monthly data)
#' @param order An integer: differencing order
#'
#' @return An object of the same dimensions as the input.
#' @export
#'
#' @examples
#' myDiff(mtcars)
#' myDiff(mtcars, lag = 2)
#' myDiff(mtcars, lag = 3, order = 2)
myDiff <- function(x, lag = 1L, order = 1L) {
  if (order == 0L) return(x)
  if (lag == 0L) stop("Differences with lag 0 are not possible. The lag should be at least 1.")
  if (length(dim(x)) == 2) {
    ret <- apply(x, 2, myDiff, lag = lag, order = order)
    if (is.data.frame(x)) ret <- as.data.frame(ret)
    if ("mts" %in% class(x)) ret <- stats::ts(ret, start = stats::start(x), freq = stats::frequency(x))
    return(ret)
  } else if (length(dim(x)) == 0) {
    xo <- c(rep(NA, lag*order), diff(x = x, lag = lag, differences = order))
    if (!is.null(names(x))) names(xo)[1] <- names(x)[1]
    if ("ts" %in% class(x)) xo <- stats::ts(xo, start = stats::start(x), freq = stats::frequency(x))
    return(xo)
  } else stop("The input should be a matrix, data frame, or a vector.")
}

#' Return the decimal date with a safeguard against floating-point .9999
#'
#' @param x A numeric vector, matrix or data frame
#' @param adjust A small numeric that will be added to the decimal time (should be less than 1 / 366 / 2)
#' @param numeric If TRUE, convert to a real numeric vector instead of a ts.
#'
#' @return An object of the same dimensions as the input.
#' @export
#'
#' @examples
#' x <- ts(2:252, start = c(2002, 2), freq = 12)
#' problem.date <- seq.Date(as.Date("2002-02-01"), to = as.Date("2022-12-01"), by = "month")
#' true.year <- rep(2002:2022, each = 12)[-1]
#' wrong.year <- floor(as.numeric(time(x)))
#' corrected.year <- myYear(x)
#' tail(data.frame(Date = as.character(problem.date), True = true.year,
#'                 Wrong = wrong.year, Adjusted = corrected.year,
#'                 HasDiscrepancy = abs(corrected.year - wrong.year)), 15)
myTime <- function(x, adjust = getOption("ts.eps"), numeric = FALSE) {
  ret <- stats::time(x) + adjust
  if (numeric) ret <- as.numeric(ret)
  return(ret)
}
#' @describeIn myTime Safely convert a ts into integer year
myYear <- function(x, adjust = getOption("ts.eps"), numeric = TRUE) {
  ret <- floor(myTime(x, adjust = adjust, numeric = FALSE))
  if (numeric) ret <- as.integer(ret)
  return(ret)
}

#' Run a command quietly
#'
#' @param x Any R expression.
#'
#' https://stackoverflow.com/a/54136863
#'
#' @return The same result as evaluated by x, only with no output (sinked into tempfile()).
#' @export
#'
#' @examples
#' quiet(print("Test"))
#' quiet(cat("Test\n"))
#' quiet(message("Test")) # Still prints a message
quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}


#' Compute mean and SD without NAs
#'
#' @param x A numeric vector.
#'
#' @return A vector of length 2.
#' @export
#'
#' @examples
#' meanSD(1:10)
meanSD <- function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE))

#' Find non-numeric gaps in the middle
#'
#' @param x Any vector, possibly with NAs, NaNs, Infs etc. in the middle
#'
#' @return Indices of middle gaps. If there are no middle gaps, returns NULL.
#' @export
#'
#' @examples
#' findMiddleGaps(c(rep(NA, 4), 1:10, rep(NA, 3))) # No middle gaps
#' findMiddleGaps(c(1:10, rep(NA, 3), 5:9)) # One middle gap
#' findMiddleGaps(c(NA, 2:4, NA, NA, 7:8, NA, 10:12, NA, NA)) # Two middle gaps
findMiddleGaps <- function(x) {
  good.seq <- is.finite(x)
  good.rle <- rle(good.seq)
  good.rle$values[1] <- TRUE # Ignoring initial bad values
  good.rle$values[length(good.rle$values)] <- TRUE # Ignoring end bad values
  y <- inverse.rle(good.rle)
  return(which(!y))
}

#' Creates periodic time series from an imported data file
#'
#' @param data The data frame, or time series, or tibble to convert
#' @param time.regex A string representing a regular expression to match the time-identifying variable
#'
#' This function tries to guess the periodicity (monthly, quarterly) from the
#' date column, as a human would do, and create regular ts series.
#'
#' @return A ts or mts object.
#' @export
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(Time = seq.Date(as.Date("2001-01-01"), as.Date("2020-01-01"), "month"),
#'            X = rnorm(229), Y = rchisq(229, df = 2), Z = runif(229))
#' dts <- makeTS(d)
#' plot(dts)
makeTS <- function(data, time.regex = "^Date|^Time|^Code") {
  if (is.null(dim(data))) { # A vector?
    if ("ts" %in% class(data)) {
      message("This is already a univariate TS, returning the input, no transformation needed.")
      return(data)
    }
  } else {
    if ("mts" %in% class(data)) {
      message("This is already a multivariate TS, returning the input, no transformation needed.")
      return(data)
    }
  }

  if (any(class(data) %in% c("tbl", "tbl_df"))) data <- as.data.frame(data)
  p <- grepl(time.regex, colnames(data), ignore.case = TRUE)
  if (!any(p)) stop("No 'Date' or 'Time' column!") else p <- which(p) # Getting the index
  if ("POSIXct" %in% class(data[, p])) data[, p] <- as.Date(data[, p])
  if ("Date" %in% class(data[, p])) {
    ind <- data[, p]
  } else if (is.character(data[, p]) | is.factor(data[, p])) { # The date was read as a character
    ind <- tryCatch(as.Date(data[, p]), error = function(e) {warning("Could not convert the candidate for time index into Date. Is the 'Date' or 'Time' column a properly encoded Date or Time object in the source? If not, is it at least in ISO YYYY-MM-DD or %Y/%m/%d format?"); stop(e)},
                    warning = function(e) {warning("Conversion of the candidate for the time index into Date threw a warning. Aborting for safety purposes. Is the 'Date' or 'Time' column a properly encoded Date or Time object in the source? If not, is it at least in ISO YYYY-MM-DD or %Y/%m/%d format?"); stop(e)})
  } else if (is.numeric(data[, p])) {
    # Is it in integer years?
    if (all(data[, p] == floor(data[, p]))) {
      ind <- as.Date(paste0(data[, p], "-01-01"))
    } else {
      ind <- lubridate::date_decimal(data[, p])
    }
  } else stop("Check the date variable; it is not a date, character, factor, or numeric.")

  data <- data[order(ind), -p]
  n <- if (is.null(dim(data))) length(data) else nrow(data)
  diff.days <- as.numeric(difftime(ind[2:n], ind[1:(n-1)], units = "days"))
  yr <- lubridate::year(ind[1])
  if (all(diff.days %in% 28:31)) { # Days in months
    freq <- 12
    first <- lubridate::month(ind[1])
  } else if (all(diff.days %in% 90:92)) { # Days in quarters
    freq <- 4
    first <- lubridate::quarter(ind[1])
  } else if (all(diff.days %in% 365:366)) {
    freq <- 1
    first <- NULL
  } else {
    tb <- table(diff.days)
    message(paste0("Differences between dates: {", paste0(names(tb), collapse = ", "), "} of multiplicity [", paste0(tb, collapse = ", "), "]."))
    stop("The dates in the series do not form regular months, quarters, or years!")
  }
  ind1 <- date2ind(ind[1], freq = freq)
  ind2 <- date2ind(ind[length(ind)], freq = freq)
  ret <- stats::ts(data, start = c(yr, first), frequency = freq)
  spn <- paste0(ind1[1], if (freq != 1) "-" else "", if (freq != 1) ind1[2] else "",
                " ... ", ind2[1], if (freq != 1) "-" else "", if (freq != 1) ind2[2] else "")
  message(paste0("Creating ", if (is.null(dim(ret))) "1" else ncol(ret), " regular ", if (freq == 4) "QUARTERLY" else if (freq == 12) "MONTHLY" else "ANNUAL", " time series (", spn, ")."))
  attr(ret, "freq") <- freq
  return(ret)
}

#' Get indices of a vector of numbers in another vector
#'
#' @param what A vector to look up.
#' @param where A vector in which to look up.
#'
#' @return A vector of indices (length 0 or above).
#' @export
#'
#' @examples findVecPos(1:3, c(7, 4, 1, 2, 3, 5))
findVecPos <- function(what, where) {
  m <- length(what)
  l <- length(where)
  if (m > l) return(NULL)
  if (m == l) return(if (all(what == where)) 1 else NULL)
  idx <- which(where == what[1])
  idx <- idx[idx <= (l-m+1)]
  if (length(idx) == 0) return(NULL)
  ret <- idx[sapply(idx, function(i) all(where[i:(i+(length(what)-1))] == what))]
  ret <- ret[!is.na(ret)] # If the first value of the match is too close to the end, discard it
}

#' Get the first and last non-missing observation and the number of NAs
#'
#' @param x A vector of any type
#'
#' @return A vector: the index of the first non-NA observation,
#' the index of the last non-NA observation, total valid observations,
#' the number of middle NA sequences, and their total count (middle NAs only)
#' @export
#'
#' @examples getRangeVec(c(NA, NA, 1:3, NA, 1:3, NA, NA, 2, NA))
getRangeVec <- function(x) {
  l <- length(x)
  n <- is.na(x)
  if (all(!n)) { # Case 1: all complete
    first <- 1
    last <- l
    na.runs <- 0
    na.count <- 0
  } else if (all(n)) { # Case 2: all missing
    first <- last <- NA
    na.runs <- 1
    na.count <- l
  } else {
    r <- rle(n) # Run-length encoding
    m <- length(r$lengths)
    begin.miss <- isTRUE(r$values[1]) # Any missingness at the beginning?
    end.miss <- isTRUE(r$values[m]) # Any missingness at the end?
    first <- if (!begin.miss) 1 else which(!n)[1] # First non-missing
    last <- if (!end.miss) l else l - r$lengths[m] # Length of the series minus the length of the last run of NAs
    mid.na <- findVecPos(c(0, 1, 0), as.numeric(r$values)) # Any NAs in the middle?
    if (length(mid.na) > 0) mid.na <- mid.na + 1 # Indices of the middle value
    na.runs <- length(mid.na)
    na.count <- sum(r$lengths[mid.na])
  }
  return(c(first = first, last = last, valid = sum(!n), na.runs = na.runs, mid.na.count = na.count))
}

#' Extract the valid range of observations
#'
#' @param data A ts or mts object
#' @param date.format Output format compatibility. EViews: 2021Q4, X13: 2021.4,
#' Date: 2021-04-01, time: 2021.25, index: integer position
#'
#' @return A vector or a data frame with the results of getRangeVec() combined
#' @export
#'
#' @examples
#' getRange(ts(cbind(x1 = c(NA, NA, 1:3, NA, 1:3, NA, NA, 2, NA), x2 = 1:13),
#'          start = c(2020, 1), freq = 12), date.format = "Date")
getRange <- function(data, date.format = c("eviews", "X13", "Date", "time", "index")) {
  date.format <- date.format[1]
  if (!any(c("ts", "mts") %in% class(data))) stop("This function can be applied to TS or MTS only.")
  f <- stats::frequency(data)
  if (!(f %in% c(4, 12))) stop("Only monthly or quarterly data are expected.")

  if (is.null(dim(data))) { # A vector?
    r <- getRangeVec(data)
    dates <- matrix(r, nrow = 1, dimnames = list("Y", names(r)))
  } else {
    dates.list <- lapply(1:ncol(data), function(i) getRangeVec(data[, i]))
    dates <- do.call(rbind, dates.list)
    dates <- as.data.frame(dates)
    colnames(dates) <- names(dates.list[[1]])
  }
  d <- stats::time(data) + getOption("ts.eps")
  yrs <- floor(as.numeric(d))
  period <- as.integer(stats::cycle(d))

  if (is.data.frame(dates)) {
    labels <- switch(date.format, eviews = paste0(yrs, if (f == 12) "M" else "Q", if (f == 12) sprintf("%02d", period) else period),
                     X13 = paste0(yrs, ".", period),
                     Date = ts2date(data),
                     time = d,
                     index = NULL)
  }
  if (!is.null(labels)) {
    dates$first <- labels[dates[, "first"]]
    dates$last  <- labels[dates[, "last"]]
  }
  return(dates)
}


#' Add labels with a white halo for readability
#'
#' @param x Passed to text().
#' @param y Passed to text().
#' @param labels Passed to text().
#' @param nhalo The number of semi-transparent copies of the object to add as a background.
#' @param halo.col Colour of halo elements (default: white semi-transparent).
#' @param hscale Horizontal radius of the halo ellipse (multiple of current plot width)
#' @param vscale Vertical radius of the halo ellipse (multiple of current plot height)
#' @param abs.h Logical: use absolute units of x instead of the current width
#' @param abs.v Logical: use absolute units of y instead of the current height
#' @param ... Passed to text(); both the main text and the halos are invoked with it.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' set.seed(1)
#' plot(runif(2000), pch = 16)
#' points(500, 0.6, pch = 16, col = 2)
#' labelsWithHalo(500, 0.6, "Point of interest", nhalo = 32, pos = 3)
labelsWithHalo <- function(x, y, labels,
                           nhalo = 16, halo.col = "#FFFFFF66",
                           hscale = 0.02, vscale = 0.02,
                           abs.h = FALSE, abs.v = FALSE, ...) {
  lims <- graphics::par("usr")
  xlim <- lims[1:2]; ylim <- lims[3:4]
  hs <- if (abs.h) hscale else diff(xlim) * hscale
  vs <- if (abs.v) vscale else diff(ylim) * vscale
  offsets <- cbind(cos(seq(0, 2*pi, length.out = nhalo + 1)) * hs,
                   sin(seq(0, 2*pi, length.out = nhalo + 1)) * vs)[-(nhalo + 1), ]
  dotargs <- list(...)
  dotargs$labels <- labels
  dotargs$col <- halo.col
  for (i in 1:nhalo) {
    dotargs$x <- x + offsets[i, 1]
    dotargs$y <- y + offsets[i, 2]
    do.call(graphics::text, dotargs)
  }
  graphics::text(x = x, y = y, labels = labels, ...)
  return(invisible(NULL))
}

#' Quickly create a formula y ~ x1 + x2 + ...
#'
#' @param y The left-hand-side variable (character scalar)
#' @param x The right-hand-side variables (character vector)
#'
#' @return A formula with the requested variables.
#' @export
#'
#' @examples wrapFormula("mpg", c("cyl", "hp", "wt"))
wrapFormula <- function(y, x) {
  stats::formula(paste0(y, " ~ ", paste0(x, collapse = " + ")))
}

#' Generated visually distinct colours
#'
#' @param n Integer: how many colours to generate?
#' @param nsamp Integer: how many colours to sample randomly in a 3D cube (YCbCr)?
#' @param range1 Luma range of the YCbCr pallette: two numbers from [0, 1]
#' @param range2 Blue difference range: two numbers from [0, 1]
#' @param range3 Red difference range: two numbers from [0, 1]
#' @param seed Random seed for uniform generation
#' @param norm Numeric scalar: 1 = Manhattan norm, 2 = Euclidean
#' @param bad.margin Exclude points within this distance from existing ones
#' @param start.corner If TRUE, picks a corner of the 3D cube as the 1st colour, otherwise the middle
#' @param reduce A function of two arguments. to reduce two colour distances.
#' @param plot If TRUE, shows the generated colours in a simple plot.
#'
#' This function generates many colours uniformly in the YCbCr space and
#' minimises a function ('reduce') of two distances (in YCbCr and RGB spaces).
#' These distances allow arbirary powers (1 = Manhattan, 2 = Euclidean)
#' \code{reduce} can be set to something like \code{function(x) 0.9*min(x) + 0.1*max(x)}
#' to use the distance closer to the smaller of the two.
#'
#' This algorithm will fail if bad.margin is too large; if many colours are needed, reduce it to prevent the space shrinkage
#'
#' @return A vector of colours that are as distinct as possible.
#' @export
#'
#' @examples
#' cl <- createDistinctCols(20, range1 = c(0, 0.6), plot = TRUE, seed = 2)
#' cl2 <- createDistinctCols(20, range1 = c(0, 0.6), norm = 1, plot = TRUE, seed = 2)
createDistinctCols <- function(n = 4, nsamp = 200,
                               range1 = c(0, 0.5), range2 = c(0.05, 0.95), range3 = c(0.05, 0.95),
                               norm = 2, reduce = min, bad.margin = 0.2,
                               start.corner = TRUE, seed = 1, plot = FALSE) {
  set.seed(seed)
  ranges <- cbind(range1, range2, range3)
  x.ycc01 <- matrix(stats::runif(nsamp*3), ncol = 3)
  if (start.corner) x.ycc01[1,] <- apply(x.ycc01, 2, function(x) range(x)[1 + stats::rbinom(1, 1, 0.5)])
  x.ycc <- sapply(1:3, function(i) x.ycc01[, i]*abs(diff(ranges[, i])) + min(ranges[, i]))
  x.rgb <- imager::as.cimg(array(x.ycc*255, dim = c(nsamp, 1, 1, 3)))
  x.rgb <- imager::YCbCrtoRGB(x.rgb)[, 1, 1, ]

  Dist <- function(x, y, hue1 = FALSE) {
    d <- x - y
    if (is.null(dim(d))) d <- matrix(d, ncol = 3)
    if (hue1)  d[d[, 1] > 0.5, 1] <- (1 - d[d[, 1] > 0.5, 1])*2
    r <- rowSums(abs(d)^norm)^(1/norm)
    if (is.null(dim(x))) r <- as.numeric(r)
    return(r)
  }

  inds <- utils::combn(seq_len(nsamp), 2)
  dist1 <- dist2 <- matrix(ncol = nsamp, nrow = nsamp)
  dist1[lower.tri(dist1)] <-  Dist(x.ycc[inds[1,], ], x.ycc[inds[2,], ], hue1 = FALSE)
  dist2[lower.tri(dist2)] <-  Dist(x.rgb[inds[1,], ], x.rgb[inds[2,], ], hue1 = FALSE)
  dist1[upper.tri(dist1)] <- t(dist1)[upper.tri(dist1)]
  dist2[upper.tri(dist2)] <- t(dist2)[upper.tri(dist2)]
  diag(dist1) <- diag(dist2) <- -Inf
  dist2[1, ] <- dist2[1, ] <- -Inf # The first element has been chosen
  dists <- array(c(dist1, dist2), dim = c(nsamp, nsamp, 2))
  dists <- apply(dists, c(1, 2), reduce)
  dists.old <- dists

  best.inds <- integer(n)
  best.inds[1] <- 1
  bad.list <- vector("list", n)
  cols <- c("#00000044", grDevices::rainbow(n-1, end = 0.8, v = 0.8, alpha = 0.8))
  col.inds <- rep(1, nsamp)

  for (i in 1:(n-1)) {
    # cx <- rep(0.8, nsamp); cx[best.inds[i]] <- 3
    # p  <- suppressWarnings(persp(c(0, 1), c(0, 1), matrix(rep(-100, 4), 2), xlim = range(x.ycc[, 1]), ylim = range(x.ycc[, 2]), zlim = range(x.ycc[, 3]), theta = 135 + i*2, scale = FALSE, ticktype = "detailed", phi = 10+i, xlab = "X", ylab = "Y", zlab = "Z"))
    # points(trans3d(x.ycc[, 1], x.ycc[, 2], x.ycc[, 3], p), col = cols[col.inds], pch = 16, cex = cx)
    # Sys.sleep(0.5)
    this.d <- dists[, best.inds[i]]
    new.bad <- is.finite(this.d) & this.d < bad.margin
    new.bad[best.inds[i]] <- TRUE
    bad.list[[i]] <- which(new.bad)
    col.inds[new.bad] <- i+1
    dists[new.bad, best.inds[i]] <- -Inf
    penal <- rowSums(dists[ , best.inds[1:i], drop = FALSE])
    furthest.i <- which.max(penal)
    best.inds[i+1] <- furthest.i
    dists[furthest.i, ] <- -Inf
    # print(rbind(Current = x[best.inds[i], ], apply(x[new.bad, ], 2, range), Furthest = x[best.inds[i+1], ]))
  }
  bad.lens <- sapply(bad.list[-n], length)
  bad.loss <- c(0, cumsum(bad.lens) / nsamp)
  res <- x.rgb[best.inds, ]

  fin <- grDevices::rgb(res[, 1], res[, 2], res[, 3])
  if (plot) plot(1:n, bad.loss, pch = 16, cex = 2, col = fin, bty = "n", xlab = "Colours", ylab = "Remaining % of points",
                 main = "Colour space shrinkage in colour picking", sub = "Last point close to 1 = the points are spaced out evenly")
  if (plot) graphics::abline(h = c(0, 1), lty = 2)
  return(fin)
}

#' Compute a weighted mean of ordered matrix rows with friendly weight names
#'
#' @param x A numeric matrix for which the ordered rows should be averaged with weights
#' @param w A numeric vector of weights or name of the weighting function: "mean", "median", "IQM", "Bartlett", "Parzen", "Tukey-Hanning", "trimean", "midhinge"
#' @param order A numeric vector of column orders or numeric matrix of orders within each row.
#'
#' By default, uses triangular (Bartlett) weights.
#' If order = NULL, uses apply(x, 1, order) values to order the values in rows of x.
#' If order is a vector, uses these to sort the columns of x.
#' If order is a matrix, uses the rows of order to sort the rows of x.
#'
#' @return A numeric vector of weighted mean of ordered x.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rt(20*9, df = 1), ncol = 9) # series with many outliers
#' matplot(x, type = "l", lty = 1)
#' xm <- cbind(orderedWMean(x, w = c("mean", "iqm", "median", "trimean", "Bartlett")),
#'             omit.ends = orderedWMean(x, w = c(0, rep(1, ncol(x)-2), 0)))
#' mycol <- rainbow(6, v = 0.7, end = 0.8, alpha = 0.8)
#' matplot(xm, type = "l", col = mycol, lty = 1, lwd = 2)
#' legend("topleft", colnames(xm), ncol = 2, lwd = 2, col = mycol)
orderedWMean <- function(x, w = NULL, order = NULL) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  k <- ncol(x)
  k0 <- floor(k/2)
  if (is.null(w)) w <- "Bartlett"
  if (is.character(w)) {
    if (length(w) > 1) { # Multiple methods at once
      out <- sapply(seq_along(w), function(i) orderedWMean(x, w = w[i], order = order))
      colnames(out) <- w
      if (stats::is.ts(x)) out <- stats::ts(out, start = stats::start(x), frequency = stats::frequency(x))
      return(out)
    }

    if (w == "mean") ws <- rep(1, k)
    if (w == "median") ws <- if (k %% 2 == 1) c(rep(0, k0), 1, rep(0, k0)) else c(rep(0, k0-1), 1, 1, rep(0, k0-1))
    if (w %in% c("iqm", "IQM")) ws <-
        if (k %% 4 == 0) rep(c(0, 1, 0), times = k/c(4, 2, 4)) else
          if (k %% 4 == 2) rep(c(0, 1, 2, 1, 0), times = c(floor(k/4), 1, k/2-1, 1, floor(k/4))) else
            rep(c(0, 1, 2, 1, 0), times = c(floor(k/4), 1, (k-1)/2 - (k %% 4 == 1), 1, floor(k/4)))
    if (w %in% c("Bartlett", "Parzen", "Tukey-Hanning")) {
      kw <- sandwich::kweights(if (k %% 2 == 1) (0:k0)/(k0+1) else (0:(k0-1))/k0, kernel = w)
      ws <- c(rev(kw), kw)
      if (k %% 2 == 1) ws <- ws[-(k0+1)] # Dropping the duplicated middle item
    }
    if (w %in% c("trimean", "midhinge")) {
      qs <- t(apply(x, 1, stats::quantile, probs = 1:3/4))
      if (w == "trimean") wm <- (qs[,1] + 2*qs[,2] + qs[,3]) / 4 else wm <- (qs[,1] + qs[,3]) / 2
      if (stats::is.ts(x)) wm <- stats::ts(wm, start = stats::start(x), freq = stats::frequency(x))
      return(wm)
    }
  } else if (is.numeric(w)) {
    if (length(w) == 1) ws <- rep(1, k)
    if (length(w) != 1 & length(w) != k) stop("Wrong weight vector length (should be ncol(x)).")
    ws <- w
  }
  if (all(ws < 0)) stop("At least one weight should be positive.")
  ws <- ws / sum(ws)
  ws <- matrix(ws, ncol = ncol(x), nrow = nrow(x), byrow = TRUE)

  if (is.null(order)) x.order <- t(apply(x, 1, order)) else
    if (is.vector(order)) {
      if (length(order) == ncol(x)) x.order <- matrix(order, ncol = ncol(x), nrow(x), byrow = TRUE) else stop("The length of order should be ncol(x).")
    } else if (is.matrix(order)) x.order <- order else stop("Wrong input for 'order' (should be matrix or vector).")

  ox <- t(sapply(1:nrow(x), function(i) x[i, x.order[i, ]]))
  wm <- rowSums(ox * ws) # The row sums of weights are 1, no need to divide
  if (stats::is.ts(x)) wm <- stats::ts(wm, start = stats::start(x), freq = stats::frequency(x))
  return(wm)
}

#' Print a table to Excel or TeX
#'
#' @param x A data frame or a matrix to print
#' @param format Character: if "excel", optimises the plain-text look; if TeX, produces a LaTeX table
#' @param digits Integer: Number of decimal positions to print
#' @param zero.replace Logical: replace true zeros with a dash?
#' @param minus.replace Logical: replace the minus dash with a proper typographic minus sign?
#' @param nozero.integers Logical: if TRUE, do not add zeros to the numbers very close to integers
#' @param nozero.tol A number is considered integer and not padded with zeros if it is closer to the nearest integer than this tolerance
#' @param include.rownames Logical: add the row names as the first column?
#' @param na String to use instead of NAs.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' # This output is ready to be copied from the console and paste into Excel
#' printTab(mtcars, digits = 2, zero.replace = FALSE, minus.replace = FALSE)
#' printTab(myDiff(datasets::airquality[1:10, ]), digits = 1,
#'          zero.replace = FALSE, format = "TeX", na = "X")
printTab <- function(x, format = c("excel", "TeX"), digits = 3,
                     zero.replace = TRUE, minus.replace = TRUE,
                     nozero.integers = TRUE, nozero.tol = 1e-6,
                     include.rownames = FALSE, na = "") {
  format <- format[1]
  xn <- rownames(x)
  sep <- if (format == "excel") "\t" else " & "
  myFormat <- function(x, zero.replace, digits, excel, minus.replace, na, nozero, tol) {
    which.zero <- which(x == 0)
    rx <- round(x)
    which.int <- abs(rx - x) <= tol
    x <- sprintf(paste0("%1.", digits, "f"), x)
    if (isTRUE(all(which.int, na.rm = TRUE))) {
      x <- as.character(rx)
    }
    if (!excel) x <- gsub("-", "$-$", x) else if (minus.replace) x <- gsub("-", "\u2212", x)
    if (zero.replace & length(which.zero) > 0) {
      x[which.zero] <- if (!excel) "--" else "\u2013"
    }
    which.na <- x == "NA" | is.na(x)
    if (any(which.na)) x[which.na] <- na
    return(x)
  }
  xd <- as.data.frame(x)
  xdf <- lapply(xd, myFormat, zero.replace = zero.replace, digits = digits,
                  excel = (format == "excel"), minus.replace = minus.replace,
                na = na, nozero = nozero.integers, tol = nozero.tol)
  x <- do.call(cbind, xdf)
  cat(paste0(c(if (include.rownames) "" else NULL, colnames(x)), collapse = sep),
      if (format != "excel") "\\\\\n" else "\n", sep = "")
  # Formatting should occur my column
  for (i in 1:nrow(x)) {
    x.good <- x[i, ]
    cat(if (include.rownames) paste0(xn[i], sep) else NULL,
        paste0(x.good, collapse = sep),
        if (format != "excel") "\\\\\n" else "\n",
        sep = "")
  }
  return(invisible(NULL))

}

#' Print stars based on the coefficients and standard errors
#'
#' @param coef A numeric vector representing coefficients
#' @param se A numeric vector representing standard errors
#' @param lev A numeric vector of p-values that determine the number of stars
#' @param add.number Logical print number with the p-values? If FALSE, returns just the stars
#' @param digits How many digits of the number to print
#'
#' @return A character vector of stars or numbers with stars
#' @export
#'
#' @examples
#' printStars(coef = seq(1, 3, 0.2), se = rep(1, 11), digits = 1)
printStars <- function(coef, se, lev = c(0.10, 0.05, 0.01),
                       add.number = TRUE, digits = 3) {
  if (length(coef) != length(se)) stop("The coefficients and the standard errors must have the same length.")
  tstat <- coef / se
  pval <- 2*stats::pnorm(-abs(tstat))
  plevs <- as.integer(as.character(cut(pval, c(Inf, sort(unique(lev)), -Inf), labels = (length(unique(lev))):0)))
  stars <- sapply(plevs, function(n) paste0(rep("*", n), collapse = ""))
  if (add.number) stars <- paste0(sprintf(paste0("%1.", digits, "f"), coef), stars)
  if (!is.null(names(coef))) names(stars) <- names(coef)
  return(stars)
}

#' Draw 2D time series with evolution made obvious
#'
#' @param x A two-dimensional TS object
#' @param i An index at which row to stop. In NULL, plots the full range
#' @param thicker The number of points at the end to make thicker
#' @param vgrid Vertical grid points
#' @param hgrid Horizontal grid points
#' @param ... Passed to plot() (axis labels etc.) before the line drawing is invoked
#'
#' This function can be run in a loop over \code{i}, and the resulting frames can be saved as an animation.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' # Visualising a co-integrating relationship
#' set.seed(1)
#' x <- cumsum(rnorm(108))
#' y <- x + rnorm(108)
#' d <- ts(cbind(x, y), start = c(1995, 1), frequency = 4)
#' plotTSAnimationFrame(d)
plotTSAnimationFrame <- function(x, i = NULL, thicker = 48, vgrid = NULL, hgrid = NULL, ...) {
  if (!stats::is.ts(x)) stop("x must be a mts object with 2 columns.")
  if (ncol(x) != 2) stop("x must contain 2 series.")
  x <- stats::na.omit(x) # Incomplete observations cannot be plotted
  if (is.null(i)) i <- nrow(x)
  xv <- as.numeric(x[, 1])
  yv <- as.numeric(x[, 2])
  dotlist <- list(...)
  if (is.null(dotlist$xlim)) dotlist$xlim <- range(xv, na.rm = TRUE)
  if (is.null(dotlist$ylim)) dotlist$ylim <- range(yv, na.rm = TRUE)
  dotlist[c("x", "y")] <- list(NULL, NULL)
  dotlist$bty <- "n"
  do.call(plot, dotlist)
  if (is.null(hgrid)) hgrid <- pretty(yv)
  if (is.null(vgrid)) vgrid <- pretty(xv)
  graphics::abline(v = vgrid, lty = 2, col = "#00000033")
  graphics::abline(h = hgrid, lty = 2, col = "#00000033")
  inds <- max(i-thicker, 1):i # These indices are plotted in a bolder and darker line
  cols <- grDevices::colorRampPalette(c("#00000055", "#000000FF"), alpha = TRUE)(length(inds)) # Darker colours
  wds <- exp(seq(log(1.5), log(4), length.out = length(cols))) # Thicker lines towards the end
  graphics::segments(x0 = xv[inds[-1]], y0 = yv[inds[-1]], x1 = xv[inds[-length(inds)]], y1 = yv[inds[-length(inds)]], col = cols[-1], lwd = wds, pch = 16)
  graphics::lines(xv[1:i], yv[1:i], col = cols[1], lwd = 1.5, pch = 16) # Lines older than this are thin and pale

  # Adding features
  yrs <- stats::time(x) + getOption("ts.eps")
  int.years <- round(yrs[abs(yrs - floor(yrs)) < 0.0001])
  these.int.years.inds <- which((abs(yrs - floor(yrs)) < 0.0001) & (yrs <= floor(yrs[i])))

  cexs <- pchs <- rep(1, length(int.years))
  cexs[int.years %% 4 == 0] <- sqrt(2) # Mark every 4th year with a larger point
  cexs[int.years %% 8 == 0] <- sqrt(3) # Mark every 8th year with an even larger point
  pchs[int.years %% 2 == 0] <- 16 # Mark every 2nd year with a filled circle
  ycols <- grDevices::rainbow(length(int.years), v = 0.7, end = 0.65) # Different colours for temporal evolution

  # Plot integer year points up until the chosen time moment
  if (length(these.int.years.inds) > 0) {
    s <- length(these.int.years.inds)
    graphics::points(xv[these.int.years.inds], yv[these.int.years.inds], cex = cexs[1:s], col = ycols[1:s], pch = pchs[1:s], lwd = 2)
  }

  # Plot marked year labels up until the chosen time moment
  marked.years <- seq(to = floor(max(int.years)/4)*4, by = 4, length.out = 7)
  marked.years.inds <- which(int.years %in% marked.years) # Marked point indices
  ycols <- grDevices::rainbow(length(int.years), v = 0.7, end = 0.65)

  # These indices = until the chosen time moment
  ryrs <- round(yrs*4) / 4 # Division by powers of 2 is lossless
  these.marked.years.inds <- which((ryrs %in% marked.years) & (ryrs <= yrs[i]))
  if (length(these.marked.years.inds > 0)) {
    s <- length(these.marked.years.inds)
    ii <- marked.years.inds[1:s]
    labelsWithHalo(x = xv[these.marked.years.inds], y = yv[these.marked.years.inds], labels = int.years[ii], pos = 3, font = 2, nhalo = 8)
  }
  return(invisible(NULL))
}

#' Compute the power spectrum of a series
#'
#' @param x A univariate ts object or a numeric vector
#' @param apply.diff If TRUE, takes the differences (standard in JDemetra+)
#' @param pad Pad the series with 0s up to this many orders of magnitude above the length of x. 0 or negative = no padding.
#' @param return.half If TRUE, return only the first half of the FFT (due to the symmetrical nature of the FFT, the second half mirrors the first one)
#' @param plot If TRUE, produces a plot
#' @param compact If TRUE, suppres the y axis and axis labels
#' @param lwd Spectrogram line width
#' @param put.legend If TRUE, adds a legend (that depends on the \code{compact} argument)
#'
#' A useful tool to visualise seasonal spikes. The presence of spikes indicates that there are periodic phenomena in the data.
#'
#' @return A list of the frequency grid and corresponding frequency intensity.
#' @export
#'
#' @examples
#' computeSpectrum(datasets::AirPassengers) # Strong seasonal components
#' # Adjust the series by LOESS and remove the seasonal fluctuations
#' ap.sa <- stl(datasets::AirPassengers, 9)
#' computeSpectrum(ap.sa$time.series[, "remainder"]) # Strong trading-day effect!
computeSpectrum <- function(x, apply.diff = TRUE, pad = 2, return.half = TRUE,
                            plot = TRUE, compact = FALSE, lwd = 1.5, put.legend = TRUE) {
  if (!stats::is.ts(x)) warning("x is not a ts object. Analysing seasonal and calendar components is harder when the frequency is unknown.")
  xs <- getRangeVec(x)
  if (xs["mid.na.count"] > 0) stop("x contains NA in the middle. Consider imputing first. Aborting.")
  x <- stats::na.omit(x)
  freq <- stats::frequency(x)
  if (!(freq %in% c(4, 12))) {
    warning("The frequency could not be inferred from the data. Cannot draw the peaks.")
    freq <- NULL
  }
  if (apply.diff) x <- diff(x)
  x <- scale(x)
  n <- length(x)
  if (pad > 0) { # Padding some orders of magnitude above
    pow2 <- floor(log(n, base = 2)) + pad + 1 # Plus one because the FFT is symmetric
    n2 <- 2^pow2
    nzeros <- n2 - n
    x <- c(x, rep(0, nzeros))
  }

  l <- length(x)
  fx  <- stats::fft(x)
  k <- if (l %% 2 == 0) l/2 + 1 else ceiling(l/2)
  yv <- as.numeric(abs(fx^2))
  yv <- yv / n
  xv <- seq(0, by = 1/length(fx), length.out = l)

  # Find a free corner for the legend -- where is the average of the peaks is higher?
  yvq1 <- yv[1:floor(k/4)]
  yvq4 <- yv[ceiling(3*k/4):k]
  yvq1m <- mean(yvq1[yvq1 >= stats::median(yvq1)])
  yvq4m  <- mean(yvq4[yvq4 >= stats::median(yvq4)])
  leg.pos <- if (yvq1m > yvq4m) "topright" else "topleft"

  # Fractional part of 365.25 / 12 / 7
  if (plot) {
    if (compact) {
      plot(xv[1:k], yv[1:k], type = "l", bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = lwd)
    } else {
      plot(xv[1:k], yv[1:k], type = "l", bty = "n", xlab = "Frequency", ylab = "Spectrum", lwd = lwd)
    }
    mycols <- c("#00000088", "#E6194BBB", "#157d22BB")
    if (!is.null(freq)) {
      calfreq <-  if (freq == 12) 0:6/12 else if (freq == 4) 0:2/4 else 0:6/12
      if (compact) {
        graphics::axis(1, at = calfreq, labels = 0:round(freq/2))
      } else {
        graphics::axis(1)
        graphics::axis(3, at = calfreq, labels = 0:round(freq/2))
        graphics::mtext("Times per year", line = 3)
      }

      if (freq == 12) { # See derivation in Cleveland & Devlin (1980)
        # Exact values in Handbook on Seasonal Adjustmeny (2018)
        td1freq <- c(0.348125, 0.431458)
        td2freq <- c(0.304, 0.220)
        legtext <- if (compact) c("Month", "TD", "Week/2") else c("Monthly freq.", "Trading-day freq.", "Semi-week. freq.")
      } else if (freq == 4) {
        td1freq <- c(0.294375, 0.205, 0.338750, 0.04464, 0.455)
        td2freq <- c(0.411)
        legtext <- if (compact) c("Qrtr", "TD", "Week/2") else c("Quarterly freq.", "Trading-day freq.", "Semi-week. freq.")
      } else {
        td1freq <- td2freq <- NULL
        legtext <- if (compact) "Month?" else "Monthly (?)"
      }
      graphics::abline(v = calfreq, lty = 2, lwd = lwd, col = mycols[1])
      graphics::abline(v = td1freq, lty = 4, lwd = lwd, col = mycols[2])
      graphics::abline(v = td2freq, lty = 3, lwd = lwd, col = mycols[3])
      if (put.legend) graphics::legend(leg.pos, legtext, lty = c(2, 4, 3), lwd = lwd, col = mycols, bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    } else {
      graphics::abline(v = 0:6/12, lty = 2, lwd = lwd)
      graphics::abline(v = c(0.3482, 0.04464), lty = 3, lwd = lwd)
      legtext <- if (compact) c("Seas. (?)", "TD (?)") else c("Possible seas. freq.", "Possible trading-day freq.")
      if (put.legend) graphics::legend(leg.pos, legtext, lty = c(2, 3), lwd = lwd)
    }

  }

  if (return.half) {
    xv <- xv[1:k]
    yv <- yv[1:k]
  }
  ret <- list(x = xv, y = yv)
  attr(ret, "legend.position") <- leg.pos

  return(invisible(ret))
}


#' Print a line of N characters
#'
#' @param x A character to print.
#' @param times How many times to print it.
#'
#' @return Nothing.
#' @export
#'
#' @examples printSym("*", 30)
printSym <- function(x = "=", times = 60) {
  xline <- paste0(rep(x, times), collapse = "")
  cat(xline, "\n", sep = "")
  return(invisible(NULL))
}


#' A set of functions for date conversion
#'
#' @param x Source vector or matrix of the type described below.
#' @param freq Frequency (4 or 12).
#'
#' time2ind: ts --> (year, cycle): 2012.4167 -> c(2012, 5)
#' ind2time: (year, cycle) --> ts: c(2012, 5) -> 2012.4167
#' date2ind: Date --> (year, cycle): "2012-03-21" -> c(2012, 3) or c(2012, 1) (quarterly)
#' ind2date: (year, cycle) --> Date: c(2012, 3) --> "2012-03-01" or "2012-07-01" (quarterly)
#' ts2date: ts --> Date
#'
#' These functions support matrix input (i.e. ind2time accepts 2-column matrices).
#'
#' @return A converted vector representing the time.
#' @export
#'
#' @examples
#' x <- ts(1:24, start = c(2001, 1), freq = 12)
#' d <- seq(as.Date("2009-01-01"), as.Date("2010-04-01"), by = "week")
#' i <- cbind(rep(2010:2014, each = 4), rep(1:4, 5))
#' time2ind(x)
#' ts2date(x)
#' date2ind(d, freq = 12)
#' ind2date(i, freq = 4)
#' ind2time(i, freq = 4)
time2ind <- function(x, freq) {
  if (stats::is.ts(x)) {
    xt <- stats::time(x) + getOption("ts.eps")
    yr <- as.numeric(floor(xt))
    cyc <- as.numeric(stats::cycle(x))
  } else if (is.numeric(x)) {
    yr <- floor(x)
    cyc <- floor((x - yr) * freq + (1/366/2/1.001)) + 1 # For numerical stability so that floor() does not clip to previous year
  } else stop("x must be ts or numeric (time).")
  ret <- cbind(year = yr, cycle = cyc)
  if (length(x) == 1) ret <- as.numeric(ret)
  return(ret)
}

#' @describeIn time2ind Pair(s) of indices to decimal time
#' @export
ind2time <- function(x, freq) {
  if (is.null(dim(x))) x[1] + (x[2]-1)/freq else x[, 1] + (x[, 2]-1)/freq
}

#' @describeIn time2ind Date to (year, period) index
#' @export
date2ind <- function(x, freq) {
  if ("Date" %in% class(x)) {
    yr <- as.numeric(substr(x, 1, 4))
    cyc <- as.numeric(substr(x, 6, 7)) # Month
    if (freq == 4) cyc <- floor((cyc-1)/3) + 1 else if (freq == 1) cyc <- NULL
  } else stop("x must be a Date.")
  ret <- if (freq != 1) cbind(year = yr, cycle = cyc) else yr
  if (length(x) == 1) ret <- as.numeric(ret)
  return(ret)
}

#' @describeIn time2ind Time series to (year, period) index
#' @export
ind2date <- function(x, freq) {
  if (is.null(dim(x))) x <- matrix(x, 1, 2)
  if (is.numeric(x)) {
    yr <- x[, 1]
    mon <- if (freq == 12) x[, 2] else (x[, 2]-1)*3 + 1
  } else stop("x must be numeric of length 2 or matrix with 2 columns.")
  ret <- paste0(yr, "-", sprintf("%02d", mon), "-", "01")
  ret <- as.Date(ret)

  return(ret)
}

#' @describeIn time2ind Time series to Date
#' @export
ts2date <- function(x) {
  freq <- stats::frequency(x)
  x.year <- floor(as.numeric(stats::time(x)) + getOption("ts.eps"))
  x.mon <- as.numeric(stats::cycle(x))
  if (freq == 4) x.mon <- (x.mon - 1)*3 + 1
  return(as.Date(paste0(x.year, "-", x.mon, "-", "01")))
}

#' Generate a linear ramp from 0 to 1 or possible reversed
#'
#' @param x Time series that has the time attribute
#' @param start Start time (or if there is no time(x) attribute, index)
#' @param end End time (or if there is no time(x) attribute, index)
#' @param rev If TRUE, starts at 1 and goes down to zero; otherwise, starts at 0
#'
#' @return A series of the same length of x where the change from 0 to 1 occurs gradually
#' @export
#'
#' @examples
#' genRamp(1:100, 40, 50, rev = TRUE)
genRamp <- function(x, start = 2001.5, end = 2002.5, rev = FALSE) {
  if (class(x)[1] == "Date") {
    x <- makeTS(data.frame(Date = x, x = rep(0, length(x))))
  }
  ramp <- x
  start.ind <- which(stats::time(x) >= start)[1]
  end.ind <- which(stats::time(x) >= end)[1]
  ramp[1:(start.ind-1)] <- 0
  ramp[(end.ind+1):length(x)] <- 1
  ramp[start.ind:end.ind] <- seq(0, 1, length.out = end.ind - start.ind + 1)
  if (rev) ramp <- 1 - ramp
  return(ramp)
}

#' Load all Excel sheets as a list of vectors
#'
#' @param xlsxFile Character: path to the Excel file.
#' @param verbose If TRUE, prints the loading progress (useful if the file is on the network drive).
#' @param vector If TRUE, returns a list of vectors; otherwise, a list of data frames.
#' @param lowercase If TRUE, converts the first column (supposedly containing the names) to lowercase.
#' @param intercept2c If TRUE, replaces all instances of `"c"` with `"(Intercept)"` in the first column.
#' @param ... Passed to `openxlsx::read.xlsx` that is applied to each detected sheet.
#'
#' To simply load the XLSX file as a list of sheets, invoke the function with
#' `vector`, `lowercase`, and `intercept2c` set to `FALSE`.
#'
#' @return A list of named vectors or (if `vector` is `FALSE`, a list of data frames from the sheets).
#' @seealso [openxlsx::read.xlsx()] for everything that can be passed to `...`.
#' @export
#'
#' @examples
#' \dontrun{
#' # To compute contributions, load an Excel file where the 1st column
#' # of each sheet contains names, and the second one, coefficient values
#' readAllSheets("myfile.xlsx")
#'
#' # Simply load a file as a list of sheets
#' readAllSheets("myfile.xlsx", vector = FALSE, lowercase = FALSE, intercept2c = FALSE)
#' }
readAllSheets <- function(xlsxFile, verbose = TRUE, vector = TRUE,
                          lowercase = TRUE, intercept2c = TRUE, ...) {
  sheet_names <- openxlsx::getSheetNames(xlsxFile)
  l <- length(sheet_names)
  sheet_list <- vector("list", l)
  names(sheet_list) <- sheet_names
  for (i in 1:l) {
    sh <- openxlsx::read.xlsx(xlsxFile, sheet = sheet_names[i], ...)
    if (lowercase) sh[, 1] <- tolower(sh[, 1])
    if (intercept2c) sh[sh[, 1] == "c", 1] <- "(Intercept)"
    if (vector) sh <- structure(sh[, 2], names = sh[, 1])
    sheet_list[[i]] <- sh
    if (verbose) cat("Read sheet ", i, "/", l, "\n", sep = "")
  }
  return(sheet_list)
}

