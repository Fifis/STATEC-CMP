#' Get the key seasonal adjustment statistics
#'
#' @param x An object of class [seasonal::seas], or a list of outputs from
#' [seasonal::seas] or [diagnoseSeasonality()].
#' @param skip.robustM7 If TRUE, does not test the significance of seasonal and
#' year dummies with a consistent VCOV matrix. Saves 0.1 seconds in most cases
#' @param sep A character used as a delimiter separating multiple outlier entries
#' (is used as the `collapse` argument of `paste0()`)
#'
#' `getSAStat()` can be applied to a [seasonal::seas] object or to the output of
#' [diagnoseSeasonality()] regardless of the number of models / classes.
#' If an adjustment is blended from multiple models, `getSAStat()` also computes the revision
#' stability in the overlapping portions.
#'
#' `getSAStatOne()` is a low-level function for one single [seasonal::seas] object.
#'
#'
#' @return A data frame with named values that are the most important for diagnostics
#' and archiving.
#' @export
#'
#' @examples
#' xs1 <- seasonal::seas(AirPassengers, x11 = "")
#' getSAStat(xs1)
#' getSAStatOne(xs1)
#' xs2 <- diagnoseSeasonality(AirPassengers, name = "AirlinePass", verbose = 1)
#' getSAStat(xs2) # This one has a name
#' # A list of multiple models, possibly mixed, coming from seas or diagnoseSeasonality
#' xs3 <- diagnoseSeasonality(AirPassengers, td = 1, leap.year = FALSE, verbose = 1)
#' set.seed(1)
#' xs4 <- diagnoseSeasonality(ts(rnorm(120), start = c(1970, 1), freq = 12),
#'          td = 1, max.length = 7, name = "WhiteNoise")
#' getSAStat(list(xs1, xs2, xs3, xs4))
getSAStat <- function(x, skip.robustM7 = FALSE, sep = " ") {
  extrDS <- function(x) {
    if (isTRUE(attr(x, "type") == "diagnoseSeasonality")) {
      if (!is.null(x$spans)) { # Extracting all spans
        x <- lapply(x$spans, "[[", "seas")
        attr(x, "spans") <- TRUE
      } else {
        cl <- x$call
        x <- x$seas # Extracting from one diagnoseSeasonality output
        attr(x, "spans") <- FALSE
        x$call <- cl # We are interested in the invokation of the external function, not how it invoked the internal one (which is always the same)
      }
    }
    return(x)
  }
  which.ds <- sapply(x, function(y) isTRUE(attr(y, "type") == "diagnoseSeasonality"))
  if (any(which.ds)) { # Some come from diagnoseSeasonality
    x[which.ds] <- lapply(x[which.ds], extrDS) # Extracting from a list of diagnoseSeasonality objects
  }
  x <- extrDS(x) # Extrcting if it is a single model

  if (class(x)[1] == "seas") {
    return(getSAStatOne(x, skip.robustM7 = skip.robustM7, sep = sep))
  } else {
    # Could be one diagnoseSeasonality output or a list
    st <- lapply(x, function(x) {
      if (isTRUE(attr(x, "spans"))) {
        r <- do.call(rbind, lapply(x, getSAStatOne, skip.robustM7 = skip.robustM7, sep = sep))
        rownames(r) <- paste0(r$name[1], "_span", 1:nrow(r))
        r
      } else getSAStatOne(x, skip.robustM7 = skip.robustM7, sep = sep)
    })
    st <- do.call(rbind, st)
    return(st)
  }
}

#' @rdname getSAStat
getSAStatOne <- function(x, skip.robustM7 = FALSE, sep = " ") {
  if (class(x)[1] != "seas") stop("The main argument must be a 'seas' object.")

  mna <- replicate(11, NA, simplify = FALSE)
  names(mna) <- paste0("M", 1:11)
  sna <- replicate(8, NA, simplify = FALSE)
  names(sna) <- c("mean.rev", "mean.abs.rev", "mean.abs.pct.rev", "sd.rev", "rel.mean.abs.rev", "sd.rev.to.sd.sa", "spans", "overlap.nyears")
  sna$spans <- 1
  sna$overlap.nyears <- NA
  out <- c(list(name = "(noname)", start = NA, end = NA, freq = NA, n.obs = NA, n.ef.obs = NA,
                custom.calendar = NA, log = NA, td = NA, easter = NA, leap.year = NA,
                arima = NA, Noutlier = NA, N.AO = NA, N.LS = NA, N.TC = NA,
                datesAO = "", datesLS = "", datesTC = "", AICc = NA, AICc.per.obs = NA),
           mna, Q_M2 = NA, robust.M7 = NA, f.seas = NA, p.seas = NA, f.stab = NA, p.stab = NA,
           sna, # Extra stability diagnostics for spans
           force.annual = NA,
           has.seas = NA, has.calend = NA, seas.rule = NA, seas.value = NA, seas.thresh = NA,
           log.rule = NA, td.rule = NA, easter.rule = NA, ly.rule = NA)

  f <- stats::frequency(x$x)
  first <- substr(as.character(ind2date(stats::start(x$x), f)), 1, 7)
  last <- substr(as.character(ind2date(stats::end(x$x), f)), 1, 7)
  adjdate <- tryCatch(as.character(as.Date(unname(x$udg["date"]), format = "%b %d, %Y")), error = function(e) unname(x$udg["date"]))
  mod <- x$model$arima$model
  outl <- x$model$regression$variables
  outl <- outl[grepl("^ao|^ls|^tc", outl)]
  if (length(outl) > 0) {
    nout <- length(outl)
    noutAO <- sum(grepl("ao", outl))
    noutLS <- sum(grepl("ls", outl))
    noutTC <- sum(grepl("tc", outl))
    outAO <- paste0(outl[grep("ao", outl)], collapse = sep)
    outLS <- paste0(outl[grep("ls", outl)], collapse = sep)
    outTC <- paste0(outl[grep("tc", outl)], collapse = sep)
  } else {
    outAO <- outLS <- outTC <- ""
    nout <- noutAO <- noutLS <- noutTC <- 0
  }
  uselog <- isTRUE(seasonal::transformfunction(x) == "log")
  m <- as.numeric(x$udg[paste0("f3.m", sprintf("%02d", 1:11))])
  names(m) <- paste0("M", 1:11)
  q <- as.numeric(x$udg["f3.qm2"])
  rtest <- if (!skip.robustM7) robustSeasTests(x) else NULL
  rm7 <- if (!skip.robustM7) rtest$robust.m7 else NA
  f.seas <- if (!skip.robustM7) sprintf("%1.2f", rtest$tests[["identifiable"]]$F[2]) else NA
  f.stab <- if (!skip.robustM7) sprintf("%1.2f", rtest$tests[["stable"]]$F[2]) else NA
  p.seas <- if (!skip.robustM7) sprintf("%1.5f", rtest$pval["identifiable"]) else NA
  p.stab <- if (!skip.robustM7) sprintf("%1.5f", rtest$pval["stable"]) else NA

  aicc <- as.numeric(x$udg["aicc"])
  n <- as.integer(x[["udg"]]["nobs"])
  nef <- as.numeric(x[["udg"]]["nefobs"])
  aiccn <- aicc / nef

  # Finding out if a custom calendar was used, and then, calendar regressors
  custom.cal <- ("User-defined Trading Day" %in% x$est$reg$group) | isTRUE(attr(x, "custom.calendar"))
  which.ea <- grep("Easter", x$est$reg$group)
  has.easter <- if (length(which.ea) > 0) x$est$reg$group[which.ea] else "FALSE" # Keeping the value character
  has.ly <- "Leap Year" %in% x$est$reg$group
  ntd <- sum(grepl("Trading Day", x$est$reg$group))

  sn <- attr(x, "seriesname")
  sncall <- as.list(x$call)$x # Series name inferred from the call
  out$name <- if (!is.null(sn)) sn else if (!is.null(sncall)) deparse(sncall) else "Series"
  out$start <- first
  out$end <- last
  out$freq <- f
  out$n.obs <- n
  out$n.ef.obs <- nef
  out$custom.calendar <- custom.cal
  out$log <- uselog
  out$td <- ntd
  out$easter <- has.easter
  out$leap.year <- has.ly
  out$arima <- mod
  out$Noutlier <- nout
  out$N.AO <- noutAO
  out$N.LS <- noutLS
  out$N.TC <- noutTC
  out$datesAO <- outAO
  out$datesLS <- outLS
  out$datesTC <- outTC
  out$AICc <- aicc
  out$AICc.per.obs <- aiccn
  out[names(mna)] <- m
  out$Q_M2 <- q
  out$robust.M7 <- rm7
  out$f.seas <- f.seas
  out$p.seas <- p.seas
  out$f.stab <- f.stab
  out$p.stab <- p.stab
  out$force.annual <- !is.null(x$series$saa)
  if (!is.null(attr(x, "stability.inds"))) {
    out[names(sna)] <- attr(x, "stability.inds")
  }
  if (!is.null(attr(x, "seasonality"))) {
    out$has.seas <- attr(x, "seasonality")["seasonal"]
    out$has.calend <- attr(x, "seasonality")["calendar"]
    out$seas.rule <- attr(x, "rule")
    out$seas.value <- attr(x, "rule.value")
    out$seas.thresh <- attr(x, "threshold")
    out$log.rule <- attr(x, "transform")
    out$td.rule <- attr(x, "td")
    out$easter.rule <- attr(x, "easter")
    out$ly.rule <- attr(x, "ly")
  } else {
    out$has.calend <- !isTRUE(all.equal(x$series$d16, x$series$d10)) # Is the seasonal adjustment the same as total?
    out$seas.rule <- "X13"
  }

  out <- as.data.frame(out)
  rownames(out) <- out$name
  return(out)
}

#' Return dates in EViews format
#'
#' @param x A ts object.
#'
#' @return A character with EViews-compatible dates.
#' @export
#'
#' @examples
#' x <- ts(1:10, start = c(2001, 1), freq = 4)
#' y <- ts(1:20, start = c(2001, 1), freq = 12)
#' dateAsEV(x)
#' dateAsEV(y)
dateAsEV <- function(x) {
  d <- stats::time(x) + getOption("ts.eps")
  f <- stats::frequency(x)
  p <- stats::cycle(d)
  if (!(f %in% c(4, 12))) stop("The data must be a quarterly or monthly TS.")
  if (f == 12) paste0(myYear(d), "M", sprintf("%02d", p)) else paste0(myYear(d), "Q", p)
}

# Helper functions to combine outlier types
# Choosing the same set of outliers that were found in at least 2 models
.OTMed2 <- function(x) { # Outlier type median
  if (!is.character(x)) x <- unlist(x)
  concordant <- FALSE
  if (all(is.na(x))) {
    ret <- NA
    concordant <- TRUE # All NA = concordance
  } else { # It at least one not NA, figure out whether there are 2, or 1 unique value
    tbx <- table(x)
    if (length(tbx) == 1) {
      if (tbx == 2) concordant <- TRUE # The same 2 outliers in a period = concordance
      # Otherwise, one is an outlier and one is NA
      ret <- names(tbx) # Table with potentially omitted entries
    } else if (length(tbx) == 2) { # Choose between the two
      if (all(sort(names(tbx)) == c("ao", "tc"))) ret <- "ao"
      if (all(sort(names(tbx)) == c("ls", "tc"))) ret <- "tc"
      if (all(sort(names(tbx)) == c("ao", "ls"))) ret <- "ao"
    }
  }
  attr(ret, "concordant") <- concordant
  return(ret)
}

# Choosing the same set of outliers that were found in at least 2 models
.OTMed3 <- function(x) { # Outlier type median
  if (!is.character(x)) x <- unlist(x)
  concordant <- FALSE
  if (all(is.na(x))) {
    ret <- NA
    concordant <- TRUE # All NA = concordance
  } else { # It at least one not NA, figure out whether there all 3, 2, or 1 unique value
    tbx <- table(x)
    if (length(tbx) == 1) {
      if (tbx == 3) concordant <- TRUE # The same 3 outliers in a period = concordance
      ret <- names(tbx) # Table with potentially omitted entries
    } else if (length(tbx) == 3) {
      ret <- "tc" # If one has to chose between an AO, LS, and TC, TC seems a reconciliation of the 3 options
    } else if (length(tbx) == 2) {
      if (any(tbx == 2)) ret <- names(tbx)[which(tbx == 2)] # The option that exists in 2 out of 3 cases
      # Otherwise, choose between the two
      if (all(names(tbx) %in% c("ao", "tc"))) ret <- "ao"
      if (all(names(tbx) %in% c("ls", "tc"))) ret <- "tc"
      if (all(names(tbx) %in% c("ao", "ls"))) ret <- "ao"
    }
  }
  attr(ret, "concordant") <- concordant
  return(ret)
}

# Print the comparable AIC stats for a `seas` model
.printSeasAIC <- function(seasmod, txt = "-- Estimated a seasonal model", ndigits = 3) {
  if (!is.null(seasmod)) {
    n <- as.integer(seasmod[["udg"]]["nobs"])
    nef <- as.integer(seasmod[["udg"]]["nefobs"])
    avgAIC <- as.numeric(seasmod[["udg"]]["aicc"]) / nef
    fullAIC <- avgAIC * n
    avgAICs <- sprintf(paste0("%1.", ndigits, "f"), avgAIC)
    fullAICs <- sprintf("%1.1f", fullAIC)
  } else {
    nef <- avgAICs <- fullAICs <- "FAIL"
  }
  cat(txt, " (obs=", nef, "). AICc/obs = ", avgAICs, ", AICcS = ", fullAICs, ".\n", sep = "")
}

#' An improved seasonal plot with extra useful information
#'
#' @param x An object of class "seas" returned by [seasonal::seas] or the output
#' of [diagnoseSeasonality()]
#' @param name A name to be put on the plot
#' @param sa.custom A `ts` of the adjusted series if they differ from the default
#' X13 ones. Useful if the user obtains their SA series from elsewhere or decides
#' that no SA adjustment (only pre-adjustment) should be done.
#' @param extend.xlim The number of years to add at the beginning and the end to extend
#' the plotting range. Useful if the plot should be produced for a larger sample
#' than the effective sample on which the adjustment was done.
#' @param skip.boxplot If `TRUE`, does not produce the box-and-whisker plot lines
#' @param skip.robustM7 If `TRUE`, does not produce the HAC-robust M7
#' @param span.index The index of the plot in a combined adjustment to show.
#' If NULL, plots the combined series in one plot. If "last", or 0, or -1, plots the last span.
#'
#' @return Nothing (NULL).
#' @export
#'
#' @examples
#' x <- diagnoseSeasonality(AirPassengers, transform = FALSE, td = 0,
#'         easter = TRUE, leap.year = TRUE, plot = NA)
#' x.force <- diagnoseSeasonality(AirPassengers, transform = FALSE, td = 0,
#'              easter = TRUE, leap.year = TRUE, plot = NA, force.annual = "regress")
#' x.alt <- seasonal::seas(AirPassengers, x11 = "")
#' plotSeas(x)
#' plotSeas(x.force)
#' plotSeas(x, sa.custom = x.alt$series$d11, skip.boxplot = TRUE)
#' # Multiple spans
#' xm <- diagnoseSeasonality(AirPassengers, max.length = 5, transform = "yes")
#' plotSeas(xm)
plotSeas <- function(x, name = NULL, sa.custom = NULL,
                     extend.xlim = c(0, 0), skip.boxplot = FALSE, skip.robustM7 = FALSE,
                     span.index = NULL) {
  blended <- FALSE
  if (isTRUE(attr(x, "type") == "diagnoseSeasonality")) {
    if (!is.null(x$spans)) { # If the model has multiple spans
      if (length(span.index) != 1) { # NULL or many
        x <- x$combinedX13
        blended <- TRUE
        if (length(span.index) > 1) warning("Requested more than 1 'span.index', ignoring them, plotting the combined result.")
      } else {
        if (span.index == "last") span.index <- length(x$spans)
        if (is.numeric(span.index)) {
          if (span.index %in% -1:0) span.index <- length(x$spans)
          if (span.index > length(x$spans)) {
            span.index <- length(x$spans)
            warning("'span.index' too large, plotting the last span instead.")
          }
          x <- x$spans[[span.index]]$seas
        } else stop("'span.index' should be NULL for the combined result, 'last', 0, or -1 for the last sample, or an integer span number.")
      }
    } else x <- x$seas
  }
  if (!("seas" %in% class(x))) stop("plotSeas takes a 'seasonal::seas' object or the output of 'diagnoseSeasonality' as input.")
  # If the supplied series coincide with the X13 output, ignore it
  if (isTRUE(all.equal(sa.custom, x$series$d11))) sa.custom <- NULL
  std.adj <- is.null(sa.custom) # Standard adjustment
  if (!std.adj & length(sa.custom) != length(x$series$d11)) stop("'sa.custom' must have the same length as the adjusted series in x.")

  s <- getSAStatOne(x, skip.robustM7 = skip.robustM7) # If (skip.robustM7), the time goes up by 0.2 s
  if (is.null(name)) name <- s$name # If no name is given

  # Extracting the series like in diagnoseSeasonalityOne
  x.orig   <- x$x
  xs <- x$series
  force <- !is.null(x$series$saa)
  if (!force) {
    x.sa     <- if (!std.adj) sa.custom else xs$d11
    x.adjfac <- xs$d16 # Seasonal and calendar effects
    x.seas   <- xs$d10
    x.cal <- if (s$log) x.adjfac / x.seas else x.adjfac - x.seas
  } else { # With forcing of annual totals, we redefine the SA component
    x.sa <- xs$saa # SA = adjusted and annual-corrected
    x.adjfac0 <- xs$d16 # AF without annual correction
    x.seas0 <- xs$d10 # Seasonal without annual correction
    x.adjfac <- if (s$log) x.orig / x.sa else x.orig - x.sa # Read adj. fac.
    x.cal  <- if (s$log) x.adjfac0 / x.seas0 else x.adjfac0 - x.seas0 # Vanilla calendar
    x.seas <- if (s$log) x.adjfac / x.cal else x.adjfac - x.cal
  }
  x.trend <- xs$d12
  x.ir <- xs$d13
  x.fc <- xs$fct
  x.bc <- xs$bct

  x.seas.factual <- x.seas # If a custom adjusted series is used, recall that orig = sa * adjfac and adjfac = seas * calend
  if (!std.adj) x.seas.factual <- if (s$log) x.orig / x.sa / x.cal else x.orig - x.sa - x.cal
  outliers <- seasonal::outlier(x) # Relies on seasonal:::extract_w_na_action, requires proper x$data
  sym.ts <- as.numeric(factor(outliers, levels = c("AO", "LS", "TC")))
  f <- stats::frequency(x.orig)
  do.calend <- !isTRUE(all.equal(diff(as.numeric(x.cal)), rep(0, length(x.cal)-1)))

  xr <- range(myTime(x.orig, numeric = TRUE)) + extend.xlim
  xlim <- c(floor(xr[1]), ceiling(xr[2])) # Integer years on the x axis
  xs5 <- (ceiling(xlim[1]/5):floor(xlim[2]/5))*5 # Getting round year numbers divisible by 5
  mycols <- c(orig = "#27ba30", sa = "#000000EE", trend = "#0056ED", seas = "#FF8512", cal = "#377EB8BB")
  lbcol <- "#FFFFFFB2" # Semi-transparent legend background
  lbxcol <- "#FFFFFF00" # Fully transparent legend box
  vxs <- setdiff(seq(round(xlim[1]), round(xlim[2])+1), xs5) # Secondary axis lines
  ylim <- range(x.sa, x.trend, x.orig, na.rm = TRUE)
  ylim <- ylim + c(-0.01, 0.01)*diff(ylim) # Extending slightly

  if (blended) { # If we have a seas modified by our chunk routine
    splits <- x$splits
    splits <- ind2time(splits, f) - 0.5/f
  } else splits <- NULL

  # Find two free corners for the legend
  xfi <- x.orig[is.finite(x.orig)]
  mdx <- stats::median(xfi)
  hdx <- mean(utils::head(xfi, f))
  tlx <- mean(utils::tail(xfi, f))
  leg.pos.l <- if (isTRUE(hdx > mdx)) "bottomleft" else "topleft"
  leg.pos.r <- if (isTRUE(tlx > mdx)) "bottomright" else "topright"

  graphics::layout(matrix(c(1, 1, 2, 2, 3, 4, 5, 6), ncol = 2, nrow = 4, byrow = TRUE), heights = c(1.5, 1, 1, 1))
  # Plot 1: series and adjusted
  withr::local_par(mar = c(0.2, 2, 0.2, 1.2))
  plot(NULL, NULL, xlim = xlim, ylim = ylim, bty = "n", xaxt = "n")
  graphics::abline(v = vxs - 0.5/f, col = "#0000003A", lty = 2)
  graphics::abline(v = xs5 - 0.5/f, col = "#0000003A", lty = 1, lwd = 1.5)
  if (blended) graphics::abline(v = splits, lty = 3, lwd = 1.5, col = "red")
  graphics::lines(x.orig, col = mycols["orig"], lwd = 2)
  if (!is.null(x.fc)) {
    graphics::lines(x.fc[, "forecast"], col = mycols["orig"], lwd = 1.5, lty = 2)
    graphics::lines(x.fc[, "upperci"], col = mycols["orig"], lwd = 1.5, lty = 3)
    graphics::lines(x.fc[, "lowerci"], col = mycols["orig"], lwd = 1.5, lty = 3)
  }
  if (!is.null(x.bc)) {
    graphics::lines(x.bc[, "backcast"], col = mycols["orig"], lwd = 1.5, lty = 2)
    graphics::lines(x.bc[, "upperci"], col = mycols["orig"], lwd = 1.5, lty = 3)
    graphics::lines(x.bc[, "lowerci"], col = mycols["orig"], lwd = 1.5, lty = 3)
  }
  graphics::lines(x.sa, col = "#FFFFFFAA", lwd = 4) # Adjusted halo
  graphics::lines(x.sa, col = mycols["sa"], lwd = 2.5)
  graphics::lines(x.trend + diff(ylim)*0.007, col = "#FFFFFFAA", lwd = 2, lty = 2) # Trend halo that does not obstruct the SA
  graphics::lines(x.trend - diff(ylim)*0.007, col = "#FFFFFFAA", lwd = 2, lty = 2)
  graphics::lines(x.trend, col = mycols["trend"], lwd = 2, lty = 2)
  graphics::legend(leg.pos.r, c("Original", "Adjusted", "Trend"), col = mycols[c("orig", "sa", "trend")],
                   lwd = c(2, 3, 2), lty = c(1, 1, 2), bg = lbcol, box.col = lbxcol)
  graphics::legend("top", paste0("Seasonality adjustment for ", name), bg = lbcol, box.col = lbxcol)
  if (blended) graphics::legend(leg.pos.l, "Split blend borders", lty = 3, lwd = 1.5, col = "red", bg = lbcol, box.col = lbxcol)
  graphics::points(x.orig, pch = as.numeric(sym.ts), col = "#FFFFFFEE", lwd = 5)
  graphics::points(x.orig, pch = as.numeric(sym.ts), col = "#FF0000", lwd = 2)
  # Labels at the proper side
  is.below.trend <- (x.sa < x.trend) & is.finite(x.sa)
  outl.pos <- !is.na(outliers)
  if (any(outl.pos)) labelsWithHalo(stats::time(x.orig)[outl.pos], y = x.orig[outl.pos], labels = outliers[outl.pos],
                               pos = ifelse(is.below.trend[outl.pos], 3, 1), cex = 0.75, offset = 0.4, hscale = 0.003, vscale = 0.01, nhalo = 16)

  # Plot 2: only seasonality
  # Find a less busy corner for the legend compared to the median
  xfi <- x.seas.factual[is.finite(x.seas.factual)]
  hdx <- utils::head(xfi, 3*f)
  mdx <- (stats::median(xfi) + if (s$log) 1 else 0) / 2
  hdl <- abs(stats::median(hdx[hdx < mdx]))
  if (!isTRUE(is.finite(hdl))) hdl <- Inf # If nothing is in the corner, good
  hdu <- abs(stats::median(hdx[hdx > mdx]))
  if (!isTRUE(is.finite(hdu))) hdu <- Inf
  leg.pos.l <- if (isTRUE(hdl > hdu)) "bottomleft" else "topleft"
  withr::local_par(mar = c(0.1, 2, 2, 1.2))
  ylim <- range(x.seas.factual, x.seas, x.cal, na.rm = TRUE)
  plot(NULL, NULL, bty = "n", xlim = xlim, ylim = ylim, xaxt = "n")
  graphics::abline(v = vxs - 0.5/f, col = "#0000003A", lty = 2)
  graphics::abline(v = xs5 - 0.5/f, col = "#0000003A", lty = 1, lwd = 1.5)
  if (blended) graphics::abline(v = splits, lty = 3, lwd = 1.5, col = "red")
  graphics::abline(h = if (s$log) 1 else 0, lty = 2)
  graphics::lines(x.seas.factual, lwd = if (force) 3.5 else 2.5, col = mycols["seas"])
  if (!std.adj) graphics::lines(x.seas, lty = 2)
  if (do.calend) {
    graphics::lines(x.cal, col = "#FFFFFFBB", lwd = 4) # Calendar halo
    graphics::lines(x.cal, col = mycols["cal"], lwd = 2)
  }
  if (force) graphics::lines(x.seas0, lty = 3)
  # Legend 1: adjustment type
  graphics::legend("topright", if (s$log) "Multiplicative" else "Additive", bg = lbcol, box.col = lbxcol)
  # Legend 2: line types
  leg.text <- c(if (std.adj) "Seasonal" else c("Factual seas.", "X13 seas."), if (force) "Seas. (no ann. force)" else NULL, if (do.calend) "Calendar" else NULL)
  leg.col <- c(mycols["seas"], if (!std.adj) "#000000" else NULL, if (force) "#000000" else NULL, if (do.calend) mycols["cal"] else NULL)
  leg.lwd <- c(if (force) 3.5 else 2.5, if (!std.adj) 1 else NULL, if (force) 1 else NULL, if (do.calend) 2 else NULL)
  leg.lty <- c(1, if (!std.adj) 2 else NULL, if (force) 3 else NULL, if (do.calend) 1 else NULL)
  graphics::legend("bottomright", legend = leg.text, ncol = length(leg.text), col = leg.col, lwd = leg.lwd, lty = leg.lty, bg = lbcol, box.col = lbxcol, text.width = NA)

  # Legend 3: quality statistics
  leg.text2 <- paste0(c("M7=", if (!skip.robustM7) "robM7=" else NULL, "Q2="), sprintf("%1.2f", c(s$M7, if (!skip.robustM7) s$robust.M7 else NULL, s$Q_M2)))
  graphics::legend(leg.pos.l, leg.text2, bg = lbcol, box.col = lbxcol, ncol = 2 + (!skip.robustM7), text.width = NA)
  if (!skip.robustM7) {
    leg.pos.l2 <- ifelse(leg.pos.l == "topleft", "bottomleft", "topleft")
    graphics::legend(leg.pos.l2, paste0(c("FS=", "FM="), c(s$f.seas, s$f.stab)), cex = 0.75, bg = lbcol, box.col = lbxcol)
  }
  axs <- xs5
  if (xlim[1] - min(xs5) > 1) axs <- c(min(axs)-5, axs) # Extending short axes
  if (xlim[2] - max(xs5) > 1) axs <- c(axs, max(axs)+5)
  graphics::axis(3, at = axs - 0.5/f, labels = axs)

  # Plot 3, 4: seasonal components
  withr::local_par(mar = c(2, 0.2, 0.2, 2))
  labs <- if (f == 4) paste0("Q", 1:4) else if (f == 12) month.abb else 1:f
  # The canvas is created in any case; the box plot lines can be hidden
  suppressWarnings(graphics::boxplot(x.orig ~ stats::cycle(x.orig), frame = FALSE, notch = TRUE, xlab = "", ylab = "", yaxt = "n",
                                     col = "#00000000", border = if (!skip.boxplot) "#858ae6" else "#00000000", names = substr(labs, 1, 2), pars = list(boxwex = 0.4, staplewex = 0.3, outwex = 0.3)))
  graphics::abline(h = stats::median(x.orig, na.rm = TRUE), lty = 2)
  stats::monthplot(x.orig, box = FALSE, add = TRUE, lwd = 2)
  graphics::legend("topleft", "Original", bg = lbcol, box.col = lbxcol)
  withr::local_par(mar = c(2, 0.2, 0.2, 0.2))
  suppressWarnings(graphics::boxplot(x.seas.factual ~ stats::cycle(x.seas.factual), frame = FALSE, notch = TRUE, xlab = "", ylab = "",
                                     col = "#00000000", border = if (!skip.boxplot) "#858ae6" else "#00000000", names = substr(labs, 1, 2), pars = list(boxwex = 0.4, staplewex = 0.3, outwex = 0.3)))
  graphics::abline(h = stats::median(x.seas.factual, na.rm = TRUE), lty = 2)
  stats::monthplot(x.seas.factual, box = FALSE, add = TRUE, lwd = 2)
  if (!std.adj) stats::monthplot(x.seas.factual, box = FALSE, add = TRUE, lty = 2, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  graphics::legend("topleft", "Seasonal", bg = lbcol, box.col = lbxcol)

  # Plots 5, 6: spectra before and after adjustment
  withr::local_par(mar = c(2, 0.2, 0.2, 0.2))
  computeSpectrum(x.orig, plot = TRUE, compact = TRUE)
  xspec <- computeSpectrum(x.sa, plot = TRUE, compact = TRUE, put.legend = FALSE)
  graphics::legend(attr(xspec, "legend.position"), c("Left: original", "Right: forced SA"), bg = lbcol, box.col = lbxcol)
  graphics::layout(matrix(1, ncol = 1))
}

.vHAC <- function(x) { # Fail-safe HAC with rule of thumb
  n <- length(x$residuals)
  rot <- round(0.75*n^(1/3))
  nw.lag <- tryCatch(sandwich::bwNeweyWest(x, prewhite = 0), error = function(e) {
    ret <- rot
    # warning("Too many dummies to proprely compute the HAC bw, using the rule-of-thumb.")
    return(ret)
  })
  if (nw.lag > length(x$residuals)/4) {
    # warning(paste0("The optimal Newey--West lag is too high (n=", n, ", lag=", round(nw.lag), "), using the rule-of-thumb (", rot, ")."))
    nw.lag <- rot
  }
  ret <- tryCatch(sandwich::kernHAC(x, prewhite = 0, kernel = "Bartlett", bw = nw.lag), error = .efv)
  if (is.null(ret)) ret <- .vHC(x)
  attr(ret, "bandwidth") <- if (isTRUE(attr(nw.lag, "fail"))) NA else nw.lag
  return(ret)
}
.vHC <- function(x) {
  ret <- tryCatch(sandwich::vcovHC(x, type = "HC0"), error = .efv) # Fallback VCOV if vHAC fails
  if (is.null(ret)) ret <- matrix(NA, nrow = length(x$coefficients), ncol = length(x$coefficients))
  return(ret)
}
.tryTest <- function(x, what, robust = TRUE) {
  failRet <- function(e) return(list("Hypothesis testing failed.", F = rep(NA, 2), `Pr(>F)` = rep(NA, 2)))
  if (robust) {
    vHA <- .vHAC(x)
    ret <- tryCatch(car::linearHypothesis(x, what, vcov. = vHA), error = .ef)
    # If it failed, either vHAC failed or vHA is not invertible
    if (is.null(ret)) ret <- tryCatch(car::linearHypothesis(x, what, vcov. = .vHC), error = .ef)
    if (is.null(ret)) ret <- tryCatch(car::linearHypothesis(x, what), error = failRet)
    attr(ret, "bandwidth") <- attr(vHA, "bandwidth")
  } else ret <- tryCatch(car::linearHypothesis(x, what), error = failRet)
  if (is.null(attr(ret, "bandwidth"))) attr(ret, "bandwidth") <- NA
  return(ret)
}

#' Robustified identifiable and stable seasonality tests for X13 M7
#'
#' @param x An object of class "seas" returned by [seasonal::seas]
#' @param type Character: "rlm"for a robust linear model (default), "lm" for a linear model,
#' "rank" for the rank regression, or "nonrobust" for the homoskedasticity-based ANOVA.
#'
#' This function provides a long-overdue improvement for the X13 ANOVA-based tests
#' due to Higgins (1975)
#' by using (1) weighted robust regression instead of the rank regression and
#' (2) heteroskedasticity-and-autocorrelation-consistent VCOV instead of the
#' unreaslistic white-noise assumption that almost never holds for the irregular
#' component.
#'
#' The Kruskal--Wallis and Friendman are non-robust to the violation of the
#' homoskedasticity assumption. Insted, we run two simple linear regressions
#' and test the joint significance of the dummies of interest.
#'
#' Define SI as the sum of the seasonal and the irregular components.
#'
#' 1. Identifiable seasonality: regress SI on the seasonal dummies
#  2. Stable seasonality: regress SI on the year dummies and period dummies
#'
#' If some periods exhibit huge seasonal effects and/or atypical values, then,
#' the rank regression (`type = "rank"`) would ignore their magnitude that could
#' potentially translate into stronger #' estimated seasonality effect.
#' If `type = "rlm"`, then, IRLS via [MASS::rlm()] reduces the impact of such
#' outliers to some extent.
#' HAC-robust ANOVA is conducted in all cases. HAC estimation may fails if
#' there are singular year dummies. In this case, such singular observations
#' are dropped. If still, the model matrix is singular, it reverts to
#' non-robust VCOV estimation with a warning.
#'
#' If IRLS does not converge, the fallback is the rank regression with HAC inference.
#' If HAC cannot be computed, the fallback is the HC0 and, should it fails as well,
#' the non-robust VCOV (although in practice, it should never be the case).
#'
#' @return A list with 4 components: `models` (the list of 2 models), `tests` (the list of 2 ANOVAs),
#' `robust.m7` (the final robust M7 statistic), and the HAC bandwidth used in VCOV estimation; if it failed
#' due to singularities, bandwidth is `NA`.
#' @export
#'
#' @examples
#' x1 <- AirPassengers
#' x2 <- window(x1, start = c(1949, 12), end = c(1960, 1)) # With singular years
#' sa1 <- seasonal::seas(x1, x11 = "")
#' sa2 <- diagnoseSeasonality(x2)
#' # getSAStat(list(sa1, sa2))
#' robustSeasTests(sa1)
#' robustSeasTests(sa2$seas)
robustSeasTests <- function(x, type = c("rlm", "lm", "rank", "nonrobust")) {
  type <- type[1]
  if ((!type %in% c("rlm", "lm", "rank", "nonrobust"))) stop("'type' must be 'rlm', 'lm', 'rank', or 'nonrobust'.")
  if (!("seas" %in% class(x))) stop("robustSeasTests takes a 'seasonal::seas' object as input.")

  do.log <- isTRUE(seasonal::transformfunction(x) == "log")
  si <- if (do.log) x$series$d10 * x$series$d13 else x$series$d10 + x$series$d13
  cyc <- factor(stats::cycle(si))
  yr <- factor(myYear(si))

  lm.seas <- switch(type, rlm = tryCatch(MASS::rlm(si ~ cyc, maxit = 200), error = .ef),
                    lm = stats::lm(si ~ cyc),
                    rank = stats::lm(rank(si) ~ cyc),
                    nonrobust = stats::lm(rank(si) ~ cyc))
  if ("rlm" %in% class(lm.seas)) lm.seas <- stats::lm(stats::formula(lm.seas), data = stats::model.frame(lm.seas), weights = lm.seas$w) # Recasting to an identical lm
  # This circumvents a bug in sandwich:::bread.rlm due to the zero derivative of Huber psi
  if (is.null(lm.seas)) lm.seas <- stats::lm(rank(si) ~ cyc)
  lm.stab <- switch(type, rlm = tryCatch(MASS::rlm(si ~ cyc + yr, maxit = 200), error = .ef),
                    lm = stats::lm(si ~ cyc + yr),
                    rank = stats::lm(rank(si) ~ cyc + yr),
                    nonrobust = stats::lm(rank(si) ~ cyc + yr))
  if ("rlm" %in% class(lm.stab)) lm.stab <- stats::lm(stats::formula(lm.stab), data = stats::model.frame(lm.stab), weights = lm.stab$w)
  if (is.null(lm.stab)) lm.stab <- stats::lm(rank(si) ~ cyc + yr)

  test.seas <- .tryTest(lm.seas, names(stats::coef(lm.seas))[-1], robust = (type != "nonrobust"))
  hyp.stab <- names(stats::coef(lm.stab))[grep("^yr", names(stats::coef(lm.stab)))]
  test.stab <- .tryTest(lm.stab, hyp.stab, robust = (type != "nonrobust"))
  attr(test.seas, "heading") <- c("Seasonalily model: Y_ij = a_i + U_ij", "(ANOVA for seasonal dummies)", attr(test.seas, "heading"))
  attr(test.stab, "heading") <- c("Evolutive seasonalily model: Y_ij = a_i + b_j + U_ij", "(ANOVA in the spirit of Higgins (1975))", attr(test.stab, "heading"))

  pval <- c(identifiable = test.seas$`Pr(>F)`[2], stable = test.stab$`Pr(>F)`[2])

  robust.m7 <- sqrt((1.5*test.stab$F[2] + 3.5) / test.seas$F[2])
  if (robust.m7 > 3) robust.m7 <- 3
  robust.m7 <- round(robust.m7, 3)
  bw <- c(attr(test.seas, "bandwidth"), attr(test.stab, "bandwidth"))
  attr(bw, "kernel") <- "Bartlett"
  return(list(models = list(identifiable = lm.seas, stable = lm.stab),
              tests = list(identifiable = test.seas, stable = test.stab),
              robust.m7 = robust.m7, pval = pval, HACbw = bw))
}

# .outlier2ind(c("ao2008.10", "tc2020.mar"))
.outlier2ind <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.character(x)) stop("'x' must be a character of X13-compatible outlier entries, e.g. 'ao1975.aug' or 'ao2008.3'.")
  for (i in 1:length(month.abb)) x <- gsub(month.abb[i], as.character(i), x, ignore.case = TRUE)
  x <- gsub("[A-Za-z]", "", x)
  if (!all(gsub("[0-9.]", "", x) == "")) stop("'Check the outlier names: it should be TYPEyear.period, e.g. 'tc2020.mar' or 'ao2008.3'.")
  xl <- strsplit(x, ".", fixed = TRUE)
  xi <- do.call(rbind, lapply(xl, as.integer))
  return(xi)
}

.p0 <- function(x, collapse = ".") paste0(x, collapse = collapse)

#' @rdname diagnoseSeasonality
#' @param nbcast A non-negative integer: how many periods to backcast before the first valid observation (only for `diagnoseSeasonalityOne`)
#' @param nfcast A non-negative integer: how many periods to forecast after the last valid observation (only for `diagnoseSeasonalityOne`)
diagnoseSeasonalityOne <- function(x, calendar = NULL, name = "(noname)",
                                 nbcast = 0, nfcast = 0,
                                 transform = c("auto", "yes", "no"),
                                 td = c("auto", "6", "1", "no"),
                                 leap.year = c("auto", "yes", "no"),
                                 easter = c("auto", "yes", "no"),
                                 forced.outliers = NULL,
                                 m7.threshold = 1,
                                 force.annual = c("no", "regress", "denton"),
                                 sa.rule = c("yes", "robustM7", "M7", "averageM7", "no"),
                                 transform.aicdiff = -2, tradingdays.aicdiff = 0,
                                 plot.file = NULL, skip.boxplot = FALSE,
                                 verbose = 2
                                ) {
  # Loading the functions that are not exported
  update.seas <- utils::getFromNamespace("update.seas", "seasonal")
  predict.seas <- utils::getFromNamespace("predict.seas", "seasonal")

  # Checking the inputs
  if (!stats::is.ts(x)) stop("'x' must be a time series with frequency 4 or 12.")
  if (!is.null(dim(x))) stop("'x' has multiple series. Use 'diagnoseSeasonality(x, ...)' instead to handle many series.")

  transform <- transform[1] # Changing the input names for X13
  if (is.null(transform) | isTRUE(is.na(transform))) transform <- "auto"
  if (transform == "log" | isTRUE(transform)) transform <- "yes"
  if (transform == "none" | isFALSE(transform)) transform <- "no"
  if (!(transform %in% c("yes", "no", "auto"))) stop("'transform' must be 'yes', 'no', or 'auto'.")

  td <- td[1]
  if (is.null(td) | isTRUE(is.na(td))) td <- "auto"
  if (isTRUE(td) | td == "yes") stop("'td' cannot be simply 'yes' or TRUE. Use 6 (for separate day effects) or 1 (for one weekday effect).")
  if (td == 0 | isFALSE(td) | td == "none") td <- "no"
  if (is.numeric(td)) td <- as.character(td)
  if (!(td %in% c("auto", "6", "1", "no"))) stop("'td' must be be 'auto', '6', '1' or 'no'.")

  leap.year <- leap.year[1]
  if (is.null(leap.year) | isTRUE(is.na(leap.year))) leap.year <- "auto"
  if (isTRUE(leap.year)) leap.year <- "yes"
  if (isFALSE(leap.year)) leap.year <- "no"
  if (!(leap.year %in% c("yes", "no", "auto"))) stop("'leap.year' must be 'yes', 'no', or 'auto'.")

  easter <- easter[1]
  if (is.null(easter) | isTRUE(is.na(easter))) easter <- "auto"
  if (isTRUE(easter)) easter <- "yes"
  if (isFALSE(easter)) easter <- "no"
  if (!(easter %in% c("yes", "no", "auto"))) stop("'easter' must be 'yes', 'no', or 'auto'.")

  force.annual <- force.annual[1]
  if (force.annual == "no") force.annual <- "none"
  force <- force.annual != "none"
  if (!(force.annual %in% c("none", "regress", "denton"))) stop("'force.annual' must be 'none', 'regress', or 'denton'.")

  sa.rule <- sa.rule[1]

  freq <- stats::frequency(x)
  if (!(freq %in% c(4, 12))) stop("The data are not monthly or quarterly, aborting.")
  if (!is.null(calendar)) {
    if (is.character(calendar)) {
      all.cals <- utils::data(package = "StatecCMP")$results[, "Item"]
      dname <- paste0(calendar, ".", if (freq == 4) "Q" else "M")
      if (!(dname %in% all.cals)) {
        all.countries <- unique(gsub("\\.[QM]", "", all.cals))
        stop(paste0("Wrong calendar name supplied (", calendar, "). Available ones: ",
                    paste0(all.cals, collapse = ", "), ".\nSupply a calendar from JDemetra+ or just remove the 'calendar' argument."))
      }
      utils::data(list = dname, package = "StatecCMP", envir = environment())
      assign("calendar", get(dname))
    }
    cnames <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "LeapYear", "WorkingDays")
    if (!all(colnames(calendar) %in% cnames)) stop(paste0("The calendar must be an 'mts' with the following series names:\n", paste0(cnames, collapse = ", "), "."))
    if (freq != stats::frequency(calendar)) stop("The frequency of the data is not equal to that of the calendar, aborting.")
    if (verbose > 1) cat("Using the user-supplied national calendar.\n")
  } else {
    warning("Using the default X13 Census calendar. Better supply a country-specific calendar.")
  }

  full.range <- getRangeVec(x)
  x.ind <- time2ind(x)
  # Not accepting NAs at the tails; trim with a warning
  if ((full.range["first"] != 1) | (full.range["last"] != length(x))) {
    warning(paste0("The series ", name, " has missing observations at the ends, dropping.\nIf you want to impute them by forecasting, use 'diagnoseSeasonality' instead."))
    x <- stats::window(x, start = x.ind[full.range["first"], ], end = x.ind[full.range["last"], ])
    full.range <- getRangeVec(x)
    x.ind <- time2ind(x)
  }
  yr <- myYear(x)
  cyc <- stats::cycle(x)

  n <- length(x)
  if (n / freq > 15) {
    wrn.msg <- paste0("The series ", name, " are very long (>", round(n / freq), " years), estimation may fail.\n! Better use 'diagnoseSeasonality(..., max.length = 15)'.")
    cat("! ", wrn.msg, "\n", sep = "")
    warning(wrn.msg)
  }
  linelen <- 60
  if (verbose > 1) printSym("=", linelen)
  if (verbose > 1) cat("Testing and adjusting seasonality in", name, "via the recommended procedure with parameter auto-selection.\n") else
    if (verbose > 0) cat("Seasonality adjustment for ", name, " via X13.\n", sep = "")
  if (verbose > 1) printSym("=", linelen)

  do.bcast <- if (nbcast > 0) "yes" else "no"
  do.fcast <- if (nfcast > 0) "yes" else "no"

  # JDemetra-like procedure with all 5 weekday dummies + Easter + leap year
  # unless the user has overridden some of the parameters
  ly <- if (leap.year != "no") "lpyear" else NULL
  ea <- if (easter != "no") "easter[8]" else NULL
  if (!is.null(calendar)) {
    # calendar <- stats::window(calendar, start = min(stats::time(x)) - head.years, end = max(stats::time(x)) + tail.years)
    td.vars <- wd.vars <- c(ea, ly)
    # These following variables must be present in the calendar for Luxembourg
    td.cal <- calendar[, c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")]
    wd.cal <- calendar[, "WorkingDays", drop = FALSE]
    td.types <- rep("td", 6)
    wd.types <- "td"
  } else {
    td.vars <- c(ea, ly, "tdnolpyear")
    wd.vars <- c(ea, ly, "td1nolpyear")
    td.cal <- wd.cal <- td.types <- wd.types <- NULL
  }
  nd.vars <- c(ea, ly)
  if (verbose > 0) cat("Start: ", .p0(x.ind[1, ]), ", end: ", .p0(x.ind[n, ]), ", (bcst, fcst) = (", nbcast, ", ", nfcast, ")",
                       if (force) " [forced annual totals]" else NULL, "\n", sep = "")

  # Checking which of the forced outliers belong to the sample
  fo.inds <- .outlier2ind(forced.outliers)
  fo.time <- ind2time(fo.inds, freq)
  impossible.outliers <- fo.time < ind2time(x.ind[1, ], freq) | fo.time > ind2time(x.ind[n, ], freq)
  forced.outliers <- forced.outliers[!impossible.outliers]
  # Deleting the zero after dot because the user may write something like ao2020.03,
  # which is conceptually the same as ao2020.3, but the string is different
  forced.outliers <- gsub("\\.0([1-9])$", ".\\1", forced.outliers)
  # Adding forced outliers as the regressors that should always be honoured
  td.vars <- c(td.vars, forced.outliers)
  wd.vars <- c(wd.vars, forced.outliers)
  nd.vars <- c(nd.vars, forced.outliers)

  miss.inds <- which(!is.finite(x)) # Used here and later to impute NAs
  if (full.range["mid.na.count"] > 0) {
    cat("The series contain internal NAs. Plugging in interpolated values and adding AOs for each point.\n")
    good.inds <- which(is.finite(x))
    x.interp <- stats::approx(x = good.inds, y = as.numeric(x)[good.inds], xout = seq_along(x))$y
    x[1:length(x)] <- x.interp
    miss.vars <- paste0("ao", yr, ".", cyc)[miss.inds]
    td.vars <- c(td.vars, miss.vars)
    wd.vars <- c(wd.vars, miss.vars)
    nd.vars <- c(nd.vars, miss.vars)
  }
  cal.names <- c("easter[8]", "lpyear", "const", "tdnolpyear", "td1coef", "td1nolpyear")
  # These are the regressor names that are not outliers and should be dropped

  # A back-up function
  fallback <- function() tryCatch(seasonal::seas(x = x, x11 = "", regression.variables = NULL,
                          automdl.maxorder = c(2, 1), transform.function = "auto",
                          x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                          x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                          forecast.save = c("bct", "fct"), spectrum.savelog = "all"), error = .efv)
  flbk.msg <- "Even the basic auto-model could not be estimated. The series cannot be adjusted with X13.\nPlot the series -- most likely is it piecewise linear, although other irregularities are possible."

  ######################################################################
  # Step 1: choosing between 2 transformations using the common outliers
  # The rest will be similar to the first call, so creating a new function
  # The default is all calendar regressors
  transf.rv <- switch(td, auto = td.vars, `6` = td.vars, `1` = wd.vars, no = nd.vars)
  transf.xreg <- switch(td, auto = td.cal, `6` = td.cal, `1` = wd.cal, no = NULL)
  transf.ut <- switch(td, auto = td.types, `6` = td.types, `1` = wd.types, no = NULL)
  estSeas1 <- function(tr) {
    dlist <- list(x = x, x11 = "", outlier.types = c("ao", "ls", "tc"),
                  regression.variables = transf.rv, regression.aictest = NULL,
                  xreg = transf.xreg, regression.usertype = transf.ut, # Depends on `td`
                  transform.function = tr, ####### The only manual argument
                  force.type = force.annual,
                  automdl.maxorder = c(3, 1), estimate.tol = 1e-7, estimate.maxiter = 5000,
                  x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                  x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                  forecast.save = c("bct", "fct"), spectrum.savelog = "all") # Default argument list
    ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    if (is.null(ret)) {
      dlist$estimate.maxiter <- 20000
      message(paste0("SA of ", name, " not converged -- increasing maxiter from 5,000 to 20,000."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-6
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-7 to 1e-6."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-5
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-6 to 1e-5."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-4
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-5 to 1e-4."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$outlier.types = c("ao", "ls")
      message(paste0("SA of ", name, " not converged -- disallowing auto-detection of TC outliers."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$outlier.types = "ao"
      message(paste0("SA of ", name, " not converged -- disallowing LS outliers (only AO left)."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$automdl.maxorder = c(2, 1)
      message(paste0("SA of ", name, " not converged -- reducing the maximum order from (3 1) to (2 1)."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    return(ret)
  }

  trans.fun <- switch(transform, auto = "auto", yes = "log", no = "none")
  if (any(x <= 0) & transform != "none") {
    if (transform == "yes") warning("The series contains non-positive values. Multplicative adjustment impossible, switching to additive.")
    trans.fun <- "none"
  }
  ntd <- switch(td, auto = 6, `6` = 6, `1` = 1, no = 0)
  oreg <- paste0(c("TD", ly, ea), collapse = " + ")
  msg.mult <- paste0("-- Model: multipl., ", ntd, " ", oreg)
  msg.add  <- paste0("-- Model: additive, ", ntd, " ", oreg)
  if (trans.fun != "auto") {
    a.td <- estSeas1(tr = trans.fun)
    if (verbose > 1) .printSeasAIC(a.td, if (trans.fun == "log") msg.mult  else msg.add)
  } else {
    if (verbose > 0) cat("Comparing transformation models...\n")
    a.td.log <- estSeas1(tr = "log")
    if (verbose > 1) .printSeasAIC(a.td.log, msg.mult)
    a.td.nolog <- estSeas1(tr = "none")
    if (verbose > 1) .printSeasAIC(a.td.nolog, msg.add)

    # Extracting the common outliers
    extrOut <- function(y) tryCatch(tolower(seasonal::outlier(y)), error = function(e) {ret <- rep(NA, length(x)); attributes(ret) <- attributes(x); ret})
    log.outlier.ts <- data.frame(log = extrOut(a.td.log), nolog = extrOut(a.td.nolog))
    log.outlier.ts.med <- lapply(1:nrow(log.outlier.ts), function(i) .OTMed2(as.matrix(log.outlier.ts)[i, ])) # Creating only one outlier per period
    log.outlier.ts.med.conc <- unlist(lapply(log.outlier.ts.med, attr, "concordant"))
    log.outlier.ts.med <- unlist(log.outlier.ts.med)
    log.outl.time <- paste0(yr, ".", cyc)
    log.outlier.tab <- paste0(log.outlier.ts.med, log.outl.time)
    log.outl.reg <- log.outlier.tab[!grepl("^NA", log.outlier.tab)]

    if (length(log.outl.reg) > 0 & any(!log.outlier.ts.med.conc)) { # There are some outlier regressors and at least one is not concordant across the two models
      if (verbose > 1) cat("Outlier regressors used in the transformation choice: ", paste0(log.outl.reg, collapse = ", "), "\n", sep = "")
      if (verbose > 1) cat("Estimation with identical outlier regressors for comparability by AICc...\n")
      # If re-estimation is successful, then, choose the best model based on auto-AIC
      a.td.auto <- tryCatch(update.seas(a.td.nolog, regression.variables = unique(c(td.vars, log.outl.reg)), outlier = NULL, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = .efv)
    } else {
      if (verbose > 1) cat("All outlier variables are identical under 2 tranformations.\n", sep = "")
      a.td.auto <- tryCatch(update.seas(a.td.nolog, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = .efv)
    }

    if (!is.null(a.td.auto)) { # If X13-based auto-selection converged
      if (length(log.outl.reg) > 0 & any(!log.outlier.ts.med.conc)) {
        if (verbose > 1) cat("Estimation with common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
      }
      if (!is.null(a.td.auto$udg["aictest.trans.aicc.nolog"])) {
        if (verbose > 1) cat("-- Model: automatic X13 transformation, AICc(log) = ",
                             sprintf("%1.1f", as.numeric(a.td.auto$udg["aictest.trans.aicc.log"])),
                             ", AICc(nolog) = ", sprintf("%1.1f", as.numeric(a.td.auto$udg["aictest.trans.aicc.nolog"])), ".\n", sep = "")
      }
      # Choosing the transformation
      trans.fun <- seasonal::transformfunction(a.td.auto)
    } else { # Something went wrong with the same set of outliers -- possibly too many -- reduce the number
      if (verbose > 0) cat("! Re-trying with fewer outliers (existing in both models).\n")
      log.outlier.tab <- table(c(a.td.log$model$regression$variables, a.td.nolog$model$regression$variables))
      log.outlier.tab <- log.outlier.tab[log.outlier.tab == 2]
      log.outl.reg <- sort(setdiff(names(log.outlier.tab), cal.names))
      if (verbose > 1) cat("Revised outlier regressors used in the transformation choice: ", paste0(log.outl.reg, collapse = ", "), "\n", sep = "")
      a.td.auto <- tryCatch(update.seas(a.td.nolog, regression.variables = unique(c(td.vars, log.outl.reg)), outlier = NULL, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = .efv)
      if (!is.null(a.td.auto)) {
        if (verbose > 1) cat("Re-estimation with fewer common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
        trans.fun <- seasonal::transformfunction(a.td.auto)
      } else {
        if (verbose > 0) cat("! Re-estimation with fewer common outliers was unsuccessful.\n! Doing selection with the original models.\n! The results might be sub-optimal because the outliers are different.\n")
        log.mod.list <- list(a.td.log, a.td.nolog)
        aiccs <- unname(sapply(log.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
        names(aiccs) <- c("log", "nolog")
        eff.len <- unname(sapply(log.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 1))
        # Adjusting the AICc based on the number of observations
        aiccs <- aiccs / eff.len * n
        if (any(is.finite(aiccs))) { # At least one model was estimated in overall
          trans.fun <- if (aiccs["nolog"] - aiccs["log"] < transform.aicdiff) "none" else "log"
        } else trans.fun <- "none" # Just a plug because nothing was estimated
      }
    }
    a.td <- if (trans.fun == "log") a.td.log else a.td.nolog
  }

  if (is.null(a.td)) {
    warning("The transformation selection failed with default (generous) settings. Falling back to the most parsimonious settings and retrying.")
    a.td <- fallback()
    if (is.null(a.td)) {
      warning(flbk.msg)
      return(NULL)
    } else {
      trans.fun <- seasonal::transformfunction(a.td)
    }
  }
  if (verbose > 0) cat("*** The adjustment model is ", if (trans.fun == "log") "MULTIPLICATIVE" else "ADDITIVE", ".\n", sep = "")
  if (verbose > 1) printSym("-", linelen)
  if (verbose == 1 & td == "auto") cat("Estimating 3 main models (6 TD, 1 TD, 0 TD)...\n")

  #####################################################################
  # Step 2: choosing calendar regressors
  estSeas2 <- function(regvar, xreg = NULL, usertype = NULL, aictest = NULL) {
    dlist <- list(x = x, x11 = "", outlier.types = c("ao", "ls", "tc"),
                  regression.variables = regvar, regression.aictest = aictest,
                  xreg = xreg, regression.usertype = usertype,
                  transform.function = trans.fun,
                  force.type = force.annual,
                  automdl.maxorder = c(3, 1), estimate.tol = 1e-7, estimate.maxiter = 5000,
                  x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                  x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                  forecast.save = c("bct", "fct"), spectrum.savelog = "all") # Default argument list
    if (is.null(usertype)) dlist$usertype <- NULL # Dropping the arguments
    if (is.null(xreg)) dlist$xreg <- NULL
    ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    if (is.null(ret)) {
      dlist$estimate.maxiter <- 20000
      message(paste0("SA of ", name, " not converged -- increasing maxiter from 5,000 to 20,000."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-6
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-7 to 1e-6."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-5
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-6 to 1e-5."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$estimate.tol <- 1e-4
      message(paste0("SA of ", name, " not converged -- relaxing the tolerance from 1e-5 to 1e-4."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$outlier.types = c("ao", "ls")
      message(paste0("SA of ", name, " not converged -- disallowing auto-detection of TC outliers."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$outlier.types = "ao"
      message(paste0("SA of ", name, " not converged -- disallowing LS outliers (only AO left)."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    if (is.null(ret)) {
      dlist$automdl.maxorder = c(2, 1)
      message(paste0("SA of ", name, " not converged -- reducing the maximum order from (3 1) to (2 1)."))
      ret <- tryCatch(seasonal::seas(list = dlist), error = .efv)
    }
    return(ret)
  }
  if (td == "auto") {
    # The previous step estimated a 6TD model--estimating the other two
    if (verbose > 1) .printSeasAIC(a.td, "-- Model: 6 TD")
    a.wd <- estSeas2(regvar = wd.vars, xreg = wd.cal, usertype = wd.types)
    if (verbose > 1) .printSeasAIC(a.wd, "-- Model: 1 TD")
    a.nd <- estSeas2(regvar = nd.vars)
    if (verbose > 1) .printSeasAIC(a.nd, "-- Model: 0 TD")

    extrOut <- function(y) tryCatch(tolower(seasonal::outlier(y)), error = function(e) {ret <- rep(NA, length(x)); attributes(ret) <- attributes(x); ret})
    outlier.ts <- data.frame(td = extrOut(a.td), wd = extrOut(a.wd), nd = extrOut(a.nd))
    outlier.ts.med <- lapply(1:nrow(outlier.ts), function(i) .OTMed3(as.matrix(outlier.ts)[i, ])) # Creating only one outlier per period
    outlier.ts.med.conc <- unlist(lapply(outlier.ts.med, attr, "concordant"))
    outlier.ts.med <- unlist(outlier.ts.med)
    outl.time <- paste0(yr, ".", cyc)
    outlier.tab <- paste0(outlier.ts.med, outl.time)
    outl.reg <- outlier.tab[!grepl("^NA", outlier.tab)]

    if (length(outl.reg) > 0 & any(!outlier.ts.med.conc)) { # There are some outlier regressors and at least one is not concordant across models
      if (verbose > 1) cat("Outlier regressors used in testing: ", paste0(outl.reg, collapse = ", "), "\n", sep = "")
      if (verbose > 1) cat("Estimating models with identical outlier regressors to make the TD regressors comparable by AICc...\n")
      a.td2 <- tryCatch(update.seas(a.td, regression.variables = unique(c(td.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.td2, "-- Refit: 6 TD + common OL")
      a.wd2 <- tryCatch(update.seas(a.wd, regression.variables = unique(c(wd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.wd2, "-- Refit: 1 TD + common OL")
      a.nd2 <- tryCatch(update.seas(a.nd, regression.variables = unique(c(nd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.nd2, "-- Refit: 0 TD + common OL")
    } else {
      if (verbose > 1) cat("All outlier variables are identical in 3 models, skipping re-estimation.\n", sep = "")
      a.td2 <- a.td
      a.wd2 <- a.wd
      a.nd2 <- a.nd
    }

    # If re-estimation is successful, then, choose the best model based on AIC among the comparable, and then, go back to the original
    if (!is.null(a.td2) & !is.null(a.wd2) & !is.null(a.nd2)) {
      if (length(outl.reg) > 0 & any(!outlier.ts.med.conc)) {
        if (verbose > 1) cat("Estimation with common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
      }
      td.mod.list <- list(a.td2, a.wd2, a.nd2)
      aiccs <- unname(sapply(td.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
      eff.len <- unname(sapply(td.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 1))
      aiccs <- aiccs / eff.len * n
    } else { # Something went wrong with the same set of outliers---possibly too many---do the quick and dirty selection with the original models
      if (verbose > 0) cat("! Estimation with common outliers failed for at least 1 model.\n")
      if (verbose > 0) cat("! Re-trying with fewer outliers (existing in at least 2 out of 3 models).\n")
      outlier.tab <- table(c(a.td$model$regression$variables, a.wd$model$regression$variables, a.nd$model$regression$variables))
      outlier.tab <- outlier.tab[outlier.tab >= 2]
      outl.reg <- sort(setdiff(names(outlier.tab), cal.names))
      if (verbose > 1) cat("Outlier regressors used in testing: ", paste0(outl.reg, collapse = ", "), "\n", sep = "")
      a.td2 <- tryCatch(update.seas(a.td, regression.variables = unique(c(td.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.td2, "-- Refit: 6 TD + common OL")
      a.wd2 <- tryCatch(update.seas(a.wd, regression.variables = unique(c(wd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.wd2, "-- Refit: 1 TD + common OL")
      a.nd2 <- tryCatch(update.seas(a.nd, regression.variables = unique(c(nd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = .efv)
      if (verbose > 1) .printSeasAIC(a.nd2, "-- Refit: 0 TD + common OL")

      if (!is.null(a.td2) & !is.null(a.wd2) & !is.null(a.nd2)) {
        if (verbose > 1) cat("Re-estimation with fewer common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
        td.mod.list <- list(a.td2, a.wd2, a.nd2)
      } else { # At least one model with fewer outliers did not converge
        if (verbose > 0) cat("! Re-estimation with fewer common outliers was unsuccessful.\n! Doing selection with the original models, but the results might be misleading because the outliers are different.\n")
        td.mod.list <- list(a.td, a.wd, a.nd)
      }
      aiccs <- unname(sapply(td.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
      eff.len <- unname(sapply(td.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 1))
      aiccs <- aiccs / eff.len * n
    }
    names(aiccs) <- c("td", "wd", "nd")
    if (any(is.finite(aiccs))) { # At least one model was estimated in overall
      best.type <- "nd" # Start with the most parsimonious one
      # If ND is unavailable or (WD is better and the comparison makes sense), choose WD
      if (isTRUE(aiccs["wd"] - aiccs["nd"] < tradingdays.aicdiff) | !is.finite(aiccs["nd"])) best.type <- "wd"
      # If the new.best is unavailable or (new.best is better and the comparison makes sense), choose TD
      if (isTRUE(aiccs["td"] - aiccs[best.type] < tradingdays.aicdiff) | !is.finite(aiccs[best.type])) best.type <- "td"
      if (verbose > 0) switch(best.type, td = cat("*** The week-day effects (6 TD) are different and important.\n"),
                              wd = cat("*** The number of week days (1 TD) is important.\n"), nd = cat("*** The week-day effect is absent in these series (0 TD).\n"))
      a.sel <- switch(best.type, td = a.td, wd = a.wd, nd = a.nd)
    } else {
      warning("No model with the recommended settings could be estimated. Trying with the simplest default X13 settings.")
      best.type <- "nd"
      a.sel <- fallback()
      if (is.null(a.sel)) {
        warning(flbk.msg)
        return(NULL)
      }
    }
  } else {
    # The model has already been estimated at the first step
    best.type <- switch(td, `6` = "td", `1` = "wd", no = "nd")
    a.sel <- a.td
    msg2 <- switch(best.type, td = "*** 6 TD (user override)", wd = "*** 1 TD (user override)", nd = "*** 0 TD (user override)")
    if (verbose > 1) .printSeasAIC(a.sel, msg2) else if (verbose > 0) cat(msg2, "\n")
  }

  xregs <- a.sel$model$regression$variables
  if (verbose > 1) cat("Outlier regressors in the model with forced Easter and LY effects: ", paste0(sort(setdiff(xregs, cal.names)), collapse = ", "), "\n", sep = "")
  rv <- switch(best.type, td = td.vars, wd = wd.vars, nd = nd.vars)
  xr <- switch(best.type, td = td.cal, wd = wd.cal, nd = NULL)
  ut <- switch(best.type, td = td.types, wd = wd.types, nd = NULL)

  #######################################################################
  # Step 3: Easter and Leap year
  # If the leap year and Easter are either "yes" or "no", no selection should be done
  # Otherwise, test whilst keeping the same outliers and calendar regressors
  # If "auto", they were included by default
  to.test <- NULL
  if (leap.year == "auto") to.test <- c(to.test, "lpyear")
  if (easter == "auto") to.test <- c(to.test, "easter")

  if (length(to.test) > 0) {
    a.ealy <- tryCatch(update.seas(a.sel, regression.aictest = to.test, transform.function = trans.fun, regression.variables = xregs, outlier = NULL), error = .efv)
    if (!is.null(a.ealy)) { # If the model without Easter/LY converges
      if (!isTRUE(all.equal(sort(names(a.ealy$est$coefficients)), sort(names(a.sel$est$coefficients))))) {
        # Easter/LY re-estimation succeeded and something changed
        if (verbose > 1) .printSeasAIC(a.ealy, paste0("-- Refit: auto-(", paste0(to.test, collapse = ", "), ") + common OL"))
        if (verbose > 1) cat("Easter and LY effect testing yielded a different model, re-estimating with auto-selection...\n", sep = "")
      } else {
        if (verbose > 1) cat("Easter and LY effect test changed nothing, re-estimating with auto-selection...\n")
      }
      # Then, we re-estimate the specification (no forced outlier xregs)
      a.ealy2 <- estSeas2(regvar = rv, aictest = to.test, xreg = xr, usertype = ut)
      if (!is.null(a.ealy2)) {
        a.sel <- a.ealy2
        if (verbose > 1) .printSeasAIC(a.ealy2, paste0("-- Refit: auto-(", paste0(c(to.test, "outliers"), collapse = ", "), ")"))
      } else {
        a.sel <- a.ealy
        if (verbose > 0) cat("! Re-estimation w/auto-Easter/LY effects failed, keeping the test-based model.\n")
      }
    } else {
      if (verbose > 0) cat("! Even testing w/auto-Easter/LY effects failed! Keeping the original model, but this might indicate that the series are problematic.\n", sep = "")
    }
  } else {
    if (verbose > 1) cat("No tests for Easter or leap year were done because those parameters were overridden.\n")
  }

  keep.easter <- "Easter[8]" %in% names(a.sel$est$coefficients)
  keep.ly <- "Leap Year" %in% names(a.sel$est$coefficients)
  if (verbose > 0) cat("*** Easter effect: ", if (!keep.easter) "NO" else "YES", ", leap-year effect: ", if (!keep.ly) "NO" else "YES", ".\n", sep = "")
  if (verbose > 0) cat("*** Outlier regressors: ", paste0(sort(setdiff(a.sel$model$regression$variables, cal.names)), collapse = ", "), "\n", sep = "")
  if (verbose > 1) printSym("-", linelen)

  stat <- getSAStatOne(a.sel)

  # Preparing the separate decomposed series
  xs <- a.sel$series
  if (!force) {
    x.sa <- xs$d11
    x.adjfac <- xs$d16 # Seasonal and calendar effects
    x.seas <- xs$d10
    x.cal <- if (stat$log) x.adjfac / x.seas else x.adjfac - x.seas
  } else { # With forcing of annual totals, we redefine the SA component
    x.sa <- xs$saa # The adjusted is adjusted and annual-corrected
    x.adjfac0 <- xs$d16 # AF without annual correction
    x.seas0 <- xs$d10 # Seasonal without annual correction
    x.adjfac <- if (stat$log) x / x.sa else x - x.sa # Read adj. fac.
    x.cal  <- if (stat$log) x.adjfac0 / x.seas0 else x.adjfac0 - x.seas0 # Vanilla calendar
    x.seas <- if (stat$log) x.adjfac / x.cal else x.adjfac - x.cal
  }
  x.trend <- xs$d12
  x.ir <- xs$d13
  do.calend <- !isTRUE(all.equal(diff(as.numeric(x.cal)), rep(0, length(x.cal)-1)))

  rob.tests <- robustSeasTests(a.sel)
  robust.m7 <- rob.tests$robust.m7

  if (verbose > 1) {
    cat("Seasonality diagnostics for ", name, ":\n", sep = "")
    print(do.call(c, stat[paste0("M", 1:11)]))
    cat("Overall quality, Q_M2: ", stat$Q_M2, ".\n", sep = "")
  } else {
    if (verbose > 0) cat("Diagnostics: M7 = ", stat$M7, ", robust M7 = ", robust.m7, ", Q_M2 = ", stat$Q_M2, " (<1 is OK).\n", sep = "")
  }

  # In case there were missing values, put the regression-based prediction for the 'predicted' series
  x.pred <- x
  if (length(miss.inds) > 0) x.pred[miss.inds] <- predict.seas(a.sel)[miss.inds]
  # Creating the back- and forecast series if such were requested
  x.bc <- if (nbcast > 0) seasonal::series(a.sel, "bct") else NULL
  x.fc <- if (nfcast > 0) seasonal::series(a.sel, "fct") else NULL
  # Prepending backcasts, appending forecasts
  x.bf <- if (nbcast > 0) stats::ts(c(x.bc[, "backcast"], x.pred), end = stats::end(x.pred), freq = freq) else x.pred
  if (nfcast > 0) x.bf <- stats::ts(c(x.bf, x.fc[, "forecast"]), start = stats::start(x.bf), freq = freq)

  # Now the choice has to be made: adjust the seasonal component or not?
  # 4 options: robust M7, original M7, average, or enforce the adjustment (assume M7=0)
  # Calendar adjustment will always be performed if RegARIMA detected it
  # If it was not detected, then, the adjusement factor is equal to just the seasonal component
  seas.strength <- switch(sa.rule, robustM7 = robust.m7, M7 = stat$M7, averageM7 = (robust.m7+stat$M7)/2, yes = -1, no = 4)
  seas.strength <- round(seas.strength, 3)
  do.seasadj <- seas.strength < m7.threshold
  if (do.calend | do.seasadj) { # There is some seasonality
    if (verbose > 0) {
      if (verbose > 1) {
        if (stat$Q_M2 <= 1) {
          cat("There seem to be no immediate seasonality adjustability problems in ", name, ".\n", sep = "")
        } else {
          cat("The adjustment could be mediocre (or the trend is flat, which is harmless).\nLook at the plot just in case.\n", sep = "")
        }
      }
    }

    # If there are only calendar but no seasonal effects, then, the adj. factor is equal to the calendar component
    if (!do.seasadj) {
      if (force) {
        warning("Forcing annual totals, but there is only calendar (no seas.) effect, which, by definition, changes the annual total. Removing calendar adjustment; X.SA = X!")
        x.adjfac <- x; x.adjfac[1:length(x)] <- as.numeric(stat$log) # Copying the TS attributes to a constant 1 or 0
        x.sa <- x
        x.ir <- if (stat$log) x / x.trend else x - x.trend
      } else {
        x.adjfac <- x.cal
        x.sa <- if (stat$log) x / x.cal else x - x.cal
        x.ir <- if (stat$log) x / x.cal / x.trend else x - x.cal - x.trend
      }
      x.seas <- x; x.seas[1:length(x)] <- as.numeric(stat$log)
    }
    multiplicative <- stat$log
  } else { # There is neither calendar nor seasonal effect
    if (verbose > 1) {
      if (stat$Q_M2 <= 1) {
        cat("The series ", name, " do not have stable seasonality and seem regular, skipping.\n", sep = "")
      } else {
        cat("The series ", name, " do not have stable seasonality and are very noisy, skipping.\n", sep = "")
      }
    }
    x.sa <- x # No adjustment is done
    x.seas <- x; x.seas[1:length(x.seas)] <- as.numeric(stat$log) # Copying the TS attributes
    x.cal <- x.adjfac <- x.seas
    x.ir <- if (stat$log) x / x.trend else x - x.trend
    multiplicative <- FALSE
  }
  seasonality <- c(seasonal = do.seasadj, calendar = do.calend)

  # Adding the metadata to be extracted by getSAStat
  attr(a.sel, "seasonality") <- seasonality
  attr(seasonality, "rule") <- attr(a.sel, "rule") <- sa.rule
  attr(seasonality, "rule.value") <- attr(a.sel, "rule.value") <- seas.strength
  attr(seasonality, "threshold") <- attr(a.sel, "threshold") <- m7.threshold
  attr(a.sel, "transform") <- transform
  attr(a.sel, "td") <- td
  attr(a.sel, "ly") <- leap.year
  attr(a.sel, "easter") <- easter

  cat("*** Calendar eff.: ", if (do.calend) "YES (AICc)" else "NO (AICc)", ", seasonal eff.: ", if (do.seasadj) "YES" else "NO", " (",
      switch(sa.rule, robustM7 = "rob. M7 = ", M7 = "M7 = ", averageM7 = "avg. M7 = ", yes = "forced", no = "forbidden"),
      if (!(sa.rule %in% c("yes", "no"))) paste0(seas.strength, if (do.seasadj) " <= " else " > ", m7.threshold) else "", ").\n", sep = "")

  out.ts <- stats::ts.union(original = x, adjusted = x.sa, adjfac = x.adjfac, seasonal = x.seas,
                            calendar = x.cal, trend = x.trend, irregular = x.ir, predicted = x.bf)
  # Restoring the middle missing values in the adjusted series
  mg <- findMiddleGaps(out.ts[, "adjusted"])
  if (length(mg) > 0) {
    out.ts[mg, "adjusted"] <- if (stat$log) out.ts[mg, "predicted"] / out.ts[mg, "seasonal"] / out.ts[mg, "calendar"] else
      out.ts[mg, "predicted"] - out.ts[mg, "seasonal"] - out.ts[mg, "calendar"]
  }

  # Plotting if plot.file is not NA
  if (!isTRUE(is.na(plot.file))) {
    if (!is.null(plot.file)) grDevices::png(plot.file, 800, 900, type = "cairo", pointsize = 20)
    ext.xlim <- c(if (nbcast > 0) -nbcast/freq else 0, if (nfcast > 0) nfcast/freq else 0)
    plotSeas(x = a.sel, sa.custom = x.sa, name = name, extend.xlim = ext.xlim, skip.boxplot = skip.boxplot)
    if (!is.null(plot.file)) grDevices::dev.off()
  }

  # Finally, trimming the residual series if it is longer
  if (length(xs$rsd) != length(xs$d11)) {
    a.sel$series$rsd <- stats::window(a.sel$series$rsd, start = stats::start(a.sel$series$d11), end = stats::end(a.sel$series$d11), extend = TRUE)
  }
  attr(seasonality, "identify") <- list(rule = sa.rule, threshold = m7.threshold, actual = seas.strength)

  ret <- list(seasonality = seasonality, multiplicative = multiplicative, quality = stat,
              series = out.ts, seas = a.sel, date = ts2date(out.ts), si.robtests = rob.tests, call = match.call())

  # Last check: positivity
  if (all(x[is.finite(x)] > 0)) { # The original data are positive
    if (verbose > 1) printSym("-", linelen)
    if (verbose > 1) cat("The input series contain only positive values in the range", range(x, na.rm = TRUE), "\n")
    if (all(x.sa[is.finite(x.sa)] > 0)) {
      if (verbose > 1) cat("The adjusted series contain values in the range", range(x.sa, na.rm = TRUE), "\n")
    } else {
      if (verbose > 0) cat("! The adjusted series contain some negative values, range", range(x.sa, na.rm = TRUE), "\n! Enforcing a log transformation may be necessary.\n")
    }
    if (trans.fun != "log" & verbose > 1) cat("! All the values of ", name, " are positive.\n", sep = "")
  }
  if (verbose > 1) printSym("=", linelen)
  if (verbose > 0) cat("\n")

  attr(ret, "type") <- "diagnoseSeasonality" # To be easily parsed by getSAStat
  attr(ret$seas, "seriesname") <- name
  attr(ret$seas, "custom.calendar") <- (!is.null(calendar))
  return(ret)
}


#' Compute discrepancy / revision summary statistics
#'
#' @param old Numeric: the series before revision
#' @param new Numeric: the series after revision
#' @param nspans Numeric: number of spans. Returned in the output, does not affect the result.
#' @param overlap.length Numeric: length of one overlap of two spans. Returned in the output, does not affect the result.
#'
#' @return A numeric vector with 5 elements: mean revision, mean absolute revision, mean absolute percentage revision
#' (known as Mean Avg. % Error in forecasting), SD of revision, relative mean absolute revision (ratio between mean
#' absolute revision and mean absolute old series), ratio of standard deviations (SD(revision) / SD(old))
#' @export
#'
#' @examples
#' before <- c(1.2, 1.2, 0.5, 1.1, 1.7, 0.6, 0.9, 1, 1.7, 0.8)
#' after <- c(1.3, 1.2, 0.6, 1.3, 2.1, 0.2, 0.8, 1.1, 1.8, 0.8)
#' diagnoseRevisions(before, after)
diagnoseRevisions <- function(old, new, nspans = NA, overlap.length = NA) {
  d <- new - old
  mean.rev <- mean(d)
  mean.abs.rev <- mean(abs(d))
  mean.abs.pct.rev <- mean(abs(d / old))*100
  sd.rev <- stats::sd(d)
  rel.mean.abs.rev <- mean.abs.rev / mean(abs(old))
  sd.rev.to.sd.sa <- sd.rev / stats::sd(old)
  stability.inds <- c(mean.rev = mean.rev, mean.abs.rev = mean.abs.rev, mean.abs.pct.rev = mean.abs.pct.rev,
                      sd.rev = sd.rev, rel.mean.abs.rev = rel.mean.abs.rev, sd.rev.to.sd.sa = sd.rev.to.sd.sa,
                      nspans = nspans, overlap.nyears = overlap.length)
  return(stability.inds)
}

#' Diagnose seasonality, adjust the series, print and visualise the results
#'
#' @param x A univariate or multivariate time series of class `ts` or `mts`.
#' `diagnoseSeasonalityOne` accepts only univariate time series.
#' @param calendar A country name (character) or an object of class `mts`
#' containing calendar regressors: Monday, ..., Saturday, LeapYear, WorkingDays.
#' If NULL, the default X13 working-day calendar is used.
#' If forecasting is required, should be longer than the series.
#' @param name A string with the name of the series for printing.
#' @param bfcast.tails Logical: if x starts or ends with NA, produce the backcasts and forecasts for those periods?
#' @param est.begin An integer vector of length 2 denoting the beginning of the estimation sample. By default, the first valid observation.
#' @param est.end An integer vector of length 2 denoting the end of the estimation sample. By default, no forecasts are produced.
#' @param bcast.begin An integer vector of length 2 denoting the first period till which backcasts should be computed. By default, no backcasts are produced.
#' @param fcast.end An integer vector of length 2 denoting the last period till which backcasts should be computed. By default, no forecasts are produced.
#' @param transform Character or logical: automatic transformation, or forced logarithmic, or none.
#' TRUE translates to "log", FALSE translates to "none", NA or NULL to "auto".
#' @param td Character or numeric: trading-day dummies (6 trading-day, 1 weekday, none, or automatic)
#' @param leap.year Character: automatic, forced, or no leap-year effect
#' @param easter Character: automatic, forced, or no Easter effect
#' @param forced.outliers A character vector of ourlier dates and types in X13 format (e.g. `c("ao2010.1", "tc2020.3")`)
#' @param force.annual A character indicating whether a correction of the SA series to add up to the annual totals should be made.
#' @param split.long If TRUE, long time series will be adjusted in chunks and then, blended together.
#' Ignored if `custom.starts` is not `NULL`
#' @param max.length Numeric: length of a chunk in years. Ignored if `split.long` is `FALSE` or if
#' `custom.[est.starts, est.ends]` is not `NULL`.
#' @param overlap.length Numeric: the end of one chunk will be blended with the beginning of another chunk over this number of years
#' @param custom.est.starts A list of length-2 vectors, a data frame or a matrix with two columns:
#' year and period denoting the starts of each chunk
#' @param custom.est.ends A list of length-2 vectors, a data frame or a matrix with two columns:
#' year and period denoting the ends of each chunk
#' @param custom.blend.starts A list of length-2 vectors, a data frame or a matrix with two columns:
#' year and period denoting the starts of blending periods between 2 estimation windows
#' @param custom.blend.ends A list of length-2 vectors, a data frame or a matrix with two columns:
#' year and period denoting the ends of blending periods between 2 estimation windows
#' @param m7.threshold Numeric between 0 and 3: the procedure will produce SA series if the M7 statistic is less than this threshold. A value of 3 means that adjustment is always done.
#' @param sa.rule Character: what kind of indicator compare to `m7.threshold`.
#' @param transform.aicdiff Numeric: the AICc difference between the additive (linear) and multiplicative (log) model; use log if AICc(nolog) - AICc(log) < aicdiff; by default, prefers multiplicative models
#' @param tradingdays.aicdiff Numeric: the AICc difference for the TD or WD dummies to be added into the model (if negative, prefers more parsimonious models)
#' @param plot.file A string containing the path to the output PNG file. If NULL, plot to the current device. Using NA prevents any plot from being created.
#' @param skip.boxplot Logical: remove boxplots in the background of the seasonal plot?
#' @param verbose Logical or integer. TRUE or 1 = basic output, 2  = detailed output.
#' @param parallel If TRUE, estimates multiple series in parallel
#' @param mc.cores The number of cores for the cluster to speed up the computations.
#'
#' To list all available countries, run `data(package = "StatecCMP")`. Each country
#' should have a monthly and a quarterly calendar.
#'
#' If prediction out of the sample is requested and the calendar is custom, then,
#' the calendar should be as long as the forecast horizon.
#'
#' Since seasonal adjustment via X13 implies estimation with TRAMO-like methods that allow missingness,
#' the returned series "predicted" will contain fitted values in places where the original data had middle gaps.
#'
#' `diagnoseSeasonalityOne` does one adjustment, `diagnoseSeasonality` does the same
#' one adjustment if the series are univariate and short, and multiple adjustments
#' otherwise: for east series, if it is shorter than `max.length`, apply `diagnoseSeasonalityOne`,
#' or split it into chunks, adjust the chunk, and cleverly glue them back without
#' losing the accessible quality statistics for the last sample and adding the extra
#' measurements of how different the revised adjusted values are.
#'
#' For plotting and parallel processing, it is recommended to adjust all the series
#' first, and then, produce the plots based on the SA objects.
#'
#' @return A list:
#' * `seasonality`: a logical indicating whether there is either seasonal or calendar effects
#' * `multiplicative`: a logical indicating whether the adjustment model is multiplicative
#' * `quality`: a list containing brief diagnostic information
#' * `series`, an mts containing the decomposition
#' * `seas`: the full output of the X13 `seas()` run for the best model)
#' * `date`: the dates of the input series
#' * `si.robtests`: robust seasonality test results (IRLS with HAC for seasonal and/or yearly dummies)
#' @export
#'
#' @examples
#' # Simulating ARIMA with trend, decaying seasonality, 1 LS and 1 AO
#' addTS <- function(x) ts(x, end = c(2022, 12), freq = 12)
#' x <- seq(0, 7.99, 1/12)
#' set.seed(1)
#' ybase <- addTS(arima.sim(n = 96, list(ar = 0.8, ma = 0.3))) # ARIMA(1, 0, 1)
#' yseas <- addTS(40/(x+20)*cos(x*pi*2))
#' oreg <- addTS(5*(x >= 6) - 4*(x == 5))
#' y1 <- 10 + 0.05*x + oreg + yseas + ybase
#' ysa <- diagnoseSeasonality(y1)
#' summary(ysa$seas)
#' ysa.mult <- diagnoseSeasonality(y1, transform = "yes")
#' # ysa.add <- diagnoseSeasonality(y1, transform = "no")
#' ysa.td6 <- diagnoseSeasonality(y1, td = 6)
#' # ysa.td1 <- diagnoseSeasonality(y1, td = 1)
#' # ysa.td0 <- diagnoseSeasonality(y1, td = 0)
#' ysa.ly <- diagnoseSeasonality(y1, leap.year = "yes")
#' # ysa.ly0 <- diagnoseSeasonality(y1, leap.year = "no")
#' ysa.ea <- diagnoseSeasonality(y1, easter = "yes")
#' ysa.le <- diagnoseSeasonality(y1, easter = "yes", leap.year = "yes")
#' # A completely manual model = several times faster
#' ysa.man <- diagnoseSeasonality(y1, td = 1, transform = "no", easter = "yes", leap.year = "no")
#' # Trim the sample, forecast tails
#' yfc <- diagnoseSeasonality(y1, transform = "no", td = "no",
#'          est.begin = c(2016, 7), est.end = c(2022, 9), bfcast.tails = TRUE)
#' plot(yfc$series[ , c("original", "predicted")], plot.type = "single", lty = 1:2)
#' # Similarly, if there are missing values at the tails, they can be easily forecast
#' # Simply add `bfcast.tails = TRUE`
#' y1m <- y1; y1m[89:96] <- NA; y1m[1:12] <- NA
#' yfc2 <- diagnoseSeasonality(y1m, bfcast.tails = TRUE, transform = "no", td = "no")
#' plot(yfc2$series[ , c("original", "predicted")], plot.type = "single", lty = 1:2)
#' # Forcing outliers
#' diagnoseSeasonality(y1, forced.outliers = "ao2016.1", transform = "no", td = "no")
#' # Missing values in the middle
#' y1g <- y1; y1g[39:40] <- NA
#' ymisssa <- diagnoseSeasonality(y1g, transform = "no", td = "no")
#' window(ts.union(ymisssa$series[, c("original", "predicted")]), start = c(2013, 1), end = c(2013, 6))
#'
#' # Only seasonality
#' y2 <- ts(10 + 1.5*cos(x*pi*2) + ybase, start = c(2010, 1), freq = 12)
#' y2sa <- diagnoseSeasonality(y2)
#'
#' # Only calendar effects---suppose that Fridays bring more sales
#' dates <- seq.Date(as.Date("2015-01-01"), as.Date("2022-12-31"), by = "day")
#' wdays   <- as.numeric(strftime(dates, "%u"))
#' yearmonth <- format(dates, "%Y-%m")
#' nfridays <- aggregate(wdays == 5, by = list(yearmonth), FUN = sum)$x
#' y3 <- addTS(10 + 0.2*x + ybase + nfridays)
#' plot(y3) # No obvious seasonality
#' y3sa <- diagnoseSeasonality(y3)
#' y3lu <- diagnoseSeasonality(y3, calendar = "Luxembourg")
#'
#' # Only weekend effects
#' set.seed(1)
#' nwe <- aggregate(wdays %in% 6:7, by = list(yearmonth), FUN = sum)$x / 5
#' y4 <- addTS(rnorm(96) + nwe)
#' y4sa <- diagnoseSeasonality(y4)
#' plot(ts.union(y4sa$series[, "adjfac"], nwe - mean(nwe)),
#'      plot.type = "single", col = 2:1)
#'
#' # Quarterly data
#' apsa <- diagnoseSeasonality(aggregate(datasets::AirPassengers, nfreq = 4))
#' summary(apsa$seas)
#'
#' # Long series without any seasonal effects
#' set.seed(1)
#' ywn <- addTS(rnorm(500))
#' ylong <- diagnoseSeasonality(ywn, verbose = 1)
#' plot(ylong$series, main = "Piecewise SA")
#' plotSeas(ylong$seas) # What X13 sees
#' plotSeas(ylong$seas, sa.custom = ylong$series[, "adjusted"]) # What is used
#' ylong2 <- diagnoseSeasonality(ywn, verbose = 1, return.only.last = FALSE)
#' plotSeas(ylong2[[1]]$seas, sa.custom = ylong2[[1]]$series[, "adjusted"])
#' plotSeas(ylong2[[3]]$seas)
#' # Stability diagnostics
#' ylong$stability.inds # Erratic because the original is white noise
#'
#' # Long series with strong but devolving seasonality
#' x <- seq(0, 40, 1/12)
#' set.seed(1)
#' ybase <- addTS(arima.sim(n = 481, list(ar = 0.8, ma = 0.3))) # ARIMA(1, 0, 1)
#' yseas <- 60/(x+10)*cos(x*pi*2)
#' ylong3 <- addTS(10 + 0.07*x + yseas + ybase)
#' plot(ts.union(ylong3, ylong3 - ybase), col = 1:2, plot.type = "single")
#' # We try windows of length 10, 15, 20
#' y10 <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'                            max.length = 10)
#' y15 <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'                            max.length = 15)
#' y15n <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'                            max.length = 15, overlap.length = 0)
#' y25 <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'                            max.length = 25, overlap.length = 1/12) # Only 1-point overlap
#' getSAStat(list(y10, y15, y15n, y25))
#' plotSeas(y10)
#' plotSeas(y15)
#' plotSeas(y15n)
#' plotSeas(y25)
#' plotSeas(y25, span.index = 1)
#' round(cbind(`10+1y` = y10$stability.inds, `15+1y` = y15$stability.inds,
#'             `15+0` = y15n$stability.inds, `25+1m` = y25$stability.inds), 3)
#' # As expected, with no overlap, there are no diagnostics, and with 1-point overlap,
#' # the diagnostics are unreliable, and no SDs are available
#'
#' # Long series with forecasting and backcasting
#' ylong.fc <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'               bfcast.tails = TRUE, bcast.begin = c(1981, 1),
#'               fcast.end = c(2000, 4), est.end = c(1999, 4), max.length = 10)
#' plotSeas(ylong.fc)
#' yshort.fc <- diagnoseSeasonality(ylong3, transform = "no", easter = "no", leap.year = "no",
#'               bfcast.tails = TRUE, est.end = c(1995, 4), bcast.begin = c(1981, 1),
#'               fcast.end = c(1995, 12))
#' plotSeas(yshort.fc)
#' # No problems when there is only one point
#' yshort.fc2 <- diagnoseSeasonality(ylong3, est.end = c(1999, 12), max.length = 10,
#'                bcast.begin = c(1982, 11), fcast.end = c(2000, 1), bfcast.tails = TRUE,
#'                transform = "no", easter = "no", leap.year = "no")
#' plotSeas(yshort.fc2)
#'
#' set.seed(1)
#' x <- ts(rnorm(500), end = c(2022, 12), freq = 12)
#' a <- diagnoseSeasonality(x, max.length = 15)
#' # Custom breaks
#' ce.st <- matrix(c(1981, 5, 1990, 1, 2005, 8, 2015, 2), ncol = 2, byrow = TRUE)
#' ce.en <- matrix(c(1991, 1, 2007, 1, 2015, 1, 2022, 12), ncol = 2, byrow = TRUE)
#' a2 <- diagnoseSeasonality(x, td = 0, transform = "no", leap.year = "no", easter = "no",
#'                           custom.est.starts = ce.st, custom.est.ends = ce.en)
#' # One break = supply the breaks as a list
#' a3 <- diagnoseSeasonality(AirPassengers, transform = "yes", td = 1,
#'                           custom.est.starts = list(c(1949, 1),  c(1953, 1)),
#'                           custom.est.ends   = list(c(1957, 12), c(1960, 12)),
#'                           custom.blend.starts = list(c(1953, 7)),
#'                           custom.blend.ends   = list(c(1953, 12))
#'                           )
#'
#' # Ill-behaved series that X13 cannot digest (failure guaranteed)
#' set.seed(1)
#' # 1. Just a constant
#' bad.series1 <- rep(100, 100)
#' # 2. Piecewise linear series (e.g. VAT data look like this)
#' bad.series2 <- c(rep(100, 90), 50, rep(125, 9))
#' # 3. A series was measured in millions at first, in original units later
#' bad.series3 <- rnorm(100) * c(rep(1, 50), rep(10^6, 50))
#' bad.series <- ts(cbind(bad.series1, bad.series2, bad.series3),
#'                  start = c(2001, 1), frequency = 12)
#' plot(bad.series)
#' a.bad <- diagnoseSeasonality(bad.series) # Series 3 takes 5 minutes
diagnoseSeasonality <- function(x, calendar = NULL, name = NULL,
                                bfcast.tails = FALSE,
                                est.begin = NULL, est.end = NULL,
                                bcast.begin = NULL, fcast.end = NULL,
                                transform = c("auto", "yes", "no"),
                                td = c("auto", "6", "1", "no"),
                                leap.year = c("auto", "yes", "no"),
                                easter = c("auto", "yes", "no"),
                                forced.outliers = NULL,
                                force.annual = c("no", "regress", "denton"),
                                m7.threshold = 1,
                                sa.rule = c("robustM7", "M7", "averageM7", "yes", "no"),
                                transform.aicdiff = -2, tradingdays.aicdiff = 0,
                                plot.file = NULL, skip.boxplot = FALSE,
                                verbose = 2,
                                split.long = TRUE, max.length = 15, overlap.length = 1,
                                custom.est.starts = NULL, custom.est.ends = NULL,
                                custom.blend.starts = NULL, custom.blend.ends = NULL,
                                parallel = FALSE, mc.cores = 4
) {
  # A wrapper for diagnoseSeasonalityOne with all the input parameters passed down except for the sample-defining and plotting ones
  diagSeas <- function(x, name = NULL, nbcast = 0, nfcast = 0, plot.file = NULL, skip.boxplot = FALSE) {
    if (is.null(name)) {
      cl <- match.call()
      name <- deparse(as.list(cl)$x)
    }
    if (is.null(name)) name <- "Series"
    diagnoseSeasonalityOne(x, calendar = calendar, name = name,
                         transform = transform, td = td, leap.year = leap.year, easter = easter,
                         forced.outliers = forced.outliers, force.annual = force.annual,
                         m7.threshold = m7.threshold, sa.rule = sa.rule,
                         transform.aicdiff = transform.aicdiff, tradingdays.aicdiff = tradingdays.aicdiff,
                         plot.file = plot.file, skip.boxplot = skip.boxplot, verbose = verbose,
                         nfcast = nfcast, nbcast = nbcast)
  }

  # Diagnose one time series, splitting it into chunks
  # Backcasting is allowed in the first chunk, forecasting in the last one
  diagSeasChunks <- function(x, name = "Series",
                             est.begin = NULL, est.end = NULL,
                             custom.est.starts = NULL, custom.est.ends = NULL,
                             custom.blend.starts = NULL, custom.blend.ends = NULL,
                             bcast.begin = NULL, fcast.end = NULL,
                             bfcast.tails = FALSE, plot.file = NULL,
                             skip.boxplot = FALSE) {
    full.range <- getRangeVec(x)
    n <- length(x)
    f <- stats::frequency(x)
    x.inds <- time2ind(x, freq = f)
    valid.begin <- x.inds[full.range["first"], ]
    valid.end <- x.inds[full.range["last"], ]

    # If there are custom splits, adjust the estimation sample
    custom.splits <- (!is.null(custom.est.starts)) & (!is.null(custom.est.ends))
    custom.blends <- (!is.null(custom.blend.starts)) & (!is.null(custom.blend.ends))
    if (custom.splits) { # Converting the splits into a list of indices
      vs <- c("custom.est.starts", "custom.est.ends", if (custom.blends) c("custom.blend.starts", "custom.blend.ends") else NULL)
      for (v in vs) {
        gv <- get(v)
        if (is.data.frame(gv)) assign(v, as.matrix(gv)) # DF to matrices all at once
        gv <- get(v)
        if (is.matrix(gv)) assign(v, split(gv, 1:nrow(gv))) # Matrices to lists all at once
      }
    }

    # Computing the maximum possible range of observations
    # Estimation start: if it is less than the first valid observation, warn
    if (is.null(est.begin)) est.begin <- if (!custom.splits) valid.begin else custom.est.starts[[1]]
    if (ind2time(est.begin, f) < ind2time(valid.begin, f)) {
      est.begin <- valid.begin
      warning(paste0("Estimation for ", name, " cannot start earlier than the first valid observation, ", valid.begin[1], "-", valid.begin[2], "."))
    }
    if (is.null(est.end)) est.end <- if (!custom.splits) valid.end else custom.est.ends[[length(custom.est.ends)]]
    if (ind2time(est.end, f) > ind2time(valid.end, f)) {
      est.end <- valid.end
      warning(paste0("Estimation for ", name, " cannot end earlier than the last valid observation, ", valid.end[1], "-", valid.end[2], "."))
    }
    if (bfcast.tails) { # If forecasts and backcasts need to be produced
      if (is.null(bcast.begin)) {
        # If there are head NAs, backcast them automatically
        bcast.begin <- if (full.range["first"] == 1) est.begin else x.inds[1, ]
      }
      if (ind2time(bcast.begin, f) > ind2time(est.begin, f)) {
        bcast.begin <- est.begin
        warning(paste0("Backcasting ('bcast.begin') for ", name, " must start earlier than the estimation sample start, ", .p0(est.begin, "-"), "; ignoring it."))
      }
      if (is.null(fcast.end)) {
        # If there are tail NAs, forecast them automatically
        fcast.end <- if (full.range["last"] == n) est.end else x.inds[n, ]
      }
      if (ind2time(fcast.end, f) < ind2time(est.end, f)) {
        fcast.end <- est.end
        warning(paste0("Forecasting ('fcast.end') for ", name, " must end later than the estimation sample end, ", .p0(est.end, "-"), "; ignoring it."))
      }
    } else {
      if (!is.null(bcast.begin) | !is.null(fcast.end)) {
        warning("diagnoseSeasonality: no forecasts were requested (bfcast.tails = FALSE), but 'bcast.begin' or 'fcast.end' are non-empty; ignoring them.")
      }
      bcast.begin <- est.begin
      fcast.end <- est.end
    }

    x.trim <- stats::window(x, start = est.begin, end = est.end) # The series for analysis without missingness at the tails
    nt <- length(x.trim)
    linelen <- 60

    head.years <- unname(ind2time(est.begin, f) - ind2time(bcast.begin, f))
    tail.years <- unname(ind2time(fcast.end, f) - ind2time(est.end, f))
    if (head.years > 10.001) stop("Cannot backcast more than 10 years (X13 limitation). Reduce the number of backcasts.")
    if (tail.years > 10.001) stop("Cannot forecast more than 10 years (X13 limitation). Reduce the number of forecasts.")
    nbcast <- round(head.years * f) # How many points are to be backcast at the beginning of the series
    nfcast <- round(tail.years * f) # How many points are to be forecast at the end of the series

    # Split the series into spans if it is too long or if custom splits were requested
    nmax <- round(max.length * f)
    nlap <- round(overlap.length * f)

    if (nt > nmax | custom.splits) {
      if (custom.splits) { # Case 1: the user requested their own splits
        # Finding the indices of the corresponding dates in the series indices
        findRow <- function(x) which(apply(x.inds, 1, function(y) all(x == y)))
        ce.st <- sapply(custom.est.starts, findRow)
        ce.en <- sapply(custom.est.ends, findRow)
        cb.st <- if (is.null(custom.blend.starts)) NULL else sapply(custom.blend.starts, findRow)
        cb.en <- if (is.null(custom.blend.ends)) NULL else sapply(custom.blend.ends, findRow)
        chunks <- splitOverlapCustom(n = nt, est.starts = ce.st, est.ends = ce.en, blend.starts = cb.st, blend.ends = cb.en)
      } else { # Case 2: the series are long -- cut them automatically
        chunks <- splitOverlapFixed(n = nt, nmax = nmax, l = nlap, plot = FALSE)
      }
      rownames(chunks$chunks) <- rownames(chunks$weights) <- as.character(ts2date(x.trim))

      p <- ncol(chunks$chunks)
      ranges <- apply(chunks$chunks, 2, range, na.rm = TRUE)
      ix <- time2ind(x.trim)
      x.sub <- lapply(1:p, function(i) stats::window(x.trim, start = ix[ranges[1, i], ], end = ix[ranges[2, i], ]))

      # Backcast only in the first chunk, forecast only the in last chunk, and do not plot yet
      par.list <- vector("list", p)
      for (i in 1:p) par.list[[i]] <- list(x = x.sub[[i]], name = paste0(name, ", chunk ", i, "/", p),
                                           nbcast = 0, nfcast = 0, plot.file = NA)
      par.list[[1]]$nbcast <- nbcast
      par.list[[p]]$nfcast <- nfcast

      if (verbose > 1) {
        printSym("#", linelen)
        cat("Diagnosing seasonality in ", name, " via with parameter auto-selection.\n", sep = "")
        if (!custom.splits) cat("The series of ", round(nt/f, 1), " years are split every ", max.length, " years with ", round(overlap.length*f), "-point overlap.\n", sep = "") else
          cat("The series of ", round(nt/f, 1), " years are split into user-defined chunks.\n", sep = "")
      }

      res.list <- lapply(1:p, function(i) do.call(diagSeas, par.list[[i]]))
      # Restoring the names without chunks
      for (i in 1:p) attr(res.list[[i]]$seas, "seriesname") <- name

      nzw <- chunks$weights > 0
      nzranges <- apply(nzw, 2, function(x) range(which(x)))
      # Overlap statistics if the overlap length is greater than 0
      # With custom chunks, overlaps <=> 0 < w < 1
      nlap.total <- if (custom.splits) sum(chunks$weights > 0 & chunks$weights < 1) else nlap
      if (custom.splits) overlap.length <- nlap.total / p / f # Recomputed to make sense: average overlap length
      if (nlap.total > 0) {
        # The old weights are always decreasing; the new ones are increasing.
        # The mixing cannot start from the 1st observation by construction (check in splitOverlapCustom) --> first FALSE
        # The weight index is relative to the first observation of the estimation sample
        getWI <- function(i, old = TRUE) {
          cond1 <- chunks$weights[, i] > 0 & chunks$weights[, i] < 1
          cond2 <- c(FALSE, if (old) diff(chunks$weights[, i]) < 0 else diff(chunks$weights[, i]) > 0)
          inds <- which(cond1 & cond2)
          first.obs <- which(!is.na(chunks$chunks[, i]))[1]
          inds <- inds - first.obs + 1
          return(inds)
        }
        old.inds <- lapply(1:p, getWI, old = TRUE)[-p] # The last sample should never have old indices (it is not blended with any new series)
        new.inds <- lapply(1:p, getWI, old = FALSE)[-1] # The first sample is never blended at start
        old.sa <- lapply(1:(p-1), function(i) if (length(old.inds[[i]]) > 0) res.list[[i]]$series[old.inds[[i]], "adjusted"] else NULL)
        new.sa <- lapply(1:(p-1), function(i) if (length(new.inds[[i]]) > 0) res.list[[i+1]]$series[new.inds[[i]], "adjusted"] else NULL)
        old.sa <- as.numeric(stats::na.omit(unlist(old.sa)))
        new.sa <- as.numeric(stats::na.omit(unlist(new.sa)))
      } else old.sa <- new.sa <- NA
      stability.inds <- diagnoseRevisions(old.sa, new.sa, nspans = p, overlap.length = overlap.length)
      attr(res.list[[p]]$seas, "stability.inds") <- stability.inds # Saving for extraction by getSAStat

      # Blending the series and the 'seas' internal series, too
      series.list <- lapply(res.list, "[[", "series")
      # Trimming off the head of the first and the tail of the last series
      # to glue them back later (ro respect the dimensionality)
      x13.list <- lapply(res.list, function(x) do.call(stats::ts.union, x$seas$series))
      x13data.list <- lapply(res.list, "[[", c("seas", "data"))

      completeSeries <- function(x) {
        cn <- c("d10", "d11", "d12", "d13", "d16", "e18", "bct.backcast", "bct.lowerci", "bct.upperci", "fct.forecast", "fct.lowerci", "fct.upperci")
        m <- matrix(NA, nrow = nrow(x), ncol = length(cn))
        colnames(m) <- cn
        cm <- intersect(colnames(x), cn)
        m[, cm] <- as.matrix(x[, cm])
        m <- stats::ts(m, start = stats::start(x), frequency = stats::frequency(x))
      }

      if (nbcast > 0) {
        series.head  <- stats::window(series.list[[1]],  end = time2ind(ind2time(ix[1, ], f) - 1/f, f))
        x13.head     <- completeSeries(stats::window(x13.list[[1]],     end = time2ind(ind2time(ix[1, ], f) - 1/f, f)))
        x13data.head <- stats::window(x13data.list[[1]], end = time2ind(ind2time(ix[1, ], f) - 1/f, f))
        series.list[[1]]  <- stats::window(series.list[[1]],  start = ix[1, ])
        x13.list[[1]]     <- completeSeries(stats::window(x13.list[[1]],  start = ix[1, ])) # Has the extra bct component
        x13data.list[[1]] <- stats::window(x13data.list[[1]], start = ix[1, ])
        # Adding empty backcast columns to other elements of X13
        for (i in 2:p) x13.list[[i]] <- completeSeries(x13.list[[i]])
      }

      if (nfcast > 0) {
        series.tail  <- stats::window(series.list[[p]],  start = time2ind(ind2time(ix[nt, ], f) + 1/f, f))
        x13.tail     <- completeSeries(stats::window(x13.list[[p]],     start = time2ind(ind2time(ix[nt, ], f) + 1/f, f)))
        x13data.tail <- stats::window(x13data.list[[p]], start = time2ind(ind2time(ix[nt, ], f) + 1/f, f))
        series.list[[p]]  <- stats::window(series.list[[p]],  end = ix[nt, ])
        x13.list[[p]]     <- completeSeries(stats::window(x13.list[[p]],     end = ix[nt, ]))
        x13data.list[[p]] <- stats::window(x13data.list[[p]], end = ix[nt, ])
        # Adding empty forecast columns to other elements of X13
        for (i in 1:(p-1)) x13.list[[i]] <- completeSeries(x13.list[[i]])
      }

      blender <- function(l, weights, to.mts = TRUE) { # Blend a list of lists or 'mts's assuming that the dimensions match
        p <- length(l)
        if (p != ncol(weights)) stop("The weights for blending must have the same number of columns as the length of input list l.")
        # Works with a list of lists, returns a blended list; if the input is a list of 'mts's, break each into a list of 'ts's
        if (stats::is.ts(l[[1]])) {
          nms <- colnames(l[[1]])
          for (i in 1:length(l)) {
            l[[i]] <- lapply(1:ncol(l[[i]]), function(j) l[[i]][, j])
            names(l[[i]]) <- nms
          }
        }
        f <- stats::frequency(l[[1]][[1]])

        nms <- names(l[[1]])

        blended.list <- lapply(1:length(l[[1]]), function(i) { # For every sub-item of all list
          ser.list <- lapply(l, "[[", i)
          ser.ts <- do.call(stats::ts.union, ser.list); colnames(ser.ts) <- paste0("X", 1:p)
          # Putting zeros where the weights are zeros (which, by construction, are NA outside the chungs)
          # to preserve the inner NAs if the series have such
          ser.ts[is.na(ser.ts) & weights == 0] <- 0
          ser.ts <- ser.ts * weights # 0*0 = 0, 0*NA = NA
          ser.s <- rowSums(ser.ts)
          ser.s <- stats::ts(ser.s, start = stats::start(ser.ts), frequency = f)
        })
        names(blended.list) <- nms
        if (to.mts) {
          blended.ser <- do.call(stats::ts.union, blended.list)
          colnames(blended.ser) <- nms
          return(blended.ser)
        } else return(blended.list)

      }
      blended.series <- blender(series.list, chunks$weights, to.mts = TRUE)
      blended.x13 <- blender(x13.list, chunks$weights, to.mts = TRUE)
      blended.x13data <- blender(x13data.list, chunks$weights, to.mts = TRUE)
      if (nbcast > 0) {
        blended.series  <- stats::ts(rbind(series.head, blended.series), end = stats::end(blended.series), frequency = f)
        blended.x13     <- stats::ts(rbind(x13.head, blended.x13), end = stats::end(blended.x13), frequency = f)
        blended.x13data <- stats::ts(rbind(x13data.head, blended.x13data), end = stats::end(blended.x13data), frequency = f)
      }
      if (nfcast > 0) {
        blended.series  <- stats::ts(rbind(blended.series, series.tail), start = stats::start(blended.series), frequency = f)
        blended.x13     <- stats::ts(rbind(blended.x13, x13.tail), start = stats::start(blended.x13), frequency = f)
        blended.x13data <- stats::ts(rbind(blended.x13data, x13data.tail), start = stats::start(blended.x13data), frequency = f)
      }

      splits <- as.vector(nzranges)
      splits <- splits[-c(1, length(splits))] # The beginning and the end are trivial
      # Without overlap, only keep the newer splits for plotting (the shift is -0.5/freq)
      if (nlap == 0) splits <- splits[seq(2, length(splits), 2)]
      splits <- x.inds[splits, ]

      out <- res.list[[p]] # Taking the last one as the base
      out$series <- blended.series
      out$spans <- res.list
      out$combined.stats <- getSAStat(res.list)
      out$stability.inds <- stability.inds

      # Adding a spoofed "seas"-like for plotting for plotSeas
      # Combining all the outliers for spoofing
      combX13 <- res.list[[p]]$seas
      combX13$splits <- splits
      combX13$series <- as.list(blended.x13)
      if ("bct.backcast" %in% names(combX13$series)) {
        cn0 <- c("backcast", "lowerci", "upperci")
        cn <- paste0("bct.", cn0)
        bct <- do.call(stats::ts.union, combX13$series[cn])
        colnames(bct) <- cn0
        combX13$series[cn] <- NULL
        combX13$series$bct <- bct
      }
      if ("fct.forecast" %in% names(combX13$series)) {
        cn0 <- c("forecast", "lowerci", "upperci")
        cn <- paste0("fct.", cn0)
        fct <- do.call(stats::ts.union, combX13$series[cn])
        colnames(fct) <- cn0
        combX13$series[cn] <- NULL
        combX13$series$fct <- fct
      }
      combX13$data <- blended.x13data
      getOL <- function(x) { # Taken from seas::outliers()
        r <- x$model$regression$variables[grepl("\\.", x$model$regression$variables)]
        r <- r[!grepl("\\/", r)]
        r
      }
      combOL <- unlist(lapply(lapply(res.list, "[[", "seas"), getOL))
      # The loop-wise nature of seasonal::outlier allows us to not worry about identical outliers
      # The earlier ones will be overwritten by the latter ones if the dates are identical
      combX13$model$regression$variables <- combOL
      combX13$x <- blended.series[, "original"]
      out$combinedX13 <- combX13 # The presence of this element is checked in plotSeas and getSAStat
      out$chunks <- chunks

      # Plotting if plot.file is not NA
      if (!isTRUE(is.na(plot.file))) {
        if (!is.null(plot.file)) grDevices::png(plot.file, 800, 900, type = "cairo", pointsize = 20)
        ext.xlim <- c(if (nbcast > 0) -nbcast/f else 0, if (nfcast > 0) nfcast/f else 0)
        plotSeas(x = out, extend.xlim = ext.xlim, skip.boxplot = skip.boxplot)
        if (!is.null(plot.file)) grDevices::dev.off()
      }

      if (verbose > 1) {
        cat("Blended the ", p, " results into one, saved the original models as the [[\"spans\"]] list.\n", sep = "")
        printSym("#", linelen)
        cat("\n")
      }
    } else {
      out <- diagSeas(x = x.trim, name = name, nbcast = nbcast, nfcast = nfcast, plot.file = plot.file)
    }
    return(out)
  }

  # If there are multiple time series, process them as a list
  if ("mts" %in% class(x)) {
    k <- ncol(x)
    x.list <- lapply(1:k, function(i) x[, i])
    x.names <- colnames(x)
    if (is.null(x.names)) x.names <- paste0("Column", 1:ncol(x))

    # A function for one series
    dSCi <- function(i) diagSeasChunks(x.list[[i]], name = x.names[i], # The names are more important than the call, and they will be used in getSAStat
                                       est.begin = est.begin, est.end = est.end,
                                       custom.est.starts = custom.est.starts, custom.est.ends = custom.est.ends,
                                       custom.blend.starts = custom.blend.starts, custom.blend.ends = custom.blend.ends,
                                       bcast.begin = bcast.begin, fcast.end = fcast.end,
                                       bfcast.tails = bfcast.tails, plot.file = plot.file)

    if (parallel & mc.cores > 1) {
      cl <- if (.Platform$OS.type != "windows") parallel::makeForkCluster(mc.cores) else parallel::makeCluster(mc.cores)
      ret <- parallel::parLapplyLB(cl = cl, X = 1:k, dSCi)
      parallel::stopCluster(cl)
    } else {
      ret <- lapply(1:k, dSCi)
    }

  } else { # x is 1-dimensional
    cl <- match.call()
    if (is.null(name)) name <- deparse(as.list(cl)$x)
    if (is.null(name)) name <- "Series" # Fallback
    ret <- diagSeasChunks(x, name = name, est.begin = est.begin, est.end = est.end,
                          custom.est.starts = custom.est.starts, custom.est.ends = custom.est.ends,
                          custom.blend.starts = custom.blend.starts, custom.blend.ends = custom.blend.ends,
                          bcast.begin = bcast.begin, fcast.end = fcast.end,
                          bfcast.tails = bfcast.tails, plot.file = plot.file)
    ret$seas$call <- cl # Replacing the trivial internal one with the call of this function
  }

  return(ret)

}

#' Compute indices and weights for sub-series with overlapping
#'
#' @param n Integer: length of the full series
#' @param nmax Integer: chunk size
#' @param l Non-negative integer: overlap length; can be 0 for no overlap
#' @param est.starts Integer vector of indices denoting the starts of each chunk; should start with 1
#' @param est.ends Integer vector of indices denoting the ends of each chunk; should end with `n`
#' Must have the same length as `est.starts`.
#' @param blend.starts Integer vector of indices denoting the starts of blending periods between 2 estimation windows.
#' Must be of length 1 shorter than est.starts (because `k` chunks = `k-1` cuts between them), must start at 2 or later.
#' Must not be less that the corresponding element of `est.starts` (excluding the 1st
#' element): one cannot blend something with nothing; blending must occur within the overlap between the estimations.
#' @param blend.ends Integer vector of indices denoting the ends of blending periods between 2 estimation windows.
#' Must be of length 1 shorter than est.ends, must end at `n-1` or earlier.
#' Must not be greater that the corresponding element of `est.ends` (excluding the last element).
#' @param plot If `TRUE`, visualises the weight matrix and the sub-samples
#'
#' When a time series is estimated in chunks, the final results should be
#' merged with ts.union. Then, they can be simply multipled by the weights
#' and added by row.
#'
#' @return A list of two matrices: sub-sample indices by columns (with `NA`s for
#' unused indices) and weights for blending.
#' @export
#'
#' @examples
#' splitOverlapFixed(n = 43, nmax = 20, l = 4, plot = TRUE)
#' # Custom: estimate on 5 and 7 years, blend in between
#' splitOverlapCustom(n = 96, est.starts = c(1, 13), est.ends = c(60, 96), plot = TRUE)
#' # Custom: estimate on years 1-5 and 2-7 years, but blend during the year 3
#' splitOverlapCustom(n = 96, est.starts = c(1, 13), est.ends = c(60, 96),
#'                    blend.starts = 25, blend.ends = 36, plot = TRUE)
#' splitOverlapCustom(n = 150, est.starts = c(NA, 45, 90), est.ends = c(50, 100, 150), plot = TRUE)
#' splitOverlapCustom(n = 200, est.starts = c(1, 15, 90, 140), est.ends = c(50, 100, 139, NA),
#'                   blend.starts = c(15, 90, NA), blend.ends = c(27, 100, NA), plot = TRUE)
#' # The ends can be trimmed
#' splitOverlapCustom(n = 200, est.starts = c(8, 15, 90, 140), est.ends = c(50, 100, 139, 190),
#'                   blend.starts = c(15, 90, NA), blend.ends = c(27, 100, NA), plot = TRUE)
#' # If there is anything wrong with the indices, it will throw a meaningflu error
#' \dontrun{
#' # Does not work: gap in the middle
#' splitOverlapCustom(n = 100, est.starts = c(1, 51), est.ends = c(49, 100), plot = TRUE)
#' # Works: no gap in the middle
#' splitOverlapCustom(n = 100, est.starts = c(1, 51), est.ends = c(50, 100), plot = TRUE)
#' # Does not work because the blending window must be within the estimation sample overlap
#' splitOverlapCustom(n = 150, est.starts = c(NA, 45, 90), est.ends = c(50, 100, NA),
#'                    blend.starts = c(44, 89), blend.ends = c(51, 101), plot = TRUE)
#' }
splitOverlapFixed <- function(n, nmax, l, plot = FALSE) {
  if (n <= nmax) stop("The chunk length must be less than the sample size.")
  if (nmax < 2*l) stop("The chunk length must greater than 2*overlap.length.")
  n1 <- nmax - l # Non-overlapping points in a chunk
  p <- 1 + ceiling((n-nmax) / n1) # Last + how many non-overlapping ones
  chunks <- matrix(rep(1:n, p), ncol = p)
  mults <- matrix(0, nrow = n, ncol = p)
  # Chunks = indices of subsamples, mults = weights for linear blending
  chunks[1:(n-nmax), p] <- NA
  blend.w <- if (l > 0) (1:l)/(l+1) else NULL
  full.w <- c(blend.w, rep(1, nmax - 2*l), rev(blend.w)) # Mixing weights
  wif <- which(is.finite(chunks[, p]))
  mults[wif, p] <- 1
  if (l > 0) mults[wif[1:l], p] <- blend.w # Non-zero mixing weights for the overlap period
  if (p > 2) {
    for (i in (p-1):2) {
      wif <- stats::na.omit(chunks[, i+1]) - n1
      chunks[setdiff(1:n, wif), i] <- NA
      mults[is.finite(chunks[, i]), i] <- full.w
    }
  }
  # Final chunk: from the first observation; take only the part not covered by other periods + overlap
  chunks[(nmax+1):n, 1] <- NA
  last.used <- which(is.finite(chunks[, 1]) & is.na(chunks[, 2]))
  if (l > 0) last.used <- c(last.used, max(last.used) + 1:l)
  first.w <- c(rep(1, nmax-l), rev(blend.w)) # Mixing weights
  mults[last.used, 1] <- utils::tail(first.w, n = length(last.used))

  if (plot) {
    withr::local_par(mar = c(0.3, 0.3, 0.3, 0.3))
    .plotChunks(mults, chunks)
  }
  return(list(chunks = chunks, weights = mults))
}

#' @rdname splitOverlapFixed
#' @export
splitOverlapCustom <- function(n, est.starts, est.ends,
                               blend.starts = NULL, blend.ends = NULL, plot = FALSE) {
  p <- length(est.starts)
  if (is.null(blend.starts)) blend.starts <- est.starts[-1]
  if (is.null(blend.ends)) blend.ends <- est.ends[-p]
  if (length(est.starts) != length(est.ends)) stop("splitOverlapCustom: 'est.starts' should have the same length as 'est.ends'.)")
  if (length(blend.starts) != length(blend.ends)) stop("splitOverlapCustom: 'blend.starts' should have the same length as 'blend.ends'.)")
  if (length(est.starts) - 1 != length(blend.ends)) stop("splitOverlapCustom: 'blend.starts' should have length 'length(est.starts) - 1'.)")
  if (!is.finite(est.starts[1])) est.starts[1] <- 1
  if (!is.finite(est.ends[p])) est.ends[p] <- n
  if (max(est.ends) > n) stop("'est.ends' must come earlier than the end of the full sample.")
  if (min(est.starts) < 1) stop("'est.starts' must come later than the start of the full sample.")
  if (any((blend.starts < est.starts[-1]) & is.finite(blend.starts))) stop("Cannot start blending the series earlier than the estimation sample starts.")
  if (any((blend.ends > est.ends[-p]) & is.finite(blend.ends))) stop("Cannot finish blending the series later than the estimation sample ends.")

  for (i in 1:(p-1)) {
    gap <- est.starts[i+1] - est.ends[i]
    if (gap == 1) blend.starts[i] <- blend.ends[i] <- NA # No gap = no blending
    if (gap > 1) stop("There is a gap between the requested estimation samples -- check the indices.")
  }
  # Chunks = indices of subsamples, mults = weights for linear blending
  chunks <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:p) {
    s <- seq(est.starts[i], est.ends[i])
    chunks[s, i] <- s
  }
  mults <- matrix(NA, nrow = n, ncol = p)
  lw <- lapply(1:(p-1), function(i) if (is.finite(blend.starts[i])) seq(blend.starts[i], blend.ends[i]) else NULL) # List of overlap periods
  if (any(duplicated(unlist(lw)))) stop("Only 2 series can overlap at the same time.")
  for (i in 1:(p-1)) { # Filling in the weights
    if (length(lw[[i]] > 0)) { # If there is a non-zero gap
      w <- seq_along(lw[[i]]); w <- w / (max(w) + 1)
      mults[lw[[i]], i] <- rev(w)
      mults[lw[[i]], i+1] <- w
    }
  }
  # Filling the first column till the first positive weight, and the last after the last positive weight
  if (length(lw[[1]]) > 0) mults[est.starts[1]:(min(lw[[1]])-1), 1] <- 1 else mults[is.finite(chunks[, i]), 1] <- 1
  if (length(lw[[p-1]]) > 0) mults[(max(lw[[p-1]])+1):est.ends[p], p] <- 1 else mults[is.finite(chunks[, p]), p] <- 1
  # Filling the remaining gaps with no values withing the blocks with 1
  if (p > 2) {
    for (i in 2:(p-1)) mults[is.finite(chunks[, i]) & is.na(mults[, i]), i] <- 1
  }
  mults[is.na(mults)] <- 0

  if (plot) {
    withr::local_par(mar = c(0.3, 0.3, 0.3, 0.3))
    .plotChunks(mults, chunks)
  }
  return(list(chunks = chunks, weights = mults))
}

.plotChunks <- function(mults, chunks, gamma = 1) {
  graphics::image(1:nrow(mults), 1:ncol(mults), mults, col = grDevices::gray.colors(40, start = 1, end = 0, gamma = 1), xaxt = "n", yaxt = "n", bty = "n")
  for (i in 1:ncol(mults)) {
    hr <- range(chunks[, i], na.rm = TRUE)
    graphics::rect(hr[1]-0.5, i-0.5, hr[2]+0.5, i + 0.5, border = 1, lty = 2)
  }
}

# pdf("S:/Projets/Modelisation/Modèles/Seasonal Adjustment/Presentations-and-User-Guide/presentation-05-practical/long-overlap.pdf", 5, 1.5)
# splitOverlapFixed(n = 460, 12*12, l = 12, plot = TRUE)
# dev.off()

#' Reconstruct missing values in time series via multiple methods
#'
#' @param x A ts object (not mts)
#' @param robust.weights # Weights assigned to the 5 or 7 methods.
# If given a numeric vector of length 5 (7 for imputePanel7), return the weighted mean of sorted forecasts (inconsistent weighting may occur).
# c(1, 3, 5, 3, 1) is robust; c(0, 1, 2, 1, 0) is very robust; c(0, 0, 1, 0, 0) is the median.
#' @param calendar Passed to diagnoseSeasonality().
#' @param sample.begin Index (e.g. `c(2010, 1)`) to shrink the estimation sample if needed.
#' @param sample.end Index (e.g. `c(2025, 4)`) to shrink the estimation sample if needed.
#' @param name Passed to diagnoseSeasonality() (series name).
#' @param force.positive # If TRUE, does not return negative values; replacing them with half of the minimal positive value.
#' @param return.multi If TRUE, return the 5 series (7 for imputePanel7) without any processing
#' @param attr.seas If TRUE, attach the best X13 seasonal model as an attribute (useful for diagnostics).
#' @param parallel If TRUE, process the series in parallel via parLApplyLB.
#' @param cores Integer: the number of cores for a multi-core cluster.
#' @param trim Trim these fractions of panel for Amelia back- and forecasting (to have a model more relevant for the sample ends)
#' @param calendar Passed to X13.
#' @param seed Seed for Amelia imputation.
#'
#' Methods used in `imputeTS5`:
#' * X13-ARIMA-SEATS forecast and backcast
#' * Kalman filter for a basic structural model (`type = "BSM"` in [stats::StructTS])
#' * Kalman filter for [forecast::auto.arima]
#' * [imputeTS::na_seadec]
#' * [imputeTS::na_seasplit]
#'
#' `imputePanel7` calls `imputeMTS5` on multi-variate time series and two more
#' panel-based methods:
#' * [Amelia::amelia]
#' * [missMDA::imputePCA]
#'
#' @return imputeTS5: if \code{return.multi}, a 5-column matrix, otherwise an imputed time series.
#' imputeMTS5/imputePanel7: if \code{return.multi}, a list of 5- or 7-column matrices, otherwise
#' an imputed MTS.
#' @export
#'
#' @examples
#' x <- seq(0, 7.99, 1/12)
#' set.seed(1) # Generating 3 series wit similar dynamics
#' y1 <- ts(10 + 0.05*x + 5*(x >= 6) - 4*(x == 5) + 40/(x+20)*cos(x*pi*2) +
#'          arima.sim(n = 96, list(ar = 0.8)), start = c(2010, 1), freq = 12)
#' y2 <- ts(3 + 0.05*x + 5*(x >= 6) - 4*(x == 5) + 40/(x+20)*sin(x*pi*2) +
#'          arima.sim(n = 96, list(ar = 0.8)), start = c(2010, 1), freq = 12)
#' y3 <- ts(3 - 0.05*x + 5*(x < 6) - 4*(x == 5) + 40/(x+20)*cos(x*pi*2) +
#'          arima.sim(n = 96, list(ar = 0.8)), start = c(2010, 1), freq = 12)
#' y2[1:11] <- NA # Creating missing observations to impute
#' y3[87:96] <- NA
#' y <- ts.union(y1, y2, y3)
#' yr <- range(y, na.rm = TRUE)
#' plot(y, plot.type = "single", ylim = yr, col = 1:3)
#' y2i <- imputeTS5(y2, return.multi = TRUE) # Takes 10--20 seconds
#' plot(y2i, plot.type = "single", col = 5:1)
#' legend("top", colnames(y2i), ncol = 5, col = 5:1, lwd = 1)
#' \donttest{
#' # This takes 20 seconds per series, i.e. 2 minutes per this section
#' yi <- imputeMTS5(y)
#' plot(yi, plot.type = "single", ylim = yr, ylab = "y", col = 1:3)
#' yp <- imputePanel7(y)
#' par(new = TRUE)
#' plot(yp, plot.type = "single", ylab = "", ylim = yr, col = 1:3, lty = 2)
#' }
imputeTS5 <- function(x, robust.weights = NULL, calendar = NULL,
                      sample.begin = NULL, sample.end = NULL, name = "",
                      force.positive = FALSE,
                      return.multi = FALSE, attr.seas = FALSE
) {
  y <- x
  freq <- stats::frequency(y)
  r <- getRangeVec(y)
  if (is.null(sample.begin)) sample.begin <- stats::start(y)
  if (is.null(sample.end)) sample.end <- stats::end(y)
  first.valid <- time2ind(stats::time(y)[r["first"]], freq) # Index of the first observed data point
  last.valid  <- time2ind(stats::time(y)[r["last"]], freq)
  y.good <- stats::window(y, start = first.valid, end = last.valid)
  # y can contain missing values, so we find the first non-missing
  npery1 <- round(ind2time(first.valid, freq)*freq) # Integer period number at the beginning
  npery2 <- round(ind2time(last.valid, freq)*freq)
  nperwanted1 <- round(ind2time(sample.begin, freq)*freq) # Desired periods at the beginning
  nperwanted2 <- round(ind2time(sample.end, freq)*freq)
  nhead <- npery1 - nperwanted1
  ntail <- nperwanted2 - npery2
  if (nhead > 0 & ntail > 0) {
    y <- stats::ts(c(rep(NA, nhead), y.good, rep(NA, ntail)), frequency = freq, start = sample.begin)
  } else if (nhead > 0) {
    y <- stats::ts(c(rep(NA, nhead), y.good), frequency = freq, start = sample.begin)
  } else if (ntail > 0) {
    y <- stats::ts(c(y.good, rep(NA, ntail)), frequency = freq, start = first.valid)
  }

  # Should any imputation be done at all?
  rnew <- getRangeVec(y)
  do.head <- rnew["first"] > 1
  do.tail <- rnew["last"] != length(y)
  do.mid  <- rnew["mid.na.count"] > 0

  # If there is nothing to do, the imputed data set is equal to the complete one
  y1.mod <- if (do.head | do.tail | do.mid | attr.seas) {
    diagnoseSeasonality(y, calendar = calendar, bfcast.tails = TRUE, name = name, verbose = 1, plot.file = NA)
  } else NULL
  if (do.head | do.tail | do.mid) {
    # 5 different methods
    y1 <- y1.mod$series[, "predicted"]
    y2 <- imputeTS::na_kalman(y, model = "StructTS", type = "BSM")
    y3 <- imputeTS::na_kalman(y, model = "auto.arima", max.p = 3, max.q = 3, max.d = 1, max.P = 1, max.Q = 1, stepwise = FALSE)
    y4 <- tryCatch(imputeTS::na_seasplit(y, algorithm = "ma"), error = function(e) imputeTS::na_seasplit(y, algorithm = "locf"))
    y5 <- imputeTS::na_seadec(y)
  } else {
    y1 <- y2 <- y3 <- y4 <- y5 <- y
    warning(paste0("The series '", name, "' does not have any missing values. Wrapping the output as expected..."))
  }

  y.imp <- stats::ts.union(y1, y2, y3, y4, y5)
  colnames(y.imp) <- c("X13", "StructTS", "auto.arima", "seasplit", "seadec")

  if (any(y.imp <= 0)) { # One can replace negative values with positive ones
    if (force.positive) {
      min.pos <- min(y.imp[y.imp > 0])
      y.imp[y.imp <= 0] <- min.pos / 2
      cat("Non-positive values replaced with the half of the smallest positive one.")
    } else warning("The series are positive, yet the imputation contains negative values. Consider using force.positive = TRUE.")
  }

  if (return.multi) { # No collapsing
    res <- y.imp
  } else {
    if (is.null(robust.weights)) res <- apply(y.imp, 1, function(x) stats::weighted.mean(x, w = c(1, 0.5, 0.5, 0.5, 0.5))) else {
      if (!is.numeric(robust.weights)) stop("robust.weights should be numeric")
      if (length(robust.weights) != 5) stop("robust.weights should have length 5; try c(1, 3, 5, 3, 1) or c(0, 1, 2, 1, 0) for tri-mean.")
      res <- t(apply(y.imp, 1, sort)) # Sorting forecasts from lowest to highest
      res <- apply(res, 1, function(x) stats::weighted.mean(x, w = robust.weights))
    }
    attributes(res) <- attributes(y1)
  }

  if (attr.seas) attr(res, "seasonal") <- y1.mod
  attr(res, "original") <- y
  return(res)
}

#' @rdname imputeTS5
#' @export
imputeMTS5 <- function(x,
                       robust.weights = NULL, calendar = NULL,
                       sample.begin = NULL, sample.end = NULL,
                       force.positive = FALSE,
                       parallel = FALSE, cores = 4, return.multi = FALSE
) {
  # Preparing the arguments to export if parallel is used
  imputei <- function(i) {
    tryCatch(imputeTS5(x = x[, i], robust.weights = robust.weights, calendar = calendar,
                       sample.begin = sample.begin, sample.end = sample.end,
                       force.positive = force.positive, name = colnames(x)[i],
                       return.multi = return.multi, attr.seas = return.multi),
             error = function(e) return(e))
  }
  if ("mts" %in% class(x)) {
    if (parallel & cores > 1) {
      cl <- parallel::makeCluster(cores)
      parallel::clusterEvalQ(cl, library(StatecCMP))
      parallel::clusterExport(cl, c("x", "robust.weights", "calendar", "sample.begin", "sample.end", "force.positive"), envir = environment())
      y.out <- parallel::parLapplyLB(cl, X = 1:ncol(x), imputei)
      parallel::stopCluster(cl); rm(cl)
    } else y.out <- lapply(1:ncol(x), imputei)

    # For easier debugging
    if (any(err <- sapply(y.out, function(x) "error" %in% class(x)))) {
      warning(paste0("The MTS routine is broken! Check the output.\nIf the series look OK, tell Michel or Andrei to add a special case.\nProblematic series: ", paste0(colnames(x)[err], collapse = ", ")))
      print(y.out[err])
      return(y.out)
    }

    if (!return.multi) { # Normal multi-column output
      y.out <- do.call(stats::ts.union, y.out)
      colnames(y.out) <- colnames(x)
    } else { # Individual estimation results
      names(y.out) <- colnames(x)
    }
    return(y.out)
  } else stop("x must be an MTS.")
}

#' @rdname imputeTS5
#' @export
imputePanel7 <- function(x, trim = c(0.25, 0.25),
                       robust.weights = NULL, calendar = NULL,
                       sample.begin = NULL, sample.end = NULL,
                       force.positive = FALSE,
                       parallel = FALSE, cores = 4, seed = 1,
                       return.multi = FALSE
) {
  if (is.null(colnames(x))) colnames(x) <- paste0("X", 1:ncol(x))
  miss.inds <- apply(x, 2, function(x) any(!is.finite(x)))
  y.imp <- imputeMTS5(x = x, return.multi = TRUE, calendar = calendar, sample.begin = sample.begin, sample.end = sample.end, force.positive = force.positive, parallel = parallel, cores = cores)
  seas.mods <- lapply(y.imp, function(x) attr(x, "seasonal"))
  y.series <- do.call(stats::ts.union, lapply(seas.mods, function(x) x$series[, "original"])) # Series of proper length
  colnames(y.series) <- colnames(x)

  # Amelia for data with seasonality
  x.pan <- expand.grid(Date = ts2date(y.series), Unit = colnames(y.series))
  x.pan$y <- NA
  y.s <- scale(y.series) # Re-scaling the variables to have reasonable imputations
  for (v in colnames(x)) x.pan$y[x.pan$Unit == v] <- y.s[, v]
  x.pan$cycle <- factor(as.numeric(substr(as.character(x.pan$Date), 6, 7)), levels = 1:12)
  # Imputing the beginning of sample trimmed data, and the end of sample only with trimmed data
  x.range <- getRange(y.series, date.format = "index")
  do.head <- !all(x.range$first == 1)
  do.mid <- !all(x.range$na.runs == 0)
  do.tail <- !all(x.range$last == nrow(x))
  d <- ts2date(y.series)

  if (do.head) {
    cond1 <- x.pan$Date <= stats::quantile(x.pan$Date, 1 - trim[2], type = 1)
    am1 <- x.pan[cond1, ]
    # Safety check
    if (!any(is.na(am1))) {
      do.head <- FALSE
    } else {
      set.seed(seed)
      imp.amel1 <- suppressWarnings(Amelia::amelia(am1,  ts = "Date", cs = "Unit", noms = "cycle",  splinetime = 1, intercs = TRUE, logs = if (force.positive) "y" else NULL, leads = "y"))
      set.seed(seed)
      imps1 <- suppressWarnings(Amelia::tscsPlot(imp.amel1, plotall = TRUE, var = "y"))
      grDevices::dev.off() # We do not really need the plot; we need the imputed values that are plotted
    }
  }
  if (do.tail) {
    cond2 <- x.pan$Date >= stats::quantile(x.pan$Date, trim[1], type = 1)
    am2 <- x.pan[cond2, ]
    if (!any(is.na(am2))) {
      do.tail <- FALSE
    } else {
      set.seed(seed)
      imp.amel2 <- suppressWarnings(Amelia::amelia(am2,  ts = "Date", cs = "Unit", noms = "cycle",  splinetime = 1, intercs = TRUE, logs = if (force.positive) "y" else NULL, lags = "y"))
      set.seed(seed)
      imps2 <- suppressWarnings(Amelia::tscsPlot(imp.amel2, plotall = TRUE, var = "y"))
      grDevices::dev.off() # We do not really need the plot; we need the imputed values that are plotted
    }
  }
  if (do.mid) {
    missall <- apply(y.series, 1, function(x) all(is.na(x))) # Fully missing end observations are useless here
    last.head <- if (missall[1]) which(diff(missall) != 0)[1] + 1 else 1
    first.tail <- if (missall[length(missall)]) nrow(y.series) - which(diff(rev(missall)) != 0)[1] else nrow(y.series)
    cond3 <- as.numeric(factor(x.pan$Date)) %in% last.head:first.tail
    am3 <- x.pan[cond3, ]
    if (!any(is.na(am2))) {
      do.mid <- FALSE
    } else {
      set.seed(seed)
      imp.amel3 <- suppressWarnings(Amelia::amelia(am3,  ts = "Date", cs = "Unit", noms = "cycle",  splinetime = 1, intercs = TRUE, logs = if (force.positive) "y" else NULL, lags = "y", leads = "y"))
      set.seed(seed)
      imps3 <- suppressWarnings(Amelia::tscsPlot(imp.amel3, plotall = TRUE, var = "y"))
      grDevices::dev.off() # We do not really need the plot; we need the imputed values that are plotted
    }
  }
  # Putting the middle values first because the NAs might be overlapping
  x.pan$yimp <- x.pan$y
  if (do.mid) x.pan$yimp[cond3]  <- apply(imps3, 1, mean)
  if (do.head) x.pan$yimp[cond1] <- apply(imps1, 1, mean)
  if (do.tail) x.pan$yimp[cond2] <- apply(imps2, 1, mean)

  imp.amelia <- matrix(x.pan$yimp, nrow = nrow(y.series), ncol = ncol(y.series))
  imp.amelia <- stats::ts(imp.amelia, start = stats::start(y.series), freq = stats::frequency(y.series))
  colnames(imp.amelia) <- colnames(x)
  # Restoring the scale
  imp.amelia <- sweep(imp.amelia, 2, attr(y.s, "scaled:scale"), "*")
  imp.amelia <- sweep(imp.amelia, 2, attr(y.s, "scaled:center"), "+")

  # PCA-based for differenced data without seasonality
  # PCA in levels usually works poorly
  # Assuming positivity
  y.sa <- do.call(stats::ts.union, lapply(seas.mods, function(x) x$series[, "adjusted"]))
  if (all(y.series[is.finite(y.series)] > 0) | force.positive) {
    y.seas <- do.call(stats::ts.union, lapply(seas.mods, function(x) x$series[, "adjfac"]))
    y.mult <- sapply(seas.mods, "[[", "multiplicative")
    y.sa[is.na(y.series)] <- NA # Restoring the middle gaps
    y.sa[y.sa <= 0] <- NA
    y.gr <- myDiff(log(y.sa)) # Seasonally adjusted growth rates
    set.seed(seed)
    k <- ncol(y.gr)
    min.ncp <- 1 + (k > 3) + (k > 7)
    max.ncp <- 1 + (k > 2) + (k > 3)*ceiling(sqrt(abs(k-3)))
    ncp <- missMDA::estim_ncpPCA(y.gr, ncp.min = min.ncp, ncp.max = max.ncp)
    # plot(names(ncp$criterion), ncp$criterion)
    gr.pca <- tryCatch(missMDA::imputePCA(y.gr, ncp = ncp$ncp, method = "EM")$completeObs, error = function(e) missMDA::imputePCA(y.gr, ncp = ncp$ncp)$completeObs)
    # par(mfrow = c(4, 4), mar = c(0.1, 0.1, 0.1, 0.1))
    # for (i in seq_len(ncol(y.series))) {
    #   plot(y.gr[, i], bty = "n", xaxt = "n", yaxt = "n")
    #   points(y.pca[, i], pch = 16, cex = 0.8, col = 2)
    #   mtext(colnames(y.series)[i], 3, line = -2)
    # }
    imp.pca <- y.sa
    for (i in 1:ncol(imp.pca)) {
      do.head <- x.range$first[i] > 1
      do.mid <- x.range$na.runs[i] > 0
      do.tail <- x.range$last[i] < nrow(x)
      do.any <- do.head | do.mid | do.tail
      if (do.head) { # Treating missingness at the beginning---computing levels
        for (j in x.range$first[i]:1) imp.pca[j, i] <- imp.pca[j+1, i] / exp(gr.pca[j+1, i]) # Today divided by today's growth is yesterday
      }
      if (do.tail) { # Treating missingness at the end---computing levels
        for (j in (x.range$last[i]+1):nrow(imp.pca)) imp.pca[j, i] <- imp.pca[j-1, i] * exp(gr.pca[j, i])
      }
      if (do.mid) { # Treating missingness in the middle by extrapolating forward
        js <- findMiddleGaps(imp.pca[, i])
        for (j in js) imp.pca[j, i] <- imp.pca[j-1, i] * exp(gr.pca[j, i])
      }
      if (do.any) {
        imp.pca[, i] <- if (y.mult[i]) imp.pca[, i] * y.seas[, i] else imp.pca[, i] + y.seas[, i]
      } else { # There is nothing to impute -- restore the original unadjusted
        imp.pca[, i] <- y.series[, i]
      }
    }
  } else { # Not applying the correction, returning NA with TS attributes
    imp.pca <- y.sa
    imp.pca[1:nrow(imp.pca), 1:ncol(imp.pca)] <- NA
    warning("imputePanel7: some of the series contain negative values. Skipping growth-rate-based imputation via PCA.")
  }

  y.imp2 <- lapply(1:length(y.imp), function(i) stats::ts.union(y.imp[[i]], imp.amelia[, i], imp.pca[, i]))
  if (return.multi) return(y.imp2)

  # Combining the results
  if (is.null(robust.weights)) {
    res <- do.call(cbind, lapply(y.imp2, function(y) apply(y, 1, function(x) stats::weighted.mean(x, w = c(1, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5), na.rm = TRUE))))
  } else {
    if (!is.numeric(robust.weights)) stop("robust.weights should be numeric")
    if (length(robust.weights) != 7) stop("robust.weights should have length 7; try c(1, 3, 5, 7, 5, 3, 1) or c(0, 1, 2, 2, 2, 1, 0).")
    res <- do.call(cbind, lapply(y.imp2, function(y) apply(y, 1, function(x) {
      x <- sort(x, na.last = TRUE)
      stats::weighted.mean(x, w = robust.weights, na.rm = TRUE)
    })))
  }
  attributes(res) <- attributes(y.series)
  return(res)
}

