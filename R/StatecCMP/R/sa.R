################################################################################
# R file containing the functions necessary for seasonal adjustment
#
# Author: Andrei V. Kostyrka
# Entity: Statec
# Date: 2022-01 -- 2022-12
################################################################################


# #' Trim a calendar to pass to X13
# #'
# #' @param ts A ts object to be adjusted
# #' @param cal A calendar ts
# #' @param file The name of the file to write. If NULL, prints to console.
# #' @param years.add Integer: how many years to append at the end. X13-RegARIMA uses 1.
# #' @param periods.add Integer: how many periods to append at the end.
# #'
# #' @return Nothing.
# #' @export
# #'
# #' @examples
# genX13cal <- function(ts, cal, file = NULL, years.add = 1, periods.add = 0) {
#   .Deprecated(msg = "This function was only used for testing with bare-bones X13.")
#   f <- stats::frequency(ts)
#   if (!(f %in% c(4, 12))) stop("Only monthly or quarterly data are expected.")
#   if (!isTRUE(f == stats::frequency(cal))) stop("The frequency of the time series and the calendar must be the same!")
#   x <- stats::ts.union(ts, cal)
#   r <- getRangeVec(x[, 1]) # The ranges of the time series in the merged data set
#   index.incr <- years.add * f + periods.add
#   t <- stats::time(x)
#   x <- stats::window(x, start = t[r["first"]], end = t[r["last"] + index.incr]) # Trimming the calendar
#   x <- x[, -1] # Droppint the time time series
#   t <- stats::time(x)
#   yr <- floor(as.numeric(t))
#   prd <- stats::cycle(t)
#   m <- as.matrix(data.frame(yr, prd, x))
#   if (is.null(file)) {
#     for (i in 1:nrow(m)) cat(paste0(m[i, ], collapse = " "), "\n")
#   } else utils::write.table(m, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE)
#   return(invisible(NULL))
# }

# #' Generate an X13 variable file
# #'
# #' @param x A matrix or a data frame.
# #' @param file The name of the file to write. If NULL, prints to console.
# #'
# #' @return Nothing.
# #' @export
# #'
# #' @examples
# genX13var <- function(x, file = NULL) {
#   .Deprecated(msg = "This function was only used for testing with bare-bones X13.")
#   if (is.null(file)) cat(x, sep = "\n") else utils::write.table(x, file = file, quote = FALSE, row.names = FALSE, col.names = FALSE)
# }

#' Get the crucial seasonal adjustment statistics
#'
#' @param x An object of class seasonal::seas.
#'
#' @return A list with named values.
#' @export
#'
#' @examples
#' xs <- seasonal::seas(datasets::AirPassengers, x11 = "")
#' getStat(xs)
getStat <- function(x) {
  if (is.null(x)) return(list(log = "", arima = "",
                              Noutliers = -1, N.AO = -1, N.LS = -1, N.TC = -1,
                              outliers = "", datesAO = "", datesLS = "", datesTC = "", M = rep(999, 11), Q_M2 = 999))
  if (class(x)[1] != "seas") stop("The main argument must be a 'seas' object.")
  mod <- x$model$arima$model
  outl <- x$model$regression$variables
  outl <- setdiff(outl, "const")
  outl <- outl[grepl("^ao|^ls|^tc", outl)]
  if (length(outl) > 0) {
    nout <- length(outl)
    noutAO <- sum(grepl("ao", outl))
    noutLS <- sum(grepl("ls", outl))
    noutTC <- sum(grepl("tc", outl))
    outAO <- paste0(outl[grep("ao", outl)], collapse = ",")
    outLS <- paste0(outl[grep("ls", outl)], collapse = ",")
    outTC <- paste0(outl[grep("tc", outl)], collapse = ",")
  } else {
    outAO <- outLS <- outTC <- ""
    nout <- noutAO <- noutLS <- noutTC <- 0
  }
  uselog <- isTRUE(seasonal::transformfunction(x) == "log")
  m <- as.numeric(x$udg[paste0("f3.m", sprintf("%02d", 1:11))])
  names(m) <- paste0("M", 1:11)
  q <- x$udg[grepl("f3\\.q", names(x$udg))]

  ic <- as.numeric(x$udg[c("aic" , "aicc", "bic")])
  names(ic) <- c("AIC", "AICc", "BIC")

  list(log = uselog, arima = mod,
       Noutliers = nout, N.AO = noutAO, N.LS = noutLS, N.TC = noutTC,
       outliers = paste0(outl), datesAO = outAO, datesLS = outLS, datesTC = outTC,
       M = m, Q_M2 = as.numeric(q[2]), IC = ic)
}


# makeStatDF <- function(x) {
#   xl <- lapply(x, getStat)
#   data.frame(
#     log = do.call(c, lapply(xl, "[[", "log")),
#     arima = do.call(c, lapply(xl, "[[", "arima")),
#     Noutliers = do.call(c, lapply(xl, "[[", "Noutliers")),
#     outliers =  do.call(c, lapply(xl, "[[", "outliers")),
#     M = do.call(rbind, lapply(xl, "[[", "M")),
#     Q_M2 = do.call(c, lapply(xl, "[[", "Q_M2"))
#   )
# }

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

#' Diagnose seasonality, adjust the series, print and visualise the results
#'
#' @param x A numeric vector of univariate time series of class "ts".
#' @param calendar An object of class "ts" containing calendar regressors: Monday, ..., Saturday, LeapYear, WorkingDays. If NULL, the default X13 working-day calendar is used.
#' @param name A string with the name of the series for printing.
#' @param bfcast.tails Logical: if x starts or ends with NA, produce the backcasts and forecasts for those periods?
#' @param sample.begin An integer vector of length 2 denoting the beginning of backcast period. By default, no backcasts are used.
#' @param sample.end An integer vector of length 2 denoting the end of forecast period. By default, predicts the seasonality adjustment factor for 1 year out of sample.
#' @param force.multiplicative If TRUE, forces multiplicative adjustment.
#' @param m7.threshold Numeric between 0 and 3: the procedure will produce SA series if the M7 statistic is less than this threshold. A value of 3 means that adjustment is always done.
#' @param q.threshold Numeric between 0 and 3 (like m7.threshold): changes the contents of the information messages depending on the acceptability threshold (purely cosmetic).
#' @param transform.aicdiff Numeric: the AICc difference between the additive (linear) and multiplicative (log) model; use log if AICc(nolog) - AICc(log) < aicdiff; by default, prefers multiplicative models
#' @param tradingdays.aicdiff Numeric: the AICc difference for the TD or WD dummies to be added into the model (if negative, prefers more parsimonious models)
#' @param plot.file A string containing the path to the output PNG file. If NULL, plot to the current device. Using NA prevents any plot from being created.
#' @param seasonal.plot Logical: replace box plot with series plot?
#' @param min.outliers An integer: the number of models in which a candidate outlier must appear to be included in the AIC test.
#' If 1, test the model for selection with a common list of all outliers (appearing at least once in the regressions with 6 TD, 1 TD, and no TD). If 2, test with outliers
#' existing in at least 2 out of 3 specifications. Cannot be greater than 3.
#' @param verbose Logical or integer. TRUE or 1 = basic output, 2  = detailed output.
#'
#' Since seasonal adjustment via X13 implies estimation with TRAMO-like methods that allow missingness,
#' the returned series "predicted" will contain fitted values in places where the original data had middle gaps.
#'
#' @return A list: logical "seasonality" indicating whether there is stable seasonality,
#' logical "multiplicative" indicating whether the adjustment model is multiplicative,
#' a list "quality" containing brief diagnostic information,
#' mts "series" containing the original, adjustment factor (with forecast), adjusted and forecast series,
#' and the output of "seas" for the best model.)
#' @export
#'
#' @examples
#' # Simulating ARIMA with trend, decaying seasonality, one level shift and one additive outlier
#' x <- seq(0, 7.99, 1/12)
#' set.seed(1)
#' y <- ts(10 + 0.05*x + 5*(x >= 6) - 4*(x == 5) + 40/(x+20)*cos(x*pi*2) +
#'         arima.sim(n = 96, list(ar = 0.8)), start = c(2010, 1), freq = 12)
#' diagnoseSeasonality(y)
diagnoseSeasonality <- function(x, calendar = NULL, name = "Y",
                                bfcast.tails = FALSE,
                                sample.begin = NULL, sample.end = NULL,
                                force.multiplicative = FALSE,
                                m7.threshold = 1, q.threshold = 1,
                                transform.aicdiff = -2, tradingdays.aicdiff = 0,
                                plot.file = NULL, seasonal.plot = TRUE,
                                min.outliers = 1, verbose = 2
                                ) {
  # Loading the functions that are not exported
  update.seas <- utils::getFromNamespace("update.seas", "seasonal")
  predict.seas <- utils::getFromNamespace("predict.seas", "seasonal")

  freq <- stats::frequency(x)
  if (!(freq %in% c(4, 12))) stop("The data are not monthly or quarterly, aborting.")
  if (!is.null(calendar)) {
    if (freq != stats::frequency(calendar)) stop("The frequency of the data is not equal to that of the calendar, aborting.")
    if (verbose > 1) cat("Using the user-supplied calendar (presumably Luxembourgish).\n")
  } else {
    if (verbose > 1) cat("Using the default X13 calendar (non-Luxembourgish).\n")
  }

  ef <- function(e) return(NULL) # For try-catching and debugging
  efv <- function(e) {print(e); return(NULL)} # For try-catching and debugging
  ################################################################
  # Helper functions to combine outlier types
  # Choosing the same set of outliers that were found in at least 2 models
  OTMed2 <- function(x) { # Outlier type median
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
  OTMed3 <- function(x) { # Outlier type median
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
  printSum <- function(seasmod, txt = "-- Estimated a seasonal model", ndigits = 3) {
    n <- as.integer(seasmod[["udg"]]["nobs"])
    nef <- as.integer(seasmod[["udg"]]["nefobs"])
    avgAIC <- as.numeric(seasmod[["udg"]]["aicc"]) / nef
    fullAIC <- avgAIC * n
    cat(txt, " (obs=", nef, "). AICc/obs = ", sprintf(paste0("%1.", ndigits, "f"), avgAIC), ", AICcS = ", sprintf("%1.1f", fullAIC), ".\n", sep = "")
  }


  ######################################
  tx <- stats::time(x)
  # Assume that the beginning of the sample if the beginning of the input series
  if (is.null(sample.begin)) sample.begin <- stats::start(tx)
  if (is.null(sample.end))   sample.end   <- stats::end(tx)

  full.range <- getRangeVec(x)
  # If sample.begin comes later than the first valid observation of this.ser, trim this.ser
  # If sample.end comes earlier than the last valid observation of this.ser, do the same
  full.dates <- time2ind(stats::time(x), freq = freq)
  good.start <- if (ind2time(sample.begin, freq) - tx[full.range["first"]] > -0.001) sample.begin else full.dates[full.range["first"], ]
  good.end   <- if (tx[full.range["last"]] - ind2time(sample.end, freq) > -0.001) sample.end else full.dates[full.range["last"], ]
  this.ser <- stats::window(x, start = good.start, end = good.end)
  # this.ser is the series for analysis without missingness at the tails
  this.range <- getRangeVec(this.ser)
  n <- length(this.ser)

  linelen <- 60
  if (verbose > 1) printSym("=", linelen)
  if (verbose > 1) cat("Testing and adjusting seasonality in", name, "via the recommended procedure with parameter auto-selection.\n") else if (verbose > 0) cat("Seasonality adjustment for ", name, " via X13.\n", sep = "")
  if (verbose > 1) printSym("=", linelen)

  head.years <- max(min(stats::time(this.ser)) - ind2time(sample.begin, freq), 0) # Years before the first observation for forecasting / adjustment factor generation
  tail.years <- max(ind2time(sample.end, freq) -  max(stats::time(this.ser)), 0) # Years for forecasting / adjustment factor generation
  nbcast <- round(head.years * freq) # How many points are missing before the beginning of the forecasting range
  nfcast <- round(tail.years * freq) # How many points are missing till the end of the forecasting range
  do.bcast <- if (nbcast > 0) "yes" else "no"
  do.fcast <- if (nfcast > 0) "yes" else "no"
  if (!bfcast.tails) {
    do.fcast <- do.bcast <- "no"
    nbcast <- nfcast <- 0
  }

  # JDemetra-like procedure with all 5 weekday dummies + Easter + leap year
  # A function without regression variables, xreg, aictest
  if (!is.null(calendar)) {
    this.cal <- stats::window(calendar, start = min(stats::time(this.ser)) - head.years, end = max(stats::time(this.ser)) + tail.years)
    td.regs <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
    wd.regs <- c("WorkingDays") # These variables are used in the calendar for Luxembourg
    td.vars <- c("easter[8]", "lpyear")
    wd.vars <- c("easter[8]", "lpyear")
    td.cal <- this.cal[, td.regs]
    wd.cal <- this.cal[, wd.regs]
    td.types <- rep("td", 6)
    wd.types <- "td"
  } else {
    td.vars <- c("easter[8]", "tdnolpyear", "lpyear")
    wd.vars <- c("easter[8]", "td1nolpyear", "lpyear")
    td.cal <- wd.cal <- td.types <- wd.types <- NULL
  }
  nod.vars <- c("easter[8]", "lpyear")

  miss.inds <- which(!is.finite(this.ser)) # Used at the end for imputed return
  if (this.range["mid.na.count"] > 0) {
    warning("The series contain internal NAs. Plugging in interpolated values and adding AOs for each point.")
    this.ser2 <- stats::approx(x = which(is.finite(this.ser)), y = stats::na.omit(as.numeric(this.ser)), xout = seq_along(this.ser))$y
    this.ser[1:length(this.ser)] <- this.ser2
    miss.years <- myYear(this.ser)
    miss.cycle <- stats::cycle(this.ser)
    miss.vars <- paste0("ao", miss.years, ".", miss.cycle)[miss.inds]
    td.vars <- c(td.vars, miss.vars)
    wd.vars <- c(wd.vars, miss.vars)
    nod.vars <- c(nod.vars, miss.vars)
  }

  # Step 1: choosing the transformation; using common outliers
  if (any(this.ser <= 0)) {
    trans.fun <- "none"
    if (force.multiplicative) {
      warning("The series contains non-positive values. Multplicative adjustment impossible, switching to additive.")
      force.multiplicative <- FALSE
    }
    a.td <- tryCatch(seasonal::seas(x = this.ser, x11 = "", outlier.types = c("ao", "ls", "tc"),
                                    regression.variables = td.vars, xreg = td.cal, regression.usertype = td.types,
                                    automdl.maxorder = c(3, 1), estimate.tol = 1e-7,
                                    regression.aictest = NULL, transform.function = trans.fun,
                                    x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                                    x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                                    forecast.save = c("bct", "fct")), error = efv)
    printSum(a.td, "-- Model: additive, 6 TD + LY + Easter")
  } else { # All positive values, selection possible
    if (verbose > 0) cat("Estimating 2 transformation models (additive / multiplicative)...\n")
    a.td.log <- tryCatch(seasonal::seas(x = this.ser, x11 = "", outlier.types = c("ao", "ls", "tc"),
                                        regression.variables = td.vars, xreg = td.cal, regression.usertype = td.types,
                                        automdl.maxorder = c(3, 1), estimate.tol = 1e-7,
                                        regression.aictest = NULL, transform.function = "log",
                                        x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                                        x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                                        forecast.save = c("bct", "fct")), error = efv)
    if (verbose > 1) printSum(a.td.log, "-- Model: multipl., 6 TD + LY + Easter")
    if (!force.multiplicative) { # If the transformation selection is automatic
      a.td.nolog <- tryCatch(seasonal::seas(x = this.ser, x11 = "", outlier.types = c("ao", "ls", "tc"),
                                            regression.variables = td.vars, xreg = td.cal, regression.usertype = td.types,
                                            automdl.maxorder = c(3, 1), estimate.tol = 1e-7,
                                            regression.aictest = NULL, transform.function = "none",
                                            x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                                            x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                                            forecast.save = c("bct", "fct")), error = efv)
      if (verbose > 1) printSum(a.td.nolog, "-- Model: additive, 6 TD + LY + Easter")
      extrOut <- function(x) tryCatch(tolower(seasonal::outlier(x)), error = function(e) {ret <- rep(NA, length(this.ser)); attributes(ret) <- attributes(this.ser); ret})
      log.outlier.ts <- data.frame(log = extrOut(a.td.log), nolog = extrOut(a.td.nolog))
      log.outlier.ts.med <- lapply(1:nrow(log.outlier.ts), function(i) OTMed2(as.matrix(log.outlier.ts)[i, ])) # Creating only one outlier per period
      log.outlier.ts.med.conc <- unlist(lapply(log.outlier.ts.med, attr, "concordant"))
      log.outlier.ts.med <- unlist(log.outlier.ts.med)
      log.outl.time <- paste0(myYear(this.ser), ".", stats::cycle(this.ser))
      log.outlier.tab <- paste0(log.outlier.ts.med, log.outl.time)
      log.outl.reg <- log.outlier.tab[!grepl("^NA", log.outlier.tab)]

      if (length(log.outl.reg) > 0 & any(!log.outlier.ts.med.conc)) { # There are some outlier regressors and at least one is not concordant across the two models
        if (verbose > 1) cat("Outlier regressors used in the transformation choice: ", paste0(log.outl.reg, collapse = ", "), "\n", sep = "")
        if (verbose > 1) cat("Estimation with identical outlier regressors for comparability by AICc...\n")
        # If re-estimation is successful, then, choose the best model based on auto-AIC
        a.td.auto <- tryCatch(update.seas(a.td.nolog, regression.variables = unique(c(td.vars, log.outl.reg)), outlier = NULL, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = efv)
        if (verbose > 1) printSum(a.td.auto, "-- Model: auto-tr., 6 TD + LY + Easter")
      } else {
        if (verbose > 1) cat("All outlier variables are identical under 2 tranformations, skipping re-estimation.\n", sep = "")
        a.td.auto <- tryCatch(update.seas(a.td.nolog, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = efv)
      }

      if (!is.null(a.td.auto)) {
        if (length(log.outl.reg) > 0 & any(!log.outlier.ts.med.conc)) {
          if (verbose > 1) cat("Estimation with common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
        }
        if (!is.null(a.td.auto$udg["aictest.trans.aicc.nolog"])) {
          if (verbose > 1) cat("With comparable regressors and sample lengths, AICc(nolog) = ",
                               sprintf("%1.1f", as.numeric(a.td.auto$udg["aictest.trans.aicc.nolog"])),
                               ", AICc(log) = ", sprintf("%1.1f", as.numeric(a.td.auto$udg["aictest.trans.aicc.log"])), ".\n", sep = "")
        }
        # Choosing the transformation and setting the default model with 6 TDs
        trans.fun <- if (isTRUE(a.td.auto$udg["aictrans"] != "No transformation") | force.multiplicative) "log" else "none"
        a.td <- if (trans.fun == "log") a.td.log else a.td.nolog
      } else { # Something went wrong with the same set of outliers -- possibly too many -- reduce the number
        if (verbose > 0) cat("! Re-trying with fewer outliers (existing in both models).\n")
        log.outlier.tab <- table(c(a.td.log$model$regression$variables, a.td.nolog$model$regression$variables))
        log.outlier.tab <- log.outlier.tab[log.outlier.tab == 2]
        log.outl.reg <- sort(setdiff(names(log.outlier.tab), c("easter[8]", "lpyear", "const", "tdnolpyear", "td1coef", "td1nolpyear")))
        if (verbose > 1) cat("Revised outlier regressors used in the transformation choice: ", paste0(log.outl.reg, collapse = ", "), "\n", sep = "")

        a.td.auto <- tryCatch(update.seas(a.td.nolog, regression.variables = unique(c(td.vars, log.outl.reg)), outlier = NULL, transform.function = "auto", transform.aicdiff = transform.aicdiff), error = efv)
        if (!is.null(a.td.auto)) {
          if (verbose > 1) printSum(a.td.auto, "-- Model: auto-tr., 6 TD + LY + Easter, fewer common outliers")
          if (verbose > 1) cat("Re-estimation with fewer common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
          trans.fun <- if (isTRUE(a.td.auto$udg["aictrans"] != "No transformation") | force.multiplicative) "log" else "none"
        } else {
          if (verbose > 0) cat("! Re-estimation with fewer common outliers was unsuccessful.\n! Doing selection with the original models.\n! The results might be sub-optimal because the outliers are different.\n")
          log.mod.list <- list(a.td.log, a.td.nolog)
          aiccs <- unname(sapply(log.mod.list, function(x) if (!is.null(x)) x[["udg"]]["aicc"] else Inf))
          names(aiccs) <- c("log", "nolog")
          eff.len <- sapply(log.mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 0)
          # Adjusting the AICc based on the number of observations
          aiccs <- tryCatch(aiccs / eff.len * n, error = function(e) return(c(log = -Inf, nolog = Inf))) # If one of the models fails, assume that log is better, treat this case latter
          trans.fun <- if (aiccs["nolog"] - aiccs["log"] < transform.aicdiff) "none" else "log"
        }
        a.td <- if (trans.fun == "log") a.td.log else a.td.nolog
      }
    } else {
      trans.fun <- "log"
      a.td <- a.td.log
    }
  }
  if (!exists("a.td")) a.td <- NULL
  if (is.null(a.td)) {
    warning("The transformation selection failed with default (generous) settings. Falling back to the most parsimonious settings and retrying.")
    a.td <- tryCatch(seasonal::seas(x = this.ser, x11 = "", regression.variables = NULL,
                                        automdl.maxorder = c(2, 1), transform.function = "auto",
                                        x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                                        x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                                        forecast.save = c("bct", "fct")), error = efv)
    if (is.null(a.td)) {
      warning("Even the most parsimonious auto-model could not be estimated. The series cannot be adjusted with X13.\nPlot the series -- most likely is it piecewise linear, although other irregularities are possible.")
      return(ret <- list(seasonality = FALSE, multiplicative = FALSE,
                         quality = list(log = NA, arima = "",
                                        Noutliers = NA, N.AO = NA, N.LS = NA, N.TC = NA,
                                        outliers = NULL, datesAO = NULL, datesLS = NULL, datesTC = NULL,
                                        M = rep(NA, 11), Q_M2 = NA, IC = rep(NA, 3)),
                         series = stats::ts.union(original = x, seasonal = 0, adjusted = x, predicted = x), seas = NULL, date = ts2date(x)))
    } else {
      trans.fun <- if (isTRUE(a.td$udg["aictrans"] != "No transformation") | force.multiplicative) "log" else "none"
    }
  }

  if (verbose > 0) cat("*** The adjustment model is ", if (trans.fun == "log") "MULTIPLICATIVE" else "ADDITIVE", ".\n", sep = "")
  if (verbose > 1) printSym("-", linelen)
  if (verbose == 1) cat("Estimating 3 main models (6 TD, 1 TD, 0 TD)...\n")
  # The rest will be similar to the first call, so creating a new function
  estSeas <- function(...) tryCatch(seasonal::seas(x = this.ser, x11 = "", outlier.types = c("ao", "ls", "tc"),
                                          transform.function = trans.fun,
                                          automdl.maxorder = c(3, 1), estimate.tol = 1e-7,
                                          x11.appendfcst = do.fcast, forecast.maxlead = nfcast,
                                          x11.appendbcst = do.bcast, forecast.maxback = nbcast,
                                          forecast.save = c("bct", "fct"),
                                          ...), error = efv)
  # 1 weekday dummy + Easter + leap year
  if (verbose > 1) printSum(a.td, "-- Model: 6 TD + LY + Easter")
  a.wd <- estSeas(regression.variables = wd.vars, regression.aictest = NULL,
                  xreg = wd.cal, regression.usertype = wd.types)
  if (verbose > 1) printSum(a.wd, "-- Model: 1 TD + LY + Easter")
  # No day effects, only Easter and leap year that are tested via AIC
  a.nod <- estSeas(regression.variables = nod.vars, regression.aictest = NULL)
  if (verbose > 1) printSum(a.nod, "-- Model: 0 TD + LY + Easter")

  extrOut <- function(x) tryCatch(tolower(seasonal::outlier(x)), error = function(e) {ret <- rep(NA, length(this.ser)); attributes(ret) <- attributes(this.ser); ret})
  outlier.ts <- data.frame(td = extrOut(a.td), wd = extrOut(a.wd), nod = extrOut(a.nod))
  outlier.ts.med <- lapply(1:nrow(outlier.ts), function(i) OTMed3(as.matrix(outlier.ts)[i, ])) # Creating only one outlier per period
  outlier.ts.med.conc <- unlist(lapply(outlier.ts.med, attr, "concordant"))
  outlier.ts.med <- unlist(outlier.ts.med)
  outl.time <- paste0(floor(stats::time(this.ser)), ".", stats::cycle(this.ser))
  outlier.tab <- paste0(outlier.ts.med, outl.time)
  outl.reg <- outlier.tab[!grepl("^NA", outlier.tab)]

  if (length(outl.reg) > 0 & any(!outlier.ts.med.conc)) { # There are some outlier regressors and at least one is not concordant across models
    if (verbose > 1) cat("Outlier regressors used in testing: ", paste0(outl.reg, collapse = ", "), "\n", sep = "")
    if (verbose > 1) cat("Estimating models with identical outlier regressors to make the TD regressors comparable by AICc...\n")

    a.td2 <- tryCatch(update.seas(a.td, regression.variables = unique(c(td.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
    if (!is.null(a.td2)) {
      if (verbose > 1) printSum(a.td2, "-- Refit: 6 TD + LY + Easter + common OL")
      a.wd2 <- tryCatch(update.seas(a.wd, regression.variables = unique(c(wd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
      if (!is.null(a.wd2)) {
        if (verbose > 1) printSum(a.wd2, "-- Refit: 1 TD + LY + Easter + common OL")
        a.nod2 <- tryCatch(update.seas(a.nod, regression.variables = unique(c(nod.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
        if (!is.null(a.nod2)) {
          if (verbose > 1) printSum(a.nod2, "-- Refit: 0 TD + LY + Easter + common OL")
        } else a.td2 <- a.wd2 <- a.nod2 <- NULL
      } else a.td2 <- a.wd2 <- a.nod2 <- NULL
    } else a.td2 <- a.wd2 <- a.nod2 <- NULL
  } else {
    if (verbose > 1) cat("All outlier variables are identical in 3 models, skipping re-estimation.\n", sep = "")
    a.td2 <- a.td
    a.wd2 <- a.wd
    a.nod2 <- a.nod
  }

  # If re-estimation is successful, then, choose the best model based on AIC among the comparables, and then, go back to the original
  if (!is.null(a.td2) & !is.null(a.wd2) & !is.null(a.nod2)) {
    if (length(outl.reg) > 0 & any(!outlier.ts.med.conc)) {
      if (verbose > 1) cat("Estimation with common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
    }
    aiccs <- unname(c(getStat(a.td2)$IC["AICc"], getStat(a.wd2)$IC["AICc"], getStat(a.nod2)$IC["AICc"]))
    # Adjusting the AICc based on the number of observations
    eff.len <- as.numeric(unname(c(a.td2[["udg"]]["nefobs"], a.wd2[["udg"]]["nefobs"], a.nod2[["udg"]]["nefobs"])))
    aiccs <- aiccs / eff.len * n
  } else { # Something went wrong with the same set of outliers---possibly too many---do the quick and dirty selection with the original models
    if (verbose > 0) cat("! Estimation with common outliers failed for at least 1 model.\n")
    if (min.outliers == 1) {
      if (verbose > 0) cat("! Re-trying with fewer outliers (existing in at least 2 out of 3 models).\n")
      outlier.tab <- table(c(a.td$model$regression$variables, a.wd$model$regression$variables, a.nod$model$regression$variables))
      outlier.tab <- outlier.tab[outlier.tab >= 2]
      outl.reg <- sort(setdiff(names(outlier.tab), c("easter[8]", "lpyear", "const", "tdnolpyear", "td1coef", "td1nolpyear")))
      if (verbose > 1) cat("Outlier regressors used in testing: ", paste0(outl.reg, collapse = ", "), "\n", sep = "")
      a.td2 <- tryCatch(update.seas(a.td, regression.variables = unique(c(td.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
      if (!is.null(a.td2)) {
        if (verbose > 1) printSum(a.td2, "-- Refit: 6 TD + LY + Easter + common OL")
        a.wd2 <- tryCatch(update.seas(a.wd, regression.variables = unique(c(wd.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
        if (!is.null(a.wd2)) {
          if (verbose > 1) printSum(a.wd2, "-- Refit: 1 TD + LY + Easter + common OL")
          a.nod2 <- tryCatch(update.seas(a.nod, regression.variables = unique(c(nod.vars, outl.reg)), outlier = NULL, transform.function = trans.fun), error = efv)
          if (!is.null(a.nod2)) {
            if (verbose > 1) printSum(a.nod2, "-- Refit: 0 TD + LY + Easter + common OL")
          } else a.td2 <- a.wd2 <- a.nod2 <- NULL
        } else a.td2 <- a.wd2 <- a.nod2 <- NULL
      } else a.td2 <- a.wd2 <- a.nod2 <- NULL
      if (!is.null(a.td2) & !is.null(a.wd2) & !is.null(a.nod2)) {
        if (verbose > 1) cat("Re-estimation with fewer common outliers succeeded, no convergence issues.\nChoosing the best model based on common outliers, using the original model with its specific outliers.\n")
        mod.list <- list(a.td2, a.wd2, a.nod2)
        aiccs   <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
        eff.len <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 0))
        aiccs <- aiccs / eff.len * n
      } else {
        if (verbose > 0) cat("! Re-estimation with fewer common outliers was unsuccessful.\n! Doing selection with the original models.\n! The results might be sub-optimal because the outliers are different.\n")
        mod.list <- list(a.td, a.wd, a.nod)
        aiccs <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
        eff.len <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 0))
        aiccs <- aiccs / eff.len * n
      }
    } else { # Min. outliers > 1
      if (verbose > 0) cat("! Re-estimation based on ", min.outliers, " common outliers (fixed number) was unsuccessful!\n! Doing selection with the original models.\n! The results might be sub-optimal because the outliers are different.\n", sep = "")
      mod.list <- list(a.td, a.wd, a.nod)
      aiccs   <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["aicc"]) else Inf))
      eff.len <- unname(sapply(mod.list, function(x) if (!is.null(x)) as.numeric(x[["udg"]]["nefobs"]) else 0))
      aiccs <- aiccs / eff.len * n
    }
  }
  names(aiccs) <- c("td", "wd", "nod")
  if (any(is.finite(aiccs))) { # At least one model was estimated in overall
    best.type <- "nod"
    if (aiccs["wd"] - aiccs["nod"] < tradingdays.aicdiff) best.type <- "wd"
    if ((aiccs["td"] - aiccs[best.type]) < tradingdays.aicdiff) best.type <- "td"
    if (verbose > 0) switch(best.type, td = cat("*** The week-day effects (6 TD) are different and important.\n"), wd = cat("*** The number of week days (1 TD) is important.\n"), nod = cat("*** The week-day effect is absent in these series (0 TD).\n"))
    best.mod <- switch(best.type, td = a.td, wd = a.wd, nod = a.nod)
  } else {
    warning("No model with the recommended settings could be estimated. Trying with the simplest default X13 settings.")
    best.type <- "nod"
    best.mod <- seasonal::seas(this.ser, x11 = "")
  }
  xregs <- best.mod$model$regression$variables
  if (verbose > 1) cat("Outlier regressors in the model with forced Easter and LY effects: ", paste0(sort(setdiff(xregs, c("easter[8]", "lpyear", "const", "td1nolpyear"))), collapse = ", "), "\n", sep = "")

  # Now testing for Easter and other effects because it was forced
  # Keeping the same outlier and regressors
  best.easter <- tryCatch(update.seas(best.mod, regression.variables = xregs, regression.aictest = c("easter", "lpyear"), transform.function = trans.fun, outlier = NULL), error = efv)
  if (!is.null(best.easter)) {
    if (!isTRUE(all.equal(sort(names(best.easter$est$coefficients)), sort(names(best.mod$est$coefficients))))) {
      # Easter/LY re-estimation succeeded and something changed
      if (verbose > 1) printSum(best.easter, "-- Refit: auto-(Easter | LY) + common OL")
      if (verbose > 1) cat("Easter and LY effect testing yielded a different model, re-estimating with auto-selection...\n", sep = "")
    } else {
      if (verbose > 1) cat("Easter and LY effect test changed nothing, re-estimating with auto-selection...\n")
    }
    # Then, we re-estimate the specification with full freedom about outliers with the set holiday effects (no forced outlier xregs)
    best.easter2 <- estSeas(regression.variables = switch(best.type, td = td.vars, wd = wd.vars, nod = nod.vars),
                            regression.aictest = c("easter", "lpyear"),
                            xreg = switch(best.type, td = td.cal, wd = wd.cal, nod = NULL),
                            regression.usertype = switch(best.type, td = td.types, wd = wd.types, nod = NULL)
    )
    if (!is.null(best.easter2)) {
      best.mod <- best.easter2
      if (verbose > 1) printSum(best.easter2, "-- Refit: auto-(Easter, LY, outl. regr.)")
    } else {
      best.mod <- best.easter
      if (verbose > 0) cat("! Re-estimation w/auto-Easter/LY effects failed, keeping the test-based model.\n")
    }
  } else {
    if (verbose > 0) cat("! Even testing w/auto-Easter/LY effects failed! Keeping the original model, but this might indicate that the series are problematic.\n", sep = "")
  }

  keep.easter <- "Easter[8]" %in% names(best.mod$est$coefficients)
  keep.ly <- "Leap Year" %in% names(best.mod$est$coefficients)
  if (verbose > 0) cat("*** Easter effect: ", if (!keep.easter) "NO" else "YES", ", leap-year effect: ", if (!keep.ly) "NO" else "YES", ".\n", sep = "")
  if (verbose > 0) cat("*** Outlier regressors: ", paste0(sort(setdiff(best.mod$model$regression$variables, c("easter[8]", "lpyear", "const"))), collapse = ", "), "\n", sep = "")
  if (verbose > 1) printSym("-", linelen)

  this.stat <- getStat(best.mod)

  # Plotting
  xs <- range(myYear(this.ser))
  if (head.years > 0) xs[1] <- xs[1] - head.years
  if (tail.years > 0) xs[2] <- xs[2] + tail.years
  xlim <- c(floor(xs[1]), ceiling(xs[2]))
  xs5 <- (ceiling(xlim[1]/5):floor(xlim[2]/5))*5 # Getting round year numbers divisible by 5
  cal.adj <- if (this.stat$log) best.mod$series$d16 / best.mod$series$d10 else best.mod$series$d16 - best.mod$series$d10 # Getting pure calendar adjustment
  mycols <- c("#CC003399", "#377EB8BB", "#27ba30", "#000000EE", "#FF4D00")
  vxs <- setdiff(seq(round(xlim[1]), round(xlim[2])+1), xs5) # Secondary axis lines

  # Outliers for plotting
  ol.ts <- seasonal::outlier(best.mod)
  sym.ts <- as.numeric(factor(ol.ts, levels = c("AO", "LS", "TC")))

  # Not plotting if plot.file is NA
  if (!isTRUE(is.na(plot.file))) {
    if (!is.null(plot.file)) grDevices::png(plot.file, 800, 900, type = "cairo", pointsize = 20)
    graphics::layout(matrix(c(1, 1, 2, 2, 3, 3, 4, 5), ncol = 2, nrow = 4, byrow = TRUE), heights = c(1.5, 1, 1, 1))
    # Plot 1: series and adjusted
    withr::local_par(mar = c(0.2, 2, 0.2, 1.2))
    ylim <- range(best.mod$series$d11, best.mod$series$d12, this.ser, na.rm = TRUE)
    ylim <- ylim + c(-0.01, 0.01)*diff(ylim) # Extending slightly
    plot(this.ser, xlim = xlim, ylim = ylim, bty = "n", col = mycols[3], lwd = 2, xaxt = "n")
    graphics::lines(best.mod$series$d11, col = "#FFFFFFAA", lwd = 4.5) # Adjusted with a halo
    graphics::lines(best.mod$series$d11, col = mycols[4], lwd = 3) # Adjusted
    graphics::lines(best.mod$series$d12, col = "#FFFFFFAA", lwd = 3.5) # Trend halo
    graphics::lines(best.mod$series$d12, col = mycols[5], lwd = 2, lty = 2) # Trend
    graphics::abline(v = vxs - 0.5/freq, col = "#0000003A", lty = 2)
    graphics::abline(v = xs5 - 0.5/freq, col = "#0000003A", lty = 1, lwd = 1.5)
    graphics::legend("bottomright", c("Original", "Adjusted", "Trend"), col = mycols[3:5], lwd = c(2, 3, 2), lty = c(1, 1, 2), bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    graphics::legend("top", paste0("Seasonality adjustment for ", name), bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    graphics::points(this.ser, pch = as.numeric(sym.ts), col = "#FFFFFFEE", lwd = 4)
    graphics::points(this.ser, pch = as.numeric(sym.ts), col = "#FF0000", lwd = 2)
    # Labels at the proper side
    is.below.trend <- (best.mod$series$d11 < best.mod$series$d12) & is.finite(best.mod$series$d11)
    if (any(!is.na(ol.ts))) labelsWithHalo(stats::time(this.ser)[!is.na(ol.ts)], y = this.ser[!is.na(ol.ts)], labels = ol.ts[!is.na(ol.ts)], pos = ifelse(is.below.trend[!is.na(ol.ts)], 3, 1), cex = 0.75, offset = 0.4, hscale = 0.003, vscale = 0.01, nhalo = 16)
    # Plot 2: only seasonality
    withr::local_par(mar = c(0.1, 2, 2, 1.2))
    ylim <- range(best.mod$series$d10, cal.adj, na.rm = TRUE)
    plot(NULL, NULL, bty = "n", xlim = xlim, ylim = ylim, xaxt = "n")
    graphics::abline(h = if (isTRUE(this.stat$log)) 1 else 0, lty = 2)
    graphics::lines(best.mod$series$d10, lwd = 3, col = mycols[1]) # Only seasonal
    do.cal <- !isTRUE(all.equal(diff(as.numeric(cal.adj)), rep(0, length(cal.adj)-1))) # Is there any deterministic calendar effect? No calendar --> d10 = d16
    if (do.cal) {
      graphics::lines(cal.adj, col = "#FFFFFFBB", lwd = 4) # Calendar halo
      graphics::lines(cal.adj, col = mycols[2],   lwd = 2)
    }
    graphics::abline(v = vxs - 0.5/freq, col = "#0000003A", lty = 2)
    graphics::abline(v = xs5 - 0.5/freq, col = "#0000003A", lty = 1, lwd = 1.5)
    graphics::legend("topright", if (isTRUE(this.stat$log)) "Multiplicative" else "Additive", bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    leg.text <- c("Seasonal", "Calendar")
    if (do.cal) graphics::legend("bottomright", leg.text, ncol = 2, col = mycols[1:2], lwd = c(3, 2), lty = 1, bg = "#FFFFFFBB", box.col = "#FFFFFFBB") else
      graphics::legend("bottomright", leg.text[1], col = mycols[1], lwd = 3, bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    graphics::legend("bottomleft", c(paste0("M7=", sprintf("%1.3f", this.stat$M[7])), paste0("Q2=", sprintf("%1.3f", this.stat$Q_M2))), bg = "#FFFFFFBB", box.col = "#FFFFFFBB")
    axs <- xs5
    if (xlim[1] - min(xs5) > 1) axs <- c(min(axs)-5, axs) # Extending short axes
    if (xlim[2] - max(xs5) > 1) axs <- c(axs, max(axs)+5)
    graphics::axis(3, at = axs - 0.5/freq, labels = axs)
    # Plot 3: seasonal components
    withr::local_par(mar = c(2, 2, 0.2, 0.2))
    if (!seasonal.plot) {
      suppressWarnings(graphics::boxplot(scale(this.ser) ~ stats::cycle(this.ser), frame = FALSE, notch = TRUE, xlab = "", ylab = "", yaxt = "n"))
      graphics::abline(h = stats::median(scale(this.ser), na.rm = TRUE), lty = 2)
    } else {
      stats::monthplot(this.ser, box = FALSE)
    }
    # Plot 3: spectra before and after adjustment
    withr::local_par(mar = c(2, 0.2, 0.2, 0.2))
    computeSpectrum(this.ser, plot = TRUE, compact = TRUE)
    xspec <- computeSpectrum(best.mod$series$d11, plot = TRUE, compact = TRUE, put.legend = FALSE)
    graphics::legend(attr(xspec, "legend.position"), c("Left: original", "Right: adjusted"), bg = "#FFFFFFBB", box.col = "#FFFFFFBB")

    if (!is.null(plot.file)) grDevices::dev.off()
  }

  if (verbose > 1) {
    cat("Seasonality diagnostics for ", name, ":\n", sep = "")
    print(this.stat$M)
    cat("Overall quality, Q_M2: ", this.stat$Q_M2, " (<1 is OK, >1 is not OK).\n", sep = "")
  } else {
    if (verbose > 0) cat("Seasonality diagnostics for ", name, ": M7 = ", this.stat$M[7], ", Q_M2 = ", this.stat$Q_M2, " (<1 is OK, >1 is not OK).\n", sep = "")
  }


  # In case there were missing values, put the regression-based prediction
  x.pred <- this.ser
  if (length(miss.inds) > 0) x.pred[miss.inds] <- predict.seas(best.mod)[miss.inds]
  # Creating the forecast series in any case
  x.bc <- if (nbcast > 0) seasonal::series(best.mod, "bct") else NULL
  x.fc <- if (nfcast > 0) seasonal::series(best.mod, "fct") else NULL
  # Prepending backcasts, appending forecasts
  x.bf <- if (nbcast > 0) stats::ts(c(x.bc[, "backcast"], x.pred), end = stats::end(x.pred), freq = freq) else x.pred
  x.bf <- if (nfcast > 0) stats::ts(c(x.bf, x.fc[, "forecast"]), start = stats::start(x.bf), freq = freq) else x.bf

  if (this.stat$M[7] <= m7.threshold) { # There is seasonality
    if (this.stat$Q_M2 <= q.threshold) {
      if (verbose > 0) cat("There seem to be no immediate seasonality identifiability and adjustability problems in ", name, ".\n", sep = "")
    } else {
      if (verbose > 0) cat("The series ", name, " have some seasonality, but the adjustment quality is mediocre (or the trend is flat, which is harmless). Use with caution!\n", sep = "")
    }

    x.sa <- best.mod$series$d11
    x.seas <- best.mod$series$d16 # Not D10 because there might be calendar effects as well
    out.ts <- stats::ts.union(original = x, seasonal = x.seas, adjusted = x.sa, predicted = x.bf, trend = best.mod$series$d12, irregular = best.mod$series$d13)
    out.date <- ts2date(out.ts)

    ret <- list(seasonality = TRUE, multiplicative = this.stat$log, quality = this.stat, series = out.ts, seas = best.mod, date = out.date)
  } else {
    if (this.stat$Q_M2 > q.threshold) {
      if (verbose > 0) cat("The series ", name, " do not have stable seasonality and are very noisy, skipping.\n", sep = "")
    } else {
      if (verbose > 0) cat("The series ", name, " do not have stable seasonality and seem regular, skipping.\n", sep = "")
    }

    x.sa <- x # No adjustment is done
    x.seas <- if (this.stat$log) rep(1, length(x)) else rep(0, length(x))
    out.ts <- stats::ts.union(original = x, seasonal = x.seas, adjusted = x.sa, predicted = x.bf, trend = best.mod$series$d12, irregular = x - best.mod$series$d12)
    # Restoring the middle missing values in the adjusted series
    mg <- findMiddleGaps(out.ts[, "adjusted"])
    if (length(mg) > 0) out.ts[mg, "adjusted"] <- out.ts[mg, "predicted"]
    out.date <- ts2date(out.ts)
    ret <- list(seasonality = FALSE, multiplicative = NA, quality = this.stat, series = out.ts, seas = best.mod, date = out.date)
  }

  # Penultimate: some simple seasonality tests
  # Not even filtering out the trend
  # The Kruskal--Wallis and Friendman are non-robust
  # We run two simple linear regressions and test joint significance
  # 1: Regress the (S+I) on seasonal dummies (robust to outliers + robust ANOVA)
  # 2: Regress the (S+I) onto the years and periods; test if something changes in years
  si <- if (this.stat$log) best.mod$series$d10 * best.mod$series$d13 else best.mod$series$d10 + best.mod$series$d13
  vHAC <- function(x) tryCatch(sandwich::kernHAC(x, prewhite = 0, kernel = "Bartlett", bw = sandwich::bwNeweyWest), error = function(e) sandwich::vcovHC(x, type = "HC0"))
  vHC <- function(x) tryCatch(sandwich::vcovHC(x, type = "HC0"))
  lm.seas <- tryCatch(MASS::rlm(si ~ factor(stats::cycle(si)), maxit = 100), error = efv)
  if (is.null(lm.seas)) lm.seas <- stats::lm(rank(si) ~ factor(stats::cycle(si)))

  test.seas <- tryCatch(car::linearHypothesis(lm.seas, names(stats::coef(lm.seas))[-1], vcov. = vHAC), error = efv)
  if (is.null(test.seas)) test.seas <- tryCatch(car::linearHypothesis(lm.seas, names(stats::coef(lm.seas))[-1], vcov. = vHC), error = efv)
  if (is.null(test.seas)) test.seas <- tryCatch(car::linearHypothesis(lm.seas, names(stats::coef(lm.seas))[-1]), error = function(e) return("Hypothesis testing failed."))
  attr(test.seas, "heading") <- c("Seasonalily model: R_ij = a_i + U_ij", "(HAC-robust ANOVA for seasonal dummies in a robust regression)", attr(test.seas, "heading"))
  si.year <- myYear(si)
  lm.stab.seas <- tryCatch(MASS::rlm(si ~ factor(stats::cycle(si)) + factor(si.year), maxit = 100), error = efv)
  if (is.null(lm.stab.seas)) lm.stab.seas <- stats::lm(si ~ factor(stats::cycle(si)) + factor(si.year))
  hyp.stab <- names(stats::coef(lm.stab.seas))[grep("factor\\(si.year\\)", names(stats::coef(lm.stab.seas)))]
  test.stab.seas <- tryCatch(car::linearHypothesis(lm.stab.seas, hyp.stab, vcov. = vHAC), error = efv)
  if (is.null(test.stab.seas)) test.stab.seas <- tryCatch(car::linearHypothesis(lm.stab.seas, hyp.stab, vcov. = vHC), error = efv)
  if (is.null(test.stab.seas)) test.stab.seas <- tryCatch(car::linearHypothesis(lm.stab.seas, hyp.stab), error = function(e) return("Hypothesis testing failed."))
  attr(test.stab.seas, "heading") <- c("Evolutive seasonalily model: Y_ij = a_i + b_j + U_ij", "(HAC-robust ANOVA in the spirit of Higgins (1975))", attr(test.stab.seas, "heading"))
  # For the Friedman test, we take full years
  # friedman.mat <- matrix(as.numeric(tail(this.ser, floor(length(this.ser) / freq)*freq)), ncol = freq, byrow = TRUE)
  # ft <- friedman.test(friedman.mat)
  robust.m7 <- sqrt((1.5*test.stab.seas$F[2] + 3.5) / test.seas$F[2])
  ret$si.seasonality.tests <- list(identifiable = test.seas, stable = test.stab.seas, robust.m7 = robust.m7)

  # Last check: positivity
  if (all(stats::na.omit(as.numeric(x)) > 0)) { # The original data are positive
    if (verbose > 1) printSym("-", linelen)
    if (verbose > 1) cat("The input series contain only positive values in the range", range(x, na.rm = TRUE), "\n")
    if (all(stats::na.omit(as.numeric(x.sa)) > 0)) {
      if (verbose > 1) cat("The adjusted series contain values in the range", range(x.sa, na.rm = TRUE), "\n")
    } else {
      if (verbose > 0) cat("! The adjusted series contain some negative values, range", range(x.sa, na.rm = TRUE), "\n! Enforcing a log transformation may be necessary.\n")
    }
    if (trans.fun != "log" & verbose > 1) cat("! All the values are positive, yet the model is additive.\n! Consider 'force.multiplicative = TRUE' and see if a multiplicative model is better.\n")
  }
  if (verbose > 1) printSym("=", linelen)
  if (verbose > 0) cat("\n")

  return(ret)
}

#' Read national calendars into the global environment.
#'
#' @param path Character: full directory path to use (ending with a slash).
#' @param country Character: file name (without the "-Q.csv") prefix
#'
#' @return Nothing; loads two MTS into the global environment.
#' @export
#'
#' @examples
#' tryCatch(getCalendars("France"), error = function(e) print("Calendar not found."))
getCalendars <- function(country = "Luxembourg", path = NULL) {
  if (is.null(path)) path <- "S:/Projets/Modelisation/Mod\u00E8les/Seasonal Adjustment/calendars/"
  if (!(substr(path, nchar(path), nchar(path)) %in% c("/", "\\"))) stop("The 'path' argument should end with / or \\.")
  cal.q <- utils::read.csv(paste0(path, country, "-Q.csv"))
  cal.m <- utils::read.csv(paste0(path, country, "-M.csv"))
  cal.q$Date <- as.Date(cal.q$Date, format = "%d/%m/%Y")
  cal.m$Date <- as.Date(cal.m$Date, format = "%d/%m/%Y")
  cal.q <- suppressMessages(makeTS(cal.q))
  cal.m <- suppressMessages(makeTS(cal.m))
  assign("cal.q", cal.q, envir = .GlobalEnv)
  assign("cal.m", cal.m, envir = .GlobalEnv)
  print("Loaded the quartely and monthly calendars into the global environment as cal.q and cal.m.")
  return(invisible(NULL))
}

#' Reconstruct missing values in time series via multiple methods
#'
#' @param x A ts object (not mts)
#' @param robust.weights # Weights assigned to the 5 or 7 methods.
# If given a numeric vector of length 5 (7 for imputePanel7), return the weighted mean of sorted forecasts (inconsistent weighting may occur).
# c(1, 3, 5, 3, 1) is robust; c(0, 1, 2, 1, 0) is very robust; c(0, 0, 1, 0, 0) is the median.
#' @param calendar Passed to diagnoseSeasonality().
#' @param sample.begin Passed to diagnoseSeasonality() (shrink the estimation sample).
#' @param sample.end Passed to diagnoseSeasonality() (shrink the estimation sample).
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
      fns <- as.character(utils::lsf.str(envir = .GlobalEnv))
      cl <- parallel::makeCluster(cores)
      parallel::clusterExport(cl, fns) # Exporting global functions
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
  if (all(stats::na.omit(as.numeric(y.series)) > 0) | force.positive) {
    y.seas <- do.call(stats::ts.union, lapply(seas.mods, function(x) x$series[, "seasonal"]))
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
        if (isTRUE(y.mult[i]))  imp.pca[, i] <- imp.pca[, i] * y.seas[, i]
        if (isFALSE(y.mult[i])) imp.pca[, i] <- imp.pca[, i] + y.seas[, i]
        # If it is NA, there is no seasonality and no need to change anything
      } else { # There is nothing to impute -- restore the original unadjusted
        imp.pca[, i] <- y.series[, i]
      }
    }
  } else { # Not applying the correction
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


