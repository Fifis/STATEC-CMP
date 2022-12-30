#' Dynamic multipliers for an ARDL model
#'
#' @param endog.lag A numeric vector of AR coefficients on Y[t-1], Y[t-1] etc.
#' @param exog.lag A list (!) of vectors of DL coefficients on X[t], X[t-1] etc.
#' @param h An integer scalar: horizon at which to terminate the iterative procedure
#' @param shock.amount A scalar or numeric vector of the same length as the number of exogenous variables: change of each variable to which Y responds; 1 yields pure dynamic multipliers.
#' @param xnames A character vector of human-readable expanded exogenous variable names.
#' @param resid.name The name for the last term of the RHS (the residual with a unit coefficient on it)
#' @param plot Logical: if TRUE, produces a plot to the current device.
#' @param yname Dependent variable names for the plot title.
#' @param samescale Logical: use the same vertical axis range to have comparable responses?
#' @param plot.horizon An integer: how many periods of the response to plot.
#' @param format.fun A function that converts a real humber into a nicely rounded number with the desired precision.
#' @param add.legend Logical: if TRUE, add an extra panel with a legend.
#' @param francais Logical: if TRUE, the labels are in French.
#' @param ... Passed to all plot functions tha produce the plots.
#'
#' @return A numeric matrix containing the dynamic multipliers for the requested number of periods (one period per line), starting with period 0 in the 1st line.
#' @export
#' @seealso [computePropag()] uses dynamic multipliers to decompose the impact of shocks
#' in a real data set.
#' [ECM2ARDL()] converts an ECM into an equivalent ARDL.
#'
#' @examples
#' rho <- 0.2
#' dx <- c(x1 = 0.9, x2 = 0.2)
#' ECT <- -0.8
#' lx <- c(x1 = -1.5, x2 = 0.5)
#' ARDL.coef <- ECM2ARDL(dx.coef.list = dx, dy.coef.vec = rho, ECT = ECT, long.run = lx)
#' dm <- genDynMult(endog.lag = ARDL.coef$AR, exog.lag = ARDL.coef$DL, add.legend = TRUE)
#' head(dm, 11)
genDynMult <- function(endog.lag, exog.lag, h = 50, shock.amount = 1,
                       xnames = NULL, resid.name = "rest", plot = TRUE,
                       yname = NULL, samescale = TRUE, plot.horizon = 12,
                       format.fun = function(x) sprintf("%1.3f", x),
                       add.legend = FALSE, francais = FALSE, ...
                       ) {
  # The presence of a constant term does not change the responses of other variables
  if (is.vector(exog.lag, mode = "numeric")) exog.lag <- list(exog.lag)
  k <- length(exog.lag[[1]]) # How many exogenous regressors
  vn <- names(exog.lag[[1]])
  p <- length(endog.lag) # AR order
  q <- length(exog.lag) - 1 # Distributed lag order: the contemporaneous value is always there, i.e. q shows how many extras are there
  if (q > 0) {
    if (!all(diff(sapply(exog.lag, length)) == 0)) stop("The list of coefficients on exogenous lags must contain elements of equal lengths.")
  }
  if (is.null(xnames)) xnames <- vn
  if (is.null(xnames)) xnames <- paste0("V", 1:k)
  if (length(xnames) != k) stop(paste0("The list of variable names must have the same length (", length(xnames), ") as the number of coefficients (", k, ")."))
  for (i in 1:(q+1)) names(exog.lag[[i]]) <- xnames

  if (!francais) {
    cat("Dynamic multiplier computation for ARDL(", p, ",", q, "):\n" , sep = "")
    cat("Y = ", paste0("c", 1:p, "*Y(-", 1:p, ")", collapse = " + "), " + ", gsub("\\(-0\\)", "", paste0("b", 0:q, "'X(-", 0:q, ")", collapse = " + ")), " + U\n" , sep = "")
    cat("Output contains the final form coefficients:\nY[t] = sum[i=0..Inf] a[i]'X[t-i] + b[i]*U[t-i]\n")
  } else {
    cat("Calcul des multiplicateurs dynamiques ARDL(", p, ",", q, "):\n" , sep = "")
    cat("Y = ", paste0("c", 1:p, "*Y(-", 1:p, ")", collapse = " + "), " + ", gsub("\\(-0\\)", "", paste0("b", 0:q, "'X(-", 0:q, ")", collapse = " + ")), " + U\n" , sep = "")
    cat("Sortie = coefficients de la forme finale :\nY[t] = sum[i=0..Inf] a[i]'X[t-i] + b[i]*U[t-i]\n")
  }

  max.order <- max(p, q)
  hh <- h + max.order + 1 # For trimming afterwards
  # Data frame structure: dep. var., X, error term
  d <- matrix(0, nrow = hh, ncol = length(exog.lag[[1]]) + 2)
  rownames(d) <- paste0("t+", rev((1:hh)-1))
  # d[, 1] is the dependent variable, d[, ncol(d)] is the error term
  # d[, -c(y, e)] is all the middle columns corresponding to regressors
  # Indices instead of names are used because many data sets have variables names like "Y", "U", "Error" etc.
  y <- 1
  e <- ncol(d)
  dl.matrix <- do.call(rbind, exog.lag)
  d[hh, e] <- 1 # Starting at the end (time 't')
  d[hh - 1:p, y]        <- endog.lag
  d[hh - 0:q, -c(y, e)] <- dl.matrix
  # At each step, substituting only one term
  # y[t] = rho1*y[t-1] + {...} + U[t]
  # y[t-1] = rho1*y[t-2] + [...] + U[t-1]
  # Therefore, y[t] = rho1*(rho1*y[t-2] + [...] + U[t-1]) + {...} + U[t]
  for (i in 1:h) {
    rho <- d[hh-i, y]
    d[hh-i, e] <- rho
    d[hh-i-(1:p), y] <- d[hh-i-(1:p), y] + endog.lag*rho # Coefficients on the further lags of Y
    d[hh-i-(0:q), -c(y, e)] <- d[hh-i-(0:q), -c(y, e)] + dl.matrix*rho
    d[hh-i, y] <- 0
  }
  pers <- sum(abs(d[1:(max.order+1), y]))
  if (pers > sqrt(.Machine$double.eps)) warning(paste0("There is too much remaining persistence (", pers, ").\nEither increase the horizon (h) or check the model (it seems to be unstable)."))
  if (!francais) {
    cat("=====\nPersistence remaining beyond lag ", h, ": ", signif(pers, 2), ".\n", sep = "")
  } else cat("=====\nPersistance restante au-del\u00E0 du retard ", h, ": ", signif(pers, 2), ".\n", sep = "")
  colnames(d) <- c("DroppedY", xnames, resid.name)
  attr(d, "endog.lag") <- endog.lag # Extra output so that we know which coefficients yielded this
  attr(d, "exog.lag") <- exog.lag
  attr(d, "(Intercept)") <- sum(d[, e])
  d <- d[nrow(d):(max.order+1), -y]

  if (length(shock.amount) == 1) shock.amount <- rep(shock.amount, k+1)
  if (!all(shock.amount == 1)) {
    if (length(shock.amount) != (k+1)) stop("The length of the shock.amount vector must equal the # of exogenous regressors + 1 (for the error)")
    d <- sweep(d, 2, shock.amount, "*")
  }

  if (plot) {
    if (is.null(yname)) yname <- "Y"
    if (plot.horizon < nrow(d)) plot.horizon <- plot.horizon+1
    colnames(d)[ncol(d)] <- resid.name
    cumd <- apply(d, 2, cumsum)
    withr::local_par(mfrow = grDevices::n2mfrow(ncol(d) + as.numeric(add.legend)), mar = c(2, 2, 2.1, 0.2))
    t <- (1:plot.horizon) - 1
    for (i in 1:ncol(d)) {
      # Change to (0.001, 0.999) to get rid of outliers
      if (samescale) yl <- stats::quantile(c(unlist(d), unlist(cumd)), c(0, 1)) else yl <- stats::quantile(c(d[, i], cumd[, i]), c(0, 1))
      plot(t[1:plot.horizon], d[1:plot.horizon, i], bty = "n", type = "b", lty = 2, ylab = "", xlab = "", main = paste0("\u0394", yname, if (francais) " d\u00FB \u00E0 \u0394" else " due to \u0394", colnames(d)[i], "=", format.fun(shock.amount[i])), ylim = yl, ...)
      graphics::abline(h = c(0, cumd[nrow(cumd), i]), lty = 3)
      graphics::lines(t[1:plot.horizon], cumd[1:plot.horizon, i], type = "b", pch = 16)
      graphics::abline(v = seq(0, 20, 5), col = "#00000011")
      obstructed <- as.character(cut(c(cumd[plot.horizon, i], 0), c(-Inf, yl[1]*0.7 + yl[2]*0.3, yl[1]*0.3 + yl[2]*0.7, Inf), labels = c("bottomright", "right", "topright"))) # Are we obstricting the view
      unobstructed <- setdiff(c("right", "topright", "bottomright"), obstructed)[1]
      graphics::legend(x = unobstructed, c(paste0(if (!francais) c("Short run = ", "Up to t+2 = ", "Long run = ") else c("Court terme = ", "Apr\u00E8s 3 ans = ", "Long terme = "), format.fun(cumd[c(1, 3, nrow(cumd)), i]))), bty = "n", bg = "#FFFFFF88")
    }
    if (add.legend) {
      plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "", bty = "n")
      graphics::legend("topleft", if (!francais) c("Dynamic multiplier", "Cumulative dynamic mult.") else c("Multiplicateur dynamique", "Mult. dynamique cumul\u00E9"), bty = "n", lty = c(3, 1), pch = c(1, 16))
      ar.round <- paste0("AR coef.: ", paste0(format.fun(endog.lag), collapse = ", "))
      dl.round <- paste0("X[t-", 0:q, "]: ", sapply(exog.lag, function(x) paste0(format.fun(x), collapse = ", ")))
      graphics::legend("bottomright", c(ar.round, dl.round), bty = "n")
    }
  }

  return(d)
}

#' Compute shock propagation and contributions with real data
#' @rdname propag
#' @param endog.lag Passed to [genDynMult()].
#' @param exog.lag Passed to [genDynMult()].
#' @param data A data frame or matrix with with the regressors (in the order in which they appear in exog.lag) and residuals in the last column.
#' @param xnames Passed to [genDynMult()].
#' @param mult A numeric scalar or vector: multiplier(s) for the units of change; 100 yields percentages. Another useful option is apply(data, 2, sd).
#' @param resid.name A label how to denote the residuals.
#' @param francais Logical: if TRUE, the column names will be in French.
#'
#' @return A numeric matrix containing the contribution by period type.
#' @export
#' @seealso [genDynMult()] is the workhorse used to compute the dynamic multipliers
#' based on the coefficients.
#'
#' @examples
#' # Imagine an error-correction model:
#' # dY = 0.2*dY(-1) + 0.9*dX1 + 0.2*dX2 - 0.8(Y(-1) - 1.5*X1(-1) + 0.5*X2(-1)) + U
#' # Its ARDL equivalent is
#' # Y = 0.4*Y(-1) - 0.2*Y(-2) + 0.9*X1 + 0.2*X2 + 0.3*X1(-1) - 0.6*X2(-1) + U
#' b <- list(AR = c(0.4, -0.2), DL = list(L0 = c(0.9, 0.2), L1 = c(0.3, 0.6)))
#' d <- simARDL(endog.lag = b$AR, exog.lag = b$DL, n = 40)
#' prpg <- computePropag(endog.lag = b$AR, exog.lag = b$DL, data = d[, -1], mult = 1)
#' prpg[[40]]
#' contrib.matrices <- quiet(printContrib(prpg, sums = FALSE, type = "twocolumn"))
#' contrib.total <- t(contrib.matrices$ST + contrib.matrices$LT)
#' tail(contrib.total)
computePropag <- function(endog.lag, exog.lag, data,
                          xnames = NULL, mult = 100,  resid.name = "reste", francais = TRUE
) {
  if (!any(class(data) %in% c("data.frame", "matrix"))) stop("'data' must be a data frame or a matrix.")
  n <- nrow(data)
  h <- n + 1 # "The dynamic multiplier horizon should be long enough to harness the full information
  dx <- rbind(0, as.matrix(data[-1, ] - data[-nrow(data), ]))

  dm <- quiet(genDynMult(endog.lag = endog.lag, exog.lag = exog.lag,
                         h = h, xnames = xnames, plot = FALSE, resid.name = resid.name))
  dm.all <- decomp <- vector("list", n)
  dm.all[2:n] <- lapply(2:n, function(i) sweep(dm, 2, dx[i, ]*mult, "*"))
  # Decomposing the changes at each period into the sum of components
  decomp[2:n] <- lapply(2:n, function(i) {
    increment.decomp <- dm[1:i, ] * dx[i:1, ] * mult # Decomposition of the present change via all dynamic multipliers
    contrib.present  <- increment.decomp[1, ] # Short-run effect, lag 0
    contrib.past     <- colSums(increment.decomp[-1, , drop = FALSE]) # Long-run contribution from all the past
    inertia.present  <- colSums(dm.all[[i]][-1, ]) # Remaining inertia from the present into the future
    # Remaining inertia from the past into the future
    # E.g. period i=5 has future inertia (t=6, ...) from the shocks in periods 2 {skipping mult. t+1..t+3}, 3 {skipping mult. t+1, t+2}, 4 {skipping mult. t+1}
    if (i > 2) {
      inertia.remain <- do.call(rbind, lapply(2:(i-1), function(j) colSums(dm.all[[j]][-(1:(i-j+1)), , drop = FALSE])))
      inertia.remain <- rbind(0, inertia.remain, matrix(0, ncol = ncol(dm), nrow = n-i))
    } else {
      inertia.remain <- matrix(0, nrow = 1, ncol = ncol(dx))
    }

    inertia.past <- colSums(inertia.remain)

    ret <- rbind(`SR part in present` = contrib.present, `LR part in present` = contrib.past, `Total change in present` = contrib.present + contrib.past,
                 `Remaining propag. of present` = inertia.present, `Remaining propag. of past` = inertia.past, `Remaining propag. total` = inertia.present + inertia.past)
    if (francais) rownames(ret) <- c("Partie CT dans pr\u00E9sent", "Partie LT dans pr\u00E9sent", "Changement pr\u00E9sent total", "Propag. restante du pr\u00E9sent", "Propag. restante du pass\u00E9", "Propag. restante totale")
    return(ret)
  })
  decomp[[1]] <- decomp[[2]]
  decomp[[1]][1:nrow(decomp[[1]]), 1:ncol(decomp[[1]])] <- NA

  return(decomp)
}

#' Simulate an ARDL process
#'
#' @param endog.lag Numeric: coefficients on the lags of Y
#' @param exog.lag A list of coefficients on X and its lags
#' @param intercept Numeric: shift for the entire series
#' @param seed Integer: seed for PRNG
#' @param n Number of observations
#'
#' @return A data frame with the dependent variable, regressors, and errors
#' @export
#'
#' @examples
#' b <- list(AR = c(0.4, -0.2), DL = list(L0 = c(0.9, 0.2), L1 = c(0.3, 0.6)))
#' d <- simARDL(endog.lag = b$AR, exog.lag = b$DL, n = 40)
#' d <- ts(d, start = c(2000, 1), freq = 4)
#' plot(d)
simARDL <- function(endog.lag, exog.lag, intercept = 20, seed = 1, n = 300) {
  set.seed(seed)
  U <- stats::rnorm(n)
  p <- length(endog.lag)
  k <- length(exog.lag[[1]])
  q <- length(exog.lag) - 1
  X <- MASS::mvrnorm(n, mu = rep(1, k), Sigma = 0.7*diag(k) + matrix(0.3, k, k))
  X <- X + (1:n)*0.01
  X <- sweep(X, 2, rep(c(1, -1), k)[1:k], "*") # To have some variability
  colnames(X) <- paste0("X", 1:k)
  Y <- numeric(n)
  b <- do.call(rbind, exog.lag)
  Y[1:p] <- sum(t(X[1:p, ]) * exog.lag[[1]]) + U[1:p]
  for (i in (p+1):n) Y[i] <- sum(Y[i - 1:p] * endog.lag) + sum(X[i - 0:q, ] * b) + U[i]
  Y <- Y + intercept
  return(data.frame(Y = Y, X, U = U))
}


#' @rdname propag
#' @param x A list of matrices each of which was returned by computePropag
#' @param time.labels A character vector of time period names (e.g. years)
#' @param var.labels A character vector of variable names
#' @param digits Integer: how many decimal places to keep (more than 1 could be false precision)
#' @param zero.replace Logical: print a dash instead of a 0 if the variable is exactly equal to 0 (to distinguish: 0.0 could be 0.004, whilst '---' is really less than 0.00000001)?
#' @param sums Logical: append row sums at the end?
#' @param rownames Logical: prepend variable names at the beginning?
#' @param excel Logical: use tabs instead of '&' for copy-pasting into Excel? If FALSE, adds TeX symbols
#' @param groups A named list containing column names of `x` to be treated added up to and treated as one variable, e.g. `groups = list(time_dum = c("c1", "c2"), unemp = c("r_ubit", "nawru"))`
#' @param type A string: `"byyear"` to output tables by year, `"bytype"` by contemporaneous remaining effect, `"twocolumn"` to output two columns (short and long term)
#'
#' @return Prints the output to console. Returns an invisible list of short-term and long-term contributions.
#' @export
printContrib <- function(x,
                  time.labels = NULL, var.labels = NULL,
                  digits = 1, zero.replace = TRUE,
                  sums = TRUE, rownames = TRUE,
                  excel = TRUE, francais = TRUE,
                  groups = NULL,
                  type = c("byyear", "bytype", "twocolumn")
                  ) {
  type <- type[1]
	sep <- if (excel) "\t" else " & "
	if (!is.list(x)) x <- list(x)
	l <- length(x)
	if (is.null(time.labels)) time.labels <- seq_len(l)
	for (i in 1:l) if (is.null(dim(x[[i]]))) x[[i]] <- matrix(x[[i]], nrow = 1)
	orig.names <- colnames(x[[1]])
	xn <- if (is.null(var.labels)) orig.names else var.labels
  # Combining variables by groups
  if (!is.null(groups)) {
    gl <- sapply(groups, length)
    gn <- names(groups)
    # if (any(gl < 2)) warning("Some groups contain fewer than two names. Renaming.")
    if (is.null(names(groups))) stop("The groups must be a named list, otherwise the reader will get lost.")
    bad.vars <- setdiff(unlist(groups), orig.names)
    if (length(bad.vars) > 0) stop(paste0("Some grouped variables are not found in the original data: ", paste0(bad.vars, collapse = ", ")))
    first.indices <- drop.indices <- NULL
    for (i in 1:length(groups)) {
      # Replacing the first variable with the group aggregate, and dropping the rest
      wh <- match(groups[[i]], orig.names)
      for (j in 1:l) x[[j]][, wh[1]] <- unname(rowSums(x[[j]][, wh, drop = FALSE]))
      first.indices <- c(first.indices, wh[1])
      gl <- length(groups[[i]])
      if (gl > 1) drop.indices <- c(drop.indices, wh[-1]) # No need to drop groups of length 1
    }
    for (j in 1:l) {
      colnames(x[[j]])[first.indices] <- gn
      x[[j]] <- x[[j]][, -drop.indices, drop = FALSE]
    }
    xn[first.indices] <- gn
    xn <- xn[-drop.indices]
  }
  words <- if (!francais) paste0(c("Short", "Long", "Total", "Pres.", "Past", "Total"), sep) else paste0(c("Court terme", "Long terme", "Totale", "du pr\u00E9sent", "du pass\u00E9", "Totale"), sep)

  if (sums) xn <- c(xn, if (francais) "Totale" else "Total")

  myFormat <- function(x, zero.replace, digits, excel, sums) {
    if (sums) x <- c(x, sum(x))
    which.zero <- which(x == 0)
    x <- sprintf(paste0("%1.", digits, "f"), x)
    if (!excel) x <- gsub("-", "$-$", x)
    if (zero.replace & length(which.zero) > 0) x[which.zero] <- if (!excel) "---" else "\u2014"
    return(x)
  }

  if (type == "byyear") { # Print the classical tables by year
    for (j in 1:l) {
      cat(if (rownames) paste0(time.labels[j], sep) else NULL, paste0(xn, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
      for (i in 1:nrow(x[[j]])) {
        this.x <- x[[j]][i, ]
        x.good <- myFormat(this.x, zero.replace, digits, excel, sums)
        cat(if (rownames) words[i] else NULL, paste0(x.good, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
        if (i %in% c(2, 5) & !excel) cat("\\midrule\n")
        if (i == 3 & excel) cat("\n")
        if (i == 6) {if (excel) cat("\n\n") else cat("%\n%\n")}
      }
    }
    return(invisible(x))
  } else if (type == "bytype") {
    for (p in c(0, 3)) { # First 3 rows for all years, last 3 rows for all years
      for (j in 1:l) {
        cat(if (rownames) paste0(time.labels[j], sep) else NULL, paste0(xn, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
        for (i in (1:3)+p) {
          this.x <- x[[j]][i, ]
          x.good <- myFormat(this.x, zero.replace, digits, excel, sums)
          cat(if (rownames) words[i] else NULL, paste0(x.good, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
          if (i %in% c(2, 5) & !excel) cat("\\midrule\n")
          if (i %in% c(3, 6) & excel) cat("\n")
        }
      }
      if (excel) cat("\n\n") else cat("%\n%\n")
      return(invisible(x))
    }
  } else if (type == "twocolumn") {
    ST <- do.call(cbind, lapply(x, function(y) y[1, ]))
    LT <- do.call(cbind, lapply(x, function(y) y[2, ]))
    if (sums) {
      ST <- rbind(ST, Total = colSums(ST))
      LT <- rbind(LT, Total = colSums(LT))
    }
    cat(if (francais) "Court terme" else "Short term", "\n", sep = "")
    cat(if (rownames) sep else NULL, paste0(time.labels, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
    for (i in 1:nrow(ST)) {
      this.x <- ST[i, ]
      x.good <- myFormat(this.x, zero.replace, digits, excel, sums = FALSE)
      cat(if (rownames) paste0(xn[i], sep) else NULL, paste0(x.good, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
    }
    cat("\n\n")
    cat(if (francais) "Long terme" else "Long term", "\n", sep = "")
    cat(if (rownames) sep else NULL, paste0(time.labels, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
    for (i in 1:nrow(LT)) {
      this.x <- LT[i, ]
      x.good <- myFormat(this.x, zero.replace, digits, excel, sums = FALSE)
      cat(if (rownames) paste0(xn[i], sep) else NULL, paste0(x.good, collapse = sep), if (!excel) "\\\\\n" else "\n", sep = "")
    }
    return(invisible(list(ST = ST, LT = LT)))
  }

}

#' Plot observed and forecast series with a distinction
#'
#' @param vars A character vector of variable names to plot as solid lines.
#' @param vars.fc A character vector of forecast variables names to plot as dashed lines.
#' @param data A aata frame containing all the \code{vars} and \code{vars.fc} variables.
#' @param expand.ylim A numeric scalar: increase the bottom white space for the legend by this share of ylim
#' @param x A numeric vector of time-dimension values (usually the time). If `NULL`, integer indices are used.
#' @param col A character vector of line colours
#' @param centre If TRUE, subtract the mean of each variable.
#' @param scale If TRUE, divide by the SD of each variable.
#' @param ... Passed to plot.
#'
#' @return Nothing (invisible NULL).
#' @export
#'
#' @examples
#' set.seed(1)
#' xfc <- matrix(cumsum(rnorm(120)), ncol = 4) # Complete data
#' x <- xfc
#' x[c(26:30, 57:60, 89:90, 111:120)] <- NA # Adding incompleteness
#' colnames(x) <- paste0("X", 1:4); #' colnames(xfc) <- paste0("fX", 1:4)
#' d <- cbind(x, xfc)
#' plotWithForecast(vars = colnames(x), vars.fc = colnames(xfc), data = d)
plotWithForecast <- function(vars, vars.fc, data,
                             expand.ylim = 0.1,
                             x = NULL,
                             col = NULL,
                             centre = TRUE, scale = TRUE,
                             ... # Passed to plot()
) {
  if (is.null(col)) col <- seq_along(vars)
  if (is.null(x)) x <- 1:nrow(data)
  if (length(vars) != length(vars.fc)) stop("The variables and their forecast counterparts must have the same lengths.")

  withr::local_par(mar = c(2, 0.7, 0.2, 0.7))
  sx <- data[, vars]
  fx <- data[, vars.fc]
  if (centre) {
    m1x <- colMeans(fx, na.rm = TRUE)
    sx <- sweep(sx, 2, m1x, "-")
    fx <- sweep(fx, 2, m1x, "-")
  }
  if (scale) {
    m2x <- apply(fx, 2, stats::sd, na.rm = TRUE)
    sx <- sweep(sx, 2, m2x, "/")
    fx <- sweep(fx, 2, m2x, "/")
  }
  xr <- range(sx, fx, na.rm = TRUE)
  # Expanding the plotting range from below to fit the legend
  xr <- c(xr[1] - expand.ylim*diff(xr), xr[2])
  plot(NULL, NULL, xlim = range(x), ylim = xr, bty = "n", ylab = "", xlab = "", yaxt = "n", ...)
  for (i in 1:ncol(sx)) {
    graphics::lines(x, sx[, i], lwd = 3, col = col[i])
    na.inds <- which(is.na(sx[, i]))
    graphics::lines(x, fx[, i], lwd = 2, lty = 2, col = col[i])
    if (length(na.inds) > 0) graphics::points(x[na.inds], fx[na.inds, i], pch = 8, lwd = 2, col = col[i])
  }
  graphics::abline(h = 0, lty = 2)
  graphics::legend("bottom", vars, ncol = 3, col = col, pt.cex = 1.5, pch = 15)
  return(invisible(NULL))
}

#' Convert ECM coefficients to an equivalent ARDL
#'
#' @param dx.coef.list A list of numeric vectors of coefficients on contemporaneous and lagged regressor differences (\code{Dx[t-0]}, \code{Dx[t-1]}, ...).
#' @param dy.coef.vec A numeric vector of coefficients on lagged dep. var. differences (\code{Dy[t-1]}, \code{Dy[t-2]}, ...).
#' @param ECT A numeric scalar: coefficient on the ECT (adjustment strength). For stable models, must be between -2 and 0 (0 means that there is no error correction).
#' @param long.run A numeric vector of coefficients on the lagged regressor levels (\code{x[t-1]}).
#' @param varnames A character vector for human-readable variable names.
#'
#' @return A list of AR coefficient vector and DL coefficient list.
#' @export
#'
#' @examples
#' # Imagine an error-correction model:
#' # dY = 0.2*dY(-1) + 0.9*dX1 + 0.2*dX2 - 0.8(Y(-1) - 1.5*X1(-1) + 0.5*X2(-1)) + U
#' # Its ARDL equivalent is
#' # Y = 0.4*Y(-1) - 0.2*Y(-2) + 0.9*X1 + 0.2*X2 + 0.3*X1(-1) - 0.6*X2(-1) + U
#' rho <- 0.2
#' dx <- c(x1 = 0.9, x2 = 0.2)
#' ECT <- -0.8
#' lx <- c(x1 = -1.5, x2 = 0.5)
#' ARDL.coef <- ECM2ARDL(dx.coef.list = dx, dy.coef.vec = rho, ECT = ECT, long.run = lx)
ECM2ARDL <- function(dx.coef.list = NULL, dy.coef.vec = NULL, ECT, long.run, varnames = NULL) {
  if (is.null(dx.coef.list)) {
    dx.coef.list <- list(numeric(length(long.run)))
    warning("There is no dX in the model, which is an atypical modelling choice.\nExogenous regressors are only in the EC term.")
  }
  if (is.null(dy.coef.vec)) {
    warning("There is no dY in the model, which is an atypical modelling choice.\nThe error-correction mechanism may be unstable.\nThe auto-regressive part is coming only from the EC term.")
  }
  if (is.vector(dx.coef.list, mode = "numeric")) dx.coef.list <- list(dx.coef.list)
  x.order <- length(dx.coef.list) # At least one because there is the EC term
  y.order <- length(dy.coef.vec)
  ar.order <- y.order + 1
  dl.order <- x.order
  k <- length(long.run)
  if (!all(sapply(dx.coef.list, length) == k)) stop("The list of coefficients on exogenous lags must contain elements of the same length as the vector of long-run elasticities.\nDid you forget the intercept as the 1st element of each vector?")

  ar.vec <- numeric(ar.order)
  ar.vec[1] <- 1 + ECT + if (ar.order > 1) dy.coef.vec[1] else 0
  if (length(ar.vec) > 1) ar.vec[length(ar.vec)] <- -dy.coef.vec[length(dy.coef.vec)]
  if (length(ar.vec) > 2) ar.vec[-c(1, length(ar.vec))] <- dy.coef.vec[-1] - dy.coef.vec[-length(dy.coef.vec)]
  names(ar.vec) <- paste0("lag", seq_along(ar.vec), "_Y")

  if (is.null(varnames)) varnames <- names(long.run)
  varnames <- gsub("^lag1_", "", varnames)
  dl.list <- vector(mode = "list", dl.order + 1) # Because the DL part starts from 0
  dl.list[[1]] <- dx.coef.list[[1]]
  dl.list[[2]] <- -dx.coef.list[[1]] + ECT*long.run
  if (dl.order > 1) {
    for (i in 2:dl.order) {
      dl.list[[i]]   <- dl.list[[i]] + dx.coef.list[[i]]
      dl.list[[i+1]] <- -dx.coef.list[[i]]
    }
  }
  dl.list <- lapply(dl.list, function(x) {names(x) <- varnames; return(x)})
  names(dl.list) <- c("X[t]", paste0("X[t-", 1:dl.order, "]"))
  cat("Converting an ECMX(", y.order, ", ", x.order, ") into an ARDL(", ar.order, ", ", dl.order, ")\n", sep = "")
  cat("Y = ", paste0("r", 1:ar.order, "*Y(-", 1:ar.order, ")", collapse = " + "), " + ", gsub("\\(-0\\)", "", paste0("b", 0:dl.order, "'X(-", 0:dl.order, ")", collapse = " + ")), " + U\n" , sep = "")

  return(list(AR = ar.vec, DL = dl.list))
}

#' Relative importance of regressors in the Elastic Net path
#'
#' @rdname varimp
#' @param x An estimated glmnet object.
#' @param max.vars An integer: how many variables to return. If NULL, returns all model variables.
#' @param sorted If TRUE, sorts the variable by relative importance, otherwise use the original name order
#'
#' @return A numeric vector with popularity values: 0 means that Elastic Net drops the variable even with small penaltiles, values close to 1 means that this variable was retained in the model for most of the path.
#' @export
#'
#' @examples
#' x <- model.matrix(mpg ~ ., data = datasets::mtcars)[, -1]
#' y <- datasets::mtcars[, "mpg"]
#' lasso <- glmnet::glmnet(x = x, y = y)
#' elnet <- glmnet::glmnet(x = x, y = y, alpha = 0.5)
#' vi1 <- varImportance(lasso); vi2 <- varImportance(elnet)
#' plotImportance(cbind(vi1, vi2), pch = 1)
#' plotCoefBetter(lasso)
varImportance <- function(x, max.vars = NULL, sorted = TRUE) {
  b <- Matrix::as.matrix(glmnet::coef.glmnet(x))
  ncoef <- rowSums(b != 0) # For how many penalties these appear non-zero
  ncoef <- ncoef[names(ncoef) != "(Intercept)"]
  if (sorted) popularity <- sort(ncoef, decreasing = TRUE) else popularity <- ncoef
  popularity <- popularity / max(popularity)
  if (!is.null(max.vars)) {
    inds <- if (sorted) match(1:max.vars, order(popularity, decreasing = TRUE)) else
      which(order(popularity, decreasing = TRUE) %in% 1:max.vars)
    popularity <- popularity[inds]
  }
  return(popularity)
}

#' @rdname varimp
#' @param vi A vector or a matrix with variable names as rows and percentages in columns.
#' @param col A vector of colours (for each columnif `vi` is a matrix).
#' @param ... Passed to plot.
#'
#' @return Nothing (invisible NULL).
#' @export
plotImportance <- function(vi, col = NULL,
                           ... # Passed to plot
) {
  x <- vi
  if (is.null(dim(x))) {
    n <- names(x)
    x <- matrix(x, ncol = 1)
    rownames(x) <- n
    single.vec <- TRUE
	} else single.vec <- FALSE
  if (is.null(rownames(x))) stop("There are no row names, the plot will be meaningless.")
  if (is.null(colnames(x)) & !single.vec) stop("There are no column names, cannot add a legend.")
  if (is.null(col)) col <- 1:ncol(x)
  ys <- rev(1:nrow(x))

  plt.args <- list(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(1, nrow(x)),
                   yaxt = "n", bty = "n", xlab = "", ylab = "", pch = 16, lty = 3,
                   xaxt = "n", sub = "% models with non-zero coef. in the regularis. path")
  dotlist <- list(...)
  if (length(dotlist) > 0) plt.args[names(dotlist)] <- dotlist
  do.call(plot, plt.args)
  rmax <- apply(x, 1, max)
  for (i in 1:nrow(x)) graphics::lines(c(0, rmax[i]), rep(ys[i], 2), col = "#00000055", lty = plt.args$lty)
  for (i in 1:ncol(x)) graphics::points(x[, i], ys, pch = plt.args$pch, col = col[i])
  graphics::axis(1, seq(0, 1, 0.2))
  graphics::axis(2, nrow(x):1, labels = rownames(x), las = 2)
  if (ncol(x) > 1) graphics::legend("bottomright", legend = colnames(x), col = col, bty = "n", pch = plt.args$pch)
  invisible(return(NULL))
}

#' Prints hypothesis testing results in a concise manner
#'
#' @param x An object of class anova (from car::linearHypothesis) or htest (from nlWaldTest::nlWaldlest)
#' @param hyp Optional: human-readable hypothesis name
#' @param digits Number of p-values digits to print.
#'
#' @return Nothing (invisible NULL).
#' @export
#'
#' @examples
printShort <- function(x, hyp = NULL, digits = 3) {
  if ("anova" %in% class(x)) {
    which.eq <- grep("=", attr(x, "heading"))
    hh <- attr(x, "heading")[which.eq]
    hh <- gsub("[ \t]+", " ", hh)
    if (length(hh) == 1) h <- paste0("Hypothesis: ", hh) else {
      h <- paste0("Hypotheses:\n", paste0("(", 1:length(hh), ") ", hh, collapse = "\n"))
    }
    df <- x[2, 2]
    f <- x[2, 3]
    p <- x[2, 4]
  } else if ("htest" %in% class(x)) {
    h <- hyp
    if (is.null(h)) h <- paste0("Parametric hypothesis test for ", x$data.name) else {
      if (length(h) == 1) h <- paste0("Hypothesis: ", h) else h <- paste0("Hypotheses:\n", paste0("(", 1:length(h), "): ", h, collapse = "\n"))
    }
    df <- unname(x$parameter)
    f <- unname(x$statistic)
    p <- unname(x$p.value)
  } else stop("This function works only for car::linearHypothesis and nlWaldTest::nlWaldtest")
  cat("\n", h, "\n", "df=", df, ", F=", sprintf("%1.2f", f), ", p=", sprintf(paste0("%1.", digits, "f"), p), "\n", sep = "")
  return(invisible(NULL))
}

#' Extract the error-correction term from an ECM
#'
#' @param x An ECM (lm or restriktor).
#' @param yname Lagged dependent variable name (e.g. `"lag1_y"`).
#' @param xnames Lagged dependent variable names (e.g. `c("lag1_x1", "lag1_x2")`).
#' @param centre If TRUE, subtract the mean. Use for unrestricted constant.
#' @param detrend If TRUE, subtract the linear trend Use for unrestricted trend
#'
#' @return A numeric vector of long-run residuals (Y[t-1] - X[t-1]'b).
#' @export
#'
#' @examples
getECT <- function(x, yname, xnames, centre = TRUE, detrend = FALSE) {
  m <- if ("restriktor" %in% class(x)) stats::model.matrix(x[["model.org"]]) else stats::model.matrix(x)
  b <- if ("restriktor" %in% class(x)) x[["b.restr"]] else stats::coef(x)
  ec <- m[, c(yname, xnames)] %*% b[c(yname, xnames)]
  ec <- as.numeric(ec) / b[yname]
  if (centre) ec <- ec - mean(ec)
  if (detrend) ec <- unname(stats::lm(ec ~ seq_along(ec))$residuals)
  names(ec) <- rownames(m)
  return(ec)
}

#' @describeIn printECM Extract the long-run elasticities
#' @export
getLR <- function(x, yname, xnames) {
  b <- if (is.numeric(x)) x else stats::coef(x)
  elast <- b[xnames] / b[yname]
  return(elast)
}

#' Print one or multiple ECMs in an easily interpretable format
#'
#' @param ... Estimated model outputs from lm, restriktor::restriktor or restriktor::conLM.
#' @param yname A string with the dependent variable name.
#' @param xnames A character vector of explanatory variable names that appear in both parts and differ by prefixed. Unique variables (like time dummies) need not be included.
#' @param d.prefix A string that should precede every differenced variables (usually d_ or U+0394 if the locale is Unicode)
#' @param lag.prefix A string that should to precede first lags.
#' @param mod.names A character vector of model names for printing.
#' @param shorten A function to format numbers (usually round, sprintf or formatC).
#' @param select A character vector showing which parts of the equation to keep for reporting.
#' @param minus Logical: multiply the OLS long-run estimates by -1?
#'
#' @return Invisibly returns the matrices that have just been shown (without rounding).
#' @export
#'
#' @examples
printECM <- function(..., yname, xnames,
                     d.prefix = "d_", lag.prefix = "lag1_", mod.names = NULL,
                     shorten = function(x) formatC(x, format = "f", digits = 2, flag = " "),
                     select = c("ec", "long", "short", "inertia", "intercept"),
                     minus = TRUE) {
  x <- list(...)
  # print(str(x, max.level = 1))
  if (length(x) == 1) { # One argument was passed
    # Possible input: list of models or one model (which is a list, but the class is "lm")
    objnames <- names(x[[1]])
    # A linear model or conLM should have at least coefficients and fitted values
    # Unless a user supplied a model list where two of the models are names
    # literally "fitted" (or "fitted.values") and "residuals"
    if (("residuals" %in% objnames) & any(c("fitted", "fitted.values") %in% objnames)) {
      x <- x[[1]]
      cat("ECM summary for one ", paste0(class(x), collapse = "-"), " model\n", sep = "")
      if (any(c("lm", "restrikor", "conLM") %in% class(x))) { # Print a single ECM
        out <- getECM(x = x, yname = yname, xnames = xnames, d.prefix = d.prefix, lag.prefix = lag.prefix)
        fx  <- cbind(Short = shorten(out$short), Long = shorten(out$long * if (minus) -1 else 1))
        # rownames(fx) <- xnames
        print(fx, quote = FALSE)
        cat("\nError correction strength: ", shorten(out$EC), ".\nAuto-regressive term: ", shorten(out$inertia), ", intercept: ", shorten(out$Intercept), ".\n", sep = "")
        return(invisible(NULL))
      } else stop("This function supports only lm, conLM and restriktor objects.")
    }
  }
  # Now that the one-model case has been treated, we check the case:
  # if a model list was supplied and redundantly wrapped in list(...),
  # we extract it; the single element is not a model (we checked above)
  if (length(x) == 1) x <- unlist(x, recursive = FALSE)
  xclasses <- sapply(x, function(x) any(c("lm", "restrikor", "conLM") %in% class(x)))
  if (!all(xclasses)) stop(paste0("This function supports only lm, conLM and restriktor objects.\nProblematic object indices: ", paste0(which(!xclasses), collapse = ", ")))
  coefs <- lapply(x, function(a) getECM(x = a, yname = yname, xnames = xnames, d.prefix = d.prefix, lag.prefix = lag.prefix))
  uniq.names.short <- unique(unlist(lapply(coefs, function(x) names(x[["short"]]))))
  uniq.names.long <- unique(unlist(lapply(coefs, function(x) names(x[["long"]]))))
  EC <- do.call(cbind, lapply(coefs, "[[", "EC"))
  Intercept <- do.call(cbind, lapply(coefs, "[[", "Intercept"))
  inertia <- do.call(cbind, lapply(coefs, "[[", "inertia"))
  short <- array(NA, dim = c(length(uniq.names.short), length(coefs)), dimnames = list(uniq.names.short, NULL))
  long <-  array(NA, dim = c(length(uniq.names.long),  length(coefs)), dimnames = list(uniq.names.long, NULL))
  for (i in 1:length(coefs)) {
    cshort <- coefs[[i]][["short"]]
    clong <- coefs[[i]][["long"]]
    short[names(cshort), i] <- cshort
    long[names(clong), i]  <- clong * if (minus) -1 else 1
  }
  # Dropping the variables that are not used in any of the model
  short.all.na <- apply(short, 1, function(x) all(!is.finite(x)))
  long.all.na  <- apply(long, 1, function(x) all(!is.finite(x)))
  short <- short[!short.all.na, , drop = FALSE]
  long <- long[!long.all.na, , drop = FALSE]

  colnames(inertia) <- colnames(short) <- colnames(long) <- colnames(EC) <- colnames(Intercept) <- mod.names
  if ("all" %in% select) select <- c("ec", "long", "short", "inertia", "intercept")
  if ("ec" %in% select) {
    cat("\nCoefficient on the ECT (correction strength):\n")
    print(shorten(EC), quote = FALSE)
  }
  if ("short" %in% select) {
    cat("\nShort-run elasticities:\n")
    print(shorten(short), quote = FALSE)
  }
  if ("long" %in% select) {
    cat("\nLong-run elasticities:\n")
    print(shorten(long), quote = FALSE)
  }
  if ("inertia" %in% select) {
    cat("\nAuto-regressive term (intertia):\n")
    print(shorten(inertia), quote = FALSE)
  }
  if ("intercept" %in% select) {
    cat("\nIntercepts:\n")
    print(shorten(Intercept), quote = FALSE)
  }

  return(invisible(list(ec = EC, short = short, long = long, inertia = inertia, intercept = Intercept)))
}


#' @describeIn printECM Internal function to sort ECM components
#' @param x An lm model or a numeric vector
#' @param na.to.zero If TRUE, replace NAs with zeros.
#' @param incorporate.intercept If TRUE, report the intercept in the short-run part.
#' @export
getECM <- function(x, yname, xnames, d.prefix, lag.prefix,
                   na.to.zero = FALSE, incorporate.intercept = FALSE) {
  b <- if (is.numeric(x)) x else stats::coef(x)
  xnames <- unique(xnames) # Safety check to ensure that setdiff() worke properly

  xx <- setdiff(names(b), "(Intercept)")

  # Common vars enter both the LR and SR part
  common.vars <- xnames[(paste0(d.prefix, xnames) %in% xx) & (paste0(lag.prefix, xnames) %in% xx)]

  sr   <- paste0(d.prefix, common.vars) # Short-run effects
  b.sr <- b[sr]; names(b.sr) <- common.vars

  ac   <- xx[grep(paste0("^", d.prefix, ".+", yname), xx)] # Auto-correlations
  b.ac <- b[ac]
  names(b.ac) <- NULL

  lr <- paste0(lag.prefix, common.vars)
  b.lr  <- getLR(x = x, yname = paste0(lag.prefix, yname), xnames = lr)
  names(b.lr) <- common.vars

  ec <- paste0(lag.prefix, yname)
  b.ec <- b[ec]
  names(b.ec) <- NULL

  # The extra variables may have a short-run prefix in unrestricted VAR
  # (unrestricted in the SR equation, seed Pesaran (2001) for a discussion)
  extra.vars <- setdiff(xx, c(sr, lr, ac, ec))
  b.rest <- b[extra.vars]

  dx <- lx <- rep(NA, length(common.vars)+length(extra.vars))
  names(dx) <- names(lx) <- c(common.vars, extra.vars)

  if (!isTRUE(is.finite(ac))) ac <- NA
  dx[c(common.vars, extra.vars)] <- c(b.sr, b.rest)
  lx[common.vars] <- b.lr

  if (incorporate.intercept) {
    dx["(Intercept)"] <- unname(b["(Intercept)"])
    lx["(Intercept)"] <- 0
  }

  if (na.to.zero) {
    b.ac[is.na(b.ac)] <- 0
    dx[is.na(dx)] <- 0
    lx[is.na(lx)] <- 0
    b.ec[is.na(b.ec)] <- 0
  }

  ret <- list(inertia = b.ac, short = dx, long = lx, EC = b.ec, Intercept = unname(b["(Intercept)"]))
  if (incorporate.intercept) ret$Intercept <- NULL

  return(ret)
}

#' Compute lags and differences of a variable for an ECM
#'
#' @param varname Character: variable name
#' @param data Data frame where this variable can be found
#' @param nlag Number of lags to compute
#'
#' @return A data frame with the lags and differences
#'
#' @examples
#' d <- cbind(X1 = rnorm(100), X2 = cumsum(rnorm(100)))
#' prepareTransform("X2", d)
prepareTransform <- function(varname, data, nlag = 2) {
  x <- data[, varname]
  lx <- sapply(0:nlag, function(i) myLag(x, lag = i))
  colnames(lx) <- c(varname, paste0("lag", 1:nlag, "_", varname))
  dx <- apply(lx, 2, myDiff)
  colnames(dx) <- paste0("d_", colnames(lx))
  return(cbind(lx[, -1, drop = FALSE], dx))
}

#' Perform Adaptive Elastic Net with optimal adaptive weights
#'
#' @param x Passed to glmnet: a matrix of predictors without the constant (all finite values).
#' @param y Passed to glmnet: a numeric vector of dependent variable values (all finite).
#' @param weights Passed to glmnet: observation weights
#' @param adaptive If TRUE, use adaptive regressor weights
#' @param adaptive.gamma Adaptive weight power in 1/|beta_initial|^gamma; a positive number, usually 1, 2 or 0.5.
#' @param upper.limits Passed to glmnet: box constraints from above.
#' @param lower.limits Passed to glmnet: box constraints from below.
#' @param alphas Numeric vector of the sequence of elastic net mixing parameters (1 = LASSO, 0 = ridge).
#' @param lambda Optional: a pre-defined sequence of elastic net penalty lambdas if one alpha is chosen.
#' @param nlambda Passed to glmnet: the number of lambdas to choose over
#' @param parallel Logical: parallelise computations?
#' @param cluster A cluster for parallel computations to be passed to parLapplyLB.
#' @param plot.cv Logical: if TRUE, produce a contour plot of the cross-validation stage.
#' @param time.series Character indicating whether the CV should be done via LOO (default), sliding window, or expanding window
#' @param ts.out.sample If 'time.series' is 'expanding' or 'sliding', the number of points at the end to use for testing
#' @param se.mult In times series cross-validation, the margin for increasing the penalty within this number of forecast MSE SDs.
#' @param best.choice Which model is considered the best: the one that minimises MSE or a more parsimonious one within the 'se.mult' limit
#' See Chen & Yang (2021) The One Standard Error Rule for Model Selection: Does It Work? for a discussion (in short, "min" should be better).
#' @param full.relax Cross-validate with respect to the MSE from the relaxed fit?
#' @param ... Passed to contour that is used for plotting.
#'
#'
#' @return A list of optimal adaptive weights, estimated Elastic Net, cross-validated cv.glmnet, optimal coefficients, optimal alpha, and optimal lambda.
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- arima.sim(list(ar = c(0.4, 0.3)), n = 100)
#' x2 <- arima.sim(list(ar = c(0.4, 0.3)), n = 100)
#' x3 <- arima.sim(list(ar = c(0.4, 0.3)), n = 100)
#' y <- 1 + x1 + x2 + rnorm(100)
#' x <- cbind(x1, x2, x3)
#' a <- AdaEN(x = x, y = y, plot.cv = TRUE, time.series = "expanding",
#'            ts.out.sample = 15, alphas = 0:5/5) # Takes 15 seconds to run
AdaEN <- function(x, y, weights = NULL,
                  adaptive = TRUE, adaptive.gamma = 1,
                  lower.limits = -Inf, upper.limits = Inf,
                  alphas = c(0, 0.01, 0.05, seq(0.1, 1, 0.1)), lambda = NULL, nlambda = 51,
                  parallel = FALSE, cluster = NULL, plot.cv = FALSE,
                  time.series = c("no", "sliding", "expanding"), ts.out.sample = 10,
                  se.mult = NULL, best.choice = c("min", "se"), full.relax = FALSE,
                  ...) {
  plot.cv.glmnet <- utils::getFromNamespace("plot.cv.glmnet", "glmnet")
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(se.mult)) se.mult <- min(0.1*sqrt(ts.out.sample), 1)
  time.series <- time.series[1]
  if (!(time.series %in% c("no", "sliding", "expanding"))) stop("AdaEN: 'time.series' sohuld be 'no', 'sliding' or 'expanding'.")
  if ((time.series != "no") & (nrow(x) - ts.out.sample < 5)) stop("AdaEN: estimation in the first sub-sample has fewer than 5 points. Reduce ts.out.sample.")
  best.choice <- best.choice[1]
  if (!(best.choice %in% c("min", "se"))) stop("AdaEN: 'best.choice' should be either 'min' or 'se'.")
  x0 <- x # Saving the series with all attributes
  y0 <- y
  # If no ts analysis is requested but there are time attributes, they must
  # be stripped before analysis
  if (stats::is.ts(x) | !is.null(stats::tsp(x))) {stats::tsp(x) <- NULL; class(x) <- setdiff(class(x), c("mts", "ts"))}
  if (stats::is.ts(y) | !is.null(stats::tsp(y))) {stats::tsp(y) <- NULL; class(y) <- setdiff(class(y), c("mts", "ts"))}

  penalty.factor <- rep(1, ncol(x))

  doCV <- function(alph) { # Ensuring that there are nlambda suitable lambdas
    # The lambdas are always decreasing in the implementation
    # A sequence of lambdas must be generated in advance because sometimes, it is dropped
    mod.for.lambda <- glmnet::glmnet(x = x, y = y, weights = weights, alpha = alph, nlambda = 5, nfolds = nrow(x), penalty.factor = penalty.factor, relax = full.relax, upper.limits = upper.limits, lower.limits = lower.limits)
    lseq <- exp(seq(log(max(mod.for.lambda$lambda)), log(min(mod.for.lambda$lambda)), length.out = nlambda))
    out  <- glmnet::cv.glmnet(x = x, y = y, weights = weights, alpha = alph, lambda = lseq, nfolds = nrow(x), grouped = FALSE, penalty.factor = penalty.factor, relax = full.relax, gamma = 0, upper.limits = upper.limits, lower.limits = lower.limits)
    if (time.series != "no") {
      # Preparing data slices for time-series OOS forecasting
      start.ind <- if (time.series == "sliding") 1:ts.out.sample else rep(1, ts.out.sample)
      end.ind <- nrow(x) - ts.out.sample:1
      oos.ind <- end.ind + 1 # Out-of-sample indices
      x.train.list <- sliceList(x, start = start.ind, end = end.ind)
      y.train.list <- sliceList(y, start = start.ind, end = end.ind)
      x.test.list  <- sliceList(x, start = oos.ind, end = oos.ind)
      y.test.list  <- sliceList(y, start = oos.ind, end = oos.ind)
      oos.stat     <- lapply(1:ts.out.sample, function(i) {
        m <- glmnet::glmnet(x = x.train.list[[i]],  y = y.train.list[[i]],
                            weights = weights[attr(x.train.list[[i]], "indices")],
                            alpha = alph, lambda = lseq, penalty.factor = penalty.factor, relax = full.relax,
                            upper.limits = upper.limits, lower.limits = lower.limits)
        yhat <- if (full.relax) glmnet::predict.relaxed(m, newx = x.test.list[[i]], gamma = 0) else glmnet::predict.glmnet(m, newx = x.test.list[[i]])
        err <- y.test.list[[i]] - as.numeric(yhat)
        return(list(yhat = yhat, err = err))
      })
      oos.fcast  <- do.call(cbind, lapply(oos.stat, function(x) as.numeric(x[["yhat"]])))
      oos.errors <- do.call(cbind, lapply(oos.stat, function(x) as.numeric(x[["err"]])))
      # Substituting the key indicators in CV
      err.sq <- oos.errors^2
      cvm <- rowMeans(err.sq)
      cvsd <- apply(err.sq, 1, stats::sd) / sqrt(ts.out.sample)
      cvsd <- cvsd * se.mult
      # The next part is taken from glmnet:::getOptcv.glmnet
      out$cvm <- cvm
      out$cvsd <- cvsd
      out$cvup <- cvm + cvsd
      out$cvlo <- cvm - cvsd
      idmin <- which.min(cvm)
      lambda.min <- lseq[idmin]
      semin <- (cvm + cvsd)[idmin]
      id1se <- cvm <= semin
      lambda.1se <- max(lseq[id1se], na.rm = TRUE)
      id1se <- match(lambda.1se, lseq)
      out$lambda.min <- lambda.min
      out$lambda.1se <- lambda.1se
      out$index <- matrix(c(idmin, id1se), 2, 1, dimnames = list(c("min", "1se"), "Lambda"))
      yhat <- if (best.choice == "min") oos.fcast[idmin, ] else oos.fcast[id1se, ]
      names(yhat) <- rownames(x)[oos.ind]
    } else {
      yhat <- as.numeric(glmnet::predict.glmnet(out$glmnet.fit, newx = x, s = if (best.choice == "min") out$lambda.min else out$lambda.1se))
      names(yhat) <- rownames(x)
    }
    # cat("dim(oos.fcast) =", dim(oos.fcast), "\n")

    attr(out, "window")   <- time.series
    attr(out, "forecast") <- yhat
    return(out)
  }

  if (adaptive) { # Ridge regression estimates to be used in adaptive weights
    CV.ridge <- doCV(0)
    b.ridge  <- abs(glmnet::coef.glmnet(glmnet::glmnet(x = x, y = y, weights = weights, alpha = 0, lambda = CV.ridge$lambda.min)))[-1]
    if (adaptive.gamma <= 0) stop("adaptive.gamma must be positive (0.5, 1, or 2 is a good choice).")
    if (adaptive.gamma != 1) b.ridge <- b.ridge^adaptive.gamma
    penalty.factor <- 1 / b.ridge * stats::median(b.ridge) # To avoid large numbers
  }

  if (is.null(lambda)) {
    if (length(alphas) > 1) { # Choose alpha by CV
      if (parallel & !is.null(cluster)) {
        parallel::clusterExport(cluster, "sliceList")
        parallel::clusterExport(cluster, c("x", "y", "weights", "penalty.factor", "full.relax", "time.series", "nlambda", "upper.limits", "lower.limits"), envir = environment())
        CV <- parallel::parLapplyLB(cl = cluster, X = alphas, fun = doCV)
      } else {
        CV <- lapply(X = alphas, FUN = doCV)
      }
      CV.MSE <- sapply(CV, function(x) x[["cvm"]])
      CV.min <- cbind(alpha = alphas, do.call(rbind, lapply(CV, function(v) c(lambda = v$lambda.min, index = nlambda+1-v$index[1, 1], CV = unname(v$cvm[v$index[1, 1]])))))
      CV.1SE <- cbind(alpha = alphas, do.call(rbind, lapply(CV, function(v) c(lambda = v$lambda.1se, index = nlambda+1-v$index[2, 1], CV = unname(v$cvm[v$index[2, 1]])))))
      best.a.ind <- which.min(CV.min[, "CV"]) # Which alpha minimises the vanilla MSE
      a.best <- alphas[best.a.ind]
      if (a.best == 0) {
        warning("The optimal (by MSE) model is Ridge (alpha=0), thus no variable selection is happening (no LASSO penalty).\nSetting alpha=0.01 (Elastic Net) to allow at least some chance of zero coefficients.")
        a.best <- 0.01
        if (a.best %in% alphas) CV.best.a <- CV[[which(alphas == a.best)]] else CV.best.a <- doCV(0.1)
      } else {
        CV.best.a <- CV[[best.a.ind]]
      }
      if (plot.cv) {
        # The lambdas are decreasing
        CV.forplot <- log(CV.MSE[nrow(CV.MSE):1, ])
        graphics::contour(x = 1:nlambda, y = alphas, z = CV.forplot,
                          levels = stats::quantile(CV.forplot, c(0.01, 0.02, 1:19/20, 0.99)),
                          labels = c(0.01, 0.02, 1:19/20, 0.99),
                          col = grDevices::rainbow(22, end = 0.7, v = 0.75, rev = TRUE),
                          xaxt = "n", yaxt = "n", bty = "n", xlab = "Lambda multiplier", ylab = "Alpha", ...)
        graphics::lines(CV.min[, "index"], alphas, lty = 2, lwd = 1.5)
        graphics::lines(CV.1SE[, "index"], alphas, lty = 3, lwd = 1.5)
        best.inds <- nlambda + 1 - as.numeric(CV.best.a$index)
        graphics::points(best.inds, rep(a.best, 2), pch = c(16, 1), cex = c(1, 1.75))
        graphics::axis(1, c(1, floor(nlambda/2), nlambda), c("0.01", "1", "100"))
        graphics::axis(2, 0:5/5, c("Ridge", 1:4/5, "LASSO"))
      }
      lambda <- if (best.choice == "min") CV.min[best.a.ind, "lambda"] else CV.1SE[best.a.ind, "lambda"]
      lambda <- unname(lambda)
    } else { # No CV for alpha
      a.best <- alphas
      CV.best.a <- doCV(a.best)
      lambda <- if (best.choice == "min") CV.best.a$lambda.min else CV.best.a$lambda.1se
      if (plot.cv) plot.cv.glmnet(CV.best.a)
    }

    best.fit <- CV.best.a$glmnet.fit

  } else { # Lambda was given
    if (length(alphas) > 1) stop("Hand-pick only 1 alpha and 1 lambda.")
    a.best <- alphas
    CV.best.a <- NULL
    best.fit <- glmnet::glmnet(x = x,  y = y, weights = weights, alpha = a.best, lambda = lambda, penalty.factor = penalty.factor, relax = full.relax, upper.limits = upper.limits, lower.limits = lower.limits)
  }

  b    <- glmnet::coef.glmnet(best.fit, s = lambda)
  b.bf <- as.numeric(b)
  names(b.bf) <- rownames(b)

  return(list(ada.weights = penalty.factor, AdaEN = best.fit, AdaCV = CV.best.a,
              coefficients = b.bf, alpha = a.best, lambda = lambda, relax = full.relax, weights = weights,
              forecast = attr(CV.best.a, "forecast")))
}

#' Generate a list of sliding-window data frames
#'
#' @param x A data frame, matrix, or numeric vector.
#' @param start Integer: indices of the sample starts.
#' @param end Integer: indices of the sample ends
#'
#' Create a list of data frames, matrices, or vectors, each trimmed to the
#' corresponding starts and ends.
#'
#' @return A list of data frames trimmed to the start and end, with index attributes.
#' @export
#'
#' @examples
#' x <- matrix(rnorm(100), ncol = 2)
#' a1 <- sliceList(x, start = c(1, 11, 21), end = c(30, 40, 50)) # Sliding window
#' lapply(a1, attr, "indices")
#' a2 <- sliceList(x, end = c(30, 40, 50)) # Expanding window
#' lapply(a2, attr, "indices")
sliceList <- function(x, start = NULL, end = NULL) {
  if (is.null(start) & is.null(end)) return(x)
  was.vec <- FALSE
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
    was.vec <- TRUE
  }
  if (is.null(start)) start <- rep(1, length(end))
  if (is.null(end)) end <- rep(nrow(x), length(start))
  if (length(start) > 1 & length(end) == 1) end <- rep(end, length(start))
  if (length(end) > 1 & length(start) == 1) start <- rep(start, length(end))
  if (length(start) != length(end)) stop("sliceList: 'start' and 'end' should have the same length.")
  lapply(1:length(start), function(i) {
    out <- x[start[i]:end[i], , drop = FALSE]
    if (was.vec) out <- as.vector(out)
    attr(out, "indices") <- start[i]:end[i]
    return(out)
  })
}

#' Plot glmnet coefficients with better visualisation parameters
#'
#' @rdname varimp
#' @param glmnet A glmnet object.
#' @param xvar Indicator to plot: penalty factor, norm, or deviance.
#' @param label If TRUE, add labels with halo.
#' @param drop.intercept If TRUE, do not plot the intercept regularisation path.
#' @param hscale Passed to [labelsWithHalo()].
#' @param vscale Passed to [labelsWithHalo()].
#' @param ... Passed to matplot.
#'
#' Based on glmnet:::plotCoef.
#'
#' @return Nothing.
#' @export
plotCoefBetter <- function(glmnet, xvar = c("lambda", "norm", "dev"),
                           label = TRUE, drop.intercept = TRUE,
                           hscale = 0.005, vscale = 0.01,
                           ...) {
  beta <- Matrix::as.matrix(glmnet::coef.glmnet(glmnet))
  if (drop.intercept) beta <- beta[rownames(beta) != "(Intercept)", , drop = FALSE]
  nzinds <- which(apply(beta, 1, function(x) any(x != 0)))
  zinds <- setdiff(1:nrow(beta), nzinds)
  nnz <- length(nzinds)
  switch(nnz, `0` = {
    warning("No plot produced since all coefficients zero")
    return(NULL)
  }, `1` = warning("1 or less nonzero coefficients; glmnet plot is not meaningful"))
  beta <- beta[nzinds, , drop = FALSE]
  xvar <- xvar[1]
  switch(xvar, norm = {
    index = if (missing(norm)) apply(abs(beta), 2, sum) else norm
    iname = "L1 Norm"
    approx.f = 1
  }, lambda = {
    index = log(glmnet$lambda)
    iname = "Log Lambda"
    approx.f = 0
  }, dev = {
    index = glmnet$dev
    iname = "Fraction Deviance Explained"
    approx.f = 1
  })
  plt.args <- list(x = index, y = t(beta), type = "l", lty = 1, xlab = iname, ylab = "Coefficient",
                   bty = "n", col = createDistinctCols(nnz, nsamp = max(nnz*2, 100), bad.margin = 0.1), lwd = 2)
  dotlist <- list(...)
  if (length(dotlist) > 0) plt.args[names(dotlist)] <- dotlist
  do.call(graphics::matplot, plt.args)
  atdf <- pretty(index)
  prettydf <- stats::approx(x = index, y = glmnet$df, xout = atdf, rule = 2, method = "constant", f = approx.f)$y
  graphics::axis(3, at = atdf, labels = prettydf)
  graphics::rug(side = 3, index[c(NA, diff(glmnet$df) > 0)], ticksize = -0.03)
  graphics::rug(side = 3, index[c(NA, diff(glmnet$df) < 0)], col = "#FF0000", ticksize = 0.03)
  if (label) {
    xpos <- if (xvar == "lambda") min(index) else max(index)
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    labelsWithHalo(xpos, ypos, nhalo = 8, labels = as.character(nzinds), halo.col = "#FFFFFF88", cex = 0.6, pos = 4, col = plt.args$col, hscale = hscale, vscale = vscale)
  }
  message(paste0(nnz, " variables are present in the model.\nYou may add customised colours, line types and width styles to better distinguish then."))
  cn <- rownames(beta)
  message(paste0(nzinds, ": ", cn, collapse = ", "))
  return(invisible(NULL))
}


#' Find and test the active constraints in a model
#'
#' @param rmod A restricted model produced by 'restriktor'.
#' @param eps Numerical difference epsilon to test the sensitivity of p-values
#' @param vcov A variance-covariance matrix for the unrestricted model or a function to apply to the unconstrained model
#' @param denominator.vars Optional; for testing ratio hypotheses, the names of variables found in the denominator
#' Should start with the lag of Y, then the lags of long-run X, e.g. c("lagY", "lagCapital"); whichever is found first
#' will be put into the denominator, and in the LR testing, the lag of Y is always in the denominator
#'
#' If no variance-covariance matrix is supplied, uses [sandwich::vcovHAC()]. If
#' the inference is non-standard (e.g. fast convergence rates in error-correction
#' models), then, the conclusions might be wrong.
#'
#' @return A numeric matrix with the test p-values and their derivatives
#' @export
#'
#' @examples
#' m <- lm(mpg ~ factor(cyl) + disp + hp + drat + wt, data = datasets::mtcars)
#' cns <- c("disp < 0", "hp < 0", "drat > 1", "wt < -6", "wt > -100")
#' r <- restriktor::restriktor(m, paste0(cns, collapse = "\n"))
#' # Some of these constraints are inactive; some are tighter than others
#' testActiveConstraints(r, vcov = sandwich::vcovHC)
testActiveConstraints <- function(rmod,
                                  eps = 1e-4, vcov = NULL,
                                  denominator.vars = NULL) {
  if (!any(c("restriktor", "conLM") %in% class(rmod))) stop("You must supply a restriktor object.")

  if (is.null(vcov)) {
    bw <- sandwich::bwNeweyWest(rmod$model.org, prewhite = FALSE, kernel = "Bartlett")
    vcov <- sandwich::kernHAC(rmod$model.org, bw = bw, prewhite = FALSE, kernel = "Bartlett")
    warning(paste0("Covariance matrix not supplied, computing the Newey-West HAC with automatic default parameters.\nChosen bandwidth: ", round(bw, 2)))
  } else if (is.function(vcov)) {
    vcov <- vcov(rmod$model.org)
  } else if (!is.matrix(vcov)) stop("'vcov' should be a matrix or a function.")


  b <- restriktor::coef.restriktor(rmod)
  cm <- rmod$constraints[sort(rmod$iact), , drop = FALSE]
  rhsm <- rmod$rhs[sort(rmod$iact)]
  colnames(cm) <- names(b)

  # Fully restricted p-value
  p.r <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm, rhs = rhsm, vcov. = vcov)[2, 4]

  hyp.list <- vector("list", nrow(cm))
  for (i in 1:length(hyp.list)) {
    cm_i <- cm[i, , drop = FALSE]
    act.ind <- which(cm_i != 0)
    hyp.vars <- names(b)[act.ind]
    # Testing if this expression has multiple parts
    # For testing ratio sensitivity, b0/b1=c to b0/b1=c+eps,
    # we must convert to b0 + -c*b1 - eps*b1 = 0
    # hyp0 <- paste0(paste0(hyp.vals, "*", hyp.vars, collapse = " + "), " = ", rhsm[i])
    hyp0 <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm_i, rhs = rhsm[i], vcov. = vcov)
    hyp0.text <- attr(hyp0, "heading")[2]
    ind.pu  <- hyp0[2, 4]
    ind.pr <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm[-i, , drop = FALSE], rhs = rhsm[-i], vcov. = vcov)[2, 4]

    if (is.character(denominator.vars)) {
      is.denom <- hyp.vars %in% denominator.vars
      if (any(is.denom)) {
        which.denom <- which(colnames(cm_i) %in% denominator.vars)[1]
        # denom <- colnames(cm_i)[which.denom]
        add.vec <- matrix(0, nrow = 1, ncol = ncol(cm))
        add.vec[, which.denom] <- eps
        add.mat <- matrix(0, nrow = nrow(cm), ncol = ncol(cm))
        add.mat[i, ] <- add.vec
        ind.pup <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm_i + add.vec, rhs = rhsm[i], vcov. = vcov)[2, 4]
        ind.pum <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm_i - add.vec, rhs = rhsm[i], vcov. = vcov)[2, 4]
        ind.prp <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm + add.mat, rhs = rhsm, vcov. = vcov)[2, 4]
        ind.prm <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm - add.mat, rhs = rhsm, vcov. = vcov)[2, 4]
      }
    } else is.denom <- rep(FALSE, length(hyp.vars))
    # Testing the restrictions and sensitivity for individual hypotheses in the unrestricted model
    if (!exists("ind.pup")) {
      eps.vec <- numeric(length(rhsm))
      eps.vec[i] <- eps
      ind.pup <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm_i, rhs = rhsm[i]+eps, vcov. = vcov)[2, 4]
      ind.pum <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm_i, rhs = rhsm[i]-eps, vcov. = vcov)[2, 4]
      ind.prp <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm, rhs = rhsm+eps.vec, vcov. = vcov)[2, 4]
      ind.prm <- car::linearHypothesis(model = rmod$model.org, hypothesis.matrix = cm, rhs = rhsm-eps.vec, vcov. = vcov)[2, 4]
    }
    hyp.list[[i]] <- c(UR = ind.pu, URplus = ind.pup, URminus = ind.pum, R = ind.pr, Rplus = ind.prp, Rminus = ind.prm)
    attr(hyp.list[[i]], "hypothesis") <- hyp0.text
    attr(hyp.list[[i]], "isRatio") <- any(is.denom)
  }

  hyp.mat <- do.call(rbind, hyp.list)

  ind.pu.grad     <- (hyp.mat[, "URplus"] - hyp.mat[, "URminus"])/(2*eps)
  ind.pu.rel.grad <- ind.pu.grad/hyp.mat[, "UR"]
  ind.pr.grad     <- (hyp.mat[, "Rplus"] - hyp.mat[, "Rminus"])/(2*eps)
  ind.pr.rel.grad <- ind.pr.grad/p.r
  out <- cbind(Single = hyp.mat[, "UR"], SingleGrad = ind.pu.grad, SingleRelGrad = ind.pu.rel.grad,
               DropOne = hyp.mat[, "R"], DropOneGrad = ind.pr.grad, DropOneRelGrad = ind.pr.rel.grad)

  attr(out, "hypothesis") <- unlist(lapply(hyp.list, attr, "hypothesis"))
  attr(out, "isRatio") <- unlist(lapply(hyp.list, attr, "isRatio"))
  attr(out, "full.restriction.pval") <- p.r
  attr(out, "active") <- sort(rmod$iact)

  return(out)
}

#' Obtain forecasts for an ECM(1)
#'
#' @param x A model or a vector of coefficients
#' @param newdata A data frame or a matrix that contains all the variables corresponding to the coefficient names.
#' Note that the dependent variable should be NA for the period for which the forecast is requested (because the forecasting starts from the last observed value).
#' @param yname A string with the dependent variable name.
#' @param resid.name A string with the name for the last term of the RHS (the residual with a unit coefficient on it)
#' @param d.prefix A string that should precede every differenced variables (usually d_)
#' @param lag.prefix A string that should to precede first lags.
#' @param linear.pred.name Name of the predicted right-hand-side variable to write
#' @param return.full.df Return the full data frame, or just the predicted levels and differences?
#'
#' @return A list of forecast differences and levels of the same length as the number of missing dependent variable values in 'newdata'.
#' @export
#'
#' @examples
predictECM <- function(x, newdata, yname,
                       d.prefix = "d_", lag.prefix = "lag1_",
                      return.full.df = FALSE,
                      linear.pred.name = "Xb",
                      resid.name = "resid"
                      ) {
  if ("restriktor" %in% class(x)) x <- restriktor::coef.restriktor(x) else
    if (!is.numeric(x)) x <- tryCatch(stats::coef(x), error = function(e) stop("'x' is not a numeric, and the coefficients could not be extracted from it via stats::coef(x). Check the input type (numeric or model)."))
  if (is.null(names(x))) stop("No names found in the coefficient vector. Make sure that the estimates are named.")
  if (!(is.matrix(newdata) | is.data.frame(newdata))) stop("'newdata' should be a matrix or a data.frame.")
  if (is.data.frame(newdata)) {
    good.classes <- c("numeric", "integer", "factor", "ordered")
    # Dropping non-numeric variables
    these.classes <- sapply(newdata, class)
    good.vars <- colnames(newdata)[these.classes %in% good.classes]
    if (length(good.vars) != ncol(newdata)) {
      message(paste0("Dropping these non-numeric variables from the supplied data frame: ", paste0(setdiff(colnames(newdata), good.vars), collapse = ", ")))
    }
    newdata <- as.matrix(newdata[, good.vars])
  }
  if (is.null(colnames(newdata))) stop("'newdata' should have column names.")
  if (is.null(rownames(newdata))) warning("'newdata' does not have row names, but it should (for the ease of analysis). Consider adding them.")
  if (!(yname %in% colnames(newdata))) stop("'newdata' should contain the dependent variable as well (because functions of the dep. var. are on the right-hand side).")

  nx <- names(x)
  names.lhs <- paste0(d.prefix, yname)
  names.ly  <- paste0(lag.prefix, yname) # The lagged level that will have to be predicted
  names.dly  <- paste0(c(d.prefix, lag.prefix), c(lag.prefix, d.prefix), yname)
  names.dly <- intersect(names.dly, nx) # Allowing any order of names
  if (length(names.dly) > 1) stop(paste0("More than 1 RHS variable look like the AR(1) term: ", paste0(names.dly, collapse = ", "), ". Check the names."))
  # names.lx <- setdiff(nx[grep(paste0("^", lag.prefix), nx)], names.ly)
  if (!(names.lhs %in% colnames(newdata))) {
    newdata <- cbind(newdata, myDiff(newdata[, yname]))
    colnames(newdata)[ncol(newdata)] <- names.lhs
    newdata <- newdata[-1, ]
  }

  if (length(names.dly) == 0) { # Model without the AR1 term, somehow
    ar1name <- paste0(c(d.prefix, lag.prefix, yname), collapse = "")
    x <- c(x, 0)
    names(x)[length(names(x))] <- ar1name
    nx <- names(x)
    warning("The model has no AR(1) term. Adding it (with zero value) for conformance.")
  }
  y.miss <- which(is.na(newdata[, yname]))
  first.miss <- y.miss[1]
  if (first.miss == 1) stop("The first observation of the level, ", yname, ", must be non-missing.")
  first.miss.d <- which(is.na(newdata[, names.lhs]))[1]
  if (is.null(first.miss)) stop(paste0("The variable of interest, ", yname, ", has no missing values. Out of the estimation sample, it should take NA values."))
  if (is.null(first.miss.d)) {
    warning("The difference of ", yname, ", ", names.lhs, ", has no missing values, but ", yname, " does. Discarding the predicted differences to match the missing levels.")
    newdata[y.miss, names.lhs] <- NA
  } else if (first.miss.d < first.miss) { # There are more missing differences than observed levels
    newdata[first.miss.d:first.miss, names.lhs] <- diff(newdata[(first.miss.d-1):first.miss, yname])
  } else if (first.miss.d > first.miss) {
    warning("The difference of ", yname, ", ", names.lhs, ", has fewer missing values than the level ", yname, ". Discarding the predicted differences to match the missing levels.")
    newdata[y.miss, names.lhs] <- NA
  }

  if (("(Intercept)" %in% nx) & !("(Intercept)" %in% colnames(newdata))) {
    newdata <- cbind(newdata, 1)
    colnames(newdata)[ncol(newdata)] <- "(Intercept)"
  }

  miss.inds <- first.miss:nrow(newdata)
  for (i in miss.inds) {
    newdata[i, names.dly] <- newdata[i-1, names.lhs] # Lagged difference
    newdata[i, names.ly]  <- newdata[i-1, yname] # Lagged levels just in case
    # All exogenous variables should be there by now!
    dyhat <- sum(newdata[i, nx] * x)
    yhat <- newdata[i-1, yname] + dyhat # Y[t] = Y[t-1] + DeltaY[t]
    newdata[i, yname] <- yhat
    newdata[i, names.lhs] <- dyhat
  }
  Xb <- as.numeric(newdata[, nx] %*% x)
  resid <- newdata[, names.lhs] - Xb
  newdata <- cbind(newdata, Xb, resid)
  colnames(newdata)[ncol(newdata) - 1:0] <- c(linear.pred.name, resid.name)

  if (return.full.df) return(as.data.frame(newdata)) else {
    dif <- newdata[miss.inds, names.lhs]
    lev <- newdata[miss.inds, yname]
    names(dif) <- names(lev) <- rownames(newdata)[miss.inds]
    return(list(difference = dif, level = lev))
  }
}

#' Extract variable names automatically from an ECM
#'
#' @param mod An lm object
#' @param d.prefix Character: a regex for identifying differenced variables.
#' @param lag.prefix Character: a regex for identifying lagged variables.
#'
#' @return A named list with separated ECM components (named vectors).
#' @export
#'
#' @examples
guessECM <- function(mod, d.prefix = "^d_", lag.prefix = "^lag1?_") {
  x <- all.vars(stats::formula(mod))
  dy <- x[1] # Dependend variable: Delta y[t]
  x <- setdiff(x, dy)
  ly <- gsub(d.prefix, "", dy)
  ar.regex <- paste0(c("^(", gsub("^\\^", "", c(d.prefix, lag.prefix)), ")|(", gsub("^\\^", "", c(lag.prefix, d.prefix)), ")", ly), collapse = "")
  dly <- x[grep(ar.regex, x)]
  if (length(dly) > 0) x <- setdiff(x, dly)
  if (length(dly) > 1) stop("At this moment, only ECM(1) is supported.")
  lagy <- x[grep(paste0(lag.prefix, ly), x)]
  if (length(lagy) != 1) stop(paste0("Not a proper ECM(1). Lags of Y on the RHS: ", length(lagy), ", expected 1."))
  x <- setdiff(x, lagy)

  # Common vars enter both the LR and SR part
  # They should appear twice in the list of all variables
  # Stripping the suffixes and comparing
  has.d <- grepl(d.prefix, x)
  has.l <- grepl(lag.prefix, x)
  strip.d <- gsub(d.prefix, "", x[has.d])
  strip.l <- gsub(lag.prefix, "", x[has.l])
  both <- c(strip.d, strip.l)
  tb <- table(both)
  common.vars <- names(tb)[tb == 2]

  # Filtering the common ones first for better ordering
  sr   <- x[grep(paste0(c(d.prefix, "(", paste0(common.vars, collapse = "|"), ")"), collapse = ""), x)] # Short-run effects
  lr   <- x[grep(paste0(c(lag.prefix, "(", paste0(common.vars, collapse = "|"), ")"), collapse = ""), x)] # Short-run effects

  x <- setdiff(x, c(sr, lr))

  # The extra variables may have a short-run or a long-run prefix in unrestricted VAR
  # (unrestricted in the SR equation, seed Pesaran (2001) for a discussion)
  sr.rest <- x[grep(d.prefix, x)]
  lr.rest <- x[grep(lag.prefix, x)]
  x <- setdiff(x, c(sr.rest, lr.rest))
  # Other variables that do not fall into any category remain in x and should
  # be added to the short run

  ret <- list(depvar = dy, deplevel = ly, AR = dly, short = c(sr, sr.rest, x), long = c(lr, lr.rest), EC = lagy, common = common.vars)

  return(ret)
}

# TODO: create a class for ECMs for easier prediction and estimation
# So far, there are too many functions not using any meta-information
