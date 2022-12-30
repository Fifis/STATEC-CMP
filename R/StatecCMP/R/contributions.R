#' Plot the positive and negative contributions above and below zero
#'
#' @param x A matrix that will be passed to \code{barplot}.
#' @param main Plot title (used separately to prevent title duplication).
#' @param xlim Passed to barplot.
#' @param ylim Passed to barplot.
#' @param width Passed to barplot.
#' @param space Passed to barplot.
#' @param panel.first Passed to plot.
#' @param ... Passed to barplot.
#'
#' @return Invisibly returns a list of two bar plots (for positive and negative values).
#' @export
#'
#' @examples
#' set.seed(1)
#' negBarPlot(matrix(rnorm(24), ncol = 6), col = rainbow(4, end = 0.7))
negBarPlot <- function(x, xlim = NULL, ylim = NULL,
                       main = "Contribuiton breakdown", width = 1, space = 0.2,
                       panel.first = NULL, ...) {
  pos <- neg <- x
  pos[pos < 0] <- 0
  neg[neg > 0] <- 0
  if (is.null(xlim)) xlim <- c(0, ncol(x)*(width + space))
  if (is.null(ylim)) ylim <- c(min(colSums(neg), na.rm = TRUE), max(colSums(pos), na.rm = TRUE))
  plot(NULL, NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = main, bty = "n", xaxt = "n", yaxt = "n", panel.first = if (is.function(panel.first)) panel.first() else NULL)
  b1 <- graphics::barplot(pos, add = TRUE, xlim = xlim, ylim = ylim, width = width, space = space, ...)
  b2 <- graphics::barplot(neg, add = TRUE, xlim = xlim, ylim = rev(ylim), width = width, space = space, ...)
  graphics::abline(h = 0, lwd = 2, lty = 2)
  return(invisible(b1)) # No need to return both plots since the x coordinates are identical
}

#' Compute contributions for a linear model
#'
#' @param model An lm object for which to compute the contributions
#' @param coef If there is no model, a named vector of coefficients obtained from any source. Requires non-empty 'newdata' with matching column names
#' @param newdata A new data frame for which to compute the contributions; by default, the training data set is used
#' @param dep.var.name If model is NULL and coef is not NULL, use this variable from newdata as the observed values
#' @param resid.name A character for residual label (cosmetic)
#' @param groups A list of character vectors determining which groups of regressors to add up
#' @param zero.residual Compute the contributions with or without the residual?
#' @param warn.on.na.depvar If TRUE and the RHS variable has NAs, wars about the zero-residual assumption.
#'
#' @return A list of two matrices: contributions to the levels and differences.
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- as.numeric(arima.sim(list(ar = c(0.4, 0.3)), n = 100))
#' x2 <- as.numeric(arima.sim(list(ar = c(0.4, 0.3)), n = 100))
#' y <- 1 + x1 + x2 + as.numeric(arima.sim(list(ma = 0.7), n = 100))
#' d <- data.frame(y, x1, x2)
#' mod <- lm(y ~ x1 + x2, data = d)
#' ctb <- computeContribLM(mod)
#' round(tail(ctb$level, 5), 3) # Contributions to Y according to the model
#' true.coef <- c(`(Intercept)` = 1, x1 = 1, x2 = 1)
#' ctb.true <- computeContribLM(coef = true.coef, newdata = d, dep.var.name = "y")
#' round(tail(ctb.true$level, 5), 3) # Real contributions to Y
#' plotContribLM(ctb.true, mar = c(0.2, 2, 2, 0.2), main = "Decomposition",
#'               left.bar.skip = 10, mtext.line = -1)
computeContribLM <- function(model = NULL,
                             coef = NULL,
                             newdata = NULL,
                             dep.var.name = NULL,
                             resid.name = "resid",
                             groups = NULL,
                             zero.residual = FALSE,
                             warn.on.na.depvar = FALSE
) {
  if (is.null(coef) & is.null(model)) stop("Supply either a model or a vector of coefficients + newdata.")

  # We need: model matrix xm, fitted values yhat, residuals res
  if (!is.null(model)) { # An R object was passed
    if (!("lm" %in% class(model))) stop("'data' must be an lm object.")
    coef <- stats::coef(model)
    if (is.null(newdata)) { # Working with the in-sample data
      xm <- stats::model.matrix(model)
      yhat <- model$fitted.values
      res <- model$residuals
      y <- yhat + res
    } else { # Working with a new data set
      # Creating a data set with the dependent variable first and all regressors next
      yx <- stats::get_all_vars(formula = stats::formula(model), data = newdata)
      y <- yx[, 1]
      # Putting any reasonable value where y is missing
      yx[!is.finite(y), 1] <- 0
      xm <- stats::model.matrix(stats::formula(model), data = yx)
      yhat <- stats::predict(model, newdata = yx)
      res <- y - yhat
    }
  } else { # A vector of coefficients was passed
    if (!is.numeric(coef)) stop("'coef' must be numeric.")
    if (is.null(names(coef))) stop("The elements of 'coef' must be named.")
    if (is.null(newdata)) stop("Since 'coef' is numeric and not a model, you must pass a data set.")
    nx <- names(coef)
    nx0 <- setdiff(nx, "(Intercept)") # Meaningful variables that must be in the data frame
    nd <- colnames(newdata)
    if (!all(nx0 %in% nd)) stop(paste0("Not all variables from 'coef' are in 'newdata': ", paste0(setdiff(nx0, nd), collapse = ", ")))

    # Ensuring that the intercept goes first if it is there
    const.ind <- which(nx == "(Intercept)")
    if (length(const.ind) == 1) {
      coef <- c(coef["(Intercept)"], coef[names(coef) != "(Intercept)"])
      nx <- names(coef)
      nx0 <- setdiff(nx, "(Intercept)")
      xm <- cbind(`(Intercept)` = 1, newdata[, nx0])
    } else {
      xm <- newdata[, nx0]
    }
    xm <- as.matrix(xm)
    yhat <- drop(xm %*% coef)

    # Computing residuals if it is possible
    if (!is.null(dep.var.name)) {
      if (!(dep.var.name %in% nd)) stop(paste0("The requested dependent variable ", dep.var.name, "cannot be found in 'newdata'."))
      y <- newdata[, dep.var.name]
      res <- y - yhat
    } else {
      if (!zero.residual) {
        warning("Contributions from residuals were requested (zero.residual = FALSE), but 'dep.var.name' was not supplied (to compute them from the coeff.). Assuming zero residuals = computing contrib. to the predicted dep. var.")
        zero.residual <- TRUE
      }
      y <- yhat
      res <- rep(0, length(y)) # Cannot compute meaningful residuals
    }
  }

  # Diagnosing the RHS matrix
  if (any(!is.finite(xm))) {
    bad.cols <- colnames(xm)[which(apply(xm, 2, function(y) any(!is.finite(y))))]
    bad.rows <- which(apply(xm, 1, function(y) any(!is.finite(y))))
    stop(paste0("The new data contain missing values for the RHS regressors.",
                if (!is.null(model)) "\nProblematic model: " else NULL,
                if (!is.null(model)) Reduce(paste, deparse(stats::formula(model))) else NULL,
                "\nProblematic variables: ", paste0(bad.cols, collapse = ", "),
                ".\nProblematic rows: ", paste0(bad.rows, collapse = ", "),
                ".\nCannot predict without any assumptions about the regressor values.\nSuggestions: extrapolate, put zeros etc."))
  }

  # Diagnosing the LHS vector
  if (any(!is.finite(y))) {
    bad.inds <- which(!is.finite(y))
    if (warn.on.na.depvar) warning("The new data contain missing values for the response variable.\nThese contributions will have no residual component.")
    y[bad.inds] <- yhat[bad.inds]
    res[bad.inds] <- 0
  }

  # cat("dim(xm) = ", dim(xm), ", length(res) = ", length(res), "\n")
  x <- cbind(xm, if (!zero.residual) res else 0) # Adding residuals
  colnames(x)[ncol(x)] <- resid.name

  # Decomposing the levels into the sum of components
  xc <- sweep(x, 2, c(coef, 1), "*")

  # Combining variables by groups
  if (!is.null(groups)) {
    gl <- sapply(groups, length)
    gn <- names(groups)
    # if (any(gl < 2)) warning("Some groups contain fewer than two names. Renaming")
    if (is.null(names(groups))) stop("The groups must be a named list, otherwise the reader will get lost.")
    orig.names <- names(coef(model))
    bad.vars <- setdiff(unlist(groups), orig.names)
    if (length(bad.vars) > 0) stop(paste0("Some grouped variables are not found in the original data: ", paste0(bad.vars, collapse = ", ")))
    first.indices <- drop.indices <- NULL
    for (i in 1:length(groups)) {
      # Replacing the first variable with the group aggregate, and dropping the rest
      wh <- match(groups[[i]], orig.names)
      xc[, wh[1]] <- unname(rowSums(xc[, wh, drop = FALSE]))
      first.indices <- c(first.indices, wh[1])
      gl <- length(groups[[i]])
      if (gl > 1) drop.indices <- c(drop.indices, wh[-1]) # No need to drop groups of length 1
    }
    colnames(xc)[first.indices] <- gn
    xc <- xc[, -drop.indices]
  }
  dxc <- myDiff(xc)

  return(list(level = xc, diff = dxc))
}

#' Compute contributions averaged across the models (with robust options)
#'
#' @param model.list A list of lm objects
#' @param coef.list A list of numerical vector coefficients
#' @param dep.var.name Passed to computeContribLM.
#' @param weights Character ("mean" or "median") or numeric of length(model.list) to use for weighting.
#' @param weight.sorted If TRUE, \code{weights} are applied to the \strong{sorted} predicted levels.
#' This adds a certain degree of robustness: \code{weights = c(0, 0.5, 1, 0.5, 0)} discards
#' the outermost predicted values and uses the middle three to compute the contributions.
#' If FALSE, then, the weights are applied to the models in their original order.
#' @param newdata Passed to computeContribLM.
#' @param resid.name Passed to computeContribLM.
#' @param zero.residual Passed to computeContribLM.
#' @param groups A list of character vectors determining which groups of regressors to add up.
#' Not passed to computeContribLM, but instead, used at the end.
#'
#' @return A list of weighted contributions to the levels, to differences, and vectors of weighted level and diference forecasts.
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- rnorm(100) + 1:100/50 + 1
#' x2 <- rnorm(100) - 1:100/50 - 1
#' x3 <- rnorm(100)
#' y <- 1 + x1 + x2 + 0.1*x3 + rnorm(100)
#' d <- data.frame(y, x1, x2, x3)
#' mod1 <- lm(y ~ x1 + x2, data = d)
#' mod2 <- lm(y ~ x1 + x3, data = d)
#' # Adding a model with strange estimates
#' mod3 <- mod2; mod3$coefficients[1:3] <- c(-1.7, -0.3, 0.8)
#' mod.list <- list(mod1, mod2, mod3)
#' d2 <- d[71:100, ]
#' # Scheme 0: average of 2 good models with 1 bad
#' ctb <- computeContribManyLM(mod.list, newdata = d2)
#' # Scheme 1: robust weights (the central prediction is weighted 10 times more)
#' ctb.w1 <- computeContribManyLM(mod.list, newdata = d2, weights = c(0.1, 1, 0.1))
#' # Scheme 2: always median prediction
#' ctb.w2 <- computeContribManyLM(mod.list, newdata = d2, weights = "median")
#' # Weights: mod1 = 1, mod2 = 1, mod3 (dubious) = 0.1
#' ctb.w3 <- computeContribManyLM(mod.list, newdata = d2, weights = c(1, 1, 0.1),
#'                               weight.sorted = FALSE)
#' round(tail(ctb$level, 5), 3) # Avg. contrib. to Y according to all models
#' # In theory, the contributions of x1 should be positive over time, and of x2, negative
#' plotContribLM(ctb, mar = c(2, 2, 2, 0.2), plot.diff = FALSE)
#' plotContribLM(ctb.w1, mar = c(2, 2, 2, 0.2), plot.diff = FALSE)
#' plotContribLM(ctb.w2, mar = c(2, 2, 2, 0.2), plot.diff = FALSE)
#' plotContribLM(ctb.w3, mar = c(2, 2, 2, 0.2), plot.diff = FALSE)
computeContribManyLM <- function(model.list = NULL, coef.list = NULL,
                                 dep.var.name = NULL,
                                 weights = c("mean", "median"), weight.sorted = TRUE,
                                 newdata = NULL, resid.name = "resid", zero.residual = FALSE,
                                 groups = NULL
) {
  if (is.null(coef.list) & is.null(model.list)) stop("Supply either a model list or a coefficient vector list + newdata.")

  if (is.null(coef.list)) { # Processing models
    m <- length(model.list) # How many models
    contrib.list <- lapply(model.list, function(x) computeContribLM(model = x, newdata = newdata, dep.var.name = dep.var.name, resid.name = resid.name, groups = NULL, zero.residual = zero.residual))
  } else {
    m <- length(coef.list)
    contrib.list <- lapply(coef.list, function(x) computeContribLM(coef = x, newdata = newdata, dep.var.name = dep.var.name, resid.name = resid.name, groups = NULL, zero.residual = zero.residual))
  }

  n <- nrow(contrib.list[[1]]$level) # How many observations
  all.var.vec  <- unique(unlist(lapply(contrib.list, function(x) colnames(x[[1]]))))
  k <- length(all.var.vec) # How many regressors in total
  contrib.arr  <- dcontrib.arr <- array(0, dim = c(n, k, m), dimnames = list(rownames(newdata), all.var.vec, paste0("Model", 1:m)))
  # Filling in the existing variables
  for (i in 1:m) {
    if (dim(contrib.arr)[1] != nrow(contrib.list[[i]][["level"]])) stop(paste0("The contributions for model ", i, " were not for some time periods. Are you using any variables with future missing values (e.g. lags of dep. var.)? If so, predict them first."))
    contrib.arr[, colnames(contrib.list[[i]][["level"]]), i]  <- contrib.list[[i]][["level"]]
    dcontrib.arr[, colnames(contrib.list[[i]][["diff"]]), i] <- contrib.list[[i]][["diff"]]
  }

  # Note that the trimmed mean of a sum is not equal to the sum of trimmed means! We need to be careful.
  # We compute the numbers of sorted forecasts
  if (is.character(weights)) {
    weights <- weights[1]
    if (weights == "median") {
      middle <- if (m %% 2 == 0) m/2+0:1 else ceiling(m/2)
      weights <- rep(0, m)
      weights[middle] <- 1
    } else if (weights == "mean") {
      weights <- rep(1, m)
    } else stop("weights should be 'mean', 'median', or numeric")
  }
  y.fc <- do.call(cbind, lapply(contrib.list, function(x) rowSums(x[["level"]])))
  dy.fc <- do.call(cbind, lapply(contrib.list, function(x) rowSums(x[["diff"]])))
  y.order <- if (weight.sorted) t(apply(y.fc, 1, order)) else matrix(rep(1:m, each = n), nrow = n)
  # For every step, we want to be able to get any weighted statistic with the weights
  # defined by the order at time i
  orderedStat <- function(x, i, FUN = mean, ...) {
    obj.order <- y.order[i, ]
    x <- x[obj.order]
    # print(x)
    FUN(x, ...)
  }
  y.wfc  <- sapply(1:n, function(i) orderedStat(x = y.fc[i, ],  i = i, FUN = stats::weighted.mean, w = weights))
  dy.wfc <- sapply(1:n, function(i) orderedStat(x = dy.fc[i, ], i = i, FUN = stats::weighted.mean, w = weights))
  wcontrib <-  wdcontrib <- matrix(0, nrow = n, ncol = k)
  dimnames(wcontrib) <- dimnames(wdcontrib) <- list(rownames(newdata), all.var.vec)
  for (i in 1:n) {
    wcontrib[i, ]  <- sapply(1:k, function(j) orderedStat(x = contrib.arr[i, j, ],  i = i, FUN = stats::weighted.mean, w = weights))
    wdcontrib[i, ] <- sapply(1:k, function(j) orderedStat(x = dcontrib.arr[i, j, ], i = i, FUN = stats::weighted.mean, w = weights))
  }

  # Combining variables by groups
  if (!is.null(groups)) {
    gl <- sapply(groups, length)
    gn <- names(groups)
    # if (any(gl < 2)) warning("Some groups contain fewer than two names.")
    if (is.null(names(groups))) stop("The groups must be a named list, otherwise the reader will get lost.")
    orig.names <- all.var.vec
    bad.vars <- setdiff(unlist(groups), orig.names)
    if (length(bad.vars) > 0) stop(paste0("Some grouped variables are not found in the original data: ", paste0(bad.vars, collapse = ", ")))
    first.indices <- drop.indices <- NULL
    for (i in 1:length(groups)) {
      # Replacing the first variable with the group aggregate, and dropping the rest
      wh <- match(groups[[i]], orig.names)
      wcontrib[, wh[1]] <- unname(rowSums(wcontrib[, wh, drop = FALSE]))
      wdcontrib[, wh[1]] <- unname(rowSums(wdcontrib[, wh, drop = FALSE]))
      first.indices <- c(first.indices, wh[1])
      gl <- length(groups[[i]])
      if (gl > 1) drop.indices <- c(drop.indices, wh[-1]) # No need to drop groups of length 1
    }
    colnames(wcontrib)[first.indices] <- colnames(wdcontrib)[first.indices] <- gn
    wcontrib <- wcontrib[, -drop.indices]
    wdcontrib <- wdcontrib[, -drop.indices]
  }
  return(list(level = wcontrib, diff = wdcontrib, yhat = y.wfc, dyhat = dy.wfc))
}


#' Visualise contributions from many linear models
#'
#' @param x A list returned by computeContribManyLM
#' @param name.replacement.mat Character matrix with two columns containing the existing and the new names in its columns
#' @param contrib.threshold Numeric: Do not plot the variables whose max. relative contribution to the levels was less than this
#' @param threshold.type Character: Drop small contributions based on minimal fraction or absolute threshold?
#' @param small.contrib Character: What to do with the contributions below the threshold: ignore, add up, or add up and check if they exceed the threshold, or do nothing
#' @param small.name Character: label for all aggregated insignificant contributions
#' @param francais If TRUE, the legends are in French
#' @param plot.diff If TRUE, adds a bottom plot with the changes
#' @param col Character vector Colours of the vars
#' @param left.bar.skip How many bars of space to skip on the left (to help the legend fit)
#' @param mar Plot margins, par("mar"); change if the x-axis labels won't fit
#' @param mtext.line A line under the level plot at which the note is placed.
#' @param ... Passed to negBarPlot.
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- rnorm(100) + 1:100/50 + 1
#' x2 <- rnorm(100) - 1:100/50 - 1
#' x3 <- rnorm(100)
#' y <- 1 + x1 + x2 + 0.5*x3 + rnorm(100)
#' d <- data.frame(y, x1, x2, x3)
#' mod1 <- lm(y ~ x1 + x2, data = d)
#' mod2 <- lm(y ~ x1 + x3, data = d)
#' d2 <- d[71:100, ]
#' ctb.one <- computeContribLM(mod1, newdata = d2)
#' ctb.many <- computeContribManyLM(list(mod1, mod2), newdata = d2)
#' plotContribLM(ctb.one)
#' plotContribLM(ctb.many, mtext.line = 2.5)
plotContribLM <- function(x, name.replacement.mat = NULL,
                          contrib.threshold = 0.03, threshold.type = c("relative", "absolute"),
                          small.contrib = c("drop", "mergekeep", "mergedrop", "preserve"),
                          small.name = "Autres", francais = FALSE, plot.diff = TRUE, col = NULL,
                          left.bar.skip = 3, mar = c(5, 2.2, 2.2, 0.2), mtext.line = -1, ...
) {
  small.contrib <- small.contrib[1]
  threshold.type <- threshold.type[1]

  if (!is.null(name.replacement.mat)) {
    nm <- colnames(x$level)
    m <- match(nm, name.replacement.mat[, 1])
    if (any(is.finite(m))) nm[!is.na(m)] <-  name.replacement.mat[stats::na.omit(m), 2]
    colnames(x$level) <- colnames(x$diff) <- nm
  }


  lev <- x$level
  dif <- x$diff
  levsum <- rowSums(lev)
  lev.rel <- switch(threshold.type, relative = abs(lev) / abs(levsum),
                    absolute = abs(lev))
  if (is.null(lev.rel)) stop("threshold.type should be 'relative' or 'absolute'.")

  small.vars <- which(apply(lev.rel, 2, max) < contrib.threshold)

  if (francais) {
    mtext.nodrop      <- "Montr\u00E9 : toutes les variables"
    mtext.merge     <- paste0("Les variables dont la contrib. absolue max. < ", if (threshold.type == "relative") round(contrib.threshold*100) else contrib.threshold, if (threshold.type == "relative") " % du total" else "", " agr\u00E9g\u00E9es comme \u00AB ", small.name, " \u00BB")
    mtext.drop      <- paste0("Montr\u00E9 : uniquement des variables dont la contrib. absolue max. est > ", if (threshold.type == "relative") round(contrib.threshold*100) else contrib.threshold, if (threshold.type == "relative") " % du total" else "")
  } else {
    mtext.nodrop      <- "Shown: all the variables"
    mtext.merge     <- paste0("The variables with max. abs. contrib. < ", if (threshold.type == "relative") round(contrib.threshold*100) else contrib.threshold, if (threshold.type == "relative") " % of total" else "", " combined as '", small.name, "'")
    mtext.drop      <- paste0("Shown: only the variables with max. abs. contrib. > ", if (threshold.type == "relative") round(contrib.threshold*100) else contrib.threshold, if (threshold.type == "relative") " % of total" else "")
  }

  mtext.this <- mtext.nodrop
  if (small.contrib == "drop") {
    if (length(small.vars) > 0) {
      mtext.this <- mtext.drop
      lev <- lev[, -small.vars, drop = FALSE]
      dif <- dif[, -small.vars, drop = FALSE]
      lev <- lev / rowSums(lev) * levsum # Normalising back
    }
  } else if (small.contrib == "mergekeep") {
    if (length(small.vars) > 1) { # Changes are needed only if there is something to merge
      small.lev.sum <- rowSums(lev[, small.vars])
      small.dif.sum <- rowSums(dif[, small.vars])
      lev <- cbind(lev[, -small.vars, drop = FALSE], small.lev.sum)
      dif <- cbind(dif[, -small.vars, drop = FALSE], small.dif.sum)
      colnames(lev)[ncol(lev)] <- colnames(dif)[ncol(dif)] <- small.name
      mtext.this <- mtext.merge
    }
  } else if (small.contrib == "mergedrop") {
    if (length(small.vars) == 1) { # If there is one small variable, drop it
      lev <- lev[, -small.vars, drop = FALSE]
      dif <- dif[, -small.vars, drop = FALSE]
      lev <- lev / rowSums(lev) * levsum # Normalising back
      mtext.this <- mtext.drop
    } else if (length(small.vars > 1)) { # If there is more than 1 small variable: merge, test, and maybe drop
      small.lev.sum <- rowSums(lev[, small.vars])
      small.dif.sum <- rowSums(dif[, small.vars])
      lev <- cbind(lev[, -small.vars, drop = FALSE], small.lev.sum)
      dif <- cbind(dif[, -small.vars, drop = FALSE], small.dif.sum)
      colnames(lev)[ncol(lev)] <- colnames(dif)[ncol(dif)] <- small.name
      mtext.this <- mtext.merge
      lev.rel <- switch(threshold.type, relative = abs(lev) / abs(levsum), # levsum is the same as before the merging
                        absolute = abs(lev))
      small.vars <- which(apply(lev.rel, 2, max) < contrib.threshold)
      if (length(small.vars) > 0) {
        lev <- lev[, -small.vars, drop = FALSE]
        dif <- dif[, -small.vars, drop = FALSE]
        lev <- lev / rowSums(lev) * levsum # Normalising back
        mtext.this <- mtext.drop
      }
    }
  } else if (small.contrib != "preserve") stop('small.contrib should be "drop", "mergekeep", "mergedrop", or "preserve".')

  if (plot.diff) withr::local_par(mar = mar, mfrow = c(2, 1)) else withr::local_par(mar = mar)
  if (is.null(col)) col <- grDevices::rainbow(ncol(lev), end = 0.8, v = 0.8)

  arg.list <- list(...)
  arg.list$xlim <- c(-left.bar.skip*1.2, nrow(lev)*1.2)
  arg.list$col <- col
  arg.list$x <- t(lev)
  b1 <- do.call(negBarPlot, arg.list)
  lsum <- rowSums(lev)
  graphics::lines(b1, lsum, lwd = 4, col = "white")
  graphics::lines(b1, lsum, lwd = 2)
  graphics::legend("topleft", colnames(lev), col = col, pch = 15, pt.cex = 1.5, bg = "#FFFFFFBB")
  # mtext("Showing only the variables whose max. abs. contribution was > 2% of the total", 1, line = 4.5, cex = 0.8)

  graphics::mtext(mtext.this, 1, line = mtext.line, cex = 0.8)

  if (plot.diff) {
    arg.list$x <- t(dif)
    arg.list$main <- if (francais) "Diff\u00E9rences des contributions ci-dessus" else "Differences of the above contributions"
    b2 <- do.call(negBarPlot, arg.list)
    graphics::legend("bottomleft", colnames(dif), col = col, pch = 15, pt.cex = 1.5, bg = "#FFFFFFBB")
  }

  return(invisible(NULL))
}


#' Compare two sets of contributions from multple models
#'
#' @param old A list returned by computeContribManyLM
#' @param new A list returned by computeContribManyLM
#'
#' @return A list with the differences of contributions
#' @export
#'
#' @examples
#' set.seed(1)
#' x1 <- rnorm(100) + 1:100/50 + 1
#' x2 <- rnorm(100) - 1:100/50 - 1
#' x3 <- rnorm(100)
#' y <- 1 + x1 + x2 + 0.5*x3 + rnorm(100)
#' d <- data.frame(y, x1, x2, x3)
#' mod1 <- lm(y ~ x1 + x2, data = d)
#' mod2 <- lm(y ~ x1 + x3, data = d)
#' mod3 <- lm(y ~ x1 + x2 + x3, data = d)
#' d2 <- d[71:100, ]
#' ctb1 <- computeContribManyLM(list(mod1, mod2), newdata = d2)
#' ctb2 <- computeContribManyLM(list(mod1, mod2, mod3), newdata = d2)
#' ctb.diff <- compareContributions(ctb1, ctb2)
#' plotContribLM(ctb.diff, main = "Contribution discrepancies")
compareContributions <- function(old, new) {
  old.vars <- colnames(old$level)
  new.vars <- colnames(new$level)
  if (is.null(old.vars) | is.null(new.vars)) stop("Difference computation impossible due to the lack of variable names.")
  old.periods <- rownames(old$level)
  new.periods <- rownames(new$level)
  if (is.null(old.periods) | is.null(new.periods)) stop("Difference computation impossible because the observations (periods) should have names. Make sure both models have observation names.")
  both.vars <- union(old.vars, new.vars)
  both.periods <- union(old.periods, new.periods)

  # Old and new level and difference arrays
  l1 <- l2 <- d1 <- d2 <- array(0, dim = c(length(both.periods), length(both.vars)), dimnames = list(both.periods, both.vars))
  # Old and new predicted value vectors
  yhat1 <- yhat2 <- dyhat1 <- dyhat2 <- array(0, dim = length(both.periods), dimnames = list(both.periods))
  l1[old.periods, old.vars] <- old$level
  l2[new.periods, new.vars] <- new$level
  d1[old.periods, old.vars] <- old$diff
  d2[new.periods, new.vars] <- new$diff
  yhat1[old.periods] <- old$yhat
  yhat2[new.periods] <- new$yhat
  dyhat1[old.periods] <- old$dyhat
  dyhat2[new.periods] <- new$dyhat

  l0 <- l2 - l1
  d0 <- d2 - d1
  yhat0 <- yhat2 - yhat1
  dyhat0 <- dyhat2 - dyhat1
  out <- list(level = l0, diff = d0, yhat = yhat0, dyhat = dyhat0)
  attr(out, "contribution.type") <- "difference"
  return(out)
}


