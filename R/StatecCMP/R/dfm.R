#' Generate EViews code with reasonable initial values
#'
#' @param signals Character vector of signal names (without first revisions)
#' @param signals.with.1strev Character vector of signal names that have a first revision (FR)
#' @param signals.1strev Character vector of vector names that are first revisions
#' @param sig.names Signal names (without FR) for human-readable output
#' @param sig.w1r.names Signal names with first revision for human-readable output
#' @param mod.name State-space model name for the code
#' @param data Data set to compute the initial values; must contain the column names identical to signal names (signals, signals.with.1strev, signals.1strev).
#' @param nfac Number of unobserved factors
#' @param err.lag Maximum AR lag of the dynamic observation error (a non-negative integer scalar)
#' @param fac.lag Maximum AR lag of the dynamic factor (a non-negative integer scalar)
#' @param rev.lag If TRUE, adds one lag of revision error
#' @param positive.first If TRUE, enforce most positive coefficients on the first principal component
#' @param trim.ar A small positive integer less than 1: if the preliminary estimates of AR coefficients are greater than this number, trim them at this level to ensure stationarity
#'
#' @return A list with 5 elements: EViews model text, tab-separated table for copying into Excel, initial-value vector, the imputed data set, and PCA done on it
#' @export
#'
#' @examples
#' # A simple model with 3 variables and 2 factors
#' genDFM()
genDFM <- function(signals = c("emp_m_g_s", "u_m_g_s", "pib_r_q_g_s"),
                   signals.with.1strev = NULL, signals.1strev = NULL,
                   sig.names = c("Total employment", "Unemployed", "Lux GDP"),
                   sig.w1r.names = NULL,
                   mod.name = "emp",
                   data = NULL,
                   nfac = 2, err.lag = 2, fac.lag = 2, rev.lag = FALSE,
                   positive.first = TRUE, trim.ar = 0.9
) {
  if (length(signals.with.1strev) != length(signals.1strev)) stop("The length of signals with 1st release and their 1st-release counterparts must be the same.")
  if (length(sig.w1r.names) != length(signals.with.1strev)) stop("The length of the name vector for signals with 1st release should be the same as the the length of the signal variables itself.")

  out <- NULL
  p0 <- function(...) paste0(..., collapse = "")
  out <- c(out, p0("sspace ", mod.name), "")

  # Interleaving the signals with their 1st-release counterparts, and then all the other vectors without 1st releases
  lens <- c(rep(length(signals.with.1strev), 2), length(signals)) # Variable type count
  types <- c(rep(c("s.with1r", "s.1r"), lens[1]),  rep("ordinary", lens[3]))
  if (length(signals.with.1strev) > 0) {
    all.signals <- c(as.vector(t(cbind(signals.with.1strev, signals.1strev))), signals)
    all.sig.names <- c(sig.w1r.names, paste0(sig.w1r.names, " 1st rel."), sig.names)
    do.1r <- TRUE
  } else {
    all.signals <- signals
    all.sig.names <- sig.names
    do.1r <- FALSE
  }

  if (nfac > 4) stop("Not more than 4 common factors are allowed.")
  coef.labels <- err.labels <-  seq_along(all.signals)
  # 1XX, 2XX, ..., 4XX --- coefficients on common factors
  # 1 [2] --- AR(1) [AR(2) if it is used] coefficients on the 1st factor
  # [3] [4] --- AR(1), AR(2) for the 2nd factor (if it exists)
  # [5] [6] --- AR(1), AR(2) for the 3rd factor (if it exists)
  # [7] [8] --- AR(1), AR(2) for the 4th factor (if it exists)
  # 1XX + 50 --- coefficients on the variance
  # 5XX + 50--- means

  # Using the same coefficients for the factor effects for 1st-releast variables
  if (do.1r) coef.labels[types == "s.1r"] <- coef.labels[which(types == "s.1r") - 1]

  coef.inds <- sprintf("%02d", coef.labels) # Name endings of coefficient vectors
  var.inds <- sprintf("%02d", err.labels) # Name endings of variance vectors
  err.inds <- sprintf("%02d", err.labels + 50) # Name endings of error vectors

  all.inds <- NULL # Vector with all coefficient indices for initial values

  out <- c(out, "' Observation (signal) equations")
  for (i in seq_along(all.signals)) {
    s <- p0(mod.name, ".append @signal ", all.signals[i], " = c(1", coef.inds[i], ")*SV1_0")
    if (types[i] != "s.1r") {
      all.inds <- c(all.inds, as.numeric(paste0(1, coef.inds[i])))
      names(all.inds)[length(all.inds)] <- paste0("SV1.sig", coef.inds[i])
    }
    if (nfac > 1) {
      for (j in 2:nfac) {
        s <- p0(s, " + c(", j, coef.inds[i], ")*SV", j, "_0")
        if (types[i] != "s.1r") {
          all.inds <- c(all.inds, as.numeric(paste0(j, coef.inds[i])))
          names(all.inds)[length(all.inds)] <- paste0("SV", j, ".sig", coef.inds[i])
        }
      }
    }
    s <- p0(s, " + U", coef.inds[i], "_0")
    if (types[i] == "s.1r") {
      s <- p0(s, " + RVER", coef.inds[i])
    }
    out <- c(out, s)
  }
  out <- c(out, "")

  out <- c(out, "' State equations, p.1: Dynamic factors (auto-regressive)")
  for (i in 1:nfac) {
    s <- p0(mod.name, ".append @state SV", i, "_0 = c(", i, "1)*SV", i, "_0(-1)")
    all.inds <- c(all.inds, as.numeric(paste0(i, 1)))
    names(all.inds)[length(all.inds)] <- paste0("AR1.SV", i)
    if (fac.lag > 1) { # State variable lags must be defined through (-1) lags---a 'feature' of EViews
      for (j in 2:fac.lag) {
        s <- paste0(s, " + c(", i, j, ")*SV", i, "_", j-1, "(-1)")
        all.inds <- c(all.inds, as.numeric(paste0(i, j)))
        names(all.inds)[length(all.inds)] <- paste0("AR", j, ".SV", i)
      }
    }
    s <- p0(s, " + ESV", i)
    out <- c(out, s)
  }
  out <- c(out, "")

  # Defining those (-1) lags as before in a separate block for readability
  if (fac.lag > 1) {
    out <- c(out, "' State equations, p.1a: State variable lags must be defined through (-1) lags---EViews not allowing higher-order lags")
    for (i in 1:nfac) {
      for (j in 2:fac.lag) {
        s <- paste0(mod.name, ".append @state SV", i, "_", j-1, " = SV", i, "_", j-2, "(-1)")
        out <- c(out, s)
      }
    }
    out <- c(out, "")
  }

  # Unit variance restriction
  out <- c(out, "' Named error variable and unit variance restriction")
  for (i in 1:nfac) {
    out <- c(out, p0(mod.name, ".append @ename ESV", i))
    out <- c(out, p0(mod.name, ".append @evar var(ESV", i, ") = 1"))
  }
  out <- c(out, "")

  # The coefficients on lagged errors have an offset of +50
  out <- c(out, "' State equations, p.2: Signal equation errors (auto-regressive)")
  for (i in seq_along(all.signals)) {
    if (types[i] == "s.1r") { # For a first-release variable, the name is different, but the structure is the same
      s <- p0(mod.name, ".append @state RVER", coef.inds[i], " = c(5", err.inds[i], ")")
      all.inds <- c(all.inds, as.numeric(paste0(5, err.inds[i])))
      names(all.inds)[length(all.inds)] <- paste0("mean.REV", coef.inds[i])
      if (rev.lag) {
        s <- p0(s, " + c(1", err.inds[i], ")*RVER", coef.inds[i], "(-1)")
        all.inds <- c(all.inds, as.numeric(paste0(1, err.inds[i])))
        names(all.inds)[length(all.inds)] <- paste0("AR1.REV", coef.inds[i])
      }
    } else {
      s <- p0(mod.name, ".append @state U", coef.inds[i], "_0 = c(1", err.inds[i], ")*U", coef.inds[i], "_0(-1)")
      all.inds <- c(all.inds, as.numeric(paste0(1, err.inds[i])))
      names(all.inds)[length(all.inds)] <- paste0("AR1.err", coef.inds[i])
      if (err.lag > 1) {
        for (j in 2:err.lag) {
          s <- p0(s, " + c(", j, err.inds[i], ")*U", coef.inds[i], "_", j-1, "(-1)")
          all.inds <- c(all.inds, as.numeric(paste0(j, err.inds[i])))
          names(all.inds)[length(all.inds)] <- paste0("AR", j, ".err", coef.inds[i])
        }
      }
    }
    s <- p0(s, " + [var = exp(c(5", var.inds[i], "))]")
    all.inds <- c(all.inds, as.numeric(paste0(5, var.inds[i])))
    names(all.inds)[length(all.inds)] <- if (types[i] != "s.1r") paste0("logvar", coef.inds[i]) else paste0("logvarREV", coef.inds[i])
    out <- c(out, s)
  }
  out <- c(out, "")
  # Defining those (-1) state error lags as before in a separate block for readability
  if (err.lag > 1) {
    out <- c(out, "' State equations, p.2a: State variable lags must be defined through (-1) lags---a 'feature' of EViews not allowing higher-order lags")
    for (i in setdiff(seq_along(all.signals), which(types == "s.1r"))) {
      for (j in 2:err.lag) {
        s <- paste0(mod.name, ".append @state U", coef.inds[i], "_", j-1, " = U", coef.inds[i], "_", j-2, "(-1)")
        out <- c(out, s)
      }
    }
    out <- c(out, "")
  }

  init.vals <- rep(0, length(all.inds))
  names(init.vals) <- names(all.inds)
  if (!is.null(data)) {
    init.vals[grep("^AR1", names(init.vals))] <- 0.5
    if (!is.null(data)) {
      data0 <- data[, all.signals]
      # First, process quarterly data and fill the gaps of size 2
      for (i in seq_len(ncol(data0))) {
        data0[, i] <- imputeTS::na_locf(data0[, i], option = "nocb", maxgap = 2) # Filling quarterly gaps
      }
      data1 <- data0[, c(signals.with.1strev, signals)] # Not using FR data because the coefficient on the factor must be the same
      if (any(is.na(data1))) { # Now, fill the rest
        nc <- missMDA::estim_ncpPCA(data1)
        # plot(as.numeric(names(nc$criterion)), nc$criterion, type = "l")
        imput <- missMDA::imputePCA(data1, ncp = max(as.numeric(names(nc$criterion)[which.min(nc$criterion)]), 3)) # Not fewer than 3 components for imputation
        data1 <- as.data.frame(imput$completeObs)
      }
      pca <- stats::prcomp(data1[stats::complete.cases(data1), ])
      rot <- pca$rotation # PCA rotation matrix
      pcs <- pca$x # Components themselves
      if ((mean(rot[, 1] < 0) > 0.5) & positive.first) {
        rot <- -rot # Getting most positive coefficients on the first PC
        pcs <- -pcs
      }
      pcs <- pcs[, 1:nfac, drop = FALSE]
      rot <- rot[, 1:nfac, drop = FALSE]
      rot <- rot / max(abs(rot)) / 3
      data1 <- cbind(data1, pcs) # Adding the signals for PC regressions
      if (do.1r) { # Adding the first-revision signal
        data1 <- cbind(data1, data[, signals.1strev])
        colnames(data1)[ncol(data1) - 1:lens[2] + 1] <- signals.1strev
      }

      # Initial values for the signal-to-state transformation based on PCA loadings
      # AR(1) coefficients from AR models
      for (i in 1:nfac) {
        init.vals[grep(paste0("SV", i, "\\.sig"), names(init.vals))] <- rot[, i]
        pc.ar <- stats::ar(pcs[, i], aic = FALSE, order.max = fac.lag)
        rho0 <- pc.ar$ar
        if (rho0[1] > trim.ar) rho0[1] <- trim.ar # Safeguarding against non-stationarity
        if (fac.lag > 1) rho0[-1] <- rho0[-1] * 0.25 # Shrinking the rest closer to zero
        init.vals[paste0("AR", 1:fac.lag, ".SV", i)] <- rho0
      }

      # Initial values for the signal variance come from the assumed AR process
      for (i in seq_along(all.signals)) {
        f <- paste0(all.signals[i], " ~ ", paste0("PC", 1:nfac, collapse = " + "))
        signal.resid <- stats::lm(stats::formula(f), data = data1) # The best prediction of the signal given the principal components
        if (types[i] != "s.1r") {
          resid.ar <- stats::ar(signal.resid$residuals, aic = FALSE, order.max = err.lag) # Error AR fit
          trimmed.resid <- resid.ar$resid[is.finite(resid.ar$resid) & resid.ar$resid > stats::quantile(resid.ar$resid, 0.025, na.rm = TRUE) & resid.ar$resid < stats::quantile(resid.ar$resid, 0.975, na.rm = TRUE)]
          res.var <- stats::var(trimmed.resid)
          rho0 <- resid.ar$ar
          if (rho0[1] > trim.ar) rho0[1] <- trim.ar # Safeguarding against non-stationarity
          if (err.lag > 1) rho0[-1] <- rho0[-1] * 0.2 # Shrinking the rest slightly more because errors are less stable than signals
          init.vals[paste0("AR", 1:err.lag, ".err", coef.inds[i])] <- rho0
          init.vals[paste0("logvar", var.inds[i])] <- log(res.var) * 0.75 # Slightly shrinking towards zero
          signal.resid.old <- signal.resid
        } else {
          # Revision error mean, using the existing prediction because 1r always goes after the main signal
          sighat <- stats::predict(signal.resid.old)
          sigdiff <- data1[, all.signals[i]] - sighat
          init.vals[paste0("mean.REV", coef.inds[i])] <- mean(sigdiff, na.rm = TRUE)
          rev.data <- cbind(data1[, all.signals[i]], data1[, c(signals.with.1strev, signals)])
          rev.imp <- missMDA::imputePCA(rev.data, ncp = 2)$completeObs[, 1]
          rev.err <- data1[, all.signals[i-1]] - rev.imp
          if (rev.lag) {
            rev.ar <- stats::ar(rev.err, aic = FALSE, order.max = 1)
            rho0 <- rev.ar$ar
            if (rho0 > trim.ar) rho0 <- trim.ar
            init.vals[paste0("AR1.REV", coef.inds[i])] <- rho0
            trimmed.resid <- rev.ar$resid[is.finite(rev.ar$resid) & rev.ar$resid > stats::quantile(rev.ar$resid, 0.025, na.rm = TRUE) & rev.ar$resid < stats::quantile(rev.ar$resid, 0.975, na.rm = TRUE)]
          } else {
            trimmed.resid <- rev.err[is.finite(rev.err) & rev.err > stats::quantile(rev.err, 0.025, na.rm = TRUE) & rev.err < stats::quantile(rev.err, 0.975, na.rm = TRUE)]
          }
          res.var <- stats::var(trimmed.resid)
          init.vals[paste0("logvarREV", coef.inds[i])] <- log(res.var) * 0.75 # Slightly shrinking towards zero
        }
      }
    }
  }

  nice.names <- nice.groups <- names(all.inds)
  if ((length(sig.names) > 0) & (length(sig.names) != length(signals))) {
    sig.names <- NULL
    warning("The length of the nice labels is not the same as that of the signal vector, omitting it.")
  }
  nice.names <- gsub("^SV([1-4])\\.sig", "CF\\1 eff. on sig", nice.names)
  nice.names <- gsub("^AR([1-4])\\.SV", "AR\\1 of CF", nice.names)
  if (length(all.sig.names) > 0) {
    for (i in seq_along(all.signals)) {
      nice.names <- gsub(paste0("sig", coef.inds[i], "$"), all.sig.names[i], nice.names)
      nice.names <- gsub(paste0("^mean\\.REV", coef.inds[i], "$"), paste0("Mean rev. err. of ", all.sig.names[i]), nice.names)
      nice.names <- gsub(paste0("^AR1\\.REV", coef.inds[i], "$"), paste0("AR1 of rev. err. of ", all.sig.names[i]), nice.names)
      nice.names <- gsub(paste0("^logvarREV", coef.inds[i], "$"), paste0("Log Variance of rev. err. of ", all.sig.names[i]), nice.names)
      nice.names <- gsub(paste0("\\.err", var.inds[i], "$"), paste0(" err. of ", all.sig.names[i]), nice.names)
      nice.names <- gsub(paste0("^logvar", var.inds[i], "$"), paste0("Log Variance of ", all.sig.names[i]), nice.names)
    }
  } else {
    for (i in seq_along(all.signals)) {
      nice.names <- gsub(paste0("mean.REV", coef.inds[i], "$"), paste0("Mean rev. err. of sig", coef.inds[i]), nice.names)
      nice.names <- gsub(paste0("^AR1\\.REV", coef.inds[i], "$"), paste0("AR1 of rev. err. of sig", coef.inds[i]), nice.names)
      nice.names <- gsub(paste0("^logvarREV", coef.inds[i], "$"), paste0("Log Variance of rev. err", coef.inds[i]), nice.names)
      nice.names <- gsub(paste0("\\.err", var.inds[i], "$"), paste0(" err. of sig", var.inds[i]), nice.names)
      nice.names <- gsub(paste0("^logvar", var.inds[i], "$"), paste0("Log Variance of err", var.inds[i]), nice.names)
    }
  }

  nice.groups <- gsub("^SV.+", "Factor effect", nice.groups)
  nice.groups <- gsub("^AR[0-9]+\\.SV.+", "Factor AR", nice.groups)
  nice.groups <- gsub("^AR[0-9]+\\.err.+", "Error AR", nice.groups)
  nice.groups <- gsub("^AR[0-9]+\\.REV.+", "Error AR", nice.groups)
  nice.groups <- gsub("^logvar.+", "Log Variance", nice.groups)
  nice.groups <- gsub("^mean.REV.+", "Revision mean", nice.groups)

  # Adding initial values for comfortable line-by-line view
  line.num <- 1
  init.lines <- names(all.inds)
  init.lines <- gsub("^AR[0-9]+\\.SV.+", line.num, init.lines)
  line.num <- line.num + 1
  for (i in 1:nfac) {
    init.lines <- gsub(paste0("^SV", i, "\\.sig.+"), line.num, init.lines)
    line.num <- line.num + 1
  }
  for (i in 1:err.lag) {
    init.lines <- gsub(paste0("^AR", i, "\\.err.+"), line.num, init.lines)
    line.num <- line.num + 1
  }
  init.lines <- gsub("^logvar.+", line.num, init.lines)
  line.num <- line.num + 1
  init.lines <- gsub(".+\\.REV.+", line.num, init.lines)
  init.lines <- as.numeric(init.lines)

  out <- c(out, p0("' Initial values ", if (is.null(data)) "at zeros (completely agnostic)" else "based on PCA and sample averages"))
  init.df <- data.frame(index = paste0("C(", all.inds, ")"), value = round(init.vals, 3))
  init.df <- split(init.df, f = init.lines)
  init.pasted <- unname(unlist(lapply(init.df, function(x) paste0(as.vector(t(as.matrix(x))), collapse = " "))))
  init.pasted[1:(length(init.pasted)-1)] <- paste0(init.pasted[1:(length(init.pasted)-1)], " _")
  s <- p0(mod.name, ".append @PARAM ", init.pasted[1])
  out <- c(out, s, paste0("        ", init.pasted[-1]))

  # Appending zeros is unnecessary, but in case somebody needs...
  # nstate <- nfac*fac.lag + (lens[1] + lens[3])*err.lag + lens[2] # Number of unobserved state variables
  # out <- c(out, p0("vector(", nstate, ") svec0"), p0(mod.name, ".append @mprior svec0"))

  vars <- data.frame(index = unname(all.inds), code = names(all.inds), name = nice.names, group = nice.groups)

  if (!is.null(data)) {
    nc <- missMDA::estim_ncpPCA(data0)
    # plot(as.numeric(names(nc$criterion)), nc$criterion, type = "l")
    imput <- missMDA::imputePCA(data0, ncp = max(as.numeric(names(nc$criterion)[which.min(nc$criterion)]), 3)) # Not fewer than 3 components for imputation
    data1 <- as.data.frame(imput$completeObs)
    pca <- stats::prcomp(data1)
  } else {
    data1 <- pca <- NULL
  }

  return(list(text = out, labels = vars[order(vars$index), ], init.val = init.vals, data = data1, pca = pca))
}


#' Print the EViews code of a DFM and a tab-separated label list
#'
#' @param x List returned by genDFM
#'
#' @return Nothing.
#' @export
#'
#' @examples
#' m <- genDFM()
#' printDFM(m)
printDFM <- function(x) {
  cat(x$text, sep = "\n")
  cat("\n\n")
  for (i in 1:nrow(x$labels)) cat(paste0(x$labels[i, 1:3], collapse = "\t"), "\n", sep = "")
  return(invisible(NULL))
}

