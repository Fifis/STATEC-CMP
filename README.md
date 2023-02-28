# STATEC-CMP
Statec CMP repository for internal and/or external use.

## TODO

These are the features / functions / changes that should be addressed before the release of version 0.1.0.

### High-priority

* BUG: ECM: The example for `printECM()` does not work because the merged names are not auto-detected. Fix the internal getECM function because with incomplete SR/LR parts, the output is broken.
* BUG: SA: example of Gabriel (inflation), outlier regressors used in testing: why is 2021.2 skipped?
* FEATURE: SA: When an mts is supplied to `diagnoseSeasonality()`, output to multiple files with series name or increment the number. Consider accepting a vector of file names for writing.
* FEATURE: SA: `diagnoseSeasonality()`: if one chunk is MULT and another is ADD, plot differently (recompute the multiplicative correction); add an optional argument `"default = 'mult'"` or `'add'`.
* FEATURE: SA: `diagnoseSeasonality()`: with mts, try-catch the errors, return, the good results, and report what went wrong
* FEATURE: SA: `getStat()`: add positivity/negativity of adj. check
* DOCUMENTATION: SA: example with outlier and missed level shift
* DOCUMENTATION: SA: example with custom chunks ending with c(NA, NA) lookup to take whatever the last date is
* TEST: SA: "S:/Projets/Conjoncture EPR/SEAS_ADJ/3_Branches/Secteur financier/Nouveaux crédits_cvs.xls"

### Medium-priority

* FEATURE: SA: save the time that it took to adjust the series (for comparabliity and reasonable time estimation; to be accessible by getSAStat).
* FUNCTION: GENERIC: add `readAllSheets0` that would simply read all sheets from an Excel file
* FEATURE: SA: fix sample shrinkage (do not trim away the original series when returning the results in `$series`)  
`plot(a$series[ , c("original", "predicted")], plot.type = "single", lty = 1:2)`
* VISUAL: SA: add a label to plots if the annual totals were forced
* VISUAL: SA: plot the non-forced annual on the plot is SAA is there
* FEATURE: SA: `diagSeas()`: return (as diagnostics) the AICcS for the log transformation and 3 TD models to enable comparison for multiple series
* FEATURE: SA: Return UDG aictest.diff.Leap Year, aictest.diff.e; print this check if `verbose = 1`
* FEATURE: SA: add QS and Shapiro-Wilk to diagnostic statistics
* DOCUMENTATION: SA: `diagnoseSeasonality()`: chunks with internal NAs example
* FEATURE: SA: getSAStat should return names if it gets a list; diagnoseSeas should return a named list for an mts
* FEATURE: SA: produce aggregate full-sample M and Q for series with spans
* VISUAL: SA: show estimation samples as brackets in plots
* CHANGE: ECM: merge getECM and printECM into a unified function that can handle models with different variables in the SR/LR parts.
* FEATURE: ECM: add a wrapper for constrained estimation within a corridor (for linear models and ECMs; the latter require the most care about the long-run part)
* FEATURE: SA: check the number of auto-detected outlier regressors; use a reasonable rule of thumb (i.e. 1 per year) to warn if there are too many.

### Low-priority

* FEATURE: SA: add a calendar for the Netherlands.
* FEATURE: ECM: create a class for ECMs for easier prediction and estimation. So far, there are too many functions not using any meta-information and doing guess-work based on variable names.
* FUNCTION: SA: re-implement the M statistics explicitly, visualise them
* FEATURE: SA: detect trading-day peaks by looking at the neighbourhood at which the trading-day frequencies are the highest (or if they are lower, but not by much, i.e. 1%); remove those peaks and compare with the rest; are they the highest?
