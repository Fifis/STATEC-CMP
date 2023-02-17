# StatecCMP 0.0.5 (2022-02-17)

* Reworked contribution functions (`ECM2ARDL()`, `computePropag()`): now they support arbitrary data sets, and compute the residuals automatically
* New function `plotContributionsECM()` to plot ECM contributions conveniently with a breakdown into short-run and long-run ones.

# StatecCMP 0.0.4 (2022-02-07)

* Bundled 4 national calendars as package data, removed the `getCalendars()` function to read them from the network drive.
* Added support for string calendars (e.g. `diagnoseSeasonality(..., calendar = "Luxembourg")`) for Belgium, France, Luxembourg, and Germany.
* New function to read all Excel sheets as a list for easier contribution computation, `readAllSheets()`.

# StatecCMP 0.0.3

* New functions: `splitOverlapFixed()` for fixed-length chunks of long series, `splitOverlapCustom()` for completely arbitrary user-defined break points.
* Renamed `getStat()` to `getSAStat()`.
* Added non-robust tests to `robustSeasTests()`.
* Improved convergence of `diagnoseSeasonalityOne()` via gradual relaxation of stopping criteria.

# StatecCMP 0.0.2

* New separate functions `diagnoseSeasonalityOne()` and `getStatOne()` to operate on single series / adjustments, and `diagnoseSeasonality()` and `getStat()` as wrappers for output lists of the former.
* New function `diagnoseRevisions()` to compute the summary of discrepancies between two series.
* New function `robustSeasTests()` to conduct regression-based seasonality tests with arbitrary heteroskedasticity and auto-correlation of errors.
* Support for adjustment of long series with splits and blending of overlapping ends.
* Improved convergence of the `diagnoseSeasonalityOne()` via safety fall-backs.

# StatecCMP 0.0.1 (2022-12-30)

* Initial commit with seasonal adjustment, contribution computation, and error-correction model estimation functions.
