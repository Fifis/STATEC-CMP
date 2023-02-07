#' National calendars for regression trading-day pre-adjustment
#'
#' These datasets contain the calendars for some European countries that were
#' obtained from JDemetra+ 2.2.3 by creating national holidays.
#' They are used in seasonal adjustment procedures (diagnoseSeasonality())
#' when the user requests a national calendar
#' (e.g. diagnoseSeasonality(x, calendar = "Luxembourg"))
#'
#' The first part of the name is the country, the second denotes frequency (monthly or quarterly).
#'
#' The source JD+ file is data-raw/All-National-Calendars-JDemetra.cfgx. It
#' requires a working
#'
#' @name calendars
#' @keywords datasets
"Belgium.M"

#' @rdname calendars
"Belgium.Q"

#' @rdname calendars
"France.M"

#' @rdname calendars
"France.Q"

#' @rdname calendars
"Germany.M"

#' @rdname calendars
"Germany.Q"

#' @rdname calendars
"Luxembourg.M"

#' @rdname calendars
"Luxembourg.Q"
