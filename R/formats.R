#' Additional Formatting Functions
#'
#' @description `r lifecycle::badge("stable")`
#'
#' This summarizes the additional Formatting Functions to work with `rtables`.
#'
#' @family formatting functions
#' @name formatting_functions
NULL

#' Formatting Fraction and Percentage
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Formats a fraction together with ratio in percent.
#'
#' @param x (`integer`)\cr with elements `num` and `denom`.
#' @param ... required for `rtables` interface.
#' @return a string in the format `num / denom (ratio %)`. If `num` is 0 the format is `num / denom`.
#'
#' @examples
#' format_fraction(x = c(num = 2L, denom = 3L))
#' format_fraction(x = c(num = 0L, denom = 3L))
#'
#' @family formatting functions
#' @export
format_fraction <- function(x, ...) {
  attr(x, "label") <- NULL

  checkmate::assert_vector(x)
  checkmate::assert_count(x["num"])
  checkmate::assert_count(x["denom"])

  result <- if (x["num"] == 0) {
    paste0(x["num"], "/", x["denom"])
  } else {
    paste0(
      x["num"], "/", x["denom"],
      " (", round(x["num"] / x["denom"] * 100, 1), "%)"
    )
  }

  return(result)
}

#' Formatting Fraction and Percentage with Fixed Single Decimal Place
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Formats a fraction together with ratio in percent with fixed single decimal place.
#' Includes trailing zero in case of whole number percentages to always keep one decimal place.
#'
#' @param x (`integer`)\cr with elements `num` and `denom`.
#' @param ... required for `rtables` interface.
#' @return a string in the format `num / denom (ratio %)`. If `num` is 0 the format is `num / denom`.
#'
#' @examples
#' format_fraction_fixed_dp(x = c(num = 1L, denom = 2L))
#' format_fraction_fixed_dp(x = c(num = 1L, denom = 4L))
#' format_fraction_fixed_dp(x = c(num = 0L, denom = 3L))
#'
#' @family formatting functions
#' @export
format_fraction_fixed_dp <- function(x, ...) {
  attr(x, "label") <- NULL
  checkmate::assert_vector(x)
  checkmate::assert_count(x["num"])
  checkmate::assert_count(x["denom"])

  result <- if (x["num"] == 0) {
    paste0(x["num"], "/", x["denom"])
  } else {
    paste0(
      x["num"], "/", x["denom"],
      " (", sprintf("%.1f", round(x["num"] / x["denom"] * 100, 1)), "%)"
    )
  }
  return(result)
}

#' Formatting Count and Fraction
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Formats a count together with fraction with special consideration when count is `0`.
#'
#' @param x (`integer`)\cr vector of length 2, count and fraction.
#' @param ... required for `rtables` interface.
#' @return a string in the format `count (fraction %)`. If `count` is 0 the format is `0`.
#'
#' @examples
#' format_count_fraction(x = c(2, 0.6667))
#' format_count_fraction(x = c(0, 0))
#'
#' @family formatting functions
#' @export
format_count_fraction <- function(x, ...) {
  attr(x, "label") <- NULL

  if (any(is.na(x))) {
    return("NA")
  }

  checkmate::assert_vector(x)
  checkmate::assert_integerish(x[1])
  assert_proportion_value(x[2], include_boundaries = TRUE)

  result <- if (x[1] == 0) {
    "0"
  } else {
    paste0(x[1], " (", round(x[2] * 100, 1), "%)")
  }

  return(result)
}

#' Formatting Count and Percentage with Fixed Single Decimal Place
#'
#' @description `r lifecycle::badge("experimental")`
#'
#' Formats a count together with fraction with special consideration when count is `0`.
#'
#' @param x (`integer`)\cr vector of length 2, count and fraction.
#' @param ... required for `rtables` interface.
#' @return a string in the format `count (fraction %)`. If `count` is 0 the format is `0`.
#'
#' @examples
#' format_count_fraction_fixed_dp(x = c(2, 0.6667))
#' format_count_fraction_fixed_dp(x = c(2, 0.5))
#' format_count_fraction_fixed_dp(x = c(0, 0))
#'
#' @family formatting functions
#' @export
format_count_fraction_fixed_dp <- function(x, ...) {
  attr(x, "label") <- NULL

  if (any(is.na(x))) {
    return("NA")
  }

  checkmate::assert_vector(x)
  checkmate::assert_integerish(x[1])
  assert_proportion_value(x[2], include_boundaries = TRUE)

  result <- if (x[1] == 0) {
    "0"
  } else if (x[2] == 1) {
    sprintf("%d (100%%)", x[1])
  } else {
    sprintf("%d (%.1f%%)", x[1], x[2] * 100)
  }

  return(result)
}

#' Formatting: XX as Formatting Function
#'
#' Translate a string where x and dots are interpreted as number place
#' holders, and others as formatting elements.
#'
#' @param str (`string`)\cr template.
#'
#' @return An `rtables` formatting function.
#'
#' @examples
#' test <- list(c(1.658, 0.5761), c(1e1, 785.6))
#'
#' z <- format_xx("xx (xx.x)")
#' sapply(test, z)
#'
#' z <- format_xx("xx.x - xx.x")
#' sapply(test, z)
#'
#' z <- format_xx("xx.x, incl. xx.x% NE")
#' sapply(test, z)
#'
#' @family formatting functions
#' @export
format_xx <- function(str) {
  # Find position in the string.
  positions <- gregexpr(pattern = "x+\\.x+|x+", text = str, perl = TRUE)
  x_positions <- regmatches(x = str, m = positions)[[1]]

  # Roundings depends on the number of x behind [.].
  roundings <- lapply(
    X = x_positions,
    function(x) {
      y <- strsplit(split = "\\.", x = x)[[1]]
      rounding <- function(x) {
        round(x, digits = ifelse(length(y) > 1, nchar(y[2]), 0))
      }
      return(rounding)
    }
  )

  rtable_format <- function(x, output) {
    values <- Map(y = x, fun = roundings, function(y, fun) fun(y))
    regmatches(x = str, m = positions)[[1]] <- values
    return(str)
  }

  return(rtable_format)
}

#' Formatting Fraction with Lower Threshold
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Formats a fraction when the second element of the input `x` is the fraction. It applies
#' a lower threshold, below which it is just stated that the fraction is smaller than that.
#'
#' @param threshold (`proportion`)\cr lower threshold.
#' @return An `rtables` Formatting Function that takes numeric input `x` where the second
#'   element is the fraction that is formatted. If the fraction is above or equal to the threshold,
#'   then it is displayed in percentage. If it is positive but below the threshold, it returns
#'   "<1" e.g. if the threshold is `0.01`. If it is zero, then just "0" is returned.
#'
#' @examples
#' format_fun <- format_fraction_threshold(0.05)
#' format_fun(x = c(20, 0.1))
#' format_fun(x = c(2, 0.01))
#' format_fun(x = c(0, 0))
#'
#' @family formatting functions
#' @export
format_fraction_threshold <- function(threshold) {
  assert_proportion_value(threshold)
  string_below_threshold <- paste0("<", round(threshold * 100))
  function(x, ...) {
    assert_proportion_value(x[2], include_boundaries = TRUE)
    ifelse(
      x[2] > 0.01,
      round(x[2] * 100),
      ifelse(
        x[2] == 0,
        "0",
        string_below_threshold
      )
    )
  }
}


#' Formatting Extreme Values
#'
#' @description `r lifecycle::badge("stable")`
#'
#' `Rtables` Formatting Functions that handle extreme values.
#'
#' @details For each input, apply a format to the specified number of `digits`. If the value is
#'    below a threshold, it returns "<0.01" e.g. if the number of `digits` is 2. If the value is
#'    above a threshold, it returns ">999.99" e.g. if the number of `digits` is 2.
#'    If it is zero, then returns "0.00".
#'
#' @param digits (`integer`)\cr number of decimal places to display.
#'
#' @family formatting functions
#' @name extreme_format
NULL

#' @describeIn extreme_format Internal helper function to calculate the threshold and create formatted strings
#'  used in Formatting Functions. Returns a list with elements `threshold` and `format_string`.
#'
#' @examples
#' h_get_format_threshold(2L)
#'
#' @export
h_get_format_threshold <- function(digits = 2L) {
  checkmate::assert_integerish(digits)

  low_threshold <- 1 / (10 ^ digits) # styler: off
  high_threshold <- 1000 - (1 / (10 ^ digits)) # styler: off

  string_below_threshold <- paste0("<", low_threshold)
  string_above_threshold <- paste0(">", high_threshold)

  list(
    "threshold" = c(low = low_threshold, high = high_threshold),
    "format_string" = c(low = string_below_threshold, high = string_above_threshold)
  )
}

#' @describeIn extreme_format Internal helper function to apply a threshold format to a value.
#'   Creates a formatted string to be used in Formatting Functions.
#'
#' @param x (`number`)\cr value to format.
#'
#' @examples
#' h_format_threshold(0.001)
#' h_format_threshold(1000)
#'
#' @export
h_format_threshold <- function(x, digits = 2L) {
  if (is.na(x)) {
    return(x)
  }

  checkmate::assert_numeric(x, lower = 0)

  l_fmt <- h_get_format_threshold(digits)

  result <- if (x < l_fmt$threshold["low"] && 0 < x) {
    l_fmt$format_string["low"]
  } else if (x > l_fmt$threshold["high"]) {
    l_fmt$format_string["high"]
  } else {
    sprintf(fmt = paste0("%.", digits, "f"), x)
  }

  unname(result)
}

#' Formatting a Single Extreme Value
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Create Formatting Function for a single extreme value.
#'
#' @inheritParams extreme_format
#' @examples
#'
#' format_fun <- format_extreme_values(2L)
#' format_fun(x = 0.127)
#' format_fun(x = Inf)
#' format_fun(x = 0)
#' format_fun(x = 0.009)
#'
#' @family formatting functions
#' @export
format_extreme_values <- function(digits = 2L) {
  function(x, ...) {
    checkmate::assert_scalar(x, na.ok = TRUE)

    h_format_threshold(x = x, digits = digits)
  }
}

#' Formatting Extreme Values Part of a Confidence Interval
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Formatting Function for extreme values part of a confidence interval. Values
#' are formatted as e.g. "(xx.xx, xx.xx)" if if the number of `digits` is 2.
#'
#' @inheritParams extreme_format
#'
#' @examples
#' format_fun <- format_extreme_values_ci(2L)
#' format_fun(x = c(0.127, Inf))
#' format_fun(x = c(0, 0.009))
#'
#' @family formatting functions
#' @export
format_extreme_values_ci <- function(digits = 2L) {
  function(x, ...) {
    checkmate::assert_vector(x, len = 2)
    l_result <- h_format_threshold(x = x[1], digits = digits)
    h_result <- h_format_threshold(x = x[2], digits = digits)

    paste0("(", l_result, ", ", h_result, ")")
  }
}
