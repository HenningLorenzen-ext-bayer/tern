#' Control Function for Descriptive Statistics
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Sets a list of parameters for summaries of descriptive statistics. Typically used internally to specify
#' details for [s_summary].
#'
#' @inheritParams argument_convention
#' @param quantiles (`numeric`)\cr of length two to specify the quantiles to calculate.
#' @param quantile_type (`numeric`)\cr between 1 and 9 selecting quantile algorithms to be used.
#'   Default is set to 2 as this matches the default quantile algorithm in SAS `proc univariate` set by `QNTLDEF=5`.
#'   This differs from R's default. See more about `type` in [stats::quantile()].
#' @param test_mean (`numeric`)\cr to test against the mean under the null hypothesis when calculating p-value.
#'
#' @return A list of components with the same names as the arguments.
#' @export
#'
control_summarize_vars <- function(conf_level = 0.95,
                                   quantiles = c(0.25, 0.75),
                                   quantile_type = 2,
                                   test_mean = 0) {
  checkmate::assert_vector(quantiles, len = 2)
  checkmate::assert_int(quantile_type, lower = 1, upper = 9)
  checkmate::assert_numeric(test_mean)
  nullo <- lapply(quantiles, assert_proportion_value)
  assert_proportion_value(conf_level)
  list(conf_level = conf_level, quantiles = quantiles, quantile_type = quantile_type, test_mean = test_mean)
}

#' Format Function for Descriptive Statistics
#'
#' @description
#'
#' Returns format patterns for descriptive statistics.
#' The format is understood by the `rtables`.
#'
#' @param type (`string`)\cr choice of a summary data type.
#' Only `counts` and `numeric` types are currently supported.
#'
#' @keywords internal
summary_formats <- function(type = "numeric") {
  if (type == "counts") {
    c(
      n = "xx.",
      count = "xx.",
      count_fraction = format_count_fraction,
      n_blq = "xx."
    )
  } else {
    c(
      n = "xx.",
      sum = "xx.x",
      mean = "xx.x",
      sd = "xx.x",
      se = "xx.x",
      mean_sd = "xx.x (xx.x)",
      mean_se = "xx.x (xx.x)",
      mean_ci = "(xx.xx, xx.xx)",
      mean_sei = "(xx.xx, xx.xx)",
      mean_sdi = "(xx.xx, xx.xx)",
      mean_pval = "xx.xx",
      median = "xx.x",
      mad = "xx.x",
      median_ci = "(xx.xx, xx.xx)",
      quantiles = "xx.x - xx.x",
      iqr = "xx.x",
      range = "xx.x - xx.x",
      cv = "xx.x",
      min = "xx.x",
      max = "xx.x",
      geom_mean = "xx.x",
      geom_cv = "xx.x"
    )
  }
}

#' Label Function for Descriptive Statistics
#'
#' @description
#'
#' Returns labels of descriptive statistics for numeric variables.
#'
#' @keywords internal
summary_labels <- function() {
  c(
    mean = "Mean",
    sum = "Sum",
    sd = "SD",
    se = "SE",
    mean_sd = "Mean (SD)",
    mean_se = "Mean (SE)",
    median = "Median",
    mad = "Median Absolute Deviation",
    iqr = "IQR",
    range = "Min - Max",
    cv = "CV (%)",
    min = "Minimum",
    max = "Maximum",
    geom_mean = "Geometric Mean",
    geom_cv = "CV % Geometric Mean",
    n = "n"
  )
}

#' Summarize Variables
#'
#' @description `r lifecycle::badge("stable")`
#'
#' We use the new S3 generic function [s_summary()] to implement summaries for
#' different `x` objects. This is used as Statistics Function in combination
#' with the new Analyze Function [summarize_vars()].
#'
#' @name summarize_variables
NULL

#' @inheritParams argument_convention
#'
#' @describeIn summarize_variables `s_summary` is a S3 generic function to produce
#'   an object description.
#'
#' @export
s_summary <- function(x,
                      na.rm = TRUE, # nolint
                      denom,
                      .N_row, # nolint
                      .N_col, # nolint
                      na_level,
                      .var,
                      ...) {
  checkmate::assert_flag(na.rm)
  UseMethod("s_summary", x)
}

#' @describeIn summarize_variables Method for numeric class. Note that,
#'   if `x` is an empty vector, `NA` is returned. This is the expected
#'   feature so as to return `rcell` content in `rtables` when the
#'   intersection of a column and a row delimits an empty data selection.
#'   Also, when the `mean` function is applied to an empty vector, `NA` will
#'   be returned instead of `NaN`, the latter being standard behavior in R.
#'
#' @param control (`list`)\cr parameters for descriptive statistics details, specified by using
#'   the helper function [control_summarize_vars()]. Some possible parameter options are:
#'   * `conf_level` (`proportion`)\cr confidence level of the interval for mean and median.
#'   * `quantiles` (`numeric`)\cr vector of length two to specify the quantiles.
#'   * `quantile_type` (`numeric`)\cr between 1 and 9 selecting quantile algorithms to be used.
#'     See more about `type` in [stats::quantile()].
#'   * `test_mean` (`numeric`)\cr value to test against the mean under the null hypothesis when calculating p-value.
#'
#' @return If `x` is of class `numeric`, returns a list with the following named `numeric` items:
#' \describe{
#'   \item{n}{the [length()] of `x`.}
#'   \item{sum}{the [sum()] of `x`.}
#'   \item{mean}{the [mean()] of `x`.}
#'   \item{sd}{the [stats::sd()] of `x`.}
#'   \item{se}{the standard error of `x` mean, i.e.: (`sd(x) / sqrt(length(x))`).}
#'   \item{mean_sd}{the [mean()] and [stats::sd()] of `x`.}
#'   \item{mean_se}{the [mean()] of `x` and its standard error (see above).}
#'   \item{mean_ci}{the CI for the mean of `x` (from [stat_mean_ci()]).}
#'   \item{mean_sei}{the SE interval for the mean of `x`, i.e.: ([mean()] -/+ [stats::sd()] / [sqrt()]).}
#'   \item{mean_sdi}{the SD interval for the mean of `x`, i.e.: ([mean()] -/+ [stats::sd()]).}
#'   \item{mean_pval}{the two-sided p-value of the mean of `x` (from [stat_mean_pval()]).}
#'   \item{median}{the [stats::median()] of `x`.}
#'   \item{mad}{the median absolute deviation of `x`, i.e.: ([stats::median()] of `xc`,
#'     where `xc` = `x` - [stats::median()]).
#'   }
#'   \item{median_ci}{the CI for the median of `x` (from [stat_median_ci()]).}
#'   \item{quantiles}{two sample quantiles of `x` (from [stats::quantile()]).}
#'   \item{iqr}{the [stats::IQR()] of `x`.}
#'   \item{range}{the [range_noinf()] of `x`.}
#'   \item{min}{the [max()] of `x`.}
#'   \item{max}{the [min()] of `x`.}
#'   \item{cv}{the coefficient of variation of `x`, i.e.: ([stats::sd()] / [mean()] * 100).}
#'   \item{geom_mean}{the geometric mean of `x`, i.e.: (`exp(mean(log(x)))`).}
#'   \item{geom_cv}{the geometric coefficient of variation of `x`, i.e.: (`sqrt(exp(sd(log(x)) ^ 2) - 1) * 100`).}
#' }
#' @method s_summary numeric
#'
#' @examples
#' # `s_summary.numeric`
#'
#' ## Basic usage: empty numeric returns NA-filled items.
#' s_summary(numeric())
#'
#' ## Management of NA values.
#' x <- c(NA_real_, 1)
#' s_summary(x, na.rm = TRUE)
#' s_summary(x, na.rm = FALSE)
#'
#' x <- c(NA_real_, 1, 2)
#' s_summary(x, stats = NULL)
#'
#' ## Benefits in `rtables` contructions:
#' require(rtables)
#' dta_test <- data.frame(
#'   Group = rep(LETTERS[1:3], each = 2),
#'   sub_group = rep(letters[1:2], each = 3),
#'   x = 1:6
#' )
#'
#' ## The summary obtained in with `rtables`:
#' basic_table() %>%
#'   split_cols_by(var = "Group") %>%
#'   split_rows_by(var = "sub_group") %>%
#'   analyze(vars = "x", afun = s_summary) %>%
#'   build_table(df = dta_test)
#'
#' ## By comparison with `lapply`:
#' X <- split(dta_test, f = with(dta_test, interaction(Group, sub_group)))
#' lapply(X, function(x) s_summary(x$x))
#'
#' @export
s_summary.numeric <- function(x, # nolint
                              na.rm = TRUE, # nolint
                              denom,
                              .N_row, # nolint
                              .N_col, # nolint
                              na_level,
                              .var,
                              control = control_summarize_vars(),
                              ...) {
  checkmate::assert_numeric(x)

  if (na.rm) {
    x <- x[!is.na(x)]
  }

  y <- list()

  y$n <- c("n" = length(x))

  y$sum <- c("sum" = ifelse(length(x) == 0, NA_real_, sum(x, na.rm = FALSE)))

  y$mean <- c("mean" = ifelse(length(x) == 0, NA_real_, mean(x, na.rm = FALSE)))

  y$sd <- c("sd" = stats::sd(x, na.rm = FALSE))

  y$se <- c("se" = stats::sd(x, na.rm = FALSE) / sqrt(length(stats::na.omit(x))))

  y$mean_sd <- c(y$mean, "sd" = stats::sd(x, na.rm = FALSE))

  y$mean_se <- c(y$mean, y$se)

  mean_ci <- stat_mean_ci(x, conf_level = control$conf_level, na.rm = FALSE, gg_helper = FALSE)
  y$mean_ci <- formatters::with_label(mean_ci, paste("Mean", f_conf_level(control$conf_level)))

  mean_sei <- y$mean[[1]] + c(-1, 1) * stats::sd(x, na.rm = FALSE) / sqrt(y$n)
  names(mean_sei) <- c("mean_sei_lwr", "mean_sei_upr")
  y$mean_sei <- formatters::with_label(mean_sei, "Mean -/+ 1xSE")

  mean_sdi <- y$mean[[1]] + c(-1, 1) * stats::sd(x, na.rm = FALSE)
  names(mean_sdi) <- c("mean_sdi_lwr", "mean_sdi_upr")
  y$mean_sdi <- formatters::with_label(mean_sdi, "Mean -/+ 1xSD")

  mean_pval <- stat_mean_pval(x, test_mean = control$test_mean, na.rm = FALSE, n_min = 2)
  y$mean_pval <- formatters::with_label(mean_pval, paste("Mean", f_pval(control$test_mean)))

  y$median <- c("median" = stats::median(x, na.rm = FALSE))

  y$mad <- c("mad" = stats::median(x - y$median, na.rm = FALSE))

  median_ci <- stat_median_ci(x, conf_level = control$conf_level, na.rm = FALSE, gg_helper = FALSE)
  y$median_ci <- formatters::with_label(median_ci, paste("Median", f_conf_level(control$conf_level)))

  q <- control$quantiles
  if (any(is.na(x))) {
    qnts <- rep(NA_real_, length(q))
  } else {
    qnts <- stats::quantile(x, probs = q, type = control$quantile_type, na.rm = FALSE)
  }
  names(qnts) <- paste("quantile", q, sep = "_")
  y$quantiles <- formatters::with_label(qnts, paste0(paste(paste0(q * 100, "%"), collapse = " and "), "-ile"))

  y$iqr <- c("iqr" = ifelse(
    any(is.na(x)),
    NA_real_,
    stats::IQR(x, na.rm = FALSE, type = control$quantile_type)
  ))

  y$range <- stats::setNames(range_noinf(x, na.rm = FALSE), c("min", "max"))
  y$min <- y$range[1]
  y$max <- y$range[2]

  y$cv <- c("cv" = unname(y$sd) / unname(y$mean) * 100)

  # Convert negative values to NA for log calculation.
  x_no_negative_vals <- x
  x_no_negative_vals[x_no_negative_vals <= 0] <- NA
  y$geom_mean <- c("geom_mean" = exp(mean(log(x_no_negative_vals), na.rm = FALSE)))
  geom_mean_ci <- stat_mean_ci(x, conf_level = control$conf_level, na.rm = FALSE, gg_helper = FALSE, geom_mean = TRUE)
  y$geom_mean_ci <- formatters::with_label(geom_mean_ci, paste("Geometric Mean", f_conf_level(control$conf_level)))

  y$geom_cv <- c("geom_cv" = sqrt(exp(stats::sd(log(x_no_negative_vals), na.rm = FALSE) ^ 2) - 1) * 100) # styler: off

  y
}

#' @describeIn summarize_variables Method for factor class. Note that,
#'   if `x` is an empty factor, then still a list is returned for `counts` with one element
#'   per factor level. If there are no levels in `x`, the function fails. If `x` contains `NA`,
#'   it is expected that `NA` have been conveyed to `na_level` appropriately beforehand with
#'   [df_explicit_na()] or [explicit_na()].
#'
#' @param denom (`string`)\cr choice of denominator for factor proportions. Options are:
#'   - `n`: number of values in this row and column intersection.
#'   - `N_row`: total number of values in this row across columns.
#'   - `N_col`: total number of values in this column across rows.
#'
#' @return If `x` is of class `factor` or converted from `character`, returns a list with named `numeric` items:
#' \describe{
#'   \item{n}{the [length()] of `x`.}
#'   \item{count}{a list with the number of cases for each level of the factor `x`.}
#'   \item{count_fraction}{similar to `count` but also includes the proportion of cases for each level of the
#'     factor `x` relative to the denominator, or `NA` if the denominator is zero.
#'   }
#' }
#' @method s_summary factor
#'
#' @examples
#' # `s_summary.factor`
#'
#' ## Basic usage:
#' s_summary(factor(c("a", "a", "b", "c", "a")))
#' # Empty factor returns NA-filled items.
#' s_summary(factor(levels = c("a", "b", "c")))
#'
#' ## Management of NA values.
#' x <- factor(c(NA, "Female"))
#' x <- explicit_na(x)
#' s_summary(x, na.rm = TRUE)
#' s_summary(x, na.rm = FALSE)
#'
#' ## Different denominators.
#' x <- factor(c("a", "a", "b", "c", "a"))
#' s_summary(x, denom = "N_row", .N_row = 10L)
#' s_summary(x, denom = "N_col", .N_col = 20L)
#'
#' @export
s_summary.factor <- function(x,
                             na.rm = TRUE, # nolint
                             denom = c("n", "N_row", "N_col"),
                             .N_row, # nolint
                             .N_col, # nolint
                             na_level = "<Missing>",
                             ...) {
  assert_valid_factor(x, any.missing = FALSE)
  denom <- match.arg(denom)

  if (na.rm) x <- fct_discard(x, na_level)

  y <- list()

  y$n <- length(x)

  y$count <- as.list(table(x, useNA = "ifany"))
  dn <- switch(denom,
    n = length(x),
    N_row = .N_row,
    N_col = .N_col
  )
  y$count_fraction <- lapply(
    y$count,
    function(x) {
      c(x, ifelse(dn > 0, x / dn, 0))
    }
  )

  y$n_blq <- sum(grepl("BLQ|LTR|<[1-9]", x)) # nolint

  y
}

#' @describeIn summarize_variables Method for character class. This makes an automatic
#'   conversion to factor (with a warning) and then forwards to the method for factors.
#' @param verbose defaults to `TRUE`. It prints out warnings and messages. It is mainly used
#'   to print out information about factor casting.
#' @note Automatic conversion of character to factor does not guarantee that the table
#'   can be generated correctly. In particular for sparse tables this very likely can fail.
#'   It is therefore better to always pre-process the dataset such that factors are manually
#'   created from character variables before passing the dataset to [rtables::build_table()].
#' @method s_summary character
#'
#' @export
#'
#' @examples
#' # `s_summary.character`
#'
#' ## Basic usage:
#' s_summary(c("a", "a", "b", "c", "a"), .var = "x", verbose = FALSE)
#' s_summary(c("a", "a", "b", "c", "a", ""), .var = "x", na.rm = FALSE, verbose = FALSE)
s_summary.character <- function(x,
                                na.rm = TRUE, # nolint
                                denom = c("n", "N_row", "N_col"),
                                .N_row, # nolint
                                .N_col, # nolint
                                na_level = "<Missing>",
                                .var,
                                verbose = TRUE,
                                ...) {
  y <- as_factor_keep_attributes(x, x_name = .var, na_level = na_level, verbose = verbose)
  s_summary(
    x = y,
    na.rm = na.rm,
    na_level = na_level,
    denom = denom,
    .N_row = .N_row,
    .N_col = .N_col,
    ...
  )
}

#' @describeIn summarize_variables Method for logical class.
#'
#' @param denom (`string`)\cr choice of denominator for proportion. Options are:
#'   - `n`: number of values in this row and column intersection.
#'   - `N_row`: total number of values in this row across columns.
#'   - `N_col`: total number of values in this column across rows.
#'
#' @return If `x` is of class `logical`, returns a list with named `numeric` items:
#' \describe{
#'   \item{n}{the [length()] of `x` (possibly after removing `NA`s).}
#'   \item{count}{count of `TRUE` in `x`.}
#'   \item{count_fraction}{count and proportion of `TRUE` in `x` relative to the denominator, or `NA` if the
#'     denominator is zero. Note that `NA`s in `x` are never counted or leading to `NA` here.
#'   }
#' }
#' @method s_summary logical
#'
#' @examples
#' # `s_summary.logical`
#'
#' ## Basic usage:
#' s_summary(c(TRUE, FALSE, TRUE, TRUE))
#'
#' ## Management of NA values.
#' x <- c(NA, TRUE, FALSE)
#' s_summary(x, na.rm = TRUE)
#' s_summary(x, na.rm = FALSE)
#'
#' ## Different denominators.
#' x <- c(TRUE, FALSE, TRUE, TRUE)
#' s_summary(x, denom = "N_row", .N_row = 10L)
#' s_summary(x, denom = "N_col", .N_col = 20L)
#'
#' @export
s_summary.logical <- function(x,
                              na.rm = TRUE, # nolint
                              denom = c("n", "N_row", "N_col"),
                              .N_row, # nolint
                              .N_col, # nolint
                              ...) {
  denom <- match.arg(denom)
  if (na.rm) x <- x[!is.na(x)]
  y <- list()
  y$n <- length(x)
  count <- sum(x, na.rm = TRUE)
  dn <- switch(denom,
    n = length(x),
    N_row = .N_row,
    N_col = .N_col
  )
  y$count <- count
  y$count_fraction <- c(count, ifelse(dn > 0, count / dn, NA))
  y$n_blq <- 0L
  y
}

#' @describeIn summarize_variables S3 generic Formatted Analysis function to produce
#'   an object description. It is used as `afun` in [rtables::analyze()].
#'
#' @export
a_summary <- function(x,
                      ...,
                      .N_row, # nolint
                      .N_col, # nolint
                      .var) {
  UseMethod("a_summary", x)
}

.a_summary_numeric_formats <- summary_formats()
.a_summary_numeric_labels <- summary_labels()

#' @describeIn summarize_variables Formatted Analysis function method for `numeric`.
#'
#' @export
#'
#' @examples
#' # `a_summary.numeric`
#' a_summary(rnorm(10), .N_col = 10, .N_row = 20, .var = "bla")
a_summary.numeric <- make_afun(
  s_summary.numeric,
  .formats = .a_summary_numeric_formats,
  .labels = .a_summary_numeric_labels
)

.a_summary_counts_formats <- summary_formats(type = "counts")

#' @describeIn summarize_variables Method for `factor`.
#'
#' @export
#'
#' @examples
#' # `a_summary.factor`
#' # We need to ungroup `count` and `count_fraction` first so that the rtables formatting
#' # functions can be applied correctly.
#' afun <- make_afun(
#'   getS3method("a_summary", "factor"),
#'   .ungroup_stats = c("count", "count_fraction")
#' )
#' afun(factor(c("a", "a", "b", "c", "a")), .N_row = 10, .N_col = 10)
a_summary.factor <- make_afun(
  s_summary.factor,
  .formats = .a_summary_counts_formats
)

#' @describeIn summarize_variables Formatted Analysis function method for `character`.
#'
#' @export
#'
#' @examples
#' # `a_summary.character`
#' afun <- make_afun(
#'   getS3method("a_summary", "character"),
#'   .ungroup_stats = c("count", "count_fraction")
#' )
#' afun(c("A", "B", "A", "C"), .var = "x", .N_col = 10, .N_row = 10, verbose = FALSE)
a_summary.character <- make_afun(
  s_summary.character,
  .formats = .a_summary_counts_formats
)

#' @describeIn summarize_variables Formatted Analysis function method for `logical`.
#'
#' @export
#'
#' @examples
#' # `a_summary.logical`
#' afun <- make_afun(
#'   getS3method("a_summary", "logical")
#' )
#' afun(c(TRUE, FALSE, FALSE, TRUE, TRUE), .N_row = 10, .N_col = 10)
a_summary.logical <- make_afun(
  s_summary.logical,
  .formats = .a_summary_counts_formats
)

#' Constructor Function for [summarize_vars()] and [summarize_colvars()]
#'
#' @description Constructor function which creates a combined Formatted
#' Analysis function for use in layout creating functions [summarize_vars()] and
#' [summarize_colvars()].
#'
#' @note Since [a_summary()] is generic and we want customization of the formatting arguments
#'   via [rtables::make_afun()], we need to create another temporary generic function, with
#'   corresponding customized methods. Then in order for the methods to be found,
#'   we need to wrap them in a combined `afun`. Since this is required by two layout creating
#'   functions (and possibly others in the future), we provide a constructor that does this:
#'   [create_afun_summary()].
#' @inheritParams argument_convention
#'
#' @export
#'
#' @examples
#' # `create_afun_summary()` to create combined `afun`
#'
#' afun <- create_afun_summary(
#'   .stats = NULL,
#'   .formats = c(median = "xx."),
#'   .labels = c(median = "My median"),
#'   .indent_mods = c(median = 1L)
#' )
#' ## Fabricated dataset.
#' dta_test <- data.frame(
#'   USUBJID = rep(1:6, each = 3),
#'   PARAMCD = rep("lab", 6 * 3),
#'   AVISIT  = rep(paste0("V", 1:3), 6),
#'   ARM     = rep(LETTERS[1:3], rep(6, 3)),
#'   AVAL    = c(9:1, rep(NA, 9))
#' )
#'
#' l <- basic_table() %>%
#'   split_cols_by(var = "ARM") %>%
#'   split_rows_by(var = "AVISIT") %>%
#'   analyze(vars = "AVAL", afun = afun)
#'
#' build_table(l, df = dta_test)
create_afun_summary <- function(.stats, .formats, .labels, .indent_mods) {
  function(x,
           ...,
           .N_row, # nolint
           .N_col, # nolint
           .var) {
    afun <- function(x, ...) {
      UseMethod("afun", x)
    }

    numeric_stats <- afun_selected_stats(
      .stats,
      all_stats = names(.a_summary_numeric_formats)
    )
    afun.numeric <- make_afun( # nolint
      a_summary.numeric,
      .stats = numeric_stats,
      .formats = extract_by_name(.formats, numeric_stats),
      .labels = extract_by_name(.labels, numeric_stats),
      .indent_mods = extract_by_name(.indent_mods, numeric_stats)
    )

    factor_stats <- afun_selected_stats(.stats, c("n", "count", "count_fraction"))
    ungroup_stats <- afun_selected_stats(.stats, c("count", "count_fraction"))
    afun.factor <- make_afun( # nolint
      a_summary.factor,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats),
      .ungroup_stats = ungroup_stats
    )

    afun.character <- make_afun( # nolint
      a_summary.character,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats),
      .ungroup_stats = ungroup_stats
    )

    afun.logical <- make_afun( # nolint
      a_summary.logical,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats)
    )

    afun(
      x = x,
      ...,
      .N_row = .N_row,
      .N_col = .N_col,
      .var = .var
    )
  }
}

#' @describeIn summarize_variables Analyze Function to add a descriptive analyze
#'   layer to `rtables` pipelines. The analysis is applied to a vector and
#'   return the summary, in `rcells`. The ellipsis (`...`) conveys arguments to
#'   [s_summary()], for instance `na.rm = FALSE` if missing data should be
#'   accounted for. When factor variables contains `NA`, it is expected that `NA`
#'   have been conveyed to `na_level` appropriately beforehand with
#'   [df_explicit_na()].
#' @inheritParams rtables::analyze
#' @param ... arguments passed to `s_summary()`.
#'
#' @template formatting_arguments
#' @export
#'
#' @examples
#' ## Fabricated dataset.
#' dta_test <- data.frame(
#'   USUBJID = rep(1:6, each = 3),
#'   PARAMCD = rep("lab", 6 * 3),
#'   AVISIT  = rep(paste0("V", 1:3), 6),
#'   ARM     = rep(LETTERS[1:3], rep(6, 3)),
#'   AVAL    = c(9:1, rep(NA, 9))
#' )
#'
#' # `summarize_vars()` in `rtables` pipelines
#' ## Default output within a `rtables` pipeline.
#' l <- basic_table() %>%
#'   split_cols_by(var = "ARM") %>%
#'   split_rows_by(var = "AVISIT") %>%
#'   summarize_vars(vars = "AVAL")
#'
#' build_table(l, df = dta_test)
#'
#' ## Select and format statistics output.
#' l <- basic_table() %>%
#'   split_cols_by(var = "ARM") %>%
#'   split_rows_by(var = "AVISIT") %>%
#'   summarize_vars(
#'     vars = "AVAL",
#'     .stats = c("n", "mean_sd", "quantiles"),
#'     .formats = c("mean_sd" = "xx.x, xx.x"),
#'     .labels = c(n = "n", mean_sd = "Mean, SD", quantiles = c("Q1 - Q3"))
#'   )
#'
#' results <- build_table(l, df = dta_test)
#' as_html(results)
#'
#' ## Use arguments interpreted by `s_summary`.
#' l <- basic_table() %>%
#'   split_cols_by(var = "ARM") %>%
#'   split_rows_by(var = "AVISIT") %>%
#'   summarize_vars(vars = "AVAL", na.rm = FALSE)
#'
#' results <- build_table(l, df = dta_test)
#'
#' ## Handle `NA` levels first when summarizing factors.
#' dta_test$AVISIT <- NA_character_
#' dta_test <- df_explicit_na(dta_test)
#' l <- basic_table() %>%
#'   split_cols_by(var = "ARM") %>%
#'   summarize_vars(vars = "AVISIT", na.rm = FALSE)
#'
#' results <- build_table(l, df = dta_test)
#' \dontrun{
#' Viewer(results)
#' }
summarize_vars <- function(lyt,
                           vars,
                           var_labels = vars,
                           nested = TRUE,
                           ...,
                           show_labels = "default",
                           table_names = vars,
                           .stats = c("n", "mean_sd", "median", "range", "count_fraction"),
                           .formats = NULL,
                           .labels = NULL,
                           .indent_mods = NULL) {
  afun <- create_afun_summary(.stats, .formats, .labels, .indent_mods)

  analyze(
    lyt = lyt,
    vars = vars,
    var_labels = var_labels,
    afun = afun,
    nested = nested,
    extra_args = list(...),
    inclNAs = TRUE,
    show_labels = show_labels,
    table_names = table_names
  )
}
