#' @include summarize_variables.R
NULL

#' Compare Variables Between Groups
#'
#' @description `r lifecycle::badge("stable")`
#'
#' We use the new S3 generic function [s_compare()] to implement comparisons for
#' different `x` objects. This is used as Statistics Function in combination
#' with the new Analyze Function [compare_vars()].
#'
#' @seealso Relevant constructor function [create_afun_compare()].
#'
#' @name compare_variables
#'
NULL

#' @inheritParams argument_convention
#' @describeIn compare_variables `s_compare` is a S3 generic function to produce
#'   an object description and comparison versus the reference column in the form
#'   of p-values.
#' @seealso [s_summary()] which is used internally for the summary part per column.
#'
#' @export
#'
s_compare <- function(x,
                      .ref_group,
                      .in_ref_col,
                      ...) {
  UseMethod("s_compare", x)
}

#' @describeIn compare_variables Method for numeric class. This uses the standard t-test
#'   to calculate the p-value.
#' @return If `x` is of class `numeric`, returns a list with named items:
#'   - all items from [s_summary.numeric()].
#'   - `pval`: the p-value.
#' @method s_compare numeric
#'
#' @export
#'
#' @examples
#' # `s_compare.numeric`
#'
#' ## Usual case where both this and the reference group vector have more than 1 value.
#' s_compare(rnorm(10, 5, 1), .ref_group = rnorm(5, -5, 1), .in_ref_col = FALSE)
#'
#' ## If one group has not more than 1 value, then p-value is not calculated.
#' s_compare(rnorm(10, 5, 1), .ref_group = 1, .in_ref_col = FALSE)
#'
#' ## Empty numeric does not fail, it returns NA-filled items and no p-value.
#' s_compare(numeric(), .ref_group = numeric(), .in_ref_col = FALSE)
s_compare.numeric <- function(x,
                              .ref_group,
                              .in_ref_col,
                              ...) {
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(.ref_group)
  checkmate::assert_flag(.in_ref_col)

  y <- s_summary.numeric(x = x, ...)

  y$pval <- if (!.in_ref_col && n_available(x) > 1 && n_available(.ref_group) > 1) {
    stats::t.test(x, .ref_group)$p.value
  } else {
    character()
  }

  y
}

#' @describeIn compare_variables Method for factor class. This uses the chi-squared test
#'   to calculate the p-value. Note that the `denom` for factor proportions can only be `n` here
#'   since the usage is for comparing proportions between columns.
#'   Therefore a row-based proportion would not make sense. Also proportion based on `N_col` would
#'   be difficult since for the chi-squared test statistic we use the counts. Therefore
#'   missing values should be accounted for explicitly as factor levels.
#' @param denom (`string`)\cr choice of denominator for factor proportions,
#'   can only be `n` (number of values in this row and column intersection).
#' @return If `x` is of class `factor` or converted from `character`, returns a list with named items:
#'   - all items from [s_summary.factor()].
#'   - `pval`: the p-value.
#' @method s_compare factor
#'
#' @export
#'
#' @examples
#' # `s_compare.factor`
#'
#' ## Basic usage:
#' x <- factor(c("a", "a", "b", "c", "a"))
#' y <- factor(c("a", "b", "c"))
#' s_compare(x = x, .ref_group = y, .in_ref_col = FALSE)
#'
#' ## Management of NA values.
#' x <- explicit_na(factor(c("a", "a", "b", "c", "a", NA, NA)))
#' y <- explicit_na(factor(c("a", "b", "c", NA)))
#' s_compare(x = x, .ref_group = y, .in_ref_col = FALSE, na.rm = TRUE)
#' s_compare(x = x, .ref_group = y, .in_ref_col = FALSE, na.rm = FALSE)
s_compare.factor <- function(x,
                             .ref_group,
                             .in_ref_col,
                             denom = "n",
                             na.rm = TRUE, # nolint
                             na_level = "<Missing>",
                             ...) {
  checkmate::assert_flag(.in_ref_col)
  assert_valid_factor(x, any.missing = FALSE)
  assert_valid_factor(.ref_group, any.missing = FALSE)
  denom <- match.arg(denom)

  y <- s_summary.factor(
    x = x,
    denom = denom,
    na.rm = na.rm,
    na_level = na_level,
    ...
  )

  if (na.rm) {
    x <- fct_discard(x, na_level)
    .ref_group <- fct_discard(.ref_group, na_level)
  }

  checkmate::assert_factor(x, levels = levels(.ref_group), min.levels = 2)

  y$pval <- if (!.in_ref_col && length(x) > 0 && length(.ref_group) > 0) {
    tab <- rbind(table(x), table(.ref_group))
    res <- suppressWarnings(stats::chisq.test(tab))
    res$p.value
  } else {
    character()
  }

  y
}

#' @describeIn compare_variables Method for character class. This makes an automatic
#'   conversion to factor (with a warning) and then forwards to the method for factors.
#' @param verbose defaults to `TRUE`. It prints out warnings and messages. It is mainly used
#'   to print out information about factor casting.
#' @note Automatic conversion of character to factor does not guarantee that the table
#'   can be generated correctly. In particular for sparse tables this very likely can fail.
#'   It is therefore better to always pre-process the dataset such that factors are manually
#'   created from character variables before passing the dataset to [rtables::build_table()].
#' @method s_compare character
#'
#' @examples
#' # `s_compare.character`
#'
#' ## Basic usage:
#' x <- c("a", "a", "b", "c", "a")
#' y <- c("a", "b", "c")
#' s_compare(x, .ref_group = y, .in_ref_col = FALSE, .var = "x", verbose = FALSE)
#'
#' ## Note that missing values handling can make a large difference:
#' x <- c("a", "a", "b", "c", "a", NA)
#' y <- c("a", "b", "c", rep(NA, 20))
#' s_compare(x,
#'   .ref_group = y, .in_ref_col = FALSE,
#'   .var = "x", verbose = FALSE
#' )
#' s_compare(x,
#'   .ref_group = y, .in_ref_col = FALSE, .var = "x",
#'   na.rm = FALSE, verbose = FALSE
#' )
#'
#' @export
s_compare.character <- function(x,
                                .ref_group,
                                .in_ref_col,
                                denom = "n",
                                na.rm = TRUE, # nolint
                                na_level = "<Missing>",
                                .var,
                                verbose = TRUE,
                                ...) {
  x <- as_factor_keep_attributes(x, x_name = .var, na_level = na_level, verbose = verbose)
  .ref_group <- as_factor_keep_attributes(.ref_group, x_name = .var, na_level = na_level, verbose = verbose)
  s_compare(
    x = x,
    .ref_group = .ref_group,
    .in_ref_col = .in_ref_col,
    denom = denom,
    na.rm = na.rm,
    na_level = na_level,
    ...
  )
}

#' @describeIn compare_variables Method for logical class. A chi-squared test
#'   is used. If missing values are not removed, then they are counted as `FALSE`.
#' @return If `x` is of class `logical`, returns a list with named items:
#'   - all items from [s_summary.logical()].
#'   - `pval`: the p-value.
#' @method s_compare logical
#'
#' @export
#'
#' @examples
#' # `s_compare.logical`
#'
#' ## Basic usage:
#' x <- c(TRUE, FALSE, TRUE, TRUE)
#' y <- c(FALSE, FALSE, TRUE)
#' s_compare(x, .ref_group = y, .in_ref_col = FALSE)
#'
#' ## Management of NA values.
#' x <- c(NA, TRUE, FALSE)
#' y <- c(NA, NA, NA, NA, FALSE)
#' s_compare(x, .ref_group = y, .in_ref_col = FALSE, na.rm = TRUE)
#' s_compare(x, .ref_group = y, .in_ref_col = FALSE, na.rm = FALSE)
s_compare.logical <- function(x,
                              .ref_group,
                              .in_ref_col,
                              na.rm = TRUE, # nolint
                              denom = "n",
                              ...) {
  denom <- match.arg(denom)

  y <- s_summary.logical(
    x = x,
    na.rm = na.rm,
    denom = denom,
    ...
  )

  if (na.rm) {
    x <- stats::na.omit(x)
    .ref_group <- stats::na.omit(.ref_group)
  } else {
    x[is.na(x)] <- FALSE
    .ref_group[is.na(.ref_group)] <- FALSE
  }

  y$pval <- if (!.in_ref_col && length(x) > 0 && length(.ref_group) > 0) {
    x <- factor(x, levels = c(TRUE, FALSE))
    .ref_group <- factor(.ref_group, levels = c(TRUE, FALSE))
    tbl <- rbind(table(x), table(.ref_group))
    suppressWarnings(prop_chisq(tbl))
  } else {
    character()
  }

  y
}

#' @describeIn compare_variables S3 generic Formatted Analysis function to produce
#'   an object description and comparison versus the reference column in the form
#'   of p-values. It is used as `afun` in [rtables::analyze()].
#'
#' @export
#'
a_compare <- function(x,
                      .ref_group,
                      .in_ref_col,
                      ...,
                      .var) {
  UseMethod("a_compare", x)
}

#' @describeIn compare_variables Formatted Analysis function method for `numeric`.
#' @export
#'
#' @examples
#' # `a_compare.numeric`
#' a_compare(
#'   rnorm(10, 5, 1),
#'   .ref_group = rnorm(20, -5, 1),
#'   .in_ref_col = FALSE,
#'   .var = "bla"
#' )
a_compare.numeric <- make_afun(
  s_compare.numeric,
  .formats = c(
    .a_summary_numeric_formats,
    pval = "x.xxxx | (<0.0001)"
  ),
  .labels = c(
    .a_summary_numeric_labels,
    pval = "p-value (t-test)"
  ),
  .null_ref_cells = FALSE
)

.a_compare_counts_formats <- c(
  .a_summary_counts_formats,
  pval = "x.xxxx | (<0.0001)"
)

.a_compare_counts_labels <- c(
  pval = "p-value (chi-squared test)"
)

#' @describeIn compare_variables Formatted Analysis function method for `factor`.
#' @export
#'
#' @examples
#' # `a_compare.factor`
#' # We need to ungroup `count` and `count_fraction` first so that the `rtables` formatting
#' # functions can be applied correctly.
#' afun <- make_afun(
#'   getS3method("a_compare", "factor"),
#'   .ungroup_stats = c("count", "count_fraction")
#' )
#' x <- factor(c("a", "a", "b", "c", "a"))
#' y <- factor(c("a", "a", "b", "c"))
#' afun(x, .ref_group = y, .in_ref_col = FALSE)
a_compare.factor <- make_afun(
  s_compare.factor,
  .formats = .a_compare_counts_formats,
  .labels = .a_compare_counts_labels,
  .null_ref_cells = FALSE
)

#' @describeIn compare_variables Formatted Analysis function method for `character`.
#' @export
#'
#' @examples
#' # `a_compare.character`
#' afun <- make_afun(
#'   getS3method("a_compare", "character"),
#'   .ungroup_stats = c("count", "count_fraction")
#' )
#' x <- c("A", "B", "A", "C")
#' y <- c("B", "A", "C")
#' afun(x, .ref_group = y, .in_ref_col = FALSE, .var = "x", verbose = FALSE)
a_compare.character <- make_afun(
  s_compare.character,
  .formats = .a_compare_counts_formats,
  .labels = .a_compare_counts_labels,
  .null_ref_cells = FALSE
)

#' @describeIn compare_variables Formatted Analysis function method for `logical`.
#' @export
#'
#' @examples
#' # `a_compare.logical`
#' afun <- make_afun(
#'   getS3method("a_compare", "logical")
#' )
#' x <- c(TRUE, FALSE, FALSE, TRUE, TRUE)
#' y <- c(TRUE, FALSE)
#' afun(x, .ref_group = y, .in_ref_col = FALSE)
a_compare.logical <- make_afun(
  s_compare.logical,
  .formats = .a_compare_counts_formats,
  .labels = .a_compare_counts_labels,
  .null_ref_cells = FALSE
)

#' Constructor Function for [compare_vars()]
#'
#' @description `r lifecycle::badge("stable")` Constructor function which creates a combined Formatted
#' Analysis function for use in layout creating functions [compare_vars()].
#'
#' @note Since [a_compare()] is generic and we want customization of the formatting arguments
#'   via [rtables::make_afun()], we need to create another temporary generic function, with
#'   corresponding customized methods. Then in order for the methods to be found,
#'   we need to wrap them in a combined `afun`. Since this is required by two layout creating
#'   functions (and possibly others in the future), we provide a constructor that does this:
#'   [create_afun_compare()].
#' @export
#' @inheritParams argument_convention
#' @seealso [compare_variables]
#'
#' @examples
#' # `create_afun_compare()` to create combined `afun`
#'
#' afun <- create_afun_compare(
#'   .stats = c("n", "count_fraction", "mean_sd", "pval"),
#'   .indent_mods = c(pval = 1L)
#' )
#'
#' lyt <- basic_table() %>%
#'   split_cols_by("ARMCD", ref_group = "ARM A") %>%
#'   analyze(
#'     "AGE",
#'     afun = afun,
#'     show_labels = "visible"
#'   )
#' build_table(lyt, df = tern_ex_adsl)
#'
#' lyt <- basic_table() %>%
#'   split_cols_by("ARMCD", ref_group = "ARM A") %>%
#'   analyze(
#'     "SEX",
#'     afun = afun,
#'     show_labels = "visible"
#'   )
#' build_table(lyt, df = tern_ex_adsl)
create_afun_compare <- function(.stats = NULL,
                                .formats = NULL,
                                .labels = NULL,
                                .indent_mods = NULL) {
  function(x,
           .ref_group,
           .in_ref_col,
           ...,
           .var) {
    afun <- function(x, ...) {
      UseMethod("afun", x)
    }

    numeric_stats <- afun_selected_stats(
      .stats,
      all_stats = c(names(.a_summary_numeric_formats), "pval")
    )
    afun.numeric <- make_afun( # nolint
      a_compare.numeric,
      .stats = numeric_stats,
      .formats = extract_by_name(.formats, numeric_stats),
      .labels = extract_by_name(.labels, numeric_stats),
      .indent_mods = extract_by_name(.indent_mods, numeric_stats),
      .null_ref_cells = FALSE
    )

    factor_stats <- afun_selected_stats(
      .stats,
      all_stats = names(.a_compare_counts_formats)
    )
    ungroup_stats <- afun_selected_stats(.stats, c("count", "count_fraction"))
    afun.factor <- make_afun( # nolint
      a_compare.factor,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats),
      .ungroup_stats = ungroup_stats,
      .null_ref_cells = FALSE
    )

    afun.character <- make_afun( # nolint
      a_compare.character,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats),
      .ungroup_stats = ungroup_stats,
      .null_ref_cells = FALSE
    )

    afun.logical <- make_afun( # nolint
      a_compare.logical,
      .stats = factor_stats,
      .formats = extract_by_name(.formats, factor_stats),
      .labels = extract_by_name(.labels, factor_stats),
      .indent_mods = extract_by_name(.indent_mods, factor_stats),
      .null_ref_cells = FALSE
    )

    afun(
      x = x,
      .ref_group = .ref_group,
      .in_ref_col = .in_ref_col,
      ...,
      .var = .var
    )
  }
}

#' @describeIn compare_variables Analyze Function to add a comparison of variables to
#'  `rtables` pipelines. The column split needs to have the reference group defined
#'   via `ref_group` so that the comparison is well defined.
#'   The ellipsis (`...`) conveys arguments to [s_compare()].
#'   When factor variables contains `NA`, it is expected that `NA`
#'   have been conveyed to `na_level` appropriately beforehand with
#'   [df_explicit_na()].
#' @inheritParams rtables::analyze
#' @param ... arguments passed to `s_compare()`.
#'
#' @template formatting_arguments
#' @export
#'
#' @examples
#' # `compare_vars()` in `rtables` pipelines
#'
#' ## Default output within a `rtables` pipeline.
#' lyt <- basic_table() %>%
#'   split_cols_by("ARMCD", ref_group = "ARM B") %>%
#'   compare_vars(c("AGE", "SEX"))
#' build_table(lyt, tern_ex_adsl)
#'
#' ## Select and format statistics output.
#' lyt <- basic_table() %>%
#'   split_cols_by("ARMCD", ref_group = "ARM C") %>%
#'   compare_vars(
#'     vars = "AGE",
#'     .stats = c("mean_sd", "pval"),
#'     .formats = c(mean_sd = "xx.x, xx.x"),
#'     .labels = c(mean_sd = "Mean, SD")
#'   )
#' build_table(lyt, df = tern_ex_adsl)
compare_vars <- function(lyt,
                         vars,
                         var_labels = vars,
                         nested = TRUE,
                         ...,
                         show_labels = "default",
                         table_names = vars,
                         .stats = c("n", "mean_sd", "count_fraction", "pval"),
                         .formats = NULL,
                         .labels = NULL,
                         .indent_mods = NULL) {
  afun <- create_afun_compare(.stats, .formats, .labels, .indent_mods)

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
