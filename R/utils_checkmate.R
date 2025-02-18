#' Additional Assertions for `checkmate`
#'
#' @description
#' Additional assertion functions which can be used together with the `checkmate` package.
#'
#' @param x (`any`)\cr object to test.
#' @param df (`data.frame`)\cr data set to test.
#' @param variables (named `list` of `character`)\cr list of variables to test.
#' @param include_boundaries (`logical`)\cr whether to include boundaries when testing
#'   for proportions.
#' @param na_level (`character`)\cr the string you have been using to represent NA or
#'   missing data. For `NA` values please consider using directly `base::is.na` or
#'   similar approaches.
#' @inheritParams checkmate::assert_factor
#' @param (`integer`)\cr minimum number of factor levels. Default is `1`.
#' @param ... a collection of objects to test.
#' @name assertions
NULL

check_list_of_variables <- function(x) {
  # drop NULL elements in list
  x <- Filter(Negate(is.null), x)

  res <- checkmate::check_list(x,
    names = "named",
    min.len = 1,
    any.missing = FALSE,
    types = "character"
  )
  # no empty strings allowed
  if (isTRUE(res)) {
    res <- checkmate::check_character(unlist(x), min.chars = 1)
  }
  return(res)
}
#' @describeIn assertions Check whether `x` is a valid list of variable names.
#'   `NULL` elements of the list `x` are dropped out with `Filter(Negate(is.null), x)`.
#'
#' @examples
#' # Check whether `x` is a valid list of variable names.
#'
#' # Internal function - assert_list_of_variables
#' \dontrun{
#' assert_list_of_variables(list(val = "a"))
#' assert_list_of_variables(list(val = c("a", "b")))
#' assert_list_of_variables(list(val = c("a", "b"), val2 = NULL))
#'
#' # The following calls fail
#' assert_list_of_variables(list(1, 2))
#' assert_list_of_variables(list("bla" = 2))
#' }
#'
#' @keywords internal
assert_list_of_variables <- checkmate::makeAssertionFunction(check_list_of_variables)

check_df_with_variables <- function(df, variables, na_level = NULL) {
  checkmate::assert_data_frame(df)
  assert_list_of_variables(variables)

  # flag for equal variables and column names
  err_flag <- all(unlist(variables) %in% colnames(df))
  checkmate::assert_flag(err_flag)

  if (isFALSE(err_flag)) {
    vars <- setdiff(unlist(variables), colnames(df))
    return(paste(
      deparse(substitute(df)),
      "does not contain all specified variables as column names. Missing from dataframe:",
      paste(vars, collapse = ", ")
    ))
  }
  # checking if na_level is present and in which column
  if (!is.null(na_level)) {
    checkmate::assert_string(na_level)
    res <- unlist(lapply(as.list(df)[unlist(variables)], function(x) any(x == na_level)))
    if (any(res)) {
      return(paste0(
        deparse(substitute(df)), " contains explicit na_level (", na_level,
        ") in the following columns: ", paste0(unlist(variables)[res],
          collapse = ", "
        )
      ))
    }
  }
  return(TRUE)
}
#' @describeIn assertions Check whether `df` is a data frame with the analysis `variables`.
#'   Please notice how this produces an error when not all variables are present in the
#'   data.frame while the opposite is not required.
#'
#' @examples
#' # Check whether `df` contains the analysis `variables`.
#'
#' # Internal function - assert_df_with_variables
#' \dontrun{
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3),
#'   variables = list(val = "a")
#' )
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3),
#'   variables = list(val = c("a", "b"))
#' )
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3),
#'   variables = list(val = c("a", "b"))
#' )
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3, e = "<Missing>"),
#'   variables = list(val = c("a", "b")), na_level = "<Missing>"
#' )
#'
#' # The following calls fail
#' assert_df_with_variables(
#'   df = matrix(1:5, ncol = 2, nrow = 3),
#'   variables = list(val = "a")
#' )
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3),
#'   variables = list(val = c("a", "b", "c"))
#' )
#' assert_df_with_variables(
#'   df = data.frame(a = 5, b = 3, e = "<Missing>"),
#'   variables = list(val = c("a", "b", "e")), na_level = "<Missing>"
#' )
#' }
#'
#' @keywords internal
assert_df_with_variables <- checkmate::makeAssertionFunction(check_df_with_variables)

check_valid_factor <- function(x,
                               min.levels = 1, # nolint
                               max.levels = NULL, # nolint
                               null.ok = TRUE, # nolint
                               any.missing = TRUE, # nolint
                               n.levels = NULL, # nolint
                               len = NULL) {
  # checks on levels insertion
  checkmate::assert_int(min.levels, lower = 1)

  # main factor check
  res <- checkmate::check_factor(x,
    min.levels = min.levels,
    null.ok = null.ok,
    max.levels = max.levels,
    any.missing = any.missing,
    n.levels = n.levels
  )

  # no empty strings allowed
  if (isTRUE(res)) {
    res <- checkmate::check_character(levels(x), min.chars = 1)
  }

  return(res)
}
#' @describeIn assertions Check whether `x` is a valid factor (has levels and no empty
#'   string levels). Note that `NULL` and `NA` elements are allowed.
#'
#' @examples
#' # Check whether `x` is a valid factor.
#'
#' # Internal function - assert_valid_factor
#' \dontrun{
#' assert_valid_factor(factor(c("a", NULL)))
#' assert_valid_factor(factor(c("a", "b")))
#' assert_valid_factor(factor(c("a", "b")), len = 2)
#' assert_valid_factor(factor(c("a", NA)), any.missing = TRUE)
#' assert_valid_factor(factor("A", levels = c("A", "B")))
#'
#' # The following calls fail
#' assert_valid_factor(-1)
#' assert_valid_factor(factor(c("a", "")))
#' assert_valid_factor(factor(c("a", NA)), any.missing = FALSE)
#' assert_valid_factor(factor(NULL))
#' assert_valid_factor(factor(c(NULL, "")))
#' assert_valid_factor(factor())
#' }
#'
#' @keywords internal
assert_valid_factor <- checkmate::makeAssertionFunction(check_valid_factor)


check_df_with_factors <- function(df,
                                  variables,
                                  min.levels = 1, # nolint
                                  max.levels = NULL, # nolint
                                  any.missing = TRUE, # nolint
                                  na_level = NULL) {
  res <- check_df_with_variables(df, variables, na_level)
  # checking if all the columns specified by variables are valid factors
  if (isTRUE(res)) {
    # searching the data.frame with selected columns (variables) as a list
    res <- lapply(
      X = as.list(df)[unlist(variables)],
      FUN = check_valid_factor,
      min.levels = min.levels,
      max.levels = max.levels,
      any.missing = any.missing
    )
    res_lo <- unlist(vapply(res, Negate(isTRUE), logical(1)))
    if (any(res_lo)) {
      return(paste0(
        deparse(substitute(df)), " does not contain only factor variables among:",
        "\n* Column `", paste0(unlist(variables)[res_lo],
          "` of the data.frame -> ", res[res_lo],
          collapse = "\n* "
        )
      ))
    } else {
      res <- TRUE
    }
  }
  return(res)
}
#' @describeIn assertions Check whether `df` is a data frame where the analysis `variables`
#'   are all factors. Note that the creation of `NA` by direct call of `factor()` will
#'   trim `NA` levels out of the vector list itself.
#'
#' @examples
#' # Check whether `df` contains all factor analysis `variables`.
#' adf <- data.frame(a = factor(c("A", "B")), b = 3)
#' bdf <- data.frame(a = factor(letters[1:3]), b = factor(c(1, 2, 3)), d = 3)
#'
#' # Internal function - assert_df_with_factors
#' \dontrun{
#' assert_df_with_factors(df = adf, variables = list(val = "a"))
#' assert_df_with_factors(df = adf, variables = list(val = "a"), min.levels = 1)
#' assert_df_with_factors(df = adf, variables = list(val = "a"), min.levels = 2, max.levels = 2)
#' assert_df_with_factors(
#'   df = data.frame(a = factor(c("A", NA, "B")), b = 3),
#'   variable = list(val = "a"),
#'   min.levels = 2,
#'   max.levels = 2
#' )
#'
#' # The following calls fail
#' assert_df_with_factors(df = adf, variables = list(val = "a"), min.levels = 1, max.levels = 1)
#' assert_df_with_factors(df = adf, variables = list(val = "a"), min.levels = 1, max.levels = 1)
#' assert_df_with_factors(df = adf, variables = list(val = "a", val = "b", val = ""))
#' assert_df_with_factors(df = adf, variables = list(val = "a", val = "b", val = "d"))
#' assert_df_with_factors(
#'   df = bdf,
#'   variables = list(val = "a", val = "b"),
#'   min.levels = 1,
#'   max.levels = 1
#' )
#' }
#'
#' @keywords internal
assert_df_with_factors <- checkmate::makeAssertionFunction(check_df_with_factors)

#' @describeIn assertions Check whether `x` is a proportion: number between 0 and 1.
#'
#' @examples
#' # Check whether `x` is between 0 and 1.
#' # Internal function - assert_proportion_value
#' \dontrun{
#' assert_proportion_value(x = 0, include_boundaries = TRUE)
#' assert_proportion_value(x = 0.3)
#'
#' # These fail
#' assert_proportion_value(x = 1.3)
#' assert_proportion_value(x = 1)
#' }
#'
#' @keywords internal
assert_proportion_value <- function(x, include_boundaries = FALSE) {
  checkmate::assert_number(x, lower = 0, upper = 1)
  checkmate::assert_flag(include_boundaries)
  if (isFALSE(include_boundaries)) {
    checkmate::assert_true(x > 0)
    checkmate::assert_true(x < 1)
  }
}
