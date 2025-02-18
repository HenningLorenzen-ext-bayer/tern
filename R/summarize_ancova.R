#' Summary for analysis of covariance (ANCOVA).
#'
#' @description `r lifecycle::badge("stable")`
#'
#' Summarize results of ANCOVA. This can be used to analyze multiple endpoints and/or
#' multiple timepoints within the same response variable `.var`.
#'
#' @name summarize_ancova
NULL

#' Helper Function to Return Results of a Linear Model
#'
#' @description `r lifecycle::badge("stable")`
#'
#' @inheritParams argument_convention
#' @param .df_row (`data.frame`)\cr data set that includes all the variables that are called
#'   in `.var` and `variables`.
#' @param variables (named `list` of `strings`)\cr list of additional analysis variables, with
#'   expected elements:
#'   - `arm` (`string`)\cr group variable, for which the covariate adjusted means of multiple
#'     groups will be summarized. Specifically, the first level of `arm` variable is taken as the
#'     reference group.
#'   - `covariates` (`character`)\cr a vector that can contain single variable names (such as
#'     `"X1"`), and/or interaction terms indicated by `"X1 * X2"`.
#' @param interaction_item (`character`)\cr name of the variable that should have interactions
#'   with arm. if the interaction is not needed, the default option is NULL.
#'
#' @examples
#' h_ancova(
#'   .var = "Sepal.Length",
#'   .df_row = iris,
#'   variables = list(arm = "Species", covariates = c("Petal.Length * Petal.Width", "Sepal.Width"))
#' )
#'
#' @export
h_ancova <- function(.var,
                     .df_row,
                     variables,
                     interaction_item = NULL) {
  checkmate::assert_string(.var)
  checkmate::assert_list(variables)
  checkmate::assert_subset(names(variables), c("arm", "covariates"))
  assert_df_with_variables(.df_row, list(rsp = .var))

  arm <- variables$arm
  covariates <- variables$covariates
  if (!is.null(covariates) && length(covariates) > 0) {
    # Get all covariate variable names in the model.
    var_list <- get_covariates(covariates)
    assert_df_with_variables(.df_row, var_list)
  }

  covariates_part <- paste(covariates, collapse = " + ")
  if (covariates_part != "") {
    formula <- stats::as.formula(paste0(.var, " ~ ", covariates_part, " + ", arm))
  } else {
    formula <- stats::as.formula(paste0(.var, " ~ ", arm))
  }

  if (is.null(interaction_item)) {
    specs <- arm
  } else {
    specs <- c(arm, interaction_item)
  }

  lm_fit <- stats::lm(
    formula = formula,
    data = .df_row
  )
  emmeans_fit <- emmeans::emmeans(
    lm_fit,
    # Specify here the group variable over which EMM are desired.
    specs = specs,
    # Pass the data again so that the factor levels of the arm variable can be inferred.
    data = .df_row
  )

  emmeans_fit
}

#' @describeIn summarize_ancova Statistics function that produces a named list of results
#'   of the investigated linear model.
#'
#' @inheritParams argument_convention
#' @inheritParams h_ancova
#' @param interaction_y (`character`)\cr a selected item inside of the interaction_item column
#'   which will be used to select the specific ANCOVA results. if the interaction is not
#'   needed, the default option is FALSE
#' @return A named list of 5 statistics:
#' \describe{
#'   \item{n}{count of complete sample size for the group.}
#'   \item{lsmean}{estimated marginal means in the group.}
#'   \item{lsmean_diff}{difference in estimated marginal means in comparison to the reference group.
#'     If working with the reference group, this will be empty.
#'   }
#'   \item{lsmean_diff_ci}{confidence level for difference in estimated marginal means in comparison to the
#'     reference group.
#'   }
#'   \item{pval}{p-value (not adjusted for multiple comparisons).}
#' }
#'
#' @examples
#' library(dplyr)
#'
#' df <- iris %>% filter(Species == "virginica")
#' .df_row <- iris
#' .var <- "Petal.Length"
#' variables <- list(arm = "Species", covariates = "Sepal.Length * Sepal.Width")
#' .ref_group <- iris %>% filter(Species == "setosa")
#' conf_level <- 0.95
#'
#' # Internal function - s_ancova
#' \dontrun{
#' s_ancova(
#'   df, .var, .df_row, variables, .ref_group,
#'   .in_ref_col = FALSE,
#'   conf_level, interaction_y = FALSE, interaction_item = NULL
#' )
#' }
#'
#' @keywords internal
s_ancova <- function(df,
                     .var,
                     .df_row,
                     variables,
                     .ref_group,
                     .in_ref_col,
                     conf_level,
                     interaction_y = FALSE,
                     interaction_item = NULL) {
  emmeans_fit <- h_ancova(.var = .var, variables = variables, .df_row = .df_row, interaction_item = interaction_item)

  sum_fit <- summary(
    emmeans_fit,
    level = conf_level
  )

  arm <- variables$arm

  sum_level <- as.character(unique(df[[arm]]))

  # Ensure that there is only one element in sum_level.
  checkmate::assert_scalar(sum_level)

  sum_fit_level <- sum_fit[sum_fit[[arm]] == sum_level, ]

  # Get the index of the ref arm
  if (interaction_y != FALSE) {
    y <- unlist(df[(df[[interaction_item]] == interaction_y), .var])
    # convert characters selected in interaction_y into the numeric order
    interaction_y <- which(sum_fit_level[[interaction_item]] == interaction_y)
    sum_fit_level <- sum_fit_level[interaction_y, ]
    # if interaction is called, reset the index
    ref_key <- seq(sum_fit[[arm]][unique(.ref_group[[arm]])])
    ref_key <- tail(ref_key, n = 1)
    ref_key <- (interaction_y - 1) * length(unique(.df_row[[arm]])) + ref_key
  } else {
    y <- df[[.var]]
    # Get the index of the ref arm when interaction is not called
    ref_key <- seq(sum_fit[[arm]][unique(.ref_group[[arm]])])
    ref_key <- tail(ref_key, n = 1)
  }

  if (.in_ref_col) {
    list(
      n = length(y[!is.na(y)]),
      lsmean = formatters::with_label(sum_fit_level$emmean, "Adjusted Mean"),
      lsmean_diff = formatters::with_label(character(), "Difference in Adjusted Means"),
      lsmean_diff_ci = formatters::with_label(character(), f_conf_level(conf_level)),
      pval = formatters::with_label(character(), "p-value")
    )
  } else {
    # Estimate the differences between the marginal means.
    emmeans_contrasts <- emmeans::contrast(
      emmeans_fit,
      # Compare all arms versus the control arm.
      method = "trt.vs.ctrl",
      # Take the arm factor from .ref_group as the control arm.
      ref = ref_key,
      level = conf_level
    )
    sum_contrasts <- summary(
      emmeans_contrasts,
      # Derive confidence intervals, t-tests and p-values.
      infer = TRUE,
      # Do not adjust the p-values for multiplicity.
      adjust = "none"
    )

    sum_contrasts_level <- sum_contrasts[grepl(sum_level, sum_contrasts$contrast), ]
    if (interaction_y != FALSE) {
      sum_contrasts_level <- sum_contrasts_level[interaction_y, ]
    }

    list(
      n = length(y[!is.na(y)]),
      lsmean = formatters::with_label(sum_fit_level$emmean, "Adjusted Mean"),
      lsmean_diff = formatters::with_label(sum_contrasts_level$estimate, "Difference in Adjusted Means"),
      lsmean_diff_ci = formatters::with_label(
        c(sum_contrasts_level$lower.CL, sum_contrasts_level$upper.CL),
        f_conf_level(conf_level)
      ),
      pval = formatters::with_label(sum_contrasts_level$p.value, "p-value")
    )
  }
}

#' @describeIn summarize_ancova Formatted Analysis function which can be further customized by calling
#'   [rtables::make_afun()] on it. It is used as `afun` in [rtables::analyze()].
#'
#' @examples
#' # Internal function - a_ancova
#' \dontrun{
#' a_ancova(
#'   df, .var, .df_row, variables, .ref_group,
#'   .in_ref_col = FALSE,
#'   interaction_y = FALSE, interaction_item = NULL, conf_level
#' )
#' }
#'
#' @keywords internal
a_ancova <- make_afun(
  s_ancova,
  .indent_mods = c("n" = 0L, "lsmean" = 0L, "lsmean_diff" = 0L, "lsmean_diff_ci" = 1L, "pval" = 1L),
  .formats = c(
    "n" = "xx",
    "lsmean" = "xx.xx",
    "lsmean_diff" = "xx.xx",
    "lsmean_diff_ci" = "(xx.xx, xx.xx)",
    "pval" = "x.xxxx | (<0.0001)"
  ),
  .null_ref_cells = FALSE
)

#' @describeIn summarize_ancova Layout creating function which can be be used for creating
#'   summary tables for analysis of covariance (ANCOVA).
#' @inheritParams argument_convention
#' @export
#'
#' @examples
#' basic_table() %>%
#'   split_cols_by("Species", ref_group = "setosa") %>%
#'   add_colcounts() %>%
#'   summarize_ancova(
#'     vars = "Petal.Length",
#'     variables = list(arm = "Species", covariates = NULL),
#'     table_names = "unadj",
#'     conf_level = 0.95, var_labels = "Unadjusted comparison",
#'     .labels = c(lsmean = "Mean", lsmean_diff = "Difference in Means")
#'   ) %>%
#'   summarize_ancova(
#'     vars = "Petal.Length",
#'     variables = list(arm = "Species", covariates = c("Sepal.Length", "Sepal.Width")),
#'     table_names = "adj",
#'     conf_level = 0.95, var_labels = "Adjusted comparison (covariates: Sepal.Length and Sepal.Width)"
#'   ) %>%
#'   build_table(iris)
summarize_ancova <- function(lyt,
                             vars,
                             var_labels,
                             ...,
                             show_labels = "visible",
                             table_names = vars,
                             .stats = NULL,
                             .formats = NULL,
                             .labels = NULL,
                             .indent_mods = NULL,
                             interaction_y = FALSE,
                             interaction_item = NULL) {
  afun <- make_afun(
    a_ancova,
    interaction_y = interaction_y,
    interaction_item = interaction_item,
    .stats = .stats,
    .formats = .formats,
    .labels = .labels,
    .indent_mods = .indent_mods
  )

  analyze(
    lyt,
    vars,
    var_labels = var_labels,
    show_labels = show_labels,
    table_names = table_names,
    afun = afun,
    extra_args = list(...)
  )
}
