% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_change.R
\name{summarize_change}
\alias{summarize_change}
\alias{s_change_from_baseline}
\alias{a_change_from_baseline}
\title{Summarize the Change from Baseline or Absolute Baseline Values}
\usage{
s_change_from_baseline(df, .var, variables, na.rm = TRUE, ...)

a_change_from_baseline(df, .var, variables, na.rm = TRUE, ...)

summarize_change(
  lyt,
  vars,
  ...,
  table_names = vars,
  .stats = c("n", "mean_sd", "median", "range"),
  .formats = NULL,
  .labels = NULL,
  .indent_mods = NULL
)
}
\arguments{
\item{df}{(\code{data.frame})\cr data set containing all analysis variables.}

\item{.var}{(\code{string})\cr single variable name that is passed by \code{rtables} when requested
by a statistics function.}

\item{variables}{(named \code{list} of \code{string})\cr list of additional analysis variables.}

\item{na.rm}{(\code{flag})\cr whether \code{NA} values should be removed from \code{x} prior to analysis.}

\item{...}{additional arguments for the lower level functions.}

\item{lyt}{(\code{layout})\cr input layout where analyses will be added to.}

\item{vars}{(\code{character})\cr variable names for the primary analysis variable to be iterated over.}

\item{table_names}{(\code{character})\cr this can be customized in case that the same \code{vars} are analyzed multiple times,
to avoid warnings from \code{rtables}.}

\item{.stats}{(\code{character})\cr statistics to select for the table.}

\item{.formats}{(named \code{character} or \code{list})\cr formats for the statistics.}

\item{.labels}{(named \code{character})\cr labels for the statistics (without indent).}

\item{.indent_mods}{(named \code{integer})\cr indent modifiers for the labels.}
}
\value{
See \code{\link[=s_summary.numeric]{s_summary.numeric()}} for the return values.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

The primary analysis variable \code{.var} indicates the numerical change from baseline results,
and additional required secondary analysis variables are \code{value} and \code{baseline_flag}.
Depending on the baseline flag, either the absolute baseline values (at baseline)
or the change from baseline values (post-baseline) are then summarized.
}
\section{Functions}{
\itemize{
\item \code{s_change_from_baseline()}: Statistics Function that summarizes baseline or post-baseline visits.

\item \code{a_change_from_baseline()}: Formatted Analysis function which can be further customized by calling
\code{\link[rtables:make_afun]{rtables::make_afun()}} on it. It is used as \code{afun} in \code{\link[rtables:analyze]{rtables::analyze()}}.

\item \code{summarize_change()}: Analyze Function for change from baseline analysis.
To be used after a split on visits in the layout, such that each data
subset only contains either baseline or post-baseline data. Allows additional
formatting options.

}}
\note{
The data in \code{df} must be either all be from baseline or post-baseline visits. Otherwise
an error will be thrown.
}
\examples{
df <- data.frame(
  chg = c(1, 2, 3),
  is_bl = c(TRUE, TRUE, TRUE),
  val = c(4, 5, 6)
)

# Internal function - s_change_from_baseline
\dontrun{
s_change_from_baseline(
  df,
  .var = "chg",
  variables = list(value = "val", baseline_flag = "is_bl")
)
}

# Internal function - a_change_from_baseline
\dontrun{
a_change_from_baseline(
  df,
  .var = "chg",
  variables = list(value = "val", baseline_flag = "is_bl")
)
}


# `summarize_change()`

## Fabricated dataset.
library(dplyr)

dta_test <- data.frame(
  USUBJID = rep(1:6, each = 3),
  AVISIT = rep(paste0("V", 1:3), 6),
  ARM = rep(LETTERS[1:3], rep(6, 3)),
  AVAL = c(9:1, rep(NA, 9))
) \%>\%
  mutate(ABLFLL = AVISIT == "V1") \%>\%
  group_by(USUBJID) \%>\%
  mutate(
    BLVAL = AVAL[ABLFLL],
    CHG = AVAL - BLVAL
  ) \%>\%
  ungroup()

results <- basic_table() \%>\%
  split_cols_by("ARM") \%>\%
  split_rows_by("AVISIT") \%>\%
  summarize_change("CHG", variables = list(value = "AVAL", baseline_flag = "ABLFLL")) \%>\%
  build_table(dta_test)
\dontrun{
Viewer(results)
}

}
\keyword{internal}
