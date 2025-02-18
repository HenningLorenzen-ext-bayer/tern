% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_patients_with_event.R
\name{count_patients_with_event}
\alias{count_patients_with_event}
\alias{s_count_patients_with_event}
\alias{a_count_patients_with_event}
\alias{s_count_patients_with_flags}
\alias{a_count_patients_with_flags}
\alias{count_patients_with_flags}
\title{Count the Number of Patients with a Particular Event}
\usage{
s_count_patients_with_event(
  df,
  .var,
  filters,
  .N_col,
  .N_row,
  denom = c("n", "N_row", "N_col")
)

a_count_patients_with_event(
  df,
  .var,
  filters,
  .N_col,
  .N_row,
  denom = c("n", "N_row", "N_col")
)

count_patients_with_event(
  lyt,
  vars,
  ...,
  table_names = vars,
  .stats = "count_fraction",
  .formats = NULL,
  .labels = NULL,
  .indent_mods = NULL
)

s_count_patients_with_flags(
  df,
  .var,
  flag_variables,
  .N_col,
  .N_row,
  denom = c("n", "N_row", "N_col")
)

a_count_patients_with_flags(
  df,
  .var,
  flag_variables,
  .N_col,
  .N_row,
  denom = c("n", "N_row", "N_col")
)

count_patients_with_flags(
  lyt,
  var,
  var_labels = var,
  show_labels = "hidden",
  ...,
  table_names = paste0("tbl_flags_", var),
  .stats = "count_fraction",
  .formats = NULL,
  .indent_mods = NULL
)
}
\arguments{
\item{df}{(\code{data.frame})\cr data set containing all analysis variables.}

\item{.var}{(\code{character})\cr name of the column that contains the unique identifier.}

\item{filters}{(\code{character})\cr a character vector specifying the column names and flag variables
to be used for counting the number of unique identifiers satisfying such conditions.
Multiple column names and flags are accepted in this format
\code{c("column_name1" = "flag1", "column_name2" = "flag2")}.
Note that only equality is being accepted as condition.}

\item{.N_col}{(\code{count})\cr row-wise N (row group count) for the group of observations being analyzed
(i.e. with no column-based subsetting) that is passed by \code{rtables}.}

\item{.N_row}{(\code{count})\cr column-wise N (column count) for the full column that is passed by \code{rtables}.}

\item{denom}{(\code{string})\cr choice of denominator for proportion. Options are:
\itemize{
\item \code{n}: number of values in this row and column intersection.
\item \code{N_row}: total number of values in this row across columns.
\item \code{N_col}: total number of values in this column across rows.
}}

\item{lyt}{(\code{layout})\cr input layout where analyses will be added to.}

\item{vars}{(\code{character})\cr variable names for the primary analysis variable to be iterated over.}

\item{...}{additional arguments for the lower level functions.}

\item{table_names}{(\code{character})\cr this can be customized in case that the same \code{vars} are analyzed multiple times,
to avoid warnings from \code{rtables}.}

\item{.stats}{(\code{character})\cr statistics to select for the table.}

\item{.formats}{(named \code{character} or \code{list})\cr formats for the statistics.}

\item{.labels}{(named \code{character})\cr labels for the statistics (without indent).}

\item{.indent_mods}{(named \code{integer})\cr indent modifiers for the labels.}

\item{flag_variables}{(\code{character})\cr a character vector specifying the names of \code{logical}
variables from analysis dataset used for counting the number of unique identifiers.}

\item{var}{(\code{string})\cr single variable name that is passed by \code{rtables} when requested
by a statistics function.}

\item{var_labels}{character for label.}

\item{show_labels}{label visibility: one of "default", "visible" and "hidden".}
}
\value{
\code{\link[=s_count_patients_with_event]{s_count_patients_with_event()}} returns the count and fraction of patients with the
defined event.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The primary analysis variable \code{.var} denotes the unique patient identifier.
}
\section{Functions}{
\itemize{
\item \code{s_count_patients_with_event()}: Statistics Function that returns the number and the fraction
of unique identifiers with a particular type of event, e.g. the number and the fraction of patients who
had treatment-emergent adverse events. Note that the user can define a new data column containing
the event of interest.

\item \code{a_count_patients_with_event()}: Formatted Analysis function which can be further
customized by calling \code{\link[rtables:make_afun]{rtables::make_afun()}} on it. It is used as \code{afun} in \code{\link[rtables:analyze]{rtables::analyze()}}.

\item \code{count_patients_with_event()}: Analyze Function which adds the count statistics
to the input layout. Note that additional formatting arguments can be used here.

\item \code{s_count_patients_with_flags()}: Statistics function that returns the number and the fraction
of unique identifiers with each particular flag. Returns a list of totals, counts, counts and
fractions with one element per flag.

\item \code{a_count_patients_with_flags()}: Formatted Analysis function which can be further customized by calling
\code{\link[rtables:make_afun]{rtables::make_afun()}} on it. It is used as \code{afun} in \code{\link[rtables:analyze]{rtables::analyze()}}.

\item \code{count_patients_with_flags()}: Analyze Function which is a modified version of \code{\link[=count_patients_with_event]{count_patients_with_event()}}.
Adds the count statistics to the input layout for multiple flag variables at once.

}}
\examples{
library(dplyr)

# `s_count_patients_with_event()`

s_count_patients_with_event(
  tern_ex_adae,
  .var = "SUBJID",
  filters = c("TRTEMFL" = "Y")
)
s_count_patients_with_event(
  tern_ex_adae,
  .var = "SUBJID",
  filters = c("TRTEMFL" = "Y", "AEOUT" = "FATAL")
)
s_count_patients_with_event(
  tern_ex_adae,
  .var = "SUBJID",
  filters = c("TRTEMFL" = "Y", "AEOUT" = "FATAL"),
  denom = "N_col",
  .N_col = 456
)
# `a_count_patients_with_event()`

a_count_patients_with_event(
  tern_ex_adae,
  .var = "SUBJID",
  filters = c("TRTEMFL" = "Y"),
  .N_col = 100,
  .N_row = 100
)
# `count_patients_with_event()`

lyt <- basic_table() \%>\%
  split_cols_by("ARM") \%>\%
  add_colcounts() \%>\%
  count_values(
    "STUDYID",
    values = "AB12345",
    .stats = "count",
    .labels = c(count = "Total AEs")
  ) \%>\%
  count_patients_with_event(
    "SUBJID",
    filters = c("TRTEMFL" = "Y"),
    .labels = c(count_fraction = "Total number of patients with at least one adverse event"),
    table_names = "tbl_all"
  ) \%>\%
  count_patients_with_event(
    "SUBJID",
    filters = c("TRTEMFL" = "Y", "AEOUT" = "FATAL"),
    .labels = c(count_fraction = "Total number of patients with fatal AEs"),
    table_names = "tbl_fatal"
  ) \%>\%
  count_patients_with_event(
    "SUBJID",
    filters = c("TRTEMFL" = "Y", "AEOUT" = "FATAL", "AEREL" = "Y"),
    .labels = c(count_fraction = "Total number of patients with related fatal AEs"),
    .indent_mods = c(count_fraction = 2L),
    table_names = "tbl_rel_fatal"
  )
build_table(lyt, tern_ex_adae, alt_counts_df = tern_ex_adsl)
# `s_count_patients_with_flags()`

# Add labelled flag variables to analysis dataset.
adae <- tern_ex_adae \%>\%
  mutate(
    fl1 = TRUE,
    fl2 = TRTEMFL == "Y",
    fl3 = TRTEMFL == "Y" & AEOUT == "FATAL",
    fl4 = TRTEMFL == "Y" & AEOUT == "FATAL" & AEREL == "Y"
  )
labels <- c(
  "fl1" = "Total AEs",
  "fl2" = "Total number of patients with at least one adverse event",
  "fl3" = "Total number of patients with fatal AEs",
  "fl4" = "Total number of patients with related fatal AEs"
)
formatters::var_labels(adae)[names(labels)] <- labels

s_count_patients_with_flags(
  adae,
  "SUBJID",
  flag_variables = c("fl1", "fl2", "fl3", "fl4"),
  denom = "N_col",
  .N_col = 1000
)
#  We need to ungroup `count_fraction` first so that the `rtables` formatting
# function `format_count_fraction()` can be applied correctly.

# `a_count_patients_with_flags()`

afun <- make_afun(a_count_patients_with_flags,
  .stats = "count_fraction",
  .ungroup_stats = "count_fraction"
)
afun(
  adae,
  .N_col = 10L, # nolint
  .N_row = 10L,
  .var = "USUBJID",
  flag_variables = c("fl1", "fl2", "fl3", "fl4")
)
# `count_patients_with_flags()`

lyt2 <- basic_table() \%>\%
  split_cols_by("ARM") \%>\%
  add_colcounts() \%>\%
  count_patients_with_flags(
    "SUBJID",
    flag_variables = formatters::var_labels(adae[, c("fl1", "fl2", "fl3", "fl4")]),
    denom = "N_col"
  )
build_table(lyt2, adae, alt_counts_df = tern_ex_adsl)
}
