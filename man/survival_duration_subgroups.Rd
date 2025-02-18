% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_duration_subgroups.R
\name{survival_duration_subgroups}
\alias{survival_duration_subgroups}
\alias{a_survival_subgroups}
\alias{tabulate_survival_subgroups}
\title{Tabulate Survival Duration by Subgroup}
\usage{
a_survival_subgroups(
  .formats = list(n = "xx", n_events = "xx", n_tot_events = "xx", median = "xx.x", n_tot
    = "xx", hr = list(format_extreme_values(2L)), ci =
    list(format_extreme_values_ci(2L)), pval = "x.xxxx | (<0.0001)")
)

tabulate_survival_subgroups(
  lyt,
  df,
  vars = c("n_tot_events", "n_events", "median", "hr", "ci"),
  time_unit = NULL
)
}
\arguments{
\item{.formats}{(named \code{character} or \code{list})\cr formats for the statistics.}

\item{lyt}{(\code{layout})\cr input layout where analyses will be added to.}

\item{df}{(\code{list})\cr of data frames containing all analysis variables. List should be
created using \code{\link[=extract_survival_subgroups]{extract_survival_subgroups()}}.}

\item{vars}{(\code{character})\cr the name of statistics to be reported among
\code{n_tot_events} (total number of events per group),
\code{n_events} (number of events per group),
\code{n_tot} (total number of observations per group),
\code{n} (number of observations per group),
\code{median} (median survival time),
\code{hr} (hazard ratio),
\code{ci} (confidence interval of hazard ratio) and
\code{pval} (p value of the effect).
Note, one of the statistics \code{n_tot} and \code{n_tot_events}, as well as both \code{hr} and \code{ci}
are required.}

\item{time_unit}{(\code{string})\cr label with unit of median survival time. Default \code{NULL} skips
displaying unit.}

\item{data}{(\code{data.frame})\cr the dataset containing the variables to summarize.}

\item{groups_lists}{(named \code{list} of \code{list})\cr optionally contains for each \code{subgroups} variable a
list, which specifies the new group levels via the names and the
levels that belong to it in the character vectors that are elements of the list.}

\item{label_all}{(\code{string})\cr label for the total population analysis.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Tabulate statistics such as median survival time and hazard ratio for population subgroups.
}
\details{
These functions create a layout starting from a data frame which contains
the required statistics. Tables typically used as part of forest plot.
}
\section{Functions}{
\itemize{
\item \code{a_survival_subgroups()}: Formatted Analysis function used to format the results of
\code{\link[=extract_survival_subgroups]{extract_survival_subgroups()}}. Returns is a list of Formatted Analysis functions with one element per statistic.

\item \code{tabulate_survival_subgroups()}: table creating function.

}}
\examples{
library(dplyr)
library(forcats)

adtte <- tern_ex_adtte

# Save variable labels before data processing steps.
adtte_labels <- formatters::var_labels(adtte)

adtte_f <- adtte \%>\%
  filter(
    PARAMCD == "OS",
    ARM \%in\% c("B: Placebo", "A: Drug X"),
    SEX \%in\% c("M", "F")
  ) \%>\%
  mutate(
    # Reorder levels of ARM to display reference arm before treatment arm.
    ARM = droplevels(fct_relevel(ARM, "B: Placebo")),
    SEX = droplevels(SEX),
    AVALU = as.character(AVALU),
    is_event = CNSR == 0
  )
labels <- c(
  "ARM" = adtte_labels[["ARM"]],
  "SEX" = adtte_labels[["SEX"]],
  "AVALU" = adtte_labels[["AVALU"]],
  "is_event" = "Event Flag"
)
formatters::var_labels(adtte_f)[names(labels)] <- labels

df <- extract_survival_subgroups(
  variables = list(
    tte = "AVAL",
    is_event = "is_event",
    arm = "ARM", subgroups = c("SEX", "BMRKR2")
  ),
  data = adtte_f
)
df
# Internal function - a_survival_subgroups
\dontrun{
a_survival_subgroups(.formats = list("n" = "xx", "median" = "xx.xx"))
}


## Table with default columns.
basic_table() \%>\%
  tabulate_survival_subgroups(df, time_unit = adtte_f$AVALU[1])

## Table with a manually chosen set of columns: adding "pval".
basic_table() \%>\%
  tabulate_survival_subgroups(
    df = df,
    vars = c("n_tot_events", "n_events", "median", "hr", "ci", "pval"),
    time_unit = adtte_f$AVALU[1]
  )
}
\seealso{
\code{\link[=extract_survival_subgroups]{extract_survival_subgroups()}}
}
\keyword{internal}
