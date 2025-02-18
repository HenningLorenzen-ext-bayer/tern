% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/survival_biomarkers_subgroups.R
\name{survival_biomarkers_subgroups}
\alias{survival_biomarkers_subgroups}
\alias{tabulate_survival_biomarkers}
\title{Tabulate Biomarker Effects on Survival by Subgroup}
\usage{
tabulate_survival_biomarkers(
  df,
  vars = c("n_tot", "n_tot_events", "median", "hr", "ci", "pval"),
  time_unit = NULL
)
}
\arguments{
\item{df}{(\code{data.frame})\cr containing all analysis variables, as returned by
\code{\link[=extract_survival_biomarkers]{extract_survival_biomarkers()}}.}

\item{vars}{(\code{character})\cr the names of statistics to be reported among:
\itemize{
\item \code{n_tot_events}: total number of events per group.
\item \code{n_tot}: total number of observations per group.
\item \code{median}: median survival time.
\item \code{hr}: hazard ratio.
\item \code{ci}: confidence interval of hazard ratio.
\item \code{pval}: p-value of the effect.
Note, one of the statistics \code{n_tot} and \code{n_tot_events}, as well as both \code{hr} and \code{ci} are required.
}}

\item{time_unit}{(\code{string})\cr label with unit of median survival time. Default \code{NULL} skips
displaying unit.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

Tabulate the estimated effects of multiple continuous biomarker variables
across population subgroups.
}
\details{
These functions create a layout starting from a data frame which contains
the required statistics. The tables are then typically used as input for forest plots.
}
\section{Functions}{
\itemize{
\item \code{tabulate_survival_biomarkers()}: table creating function.

}}
\note{
In contrast to \code{\link[=tabulate_survival_subgroups]{tabulate_survival_subgroups()}} this tabulation function does
not start from an input layout \code{lyt}. This is because internally the table is
created by combining multiple subtables.
}
\examples{
library(dplyr)

adtte <- tern_ex_adtte

# Save variable labels before data processing steps.
adtte_labels <- formatters::var_labels(adtte)

adtte_f <- adtte \%>\%
  filter(PARAMCD == "OS") \%>\%
  mutate(
    AVALU = as.character(AVALU),
    is_event = CNSR == 0
  )
labels <- c("AVALU" = adtte_labels[["AVALU"]], "is_event" = "Event Flag")
formatters::var_labels(adtte_f)[names(labels)] <- labels

df <- extract_survival_biomarkers(
  variables = list(
    tte = "AVAL",
    is_event = "is_event",
    biomarkers = c("BMRKR1", "AGE"),
    strata = "STRATA1",
    covariates = "SEX",
    subgroups = "BMRKR2"
  ),
  data = adtte_f
)
df
## Table with default columns.
tabulate_survival_biomarkers(df)

## Table with a manually chosen set of columns: leave out "pval", reorder.
tab <- tabulate_survival_biomarkers(
  df = df,
  vars = c("n_tot_events", "ci", "n_tot", "median", "hr"),
  time_unit = as.character(adtte_f$AVALU[1])
)

## Finally produce the forest plot.
\dontrun{
g_forest(tab, xlim = c(0.8, 1.2))
}
}
\seealso{
\code{\link[=h_tab_surv_one_biomarker]{h_tab_surv_one_biomarker()}} which is used internally, \code{\link[=extract_survival_biomarkers]{extract_survival_biomarkers()}}.
}
