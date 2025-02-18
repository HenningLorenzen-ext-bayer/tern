% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/h_pkparam_sort.R
\name{h_pkparam_sort}
\alias{h_pkparam_sort}
\title{Sort \verb{PK PARAM} variable}
\usage{
h_pkparam_sort(pk_data, key_var = "PARAMCD")
}
\arguments{
\item{pk_data}{(\code{data.frame})\cr Pharmacokinetics dataframe}

\item{key_var}{(\code{character})\cr key variable used to merge pk_data and metadata created by \code{d_pkparam()}}
}
\value{
dataframe with \verb{PK PARAM} variable ordered
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
library(dplyr)

adpp <- tern_ex_adpp \%>\% mutate(PKPARAM = factor(paste0(PARAM, " (", AVALU, ")")))
pk_ordered_data <- h_pkparam_sort(adpp)
}
