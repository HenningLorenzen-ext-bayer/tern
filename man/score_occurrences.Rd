% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/score_occurrences.R
\name{score_occurrences}
\alias{score_occurrences}
\alias{score_occurrences_cols}
\alias{score_occurrences_subtable}
\alias{score_occurrences_cont_cols}
\title{Occurrence Table Sorting}
\usage{
score_occurrences(table_row)

score_occurrences_cols(...)

score_occurrences_subtable(...)

score_occurrences_cont_cols(...)
}
\arguments{
\item{table_row}{(\code{TableRow})\cr an analysis row in a occurrence table.}

\item{...}{arguments for row or column access, see \code{\link{rtables_access}}: either
\code{col_names} (\code{character}) including the names of the columns which should be used,
or alternatively \code{col_indices} (\code{integer}) giving the indices directly instead.}
}
\value{
\code{\link[=score_occurrences]{score_occurrences()}} returns the sum of counts across all columns of a table row.

\code{\link[=score_occurrences_cols]{score_occurrences_cols()}} returns a function that sums counts across all specified columns
of a table row.

\code{\link[=score_occurrences_subtable]{score_occurrences_subtable()}} returns a function that sums counts in each subtable
across all specified columns.
}
\description{
Functions to score occurrence table subtables and rows which can be used in the
sorting of occurrence tables.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\section{Functions}{
\itemize{
\item \code{score_occurrences()}: Scoring function which sums the counts across all
columns. It will fail if anything else but counts are used.

\item \code{score_occurrences_cols()}: Scoring functions can be produced by this constructor to only include
specific columns in the scoring. See \code{\link[=h_row_counts]{h_row_counts()}} for further information.

\item \code{score_occurrences_subtable()}: Scoring functions produced by this constructor can be used on
subtables: They sum up all specified column counts in the subtable. This is useful when
there is no available content row summing up these counts.

\item \code{score_occurrences_cont_cols()}: Produce score function for sorting table by summing the first content row in
specified columns. Note that this is extending \code{\link[rtables:score_funs]{rtables::cont_n_onecol()}} and \code{\link[rtables:score_funs]{rtables::cont_n_allcols()}}.

}}
\examples{
lyt <- basic_table() \%>\%
  split_cols_by("ARM") \%>\%
  add_colcounts() \%>\%
  analyze_num_patients(
    vars = "USUBJID",
    .stats = c("unique"),
    .labels = c("Total number of patients with at least one event")
  ) \%>\%
  split_rows_by("AEBODSYS", child_labels = "visible", nested = FALSE) \%>\%
  summarize_num_patients(
    var = "USUBJID",
    .stats = c("unique", "nonunique"),
    .labels = c(
      "Total number of patients with at least one event",
      "Total number of events"
    )
  ) \%>\%
  count_occurrences(vars = "AEDECOD")

tbl <- build_table(lyt, tern_ex_adae, alt_counts_df = tern_ex_adsl) \%>\%
  prune_table()

tbl_sorted <- tbl \%>\%
  sort_at_path(path = c("AEBODSYS", "*", "AEDECOD"), scorefun = score_occurrences)

tbl_sorted

score_cols_a_and_b <- score_occurrences_cols(col_names = c("A: Drug X", "B: Placebo"))

# Note that this here just sorts the AEDECOD inside the AEBODSYS. The AEBODSYS are not sorted.
# That would require a second pass of `sort_at_path`.
tbl_sorted <- tbl \%>\%
  sort_at_path(path = c("AEBODSYS", "*", "AEDECOD"), scorefun = score_cols_a_and_b)

tbl_sorted

score_subtable_all <- score_occurrences_subtable(col_names = names(tbl))

# Note that this code just sorts the AEBODSYS, not the AEDECOD within AEBODSYS. That
# would require a second pass of `sort_at_path`.
tbl_sorted <- tbl \%>\%
  sort_at_path(path = c("AEBODSYS"), scorefun = score_subtable_all, decreasing = FALSE)

tbl_sorted
}
\seealso{
\code{\link[=h_row_first_values]{h_row_first_values()}}

\code{\link[=h_row_counts]{h_row_counts()}}
}
