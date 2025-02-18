% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/footnotes.R
\name{footnotes<-}
\alias{footnotes<-}
\title{Assign Footnotes}
\usage{
footnotes(x) <- value
}
\arguments{
\item{x}{an object}

\item{value}{character vector}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}}

Assign value to attribute footnote of object \code{x}.

This function is deprecated and will be removed from \code{tern} in version 0.7.11.
Please use the equivalent implementation from \code{rtables} or \code{attr(x, "footnote")} directly.
}
\examples{
x <- table(iris$Species)
footnotes(x) <- "Species are equally distributed"
attributes(x)
}
