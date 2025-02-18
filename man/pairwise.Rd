% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coxph.R
\name{pairwise}
\alias{pairwise}
\title{Pairwise Formula Special Term}
\usage{
pairwise(x)
}
\arguments{
\item{x}{the variable for which pairwise result is expected}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

The special term \code{pairwise} indicate that the model should be fitted individually for
every tested level in comparison to the reference level.
}
\details{
Let's \code{ARM} being a factor with level A, B, C; let's be B the reference level,
a model calling the formula including \code{pairwise(ARM)} will result in two models:
\itemize{
\item A model including only levels A and B, and effect of A estimated in reference to B.
\item A model including only levels C and B, the effect of C estimated in reference to B.
}
}
