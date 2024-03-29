% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/select.R
\name{select.l.u}
\alias{select.l.u}
\title{Associated taxa selected by combining lower- and upper-tail p-values}
\usage{
select.l.u(lowerlist, upperlist, com, label)
}
\arguments{
\item{lowerlist}{The index of selected taxa by combining lower-tail p-values.}

\item{upperlist}{The index of selected taxa by combining upper-tail p-values.}

\item{com}{The index of selected one-sided test after combining two one-sided community-level tests.
1: the one-sided test using lower-tail p-values.
2: the one-sided test using upper-tail p-values.}

\item{label}{The lables of all non-rooted nodes.}
}
\value{
A list of selected associated taxa after combining two one-sided community-level tests:
        \item{lower}{Under-presented taxa based on lower-tail p-values.}
        \item{upper}{Over-presented taxa based on upper-tail p-values.}
}
\description{
Associated taxa selected by combining lower- and upper-tail p-values of a test
             using one of abundance representations.
}
