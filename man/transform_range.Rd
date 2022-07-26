% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lp_ideal_rule_benefit_bounds.R
\name{transform_range}
\alias{transform_range}
\title{Helper Function to Transform Outcome}
\usage{
transform_range(y_range, mean_worse, mean_better, scale)
}
\arguments{
\item{y_range}{the range of outcomes}

\item{mean_worse}{mean of the worse treatment arm}

\item{mean_better}{mean of the better treatment arm}

\item{scale}{whether "higher" or "lower" outcomes are better}
}
\value{
a list of the new range, new mean of the better arm, and new mean of the worse arm
}
\description{
Transforms outcome range and means. First, if original scale means lower outcomes are beneficial,
outcome scale is multiplied by -1. Then, once the scale is such that higher outcomes are beneficial,
the range is shifted to start at 0.
The intervals between possible values in the outcome range are not affected.
}
\examples{
#this is an internal function
}
