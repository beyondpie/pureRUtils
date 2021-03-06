% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/string.R
\name{splitSelect}
\alias{splitSelect}
\title{Get the field-th column of each string splitted by given symbols.}
\usage{
splitSelect(strs, split = "_", fixed = TRUE, field = 1)
}
\arguments{
\item{strs, }{1d vector of characters}

\item{split}{characters, default is "_"}

\item{fixed}{bool, default is TRUE}

\item{field}{integer, default is 1}
}
\value{
1d vector of characters
}
\description{
Get the field-th column of each string splitted by given symbols.
}
