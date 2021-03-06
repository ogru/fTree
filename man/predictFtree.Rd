% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{predictFtree}
\alias{predictFtree}
\title{Predict function for ftree}
\usage{
predictFtree(ftreeObj = NULL, .Xnew = NULL)
}
\arguments{
\item{ftreeObj}{- An object of class \code{ftree} as produced by the \code{ftree} function}

\item{.Xnew}{- A data frame with a new set of covariates for which to make a prediction(s).}
}
\details{
The code returns a list of length \code{nrow(.Xnew)}. Each element of the list is a matrix whose
         number of rows equals the number of the original argument values (discretizations) of the
         functions contained in \code{ftreeObj}, while the number of columns is equal to the number of
         trees contained in the \code{ftreeObj} (i.e. it's a vector for the \code{'single'} mode ).
         This function does not work with trees of type \code{"wss"}.
}
\author{
Ogy Grujic \email{ognjengr@gmail.com}
}
