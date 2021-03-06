% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SnapATACUtils.R
\name{runKNN}
\alias{runKNN}
\title{Run SnapATAC KNN.}
\usage{
runKNN(
  snapAll = NULL,
  snapAllFile = NULL,
  snapLandmark = NULL,
  snapLandmarkFile = NULL,
  snapQuery = NULL,
  snapQueryFile = NULL,
  removeBmat = TRUE,
  removeJmat = TRUE,
  k = 50,
  dims = 1:30,
  method = "RANN",
  runUMAP = TRUE,
  umapNcores = 1,
  outmmtxFile = NULL,
  outSnapFile = NULL
)
}
\arguments{
\item{snapAll}{SnapObject, default NULL}

\item{snapAllFile}{characters, default NULL}

\item{snapLandmark}{SnapObject, default NULL}

\item{snapLandmarkFile}{characters, default NULL}

\item{snapQuery}{SnapObject, default NULL}

\item{snapQueryFile}{SnapObject, default NULL}

\item{removeBmat}{bool, default TRUE}

\item{removeJmat}{bool, default TRUE}

\item{k}{integer, default 50}

\item{dims}{integer or vector, default is 1:30}

\item{method}{characters, method for KNN (either "RANN" or "Annoy")
default "RANN"}

\item{runUMAP}{bool, default TRUE}

\item{umapNcores}{integer, default 1}

\item{outmmtxFile}{characters, output file of mmtx, default NULL}

\item{outSnapFile}{characters, default NULL}
}
\value{
SnapObject
}
\description{
It will check snapLandmark and snapQuery firstly,
if snapLandmark is not NULL, will use it.
if snapQuery then is also NULL, will merge them.
if both of them are NULL, then use snapAllFile, snapAll in order.
}
