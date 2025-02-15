% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{DataNorm}
\alias{DataNorm}
\title{Data Normalisation using SCTransform
input- raw data matrix - gene counts x cells
ouput- normalised (default including log transfrom) and corrected data matrix: genes x cells}
\usage{
DataNorm(data, log_transform = c(TRUE, FALSE))
}
\description{
Data Normalisation using SCTransform
input- raw data matrix - gene counts x cells
ouput- normalised (default including log transfrom) and corrected data matrix: genes x cells
}
