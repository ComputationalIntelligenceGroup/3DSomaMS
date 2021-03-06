% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/2_extract_soma.R
\name{segment.soma3d}
\alias{segment.soma3d}
\title{Extract single soma from data}
\usage{
segment.soma3d(soma, sdf.path = tempdir(), broken_sdf.path = tempdir(),
  reconstruction.path = tempdir(), name = "soma")
}
\arguments{
\item{soma}{A filepath or the structure returned by read_VRML or read_PLY representing damaged soma mesh}

\item{sdf.path}{folder where intermediate files with shape diameter function data are saved}

\item{broken_sdf.path}{Folder where intermediate files with broken somas after GMM threshold are saved}

\item{reconstruction.path}{Folder where reconstructed file is stored}

\item{name}{Name of the soma file}
}
\value{
Extracted soma
}
\description{
Apply Gaussian mixture models to segment soma from dendrites
}
\seealso{
\code{read_VRML}
}

