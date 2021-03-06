% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/0_convert_somas_to_PLY.R
\name{convert_somas_to_PLY}
\alias{convert_somas_to_PLY}
\title{Read VRMLs of somas}
\usage{
convert_somas_to_PLY(read_directory, write_directory, parallel = TRUE)
}
\arguments{
\item{read_directory}{path to the folder where VRML files are placed}

\item{write_directory}{path to the folder where PLY files will be saved}
}
\value{
None
}
\description{
Read all the VRMLs of somas in a folder and export them to PLY format in other folder
}
\examples{
convert_somas_to_PLY(read_directory = system.file("test/VRMLs",package="SomaMS"), write_directory = file.path(tempdir(),"PLYs"), parallel = TRUE)
}

