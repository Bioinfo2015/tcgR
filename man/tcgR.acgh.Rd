\name{tcgR.acgh}
\alias{tcgR.acgh}
\title{tcgR acgh}
\description{Create a copynumbR.eset input file from Level 2 aCGH data}
\usage{tcgR.acgh(path = ".", file = "tcga_acgh.txt", filter = NULL, 
    adf.file = NULL, adf.col = "Reporter_ID", tumor.tcga.only = TRUE, 
    sdrf.file, sep = "\\t", ...)}
\arguments{
  \item{path}{The path containing the Level_2 folder.}
  \item{file}{The output file}
  \item{filter}{Consider only files passing the specified grep filter. Default all files
in the Level_2 folder.}
  \item{adf.file}{Optionally, the ADF file containing probe information.}
  \item{adf.col}{The column with the probe ids.}
  \item{tumor.tcga.only}{Filter normal and control TCGA samples.}
  \item{sdrf.file}{Translate filename to TCGA barcode with specified sdrf file.}
  \item{sep}{The field separator character. See write.table().}
  \item{\dots}{Additional arguments passed to tcgR.segmented().}
}



\author{Markus Riester}





