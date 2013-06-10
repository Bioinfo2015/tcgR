\name{tcgR.affybatch}
\alias{tcgR.affybatch}
\title{tcgR affybatch}
\description{Create a tab separated file for the phenoData slot. Used for Level 1
Affymetrix data.}
\usage{tcgR.affybatch(clinical.biotab.file, clinical.id.col = 1, sdrf.file, 
    sdrf.barcode.col = 2, celfile.path, ...)}
\arguments{
  \item{clinical.biotab.file}{The filename with the TCGA clinical data in Biotab format.}
  \item{clinical.id.col}{Column with the TCGA patient ID.}
  \item{sdrf.file}{The Level 1 CEL sdrf file.}
  \item{sdrf.barcode.col}{Column with the barcode.}
  \item{celfile.path}{The path to the CEL files.}
  \item{\dots}{alternative arguments passed to ReadAffy}
}



\author{Markus Riester}





