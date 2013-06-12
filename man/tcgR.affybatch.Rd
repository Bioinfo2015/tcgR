\name{tcgR.affybatch}
\alias{tcgR.affybatch}
\title{tcgR affybatch}
\description{Create an affybatch object out of TCGA clinical and Level 1
Affymetrix data.}
\usage{tcgR.affybatch(celfile.path = ".", clinical.biotab.file, clinical.id.col = 1, 
    sdrf.file, sdrf.barcode.col = 2, ...)}
\arguments{
  \item{celfile.path}{The path to the CEL files.}
  \item{clinical.biotab.file}{The filename with the TCGA clinical data in Biotab format.}
  \item{clinical.id.col}{Column with the TCGA patient ID.}
  \item{sdrf.file}{The Level 1 CEL sdrf file.}
  \item{sdrf.barcode.col}{Column with the barcode.}
  \item{\dots}{alternative arguments passed to ReadAffy}
}



\author{Markus Riester}





