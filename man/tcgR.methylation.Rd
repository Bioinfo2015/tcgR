\name{tcgR.methylation}
\alias{tcgR.methylation}
\title{tcgR methylation}
\description{Create a copynumbR.eset input file from Level 3 methylation data}
\usage{tcgR.methylation(path = ".", file = "tcga_methylation.txt", filter = NULL, 
    sep = "\\t", gene.level = TRUE, ...)}
\arguments{
  \item{path}{The path containing the Level_3 folder.}
  \item{file}{The output file}
  \item{filter}{Consider only files passing the specified grep filter. Default all files
in the Level_3 folder.}
  \item{sep}{The field separator character. See write.table().}
  \item{gene.level}{Average beta values for each gene.}
  \item{\dots}{Additional arguments passed to write.table().}
}



\author{Markus Riester}





