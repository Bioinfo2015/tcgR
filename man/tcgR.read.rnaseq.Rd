\name{tcgR.read.rnaseq}
\alias{tcgR.read.rnaseq}
\title{tcgR read rnaseq}
\description{Create a copynumbR.eset input file from Level 3 methylation data}
\usage{tcgR.read.rnaseq(path = ".", file = "tcga_rnaseq.txt", sep = "\\t", 
    filter = "rsem_gene_normalized", ...)}
\arguments{
  \item{path}{The path containing the Level_3 folder.}
  \item{file}{The output file}
  \item{sep}{The field separator character. See write.table().}
  \item{filter}{Only consider the normalized gene level data by default.}
  \item{\dots}{Additional arguments passed to write.table().}
}



\author{Markus Riester}





