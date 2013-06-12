\name{tcgR.rnaseq}
\alias{tcgR.rnaseq}
\title{tcgR rnaseq}
\description{Create a copynumbR.eset input file from Level 3 methylation data}
\usage{tcgR.rnaseq(path = ".", file = "tcga_rnaseq.txt", filter = "rsem_gene_normalized", 
    sep = "\\t", ...)}
\arguments{
  \item{path}{The path containing the Level_3 folder.}
  \item{file}{The output file}
  \item{filter}{Only consider the normalized gene level data by default.}
  \item{sep}{The field separator character. See write.table().}
  \item{\dots}{Additional arguments passed to write.table().}
}



\author{Markus Riester}





