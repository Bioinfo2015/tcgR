\name{tcgR.maf}
\alias{tcgR.maf}
\title{tcgR maf}
\description{Parse a MAF file.}
\usage{tcgR.maf(filename, file = "tcga_somatic.txt", coding.fun = function(x) ifelse(x == 
    "Silent", 1, 2), sep = "\\t", verbose = TRUE)}
\arguments{
  \item{filename}{The filename of the MAF file. Can be a vector of filenames.}
  \item{file}{The output file}
  \item{coding.fun}{This function is used to code mutations numerically. 0 means no
mutation. Useful for loading the data as a numeric ExpressionSet.}
  \item{sep}{The field separator character. See write.table().}
  \item{verbose}{Print some additional progress information.}
}



\author{Markus Riester}





