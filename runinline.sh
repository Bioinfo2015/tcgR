#!/bin/sh

cd ..
R --vanilla <<RSCRIPT
library(inlinedocs);
package.skeleton.dx("tcgR")
RSCRIPT

R CMD install --build tcgR
cd tcgR
