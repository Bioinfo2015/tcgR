tcgR.affybatch <- function
### Create an affybatch object out of TCGA clinical and Level 1
### Affymetrix data.
(
celfile.path=".",
### The path to the CEL files.
clinical.biotab.file, 
### The filename with the TCGA clinical data in Biotab format.
clinical.id.col=1,
### Column with the TCGA patient ID.
sdrf.file, 
### The Level 1 CEL sdrf file.
sdrf.barcode.col=2,
### Column with the barcode.
...
### alternative arguments passed to ReadAffy
) {
    clinical <- read.delim(clinical.biotab.file,as.is=TRUE)
    sdrf <- read.delim(sdrf.file, as.is=TRUE)
    pd <- cbind(ID=paste(sdrf[,sdrf.barcode.col],".CEL", sep=""),
        clinical[match(substr(sdrf[,sdrf.barcode.col],1,12),
        clinical[,clinical.id.col]),])
    rownames(pd) <- pd[,1]
    pd <- pd[,-1]
    pd <- new("AnnotatedDataFrame", pd)

    celfiles <- list.celfiles(celfile.path,full.names=TRUE)
    isc <- intersect(sapply(celfiles, basename), sampleNames(pd))
    celfiles <- celfiles[isc %in% sapply(celfiles, basename)]
    
    pd <- pd[isc,]

    ReadAffy(filenames=celfiles, phenoData=pd,...)
}

tcgR.segmented <- function
### Create a segmented copy number file from TCGA Level 3 data
(path=".", 
### The path containing the Level_3 folder.
file="tcga.seg",
### The segemented data output file.
filter="\\.hg18", 
### Consider only files passing the specified grep filter. Default the genome version.
sdrf.file=NULL,
### Translate filename to TCGA barcode with specified sdrf file.
sdrf.barcode.col=1,
### Column with the barcode.
sdrf.filename.col="Derived.Array.Data.Matrix.File",
### Column with the data matrix.
verbose=TRUE,
### Print some additional progress information.
...
### Additional arguments passed to read.delim() for reading the files.
) {
    files <- dir(path, full.names=TRUE)

    if (!is.null(filter)) {
        files <- files[grep(filter, files)]
    }

    if (verbose) cat("Reading", files,sep="\n") 
    data <- lapply(files, read.delim, stringsAsFactors=FALSE,...)

    if (!is.null(sdrf.file)) {
        sdrf <- read.delim(sdrf.file, as.is=TRUE)
        barcode <- sdrf[match(sapply(files, basename), sdrf[,sdrf.filename.col]),
            sdrf.barcode.col]
        data <- lapply(1:length(data), function(i) cbind(barcode=barcode[i],
            data[[i]]))
    }
    
    if (!is.null(file)) {
        for (i in 1:length(data)) {
            write.table(data[[i]], file=file, append=i!=1, col.names=i==1,
            row.names=FALSE, quote=FALSE, sep="\t") 
        }
    }
    do.call(rbind, data)
### A data.frame containing the segmented data.   
}

tcgR.methylation <- function
### Create a copynumbR.eset input file from Level 3 methylation data
(path=".", 
### The path containing the Level_3 folder.
file="tcga_methylation.txt",
### The output file
sep="\t",
### The field separator character. See write.table().
gene.level=TRUE,
### average beta values for each gene
...
### Additional arguments passed to write.table().
) {
    x <- tcgR.segmented(path=path,file=NULL, filter=NULL)
    if (gene.level) {
    data <-
    dcast(gene.symbol~barcode,value.var="beta.value",data=x,
    fun.aggregate=mean, na.rm=TRUE)
    } else {
        data <-
        dcast(probe.name~barcode,value.var="beta.value",data=x)
    }

    data <- data[!is.na(data[,1]),]
    data <- data[data[,1] != "",]

    write.table(data, file=file, sep=sep, row.names=FALSE,...)
}

tcgR.rnaseq <- function
### Create a copynumbR.eset input file from Level 3 methylation data
(path=".", 
### The path containing the Level_3 folder.
file="tcga_rnaseq.txt",
### The output file
sep="\t",
### The field separator character. See write.table().
filter="rsem_gene_normalized",
### Only consider the normalized gene level data by default.
...
### Additional arguments passed to write.table().
) {
    x <- tcgR.segmented(path=path,file=NULL, filter=filter)
    data <-
    dcast(gene_id~barcode,value.var="normalized_count",data=x,
    fun.aggregate=mean, na.rm=TRUE)

    write.table(data, file=file, sep=sep, row.names=FALSE,...)
}

tcgR.acgh <- function
### Create a copynumbR.eset input file from Level 2 aCGH data
(path=".", 
### The path containing the Level_2 folder.
file="tcga_acgh.txt",
### The output file
sep="\t",
### The field separator character. See write.table().
adf.file=NULL,
### Optionally, the ADF file containing probe information.
adf.col="Reporter_ID",
### The column with the probe ids.
tumor.tcga.only=TRUE,
### Filter normal and control TCGA samples.
sdrf.file,
### Translate filename to TCGA barcode with specified sdrf file.
...
### Additional arguments passed to tcgR.segmented().
) {
    x <- tcgR.segmented(path=path,file=NULL, filter=NULL,
    skip=1, sdrf.file=sdrf.file, ...)

    if (tumor.tcga.only) {
        x <- x[grepl("TCGA-\\d\\d-\\d\\d\\d\\d-0", x[,1]),]
    }

    data <-
    dcast(Composite.Element.REF~barcode,value.var="normalizedLog2Ratio",data=x,
    fun.aggregate=mean)

    if (!is.null(adf.file)) {
        adf <- read.delim(adf.file, as.is=TRUE)
        idx <- match(data[,1], adf[, adf.col])
        data <- cbind(data[,1,drop=FALSE], adf[idx,], data[,2:ncol(data)])
    }

    write.table(data, file=file, sep=sep, row.names=FALSE)
}

