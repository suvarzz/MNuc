#' Get lengths of genes from SGS txt file
#' 
#' @details Read sgd file, where name of genes (name), start of genes (txStart)
#'   and end of genes (txEnd) are specified. Write output text file with
#'   name~length columns.
#'
#' @param file SGS text file.
#' @param outdir Output directory.
#' @param filename Optional argument of file name.
#' @param from Shortest length to filter.
#' @param to Longest length to filter.
#' 
#' @return None
#' @export

genelength.from.sgd <- function(file, outdir, filename="genelength", from=NULL, to=NULL) {
    if (is.null(file)) stop("Specify location of sgd text file.")
    if (is.null(outdir)) stop("Specify output directory.")

    dt <- fread(file, header=TRUE, sep="\t", quote="", na.strings="NA")
    # Get sorted lengh of genes
    gl <- dt[, .(name, length=(txEnd-txStart))][order(length)]

    # Filter by length
    if (is.null(from)) from = min(gl[,length])
    if (is.null(to)) to = max(gl[,length])
    gl <- gl[length >= from][length <= to]
    
    # Write data
    dir.create(outdir, recursive=T)
    write.table(gl, file=paste(outdir, filename, ".txt", sep=""), append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
}