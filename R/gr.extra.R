#' Extend Genomic Range depending on strand
#' 
#' @param upstream Number of nucleotides to extend upstream.
#' @param downstream Number of nucleotides to extend downstream.

#' @return Genomic Range
#' @export

# Source: from Hervé Pagès advice (https://support.bioconductor.org/p/78652)
# Extend Grange depending on strand
gr.extend <- function(x, upstream=0, downstream=0)     
{
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- strand(x) == "+" | strand(x) == "*"
    new_start <- start(x) - ifelse(on_plus, upstream, downstream)
    new_end <- end(x) + ifelse(on_plus, downstream, upstream)
    ranges(x) <- IRanges(new_start, new_end)
    trim(x)
}

################################################################################
#' Extend GRange and trim
#' @description Extend GRange to the certain length and trim till the trim length if necessary.
#' @return Genomic Ranges
#' @export

gr.extend.reads <- function(x, fragmentLen=75, trim=0)
{
    if (!(var(width(x))==0))
        warning("The range width is different. Reads must have identical size") 
    if (fragmentLen < width(x)[1])
        stop("The length of fragments can not be smaller than width of reads") 
    if (trim > fragmentLen)
        stop("The trim length is greater than fragment length")
    down <- fragmentLen - width(x)[1]
    y <- gr.extend(x, downstream=down)
    y <- y - (fragmentLen-trim)/2
}

###############################################################################

#' Save Genomic Ranges List as a set of csv files.
#' 
#' @param grl Named genomic ranges list
#' @param oudir Output directory.
#'
#' @return None
#' 
#' @export

saveGRlist <- function(grl, outdir=NULL) {
  gr_names <- names(grl)
  if (is.null(gr_names)) {
    stop("genomic ranges in a list must be named") }
  if (is.null(outdir)) {
    stop("output directory must be specified") }
  if (!file.exists(outdir)) {
    dir.create(file.path(outdir)) }
  silent <- lapply(names(grl), function(x) {
    write.table(as(grl[[x]], "data.frame"), file = paste(outdir, x, ".csv", sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
    system2("gzip", args=paste(outdir, x, ".csv", sep=""))
  })}

################################################################################

#' Get seqlevels from chromosome sizes txt file
#'
#' @description Transform chromosome sizes in a text file in genomic ranges seqlevels
#'
#' @param chromsizes Tab separated text file containing sizes of all chromosomes in columns (chromosome ~ size).
#'
#' @return seql (named vector)
#'
#' @export

seqlevels.from.chrsizes <- function(chromsizes) {

    if (is.null(chromsizes))
        stop("Please specify \"chromsizes\" chromosome sizes file (chromosome ~ size).\n")

    # Set chromosomes length parameter - seql
    chrom.sizes <- read.table(chromsizes, sep="\t", header=FALSE, quote="", col.names=c("chromosome", "length"))
    seql = as.numeric(chrom.sizes$length)
    names(seql)=as.character(chrom.sizes$chromosome)
	return(seql)
}

################################################################################

#' Get signal
#'
#' @description Get signal from metacolumn of genomic ranges based on given chromosome, start and end genomic coordinates.
#' Return vector containing signal per nucleotide within a certain region.
#' 
#' @param gr Genomic Range containing metacolumn.
#' @param chr Chromosome (e.g. "chrI").
#' @param start Genomic coordinate to start get signal.
#' @param end Genomic coordinate to end get signal.
#'
#' @seealso \code{\link{coverage.plot}}
#'
#' @export

get.signal <- function(gr, chr, start, end) {
	if (!is.character(chr)) stop("\"end\" must be a character.")
	if (!is.numeric(start) | !is.numeric(end)) stop("\"start\" and \"end\" genomic coordinates must be numeric.")
	if (end < start) stop("\"end\" must be greater than \"start\" coordinate.")

    region <- GRanges(chr,IRanges(start,end))
    ovlp <- findOverlaps(region, gr)
    gvalue <- gr[subjectHits(ovlp)]
    border <- region[queryHits(ovlp)]
    start(gvalue[1]) = start(border[1])
    end(gvalue[length(gvalue)]) = end(border[length(gvalue)])
    vec <- rep(gvalue$score, width(gvalue))
    return(vec)
}
