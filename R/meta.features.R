#############################################################################################
### META_TSS

#' TSSs regions for metaplot
#' 
#' @details Create a csv.gz file containing saved Genomic Range of regions around TSSs (transcriptional start sites).
#'
#' @param outdir Output directory.
#' @param upstream Number of nucleotides upstream from TSSs.
#' @param downstream Number of nucleotides downstream from TSSs (transcriptional start sites).
#' @param chrom_sizes File with sizes of chromosomes (.chrom.sizes).
#' 
#' @return None
#'
#' @seealso \code{\link{meta.ars}}
#'
#' @export

meta.tss <- function (outdir, upstream, downstream, chrom_sizes, exclude_chrom = NULL)
{
    if (is.null(outdir))
        stop("Please specify \"outdir\" output directory.\n")
    if (is.null(chrom_sizes))
        stop("Please specify \"chrom_sizes\" chromosome sizes file.\n")
    # Set chromosomes
    chrom.sizes <- read.table(chrom_sizes, sep="\t", header=FALSE, quote="", col.names=c("chromosome", "length"))
    chroms <- as.character(chrom.sizes$chromosome) # list of chromosomes from chrom.sizes file
    
    if (!all(exclude_chrom %in% chroms))
        stop("Specified \"exclude_chrom\" chromosome does not exists. Check names of chromosomes.\n")
    chroms <- chroms[! chroms %in% exclude_chrom] # exclude chromosome
    
    stopifnot(is.numeric(upstream) & is.numeric(downstream))
    
win = upstream+downstream
txdb = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
seqlevels(txdb) <- chroms
## Choose promoter regions and out-of-bound filtering
csp <- suppressWarnings(trim(promoters(txdb, upstream, downstream)))
## Discard trancated promoter regions at the ends of chromosomes (csp < win)
csp <- csp[width(csp) == win]

# Add csp into genomic ranges list
grl <- GenomicRanges::GRangesList(csp)
names(grl) <- paste("TSS_meta-", upstream, "-", downstream, sep="") 

# save genomic ranges list as tsv file
MNuc::saveGRlist(grl, outdir=outdir)
}

#############################################################################################
### META_ARS

#' ARSs regions for metaplot
#' 
#' @details Create a csv.gz file containing saved Genomic Range of regions around ARSs.
#'
#' @param ars_csv Text file containing ranges.
#' @param outdir Output directory.
#' @param win Windows around ARS's centers.
#' @param chrom_sizes File with sizes of chromosomes (.chrom.sizes).
#'
#' @usage 
#' meta.ars <- function (ars_csv = "/home/suvar/Projects/Sources/Cs_genomic_features_v1/ARS.csv.gz",
#'     outdir = "/home/suvar/Projects/Sources/Cs_genomic_features_v1/meta/",
#'     win = 4000,
#'     chrom_sizes = "/home/suvar/Projects/Sources/sacCer3.chrom.sizes")
#' 
#' @return None
#'
#' @seealso \code{\link{meta.tss}}
#'
#' @export

meta.ars <- function (ars_csv,
                      outdir,
                      win = 2000,
                      chrom_sizes)
{

    ## SET PARAMETERS
    chrom.sizes <- read.table(chrom_sizes, sep="\t", header=FALSE, quote="", col.names=c("chromosome", "length"))
    seql = as.numeric(chrom.sizes$length)
    names(seql)=as.character(chrom.sizes$chromosome)

    ars <- makeGRangesFromDataFrame(read.table(ars_csv, sep="\t", header=TRUE, quote =""), keep.extra.columns=TRUE, seqinfo=seql)
    # find the center of ARSs
    start(ars) <- (start(ars) + (end(ars)-start(ars))/2)
    end(ars) <- start(ars)
    # Make a rages around the centers of ARSs
    ars <- suppressWarnings(trim(ars + win/2))
    # Remove one nucleotide from the left to get width equal to win
    ars <- narrow(ars, start=2) 
    # Filter for ARSs with the length of win
    ars <- ars[width(ars) == win]

    # Add csp into genomic ranges list
    grl <- GenomicRanges::GRangesList(ars)

    names(grl) <- paste("ARS_meta-win-", win, sep="") 

    # save genomic ranges list as csv file
    MNuc::saveGRlist(grl, outdir=outdir)
}
