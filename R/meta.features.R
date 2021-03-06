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
names(grl) <- paste("TSS_meta_", upstream, "_", downstream, sep="") 

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
                      filename = "ARS_meta",
                      win = 2000,
                      chromsizes = 'sc')
{
    seql = MNuc::seqlevels.from.chrsizes(chromsizes=chromsizes) 

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

    # save genomic range as a csv file
    MNuc::saveGR(gr=ars, outdir=outdir, filename=paste(filename, "_win_", win, sep=""))
}

###############################################################################
### Split intergenic

#' Split intergenic regions and save as csv genomic ranges:
#'     intergenic - all intergenics
#'     convergent - convergent intergenics
#'     devergent - devergent intergenics
#'     direct - direct intergenics (strand '+')
#'     reverse - reverse intergenics (strand '-')
#'     direct_reverse - combined direct and reverse intergenics
#' 
#' @details Create a csv.gz file containing saved Genomic Range of split intergenic regions.
#'
#' @param outdir Output directory.
#'
#' @return None
#'
#' @export
split_intergenic <- function(outdir) {
    # Get data from TxDb and BSgenome databases
    txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
    
    # Get intergenic regions
    genes <- unlist(range(transcriptsBy(txdb, by='gene')))
    genic <- reduce(genes, ignore.strand=TRUE)
    ir <- gaps(genic)
    ir <- ir[strand(ir) == '*'] # remove 2 of each whole chromosome ('+' & '-') ranges
    
    # Add names intergenics
    internames <- sapply(seq(1:length(ir)), function(x) {paste("intergenic_", x, sep="")} )
    mcols(ir)$names <- internames
    
    # Create logical vectors
    convergent <- divergent <- direct <- reverse <- logical(length = length(ir))
    
    cat("Split intergenic started. It takes time.\n")
    for (i in 1:length(ir)) {
        sn <- as.character(seqnames(ir[i]))
        st <- start(ir[i])-1
        ed <- end(ir[i])+1
        #left gene strand
        lgs <- as.vector(strand(genes[seqnames(genes) == sn & end(genes) == st]))
        # right gene strand
        rgs <- as.vector(strand(genes[seqnames(genes) == sn & start(genes) == ed]))
        
        if(length(lgs) !=0 & length(rgs) !=0) {  
            
            if (lgs == '+' & rgs == '+') { direct[i] <- 'TRUE'
            } else if (lgs == '-' & rgs == '-') { reverse[i] <- 'TRUE'
            } else if (lgs == '+' & rgs == '-') { convergent[i] <- 'TRUE'
            } else if (lgs == '-' & rgs == '+') { divergent[i] <- 'TRUE' } }
        
        # remove variables befor the next cycle
        rm(sn, st, ed, lgs, rgs)
    }
    
    # Select intergenic regions
    gr_divergent <- ir[divergent=='TRUE']
    gr_convergent <- ir[convergent=='TRUE']
    gr_direct <- ir[direct=='TRUE']
    gr_reverse <- ir[reverse=='TRUE']
    
    strand(gr_direct) <- "+"
    strand(gr_reverse) <- "-"
    
    gr_direct_reverse <- c(gr_direct, gr_reverse)
    gr_direct_reverse <- sortSeqlevels(gr_direct_reverse)
    gr_direct_reverse <- sort(gr_direct_reverse, ignore.strand=TRUE)
    
    # WRITE DATA INTO FILES
    # list of genomicranges
    grl <- GenomicRanges::GRangesList(INTERGENIC=ir,
                                      INTERGENIC_CONVERGENT=gr_convergent, 
                                      INTERGENIC_DIVERGENT=gr_divergent,
                                      INTERGENIC_DIRECT=gr_direct,
                                      INTERGENIC_REVERSE=gr_reverse,
                                      INTERGENIC_DIRECT_REVERSE=gr_direct_reverse)
    
    MNuc::saveGRlist(grl, outdir=outdir)
}
