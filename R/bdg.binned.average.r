#' Binned average from bedgraph file
#' 
#' @details Calculate binned average from bedgraph files
#' 
#' @param indir Directory containing bedgraph files
#' @param outdir Output directory
#' @param binsize Length of one bin (default - 100 bp)
#' 
#' @return NULL
#' @export

bdg.binned.average <- function(indir,
                               outdir,
                               binsize = 100)
{
    # Create output directory
    dir.create(file.path(outdir), recursive = TRUE)

    bins <- tileGenome(seqinfo(BSgenome.Scerevisiae.UCSC.sacCer3), tilewidth=binsize, cut.last.tile.in.chrom=TRUE)

    files <- list.files(path=indir, pattern="\\.bdg.gz$", full.names=TRUE, recursive=FALSE)
    invisible <- lapply(files, function(f) {

        # Import bedGraph to GRanges
	    signal <- import(f, format="bedGraph")

        # Sort GRange
        signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
    
        # Find coverage of GRange
        score <- GenomicRanges::coverage(signal, weight="score")
    
        # Calculate binned average
        binned_data <- GenomicRanges::binnedAverage(bins, score, "average_score")
    
        # Vector of average score
        bin_avg <- mcols(binned_data)$average_score

        # Write average score vector into a file
        wfile <- gzfile(paste(outdir, tools::file_path_sans_ext(basename(f)), ".bins.gz", sep=""), "w")
        write.table(bin_avg, file=wfile, col.names = F, row.names = F, append = F, na= "NA")
        close(wfile)
        rm(bin_avg, binned_data, score, signal)
	})
}