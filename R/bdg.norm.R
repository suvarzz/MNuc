#' Normalize bed graph files
#' 
#' @details Normalize bed graph files
#' 
#' @param indir Directory containing bedgraph files with peaks
#' @param outdir Output directory to save bedgraph files after normalization
#' @param normk Vector of lengths of files containing multipliers for normalization
#' 
#' @return NULL
#' @export

bdg.norm <- function(indir,
                     outdir,
                     normk=NULL)
{
    files <- list.files(path=indir, pattern="\\.bdg.gz$", full.names=TRUE, recursive=FALSE)
    
    if (length(files)!=length(normk)) stop("Length of normk must be equal to number of files to normalize")
    if (is.null(normk)) stop("Please specify nomralizing multipliers.")
    
    silent <- mclapply(1:length(files), function(i) {
        # Import bedGraph to GRanges
        gr <- rtracklayer::import(files[i], format="bedGraph")
        gr$score <- (gr$score)*normk[i]
        
        # Export bedgraph file
        rtracklayer::export(gr, con=paste(outdir, basename(files[i]), sep=""), format="bedGraph" )
    })
}
    