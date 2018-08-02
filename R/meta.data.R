#' Meta data generator
#' @description Generate meta data for metaplots.
#' 
#' @param meta.csv csv file with meta features.
#' @param bdg.data Directory with bedgraph files of signal peaks (e.g. after MACS2 peak calling).
#' @param upstream Upstream nucleotides from feature.
#' @param downstream Downstrem nucleotides from feature.
#' @param chromsizes Text file containing sizes of all chromosomes in columns (chromosome ~ size).
#' @param threads Optional parameter of number of working processors. If not set, all-1 will be taken.
#' 
#' @return TSS regions as CSV file dat.gz
#' @export

meta.data <- function(meta.csv, bdg.data, outdir, upstream, downstream, chromsizes, threads=NULL)
{
    if (is.null(bdg.data))
        stop("Please specify \"bdg.data\" bedgraph data containing signals")
    
	files <- list.files(path=bdg.data, pattern="*.bdg.gz", full.names=TRUE, recursive=FALSE)
    
    if (is.null(outdir))
        stop("Please specify \"outdir\" output directory.\n")

    # Create output directory
    dir.create(file.path(outdir), recursive=T)

	seql <- MNuc::seqlevels.from.chrsizes(chromsizes)

    if (!is.numeric(upstream) & !is.numeric(downstream))
        stop("Upstream region \"upstream\" and downstream region \"downstream\" must be numeric.\n")

    win = upstream+downstream

    if (is.null(meta.csv))
        stop("Please specify \"meta.csv\" output directory.\n")

    ftr <- makeGRangesFromDataFrame(read.table(meta.csv, sep="\t", header=TRUE, quote =""), keep.extra.columns=TRUE, seqinfo=seql)

	if (!all(width(ftr)==win))
        stop("Width of ranges from \"meta.csv\" is not equal to window width \"upstream\" + \"downtream\".")

	if (!is.null(threads) & is.numeric(threads)) {
        numcores = threads
    } else { numcores=parallel::detectCores() - 1 }

    # Write dat files and return NULL
    silent <- mclapply(files, function(f) {
        gr <- import(f, format="bedGraph")
	    ## Create a matrix to store score data for all windows
	    mx <- matrix(NA, win, length(ftr)+1) 
	    # Set names of mx columns 
	    colnames(mx) <- c(1:(length(ftr)+1))
	    # first column of mx contains genomic coordinates
	    colnames(mx)[1] <- "Coord"
	    # Set coordinates
	    mx[,1] <- seq(from = -1*upstream+1, to = downstream, by = 1)
        for (i in 1:length(ftr)) {
	        ftr.i <- ftr[i]
	        ovlp <- findOverlaps(ftr.i, gr)
	  
	        ## Remove no overlaps (when reads at the end of chromosomes do not exist)
	        if (length(ovlp) > 0 ) {
	            gvalue <- gr[subjectHits(ovlp)]
	            border <- ftr.i[queryHits(ovlp)]
	            start(gvalue[1]) = start(border[1])
	            end(gvalue[length(gvalue)]) = end(border[length(gvalue)])
	            vec <- rep(gvalue$score, width(gvalue))
	    
	            if (length(vec) == win) {
                    if (as.logical(strand(ftr[i]) == "-")) {vec <- rev(vec)}
	                mx[,i+1] <- vec
	                # the last column in feature genomic range is a name of gene/ars
	                colnames(mx)[i+1] <- as.character(mcols(ftr)[[length(mcols(ftr))]][i])
	            }
            }
        }
	    # Write dat file
	    name.1 <- sub('\\.bdg$', '', tools::file_path_sans_ext(basename(f)))
	    name.2 <- sub('\\.csv$', '', tools::file_path_sans_ext(basename(meta.csv)))
	    wfile <- gzfile(paste(outdir, name.1, "_", name.2, ".dat.gz", sep=""), "w")
	    write.table(mx, file=wfile, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
	    close(wfile)
	    rm(mx, gr)
	}, mc.cores=numcores)
}
