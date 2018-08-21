#' All features averages plot
#'
#' @param indir Input directory containing bedgraph files.
#' @param outdir Output directory.
#' @param fdir Directory containing different features in files.
#' @param chromsizes Text file containing sizes of all chromosomes in columns (chromosome ~ size).
#' @param title Main title of the plot
#' @param filename Name of the output plot file
#' @param notes Optional notes, which is displayed in the bottom of the plot.
#' @param exclude_seq Optional argument to exclude chromosomes.
#'
#' @export

features.box.plot <- function(indir,
							  outdir,
							  fdir,
							  chromsizes='sc',
							  lables = NULL,
							  title="All features averages plot",
							  filename="All_features_avg_plot",
							  exclude_seq="chrM",
							  ylim=NULL)
{
	seql <- MNuc::seqlevels.from.chrsizes(chromsizes=chromsizes)

	### Import Features > GR
	features <- list.files(path=fdir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
	feature_list <- lapply(features, function(f) {
    gr <- makeGRangesFromDataFrame(read.table(f, sep="\t", header=TRUE, quote =""), keep.extra.columns=TRUE, 
                                   seqinfo=seql)
    #seqlengths(gr) <- seql
    gr <- dropSeqlevels(gr, value=exclude_seq, pruning.mode="coarse")
	} )

	feature_names <- lapply(features, function(f) { tools::file_path_sans_ext(basename(f)) })
	cat("All features imported as genomicranges")

	### Import Signals > GR
	cat("Importing signals started")
	signals <- list.files(path=indir, pattern="*.bdg.gz", full.names=TRUE, recursive=FALSE)
		
	scores_list <- lapply(signals, function(x) {
    	signal <- sort(GenomeInfoDb::sortSeqlevels(import(x, format="bedGraph")))
    	# Set seqlevels to make a correct coverage across all chromosome lengths
    	seqlengths(signal) <- seql
    	signal <- dropSeqlevels(signal, value=exclude_seq, pruning.mode="coarse")
    	# Find coverage of GRange
    	GenomicRanges::coverage(signal, weight="score")
	} )

	cat("Importing signals ended, scores list created without errors")

	### Get Signal Names
	if (is.null(lables)) {
    	signal_names <- sapply(strsplit(tools::file_path_sans_ext(basename(signals)), "_"), function(x) x[2])
	} else { signal_names <- lables}   
	### START PLOT
	cat("Plot started")
	dir.create(file.path(outdir), recursive = TRUE)
	
	pdf(paste(outdir, filename, ".pdf", sep=""), width=7, height=13, pointsize=5)
	par(mfrow=c(6,3)) # how many diagrams on one plot?
	par(mar=c(4,4,2,1.5)) # margins size
	par(oma=c(4,1,4,1)) # outer margins in lines
	col=c('darkslateblue', 'brown4')

	invisible(lapply(seq_along(feature_list), function(idx) {
    	vec_data <- sapply(scores_list, function(sc) {
    	    data <- GenomicRanges::binnedAverage(feature_list[[idx]], sc, "avg_score")
            m <- mean(data$avg_score)
    	})
    	
    	lim_min <- min(vec_data)
    	lim_max <- max(vec_data)
    	ylim = c(lim_min, lim_max)
    	
    	mx <- matrix(vec_data, length(indir), byrow=TRUE)
    	x <- c(1:ncol(mx))
    	plot(0,type='n',
    	     axes=T,
    	     xaxt="n",
    	     ann=T, 
    	     ylim=ylim, 
    	     xlim=c(1, ncol(mx)),
    	     main=feature_names[idx])
    	
        for (i in 1:nrow(mx)) {
            y <- mx[i,] 
            points(x, y, col=col[i])
            smoothingSpline = smooth.spline(x, y, spar=0.1)
            xl <- seq(min(x), max(x), length=length(x)*10)
            lines(predict(smoothingSpline, xl), col=col[i])
        }
        axis(1, at=1:length(signal_names[x]), labels=signal_names[x], srt = 45, xpd=TRUE)
	}))
	
	mtext(title, side=3, line=1, outer=TRUE, cex=1.5)
	mtext(Sys.Date(), side=1, line=1, outer=TRUE, adj=1)
	dev.off()
}
