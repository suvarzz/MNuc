#' All features averages plot
#'
#' @param indir Input directory containing bedgraph files.
#' @param outdir Output directory.
#' @param fdir Directory containing different features in files.
#' @param chromsizes Text file containing sizes of all chromosomes in columns (chromosome ~ size).
#' @param title Main title of the plot.
#' @param smooth If lines between points are smoothed.
#' @param log Boolean, if log difference of signal is needed.
#' @param nlog Id of sample in files in 'indir' which is used for log difference (e.g. 1 - first sample) 'log' must be TRUE, log.num must be NULL.
#' @param log.num Numeric to calculate log difference (e.g. 21.5). 'log' must be TRUE, nlog must be NULL.
#' @param filename Name of the output plot file
#' @param exclude_seq Optional argument to exclude chromosomes.
#'
#' @export

features.avg.plot <- function(indir,
							  outdir,
							  fdir,
							  chromsizes = 'sc',
							  lables = NULL,
							  title = "All features averages plot",
							  smooth = FALSE,
							  log = FALSE,
							  nlog = NULL,
							  log.num = NULL,
							  xlab = "Time",
							  ylab = "Signal intensity",
							  filename = "All_features_avg_plot",
							  exclude_seq = "chrM",
							  ylim = NULL)
{
  if (log & !is.null(nlog) & !is.null(log.num))
    stop("Please choose between \"nlog\"  and \"log.num\" arguments.\n")	
  
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
	cat("All features imported as genomicranges.\n")

	### Import Signals > GR
	cat("Importing signals started. It takes time.\n")
	signals <- list.files(path=indir, pattern="*.bdg.gz", full.names=TRUE, recursive=FALSE)
		
	scores_list <- lapply(signals, function(x) {
    	signal <- sort(GenomeInfoDb::sortSeqlevels(import(x, format="bedGraph")))
    	# Set seqlevels to make a correct coverage across all chromosome lengths
    	seqlengths(signal) <- seql
    	signal <- dropSeqlevels(signal, value=exclude_seq, pruning.mode="coarse")
    	# Find coverage of GRange
    	GenomicRanges::coverage(signal, weight="score")
	} )

	cat("Importing signals ended, scores list created without errors.\n")

	### Get Signal Names
	if (is.null(lables)) {
    	signal_names <- sapply(strsplit(tools::file_path_sans_ext(basename(signals)), "_"), function(x) x[2])
	} else { signal_names <- lables}   
	### START PLOT
	cat("Plot started")
	dir.create(file.path(outdir), recursive = TRUE)
	
	pdf(paste(outdir, filename, ".pdf", sep=""), width=7, height=13, pointsize=5)
	par(mfrow=c(6,4)) # how many diagrams on one plot?
	par(mar=c(4,4,2,1.5)) # margins size
	par(oma=c(4,1,4,1)) # outer margins in lines
	col=c('darkslateblue', 'brown4', 'grey')

	invisible(lapply(seq_along(feature_list), function(idx) {
    	vec_data <- sapply(scores_list, function(sc) {
    	    data <- GenomicRanges::binnedAverage(feature_list[[idx]], sc, "avg_score")
            m <- mean(data$avg_score)
    	})
    	
    	# ylim calculation
    	if (!log & is.null(ylim)) {
    	    lim_min <- min(vec_data)
    	    lim_max <- max(vec_data)
    	    ylim = c(lim_min, lim_max)
    	}
    	
    	mx <- matrix(vec_data, nrow=length(indir), byrow=TRUE)
    	x <- c(1:ncol(mx))
    	if (log & is.null(ylim)) {
    	    if (!is.null(nlog)) {
    	      limits <- apply(mx, 1, function(x) log2(x/x[nlog]))
    	      ylim <- c(min(limits), max(limits))
    	    }
    	    if (!is.null(log.num)) {
    	      limits <- apply(mx, 1, function(x) log2(x/log.num))
    	      ylim <- c(min(limits), max(limits))
    	    }

    	}
    	
    	plot(1,type = 'n',
    	     axes = T,
    	     xaxt = "n",
    	     xlab = xlab,
    	     ylab = ylab,
    	     ann = T, 
    	     ylim = ylim,
    	     las=1,
    	     xlim = c(1, ncol(mx)),
    	     main = feature_names[idx])
    	
        for (i in 1:nrow(mx)) {
            y <- mx[i,]
            if (log & !is.null(nlog)) y = log2(y/y[nlog])
            if (log & !is.null(log.num)) y = log2(y/log.num)
            points(x, y, col=col[i])
            if (smooth) {
            smoothingSpline = smooth.spline(x, y, spar=0.1)
            xl <- seq(min(x), max(x), length=length(x)*10)
            lines(predict(smoothingSpline, xl), col=col[i])
            } else {
                lines(x, y, col=col[i])
            }
        }
        axis(1, at=1:length(signal_names), labels=signal_names, srt = 45, xpd=TRUE)
	}))
	
	mtext(title, side=3, line=1, outer=TRUE, cex=1.5)
	mtext(Sys.Date(), side=1, line=1, outer=TRUE, adj=1)
	dev.off()
}
