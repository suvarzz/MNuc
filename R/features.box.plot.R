#' All features box plot
#'
#' @param indir Input directory containing bedgraph files.
#' @param outdir Output directory.
#' @param fdir Directory containing different features in files.
#' @param chromsizes Text file containing sizes of all chromosomes in columns (chromosome ~ size).
#' @param title Main title of the plot
#' @param filename Name of the output plot file
#' @param notes Optional notes, which is displayed in the bottom fo the plot.
#' @param exclude_seq Optional argument to exclude chromosomes.
#'
#' @export

features.box.plot <- function(indir,
							  outdir,
							  fdir,
							  chromsizes,
							  title="Average log2 difference of H4K16Ac vs 1.5h sample",
							  filename="All_features_box_plot",
							  notes="",
							  exclude_seq="chrM")
{
	seql <- MNuc::seqlevels.from.chrsizes(chromsizes)

	### Import Features > GR
	features <- list.files(path=fdir, pattern="*.csv", full.names=TRUE, recursive=FALSE)
	feature_list <- lapply(features, function(f) {
    gr <- makeGRangesFromDataFrame(read.table(f, sep="\t", header=TRUE, quote =""), keep.extra.columns=TRUE, 
                                   seqinfo=seql)
    #seqlengths(gr) <- seql
    gr <- dropSeqlevels(gr, value=exclude_seq, pruning.mode="coarse")
	} )

	feature_names <- lapply(features, function(f) { tools::file_path_sans_ext(basename(f)) })
	print("All features imported as genomicranges")

	### Import Signals > GR
	print("Importing signals started")
	signals <- list.files(path=indir, pattern="*.bdg.gz", full.names=TRUE, recursive=FALSE)
	scores_list <- lapply(signals, function(x) {
    	signal <- sort(GenomeInfoDb::sortSeqlevels(import(x, format="bedGraph")))
    	# Set seqlevels to make a correct coverage across all chromosome lengths
    	seqlengths(signal) <- seql
    	signal <- dropSeqlevels(signal, value=exclude_seq, pruning.mode="coarse")
    	# Find coverage of GRange
    	GenomicRanges::coverage(signal, weight="score")
	} )

	print("Importing signals ended, scores list created without errors")

	### Get Signal Names
	signal_names <- sapply(strsplit(tools::file_path_sans_ext(basename(signals)), "_"), function(x) x[2])

	### START PLOT
	print("Plot started")
	dir.create(file.path(outdir), recursive = TRUE)
	pdf(paste(outdir, filename, ".pdf", sep=""), width=7, height=13, pointsize=5)
	par(mfrow=c(6,3)) # how many diagrams on one plot?
	par(mar=c(4,4,2,1.5)) # margins size
	par(oma=c(4,1,4,1)) # outer margins in lines
	col=rep(c('aliceblue','mistyrose'), 11)
	bor=rep(c('darkslateblue', 'brown4'), 11)
	lim_min <- -0.5    # min(unlist(list_data))
	lim_max <- 0.5    # max(unlist(list_data))

	invisible(lapply(seq_along(feature_list), function(idx) {
    	list_data <- lapply(scores_list, function(sc) {
    	    data <- GenomicRanges::binnedAverage(feature_list[[idx]], sc, "avg_score")
    	    data$avg_score
    	})
    	#lim_min <- min(unlist(list_data))
    	#lim_max <- max(unlist(list_data))
    	boxplot(list_data, 
    	        ylim=c(lim_min, lim_max), 
    	        main=feature_names[idx],
    	        las=1,
    	        xaxt="n",
    	        col=col,
    	        border=bor,
    	        outline = FALSE) # do not draw outliers
    	axis(1, at=1:length(signal_names), labels=signal_names, srt = 45, xpd=TRUE)
	}))
	mtext(title, side=3, line=1, outer=TRUE, cex=1.5)
	mtext(notes, side=1, line=1, outer=TRUE, adj=0)
	mtext(Sys.Date(), side=1, line=1, outer=TRUE, adj=1)
	dev.off()
}
