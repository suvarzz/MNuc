#' Coverage plot
#'
#' @param indir Directory containing bedgraph files with signals (e.g. after MACS2).
#' @param outdir Output directory.
#'
#' @export

coverage.plot <- function(indir,
						  outdir,
						  names = c("HU 0 min","","","","","","","","","HU 90 min","NOC 0 min","","","","","","","","","NOC 90 min"),
						  num_files = c(1,10,11,20),
						  y_lim = 30,
						  chromosome="chrIII",
						  start=1,
						  end=30000)
{
width=end-start

## Read files
files <- list.files(path=indir, pattern="*.bdg", full.names=T, recursive=FALSE)	# pattern="HU.*\\.bdg"

## Get genes for plot
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene
filter <- list(tx_chrom = c(chromosome))
cds <- cds(txdb, columns=c("cds_id","gene_id"), filter=filter)
g <- GRanges(seqname=c(chromosome), ranges=IRanges(start=start, width=width))
o <- findOverlaps(g, cds)
    gval <- cds[subjectHits(o)]
    board <- g[queryHits(o)]
	if (start(gval[1]) < start) { start(gval[1]) = start }
	if (end(gval[length(gval)]) > end ) { end(gval[length(gval)]) = end }

## How many plots will be drawn
pl <- length(num_files)+1

## Draw plot for all files
svg(paste(outdir, "plot.svg", sep=""), width=10, height=pl*1)	# height is 1 cm per plot

par(mfrow=c(pl,1),
	bty = 'n',	# remove box from plot
	mai=c(0.1, 0.3, 0.1, 1))	# margins

## Draw gene models
	plot(0, col="white", xaxt = "n", yaxt = "n", ylim=c(0,y_lim), xlim=c(0,width), xlab="", ylab="" )
	abline(h = y_lim/2, col="black")
	for(gene in 1:length(gval)) {
	rect(start(gval[gene])-start, y_lim/2-2, end(gval[gene])-start, y_lim/2+2, col="black")
	}
## Draw scaling bar
	segments(5000, 0.5, 10000, 0.5)	# (x0,y0,x1,y1)
	segments(5000, 0,5000,1)
	segments(10000,0, 10000,1)
	text(7500, 3, labels="5 kb")
	

for (f in 1:length(files)) {

    if (f %in% num_files) {

	gr <- import(files[f], format="bedGraph")
	vec <- MNuc::get.signal(gr=gr, chr=chromosome, start=start, end=end) # maybe set score is necessary.

	plot(vec,
		xaxt = "n",	# remove x axis
		yaxt = "n",	# remove y axis
		ylim=c(0,y_lim),	# custom y-limit
		xlab="",
		ylab="",
		type="l",
		col="black",
		lwd=0.5,	# line width
		)	# lables orientation, all horizontal
	segments(0,0,width,0)
	axis(side=2,
		at=c(0,y_lim),
		tick=T,
		las=1)
	text(x=width*1.05, y=15, labels=names[f], xpd=NA)	# names of plots on the rigt sides

	}
 }

dev.off()

}
