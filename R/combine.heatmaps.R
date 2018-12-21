#' Combine heatmaps png into one image
#' 
#' @param indir Directory containing png heatmaps.
#' @param outdir Output directory.
#' 
#' @return NULL
#' @export

combine.heatmaps <- function(indir, outdir)
{
    files <- list.files(path=indir, pattern="\\.png$", full.names=TRUE, recursive=FALSE)

    ims <- lapply(files, function(f) { readPNG(f) })



    plotheight = 1200*length(files)
    plotwidth = max(unlist(lapply(ims, function(im) dim(im)[2])))

    svg(paste(outdir, "combined_heatmap.svg", sep=""), width=100, height=length(files)*10)

    plot(1,
         bty="n",
         type="n",
         xaxt = "n",	# remove x axis
         yaxt = "n",	# remove y axis
         xlab="",
         ylab="",
         xlim=c(1, plotwidth), 
         ylim=c(1, plotheight))

    ytop=plotheight

    for (i in 1:length(ims)) {
        rasterImage(ims[[i]], xleft=1, xright=dim(ims[[i]])[2], ybottom=ytop-1000, ytop=ytop)
        ytop = ytop - 1200
    }
    dev.off()
}
