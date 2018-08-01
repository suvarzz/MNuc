#' @title TSS & ARS regions for metaplot
#' 
#' These functions allows to create txt file with a list of TSS or ARS with a given size.
#' 
#' @rdname meta.plot
#' 
#' @param indir Input directory with means DAT files (mean.dat.gz).
#' @param outdir Output directory.
#' @param col.file One column text file containing colors for curves.
#' @param filename Name of output file.
#' @param legend_names A vector containing legend names.
#' @param ncur Numeric vector of lines to draw (e.g  c(2,4,6)). The
#'   first element in \code{ncur} is a denominator if \code{log} argument is TRUE.
#' @param log Log2 of difference. The denominator is the first element in \code{ncur} argument.  
#' @param main Main header of the plot.
#' @param xlab downstrem nucleotides from feature.
#' @param col.file File containing list of colors in the first column
#' @param ... Optional graphics arguments
#' 
#' @return None
#' 
#' @seealso \code{\link{meta.data}}, \code{\link{meandata}}
#' 
#' @export
#' 

### META_PLOT
# TODO ncur subset names from legend names for plot. And legend_names set as in mean dat file. rename legend_names to curnames
# TODO subset colors from the file in the same way using ncur
# TODO make legend_names (curnames) by default the names of columns in dat files. in each cycle? or check if identical in all files?
# FIXME fix ylim if log is TRUE
meta.plot <- function (indir,
                       outdir,
                       col.file,
                       filename="Meta_plot",
                       legend_names,
                       ncur=NULL,
                       log=FALSE,
                       main="Meta Plot",
                       xlab="Coordinates",
                       ylab="Signal",
                       xlim=NULL,
                       ylim=NULL,
                       legendpos=c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center"),
                       las=1,		 # labels orientations (all horizontal)
                       cex.main=1, # font size for titel
                       cex.lab=1,	 # font size for labels
                       cex.axis=1, # font size for axis
                       cex=1)     # font size legend
{
    if (is.null(indir) | is.null(outdir)) stop("Specify \"indir\" and \"outdir\" input and output directories.")
    if (is.null(legend_names)) stop("Specify legend names")
    legendpos <- match.arg(legendpos)
    ## READ FILES, get data and keep the first (coordinates) and the last (averages) of data
    
    files <- list.files(path=indir, pattern="*.dat.gz", full.names=T, recursive=FALSE)
    data <- lapply(1:length(files), function(f) {
        df <- read.table(files[f], header=TRUE, sep="\t", na.strings="NA")
    })
    if (!all(apply(sapply(data, dim), 1, function(x) length(unique(x)) == 1) == TRUE)) stop("Dimentions of data in files are different.")
    if (is.null(ncur)) ncur <- 1:(ncol(data[[1]])-1)
        
    ## Read colors from file and take the last colors
    colors <- as.vector(read.table(col.file)[,1])
    colors <- tail(colors, n=(ncol(data[[1]])-1))
    
    ## GRAPHICAL PARAMETERS
    ## get coordinates and values limits
    if(is.null(xlim)) {
        all_coord <- unlist(lapply(data, '[', 1 ))
        xlim <- c(min(all_coord), max(all_coord))
    } 
    
    if (is.null(ylim) & log==FALSE) {
        all_values <- unlist(lapply(data, '[', 2:ncol(data[[1]])))
        ylim <- c(min(all_values), max(all_values))
    }
    
    ##### Draw plot for all files
    dir.create(outdir, recursive=T)
    rows=ceiling(length(files)/2)
    cols=ceiling(length(files)/rows)
    pdf(paste(outdir, filename, ".pdf", sep=""), width=3*cols, height=3*rows, pointsize=5)
    par(mfrow=c(rows, cols))
    
    for (f in 1:length(data)) {
        ## Empty plot
        plot(1, type = "n",
             xlim=xlim,
             ylim=ylim,
             main=main,
             xlab=xlab,
             ylab=ylab,
             las=las,
             cex.main=cex.main,
             cex.lab=cex.lab,
             cex.axis=cex.axis)
    
        for (l in ncur) {
            if(log==FALSE) {
                lines(data[[f]][,1], data[[f]][,(l+1)], col=colors[l], ...) }
            if(log==TRUE ) {
                # log2 of coordinates, draw log difference (line(l) - line(1))
                lines(data[[f]][,1], log2(data[[f]][,(l+1)] / data[[f]][,(ncur[1]+1)]), col=colors[l]) }
        }
        
        ##### LEGEND
        legend(legendpos, legend_names,
               ncol=ifelse(length(ncur) > 5, 2, 1),
               lty=rep(1, times=length(ncur)),	# line type
               lwd=rep(2, times=length(ncur)),	# line width
               col=colors[ncur],
               cex=cex)	# text size
    
        dev.off()
    }
}
