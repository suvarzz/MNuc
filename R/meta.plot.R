#' @title TSS & ARS regions for metaplot
#' 
#' These functions allows to create txt file with a list of TSS or ARS with a given size.
#' 
#' @rdname meta.plot
#' 
#' @param dir Input directory with DAT files (dat.gz).
#' @param legend_names A vector with legend names.
#' @param main Main header of a plot.
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

meta.plot <- function (indir,
                       outdir,
                       col.file,
                       file_name="Meta_plot",
                       legend_names,
                       lines_to_draw=NULL,
                       main="Meta Plot",
                       xlab="Coordinates",
                       ylab="Signal",
                       xlim=NULL,
                       ylim=NULL,
                       las=1,		 # labels orientations (all horizontal)
                       cex.main=1, # font size for titel
                       cex.lab=1,	 # font size for labels
                       cex.axis=1, # font size for axis
                       cex=1) # font size legend
{
    if (is.null(indir) | is.null(outdir)) stop("Specify \"indir\" and \"outdir\" input and output directories.")
    if (is.null(legend_names)) stop("Specify legend names")
    
    ## READ FILES, get data and keep the first (coordinates) and the last (averages) of data
    
    files <- list.files(path=indir, pattern="*.dat.gz", full.names=T, recursive=FALSE)
    data <- lapply(1:length(files), function(f) {
        df <- read.table(files[f], header=TRUE, sep="\t", na.strings="NA")
    })
    if (!all(apply(sapply(data, dim), 1, function(x) length(unique(x)) == 1) == TRUE)) stop("Dimentions of data in files are different.")
    if (is.null(lines_to_draw)) lines_to_draw <- 1:(ncol(data[[1]])-1)
        
    ## Read colors from file and take the last colors
    colors <- as.vector(read.table(col.file)[,1])
    colors <- tail(colors, n=(ncol(data[[1]])-1))
    
    ## GRAPHICAL PARAMETERS
    ## get coordinates and values limits
    if(is.null(xlim) & is.null(ylim)) {
        all_values <- unlist(lapply(data, '[', 2:ncol(data[[1]])))
        all_coord <- unlist(lapply(data, '[', 1 ))
        xlim <- c(min(all_coord), max(all_coord))
        ylim <- c(min(all_values), max(all_values))
    }
    
    ##### Draw plot for all files
    dir.create(outdir, recursive=T)
    rows=ceiling(length(files)/2)
    cols=ceiling(length(files)/rows)
    pdf(paste(outdir, file_name, ".pdf", sep=""), width=3*cols, height=3*rows, pointsize=5)
    par(mfrow=c(rows, cols))
    
    for (f in 1:length(data)) {
        ## Empty plot
        plot(1, type = "n",
             xlim = xlim,
             ylim = ylim,
             main=main,
             xlab=xlab,
             ylab=ylab,
             las=las,
             cex.main=cex.main,
             cex.lab=cex.lab,
             cex.axis=cex.axis)
    
        for (l in lines_to_draw) {
            lines(data[[f]][,1], data[[f]][,(l+1)], col=colors[l])
        }
        
        ##### LEGEND
        legend("bottomright", legend_names,
               ncol=ifelse(length(lines_to_draw) > 5, 2, 1),
               lty=rep(1, times=length(lines_to_draw)),	# line type
               lwd=rep(2, times=length(lines_to_draw)),	# line width
               col=colors[lines_to_draw],
               cex=cex)	# text size
    
        dev.off()
    }
}
