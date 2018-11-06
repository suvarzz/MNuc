#' @title Create list of matrices from several experiments.
#'
#' Create list of matrices, where length of list is number of files (splits) in directories.
#' For each matrix columns contains timepoints and rows contain experiments.
#'
#' @param indir Input directories. Each directory correspont to one experiment.
#' @param expnames Names of experiments. Names will be detected from names of directories if not specified.
#' @param timepoints Datasets in columns of each file. Names of datapoints will be detected from column names of the first file.
#' @param splitnames Names of files in directories. Splitnames will be detected from names of files in the first directory.
#' 
#' @author Mark Boltengagen \email{m.boltengagen@@gmail.com}
#' 
#' @return List of matrices
#' 
#' @export
mxl <- function(indir,
               expnames = NULL,
               timepoints = NULL,
               splitnames = NULL) {
    
    if (is.null(indir)) stop("Specify \"indir\" input directories.")
    nl <- length(indir) # length of list
    
    if (is.null(timepoints)) {
        # TODO detect names of columns
        # TODO check if all the same in all directories -> warning if different
        
    }
    
    if (is.null(expnames)) {
        # TODO detect names from the first file in a directory
        # TODO compair names in all directories -> warning if different
    }
    
    if (length(expnames) != length(indir)) stop("Number of experiments must be equal to number of input directories.")
    
    if (is.null(splitnames)) {
        # TODO detect names from names of files in directories
        # TODO compair names in all directories -> warning if different
    }
    
    ncol <- length(timepoints)
    nrow <- length(expnames)
    mx <- matrix(0, nrow = nrow, ncol = ncol)
    colnames(mx) <- timepoints
    rownames(mx) <- expnames
    
    mxl <- lapply(1:length(splitnames), function(x) mx)
    names(mxl) <- splitnames
    return(mxl)
    
}

#' @title Dot plot from list of matrices.
#'
#' Take list of matrices and plot each matrix on separate figure. Legends are in name of matrices.
#'
#' @param mxl List of matrices.
#' @param outdir Output directory.
#' @param filename Name of output file.
#' @param ... Optional graphical parameters.
#' 
#' @author Mark Boltengagen \email{m.boltengagen@@gmail.com}
#' 
#' @return None
#' 
#' @export

mxl.dotplot <- function(mxl,
                        outdir, 
                        filename, 
                        title, ...) {

    dir.create(file.path(outdir), recursive = TRUE)
    
    ylim <- c(min(unlist(mxl)), max(unlist(mxl)))
    
    pdf(paste(outdir, filename, ".pdf", sep=""), width=12, height=3, pointsize=5)
    col=c('black', 'red', 'blue')
    par(mfrow=c(1,length(mxl)))
    par(mar=c(4,4,2,1.5)) # margins size
    par(oma=c(4,1,4,1)) # outer margins in lines
    
    for(f in 1:length(mxl)) {
        mx <- mxl[[f]]
        x <- 1:ncol(mx)
        # Empty plot
        plot(1,type = 'n',
             axes = T,
             xaxt = "n",
             ylim = ylim, 
             xlim = c(1, ncol(mx)),
             ann = T,
             main = names(mxl[f]))
        
        # Draw 
        for (i in 1:nrow(mx)) {
            y <- mx[i,]
            # Points
            points(x, y, col=col[i])
            
            # Smoothed lines
            smoothingSpline = smooth.spline(x, y, spar=0.1)
            xl <- seq(min(x), max(x), length=length(x)*10)
            lines(predict(smoothingSpline, xl), col=col[i])
        }
        axis(1, at=1:ncol(mx), labels=colnames(mx))
    }
    
    mtext(title, side=3, line=1, outer=TRUE, cex=1.5)
    dev.off()
}