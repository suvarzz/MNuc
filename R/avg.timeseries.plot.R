#' Calculate matrix of averages from meandat files and plot in time series
#'
#' @param indir Input director-y/-ies containing meandat files (xxx.mean.dat.gz). Files within one directory may be split by groups/quantiles etc.
#' @param outdir Output directory.
#' @param expnames Names of experiments using for plot (must be the same as number of entry directories in indir).
#' @param timepoints Names of timepoints for plot (must be in the same order as in meandat files).
#' @param splitnames Names of split meandat files (e.g by groups or quantiles of replication time or expression levels).
#' @param title Main title of the plot
#' @param filename Name of the output plot file
#' @param frommin How many bp skeep from minimal value. Useful if necessary to exclude a region wihin nucleosome free regions (NFR), e.g. near TSS or ARS.
#' @param savemx Save data matrix as csv for plots
#'
#' @author Mark Boltengagen \email{m.boltengagen@@gmail.com}
#'
#' @export

avg.timeseries.plot <- function(indir,
                                expnames,
                                timepoints,
                                splitnames,
							                  outdir,
							                  title = "Average signal in time series plot",
        						            xlab = "Time",
							                  ylab = "Signal intensity",
							                  filename = "avg_timeseries_plot",
							                  frommin = 0,
							                  savemx = FALSE)
{
  if (is.null(indir))
    stop("Please specify \"indir\" input directory containing meandat files matrices (xxx.mean.dat.gz).")
  if (is.null(expnames))
    stop("Please specify \"expnames\" names of experiments.")
  if (is.null(timepoints))
    stop("Please specify \"timepoints\" names of time series.")
  if (is.null(splitnames))
    stop("Please specify \"splitnames\" names splits. Et least one must be given.")
  if (is.null(outdir))
    stop("Please specify \"outdir\" output directory.")
  
  # Create list of matrices
  mxl <- MNuc::mxl(indir=indir, expnames=expnames, timepoints=timepoints, splitnames=splitnames)
    
  # Fill mxl matrices with mean data
  for(d in 1:length(indir)) {
    # TODO make function to open mean.dat.gz files in a directory. Use fetch in data.table?
    files <- list.files(path=indir[d], pattern="*.mean.dat.gz", full.names=T, recursive=FALSE)
    data <- lapply(1:length(files), function(f) {
      # first column turns to rownames!
      df <- read.table(files[f], row.names = 1, header=TRUE, sep="\t", na.strings="NA", check.names=TRUE)
      })
        
    # detect minimal coordinate
    df <- data[[1]]
    minvalue <- min(df[1:ncol(df)])
    nr <- which(df == minvalue, arr.ind=TRUE)[1]
        
    # subset of data in columns 250 bp from minimal value and save mean in mxl
    for(nl in 1:length(data)) {
      curdf <- data[[nl]]
      for(nc in 1:ncol(curdf)) {
        curcol <- curdf[nc]
        datvec <- c(curcol[1:(nr-frommin),], curcol[(nr+frommin):nrow(curcol),])
        mxl[[nl]][d,nc] <- mean(datvec)
      }}
    }
    
    # Save mxl
    if (savemx) {
      for (m in 1:length(mxl)) {
        dir.create(file.path(outdir), recursive = TRUE, showWarnings = FALSE)
        write.table(mxl[m], file = paste(outdir, filename, "_", splitnames[m], "_matrix.csv", sep=""), quote = FALSE, sep = "\t", col.names=NA)
      }
    }
  
    # Plot mxl
    MNuc::mxl.dotplot(mxl, outdir=outdir, filename=filename, title=title, xlab=xlab, ylab=ylab)
}
