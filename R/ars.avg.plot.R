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

ars.avg.plot <- function(indir = c('/home/suvar/Projects/007_ChIPseq_HU-SR/output/ars_dat/quantiles/Et',
                                       '/home/suvar/Projects/007_ChIPseq_HU-SR/output/ars_dat/quantiles/IAA'),
                              expnames = c('HU-Et','HU-IAA'),
                              timepoints = c('2h', '1.5h', '0', '10', '20', '30', '40', '50', '60', '70', '80', '90'),
                              splitnames = c('Q1', 'Q2', 'Q3', 'Q4'),
							  outdir = '/home/suvar/Desktop/',
							  setmink = TRUE,
							  title = "All features averages plot",
							  xlab = "Time",
							  ylab = "Signal intensity",
							  filename = "ars_avg_plot")
{
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
        
        if(setmink) {
            data <- lapply(data, FUN=MNuc::setmink2.dat)
        }
        
        # detect minimal coordinate
        df <- data[[1]]
        minvalue <- min(df[1:ncol(df)])
        nr <- which(df == minvalue, arr.ind=TRUE)[1]
        
        # subsett of data in columns 250 bp from minimal value and save mean in mxl
        for(nl in 1:length(data)) {
            curdf <- data[[nl]]
            for(nc in 1:ncol(curdf)) {
                curcol <- curdf[nc]
                datvec <- c(curcol[1:(nr-250),], curcol[(nr+250):nrow(curcol),])
                mxl[[nl]][d,nc] <- mean(datvec)
            }}
    }
    
    # Plot mxl
    MNuc::mxl.dotplot(mxl, outdir=outdir, filename=filename, title=title)
}

