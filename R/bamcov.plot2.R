#' Bam Coverage Statter Plot
#'
#' @param indir Directory containing bam coverage binned bedgraph files
#' @param outdir Output directory
#' 
#' @export

bamcov.plot <- function(indir1 = "/home/suvar/Projects/007_ChIPseq_HU-SR/output/bam_cov/Et/IP",
                        indir2 = "/home/suvar/Projects/007_ChIPseq_HU-SR/output/bam_cov/IAA/IP",
						outdir = "/home/suvar/Desktop/",
						lables = c("HU-SR 2h", "HU-SR 0'", "HU-SR 10'", "HU-SR 20'", "HU-SR 30'",
						           "HU-SR 40'", "HU-SR 50'", "HU-SR 60'", "HU-SR 70'", "HU-SR 80'",
						           "HU-SR 90'"),
                        filename = "BamCoverage_plot.png",
						log=FALSE)

{

files1 <- list.files(path=indir1, pattern=".bg.gz", full.names=TRUE, recursive=FALSE)
files2 <- list.files(path=indir2, pattern=".bg.gz", full.names=TRUE, recursive=FALSE)

dir.create(file.path(outdir), recursive = TRUE)

# PLOT
png(paste(outdir, filename, sep=""), width=150, height=200, units="mm", res=600, pointsize=7)
# Figure parameters
par(mfrow=c(4,3)) # how many diagram on one plot?
par(mar=c(4,4,2,1.5)) # margins size

xlim <- c(0, 170)
ylim <- c(0, 170)

silent <- lapply(seq_along(files1), function(i) {
    # Generate data for plot 
    x <- read.table(files1[i])
    y <- read.table(files2[i])
    
    if (log) {
        x$V4 <- log2(x$V4)
        y$V4 <- log2(y$V4)
        xlim <- c(0, log2(xlim[2]))
        ylim <- c(0, log2(ylim[2]))
    }
    
    plot(x$V4, y$V4,
         mgp=c(2.5,0.5,0), # margins of lables
         cex=0.3, # dot size
         bty="n", # remove box
         main=lables[i],
         col=rgb(red = 0, green = 0, blue = 0, alpha = 0.1),
         lwd=0.05,
         pch=19, 
         xlim=xlim, 
         ylim=ylim,
         xlab = "Control (Et)",
         ylab = "Test (IAA)",
         las=1,
         axes=FALSE)
    box(lwd=0.5)
    axis(side=1, lwd=0.5, las=1)
    axis(side = 2, lwd = 0.5, las=1)
    
    abline(0, 1, col="red", lwd=0.25)
})

dev.off()
}

