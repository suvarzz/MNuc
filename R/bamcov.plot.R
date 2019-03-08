#' Bam Coverage Statter Plot
#'
#' @param indir Directory containing bam coverage binned bedgraph files
#' @param outdir Output directory
#' 
#' @export

bamcov.plot <- function(indir,
						outdir)

{

files <- list.files(path=inpdir, pattern="_IP_bins.bg", full.names=TRUE, recursive=FALSE)
dir.create(file.path(outdir), recursive = TRUE, showWarnings = FALSE)
png(paste(outdir, "BamCoverage_plot.png", sep=""), width=135, height=90, units="mm", res=600, pointsize=7)
# Figure parameters
par(mfrow=c(2,3)) # how many diagram on one plot?
par(mar=c(4,4,2,1.5)) # margins size

xlim <- c(0, 160)
ylim <- c(0, 160)

for(ip_file in files) {

  # Basename of current ip_file
  bn <- tools::file_path_sans_ext(basename(ip_file))
  # Split name of ip_file
  sn <- strsplit(bn, "_")
  strain_id <- sn[[1]][2]
  #Take the in_file using pattern from ip_file name
  in_file <- list.files(path=inpdir, pattern=paste(strain_id, "IN", sep="_"), full.names=TRUE, recursive=FALSE)

  # Generate data for plot 
  x <- read.table(in_file)
  y <- read.table(ip_file)

 
  plot(x$V4, y$V4,
       mgp=c(2.5,0.5,0), # margins of lables
       cex=0.3, # dot size
       bty="n", # remove box
       main=strain_id,
       col=rgb(red = 0, green = 0, blue = 0, alpha = 0.15),
       lwd=0.05,
       pch=19, 
       xlim=xlim, 
       ylim=ylim,
       xlab = "Input coverage",
       ylab = "IP coverage",
       las=1,
       axes=FALSE)
  box(lwd=0.5)
  axis(side=1, lwd=0.5, las=1)
  axis(side = 2, lwd = 0.5, las=1)
  
  #abline(lm(y$V4 ~ x$V4), col="red", lwd=0.5)
}

dev.off()
}
