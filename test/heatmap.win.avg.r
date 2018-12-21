#' Heatmap of binned averages from bedgraph files
#' 
#' @details Heatmap plot of binned averages from bedgraph files
#' 
#' @param indir Directory containing bedgraph files with peaks
#' @param outdir Output directory
#' @param chromosome Chromosome to draw
#' @param fname File name
#' @param norm Vector of normalization coeficients
#' @param log Log difference (to common average)
#' 
#' @return NULL
#' @export

heatmap.win.avg <- function(indir="/home/suvar/Desktop/test_heatmap/IAA/",
                            outdir="/home/suvar/Desktop/test_heatmap/IAA/",
                            chromosome = 7,
                            # HU-Et
                            #norm = c(2.42149621199355, 2.86362316624363, 2.93820113964696, 3.08496113754195, 3.20350471237092, 3.3, 3.43696810741292, 3.425, 3.41924243158332, 4.32304105965692, 4.69223264170933, 3.59356408750295),
                            # HU-IAA
                            norm = c(2.42149621199355, 2.11643548730179, 2.91244993294935, 2.26411435822311, 2.76278946940894, 2.46749312428193, 2.76310762361151, 2.0662504878982, 3.29877868837283, 2.13639926249044, 2.75545300546037, 2.5922193133499),
                            fname = "HeatMap_win_IAA_mean",
                            log = TRUE)
{
    # Create output directory
    dir.create(file.path(outdir), recursive = TRUE)

    # GR from the whole genome
    sc <- BSgenome.Scerevisiae.UCSC.sacCer3
    whole_genome_gr <- GRanges(seqnames(sc), IRanges(start=1, width=seqlengths(sc)), strand="*")
    
    # sliding windows for the whole genome
    wins <- GenomicRanges::slidingWindows(whole_genome_gr, width=15000, step=2000)

    # Input bedgraph files
    files <- list.files(path=indir, pattern="\\.bdg.gz$", full.names=TRUE, recursive=FALSE)
    
    cat("Calculation of matrix for files started and takes time.")
    list_data <- mclapply(files, function(f) {

        # Import bedGraph to GRanges
        signal <- import(f, format="bedGraph")

        # Sort GRange
        signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
        
        # Find coverage of GRange
        score <- GenomicRanges::coverage(signal, weight="score")
    
        # Calculate binned average
        win_data <- GenomicRanges::binnedAverage(wins[[chromosome]], score, "avg_score")
        mcols(win_data)$avg_score
    
    })

    mx <- matrix(unlist(list_data), nrow = length(list_data), byrow = TRUE)
    
    # Normalizing
    mx <- mx*norm
    
    if (log) {
    # log difference with whole mean
    mx <- log2(mx/mean(mx)) # log2(mx/20.73039) 
    # log difference with the first sample
    #mx <- log2(t(t(mx)/mx[1,]))
    # remove extreme values
    mx[!is.finite(mx)] <- 0
    }

    # set limits for better color usage
    mx[mx < -1] <- -1
    mx[mx > 1] <- 1
    
    # get seqlevels to adjust figure size to chromosome size
    seql <- MNuc::seqlevels.from.chrsizes()
    width <- as.numeric(round(0.00347*seql[chromosome], digits=0))
    
    png(file=paste(outdir, fname, "_", chromosome, ".png", sep=""), width=width, height=1000)
    par(mai=c(0.1,0.1,0.1,0.1))

    # set color scheme
    my_palette <- colorRampPalette(c("blue", "light blue", "white", "orange", "brown"))(n = 299) 

    # draw heatmap
    heatmap.2(mx,
              # heatmap parameters
              col=my_palette,
              Rowv=F,
              Colv=F,
              dendrogram = "none",
              trace="none", 
              labRow=F,
              labCol=F,
              
              # remove margins from the plot 
              lwid=c(0.01,5),
              lhei=c(0.01,5),
              margins=c(0.3,0.3),
             
              # color key parameters
              density.info="none",
              key = F) # show color key?
    dev.off()
    
}
