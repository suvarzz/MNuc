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

heatmap.win.avg <- function(indir,
                            outdir,
                            chromosome = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                            norm = NULL,
                            fname = "HeatMap",
                            log = TRUE,
                            width = NULL,
                            step = NULL,
                            comb=TRUE)
{
    # Create output directory
    dir.create(file.path(outdir), recursive = TRUE)

    # GR from the whole genome
    sc <- BSgenome.Scerevisiae.UCSC.sacCer3
    whole_genome_gr <- GRanges(seqnames(sc), IRanges(start=1, width=seqlengths(sc)), strand="*")
    
    # Input bedgraph files
    files <- list.files(path=indir, pattern="\\.bdg.gz$", full.names=TRUE, recursive=FALSE)
    
    scores <- mclapply(files, function(f) {
        # Import bedGraph to GRanges
        signal <- import(f, format="bedGraph")
        
        # Sort GRange
        signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
        
        # Find coverage of GRange
        score <- GenomicRanges::coverage(signal, weight="score")
    })
    
    
    # sliding windows for the whole genome
    wins <- GenomicRanges::slidingWindows(whole_genome_gr, width=width, step=step)
    
    # get chromosome sizes
    seql <- MNuc::seqlevels.from.chrsizes()
    
    for (chr in chromosome) {
    wins.chr <- wins[[chr]]
    
    list_data <- mclapply(scores, function(score) {
        # Calculate binned average
        win_data <- GenomicRanges::binnedAverage(wins.chr, score, "avg_score")
        mcols(win_data)$avg_score
    })

    mx <- matrix(unlist(list_data), nrow = length(list_data), byrow = TRUE)
    
    # Normalizing
    if (!is.null(norm)) {
        mx <- mx*norm
    }
    
    if (log) {
    # log difference with whole mean
   mx <- log2(mx/18.10279)#log2(mx/mean(mx)) # log2(mx/20.58096)  # 
    # log difference with the first sample
    #mx <- log2(t(t(mx)/mx[1,]))
    # remove extreme values
    mx[!is.finite(mx)] <- 0
    }

    # set limits for better color usage
    mx[mx < -1] <- -1
    mx[mx > 1] <- 1
    
    # width of plot addjusted to chromosome size
    width <- as.numeric(round(0.00347*seql[chr], digits=0))
    
    png(file=paste(outdir, fname, "_", chr, ".png", sep=""), width=width, height=1000)
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
    
    if(comb==TRUE) {
        MNuc::combine.heatmaps(indir=outdir,outdir=outdir)
    }
}
