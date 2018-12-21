                              indirs=c("/home/suvar/Projects/007_ChIPseq_HU-SR/output/bdg/Et/",
                                       "/home/suvar/Projects/007_ChIPseq_HU-SR/output/bdg/IAA/",
                                       "/home/suvar/Projects/001_ChIPseq_HU/output/bdg/",
                                       "/home/suvar/Projects/008_ChIPseq_NOC-SR/output/bdg/Et/",
                                       "/home/suvar/Projects/008_ChIPseq_NOC-SR/output/bdg/IAA/",
                                       "/home/suvar/Projects/002_ChIPseq_NOC/output/bdg/")
                              normk="/home/suvar/Projects/Normalization/Norm_k_spiking_ok.csv"
                              samples=""
                              exclude_seq = "chrM"

    # build normalizing matrix: k = (IPsc/IPsp)*(INsp/INsc), where IPsc, IPsp, INsp, INsc - reads in IP and IN in spiking (s. pombe) and sample (s. cerevisiae)
    norm.mx <- as.matrix(read.table(normk, sep="\t", row.names = 1, header=TRUE, quote =""))
    
    # make matrix filled with averages
    mx <- norm.mx
    mx[1:6,1:12] <- NA
    
    for (i in 1:nrow(mx)) {
        indir <- indirs[i]
    
    # get chromosome sizes
    seql <- MNuc::seqlevels.from.chrsizes()

    # GR from the whole genome
    sc <- BSgenome.Scerevisiae.UCSC.sacCer3
    whole_genome_gr <- GRanges(seqnames(sc), IRanges(start=1, width=seqlengths(sc)), strand="*")
    whole_genome_gr <- dropSeqlevels(whole_genome_gr, value=exclude_seq, pruning.mode="coarse")
    
    # Input bedgraph files
    files <- list.files(path=indir, pattern="\\.bdg.gz$", full.names=TRUE, recursive=FALSE)
    
    scores <- mclapply(files, function(f) {
        # Import bedGraph to GRanges
        signal <- import(f, format="bedGraph")
        
        # Sort GRange
        signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
        
        # Set seqlevels to make a correct coverage across all chromosome lengths
        seqlengths(signal) <- seql
        signal <- dropSeqlevels(signal, value=exclude_seq, pruning.mode="coarse")
        
        # Find coverage of GRange
        score <- GenomicRanges::coverage(signal, weight="score")
    })
        
    vec_data <- sapply(scores, function(score) {
        data <- GenomicRanges::binnedAverage(whole_genome_gr, score, "avg_score")
        m <- mean(data$avg_score)
        })
    
    if (length(vec_data)==10) {
        mx[i,3:12] <- vec_data
    } else { mx[i,] <- vec_data }
        }
        
        mx
        
        nmx <- mx*norm.mx
        
        # normalize amount of spiking in NR and released samples
        # 0.6667 less was added to released samples
        nmx[3,] <- nmx[3,]*0.666667
        nmx[6,] <- nmx[6,]*0.666667
        
        par(mfrow=c(2, 2))
        
        # HU not normalized
        plot(mx[1,], type="o", col="black", xlim=c(1,12), ylim=c(5.90,6.0), main="HU (not normalized)",xaxt="n")
        lines(mx[2,], type="o", col="red")
        lines(mx[3,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60", "70", "80", "90"))
        
        # NOC not normalized
        plot(mx[4,], type="o", col="black", xlim=c(1,12), ylim=c(5.90,6.0), main="NOC (not normalized)",xaxt="n")
        lines(mx[5,], type="o", col="red")
        lines(mx[6,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60", "70", "80", "90"))
        
        # HU normalized - spiking
        plot(nmx[1,], type="o", col="black", xlim=c(1,12), ylim=c(10,30), main="HU (not normalized)",xaxt="n")
        lines(nmx[2,], type="o", col="red")
        lines(nmx[3,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60", "70", "80", "90"))
        
        # NOC normalized - spiking
        plot(nmx[4,], type="o", col="black", xlim=c(1,12), ylim=c(10,30), main="NOC (not normalized)",xaxt="n")
        lines(nmx[5,], type="o", col="red")
        lines(nmx[6,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9,10,11,12), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60", "70", "80", "90"))
        
        
        plot.new()