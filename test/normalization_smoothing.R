                              indirs=c("/home/suvar/Projects/007_ChIPseq_HU-SR/output/bdg/Et/",
                                       "/home/suvar/Projects/007_ChIPseq_HU-SR/output/bdg/IAA/",
                                       "/home/suvar/Projects/001_ChIPseq_HU/output/bdg/",
                                       "/home/suvar/Projects/008_ChIPseq_NOC-SR/output/bdg/Et/",
                                       "/home/suvar/Projects/008_ChIPseq_NOC-SR/output/bdg/IAA/",
                                       "/home/suvar/Projects/002_ChIPseq_NOC/output/bdg/")
                              normk="/home/suvar/Projects/Normalization/Norm_k_spiking_ok.csv"
                              exclude_seq = "chrM"

    # build normalizing matrix: k = (IPsc/IPsp)*(INsp/INsc), where IPsc, IPsp, INsp, INsc - reads in IP and IN in spiking (s. pombe) and sample (s. cerevisiae)
    norm.mx <- as.matrix(read.table(normk, sep="\t", row.names = 1, header=TRUE, quote =""))
    
    # omit last 3 samples
    norm.mx <- norm.mx[,1:9]
    
    # make matrix filled with averages
    mx <- norm.mx
    mx[1:6,1:9] <- NA
    
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
    
    if (length(vec_data)==7) {
        mx[i,3:9] <- vec_data
    } else { mx[i,] <- vec_data }
        }
        
        # mean all chip-seq data
        mx
        
        # normalize to spiking
        nmx <- mx*norm.mx
        
        # correction of not-released HU (less spiking added 30/45 ul) 0.6667 less was added to released samples
        nmx[3,] <- nmx[3,]*0.66666667
        
        # correction for NOC-SR-Et outlier
        nmx[4,7] <- mean(nmx[4,6], nmx[4,8])

        par(mfrow=c(1, 2), oma = c(0, 0, 2, 0))
        
        # HU normalized - spiking
        plot(nmx[1,], type="o", col="black", xlim=c(1,9), ylim=c(10,30), 
             main="HU normalization",
             xaxt="n",
             xlab="Time, min",
             ylab="H4K16Ac genomic mean")
        lines(nmx[2,], type="o", col="red")
        lines(nmx[3,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60"))
        legend("topright", c("HU-SR-Et","HU-SR-IAA", "HU-NR"), lwd=c(2,2,2), bty='n', col=c("black", "red", "green"))
        mtext("Average H4K16Ac signals per genomes (spiking normalization)", outer = TRUE, cex = 1.5)
        
        # Predict normalizing coefficients for HU-SR-Et
        x <- 1:9
        hu.et <- as.vector(nmx[1,1:9])
        hu.et.fit <- lm(hu.et ~ poly(x,2))
        pr.hu.et <- predict(hu.et.fit)
        lines(pr.hu.et, lty=2, col="black")
        
        # Predict normalizing coefficients for HU-SR-IAA
        x <- 1:9
        hu.iaa <- as.vector(nmx[2,1:9])
        hu.iaa.fit <- lm(hu.iaa ~ poly(x, 2))
        pr.hu.iaa <- predict(hu.iaa.fit)
        lines(pr.hu.iaa, lty=2, col="red")
        
        # Predict normalizing coefficients for HU-NR
        x <- 3:9
        hu.nr <- as.vector(nmx[3,3:9])
        hu.nr.fit <- lm(hu.nr ~ poly(x, 2))
        pr.hu.nr <- predict(hu.nr.fit)
        lines(pr.hu.nr~x, lty=2, col="green")
        

        # NOC normalized - spiking
        plot(nmx[4,], type="o", col="black", xlim=c(1,9), ylim=c(10,30), 
             main="NOC normalization", 
             xaxt="n",
             xlab="Time, min",
             ylab="H4K16Ac genomic mean")
        lines(nmx[5,], type="o", col="red")
        lines(nmx[6,], type="o", col="green")
        axis(1, at=c(1,2,3,4,5,6,7,8,9), labels=c("1.5h", "2h", "0", "10", "20", "30", "40", "50", "60"))
        legend("topright", c("NOC-SR-Et","NOC-SR-IAA", "NOC-NR"), lwd=c(2,2,2), bty='n', col=c("black", "red", "green"))

        # Predict normalizing coefficients for NOC-SR-Et
        x <- 1:9
        noc.et <- as.vector(nmx[4,1:9])
        noc.et.fit <- lm(noc.et ~ poly(x,2))
        pr.noc.et <- predict(noc.et.fit)
        lines(pr.noc.et, lty=2, col="black")
        
        # Predict normalizing coefficients for NOC-SR-IAA
        x <- 1:9
        noc.iaa <- as.vector(nmx[5,1:9])
        noc.iaa.fit <- lm(noc.iaa ~ poly(x, 2))
        pr.noc.iaa <- predict(noc.iaa.fit)
        lines(pr.noc.iaa, lty=2, col="red")
        
        # Predict normalizing coefficients for NOC-NR
        x <- 3:9
        noc.nr <- as.vector(nmx[6,3:9])
        noc.nr.fit <- lm(noc.nr ~ poly(x, 2))
        pr.noc.nr <- predict(noc.nr.fit)
        lines(pr.noc.nr~x, lty=2, col="green")
        
        # matrix with predicted averages
        pr.mx <- mx
        pr.mx[1:6,1:9] <- NA
        
        pr.mx[1,] <- pr.hu.et
        pr.mx[2,] <- pr.hu.iaa
        pr.mx[3,3:9] <- pr.hu.nr
        pr.mx[4,] <- pr.noc.et
        pr.mx[5,] <- pr.noc.iaa
        pr.mx[6,3:9] <- pr.noc.nr
        
        # matrix of coefficients
        # predicted normalizing coefficients
        pr.k.mx <- pr.mx/mx
        
        pr.k.mx
        
        # Normalize bedgraph files
        MNuc::bdg.norm(indir=indirs[1], 
                       outdir="/home/suvar/Projects/010_Spiking/HU_SR_Et/bdg/", 
                       normk=pr.mx[1,])
        