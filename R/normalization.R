#' Calculate normalizing coefficients
#' 
#' @param fbam Directory containing filtered bam files for peak calculaton.
#' @param fbams  Directory containing filtered bam files of spiking using for normalizing.
#' @param outdir Output directory.
#' @param filename Name of output file. 
#' 
#' @return NULL
#' @export

find.norm.k <- function(fbam,
                        fbams,
                        outdir,
                        filename = "normfactors.txt")
{
    bam_files <- list.files(path=fbam, pattern="\\.bam$", full.names=TRUE, recursive=FALSE)
    bam_spikes_files <- list.files(path=fbams, pattern="\\.bam$", full.names=TRUE, recursive=FALSE)
    cat("Calculation of factors using number of reads in bam files takes time\n")
    factors <- sapply(seq_along(bam_files), function(i) {
               reads <- as.numeric(system(paste("samtools view -c", bam_files[i]), intern = TRUE))
               reads_spikes <- as.numeric(system(paste("samtools view -c", bam_spikes_files[i]), intern = TRUE))
               k <- reads/reads_spikes
    })
    dt <- data.table(name=basename(bam_files), factor=factors)
    dir.create(file.path(outdir), recursive = TRUE, showWarnings = FALSE)
    fwrite(dt, file=paste(outdir, "/", filename, sep=""), sep="\t")
}

#' Split files in a directory into groups
#' 
#' @param norm File containing normalizing coefficients (use find.norm.k function).
#' @param mdat Directory containing mean dat files which should be normalized.
#' @param outdir Output directory for normalized mean dat file.
#' 
#' @return NULL
#' @export

norm.mean.dat <- function(norm,
                          mdat,
                          outdir) 
{
    normk <- read.table(norm, header=TRUE, sep="\t", na.strings="NA", check.names=TRUE)
    files <- list.files(path=mdat, pattern="*.dat.gz", full.names=T, recursive=FALSE)
    dir.create(file.path(outdir), recursive=T, showWarnings = FALSE)
    silent <- lapply(1:length(files), function(f) {
        dat <- read.table(files[f], header=TRUE, sep="\t", na.strings="NA", check.names=TRUE)
        for (cln in 2:ncol(dat)) {
            dat[, cln] <- dat[,cln]*normk[cln-1, ]$factor
        }
        filename=paste(outdir, tools::file_path_sans_ext(basename(files[f])), sep="")
        write.table(dat, file=filename, quote=FALSE, append=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
        system2("gzip", args=filename)
    })
}

#' Normalize data in data frame to minimal value (Null value) by subtraction of difference.
#' 
#' @param df Data frame from mean data file.
#' 
#' @return data frame (mean data)
#' @export

setmin.dat <- function(df)
{
        minvalue <- min(df[2:ncol(df)])
        nr <- which(df == minvalue, arr.ind=TRUE)[1]
        for (i in 2:ncol(df)) {
            dif <- df[nr, i]-minvalue
            df[,i] <- df[,i]-dif
        }
        return(df)
}

#' Normalize data in data frame to minimal value (Null value) by coefficient.
#' 
#' @param df Data frame from mean data file.
#' 
#' @return data frame (mean data)
#' @export

setmink.dat <- function(df)
{
    minvalue <- min(df[2:ncol(df)])
    nr <- which(df == minvalue, arr.ind=TRUE)[1]
    for (i in 2:ncol(df)) {
        k <- minvalue/df[nr, i]
        df[,i] <- df[,i]*k
    }
    return(df)
}

#' Normalize data in data frame to minimal value (Null value) by coefficient.
#' For data frames where rownames contain coordinates, and data frame contains data only.
#' @param df Data frame from mean data file.
#' 
#' @return data frame (mean data)
#' @export

setmink2.dat <- function(df)
{
    minvalue <- min(df[1:ncol(df)])
    nr <- which(df == minvalue, arr.ind=TRUE)[1]
    for (i in 1:ncol(df)) {
        k <- minvalue/df[nr, i]
        df[,i] <- df[,i]*k
    }
    return(df)
}

#' Find normalizing coefficients based on Nucleosome Free Regions (NFR).
#' For data frames where rownames contain coordinates, and data frame contains data only.
#' @param dir Directory containing vector of meandata for all features containing Nucleosome Free Regions.
#' e.g. TSS, ARS
#' 
#' @param snames Sample names
#' @param dirnames Vector containing names of features (one feature in the input directory)
#' @param weights Vector containings number of every feature in a genome. Need for calculation of weighted means.  
#' @return matrix of coefficients for each sample in each given feature
#' @export

mxkmink <- function(indir=c("/home/suvar/Projects/007_ChIPseq_HU-SR/output/tss_data/Et/means/all/Meandata-All_1.mean.dat.gz",
                          "/home/suvar/Projects/007_ChIPseq_HU-SR/output/ars_dat/all/Et/HU_SR_Et_ars_meandata_All_1.mean.dat.gz"),
                    snames=c("2h", "1.5h", "0'", "10'", "20'", "30'", "40'" , "50'", "60'", "70'", "80'" , "90'"),
                    dirnames=c("et_tss", "et_ars"),
                    weights=c(4547, 352)) # tss, ars
{
    # Read files 
    data <- lapply(indir, function(f) {
        df <- read.table(f, header=TRUE, sep="\t", na.strings="NA", check.names=TRUE)
    })
    ## create matrix of normalizing coefficients
    mx = matrix(data=0, nrow=length(indir), ncol=ncol(data[[1]])-1)
    colnames(mx) <- snames
    rownames(mx) <- dirnames
    
   ## Calculate normalizing coefficients
   for (l in 1:length(data)) {
       df <- data[[l]]
       minvalue <- min(df[2:ncol(df)])
       nr <- which(df == minvalue, arr.ind=TRUE)[1]
       for (i in 2:ncol(df)) {
           k <- minvalue/df[nr, i]
           mx[l,i-1] <- k
       }
   }
    
   # Weighted mean
   mx <- rbind(mx, wmean_norm_k=apply(mx, 2, function(x) weighted.mean(x, w=weights)))
   mx
   return(mx)
}
