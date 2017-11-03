# The idea is to take the gene coordinates 
#   > set both start and end of the gene to the position in the middle of the gene
#   > split the gene length in half and add 10kb - this is the range to supply to the profileSites function
#   > This way, the profileSites will calculate the profile over the region that is half length of the gene + 10 kb either side of the centre of the gene. This should cover all the gene + 10 kb up and down-stream
# The genes coordinates are based on the subread in-built RefSeq annotation and are gene-based, i.e. they span the length of the longest possible transcript for any given gene.

## In its simplest form, the script takes as arguments the list of bam files and the gene symbol + genome or genomic region in a form of GRanges

# rm(list=ls())


plot.coverage.profile <- function(target, bams, lib.size = NULL, fragment.size = 180, genome = NULL, expt.name="", flank.region = 5000L, is.PE = F) {
    require(GenomicRanges)
    require(csaw)
    require(edgeR)
    require(RColorBrewer)
    clrs <- brewer.pal(length(bams), "Dark2")

    print(target)
    #   Make a folder to place the graphs in
    if (expt.name != "") {
        system(paste0("mkdir ./graphs/", expt.name))
    }

    # Select the coordinates for the gene being plotted. If there's more than 1 transcript after reducing, use the first one.
    # target.region <- anno[ anno$Symbols %in% target][1]
    # target.region
    if(class(target) == "GRanges") {
        target.region <- target
    }

    # Set the flanking region by dividing gene length in half and adding flanking distance
    flank.dist <- round(width(target)/2, 0)+flank.region
    # Set target coordinates to the centre of the target region
    start(target.region) <- end(target.region) <- round((start(target.region) + end(target.region))/2, 0)
    # Make a matrix to hold the coverage at each position
    target.profile <- matrix(NA, length(bams), flank.dist*2+1)
    dim(target.profile)
    rownames(target.profile) <- bams
    colnames(target.profile) <- -flank.dist:flank.dist
    target.profile[,1:10]


    ### Extract coverage for the region from each of the bam files
    ## Set extraction parameters
    if (is.PE) {
        dedup.on <- readParam(dedup=TRUE, minq=50, pe = "both", max.frag = 1000)
    } else {
        dedup.on <- readParam(dedup=TRUE, minq=50, pe = "none")
    }
    for (i in bams){
        print(paste0("Processing ", i))
        target.profile[i,] <- profileSites(bam.files = i, regions = target, range = flank.dist, ext = fragment.size, param = dedup.on)
    }

    ## If the library size not specified - plot the coverage directly
    if (is.null(lib.size)) {
        target.profile.norm <- target.profile
    } else {
        # Normalise the coverage by given library size
        target.profile.norm <- target.profile/lib.size*10^6
    }

    ## Plot individual profiles
    plot(colnames(target.profile.norm), target.profile.norm[1,], pch=16, cex=0.5, type="n", main=target, ylab="Normalised coverage", xlab="Position relative to TSS", ylim=c(0, max(target.profile.norm)))
    for (i in 1:nrow(target.profile.norm)) {
        points(colnames(target.profile.norm), target.profile.norm[i,], pch=16, cex=0.5, type="l", col=clrs[i])
    }
    legend("topright", legend = basename(bams), col = clrs, lwd = 2, bty = "n")
    if (class(target) == "GRanges") {
        abline(v = c(-width(target)/2, width(target)/2), lty = 2, lwd  = 2)
    }


    # # Plot average profiles
    # #   h3 <- target.profile.norm[1, ]
    # t2 <- colMeans(target.profile.norm[2:4, ])
    # d5 <- colMeans(target.profile.norm[5:7, ])
    # colnames(target.profile.norm) <- as.numeric(colnames(target.profile)) + flank.dist-flank.region
    # pdf(paste("./graphs/", expt.name, "/", target$Symbols, ".pdf", sep=""), 6, 4)
    # plot(colnames(target.profile.norm), h3, pch=16, cex=0.5, type="n", main=target$Symbols, ylab="Normalised coverage", xlab="Position relative to TSS", col="forestgreen", ylim=c(0, max(target.profile.norm)))
    # points(colnames(target.profile.norm), t2, type ="l", col="dodgerblue")
    # points(colnames(target.profile.norm), d5, type ="l", col="maroon")
    # abline(v=c(0, (flank.dist-flank.region)*2), lty=2, lwd=2)
    # legend("topright", legend = c("H3 control", "WT H3K36me2", "d5KI H3K36me"), fill=c("forestgreen", "dodgerblue", "maroon"), bty="n")
    # dev.off()
}
# 
# plot.gene.profile.wo.control <- function(target, bams, eff.lib.size, fragment.size, anno, expt.name="", flank.region = 15000L, bin.coords=NULL) {
#     require(csaw)
#     print(target)
#     #   Make a folder to place the graphs in
#     if (expt.name != "") {
#         system(paste0("mkdir ./graphs/", expt.name))
#     }
#     # Select the coordinates for the gene being plotted.
#     target <- anno[ anno$Symbols %in% target]
#     target
#     target.reduced <- reduce(target)[1]
#     # Set the flanking region by dividing gene length in half and adding 10kb
#     flank.dist <- round(width(target.reduced)/2, 0)+flank.region
#     # Set gene coordinates to the centre of the gene
#     start(target.reduced) <- end(target.reduced) <- round((start(target.reduced) + end(target.reduced))/2, 0)
#     # Make a matrix to hold the coverage at each position
#     target.profile <- matrix(NA, length(bams), flank.dist*2+1)
#     rownames(target.profile) <- bams
#     colnames(target.profile) <- -flank.dist:flank.dist
#     # Calculate coverage for the gene from each of the bam files
#     for (i in bams){
#         print(paste0("Processing ", i))
#         target.profile[i,] <- profileSites(i, target.reduced, range = flank.dist, ext = fragment.size, use.strand = T)
#     }
#     # Normalise the coverage by efective library size
#     target.profile.norm <- target.profile/eff.lib.size*10^6
#     clrs <- rep(c("dodgerblue", "maroon"), each=3)
# 
#     # Plot average profiles
#     wt <- colMeans(target.profile.norm[grep("wt", bams), ])
#     d5 <- colMeans(target.profile.norm[grep("d5", bams), ])
#     colnames(target.profile.norm) <- as.numeric(colnames(target.profile)) + flank.dist-flank.region
#     pdf(paste("./graphs/", expt.name, "/", unique(target$Symbols), ".pdf", sep=""), 12, 8)
#     plot(colnames(target.profile.norm), wt, type ="l", col="dodgerblue", ylim=c(0, max(c(wt, d5))), main=unique(target$Symbols), ylab="Normalised coverage, cpm", xlab="Position relative to TSS")
#     points(colnames(target.profile.norm), d5, type ="l", col="maroon")
#     # Plot TES
#     abline(v=(flank.dist-flank.region)*2, lty=2, lwd=3)
#     # Plot all possible TSSs
#     if (unique(strand(target)) == "+") {
#         abline(v=start(target)-min(start(target)), col="darkgreen", lty=2, lwd=3)
#     } else {
#         abline(v=-(end(target)-max(end(target))), col="forestgreen", lty=2, lwd=3)
#     }
#     legend("topright", legend = c("", "", "Transcription Start Sites", "Transcription End Site", "Differentially Bound Region"), col=c("#FFFFFF00", "#FFFFFF00", "darkgreen", "black", "black"), bty="n", lty=c(1, 1, 2, 2, 1), lwd=c(0, 0, 3, 3, 5))
#     legend("topright", legend = c("WT H3K36me2", "d5KI H3K36me"), fill=c("dodgerblue", "maroon"), bty="n")
# 
#     # If bin coordinates are present, plot them
#     if (ncol(bin.coords) == 3) {
#         bin.coord <- bin.coords[bin.coords$Symbols == unique(target$Symbols), 2:3]
#         if (unique(strand(target)) == "+") {
#             bin.coord <- as.numeric(bin.coord - (start(target.reduced) - flank.dist + flank.region))
#         } else {
#             bin.coord <- as.numeric(start(target.reduced) + flank.dist - flank.region - bin.coord)
#         }
#         segments(x0 = bin.coord[1], y0 = 0, x1 = bin.coord[2], y1 = 0, lwd=5, col="black")
#     }
#     dev.off()
# }
