### Script to plot average coverage over multiple TSSs

plot.tss.coverage.profile <- function(tss.annotation, bams, lib.size = NULL, fragment.size = 180, expt.name="", flank.region = 5000L, is.PE = F) {
    require(GenomicRanges)
    require(csaw)
    require(edgeR)
    require(RColorBrewer)
    clrs <- brewer.pal(length(bams), "Dark2")
    
    # Make a matrix to hold the coverage at each position
    target.profile <- matrix(NA, length(bams), flank.region*2+1)
    dim(target.profile)
    rownames(target.profile) <- bams
    colnames(target.profile) <- -flank.region:flank.region
    
    ### Extract coverage for the region from each of the bam files
    ## Set extraction parameters
    if (is.PE) {
        dedup.on <- readParam(dedup=TRUE, minq=50, pe = "both", max.frag = 1000)
    } else {
        dedup.on <- readParam(dedup=TRUE, minq=50, pe = "none")
    }
    ## Count reads
    for (i in bams){
        print(paste0("Processing ", i))
        target.profile[i,] <- profileSites(bam.files = i, regions = tss.annotation, range = flank.region, ext = fragment.size, param = dedup.on, average = T, normalize = "none")
    }
    
    ## If the library size not specified - plot the coverage directly
    if (is.null(lib.size)) {
        target.profile.norm <- target.profile
    } else {
        # If library size is specified, normalise the coverage by given library size
        target.profile.norm <- target.profile/lib.size*10^6
    }
    
    plot(colnames(target.profile.norm), target.profile.norm[1,], type="n", main="Average coverage profiles over TSS", ylab="Normalised coverage", xlab="Position relative to TSS", ylim=c(0, max(target.profile.norm)))
    for (bam in 1:nrow(target.profile.norm)) {
        points(colnames(target.profile.norm), target.profile.norm[bam,], pch=16, cex=0.5, type="l", col=clrs[bam])
    }
    abline(v = seq(-flank.region, flank.region, 250), lty = 2)
    legend("topright", legend = basename(bams), col = clrs, lwd = 2)

}