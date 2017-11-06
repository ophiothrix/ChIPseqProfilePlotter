# rm(list=ls())


plot.coverage.profile <- function(target, bams, lib.size = NULL, fragment.size = 180, genome = NULL, flank.region = 5000L, is.PE = F, groups = NULL, groups.only = F) {
    require(GenomicRanges)
    require(csaw)
    require(edgeR)
    require(RColorBrewer)
    if(is.null(groups)) {
        clrs <- brewer.pal(length(bams), "Dark2")
    } else {
        clrs <- brewer.pal(length(bams)+length(groups), "Dark2")
    }
    
    ## group plotting sanity check
    if (is.null(groups) & groups.only) {
        warning("'groups.only' parameter is set to TRUE, but 'groups' vector is not specified. Resetting 'groups.only' parameter to FALSE")
        groups.only <- F
    }
    
    print(paste0("Plotting coverage profile for ", target))
    
    ### Set target region for when target is specified as a gene symbol
    if (is.character(target)) {
        ## Getting target region from gene requires the genome to be specified
        if (is.null(genome)) {
            stop("In order to specify target region as a gene symbol, you need to specify genome, e.g. mm10, hg19, etc")
        } else {
            require(paste0("org.", toupper(substr(genome, 1, 1)), substr(genome, 2, 2), ".eg.db"), character.only = T)
            symbs <- toTable(get(paste0("org.", toupper(substr(genome, 1, 1)), substr(genome, 2, 2), ".egSYMBOL")))
            target.id <- symbs$gene_id[match(target, symbs$symbol)]
            genome <- paste0("TxDb.Mmusculus.UCSC.", genome, ".knownGene")
            require(genome, character.only = T)
            anno <- transcriptsBy(x = get(genome), by = "gene")
            target.regions <- anno[[target.id]]
            ## In case there are multiple transcripts corresponding to a single gene, collapse them into a single range
            target.region <- reduce(target.regions)
        }
    } else {
        ### Set target region for when target is specified as a GRanges object
        if(class(target) == "GRanges") {
            target.region <- target
        } else {
            stop(paste0("You can specify the target region to be plotted either as an official gene symbol or the GRanges object. The target region specified doesn't seem to match either of those formats. Specified target region is of class '", class(target), "'"))
        }
    }
    

    # Set the flanking region by dividing gene length in half and adding flanking distance
    flank.dist <- round(width(target.region)/2, 0)+flank.region
    # Set target coordinates to the centre of the target region
    start(target.region) <- end(target.region) <- round((start(target.region) + end(target.region))/2, 0)
    # Make a matrix to hold the coverage at each position
    target.profile <- matrix(NA, length(bams), flank.dist*2+1)
    dim(target.profile)
    rownames(target.profile) <- bams
    colnames(target.profile) <- (start(target.region)-flank.dist):(end(target.region)+flank.dist)
    

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
        target.profile[i,] <- profileSites(bam.files = i, regions = target.region, range = flank.dist, ext = fragment.size, param = dedup.on)
    }

    ## If the library size not specified - plot the coverage directly
    if (is.null(lib.size)) {
        target.profile.norm <- target.profile
        y.label <- "Raw coverage"
    } else {
        # If library size is specified, normalise the coverage by given library size
        target.profile.norm <- target.profile/lib.size*10^6
        y.label <- "Normalised coverage, CPM"
    }
    
    ### If group vector is specified, take average of all samples in each group
    if (!is.null(groups)) {
        groups <- as.factor(groups)
        ## Create a matrix to hold group averages
        group.avg <- matrix(NA, nrow = length(levels(groups)), ncol = ncol(target.profile.norm))
        colnames(group.avg) <- colnames(target.profile.norm)
        rownames(group.avg) <- levels(groups)
        ## For each group compute average coverage between samples
        for (grp in levels(groups)) {
            group.avg[grp,] <- base::colMeans(target.profile.norm[groups == grp,], na.rm = T)
        }
    }

    ## Initialise the plot
    plot(colnames(target.profile.norm), target.profile.norm[1,], type="n", main=target, ylab="Normalised coverage", xlab=paste0("Position on ", seqnames(target.region)), ylim=c(0, max(target.profile.norm)))
    
    
    ## Plot individual profiles
    ## if groups.only is set to T, do not plot individual profiles
    if (!groups.only) {
        for (bam in 1:nrow(target.profile.norm)) {
            points(colnames(target.profile.norm), target.profile.norm[bam,], pch=16, lwd=0.5, type="l", col=clrs[bam])
        }
    }
    
    
    ## Plot group averages, if groups vector is supplied
    if (!is.null(groups)) {
        for (grp in 1:nrow(group.avg)) {
            points(colnames(group.avg), group.avg[grp,], pch=16, cex=0.5, type="l", col=clrs[grp+length(bams)], lty = 3, lwd = 2)
        }
    }
    
    ### Plot legends and annotations
    ## For plotting genomic region
    if (class(target) == "GRanges") {
        abline(v = c(start(target), end(target)), lty = 2, lwd  = 2)
        if (groups.only) {
            legend("topright", legend = "Target region", col = "black", lwd = 2, lty = 2, bty = "n")
        } else {
            legend("topright", legend = c(basename(bams), "Target region"), col = c(clrs[1:length(bams)], "black"), lwd = 2, lty = c(rep(1, length(bams)), 2), bty = "n")
        }
    }
    
    ## For plotting a gene
    if (is.character(target)) {
        if (unique(as.character(strand(target.regions))) == "+") {
            tss.coord <- unique(start(target.regions))
            tes.coord <- unique(end(target.regions))
        } else {
            tes.coord <- unique(start(target.regions))
            tss.coord <- unique(end(target.regions))
        }
        abline(v = tss.coord, lty = 2, lwd  = 2, col = "forestgreen")
        abline(v = tes.coord, lty = 2, lwd  = 2, col = "dodgerblue")
        if (groups.only) {
            legend("topright", legend = c("TSSs", "TESs"), col = c("forestgreen", "dodgerblue"), lwd = 2, lty = 2, bty = "n")
        } else {
            legend("topright", legend = c(basename(bams), "TSSs", "TESs"), col = c(clrs[1:length(bams)], "forestgreen", "dodgerblue"), lwd = 2, lty = c(rep(1, length(bams)), 2, 2), bty = "n")
        }
    }
    
    ## For group averages
    if (!is.null(groups)) {
        legend("topleft", levels(groups), bty = "n", col = clrs[length(bams) + unique(as.numeric(groups))], lty = 3, lwd = 2)
    }
}