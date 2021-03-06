# rm(list = ls())
# gc()
require(GenomicRanges)
#### Use case 1 ####
### Plot coverage from a set of bam files over given genomic range ###
source("./lib/plot.coverage.profile.R")

### Required parameters ###
## Set path to bam files to be plotted
## Bam files need to be sorted by genomic position and indexed
path.to.bams <- list.files("./data", "bam$", full.names = T)

## Set genomic range to plot
## Test bam files are restricted to chromosome 19
target.range <- GRanges("chr19", IRanges(34921751, 34922500))

### Optional parameters ###
## Flanking region size - specifies the size of the region flanking the target region
## Defaults to 5kb
flank.region = 3000L

## Library sizes - if specified, the coverages are converted to CPM
library.sizes <- c(15489132, 14496836, 23424741, 20148416)

## Fragment size - corresponds to average library fragment size. Reads are extended to this length
frag.size <- 200

## Paired ended - used to specify if bam file contains paired reads
## If set as TRUE, the fragment.size parameter is ignored and read pairs are treated as fragments instead
## If Bam file contains paired reads, but is.PE is set to FALSE (default) the mates are (incorrectly) considered as independent reads.
is.PE = T

## genome - only required when target region is specified as gene symbol

## groups - character or factor vector of the same length as the vector of bam files, corresponding to the sample groups
## If vector of groups is specified, additional lines are plotted for the average coverage between samples of each group.
groups <- c("ctrl", "ctrl", "expt", "expt")

## groups only - if "groups" is not null, a scalar prescribing whether only group average should be plotted, or whether it should be plotted alongside with individual coverages (default).

## Calculate and plot the coverage
png("./graphs/use.case.1.png", 800, 600)
plot.coverage.profile(target = target.range, bams = path.to.bams, lib.size = library.sizes, is.PE = is.PE, flank.region = flank.region, fragment.size = frag.size, groups = groups)
dev.off()

#### Use case 2 ####
### Plot coverage from a set of bam files over a gene of interest ###
source("./lib/plot.coverage.profile.R")

### Required parameters ###
## Set path to bam files to be plotted
## Bam files need to be sorted by genomic position and indexed
path.to.bams <- list.files("./data", "bam$", full.names = T)

## Set genomic range to plot
## Test bam files are restricted to chromosome 19
target.gene <- "Sfxn4"

## Library sizes - if specified, the coverages are converted to CPM
library.sizes <- c(15489132, 14496836, 23424741, 20148416)

## Paired ended - used to specify if bam file contains paired reads
## If set as TRUE, the fragment.size parameter is ignored and read pairs are treated as fragments instead
## If Bam file contains paired reads, but is.PE is set to FALSE (default) the mates are (incorrectly) considered as independent reads.
is.PE = T

## genome - only required when target region is specified as gene symbol
genome <- "mm10"

## groups - character or factor vector of the same length as the vector of bam files, corresponding to the sample groups
## If vector of groups is specified, additional lines are plotted for the average coverage between samples of each group.
groups <- c("ctrl", "ctrl", "expt", "expt")


## Calculate and plot the coverage
png("./graphs/use.case.2.png", 800, 600)
plot.coverage.profile(target = target.gene, bams = path.to.bams, lib.size = library.sizes, is.PE = is.PE, genome = genome, groups = groups, groups.only = T)
dev.off()

#### Use case 3 ####
### Plot averaged coverage over multiple genomic regions ###
## The script is only tested with genomic regions of equal size
source("./lib/plot.tss.coverage.profile.R")

### Required parameters ###
## Set path to bam files to be plotted
## Bam files need to be sorted by genomic position and indexed
path.to.bams <- list.files("./data", "bam$", full.names = T)

## Generate the list of genomic ranges over which to calculate the coverage
## Test bam files are restricted to chromosome 19
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
tss.anno <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, 1, 1)
tss.anno <- tss.anno[seqnames(tss.anno) == "chr19"]

### Optional parameters ###
### Optional parameters ###
## Flanking region size - specifies the size of the region flanking the target region
## Defaults to 5kb
flank.region = 3000L

## Library sizes - if specified, the coverages are converted to CPM
library.sizes <- c(15489132, 14496836, 23424741, 20148416)

## Paired ended - used to specify if bam file contains paired reads
## If set as TRUE, the fragment.size parameter is ignored and read pairs are treated as fragments instead
## If Bam file contains paired reads, but is.PE is set to FALSE (default) the mates are (incorrectly) considered as independent reads.
is.PE = T

## Calculate and plot the coverage
png("./graphs/use.case.3.png", 800, 600)
plot.tss.coverage.profile(target.ranges = tss.anno, bams = path.to.bams, flank.region = flank.region, lib.size = library.sizes, is.PE = is.PE)
dev.off()
