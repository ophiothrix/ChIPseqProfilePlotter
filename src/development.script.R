target <- GRanges("chr19", IRanges(6502001, 6502250))
target <- "Nrxn2"
bams <- list.files("./data", "bam$", full.names = T)
windowed.counts <- readRDS("../Epigenetics.of.EE/unified.counting.script/atac.seq.cortex/cache/atac.seq.cortex.windowed.250bp.counts.rds")
lib.size <- windowed.counts$totals
fragment.size = 180
genome = "mm10"
expt.name=""
flank.region = 5000L
is.PE = T

source("./lib/plot.coverage.profile.R")
plot.coverage.profile(target, bams, lib.size = lib.size, is.PE = T)
plot.coverage.profile("Sfxn4", bams, genome = "mm10", lib.size = lib.size, is.PE = T)

gene.anno <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene)
gene.anno[seqnames(gene.anno) == "chr19"]
symbs[symbs$gene_id == "94281",]

tss.anno <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, 1, 1)
tss.anno <- tss.anno[seqnames(tss.anno) == "chr19"]

source("./lib/plot.tss.coverage.profile.R")
plot.tss.coverage.profile(tss.annotation = tss.anno, bams = bams, flank.region = 3000L, lib.size = lib.size)
