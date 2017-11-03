target <- GRanges("chr19", IRanges(6502001, 6502250))
bams <- list.files("./data", "bam$", full.names = T)
windowed.counts <- readRDS("../Epigenetics.of.EE/unified.counting.script/atac.seq.cortex/cache/atac.seq.cortex.windowed.250bp.counts.rds")
lib.size <- windowed.counts$totals

source("./lib/plot.coverage.profile.R")
plot.coverage.profile(target, bams, lib.size = lib.size, is.PE = T)
