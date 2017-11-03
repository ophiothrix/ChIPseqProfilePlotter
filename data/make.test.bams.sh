#!/bin/bash
for fls in `ls /no_backup/GD/projects/Epigenome/TgDyrk1a_mouse/epigenome/Chromatin_interaction_prediction/our_data_mm10/ATACseq/*.dedup.sort.bam`
do
	echo $fls
	outf=$(basename $fls)
	echo $outf
	samtools view -bh -o $outf $fls chr19 &&
	samtools index $outf
done
wait
