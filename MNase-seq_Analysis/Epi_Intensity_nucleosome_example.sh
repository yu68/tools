# change ChIP-seq fragment length as 300

folder=~/ChIPseq_map_BAM/

nucleosome_file=/home/yu68/MNase-seq/Danpos/whole_MNase/result/pooled/Pn_E14_whole_mm9.sort.Fnor.ajClonal.smooth.peaks.xls

python Epi_Intensity_nucleosome.py -v -l 500 -N $nucleosome_file -b ${folder}mouse_H3K4me3_d0.sort.bam ${folder}mouse_H3K4me2_d0.sort.bam ${folder}mouse_H3K4me_d0.sort.bam ${folder}mouse_H3K27me3_d0.sort.bam ${folder}mouse_H3K36me3_d0.sort.bam ${folder}mouse_H3K27Ac_d0.sort.bam ${folder}mouse_H2A.Z_d0.sort.bam ${folder}mouse_H3K9me3_d0.sort.bam -n H3K4me3 H3K4me2 H3K4me1 H3K27me3 H3K36me3 H3K27ac H2AZ H3K9me3 -w 75 125 -o Epi_Intensity_nucleosome-300bp_75-125.txt


