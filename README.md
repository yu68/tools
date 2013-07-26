tools
=====

### miscellaneous tools for bioinformatics ###

---

#### count\_bam.py

given a list of genomic regions (-i) and a list of BAM files (-b), output the count of reads in each BAM file within these genomic regions.

```
usage: count_bam.py [-h] [-i INTERVAL] [-b BAMS [BAMS ...]] [-o OUTPUT]
                    [-l LEN] [-n NAME [NAME ...]]

count of reads in a list of invervals

optional arguments:
  -h, --help            show this help message and exit
  -i INTERVAL, --interval INTERVAL
                        the bed file contains location information of
                        intervals
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        the list of bam files containing mapped reads for
                        MNase-seq
  -o OUTPUT, --output OUTPUT
                        name prefix of output files: *_count.txt
  -l LEN, --length LEN  choose the center region of this length for each
                        interval to count
  -n NAME [NAME ...], --name NAME [NAME ...]
                        name of each bam sample (to be wrote on the header)

Library dependency: pysam, bedtools,numpy,math
```

#### conunt\_hits.py
given a list of genomic regions (-i) and a list of BED files with aligned reads (-f & -s), output the count of reads in each BED file within these genomic regions.

```
Usage: count_hits.py [-h] [-i interval_file] [-f data_folder] [-p ovlp_pct]
                     [-s suffix] mark1 mark2 ...
Example:counts.py -i Pn_E14_mm9_nucleosome_peak.bed -f ~/ChIPseq_map_data/
        -s _d0_extend_sort.bed mouse_H3K4me3 mouse_H3K4me2 mouse_H3K4me1
        > count_epi_nucleosome.txt
Arguments:
  -h, --help           Show this help.
  -i, --interval       file contains intervals to be counted
  -f, --folder         folder for bed data
  -p, --ovlp_pct       minimum overlap_percentage to be included as a count
  -s, --suffix         uniform suffix for bed data files within the folder
Library dependency: bedtools, getopt

```


#### sortbedTOwig.py
convert sorted BED file into WIG file ,which can be potentially single base resolution depending on read depth.  
BED file can be sorted using linux commend:
sort -k1,1 -k2,2n <unsorted.bed> > <sorted.bed>

```
Version:1.0

Library dependency: csv


Usage: sortbed2wig.py <options> -i input.bed -n name_of_output -e -l extended_read_length -s column_num_for_strand
       sortbed2wig.py <options> -i input.bed -n name_of_output -e
       sortbed2wig.py <options> -i input.bed -n name_of_output
Example:  sortbed2wig.py <options> -i mm9_H3K9me3.bed -n mm9_H3K9me3 -e -l 150 -s 4
Options:
   -h,--help          show help information
   -i,--inputfile     input bed file (with strand information for extend option)
   -o,--outputFolder  folder for output wid file (default: /home/GenomeBrowser/lab_tracks/
   -n,--wigname       name of the output wig file
   -e,--extend        extend read in bed file or not (default: false)
   -l,--readlength    the extended length of each read (default: 150, effective only when extend=True)
   -s,--strandLoc     the column # for strand information in bed (default: 4, effective only when extend=True)
```

#### pairend\_fragmentLen.py
Draw distribution of fragment length from a pairend dataset (BAM file, -i)

```
pairend_fragmentLen.py: draw distribution of fragment lenghs from pair end NGS data and fit with Gaussian Kernel Density Estimation (KDE)
Version:1.0

Library dependency: matplotlib, numpy, scipy, pysam

Usage:
    python pairend_fragmentLen.py -i [NGS_pairend_mapped_bam] -x min_x,max_x -n 100000 -o [output_figure]
    python pairend_fragmentLen.py -i [NGS_pairend_mapped_bam] -o [output_figure]
Example:
    python pairend_fragmentLen.py -i H209_pairend_5mark.sort.bam -x 0,500 -n 100000 -o H209_pairend_5mark_fragmentLen.png
Options:
    -h,--help          show help information
    -i,--inputbam      input bam file(with correspnding bai file in same folder
    -x,--xlim          range for x axis: min_x, the left bound (default 0); max_x, the right bound (default 350)
    -l,--lambda        covariance_factor lambda for KDE (default 0.25)
    -n,--num           number of fragments to be processed for plotting
    -o,--output        the output figure file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz

```

#### random\_seq\_generator.py
generate random sequences from genome specified (not exceeding the chromosome size boundary). One can adjust the mean and SD for size of  random sequences.
__probability to choose each chrom based on the size distribution.__

```
usage: random_seq_generator.py [-h] [-g GENOME] [-m MEAN] [-s SD] [-n NUM]

generate random sequences with customized length and number (for random peaks
et...)
probability to choose each chrom based on the size distribution

optional arguments:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        specify genome name to get chromosome info from
                        UCSCGB, default: mm9
  -m MEAN, --mean MEAN  mean length of each random sequence,default:200
  -s SD, --sd SD        sd of random sequence lengths,default:20
  -n NUM, --num NUM     number of sequences to be randomly
                        sampled,default:10000

library dependency: cruzdb (https://github.com/brentp/cruzdb),sqlalchemy

```
