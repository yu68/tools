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
