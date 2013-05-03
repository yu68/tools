#/usr/bin/env python
import sys
from bedtools import IntervalFile,Interval
from getopt import getopt

def show_help():
  print >>sys.stderr,""
  print >>sys.stderr,"Program: count_hits.py"
  print >>sys.stderr,"Usage: count_hits.py [-h] [-i interval_file] [-f data_folder] [-p ovlp_pct]"
  print >>sys.stderr,"                     [-s suffix] mark1 mark2 ... "
  print >>sys.stderr,"Example:counts.py -i Pn_E14_mm9_nucleosome_peak.bed -f ~/ChIPseq_map_data/ "
  print >>sys.stderr,"        -s _d0_extend_sort.bed mouse_H3K4me3 mouse_H3K4me2 mouse_H3K4me1"
  print >>sys.stderr,"        > count_epi_nucleosome.txt "
  print >>sys.stderr,"Arguments:  "
  print >>sys.stderr,"  -h, --help           Show this help."
  print >>sys.stderr,"  -i, --interval       file contains intervals to be counted"
  print >>sys.stderr,"  -f, --folder         folder for bed data"
  print >>sys.stderr,"  -p, --ovlp_pct       minimum overlap_percentage to be included as a count"
  print >>sys.stderr,"  -s, --suffix         uniform suffix for bed data files within the folder"
  print >>sys.stderr,"Library dependency: bedtools, getopt"
  print >>sys.stderr,""

Version="1.0"
if len(sys.argv)<2:
  show_help()
  exit(0)

opts,restlist=getopt(sys.argv[1:],"hi:f:p:s:",["help","interval","folder","ovlp_pct","suffix"])
folder="./"
percent=0.5
for o,a in opts:
  if o in ("-h","--help"):
    show_help()
    exit(0)
  elif o in ("-i","--interval"):
    interval_file=a
  elif o in ("-f","--folder"):
    folder=a
  elif o in ("-p","--ovlp_pct"):
    percent=float(a)
  elif o in ("-s","--suffix"):
    suffix=a

bedlist=[]
name='chr\tstart\tend'
for i in restlist:
  i=i.strip()
  bedlist.append(IntervalFile(folder+i+suffix))
  name=name+'\t'+i.split("_")[-1]
print name
intervals=IntervalFile(interval_file)
for i in intervals:
  line=i.chrom+'\t'+str(i.start)+'\t'+str(i.end)
  for j in bedlist:
    num=j.count_hits(i,ovlp_pct=percent)
    line=line+'\t'+str(num)
  print line
