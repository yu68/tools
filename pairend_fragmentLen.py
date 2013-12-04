import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import scipy
import pysam
from getopt import getopt

plt.ioff()

def show_help():
  print >>sys.stderr,"\n"
  print >>sys.stderr,"pairend_fragmentLen.py: draw distribution of fragment lenghs from pair end NGS data and fit with Gaussian Kernel Density Estimation (KDE)"
  print >>sys.stderr,"Version:"+Version+"\n"
  print >>sys.stderr,"Library dependency: matplotlib, numpy, scipy, pysam\n"
  print >>sys.stderr,"Usage:"                      
  print >>sys.stderr,"    python pairend_fragmentLen.py -i [NGS_pairend_mapped_bam] -x min_x,max_x -n 100000 -o [output_figure]"
  print >>sys.stderr,"    python pairend_fragmentLen.py -i [NGS_pairend_mapped_bam] -o [output_figure]"
  print >>sys.stderr,"Example:"
  print >>sys.stderr,"    python pairend_fragmentLen.py -i H209_pairend_5mark.sort.bam -x 0,500 -n 100000 -o H209_pairend_5mark_fragmentLen.png"
  print >>sys.stderr,"Options:"
  print >>sys.stderr,"    -h,--help          show help information"
  print >>sys.stderr,"    -i,--inputbam      input bam file(with correspnding bai file in same folder"
  print >>sys.stderr,"    -x,--xlim          range for x axis: min_x, the left bound (default 0); max_x, the right bound (default 350)"
  print >>sys.stderr,"    -l,--lambda        covariance_factor lambda for KDE (default 0.25)"
  print >>sys.stderr,"    -n,--num           number of fragments to be processed for plotting"
  print >>sys.stderr,"    -o,--output        the output figure file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz"
  print 

Version="1.0"
if len(sys.argv)<2:
    show_help()
    exit(0)
opts,restlist=getopt(sys.argv[1:],"hi:x:n:l:o:",["help","inputbam=","xlim=","num","lambda=","output="])
min_x=0
max_x=350
lambd=.25
total_num=99999999999
for o,a in opts:
  if o in ("-h","--help"):
    show_help()
    exit(0)
  elif o in ("-i","--inputbam"):
    inputbam=a
  elif o in ("-x","--xlim"):
    min_x=int(a.split(",")[0])
    max_x=int(a.split(",")[1])
  elif o in ("-l","--lambda"):
    lambd=float(a)
  elif o in ("-n","--num"):
    total_num=int(a)
  elif o in ("-o","--output"):
    output=a

ChIPseq_bam=pysam.Samfile(inputbam,'rb')
fragmentLen=[]
num=0
for i in ChIPseq_bam:
  if i.is_proper_pair and i.is_read1 and not i.is_duplicate:
    pointer=ChIPseq_bam.tell()
    try:
      j=ChIPseq_bam.mate(i)
      span=max(i.aend+j.alen-j.aend,j.aend+i.alen-i.aend)
      fragmentLen.append(span)
      num=num+1
      if (num-1)%10000==0:
        print >>sys.stderr,"Processing:%dth~%dth read pairs\r"%(num,min(num+10000-1,total_num)),
      if num==total_num:
        break
    except ValueError:
      continue
    finally:
      ChIPseq_bam.seek(pointer)

print "start drawing..."
plt.hist(fragmentLen,bins=100,normed=1,range=(min_x,max_x))
density=gaussian_kde(scipy.array(fragmentLen,float))
xs = np.linspace(min_x,max_x,200)
density.covariance_factor = lambda : lambd
density._compute_covariance()
print "KDE density computated"
plt.plot(xs,density(xs),label='KDE',color='r')
plt.xlabel("Length")
plt.ylabel("Density")
plt.title("Histogram of "+inputbam.split(".")[0])
plt.legend()
plt.savefig(output)
