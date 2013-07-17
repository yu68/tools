from time import time

import sys,argparse,os
from xplib import TableIO
from xplib import DBI
from xplib.Annotation import Bed
import numpy as np

def ParseArg():
    p=argparse.ArgumentParser(description="get intensity of epi-modification on single nucleosome",epilog="library dependency: xplib")
    p.add_argument("-N","--Nucleosome",dest="nucleosome",type=str,help="xls file containing nucleosome location information from Danpos output")
    p.add_argument("-b","--bams",nargs='+',dest="bams",type=str,help="bed/bam files for epigenetic data")
    p.add_argument('-f','--fmt',type=str,default='bam',help='format: bed/bam,default:bam')
    p.add_argument("-l","--length",dest="len",type=int,default=200,help="average length of ChIP-seq fragment,default:200")
    p.add_argument("-n","--name",nargs='+',dest='name',type=str,help='name of each bed sample (to be wrote on the header)')
    p.add_argument('-r','--rangeS',type=int,default=100,help='search range to find the maximum Epi-intensity in each nucleosome location (default: 100bp)')
    p.add_argument("-w","--weightP",nargs='+',dest='weightP',type=int,default=[75,125],help='parameters to calculate the weight for each read,[half_len of core nucleosome region and half_len of whole regions],default: [75,125]')
    p.add_argument("-o","--output",dest="output",type=str,help="output file name (can be .txt)")
    p.add_argument("-v",'--verbose',action='store_true',help='set to output shifted nucleosome centers for each histone mark, one bed file per mark')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def find_center(read,fmt):

    # read.id in strand for bed with for columns
    if fmt=='bed':
        strand=read.id
    elif fmt=='bam':
        strand=read.strand
    if strand=="+":
        center=read.start+half_len
    elif strand=="-":
        center=read.stop-half_len
    else:
        print >>sys.stderr,read.id,"No strand"
    return center


def getCount(read_centers,center):
    # get the count of each epi-mark in each nucleosome by searching from [nuc_center-rangeS, nuc_center+rangeS]
    # stepsize: 5bp
    max_count=0
    offset=-rangeS
    for off in range(-rangeS,rangeS,5):
        count=0
        for j_center in read_centers:
            weight = max(min(1,(ma-abs(j_center-(center+off)))*1.0/(ma-mi)),0)
            count+=weight
        if (count>max_count) or ((count>=max_count) and (abs(off)<abs(offset))): # choose the offset with largest count and smallest offset
            max_count=count
            offset=off
    return [offset,max_count]



'''
      mi
    -----
    |    \
    |     \
    |      \
    |       \
    --------
  Mid  ma

'''

global mi,ma,args,rangeS,half_len
args=ParseArg()
mi=args.weightP[0]
ma=args.weightP[1]
rangeS=args.rangeS
half_len=int(args.len/2)



def Main():
    args=ParseArg()

    #store bed files with indexing and count information:
    bam={}

    print >>sys.stderr,"Starting index bed files:"
    for i in range(len(args.bams)):
        temp_name=args.name[i]
        print >>sys.stderr,"  #Indexing for bed file of",temp_name,"\r",
        bam[temp_name]=DBI.init(args.bams[i],args.fmt)

    print >>sys.stderr
    print >>sys.stderr,"Reading nucleosome peak xls file from Danpos."
    nucleosomes=TableIO.parse(args.nucleosome,'metabed',header=True)

    print >>sys.stderr,"Initial output files..."

    out=open(args.output,"w")
    # -- for verbose ---
    if args.verbose:
        out_mark=[]
        for n in args.name:
            out_mark.append(open(n+'_shift_nucleosomes.bed','w'))
    # ------------------ 
    line_head=open(args.nucleosome,'r').readline().strip()
    line_head=line_head+"\t"+"\t".join(str(f) for f in args.name)+'\t'+"\t".join(str(f)+'_off' for f in args.name)
    print >>out,line_head
    
    print >>sys.stderr,"Start Counting..."
    num=0
    t0 = time()
    for i in nucleosomes:
        chrom=i.chr
        if i.smt_pval>0.01 or i.fuzziness_pval>0.01: continue # only choose nucleosomes with high value and low fuzziness   
        if chrom == 'chrY' or chrom == 'chrX' or chrom == 'chrM':
            continue
        num=num+1
        center=int(i.start+i.end)/2
        count=np.zeros(len(args.bams),dtype="float")
        offset=np.zeros(len(args.bams),dtype='int')
        line=str(i)
        for k,name in enumerate(bam.keys()):
            if args.fmt=='bam':
                query=bam[name].query(Bed([chrom,center-ma-(half_len-75)-rangeS,center+ma+(half_len-75)+rangeS]),method='fetch')
            else:
                query=bam[name].query(Bed([chrom,center-ma-(half_len-75)-rangeS,center+ma+(half_len-75)+rangeS]))
            read_centers=[]
            for j in query:
                read_centers.append(find_center(j,args.fmt))
            [o,c]=getCount(read_centers,center)
            count[k]=c
            offset[k]=o
            # -- for verbose ---
            if args.verbose:
                print >>out_mark[k],chrom+'\t%d\t%d'%(i.start+o,i.end+o)
            # ------------------
        line = line + "\t" + "\t".join(str(f) for f in count) + '\t' + "\t".join(str(f) for f in offset)
        if num%20000==0:
            t1 = time()
            print >>sys.stderr,"processing %dth nucleosome..., time: %.2fs."%(num,t1-t0),'\r',
            t0 = time()    
        print >>out,line
    print
    out.close()
    
    # -- for verbose ---
    if args.verbose:
        for k in out_mark:
            k.close()
    # ------------------
if __name__=="__main__":
    Main() 


        
                        

'''
TableIO.parse("Pn_E14_merge_mm9.sorted.Fnor.ajClonal.smooth.peaks.xls",'metabed',header=True)
'''



