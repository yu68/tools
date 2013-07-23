#-------------------------------------------------------------------------------
# Name:        count reads
# Purpose:
#
# Author:      Pengfei
#
# Created:     13/09/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

#!/usr/bin/python2.7
import sys,os,argparse
from bedtools import IntervalFile
import pysam
import numpy as np
import math

def ParseArg():
    p=argparse.ArgumentParser( description = "count of reads in a list of invervals", epilog="Library dependency: pysam, bedtools,numpy,math")
    p.add_argument("-i","--interval",dest="interval",type=str,help="the bed file contains location information of intervals")
    p.add_argument("-b","--bams",nargs='+',dest="bams",type=str,help="the list of bam files containing mapped reads for MNase-seq")
    p.add_argument("-o","--output",dest="output",type=str,help="name prefix of output files: *_count.txt")
    p.add_argument("-l","--length",dest="len",type=int,default=0,help='choose the center region of this length for each interval to count, default: whole interval')
    p.add_argument("-f","--fragmentL",type=int,default=200,help="fragment length of the samples which the bam files belong to, default: 200bp")
    p.add_argument("-n","--name",nargs='+',dest='name',type=str,help='name of each bam sample (to be wrote on the header)')
    p.add_argument("-N","--normalize",action='store_true',help="specify to normalize to RPKM")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()




def Count_num(bam,region,length,nu_l,total_n):
    count=0
    if length!=0:
        start=(region.start+region.end-length)/2
        end=(region.start+region.end+length)/2
    else:
        start=region.start
        end=region.end
    for read in bam.fetch(region.chrom,start,end):
        if (not read.is_unmapped):
            if (not read.is_reverse):
                center=read.pos+nu_l/2
            elif (read.is_reverse):
                center=read.aend-nu_l/2
            if (start)<center<(end):
                count+=1
    # if normalize option
    if total_n>0:
        count=count*1000000*1000/total_n/(end-start+1)
    return count


def Main():
    args=ParseArg()


    #store bam files and count information:
    bams={}
    total_reads=np.zeros(len(args.bams))
    for i in range(len(args.bams)):
        temp_name=args.name[i]
        print >> sys.stderr, "\nReading bam file:"+temp_name+"..."
        bams[temp_name]=pysam.Samfile(args.bams[i],'rb')
        if args.normalize:
            for b in bams[temp_name]:
                if not b.is_unmapped:
                    total_reads[i]+=1
                if total_reads[i]%10000==0:
                    print >> sys.stderr, "  reading %d reads..\r"%(total_reads[i]),



    output=open(args.output+"_count.txt",'w')
    #read interval regions:
    intervals=IntervalFile(args.interval)
    header='\t'.join (str(f) for f in ['chr','start','end','name','score']) + '\t' + '\t'.join(str(f) for f in args.name )

    output.write(header+'\n')

    print >> sys.stderr,"\n\n Start counting reads for intervals..."
    for interval in intervals:
        if 'random' in interval.chrom: continue
        print_line='\t'.join (str(f) for f in [interval.chrom,interval.start,interval.end,interval.name,interval.score])
        for i in range(len(args.bams)):
            name=args.name[i]
            count=Count_num(bams[name],interval,args.len,args.fragmentL,total_reads[i])
            print_line=print_line+'\t'+str(count)
        output.write(print_line+'\n')

    #close files



    output.close()


if __name__ == '__main__':
    Main()

