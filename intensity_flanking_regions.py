#-------------------------------------------------------------------------------
# Name:        intensity_flaking_regions
# Purpose:
#
# Author:      Pengfei
#
# Created:     03/06/2012
# Copyright:   (c) Pengfei 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import sys,os,argparse
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from xplib import TableIO
import numpy as np
import scipy
from scipy.interpolate import spline

plt.ioff()

#for re-sampling of data
def find_peak(x,y,mi,ma):
    sx=x[(x>mi)&(x<ma)]
    sy=y[(x>mi)&(x<ma)]
    return sx[np.argmax(sy)]

def detect_peak(x,y):
    n=np.arange(0,5,dtype="float")
    l=np.arange(0,5,dtype="float")
    for i in range(0,5):
        l[i]=find_peak(x,y,185*i+120,185*i+160)
    print l
    n=np.array([n,np.ones(5)])
    slope=np.linalg.lstsq(n.T,l)[0]
    return slope[0]
    
    



def ParseArg():
    p=argparse.ArgumentParser( description = "draw intensity flanking specific regions", epilog="Library dependency: xplib,pysam, matplotlib, scipy, numpy")
    p.add_argument("input",type=str,metavar='input_rmdup',help='duplication removed input bam file for mapped sequencing data')
    p.add_argument("-r","--regions",type=str,dest="regions",help="regions whose flanking region to be drawn")
    p.add_argument("-o","--output",type=str,dest="output",help="the output figure file, can be format of emf, eps, pdf, png, ps, raw, rgba, svg, svgz")
    #p.add_argument("-l","--lambda",type=float,dest="lambd",default=0.05,help="covariance_factor lambda for KDE (default 0.05)")
    p.add_argument("-l","--len",type=int,dest="len",default=74,help="length of each read covered for each fragment,default:74")
    p.add_argument("-x","--xlim",type=int, nargs='+',default=[-1000,1000],dest="xlim",help="range for x axis: min_x, the left bound; max_x, the right bound. (default: -1000,1000)")
    p.add_argument("-n","--num",type=int,dest="num",default=5000000,help="number of distances included to draw density function (default: 5000000)")
    p.add_argument("-p","--phase",action='store_true',help="whether cauculate phase in a '*_phase_distances.txt' file")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)

    return p.parse_args()

def Main():
    args=ParseArg()
    min_x=int(args.xlim[0])
    max_x=int(args.xlim[1])
    print  >>sys.stderr,"Reading dup-removed bam file..."
    reads=pysam.Samfile(args.input,'rb')
    counts=np.zeros(max_x-min_x)

    half_l=int(args.len/2)
    num=0
    n=0

    #input bed file for regions
    regions=TableIO.parse(args.regions,"bed")

    #for resampling purpose
    slopes=[]
    count_temp=np.zeros(max_x-min_x)

    for i in regions:
        n=n+1
        #print "%d reads included"%(num),"\r",
        chrom=i.chr
        center=int((i.start+i.stop)/2)
        # Skip chrM            
        if (chrom=='chrM' or len(chrom)>6):
            continue

        for j in reads.fetch(chrom,center+min_x-1,center+max_x):
            # loc is middle points of each reads
            if (not j.is_unmapped) and (j.is_reverse):
                loc=j.aend-74
            elif (not j.is_unmapped) and (not j.is_reverse):
                loc=j.pos+74
            if (min_x+half_l<=loc-center<=max_x-half_l):
                counts[(loc-center-half_l-min_x):(loc-center-half_l+args.len-min_x)]+=1
                count_temp[(loc-center-half_l-min_x):(loc-center-half_l+args.len-min_x)]+=1
                num=num+1
                if num%10000==0:
                    print >>sys.stderr,"Processing:",num,"reads",chrom,center,"region#",n,"\r",

                #resampling purpose, find standard deviation of phase distance.
                if (args.phase) and num%50000==0:
                    xnew=np.linspace(0,max_x,200)
                    smooth=spline(range(0,max_x),count_temp[(0-min_x):],xnew)
                    slope=detect_peak(xnew,smooth)
                    print >>sys.stderr,"In process of %dth 50000 distances,slope is %f"%(num/50000,slope)
                    count_temp=np.zeros(max_x-min_x)
                    slopes.append(slope)
                    
        if num>=args.num:
            break
    reads.close()
    
    print
    print np.sum(counts)
    if np.sum(counts)==0:
        print "EXIT"
        sys.exit(0)
    print
    print >>sys.stderr,"start drawing..."

    #plt.hist(distancePool,bins=(max_x-min_x),normed=0,range=(min_x,max_x),alpha=0.4)
    #hist,bins=np.histogram(distancePool,bins=(max_x-min_x),range=(min_x,max_x),density=False)
    np.savetxt(args.output.split(".")[0]+"_%d~%dbp.txt"%(min_x,max_x),np.column_stack((np.array(range(min_x,max_x),int),counts)),delimiter='\t')
    
    plt.plot(range(min_x,max_x),counts,color='r')

    # smoth the plot
    xnew=np.linspace(min_x,max_x,400)
    smooth=spline(range(min_x,max_x),counts,xnew)
     
    #resampling purpose
    if args.phase:
        np.savetxt(args.output.split(".")[0]+"phase_distances.txt",slopes,delimiter="\t")
        slope=detect_peak(xnew,smooth)
        print >>sys.stderr,"final slope",slope
    plt.plot(xnew,smooth,color='g')

    plt.xlabel("Length")
    plt.ylabel("Density")
    plt.xlim(min_x,max_x)

    plt.title("intensity of "+os.path.basename(args.input).split(".")[0])
    plt.legend()
    plt.savefig(args.output)
    print >>sys.stderr,"output figure file generated!!"

if __name__ == '__main__':
    Main()





