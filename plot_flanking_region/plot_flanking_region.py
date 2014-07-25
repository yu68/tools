import matplotlib
matplotlib.use('Agg')
import numpy as np
from xplib import TableIO
import pysam
from xplib.Annotation import Bed
import sys,os,argparse
import scipy.cluster.hierarchy as sch
from collections import Counter
import pylab
from Bio.Cluster import somcluster
from scipy.cluster.vq import kmeans,vq
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
from collections import Counter

# font size
from matplotlib.font_manager import FontProperties
fontP = FontProperties()
fontP.set_size('small')

plt.ioff()

def ParseArg():
    ''' This Function Parse the Argument '''
    p=argparse.ArgumentParser( description = 'plot heatmap or average patterns for intensities around flanking regions of interested locations', epilog='Library dependency : xplib, pylab, scipy, numpy, pysam, matplotlib')
    group=p.add_mutually_exclusive_group()
    group.add_argument("-H","--Heatmap",action='store_true',help="draw heatmap for intensities of all regions")
    group.add_argument("-A","--Average",action='store_true',help="draw average pattern for intensities of all regions")
    p.add_argument('-b','--bams',nargs='+',dest="bams",type=str,help="the list of bam files containing raw reads information")
    p.add_argument('-i','--intervals',nargs='+',dest='intervals',type=str,help="the list of bed files containing location information of intervals")
    p.add_argument('-d','--direction',action='store_true',help='consider direction of each interval, useful for TSS around regions (default: False)')
    p.add_argument('-r','--resolution',type=int,default=5,help='the resolution for counting reads (default: 5)')
    p.add_argument('-l','--length',type=int,default=1000,help='the length extending from the center to both directions to be drawn (default: 1000)')
    p.add_argument('-f','--frag_l',nargs='+',type=int,default=[300],help='the average length of sequencing fragments, used to determine center location of each read (default: 300, 150 for MNase-seq)')
    p.add_argument('-n','--bamname',nargs='+',dest='bamnames',type=str,help="names for the bam files to be shown on the figure")
    p.add_argument('-N','--intervalname',nargs='+',dest='intervalnames',type=str,help='names for the bam files to be shown on the figure')
    p.add_argument('-w','--win_l',type=int,default=3,help='smooth window length for counts in each interval, (default:3, no smooth)')
    p.add_argument('-m','--method_c',type=str,default='kmeans',help='clustering method for first heatmap, can be: kmeans, somcuster, hcluster[hcluster only for small size intervals]. default: kmeans')
    p.add_argument('-o','--output',type=str,default='test',help='suffix of output figure file,can be (pdf, eps, png,jpg,...) final will be average_* or heatmap_*')
    if len(sys.argv)==1:
        print >> sys.stderr, p.print_help()
        sys.exit(0)
    return p.parse_args()

def find_center(read,half_len):
    '''
    find center of each fragment from read
    '''
    if not read.is_reverse:
        center=read.pos+half_len
    elif read.is_reverse:
        center=read.aend-half_len
    else:
        print >>sys.stderr,read.id,"No strand"
    return center

def smooth(x,window_len=11,window='hanning'):
        '''
        smooth the count for each interval.
        from:
        http://stackoverflow.com/questions/5515720/python-smooth-time-series-data
        '''
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=mp.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


# use center 50bp of the reads
# use center 50bp of the reads
def get_count(interval,bam,r,l,fl,direction=False,win_l=3):
    '''
    generate counts for +-<l> length regions from center of interval
    with <r> resolution. For each read in bam we only use the center
    of corresponding fragment to count (fragment size <fl>)
    '''
    half_len=int(fl/2)
    counts=[]
    for i in interval:
        count=np.zeros(int(2*l/r))
        center=int((i.start+i.stop)/2)
        if 'random' in i.chr: continue
        if i.chr=='chrM': continue
        for j in bam.fetch(i.chr,center-l-(half_len),center+l+(half_len)):
            j_center=find_center(j,half_len)
            bin=int((j_center-center+l)/r)
            low_bound=min(int(2*l/r),max(0,bin-int(50/r)/2))
            high_bound=max(0,min(int(2*l/r),bin+int(50/r)/2))
            count[low_bound:high_bound]+=1
        if direction and i.strand=='-':
            count=count[::-1]
        elif direction and i.strand!='+':
            print >> sys.stderr, '##ERROR: Strand information not there, please omit -d option'
            sys.exit(0)
        count=smooth(count,window_len=win_l)
        counts.append(count)
        if len(counts)%1000 ==0:
            print >> sys.stderr, "    ##  Counting for interval:",len(counts),"\r",
    print
    return np.array(counts)

global cdict
cdict  = {'red':  ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 1.0, 1.0),
                   (1.0, 0.8, 0.8)),

         'green': ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.8, 0.8),
                   (0.1, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

def feature_scale(matrix):
    '''
    Feature scaling and standardization
    x=(x-mean(X))/(max(X)-min(X))
    '''
    matrix=(matrix-np.mean(matrix))/(np.max(matrix)-np.min(matrix))
    
    return matrix


def heatmap_oneBam(collect,fig,xstart,xlen,cum_interval_n,leng,name):
    axmatrix=fig.add_axes([xstart,0.1,xlen,0.8])
    plt.register_cmap(name='cust_cmap', data=cdict)
    collect=np.minimum(collect,np.mean(collect)+3*np.std(collect)) # remove outlier with are extremely large
    im=axmatrix.matshow(collect,aspect='auto',origin='lower',cmap=plt.get_cmap('Reds'))
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    axmatrix.set_title(name)
    
    for n in cum_interval_n:
        axmatrix.axhline(y=n,color='black')
        
    
    #xticks
    axlabel=fig.add_axes([xstart,0.08,xlen,0.0])
    axlabel.set_xticks([-leng,0,leng])
    axlabel.set_yticks([])
    #plot colorbar
    axcolor=fig.add_axes([xstart+xlen*1.01,0.1,xlen*0.04,0.8])
    pylab.colorbar(im,cax=axcolor)




def main():
    args=ParseArg()
    frag_l=args.frag_l
    # determin if intervals or bams is multiple
    if len(args.intervals)==1 and not args.intervalnames: 
        interval_names=['interval']
    else: 
        interval_names=args.intervalnames

    if len(args.bams)==1 and not args.bamnames: 
        bam_names=['bam']
    else: 
        bam_names=args.bamnames
    if len(args.bams)!=len(bam_names) or len(args.intervals)!=len(interval_names):
        print >>sys.stderr,"length of names are not matching length of files, please check your command"
        sys.exit(0)
    if len(args.bams)!=len(frag_l):
        print >>sys.stderr,"number of BAMs are not matching number of frag_l, assign all to 300bp"
        frag_l=np.array([300]*len(args.bams))


    # store bam files with indexing and count information
    bams={}
    read_numbers={}
    print >> sys.stderr, "##  Starting read and index BAM files: "
    for i in range(len(args.bams)):
        temp_name=bam_names[i]
        print >> sys.stderr,"  ##  Indexing for bam file of '"+temp_name+"'"
        bams[temp_name]=pysam.Samfile(args.bams[i],'rb')
        print >> sys.stderr,"    ##  counting total reads number <slow>"
        ss=0
        for chr in bams[temp_name].references:
            ss+=bams[temp_name].count(chr)
        read_numbers[temp_name]=ss
        print >> sys.stderr

    # store interval files
    intervals={}
    print >> sys.stderr, "##  Starting reading intervals:"
    for i in range(len(args.intervals)):
        temp_name=interval_names[i]
        print >> sys.stderr,"  ##  Reading for interval file of '"+temp_name+"'\r",
        intervals[temp_name]=TableIO.parse(args.intervals[i],'bed')
        print >> sys.stderr

    resol=args.resolution
    leng=args.length
    # draw heatmap
    if args.Heatmap:
        print >> sys.stderr,"## Start count reads"
        collects={}
        interval_n=[0]
        order={}
        for k,nab in enumerate(bam_names):
            collect=[]
            for l,name in enumerate(interval_names):
                if k>=1:
                    intervals[name]=TableIO.parse(args.intervals[l],'bed')
                print >> sys.stderr, "  ## counting for bam["+nab+"] - interval["+name+"]"
                H_counts=get_count(intervals[name],bams[nab],resol,leng,frag_l[k],args.direction,args.win_l)
                H_counts=H_counts*5E7/read_numbers[nab]
                H_counts=np.log(H_counts+1)
                #H_counts=feature_scale(H_counts)
                if k==0:
                    if args.method_c=='kmeans':
                        centroids,_=kmeans(H_counts,5)
                        idx,_=vq(H_counts,centroids)
                        print >> sys.stderr,"  ## size of clusters using kmeansfor bam["+nab+"] - interval["+name+"]: "
                        cluster_size=Counter(idx)
                        for c in cluster_size: 
                            print >> sys.stderr,"   Cluster[%d]:%d"%(c,cluster_size[c])
                        order[name]=[i[0] for i in sorted(enumerate(idx), key=lambda x:x[1])]
                    elif args.method_c=='somcluster':
                        clusterid,_=somcluster(data=H_counts,nxgrid=5,nygrid=5)
                        order[name] = np.lexsort((clusterid[:,1],clusterid[:,0]))
                    elif args.method_c=='hcluster':
                        distMatrix = dist.pdist(H_counts)
                        distSquareMatrix = dist.squareform(distMatrix)
                        linkageMatrix = hier.linkage(distSquareMatrix)
                        dendro = hier.dendrogram(linkageMatrix)
                        order[name] = dendro['leaves']
                    interval_n.append(H_counts.shape[0])
                H_counts=H_counts[order[name],:]
                collect.append(H_counts)                
            collect=np.vstack(collect)
            collects[nab]=collect
        cum_interval_n=np.cumsum(interval_n)

        
        fig=pylab.figure(figsize=(3*len(bam_names),8))
        print >> sys.stderr,"## Start draw heatmap for intereval"
        #heatmap
            #yticks
        aylabel=fig.add_axes([0.1,0.1,0.0,0.8])
        aylabel.set_yticks(cum_interval_n)
        aylabel.set_xticks([])
        
        j=0
        width=0.8/len(bam_names)
        for name in bam_names:
            heatmap_oneBam(collects[name],fig,0.15+j*width,0.7*width,cum_interval_n,leng,name)
            j=j+1
        fig.savefig('heatmap_'+args.output)

        
        # print clustering information for first bam
        cluster_info=open("interval_with_cluster.txt",'w')
        name=interval_names[-1]
        intervals[name]=TableIO.parse(args.intervals[l],'bed')
        n=0
        for l in intervals[name]:
            if 'random' in l.chr: continue
            try:
                print >>cluster_info,'\t'.join(str(f) for f in [l.chr,l.start,l.stop,idx[n]])
            except:
                print "Error: the number of intervals are not consistent" 
            n=n+1
                
    # draw averge patterns
    if args.Average:
        print >> sys.stderr,"##  Start count reads"
        collect={}
        y_max=0
        for k,nab in enumerate(bam_names):
            for j,name in enumerate(interval_names):
                if k>=1:
                    intervals[name]=TableIO.parse(args.intervals[j],'bed')
                print >> sys.stderr, "  ## counting for bam["+nab+"] - interval["+name+"]"
                H_counts=get_count(intervals[name],bams[nab],resol,leng,frag_l[k],args.direction,args.win_l)
                print H_counts.shape[0]
                if name=='interval':
                    collect[nab]=np.sum(H_counts,axis=0)/H_counts.shape[0]*5E7/read_numbers[nab]
                    y_max=max(y_max,max(collect[nab]))
                else:
                    collect[(nab,name)]=np.sum(H_counts,axis=0)/H_counts.shape[0]*5E7/read_numbers[nab]
                    y_max=max(y_max,max(collect[nab,name]))
        
        fig=plt.figure(figsize=(3*len(bam_names),4))
        for i,nab in enumerate(bam_names):
            if len(interval_names)>1:
                ax = plt.subplot2grid((1,len(bam_names)),(0,i))
                for j,name in enumerate(interval_names):
                    col=matplotlib.cm.Paired((j+1)*1.0/(len(interval_names)+2),1)
                    ax.plot(np.array(range(-leng,leng,resol))+resol/2.0,collect[(nab,name)],color=col)
                ax.legend(interval_names,loc='upper right')
                ax.set_ylim(0,y_max+1)
                ax.set_title(nab)
            else:
                col=matplotlib.cm.Paired((i+1)*1.0/(len(bam_names)+2),1)
                plt.plot(np.array(range(-leng,leng,resol))+resol/2.0,collect[nab],color=col)
        if len(interval_names)==1:
            plt.legend(bam_names,bbox_to_anchor=(0., 0.95, 1., .100), loc=3, ncol=len(bam_names), mode="expand", borderaxespad=0.,fontsize=5)#,prop={'size':15},fontsize=10)
            plt.ylim(0,y_max+1)
        plt.xlabel('Distance to center')
        plt.ylabel('Average coverage for 5E7 reads')
        plt.tight_layout()
        fig.savefig('average_'+args.output)
    



if __name__=="__main__":
    main()









