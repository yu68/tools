import os,sys,argparse
from random import choice,randrange,gauss
from cruzdb import Genome

def ArgParse():
    p=argparse.ArgumentParser( description = "generate random sequences with customized length and number (for random peaks et...)", epilog="library dependency: cruzdb (https://github.com/brentp/cruzdb),sqlalchemy")
    p.add_argument("-g",'--genome',dest='genome',default='mm9',type=str,help='specify genome name to get chromosome info from UCSCGB, default: mm9 ')
    p.add_argument('-m','--mean',dest='mean',type=float,default=200,help='mean length of each random sequence,default:200')
    p.add_argument('-s','--sd',dest='sd',type=float,default=20,help='sd of random sequence lengths,default:20')
    p.add_argument('-n','--num',dest='num',type=int,default=10000,help='number of sequences to be randomly sampled,default:10000')
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def Main():
    args=ArgParse()

    # initiate genome object from UCSC genome browser
    chromInfo=Genome(db=args.genome).chromInfo
    
    chroms={}
    #store chromInfo 
    for i in range(40):
        try: 
            chroms[chromInfo[i].chrom]=chromInfo[i].size
        except:
            break
    print >>sys.stderr, "Chromosome information readed, %d chromosomes"%(len(chroms))
    i=0
    while i < args.num:
        # randomly select one chromosome and a region in this chromosome
        chrom=choice(chroms.keys())
        size=chroms[chrom]
        if 'random' in chrom: continue
        length=int(gauss(args.mean,args.sd))
        start=randrange(1,size-length-1)
        end=start+length-1
        strand=choice(['+','-'])
        print "\t".join(str(f) for f in [chrom,start,end,i+1,0,strand])
        i=i+1

if __name__=="__main__":
    Main()
