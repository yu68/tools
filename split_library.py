import sys,os,argparse
from Bio import SeqIO

def ParseArg():
    '''Parse the argument'''
    p=argparse.ArgumentParser( description = 'Example: %(prog)s -q -i example.fastq -b barcode.txt', epilog = 'Library dependency: Bio')
    group=p.add_mutually_exclusive_group()
    group.add_argument("-f","--fasta",action='store_true',help='add this option for fasta input file')
    group.add_argument("-q","--fastq",action='store_true',help='add this option for fastq input file')
    p.add_argument('-v','--version',action='version',version='%(prog)s 0.1')
    p.add_argument('-i','--input',dest='input',type=str,help='input fastq/fasta file')
    p.add_argument('-b','--barcode',dest='barcode',type=str,help='barcode file')
    p.add_argument('-t','--trim',action='store_true',help='trim sequence before and within barcode')
    p.add_argument('-m','--max_score',dest='max_score',type=int,default=2, help="max(mismatch+indel) allowed for barcode match, otherwise move reads into 'unassigned' file") 
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()



def fuzzy_substring(needle, haystack):
    """Calculates the fuzzy match of needle in haystack,
    using a modified version of the Levenshtein distance
    algorithm.
    The function is modified from the levenshtein function
    in the bktree module by Adam Hupp
    http://ginstrom.com/scribbles/2007/12/01/fuzzy-substring-
    matching-with-levenshtein-distance-in-python/"""
    m, n = len(needle), len(haystack)

    # base cases
    if m == 1:
        return not needle in haystack
    if not n:
        return m

    row1 = [0] * (n+1)
    minS=m
    for i in range(0,m):
        row2 = [i+1]
        for j in range(0,n):
            cost = ( needle[i] != haystack[j] )

            row2.append(   min(row1[j+1]+1, # deletion
                               row2[j]+1, #insertion
                               row1[j]+cost) #substitution
                          )
            if i == m-1:
                if row2[j+1] <= minS:
                    minS=row2[j+1]
                    end=j+1
        row1 = row2
    return minS, end
'''
TEST:
print (fuzzy_substring('ACTC', 'C_ ATCG'))
print (fuzzy_substring('ACTC', 'C_ ACTGG'))
print (fuzzy_substring("ACTAAC", "ACTAACTAGCCATGCAATGGCTAG"))
'''

def Main():
    args=ParseArg()
    if args.fastq:
       type="fastq"
    elif args.fasta:
       type="fasta"
    seq={}
    #----------- read barcode ----------
    barcodes=[]
    for i in open(args.barcode,'r'):
        i=i.strip()
        barcodes.append(i)
        barcode_len=len(i)
        seq[i]=[]
    seq['unassign']=[]
    #-----------------------------------
    
    print "start to assign sequence to different barcodes..."    
    print "----------"
    for record in SeqIO.parse(args.input,type):
        miScore=barcode_len
        print record
        for i in barcodes:
            score,j=fuzzy_substring(i,record.seq)
            if score<miScore:
                barcode=i
                end=j
                miScore=score
            print score,i
        if miScore>args.max_score:
            seq['unassign'].append(record)
        else:
          if args.trim:
              seq[barcode].append(record[end:])
          else:
              seq[barcode].append(record)
    barcodes.append('unassign')
    for barcode in barcodes:
        print "#Writing "+type+" file for barcode: "+barcode+". # of seq: %i ."%len(seq[barcode])
        output=open(barcode+args.input,'w')
        SeqIO.write(seq[barcode],output,type)

        
if __name__=="__main__":
    Main()

