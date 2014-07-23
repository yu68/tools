#!/usr/bin/env python
# Pengfei Yu yupf05@gmail.com

import sys,os,argparse
from Bio.Blast import NCBIXML
import itertools
from Bio import SeqIO
from time import time


'''
dir(hsp):
['__class__', '__delattr__', '__dict__', '__doc__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', 'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
'''


def ParseArg():
    
    p=argparse.ArgumentParser( description = 'DESCRIPTION: Run BLAST and select positive results from a large set of fastq or fasta sequences to find existence of specific sequences', epilog='')
    p.add_argument("input",type=str,help="the input fastq/fasta file containing reads sequences")
    p.add_argument("-e","--evalue",dest="evalue",type=float,default=0.00001,help="cutoff evalues, only choose alignment with evalue less than this cutoffs (default: 1e-5).")
    p.add_argument("--seq_db",dest="seq_db",type=str,help="BLAST database of specific sequences, or fasta sequences",default="~/Stitch-seq/blast_db/linker.fa")
    p.add_argument("--blast_path",dest="blast_path",type=str,help="path for the local blast program",default="/home/yu68/Softwares/ncbi-blast-2.2.27+/bin/blastn")
    p.add_argument("-o","--output",dest="output",type=str,help="output file contains sequences with both linker alignment")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        exit(0)
    return p.parse_args()

def blast_align(fasta,blast_path,linker_db):
    os.system(blast_path+" -task blastn -outfmt 5 -num_threads 6 -evalue 0.1 -db "+linker_db+" -query ./temp/"+fasta+" > "+name+"temp_blast_linker.xml")
    linker_records=NCBIXML.parse(open(name+"temp_blast_linker.xml"))
    os.system("rm ./temp/"+fasta)
    return (linker_records)

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch



def main():
    #initialization
    n=0 # total number of query seq
    align_linker=0 # aligned_reads
    counts={}  # store number of alignment to different linkers

    args=ParseArg()
    output=open(args.output,'w')
 
    os.system("mkdir temp") # create temp folder
    
    global name
    name=args.output.split("_")[0]
  
    blast_path=args.blast_path
    linker_db=args.seq_db

    # create blast database if not yet
    if not os.path.isfile(linker_db+".nhr"):
        print >> sys.stderr, "  ## blast database not detected, making blast db for seq"
        blastdir,_ = os.path.split(blast_path)
        os.system(blastdir+"/makeblastdb -in "+linker_db+" -dbtype 'nucl' -title "+os.path.splitext(linker_db)[0])
    
    # E-values
    evalue=float(args.evalue)
    
    # determine input sequence file type
    types="fastq"
    if args.input.split(".")[-1] in ["fa","fasta"]:
        types="fasta"

    seq_file=SeqIO.parse(args.input,types)
    
    ###################################
    ##    start editing from here    ## 
    ###################################
    
    for i, batch in enumerate(batch_iterator(seq_file, 100000)):
        t0=time()
        filename=name+"group_%i.fasta" % (i+1)
        handle=open("./temp/"+filename, "w")
        count=SeqIO.write(batch,handle,"fasta")
        handle.close()
        print "Wrote %i records to %s" % (count,filename)
        
        linker_records = blast_align(filename,blast_path,linker_db)
        print "BLAST aligned for %s." % (filename)
        
        
        print "Start to parse BLAST results for %s" %(filename)
        for linker_record in linker_records:
            temp_output=''
            n=n+1
            line=''
            align=False
            temp_output=linker_record.query+'\t'
            for alignment in linker_record.alignments:
                expect=evalue
                temp_output=linker_record.query+'\t'
                for hsp in alignment.hsps:
                    if hsp.expect < expect:
                        expect=hsp.expect
                        line="\t".join (str(f) for f in [hsp.query_start,hsp.query_end,alignment.title,hsp.sbjct,hsp.sbjct_start,hsp.sbjct_end,hsp.expect,hsp.score])
                temp_output=temp_output+line+'\n'
                
                if expect<evalue:
                    output.write(temp_output)
                    if alignment.hit_def in counts.keys():
                        counts[alignment.hit_def]+=1
                    else:
                        counts[alignment.hit_def]=1
                    align=True
            if align==True: 
                align_linker+=1
        stat=''
        for i in counts.keys():
            stat=stat+'\n'+i+':%d'%(counts[i])
        t1=time()
        print "After %s, got %i sequences, %i align to linkers.\n %s. \nThis batch takes %.2f min."%(filename,n,align_linker,stat,(t1-t0)/60)

    output.close()
    os.system("rm *.xml")
    os.system("rm -r temp")

if __name__=="__main__":
    main()
