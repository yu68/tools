#normalized by total_read_no/5000001
 
import sys
import csv
from getopt import getopt

def show_help():
    print >>sys.stderr,"sortbed2wig.py: convert sorted BED file into WIG file with 1bp resolution"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Library dependency: csv\n\n"
    print >>sys.stderr,"Usage: sortbed2wig.py <options> -i input.bed -n name_of_output -e -l extended_read_length -s column_num_for_strand"
    print >>sys.stderr,"       sortbed2wig.py <options> -i input.bed -n name_of_output -e"
    print >>sys.stderr,"       sortbed2wig.py <options> -i input.bed -n name_of_output"
    print >>sys.stderr,"Example:  sortbed2wig.py <options> -i mm9_H3K9me3.bed -n mm9_H3K9me3 -e -l 150 -s 4"
    print >>sys.stderr,"Options:"
    print >>sys.stderr,"   -h,--help          show help information"
    print >>sys.stderr,"   -i,--inputfile     input bed file (with strand information for extend option)" 
    print >>sys.stderr,"   -o,--outputFolder  folder for output wid file (default: /home/GenomeBrowser/lab_tracks/"
    print >>sys.stderr,"   -n,--wigname       name of the output wig file"
    print >>sys.stderr,"   -e,--extend        extend read in bed file or not (default: false) "
    print >>sys.stderr,"   -l,--readlength    the extended length of each read (default: 150, effective only when extend=True) "
    print >>sys.stderr,"   -s,--strandLoc     the column # for strand information in bed (default: 4, effective only when extend=True)"

Version="1.0"
if len(sys.argv)<2:
    show_help()
    exit(0)
opts,restlist=getopt(sys.argv[1:],"hi:o:n:el:s:",["help","inputfile=","outputFolder=","wigname=","extend","readlength=","strandLoc="])
extend=False
read_len=150
strandLoc=4
outputFolder="/home/GenomeBrowser/lab_tracks/"
for o,a in opts:
  if o in ("-h","--help"):
      show_help()
      exit(0)
  elif o in ("-i","--inputfile"):
      file_name=a
  elif o in ("-n","--wigname"):
      name=a
  elif o in ("-o","--outputFolder"):
      outputFolder=a
  elif o in ("-e","--extend"):
      extend=True
  elif o in ("-l","--readlength"):
      read_len=int(a)
  elif o in ("-s","--strandLoc"):
      strandLoc=int(a)
output = csv.writer(open(outputFolder+name+".wig", 'wb'), delimiter='\t')
rowtotal=[]
rowtotal.append('track type=wiggle_0 name='+name+' maxHeightPixels=64:64:11 color=120,31,180 visibility=full')
output.writerow(rowtotal)
file_obj = open(file_name,'r')
chr_pos_st =file_obj.read().split('\n')
file_obj.close()
total_no =len(chr_pos_st) # total number of rows
factor =total_no/5000000.0


#generate wig file for each chromosome
def wigonechro(start,end):
  chr_pos_st1=chr_pos_st[start:end]
  a=[['chr1','0','0','+']]
  for i in chr_pos_st1:
    a.append(i.split('\t'))
  d=range(0,len(a))
  if extend:
    for i in range(1,len(a)):
      if a[i][strandLoc-1]=='-':
        d[i]=max(int(a[i][2])-read_len+1,1)
      else:
        d[i]=max(int(a[i][1]),1)
      e=d[i]
      f=1
      while e<d[i-f]:
        d[i+1-f]=d[i-f]
        d[i-f]=e
        f=f+1
  else:
    for i in range(1,len(a)): 
      d[i]=max(int(a[i][1]),1)
  n=1
  end=d[1]+read_len-1
  d[len(a)-1]=1000000000
  for i in range (1,len(a)-1):
    b=d[i-1]
    c=d[i]
    if c<d[i+1]:
      reads_no=1
      while c>=b>=c-read_len+1:
        reads_no=reads_no+1
        b=d[i-reads_no]
      rowtotal=[]
      rowtotal.append(c)
      rowtotal.append("%.2f"%(reads_no/factor))
      output.writerow(rowtotal)
      while c<end<=d[i+1]:
        if (d[n]<d[n+1])&(i!=n)&(end!=d[i+1]):
          rowtotal=[]
          rowtotal.append(end)
          rowtotal.append("%.2f"%((i-n)/factor))
          output.writerow(rowtotal)
        n=n+1
        end=d[n]+read_len-1
  return
   
start=0
for i in range(1,total_no):
  chro2=chr_pos_st[i].split('\t')[0]
  chro1=chr_pos_st[i-1].split('\t')[0]
  if chro1!=chro2:
    end=i
    rowtotal=[]
    rowtotal.append('variableStep chrom='+chro1+' span=%d'%(read_len))
    output.writerow(rowtotal)
    wigonechro(start,end)
    start=i
  if i%100000==0:
    print 'total number:', total_no, 'processed:', i



