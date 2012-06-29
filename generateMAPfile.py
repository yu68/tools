'''
Usage:
  python generateMAPfile.py chromInfo_hg18.txt cancer_indel_dup_inv.vcf > map_chain.txt 

'''

import sys
import vcf
import numpy as np
import copy
from bedtools import IntervalFile
from getopt import getopt

def show_help():
    print >>sys.stderr,"generateMAPfile.py: generation of MAP(chain) file for conversion between cancer and ref genome"
    print >>sys.stderr,"Version:"+Version+"\n"
    print >>sys.stderr,"Library dependency: vcf, numpy, copy, bedtools"
    print >>sys.stderr,"Usage:   generateMAPfile.py -l [chromLength_file] -i [cancer_vcf_file(indel_dup_inv)] > output_chain_file"
    print >>sys.stderr,"Example:  generateMAPfile.py -l chromInfo_hg18.txt -i cancer_indel_dup_inv.vcf > map_chain.txt"
    print >>sys.stderr,"Option:"
    print >>sys.stderr,"    -h,--help"
    print >>sys.stderr,"    -l,--chromLen   file containing information of chromosome length for specific genome"
    print >>sys.stderr,"                    format:"
    print >>sys.stderr,"                        chr1	1	247249719"
    print >>sys.stderr,"                        chr10	1	135374737"
    print >>sys.stderr,"                        ......"
    print >>sys.stderr,"    -i,--inputVcf   standard vcf file containing information of indels, duplications and inversions for cancer"

Version="1.0"
if len(sys.argv)<2:
    show_help()
    exit(0)
opts,restlist=getopt(sys.argv[1:],"hl:i:",["help","chromLen","inputVcf="])
for o,a in opts:
  if o in ("-h","--help"):
      show_help()
      exit(0)
  elif o in ("-l","chromLen"):
      chro_length=open(a)
  elif o in ("-i","inputVcf"):
      vcf_name=a


#chro_length=open(sys.argv[1])
#vcf_name=sys.argv[2]
vcf_data=vcf.VCFReader(open(vcf_name))
open(vcf_name).close()
#collect and store indels and tandom duplications information
chro_length=chro_length.read().split('\n')

indicator=[]
for i in vcf_data:
    if len(i.CHROM)>2:
        sv_chro=i.CHROM
    else:
        sv_chro="chr"+i.CHROM
    sv_pos=i.POS
    sv_ref=i.REF
    alt=i.ALT[0]
    if alt.startswith('<'):
        svtype=i.INFO['SVTYPE']
        length=i.INFO['SVLEN']
        end=i.INFO['END']
    else:
        svtype='none'
        length=len(i.ALT)-len(sv_ref)

    if (svtype!='none')&(svtype!='INV'):
        indicator.append((sv_chro,sv_pos,length,svtype))
    elif svtype=='INV':
        indicator.append((sv_chro,end,end-sv_pos+1,svtype))
    elif (svtype=='none')&(length>0):
        indicator.append((sv_chro,sv_pos,length,"INS"))
    elif (svtype=='none')&(length<0):
        indicator.append((sv_chro,sv_pos,length,"DEL"))

for j in chro_length:
    if j!='':
        k=j.split('\t')
        indicator.append((k[0],int(k[1]),0,'STA'))
        indicator.append((k[0],int(k[2]),0,'END'))

dtype=[('chro','S8'),('start','i4'),('length','i4'),('type','S3')]
indi_array=np.array(indicator,dtype=dtype)
indi_array=np.sort(indi_array,order=['chro','start'])

'''
for i in indi_array:
    print "\t".join(str(f) for f in i)
'''

#function for inversion covers previous segments

def changeUponInversion(window_num,pos,length,pos2,start,start2,map_file,map_window):
    k=len(map_file)-1
    while map_file[k][1]>pos-length:
        k=k-1
    start_inv=pos-length+1
    if map_file[k][7]=='-':
        window_num=window_num+1
        map_window.append([chro,start_inv-99,start_inv-1,"+","chrW"+str(window_num),1,99,"+"])
        map_window.append([chro,pos-98,pos,"+","chrW"+str(window_num),100,198,"-"])
        map_window.append([chro,start_inv,start_inv+98,"+","chrW"+str(window_num+1),1,99,"-"])
        map_window.append([chro,pos+1,pos+99,"+","chrW"+str(window_num+1),100,198,"+"])
        window_num=window_num+1
    else:
        start2_inv=start_inv-map_file[k][2]+map_file[k][6]
        insert=copy.deepcopy(map_file[k])
        insert[1]=start_inv
        insert[5]=start2_inv
        map_file.insert(k+1,insert)
        map_file[k][2]=start_inv-1
        map_file[k][6]=start2_inv-1
        map_file.append([chro,start,pos,"+",chro,start2,pos2,"+"])
        for i in range(k+1,len(map_file)):
            temp=map_file[i][5]
            map_file[i][5]=start2_inv+pos2-map_file[i][6]
            map_file[i][6]=start2_inv+pos2-temp
            temp=map_file[i][7]
            if temp=='+':
                map_file[i][7]='-'
            else:
                map_file[i][7]='+'
        start=pos+1
        start2=pos2+1
    return (window_num,start,start2,map_file,map_window)
    



map_file=[]
start=start2=1
map_window=[]
window_num=0
for i in range(1,len(indi_array)):
    info=indi_array[i]
    chro=info[0]
    pos=info[1]
    length=info[2]
    typ=info[3]
    pos2=pos-start+start2
    if typ=='STA':
        start=start2=1
        continue
    elif typ=='INS':
        map_file.append([chro,start,pos,'+',chro,start2,pos2,'+'])
        map_file.append([chro,pos+1,pos+1,'+',chro,pos2+1,pos2+1+length,'+'])
        start=pos+2
        start2=pos2+2+length
    elif typ=='DEL':
        map_file.append([chro,start,pos,'+',chro,start2,pos2,'+'])
        map_file.append([chro,pos+1,pos+1-length,'+',chro,pos2+1,pos2+1,'+'])
        start=pos+2-length
        start2=pos2+2
    elif typ=='DUP':
        map_file.append([chro,start,pos+length,'+',chro,start2,pos2+length,'+'])
        map_file.append([chro,pos+1,pos+length,'+',chro,pos2+length+1,pos2+2*length,'+'])
        start=pos+length+1
        start2=pos2+2*length+1
    elif typ=='END':
        map_file.append([chro,start,pos,'+',chro,start2,pos2,'+'])
    elif typ=='INV':
        if start<=pos-length:
            map_file.append([chro,start,pos-length,'+',chro,start2,pos2-length,'+'])
            map_file.append([chro,pos-length+1,pos,'+',chro,pos-length+1,pos2,'-'])
            start=pos+1
            start2=pos2+1
        else:
            window_num,start,start2,map_file,map_window=changeUponInversion(window_num,pos,length,pos2,start,start2,map_file,map_window)

for i in map_file:
    print "\t".join(str(f) for f in i)

for i in map_window:
    print "\t".join(str(f) for f in i)


