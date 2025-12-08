#!/usr/bin/python3 readsam2.py
'''updated 051224, to be used where ISS_plasmid and ISS_unintegrated not generated during mapping'''
'''updated 061424, now workpath is defined by input parameter'''


import re
import sys
import glob
import os
from pathlib import Path


def main():
  tag = sys.argv[1]
  folder = sys.argv[2]

  workpath = '/u/scratch/y/yuanshi/{}'.format(folder)

  for subfolder in ['ISS_archive','ISS_linkage']:
      Path(os.path.join(workpath,subfolder)).mkdir(parents=True, exist_ok=True)
  archfile  = open(os.path.join(workpath,'ISS_archive','record_{}.txt'.format(tag)),'w')
  linkfile  = open(os.path.join(workpath,'ISS_linkage','linkage_{}.txt'.format(tag)),'w')
  sumfile  = open(os.path.join(workpath,'ISS_archive','summary_{}.txt'.format(tag)),'w')

  #hg38files = sorted(glob.glob(workpath+'ISS_bowtie/hg38/'+tag+'.sam'))
  linkdict = {}; logdict = {'unmapped':0,'unpaired':0,'positive_strand':0,'negative_strand':0}
  infile = os.path.join(workpath,'ISS_bowtie','hg38','{}.sam'.format(tag))
  linkdict,logdict = readsam2(linkdict,logdict,infile,'hg38',archfile)
  infile = os.path.join(workpath,'ISS_bowtie','HIV','{}.sam'.format(tag))
  linkdict,logdict = readsam2(linkdict,logdict,infile,'HIV', archfile)

  #infile = os.path.join(workpath,'ISS_unintegrated','{}.txt'.format(tag))
  #linkdict,logdict = readtxt(linkdict,logdict,infile,'free',archfile)
  #infile = os.path.join(workpath,'ISS_plasmid','{}.txt'.format(tag))
  #linkdict,logdict = readtxt(linkdict,logdict,infile,'plasmid',archfile)
  
  infiles = glob.glob(os.path.join(workpath,'ISS_log','{}*log'.format(tag)))
  linkdict,logdict = readlog(linkdict,logdict,infiles,archfile)

  bclist = []; annotdict = {}; depth = 0; depdict = {}
  for UMI in linkdict:
    topbc,bctopfreq,bcsecfreq  = find_most([x[0] for x in linkdict[UMI]])
    topannot,atopfreq,asecfreq = find_most([x[1] for x in linkdict[UMI]])
    atype = topannot.rsplit(':')[0]
    topa  = topannot.rsplit('|')[0].rsplit(':')[1:]
    topa  = ':'.join(topa)
    mlength = topannot.rsplit('|')[-1]
    linkfile.write(UMI+'\t'+str(len(linkdict[UMI]))+'\t'+topbc+'\t'+str(bctopfreq)+'\t'+str(bcsecfreq)+'\t'+atype+'\t'+topa+'\t'+str(atopfreq)+'\t'+str(asecfreq)+'\t'+mlength+'\n')
    bclist.append(topbc)
    Udep  = len(linkdict[UMI])
    depth += Udep
    if atype not in annotdict: annotdict[atype] = 0
    annotdict[atype] += 1
    if Udep not in depdict: depdict[Udep] = 0
    depdict[Udep] += 1

  sumfile.write('total read: '+str(depth)+'\n')
  sumfile.write('total UMI: '+str(len(linkdict))+'\n')
  sumfile.write('unique barcodes: '+str(len(set(bclist)))+'\n')
  if 'hg38' in annotdict:
    sumfile.write('hg38 UMI: '+str(annotdict['hg38'])+'\n')
  if 'HIV' in annotdict:
    sumfile.write('HIV UMI: '+str(annotdict['HIV'])+'\n')
  if 'free' in annotdict:
    sumfile.write('unintegrated UMI: '+str(annotdict['free'])+'\n')
  if 'plasmid' in annotdict:
    sumfile.write('plasmid UMI: '+str(annotdict['plasmid'])+'\n')

  sumfile.write('unpaired reads:'+str(logdict['unpaired'])+'\n')
  sumfile.write('positive strand mapped reads:'+str(logdict['positive_strand'])+'\n')
  sumfile.write('negative strand mapped reads:'+str(logdict['negative_strand'])+'\n')
  sumfile.write('========UMI depth distribution=======\n')
  sumfile.write('UMI_read_depth\tunique_UMI_count\n')

  for Udep in depdict:
    sumfile.write(str(Udep)+'\t'+str(depdict[Udep])+'\n')
  sumfile.close()


def find_most(bclist):
  bcdict = {}
  for bc in bclist:
    if bc not in bcdict: bcdict[bc] = 0
    bcdict[bc] += 1
  sortbcdict = sorted((v,k) for (k,v) in bcdict.items())
  if len(sortbcdict) == 1:
    return sortbcdict[0][1], 1, 'NA'
  else:
    return sortbcdict[-1][1], float(sortbcdict[-1][0])/len(bclist), float(sortbcdict[-2][0])/len(bclist)



def readsam2(linkdict,logdict,infile,samtype,archfile):
  #This is for paired end read
  lc = 0
  samfile = open(infile)
  for line in samfile:
    if line[0] == '@': continue
    lc += 1
    if lc % 2 == 1:
      line1 = line
      continue
    line1   = line1.rstrip().rsplit('\t')
    line2   = line.rstrip().rsplit('\t')
    if line1[0] != line2[0]:
      line1 = '\t'.join(line[2])
      lc -= 1
      continue
    UMI     = line1[0].rsplit(':')[-2]
    bc      = line1[0].rsplit(':')[-1]
    #In the new design, R2 is forward read after LTR.
    flag1   = bin(int(line1[1]))
    flag2   = bin(int(line2[1]))
    if len(flag1) > 7 and flag1[-8] == '1':
      read2 = line1; flag2 = bin(int(line1[1]))
      read1 = line2; flag1 = bin(int(line2[1]))
    else:
      read1 = line1
      read2 = line2
    qual1   = int(line1[4])
    qual2   = int(line2[4])
    if qual2 >= qual1:
      if flag2[-3] == '1' or qual2 == 1:
        logdict['unmapped'] += 1
        continue
      strand = flag2strand(flag2)
      if strand == '+':
        pos = read2[3]
        cigar = read2[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
      else:
        cigar = read2[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
        pos = str(int(read2[3])+mlength)
    else:
      if flag1[-3] == '1' or qual1 == 1:
        logdict['unmapped'] += 1
        continue
      strand = flag2strand(flag1)
      if strand == '+':
        strand = '-'
        cigar = read1[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
        pos = str(int(read1[3])+mlength)
      else:
        strand = '+'
        pos = read1[3]
        cigar = read1[5]
        mlength = re.findall(r"(\d+)M", cigar)
        mlength = sum([int(x) for x in mlength])
    if strand == '+':
      logdict['positive_strand'] += 1
    else:
      logdict['negative_strand'] += 1
    chrom1 = read1[2]
    chrom2 = read2[2]
    if chrom1 != chrom2 and chrom1 != '*' and chrom2 != '*':
      logdict['unpaired'] += 1
      continue
    else:
      chrom = chrom2
    annot = chrom+':'+pos+':'+strand
    archfile.write(UMI+'\t'+bc+'\t'+samtype+'\t'+annot+'\n')
    if UMI not in linkdict:
      linkdict[UMI] = []
    linkdict[UMI].append((bc,samtype+':'+annot+'|'+str(mlength)))
  samfile.close()
  return linkdict,logdict


def flag2strand(flag):
  if len(flag) > 5 and flag[-5] == '1':
    return '-'
  else:
    return '+'

def readtxt(linkdict,logdict,infile,txttype,archfile):
  txtfile = open(infile)
  for line in txtfile:
    line    = line.rstrip('\n').rsplit('\t')
    bc      = line[1]
    UMI     = line[0]
    annot   = line[2]
    archfile.write(UMI+'\t'+bc+'\t'+txttype+'\t'+annot+'\n')
    if UMI not in linkdict:
      linkdict[UMI] = []
      linkdict[UMI].append((bc,txttype+':'+annot+'|'+str(len(annot))))
  txtfile.close()
  return linkdict,logdict

def readlog(linkdict,logdict,infiles,archfile):
  for infile in infiles:
    with open(infile) as f:
            next(f)
            for line in f:
              try:
                _,flag1,_,flag2,bc,_,UMI,_,_,annot=line.split()
              except:
                continue
              if (flag1.startswith('221') or flag1.startswith('222')) and flag2[0]=='2' and flag2[1]!='0': txttype = 'unintegrated'
              elif (flag1.startswith('221') or flag1.startswith('222')) and flag2[0]=='1' and flag2[1]!='0': txttype = 'plasmid'
              else: continue
              archfile.write(UMI+'\t'+bc+'\t'+txttype+'\t'+annot+'\n')
              if UMI not in linkdict:
                linkdict[UMI] = []
                linkdict[UMI].append((bc,txttype+':'+annot+'|'+str(len(annot))))
  return linkdict,logdict

if __name__ == '__main__':
    main()
