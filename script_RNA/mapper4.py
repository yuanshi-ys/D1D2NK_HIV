#!/usr/bin/python3 mapper4.py
'''mapper updated on 062823, based on mapper3.py which is designed to map amplicon longer than sequenced, 3' and 5' on R1 and R5 has overhangs
   new added functions include 
   1) decorator to time how long the function runs
   2) feed in flag parameters to let user fine-tune
   
   updated on 032124: output foldername change add HIV'''
'''modified 10022024, add a binary flag, if set as true, the input is in .gz format'''


import os, math, json, operator, sys
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import time, argparse
import argparse
from pathlib import Path
import gzip

def timeit_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper

def hamming(seq1,seq2):
    if len(seq1)!=len(seq2):
        return -1
    return sum(pos1 != pos2 for pos1, pos2 in zip(seq1, seq2))

def error_correction(seq1,seq2,qual=30):
    '''error correction with two sequence, for each position if nucleotide the same, pick the concensus, if different compare with quality score:
if one is higher than or equal to 30, the other is lower than 30, pick the higher one
otherwise is N (discard later)'''
    if seq1.seq==seq2.seq:
        return str(seq1.seq)
    result=''
    for i,j,k,l in zip(seq1,seq2,seq1.letter_annotations['phred_quality'],seq2.letter_annotations['phred_quality']):
        if i==j:
            result+=i
            continue
        if k>30 and l<=qual:
            result +=i
            continue
        if l>30 and k<=qual:
            result +=j
            continue
        result += 'N'
    return result

def offset_check (source, sequence,trial_start=0,**kwargs): #need to set win, cutoff, verbose
    '''scanning from -win to +win around trial_start, map source with sequence, can define cutoff'''
    for mis in range(-args.win,args.win,1):#allowed scanning window -6 to +6
        if trial_start+mis<0:
            continue
        if trial_start+mis+len(source)>len(sequence): 
            continue
        target=sequence[trial_start+mis:trial_start+mis+len(source)]
        if args.verbose:
            print (source,target)
        dis = hamming(source,target)
        cutoff=kwargs.get('cutoff','NA')
        if cutoff == 'NA': cutoff=int(len(source)*0.2)
        if cutoff>=0.5*len(source): #too flexible
            cutoff=int(len(source)*0.4)
            #print ('Warning: mapping sequence used is too short, cutoff adjusted.',file=sys.stderr)
        if dis<=cutoff and dis>=0:
            if args.verbose:
                print('mapped',source,target,dis)
            return (mis) #return relative position from trial_start that got feed into the function
    return ('NA')

def check_strand(R1,R2,primers,**kwargs):
    trial_start=0
    offsetF=offset_check(primers['concensus5'][:15],R1,trial_start=trial_start,**vars(args))#use first 15nt for mapping
    offsetR=offset_check(primers['concensus5'][:15],R2,trial_start=trial_start,**vars(args))
    if offsetF == 'NA' and offsetR =='NA':
        return ('N') #not mapped
    if offsetF !='NA' and offsetR != 'NA':
        if args.verbose:
            print (offsetF,offsetR)
        return ('A') #ambiguous
    if offsetF != 'NA':
        return ('F')
    if offsetR != 'NA':
        return ('R')
    
def mapping(primers,sequence,overhang='5',bc_length=21,UMI_length=15,**kwargs):
    
    if overhang == '3': 
        seq5=primers['concensus5']
        if args.short3:
            seq3=primers['concensus3'][:args.short3]
        else:
            seq3=primers['concensus3']
    if overhang == '5': 
        seq5=primers['concensus5'][-args.short5:]
        seq3=primers['concensus3']
    
    use5=15 if len(seq5)>15 else len(seq5)
    usem=15 if len(primers['midconcensus'])>15 else len(primers['midconcensus'])
    use3=15 if len(seq3)>15 else len(seq3)
        
    trial_start=len(seq5)-use5
    if args.verbose:
        print ('Mapping 5\' concensus {}, using {}'.format(seq5,seq5[-use5:]) )
        print ('Target {}'.format(sequence[trial_start:]))
        print ('offset1 start position is {}'.format(trial_start))
    offset1=offset_check(seq5[-use5:],sequence,trial_start=trial_start,**vars(args))
    if args.verbose:
        print('offset1 relative position is {}'.format(offset1))
    if offset1 == 'NA':
        return ('P',1,use5,len(seq5)-use5,seq5,str(sequence),sequence[len(seq5)-use5:])
    offset1=offset1+trial_start+use5
    if args.verbose:
        print('offset1 real position is {}'.format(offset1))
    
    trial_start=offset1+bc_length
    if args.verbose:
        print ('Mapping middle 5\' concensus {}, using {}'.format(primers['midconcensus'],primers['midconcensus'][:usem]) )
        print ('Target {}'.format(sequence[trial_start:]))
        print ('offset2 start position is {}'.format(trial_start))
    offset2=offset_check(primers['midconcensus'][:usem],sequence,trial_start=trial_start,**vars(args))
    if args.verbose:
        print ('offset2 relative position is {}'.format(offset2))
    if offset2 == 'NA':
        return ('P',2)
    offset2=offset2+trial_start
    if args.verbose:
        print('offset2 real position is {}'.format(offset2))
    
    
    trial_start=offset1+bc_length+len(primers['midconcensus'])-usem
    if args.verbose:
        print ('Mapping middle 3\' concensus {}, using {}'.format(primers['midconcensus'],primers['midconcensus'][-usem:]) )
        print ('Target {}'.format(sequence[trial_start:]))
        print ('offset3 start position is {}'.format(trial_start))
    offset3=offset_check(primers['midconcensus'][-usem:],sequence,trial_start=trial_start,**vars(args))
    if args.verbose:
        print ('offset3 relative position is {}'.format(offset3))
    if offset3 == 'NA':
        return ('P',3)
    offset3=offset3+trial_start+usem
    if args.verbose:
        print('offset3 real position is {}'.format(offset3))

        
    trial_start=offset3+UMI_length
    if args.verbose:
        print ('Mapping 3\' concensus {}, using {}'.format(seq3,seq3[:use3]) )
        print ('Target {}'.format(sequence[trial_start:]))
        print ('offset4 start position is {}'.format(trial_start))
    offset4=offset_check(seq3[:use3],sequence,trial_start=trial_start,**vars(args)) 
    if args.verbose:
        print ('offset4 relative position is {}'.format(offset4))
    if offset4 == 'NA':
        return ('P',4,use3)
    offset4=offset4+trial_start
    if args.verbose:
        print('offset4 real position is {}'.format(offset4))
        
    return('G',(offset1,offset2,offset3,offset4))

def check_qual (seq,win=10,cutoff=25,**kwargs):
    '''scan through seq with quality with window size of win, if average quality score lower than cutoff,
identify as bad quality (return value 0, discard later), otherwise identify as good (return value 1)'''
    for i in range(len(seq)-win):
        if sum(seq.letter_annotations['phred_quality'][i:i+win])/win<cutoff:
            return 0
    return 1

@timeit_decorator
def wrapper(primers,infile1,infile2,outfile,archfile1,archfile2,**kwargs):
    inhandle1 = SeqIO.parse(infile1,'fastq')
    inhandle2 = SeqIO.parse(infile2,'fastq')
    
    readcount = 0; goodcount = 0; 
    errorcounts={'nomap':0,'ambiguous':0,'partial':0} #nocount, ambiguous count, bad quality count, partial count
    
    for record1, record2 in zip(inhandle1,inhandle2):

        assert record1.id == record2.id
        readcount += 1
        
        if readcount % 1000 == 0: 
            print('readline:'+str(readcount)+' Good count:'+str(goodcount)+' Error count:',errorcounts)
        
        if args.verbose:
            print ('Mapping {}...'.format(readcount))
        
        #if check_qual (record1) == 0 or check_qual (record2) == 0:
        #    bcount+=1
        #    continue
            
        if args.verbose:
            print ('Mapping strand: ',end='\t')
        strand=check_strand(record1.seq,record2.seq,primers,**vars(args))
        if args.verbose:
            print ('strand is {}'.format(strand))
            print ('*'*10)
        
        if strand == 'A':
            errorcounts['ambiguous']+=1
            continue
        if strand == 'N':
            errorcounts['nomap']+=1
            continue
            
        if strand == 'F':
            if args.verbose:
                print ("Mapping result1")
            result1 = mapping(primers,record1.seq,overhang='3',bc_length=21,UMI_length=15,**vars(args))
            if args.verbose:
                print ("Mapping result2")
            result2 = mapping(primers,record2.seq.reverse_complement(),overhang='5',bc_length=21,UMI_length=15,**vars(args))
            if result1[0]=='G' and result2[0]=='G':
                goodcount+=1
            else:
                errorcounts['partial']+=1
                continue
        
        if strand == 'R':
            record1,record2=record2,record1
            if args.verbose:
                print ("Mapping result1")
            result1 = mapping(primers,record1.seq,overhang='3',bc_length=21,UMI_length=15,**vars(args))
            if args.verbose:
                print ("Mapping result2")
            result2 = mapping(primers,record2.seq.reverse_complement(),overhang='5',bc_length=21,UMI_length=15,**vars(args))
            if result1[0]=='G' and result2[0]=='G':
                goodcount+=1
            else:
                errorcounts['partial']+=1
                continue
                
        #error correction
        offsets1=result1[1]
        offsets2=result2[1]
        BC1=record1[offsets1[0]:offsets1[1]]
        BC2=record2.reverse_complement()[offsets2[0]:offsets2[1]]
    
        UMI1=record1[offsets1[2]:offsets1[3]]
        UMI2=record2.reverse_complement()[offsets2[2]:offsets2[3]]
                
        BC=error_correction(BC1,BC2)
        UMI=error_correction(UMI1,UMI2)
    
        if len(BC) == 0: BC='-'
        if len(UMI) == 0: UMI='-'
        outfile.write(BC+'\t'+UMI+'\n')
    
        if args.verbose: print (BC,UMI)
    
        SeqIO.write(record1,archfile1,'fastq')
        SeqIO.write(record2,archfile2,'fastq')
        
    print('readline:'+str(readcount)+' Good count:'+str(goodcount)+' Error count:',errorcounts)
    
    return readcount,goodcount,errorcounts

def main(args):
    #change info here
    foldername=args.foldername
    primers = {'concensus5':'TACAAGAATAAGACAGGGCTTGGAAAGGGCTTTGTTGTAGG','midconcensus':'TATAAGATGGGTGGCAAGTGGtCAAAACGTAGTGTGCCTGGATGGTCTACTGTAAGGGAAAGAATG','concensus3':'CACTGGCTGACACTCATCCT'}
    workpath = '/u/scratch/y/yuanshi/{}/'.format(foldername)
    inpath = '/u/scratch/y/yuanshi/{}/'.format(foldername)
    logpath='/u/home/y/yuanshi/JCK/joblog/{}_HIV'.format(foldername) #update on 040324, put logfiles into _HIV folders, this is in compatible with cleanup2.py    
    outpath=os.path.join(workpath,'HIV_barcode')
    archpath=os.path.join(workpath,'HIV_archive')
    Path(outpath).mkdir(parents=True,exist_ok=True)
    Path(archpath).mkdir(parents=True,exist_ok=True)

    #setup flag
    win=args.win
    cutoff=args.cutoff
    short3=args.short3
    short5=args.short5
    filename=args.filename
    
    #setup filehandle #sys.argv only contain the R1 filename
    if args.binary == 'text':
        infile1 = open(os.path.join(inpath,'split',filename))
        infile2 = open(os.path.join(inpath,'split',filename.replace('R1','R2'))) #R2 file
 
    if args.binary == 'gz':
        infile1 = gzip.open(os.path.join(workpath,'split',filename),'rt')
        infile2 = gzip.open(os.path.join(workpath,'split',filename.replace('R1','R2')),'rt') #R2 file
   
    outfile = open(os.path.join(workpath,'HIV_barcode',filename.replace('R1','barcode')),'w')

    archfile1 = open(os.path.join(workpath,'HIV_archive',filename),'w')
    archfile2 = open(os.path.join(workpath,'HIV_archive',filename.replace('R1','R2')),'w')
    
    #run
    results=wrapper(primers,infile1,infile2,outfile,archfile1,archfile2,**vars(args)),
    Path(logpath).mkdir(parents=True,exist_ok=True)
    
    with open(os.path.join(logpath,filename.replace('R1','count')),'w') as f:
        json.dump(results,f)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Define the verbose flag
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode.')
    
    # Define cutoff
    parser.add_argument('-c', '--cutoff', type=int, default=3, help='Setting up cutoff for hamming distance measurement.')
    parser.add_argument('-w', '--win', type=int, default=6, help='Setting up sliding window size for offset mapping.')

    
    # Define overhang length
    parser.add_argument('-3', '--short3', type=int, default=0, help='Set up 3\' overhang.')
    parser.add_argument('-5', '--short5', type=int, default=0, help='Set up 5\' overhang.')
    parser.add_argument('-f','--filename',type=str,help='Get R1 filename.')    
    parser.add_argument('-d','--foldername',type=str,help='Get folder name in scratch.')    

    
    #add binary flag
    parser.add_argument('-b','--binary',type=str, default='gz', help='if the input fastq is in gz format')

    args = parser.parse_args()

    main(args)
