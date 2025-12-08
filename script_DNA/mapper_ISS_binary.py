#!/usr/bin/python3 mapper_ISS.py
'''written 01202024, based on Tianhao's ISS mapper and barcode-primerID mapper. map HIV integration site, do not do error correction here
extract all extractable information from R1 and R2, store in log folder as table. genome is defined as IS longer than 10nt in R2 mapping and 
no plasmid sequence mapped. original fastq files for partial file which can still somehow extract information are stored. 
effective map only when IS is potentially genome is stored as fastq, barcode and UMI are stored as its id. For plasmid and short reads,
barcode and UMI can extracted from the log txt file '''
'''modified 10022024, add a binary flag, if set as true, the input is in .gz format'''

import os, math, json, operator, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
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

def hamming(s1,s2):
    assert len(s1) == len(s2)
    s1 = s1.upper()
    s2 = s2.upper()
    return sum([0 if i==j else 1 for i,j in zip(s1,s2)])

def rc(seq):
    seq = seq.upper()
    assert set(seq).issubset(set ('ATGCN'))
    tb = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    seq = ''.join((tb[i] for i in seq))
    return seq[::-1]

def offset_check (primer, sequence,end=5, mapsize=15,**kwargs): 
    '''scanning from -win to +win, use concensus to map sequence, can define cutoff'''
    # setup parameters
    win5 = kwargs.get('win5',kwargs['win'])
    win3 = kwargs.get('win3',kwargs['win']) #if not set, scanning set to -6 to +6
    verbose = kwargs.get('verbose',False) # if not set, not verbose

    if end == 5: #use primer's 5'end for mapping, will return the 5' coordinate
        c = primer[:mapsize] if len(primer) >= mapsize else primer #use the beginning 15nt of concensus for mapping
        start = 0
    if end == 3: #use primer's 3'end for mapping, will return the 3' coordinate
        c = primer[-mapsize:] if len(primer) >= mapsize else primer #use the last 15nt of concensus for mapping
        start = len(primer)-len(c) #start >= 0

    cutoff = kwargs.get('cutoff',mapsize*0.2) #allow 20% mismatch
    if cutoff >= 0.5*mapsize: 
        cutoff = mapsize*0.4 #too flexible, change to 40% mismatch
        
    for mis in range(-win5,win3,1):#allowed scanning window -win to +win
        if start+mis+len(c)>len(sequence):
            break #cannot map
        if start+mis<0:
            target = 'N'*-(start+mis)+sequence[:start+mis+len(c)]
        else:
            target = sequence[start+mis:start+mis+len(c)] #extract target sequence for mapping
        if verbose:
            print (target,c)
        dis = hamming(target,c)
        if dis <= cutoff and dis >= 0:
            if verbose:
                print('mapped',c,target,dis)
            if end == 5:
                return (mis+start)
            if end == 3:
                return (mis+start+len(c))
    return ('NA')

def R1_map(concensus,extract,sequence,**kwargs):
    '''R1 read configuration:
       Nef: `5'-TTGACCACTTGCCACCCATCTTATA-3'` 25nt
       barcode: `5'-GNNGNNGNNGNNGNNGNNGNN-3'` 21nt
       Env+EcoRI: `5'-CCTACAACAAAGCCCTTTCCAAGCCCTGAATTC-3'` 33nt
       UMI: `5'-NNNNNNNNNNNNNN-3'`14nt
       Adatpor: `5'-CAACTGCAAACCTCAAT-3'` 17nt
       Integration Site: `5'-NNNNNNNNNNNNNN-3'` varies
       3UTR: `5'-TGCTAGAGATTTTCCACACTGACTAA-3'` 26nt''' 
    
    verbose = kwargs.get('verbose',False)
    result = {'barcode':'','UMI':'','IS':''}
    flag = [0,0,0] #1st digit barcode, 2nd digit UMI, 3rd digit IS; 0:non can be anchored, 1: 5'end anchored, 2: 3'end anchored
    # map barcode
    offset1 = offset_check(concensus['Nef'],sequence,end=3,**kwargs)
    if offset1 == 'NA':
        return 'bad',flag,result
    
    offset2 = offset_check(concensus['Env_EcoRI'],sequence[offset1:],win5=0,win3=extract['barcode']+6,end=5,**kwargs)
    if offset2 == 'NA':
        flag[0] = 1
        result['barcode'] = sequence[offset1:]
        return 'partial',flag, result #Only Nef mapped, barcode 3'end is not anchored
    offset2 += offset1
    flag[0] = 2
    result['barcode'] = sequence[offset1:offset2]
    
    # map UMI
    offset3 = offset_check(concensus['Env_EcoRI'],sequence[offset2:],end=3,**kwargs)
    if offset3 == 'NA':
        return 'partial',flag,result #Only Nef+Env5' mapped, barcode can be extracted
    offset3 += offset2
    
    offset4 = offset_check(concensus['Adaptor'],sequence[offset3:],win3=extract['UMI']+6,end=5,**kwargs)
    if offset4 == 'NA':
        flag[1] = 1
        result['UMI'] = sequence[offset3:]
        return 'partial',flag,result #Only Nef+Env mapped, UMI 3'end is not anchored
    offset4 += offset3
    flag[1] = 2
    result['UMI'] = sequence[offset3:offset4]
    
    # map IS
    offset5 = offset_check(concensus['Adaptor'],sequence[offset4:],end=3,**kwargs)
    if offset5 == 'NA':
        return 'partial',flag,result #Only Nef+Env+Adaptor mapped, barcode and UMI can be extracted
    offset5 += offset4

    offset6 = end_mapping(concensus['3UTR'],sequence,**kwargs)
    if offset6 == 'NA':
        flag[2] = 1
        result['IS'] = sequence[offset5:]
        return 'mapped',flag,result
    else:
        flag[2] = 2
        result['IS'] = sequence[offset5:offset6]
    return 'mapped',flag,result
        

    
def R2_map(concensus,extract,sequence,**kwargs):
    '''R2 read configuration:
       3UTR: `5'-TGCTAGAGATTTTCCACACTGACTAA-3'` 26nt
       Integration Site: `5'-NNNNNNNNNNNNNN-3'` varies
       Adatpor: `5'-CAACTGCAAACCTCAAT-3'` 17nt
       UMI: `5'-NNNNNNNNNNNNNN-3'`14nt
       Env+EcoRI: `5'-CCTACAACAAAGCCCTTTCCAAGCCCTGAATTC-3'` 33nt
       barcode: `5'-GNNGNNGNNGNNGNNGNNGNN-3'` 21nt
       Nef: `5'-TTGACCACTTGCCACCCATCTTATA-3'` 25nt
       '''
    #concensus = {'Nef':'TTGACCACTTGCCACCCATCTTATA','barcode':'GNNGNNGNNGNNGNNGNNGNN','Env_EcoRI':'CCTACAACAAAGCCCTTTCCAAGCCCTGAATTC','UMI':'NNNNNNNNNNNNNN','Adaptor':'CAACTGCAAACCTCAAT','3UTR':'TGCTAGAGATTTTCCACACTGACTAA'}
    #extract = {'barcode':21,'UMI':14} #
    verbose = kwargs.get('verbose',False)
    concensus_RC ={k:rc(v) for k,v in concensus.items()}
    result = {'barcode':'','UMI':'','IS':''}
    flag = [0,0,0,0] #1st digit IS and type, 3rd digit UMI, 4th digit barcode; for 2-4th place: 0 is not found, 1 is 5' end anchored, 2 is 3' end anchored; for the 1st place: 1 is plasmid, 2 is short, 0 is default(genome)
    # map IS
    offset1 = offset_check(concensus_RC['3UTR'],sequence,end=3,**kwargs)
    if verbose: print ('offset1:',offset1)
    if offset1 == 'NA':
        return 'bad',flag,result
    # scan till the end
    plasmid_offset = offset_check(concensus_RC['plasmid'],sequence[offset1:],win3=10,end=5,**kwargs)
    if plasmid_offset != 'NA':
        flag[0] = 1 #plasmid
    offset2 = offset_check(concensus_RC['Adaptor'],sequence[offset1:],win3=len(sequence)-offset1,end=5,**kwargs)
    if offset2 == 'NA':
        offset2 = end_mapping(concensus_RC['Adaptor'],sequence,**kwargs)
        if verbose: print ('end map, offset2:',offset2)
        if offset2 == 'NA':
            flag[1] = 1
            result['IS'] = sequence[offset1:]
        else:
            flag[1] = 2
            result['IS'] = sequence[offset1:offset2]
        if len(result['IS'])<=10: flag[3] = 2 #short
        return 'mapped',flag, result #only IS is retrieved, pn (plasmid, not complete)
    offset2 += offset1
    flag[1] = 2
    result['IS'] = sequence[offset1:offset2]
    if len(result['IS'])<=10: flag[0] = 2 #short
    if verbose:print ('offset2:',offset2,result)
    
    # map UMI
    offset3 = offset_check(concensus_RC['Adaptor'],sequence[offset2:],end=3,**kwargs)
    if offset3 == 'NA':
        return 'mapped',flag,result #only IS is retrieved
    offset3 += offset2
    if verbose: print ('offset3:',offset3)
    offset4 = offset_check(concensus_RC['Env_EcoRI'],sequence[offset3:],win5=0,win3=extract['UMI']+6,end=5,**kwargs)
    if offset4 == 'NA':
        offset4 = end_mapping(concensus_RC['Env_EcoRI'],sequence,**kwargs)
        if verbose: print ('endmap,offset4:',offset4)
        if offset4 == 'NA':
            flag[2] = 1
            result['UMI'] = sequence[offset3:]
        else:
            flag[2] = 2
            result['UMI'] = sequence[offset3:offset4]
        return 'mapped',flag,result #IS and UMI can be retrieved
    offset4 += offset3
    flag[2] = 2
    result['UMI'] = sequence[offset3:offset4]
    if verbose:print('offset4:',offset4,result)

    # map barcode
    offset5 = offset_check(concensus_RC['Env_EcoRI'],sequence[offset4:],end=3,**kwargs)
    if offset5 == 'NA':
        return 'mapped',flag,result #IS and UMI can be retrieved
    offset5 += offset4
    if verbose:print('offset5:',offset5,result)
    offset6 = offset_check(concensus_RC['Nef'],sequence[offset5:],win5=0,win3=extract['barcode']+6,end=5,**kwargs)
    if verbose: print ('offset6:',offset6,offset5)
    if offset6 == 'NA':
        offset6 = end_mapping(concensus_RC['Nef'],sequence,**kwargs)
        if verbose: print ('end map,offset6:',offset6)
        if offset6 == 'NA':
            flag[3] = 1
            result['barcode'] = sequence[offset5:]
            return 'mapped',flag,result #IS, UMI and barcode can be retrieved
        else:
            flag[3] = 2
            result['barcode'] = sequence[offset5:offset6]
    else:
        flag[3] = 2
        offset6+=offset5
        result['barcode'] = sequence[offset5:offset6]
    return 'mapped',flag,result

def end_mapping(primer,sequence,**kwargs):
    for i in range(len(primer),0,-1):
        dis = hamming(sequence[-i:],primer[:i])
        if dis <= 0.2*i:
            return -i
    return 'NA'

@timeit_decorator
def map_record(infile1, infile2, concensus,extract,out_handlers,**kwargs):
    inhandle1 = SeqIO.parse(infile1,'fastq')
    inhandle2 = SeqIO.parse(infile2,'fastq')

    #good: mapped result good from both R1 and R2, extract barcode UMI from R1, IS from R1 and R2, save fastq file as well as the extracted genome file
    #partial: only R1 or R2 can be used, but can extract barcode UMI and IS from that read, save fastq only for later examination
    #problem: cannot extract all information, but not all bad, save fastq only for later examination
    #no map: cannot map in either R1 or R2, dump, only count
    counts = {'total':0,'genome':0,'nomap':0,'problem':0,'partial':0,'short':0,'plasmid':0} #nocount, ambiguous count, bad quality count, partial count    

    print('\t'.join(['map_R1','R1_flag','map_R2','R2_flag','R1_barcode','R2_barcode','R1_UMI','R2_UMI','R1_IS','R2_IS']),file=out_handlers['log'])

    for record1,record2 in zip(inhandle1,inhandle2):
        assert record1.id == record2.id
        counts['total'] += 1
        if counts['total'] % 1000 == 0:
            for k,v in counts.items():
                print ('{}:{}'.format(k,v),end='\t')
            print('\n')

        check1,flag1,result1 = R1_map(concensus,extract,record1,**kwargs)
        flag1 = ''.join(str(i) for i in flag1)
        check2,flag2,result2 = R2_map(concensus,extract,record2,**kwargs)
        flag2 = ''.join(str(i) for i in flag2)    
    
        print('\t'.join([check1,flag1,check2,flag2]),end='\t',file=out_handlers['log'])
        for (v1),(v2) in zip(result1.values(),result2.values()):
            s1=str(v1.seq) if v1 else '-'
            s2=str(v2.seq.reverse_complement()) if v2 else '-'
            print ('\t'.join([s1,s2]),end='\t',file=out_handlers['log'])
        print(file=out_handlers['log'])

        if check1 == 'bad' and check2 == 'bad':
            counts['nomap'] += 1
            continue

        if check1 != 'mapped' and (flag2.endswith('222') or flag2.endswith('221')) : #everything extract from R2
            counts['partial'] += 1
            out_handlers['problem_R1'].write(record1.format('fastq'))
            out_handlers['problem_R2'].write(record2.format('fastq'))
            continue

        if check2 != 'mapped' and (flag1=='221' or flag1=='222'): #everything extract from R1
            counts['partial'] += 1
            out_handlers['problem_R1'].write(record1.format('fastq'))
            out_handlers['problem_R2'].write(record2.format('fastq'))
            continue

        if check1 == 'mapped' and flag2.startswith('2'): #short
            counts['short'] += 1
            continue

        if check1 == 'mapped' and flag2.startswith('1'): #plasmid
            counts['plasmid']+=1
            continue
    
        if flag1.startswith('22') and (flag2.startswith('02') or flag2.startswith('01')): #good
            counts['genome']+=1
            barcode = str(result1['barcode'].seq)
            UMI = str(result1['UMI'].seq)
            r1 = result1['IS']
            r2 = result2['IS']
            if not r1:
                r1 = SeqRecord(Seq(""),id=record1.id,letter_annotations={"phred_quality": []})
            if not r2:
                r2 = SeqRecord(Seq(""),id=record2.id,letter_annotations={"phred_quality": []})
            r1.id += (':'+UMI+':'+barcode)
            r2.id += (':'+UMI+':'+barcode)
            r1.description = ''; r2.description = ''
            SeqIO.write(r1, out_handlers['genome_R1'], "fastq")
            SeqIO.write(r2, out_handlers['genome_R2'], "fastq")

    return counts


def main(args):
    #change info here
    concensus = {'Nef':'TTGACCACTTGCCACCCATCTTATA','barcode':'GNNGNNGNNGNNGNNGNNGNN','Env_EcoRI':'CCTACAACAAAGCCCTTTCCAAGCCCTGAATTC',
             'UMI':'NNNNNNNNNNNNNN','Adaptor':'CAACTGCAAACCTCAAT','3UTR':'TGCTAGAGATTTTCCACACTGACTAA',
             'plasmid':'TTGGCTCACTGCAACCTCTACCTCCTGGG'}
    extract = {'barcode':21,'UMI':14} #

    filename=args.filename
    
    foldername=args.foldername
    workpath = '/u/scratch/y/yuanshi/{}/'.format(foldername)
    logpath='/u/home/y/yuanshi/JCK/joblog/{}_ISS'.format(foldername)
    Path(logpath).mkdir(parents=True,exist_ok=True)

    subfolders = ['ISS_genome','ISS_log','ISS_problem']
    
    for subfolder in subfolders:
        Path(os.path.join(workpath,subfolder)).mkdir(parents=True, exist_ok=True)
    infile1 = filename
    infile2 = infile1.replace('R1','R2')

    name = infile1.replace('R1_','')
    logfile = open(os.path.join(workpath,'ISS_log','{}_log'.format(name)),'w')
    p1 = open(os.path.join(workpath,'ISS_problem','{}_R1.fastq'.format(name)),'w')
    p2 = open(os.path.join(workpath,'ISS_problem','{}_R2.fastq'.format(name)),'w')
    g1 = open(os.path.join(workpath,'ISS_genome','{}_R1.fastq'.format(name)),'w')
    g2 = open(os.path.join(workpath,'ISS_genome','{}_R2.fastq'.format(name)),'w')

    out_handlers={'log':logfile,'problem_R1':p1,'problem_R2':p2,'genome_R1':g1,'genome_R2':g2}


    #setup filehandle #sys.argv only contain the R1 filename

    if args.binary == 'text':
        infile1 = open(os.path.join(workpath,'split',filename))
        infile2 = open(os.path.join(workpath,'split',filename.replace('R1','R2'))) #R2 file

    if args.binary == 'gz':
        infile1 = gzip.open(os.path.join(workpath,'split',filename),'rt')
        infile2 = gzip.open(os.path.join(workpath,'split',filename.replace('R1','R2')),'rt') #R2 file

    #run
    results=map_record(infile1, infile2, concensus,extract,out_handlers,**vars(args))

    with open(os.path.join(logpath,filename.replace('R1','count')),'w') as f:
        json.dump(results,f)

    for hd in out_handlers.values():
        hd.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Define the verbose flag
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose mode.')

    # Define cutoff
    parser.add_argument('-c', '--cutoff', type=int, default=3, help='Setting up cutoff for hamming distance measurement.')
    parser.add_argument('-5', '--win', type=int, default=6, help='Setting up sliding window size for offset mapping.')

    parser.add_argument('-f','--filename',type=str,help='Get R1 filename.')
    parser.add_argument('-d','--foldername',type=str,help='Get folder name in scratch.')

    # add binary flag
    parser.add_argument('-b','--binary',type=str, default='gz', help='if the input fastq is in gz format')

    args = parser.parse_args()

    main(args)
