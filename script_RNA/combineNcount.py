#!/usr/bin/python3 combineNcount.py
'''modified on 022924, now folder is taken from user, in the second argument'''
'''modified on 101824, take argument from user'''
'''modified on 042925, check git log, used for NovaSeq041725'''


import sys,os
from collections import Counter
from pathlib import Path
import time, argparse

def timeit_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper


def count(path,folder,files,outputfile):
    counter = Counter()
    for inputfile in files:
        # Open the input file
        with open(inputfile, 'r') as f:
        # Read through the file line by line
            for line in f:
            # Remove the newline character at the end of the line and split the line into barcode and UMI
                barcode, umi = line.rstrip('\n').split('\t')
            # Update the counter
                counter[(barcode, umi)] += 1

    # Print out the barcode, UMI, and count for each unique barcode-UMI combination
    with open(outputfile,'w') as f:
        for (barcode, umi), count in counter.items():
            print(f'{barcode}\t{umi}\t{count}',file=f)

@timeit_decorator
def main(args):
    path = '/u/scratch/y/yuanshi/'
    foldername = args.foldername
    project = args.project
    index = args.index
    print (index)

    #inputfiles=[os.path.join(path,foldername,'{}_barcode'.format(project),i) for i in os.listdir(os.path.join(path,foldername,'{}_barcode'.format(project))) if i.startswith(index)] #check this
    print (path)
    print (foldername)
    print (os.path.join(path,foldername,"Barcode"))
    print (os.listdir(os.path.join(path,foldername,f"{project}_barcode")))

    inputfiles=[os.path.join(path,foldername,f"{project}_barcode",i) for i in os.listdir(os.path.join(path,foldername,f"{project}_barcode")) if i.startswith(index)]
    Path(os.path.join(path,foldername,'{}_output'.format(project))).mkdir(parents=True,exist_ok=True)
    outputfile=os.path.join(path,foldername,'{}_output'.format(project),'{}_count.txt'.format(index))

    count(path,foldername,inputfiles,outputfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--foldername',type=str,help='Get folder name in scratch.')
    parser.add_argument('-p','--project',type=str,help='Get project name.')
    parser.add_argument('-i','--index',type=str,help='Get index name.')
    args = parser.parse_args()

    main(args)
