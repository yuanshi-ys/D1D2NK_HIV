#! /bin/python3
'''changed on 4/28/2025, all .fastq.gz files are at the same folder, also the library names is no longer UDP***, but the indicated sample names
    perform one library split for each task'''
'''updated again, to allow system IO sync after splitting, also ensure gzip finish before the job close'''
'''updated 05/09/2025, remove the gzip step that causes trouble all the time, do altogether afterwards, requires a large buffer space'''

from pathlib import Path
import glob, re, os, subprocess, time
from subprocess import PIPE
import sys,json
import argparse
import logging
from datetime import datetime

def timeit_decorator(func):
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} executed in {end_time - start_time} seconds")
        return result
    return wrapper

def exe(args):
    result = subprocess.run(args, shell=True, stdout=PIPE, stderr=PIPE)
    if result.returncode != 0:
        print(f"Command failed with exit code {result.returncode}")
        print(f"Stdout: {result.stdout.decode('utf-8')}")
        print(f"Stderr: {result.stderr.decode('utf-8')}")
        raise RuntimeError(f"Command failed: {args}")
    return result.stdout.decode('utf-8')

def gen_args(filenames, outpath):
    for n in filenames:
        index = os.path.basename(n).split('_')[0]
        R = os.path.basename(n).split('_')[-2]
        outname = os.path.join(outpath, f'{index}_{R}_')
        yield 'zcat "{n}" | split -l 4000000 - "{outname}"'.format(n=n, outname=outname)

def setup_logger(name, log_dir='logs'):
    # Create logs directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)

    # Generate unique log filename
    job_id = os.getenv('JOB_ID', 'local')
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = f'{log_dir}/{name}_{job_id}_{timestamp}.log'

    # Configure logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)

    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.propagate = False
    return logger

@timeit_decorator
def main(flags):
    logger = setup_logger('splitjob')
    # load indexfile
    indexfile = os.path.join('/u/home/y/yuanshi/JCK/Fasta',flags.indexfile)
    if not os.path.exists(indexfile):
        raise FileNotFoundError(f"Index file not found: {indexfile}")

    if indexfile.endswith('.json'):
        with open(indexfile) as f:
            UDPs = json.load(f)
    else:
        UDPs = open(indexfile, 'r').read().splitlines()
    if flags.verbose:
        print (f'Index file used {flags.indexfile}, indexes are:', *UDPs)
    logger.info(f'Index file: {flags.indexfile}')

    # load folder
    path = "/u/scratch/y/yuanshi/{}/data/".format(flags.foldername)
    outpath = "/u/scratch/y/yuanshi/{}/split".format(flags.foldername)
    Path(outpath).mkdir(parents=True, exist_ok=True)

    task_id = os.environ.get('SGE_TASK_ID','1')
    
    token = UDPs[int(task_id)-1]

    logger.info(f'Index: {token}')
    search_path = "{}{}{}*.gz".format(path,"*/"*flags.subfolders,token) #allow any among of subfolders
    inputnames = glob.glob(search_path) 
    if flags.verbose:
        print ('Current index is {}, invloving these files {}'.format(token,'\n'.join([os.path.basename(name) for name in inputnames])))
        # name is like this UDP0100_S100_L006_R1_001.fastq.gz
    args = gen_args(inputnames, outpath)

    for arg in args:
        if flags.verbose:
            print (f'executables is {arg}')
        logger.info(f'{arg}')
        exe(arg)

    logger.info('Finished splitting. Forcing file system sync..')
    os.sync()
    time.sleep(2)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # Define the verbose flag
    parser.add_argument('-v', '--verbose', action='store_true', default=False,help='Enable verbose mode.')
    # Define
    parser.add_argument('-d','--foldername',type=str,default='NovaSeq121424',help='Get folder name in scratch.')
    parser.add_argument('-s','--subfolders',type=int,default=0)
    parser.add_argument('-i','--indexfile',type=str,default='UDPs.txt',help='Get index file name, default is UDP.txt with all 96 UDPs.')
    parser.add_argument('-b','--binary', action='store_true',default=True,help='Set binary flag, if gzip type is needed')
    parser.add_argument('--no-binary', action='store_false', dest='binary',help='Disable gzip compression')
    flags = parser.parse_args()

    main(flags)
