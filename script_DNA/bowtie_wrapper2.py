#!/usr/bin/env python bowtie_wrapper.py
'''changed based on Tianhao's bowtie_wrapper_R5.py'''
import os
import glob
import sys
from pathlib import Path

def main():
  workpath = '/u/scratch/y/yuanshi/NovaSeq101625/'
  homepath = '/u/home/y/yuanshi/JCK/script'
  #prepare path
  subfolders = ['ISS_bowtie','ISS_bowtie/hg38','ISS_bowtie/HIV']
  for subfolder in subfolders:
      Path(os.path.join(workpath,subfolder)).mkdir(parents=True,exist_ok=True)

  infiles = sorted(glob.glob(workpath+'ISS_genome/UDP*R1*.fastq'))
  samplelist = [x.rsplit('/')[-1].rsplit('_')[0] for x in infiles] 
  samplelist = list(set(samplelist)) #UDP indexes
  for sample in samplelist:
      os.system('cat '+os.path.join(workpath,'ISS_genome','{}*R1*fastq'.format(sample))+' > '+os.path.join(workpath,'ISS_genome','combined_{}_R1.fastq'.format(sample)))
      os.system('cat '+os.path.join(workpath,'ISS_genome','{}*R2*fastq'.format(sample))+' > '+os.path.join(workpath,'ISS_genome','combined_{}_R2.fastq'.format(sample)))
      infile1 = workpath+'ISS_genome/combined_'+sample+'_R1.fastq'
      infile2 = infile1.replace('_R1','_R2')
      for j in ['hg38','HIV']:
          outfile = os.path.join(workpath,'ISS_bowtie',j,'{}.sam'.format(sample))
          metfile = os.path.join(workpath,'ISS_log','bowtie_mapping_{}_{}.txt'.format(j,sample))
          bashfile = open(os.path.join(homepath,'tmp','bowtie_wrapper_{}_{}.sh'.format(j,sample)),'w')
          bashfile.write('#!/bin/bash\n')
          bashfile.write('source /u/local/Modules/default/init/modules.sh\n')
          bashfile.write('module load bowtie2\n')
          if j == 'hg38': reffile = '~/ref/bowtie2/hg38'
          #if j == 'HIV': reffile = '~/ISS/ref/HIV'
          if j == 'HIV': reffile = '~/ref/bowtie2/HIV_nfnsx'
          bashfile.write('bowtie2 --very-sensitive -x '+reffile+' -1 '+infile1+' -2 '+infile2+' -S '+outfile+' --met-file '+metfile)
          bashfile.close()
          os.system('chmod 777 '+os.path.join(homepath,'tmp','bowtie_wrapper_{}_{}.sh'.format(j,sample)))
          os.system('qsub -cwd -V -N PJbowtie -l h_data=6144M,h_rt=3:00:00 '+os.path.join(homepath,'tmp','bowtie_wrapper_{}_{}.sh'.format(j,sample)))


if __name__ == '__main__':
    main()
