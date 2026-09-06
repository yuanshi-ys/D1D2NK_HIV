#!/bin/bash
#$ -cwd
#$ -j y
#$ -N LINKAGE_GSE335409
#$ -l h_data=4G,h_rt=1:00:00
#$ -t 1-117

echo $SGE_TASK_ID
python3 /u/home/y/yuanshi/JCK/script/readsam2.py `sed -n ${SGE_TASK_ID}p /u/home/y/yuanshi/JCK/Fasta/GSE335409_ISS.txt` GSE335409
