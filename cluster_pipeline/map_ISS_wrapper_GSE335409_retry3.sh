#!/bin/bash
#$ -cwd
#$ -j y
#$ -N MAP_GSE335409_RETRY3
#$ -l h_data=4G,h_rt=12:00:00

echo $SGE_TASK_ID
python3 /u/home/y/yuanshi/JCK/script/mapper_ISS_binary.py -f `sed -n ${SGE_TASK_ID}p /u/home/y/yuanshi/JCK/Fasta/Filename_GSE335409_ISS_retry_final2.txt` -d GSE335409 -b gz
