#!/bin/bash
#$ -cwd
#$ -j y
#$ -N GZIP1_GSE335409
#$ -l h_data=4G,h_rt=12:00:00

INDEXFILE="/u/home/y/yuanshi/JCK/Fasta/GSE335409_gzip_first.json"
TOKEN=$(python3 -c "import sys, json; print(json.load(open(\"$INDEXFILE\"))[${SGE_TASK_ID}-1])")

echo "Task ${SGE_TASK_ID}: Gzipping ${TOKEN}..."
for f in /u/scratch/y/yuanshi/GSE335409/split/${TOKEN}_*; do
    if [[ -f "$f" ]] && [[ "$f" != *.gz ]]; then
        gzip "$f"
    fi
done
echo "Task ${SGE_TASK_ID}: Done ${TOKEN}"
