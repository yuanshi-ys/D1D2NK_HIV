#!/bin/bash
#$ -cwd
#$ -j y
#$ -N SPLITGZIP_GSE335409_B
#$ -l h_data=4G,h_rt=6:00:00

INDEXFILE="/u/home/y/yuanshi/JCK/Fasta/GSE335409_retry_others.json"
TOKEN=$(python3 -c "import sys, json; print(json.load(open(\"$INDEXFILE\"))[${SGE_TASK_ID}-1])")

echo "Task ${SGE_TASK_ID}: splitting ${TOKEN}..."
python3 /u/home/y/yuanshi/JCK/script/split.py -s 0 -d GSE335409 -i GSE335409_retry_others.json -v

echo "Task ${SGE_TASK_ID}: gzipping ${TOKEN}..."
for f in /u/scratch/y/yuanshi/GSE335409/split/${TOKEN}_*; do
    if [[ -f "$f" ]] && [[ "$f" != *.gz ]]; then
        gzip "$f"
    fi
done
echo "Task ${SGE_TASK_ID}: done ${TOKEN}"
