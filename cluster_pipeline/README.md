Cluster-side pipeline that produces `ISS_linkage/linkage_UDP*.txt` (the input to `../get_perIS_lite.py`). Runs on Hoffman2 (SGE), not locally -- paths are hardcoded to `/u/home/y/yuanshi/JCK` and `/u/scratch/y/yuanshi/D1D2_BIseq`. Thin SGE array-job launcher scripts that just qsub one of these five scripts per sample/UDP are omitted here.

Run order:

1. `split.py` -- split raw per-lane `fastq.gz` into per-UDP chunks under `D1D2_BIseq/split/`.
2. `gzip_split_D1D2_BIseq_first.sh` -- gzip any leftover unzipped split chunks.
3. `mapper_ISS_binary.py` -- extract barcode/UMI/integration-site from each read pair, write genome-mapped fastq to `ISS_genome/`.
4. `bowtie_wrapper_D1D2_BIseq.py` -- concatenate per-chunk fastqs by UDP, align with `bowtie2 --very-sensitive` against hg38 and HIV_nfnsx, write `ISS_bowtie/{hg38,HIV}/UDP*.sam`.
5. `readsam2.py` -- parse both SAM files per UDP, majority-vote barcode/site per UMI, write `ISS_linkage/linkage_UDP*.txt`.

Source: `~/JCK/script/` on Hoffman2.
