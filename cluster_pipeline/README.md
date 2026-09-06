Cluster-side pipeline that produces `ISS_linkage/linkage_UDP*.txt` (the input to `../get_perIS_lite.py`). Runs on Hoffman2 (SGE), not locally — paths are hardcoded to `/u/home/y/yuanshi/JCK` and `/u/scratch/y/yuanshi/GSE335409`.

Run order:

1. `split_gzip_interleave_GSE335409.sh` / `_others.sh` (SGE array jobs) -> `split.py` -- split raw per-lane `fastq.gz` into per-UDP chunks under `GSE335409/split/`.
2. `gzip_split_GSE335409_first.sh` -- gzip any leftover unzipped split chunks.
3. `map_ISS_wrapper_GSE335409.sh` (+ `_retry`, `_retry2`, `_retry3`) -> `mapper_ISS_binary.py` -- extract barcode/UMI/integration-site from each read pair, write genome-mapped fastq to `ISS_genome/`.
4. `bowtie_wrapper_GSE335409.py` -- concatenate per-chunk fastqs by UDP, align with `bowtie2 --very-sensitive` against hg38 and HIV_nfnsx, write `ISS_bowtie/{hg38,HIV}/UDP*.sam`.
5. `linkage_wrapper_GSE335409.sh` (SGE array, 1-117) -> `readsam2.py` -- parse both SAM files per UDP, majority-vote barcode/site per UMI, write `ISS_linkage/linkage_UDP*.txt`.

Source: `~/JCK/script/` and `~/JCK/Fasta/` on Hoffman2.
