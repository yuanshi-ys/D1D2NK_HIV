# D1D2NK_HIV
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22476665.svg)](https://doi.org/10.5281/zenodo.22476665)

Code to analyze DNA barcode and integration sites for the manuscript "NK Cells Engineered with a Chimeric Antigen Receptor Delay HIV Rebound and Reshape HIV Reservoir Composition".

Two stages:

- `cluster_pipeline/` -- runs on Hoffman2 (SGE). Raw paired-end FASTQ (D1D2_BIseq) -> barcode/UMI/integration-site extraction -> bowtie2 alignment to hg38 + HIV -> `linkage_UDP*.txt` per sample. See `cluster_pipeline/README.md`.
- `source_data_pipeline/` -- runs locally. `get_perIS_lite.py` turns the linkage files into a per-integration-site table; `map_to_database.py` annotates it against gene/regulatory/histone-mark databases (bedtools); `generate_source_data.py` produces the Source Data Excel files for Fig. 5, 6, S7, S8.
