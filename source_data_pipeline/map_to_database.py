import os
import re
import subprocess
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.dirname(SCRIPT_DIR)

DB_DIR = os.path.join(UPLOAD_DIR, "database")
BED_DIR = os.path.join(UPLOAD_DIR, "bedfiles")
OUT_DIR = os.path.join(UPLOAD_DIR, "output")

PERIS_LITE_PATH = os.path.join(OUT_DIR, "perIS_lite.csv")
INTERSECT_OUT = os.path.join(OUT_DIR, "perIS_intersect.csv")
DISTANCE_OUT = os.path.join(OUT_DIR, "perIS_distance.csv")

GENE_GFF = os.path.join(DB_DIR, "hg38_Ensembl_genes.gff3")
REGU_GFF = os.path.join(DB_DIR, "hg38_Ensembl_release108_Regulatory_activity_CD4_ab.gff")
CHROM_SIZES = os.path.join(DB_DIR, "hg38.chrom.sizes")
CHROM_SIZES_CHR = os.path.join(DB_DIR, "hg38.chrom.sizes.chr")
GENE_SORTED_BED = os.path.join(DB_DIR, "hg38_Ensembl_genes_sorted.bed")

GFF3_COL = ["chrom2", "source", "type", "start", "end", "score", "strand2", "phase", "attributes"]
HISTONE_MARKS = ["H3K36me3", "H3K27ac", "H3K4me1", "H3K9me3", "H3K27me3", "H3K4me3"]
REGULIST = ["CTCF_binding_site", "TF_binding_site", "enhancer", "open_chromatin_region", "promoter"]

RANK_MAP = {
    "protein_coding": 1,
    "TR_V_gene": 2, "TR_C_gene": 2, "IG_V_gene": 2,
    "lncRNA": 3, "ncRNA_gene": 3, "snRNA": 3, "snoRNA": 3,
    "miRNA": 4, "misc_RNA": 4, "rRNA": 4,
    "transcribed_processed_pseudogene": 5, "transcribed_unprocessed_pseudogene": 5,
    "transcribed_unitary_pseudogene": 5,
    "processed_pseudogene": 6, "unprocessed_pseudogene": 6, "unitary_pseudogene": 6,
    "rRNA_pseudogene": 6,
    "TEC": 7, "pseudogene": 7, "": 7,
}


def load_perIS():
    perIS = pd.read_csv(PERIS_LITE_PATH)
    perIS = perIS.rename(columns={"AnimalKey": "AnimalID"})
    perIS["SampleID"] = perIS["AnimalID"] + "-" + perIS["Organ"]
    perIS["Mouse_Cohort"] = perIS["AnimalID"]
    perIS["MouseID"] = "placeholder"
    return perIS


def prepare_regulatory_files():
    missing = [e for e in REGULIST if not os.path.exists(os.path.join(DB_DIR, f"regu_{e}.gff3"))]
    if missing:
        reguref = pd.read_csv(REGU_GFF, sep="\t", comment="#",
                               names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
        reguref = reguref[reguref["seqid"].isin([str(x) for x in range(1, 23)] + ["X", "Y"])]
        for element in REGULIST:
            tdata = reguref[reguref["type"] == element]
            tdata.to_csv(os.path.join(DB_DIR, f"regu_{element}.gff3"), index=False, header=False, sep="\t")


def parse_attr(attribute):
    search = re.search(r";biotype=(.*?);", attribute)
    biotype = search.group(1) if search else ""
    search = re.search(r";gene_id=(.*?);", attribute)
    ensg = search.group(1) if search else ""
    search = re.search(r";Name=(.*?);", attribute)
    gene = search.group(1) if search else ensg
    return biotype, ensg, gene


def write_bed(df, perIS_col, out_path):
    with open(out_path, "w") as f:
        for _, row in df.iterrows():
            site = row["site"]
            if "chr" in site:
                record = site.rsplit(":")
                cols = [record[0][3:], record[1], record[1], record[2]]
                for field in perIS_col:
                    cols.append(str(row[field]))
                f.write("\t".join(cols) + "\n")


def annotate(perIS):
    perIS = perIS.reset_index(drop=True).copy()
    perIS_col = list(perIS.columns)

    bed_path = os.path.join(BED_DIR, "perIS.bed")
    write_bed(perIS, perIS_col, bed_path)

    genes_bed = os.path.join(BED_DIR, "perIS_genes.bed")
    subprocess.run(f'bedtools intersect -wb -a "{bed_path}" -b "{GENE_GFF}" > "{genes_bed}"', shell=True, check=True)

    bed_cols = ["chrom", "pos", "pos1", "strand"] + perIS_col
    interdata = pd.read_csv(genes_bed, sep="\t", names=bed_cols + GFF3_COL)
    interdata["Gene_Name"] = interdata["attributes"].apply(lambda x: parse_attr(x)[2])
    interdata["ENSG"] = interdata["attributes"].apply(lambda x: parse_attr(x)[1])
    interdata["Gene_Type"] = interdata["attributes"].apply(lambda x: parse_attr(x)[0])

    multiple = interdata.groupby(bed_cols).count()[["chrom2"]]
    multiple.columns = ["multiple genes"]
    multiple = multiple.reset_index()
    interdata = pd.merge(interdata, multiple, on=bed_cols)

    interdata["rank"] = interdata["Gene_Type"].map(RANK_MAP)
    interdata = interdata.sort_values("rank").reset_index(drop=True)
    interdata = interdata.drop_duplicates(perIS_col)

    fields = ["type", "Gene_Name", "ENSG", "Gene_Type", "multiple genes"]
    interdata = interdata[perIS_col + fields]

    perIS = pd.merge(interdata, perIS, on=perIS_col, how="outer")
    perIS["In_Gene"] = perIS["type"].isin(["gene"])
    perIS = perIS.drop(columns=["type"])

    prepare_regulatory_files()
    for element in REGULIST:
        regu_bed = os.path.join(DB_DIR, f"regu_{element}.gff3")
        out_bed = os.path.join(BED_DIR, f"perIS_{element}.bed")
        subprocess.run(f'bedtools intersect -a "{bed_path}" -b "{regu_bed}" > "{out_bed}"', shell=True, check=True)
        interdata = pd.read_csv(out_bed, names=bed_cols, sep="\t")[perIS_col]
        interdata = interdata.drop_duplicates()
        interdata["type"] = element
        perIS = pd.merge(perIS, interdata, on=perIS_col, how="outer")
        perIS[f"In_{element}"] = ~perIS["type"].isna()
        perIS = perIS.drop(columns="type")

    return perIS


def sort_bed(in_path, out_path):
    subprocess.run(f'bedtools sort -i "{in_path}" -g "{CHROM_SIZES}" > "{out_path}"', shell=True, check=True)


def _load_valid_chroms(path):
    chroms = set()
    with open(path) as f:
        for line in f:
            chroms.add(line.strip().split("\t")[0])
    return chroms


def add_gene_distance(perIS):
    perIS_col = list(perIS.columns)
    bed_cols = ["chrom", "pos", "pos1", "strand"] + perIS_col
    valid_chroms = _load_valid_chroms(CHROM_SIZES)
    perIS_for_bed = perIS[perIS["site"].apply(
        lambda s: s.startswith("chr") and s.split(":")[0][3:] in valid_chroms
    )].reset_index(drop=True)
    raw_bed = os.path.join(BED_DIR, "perIS_distance.bed")
    write_bed(perIS_for_bed, perIS_col, raw_bed)
    sorted_bed = os.path.join(BED_DIR, "perIS_distance_sorted.bed")
    sort_bed(raw_bed, sorted_bed)

    out_bed = os.path.join(BED_DIR, "perIS_distance_to_gene.bed")
    subprocess.run(
        f'bedtools closest -d -a "{sorted_bed}" -b "{GENE_SORTED_BED}" -g "{CHROM_SIZES}" > "{out_bed}"',
        shell=True, check=True,
    )
    sorted_bed_col = ["chrom2", "pos2", "pos21", "attributes", "score2", "strand2", "distance"]
    closedata = pd.read_csv(out_bed, sep="\t", names=bed_cols + sorted_bed_col)
    closedata["Gene_Name"] = closedata["attributes"].apply(lambda x: parse_attr(x)[2])
    closedata["ENSG"] = closedata["attributes"].apply(lambda x: parse_attr(x)[1])
    closedata["Gene_Type"] = closedata["attributes"].apply(lambda x: parse_attr(x)[0])

    multiple = closedata.groupby(bed_cols).size().reset_index(name="multiple genes")
    closedata = pd.merge(multiple, closedata, on=bed_cols)
    closedata["rank"] = closedata["Gene_Type"].map(RANK_MAP)
    closedata = closedata.sort_values("rank").reset_index(drop=True)
    closedata = closedata.drop_duplicates(perIS_col)

    perIS = pd.merge(perIS, closedata, on=perIS_col, how="outer")
    tlist = perIS_col + ["multiple genes", "Gene_Name", "ENSG", "Gene_Type", "distance"]
    perIS = perIS[tlist].rename(columns={"distance": "distance_gene"})
    return perIS, sorted_bed


def add_regulatory_distance(perIS, sorted_bed):
    perIS_col_base = [c for c in perIS.columns if c not in (
        "multiple genes", "Gene_Name", "ENSG", "Gene_Type", "distance_gene")]
    for element in REGULIST:
        regu_gff = os.path.join(DB_DIR, f"regu_{element}.gff3")
        regu_sorted = os.path.join(DB_DIR, f"regu_{element}.sort.gff3")
        if not os.path.exists(regu_sorted):
            sort_bed(regu_gff, regu_sorted)
        out_bed = os.path.join(BED_DIR, f"perIS_distance_{element}.bed")
        subprocess.run(
            f'bedtools closest -D b -a "{sorted_bed}" -b "{regu_sorted}" > "{out_bed}"',
            shell=True, check=True,
        )
        names = ["chrom", "pos", "pos1", "strand1"] + perIS_col_base + \
                ["chrom2", "source", "type", "start", "end", "score", "strand2", "phase", "attributes",
                 f"distance_{element}"]
        df = pd.read_csv(out_bed, sep="\t", names=names)
        df = df.groupby(perIS_col_base)[[f"distance_{element}"]].first().reset_index()
        perIS = pd.merge(perIS, df, on=perIS_col_base, how="outer")
    return perIS


def _make_chr_prefixed_bed(sorted_bed):
    valid_chr = _load_valid_chroms(CHROM_SIZES_CHR)
    unsorted_chr_bed = os.path.join(BED_DIR, "perIS_distance_unsorted_chr.bed")
    with open(sorted_bed) as f_in, open(unsorted_chr_bed, "w") as f_out:
        for line in f_in:
            chrom = line.split("\t", 1)[0]
            chrom_chr = chrom if chrom.startswith("chr") else "chr" + chrom
            if chrom_chr in valid_chr:
                f_out.write((chrom_chr + line[len(chrom):]) if not line.startswith("chr") else line)
    chr_bed = os.path.join(BED_DIR, "perIS_distance_sorted_chr.bed")
    subprocess.run(f'bedtools sort -i "{unsorted_chr_bed}" -g "{CHROM_SIZES_CHR}" > "{chr_bed}"', shell=True, check=True)
    return chr_bed


def add_histone_distance(perIS):
    perIS_col_base = list(perIS.columns)
    chip_col = ["chr", "start", "end", "peak_name", "peak_score", "strand2",
                "signalValue", "pValue", "qValue", "summit"]

    valid_chroms = _load_valid_chroms(CHROM_SIZES)
    perIS_for_bed = perIS[perIS["site"].apply(
        lambda s: s.startswith("chr") and s.split(":")[0][3:] in valid_chroms
    )].reset_index(drop=True)
    raw_bed = os.path.join(BED_DIR, "perIS_distance_histone_query.bed")
    write_bed(perIS_for_bed, perIS_col_base, raw_bed)
    no_prefix_sorted = os.path.join(BED_DIR, "perIS_distance_histone_query_sorted.bed")
    sort_bed(raw_bed, no_prefix_sorted)
    chr_sorted_bed = _make_chr_prefixed_bed(no_prefix_sorted)

    for histone in HISTONE_MARKS:
        histone_bed = os.path.join(DB_DIR, f"Tcell_{histone}.sort.bed")
        out_bed = os.path.join(BED_DIR, f"IS_distance_{histone}.bed")
        subprocess.run(
            f'bedtools closest -D a -a "{chr_sorted_bed}" -b "{histone_bed}" -g "{CHROM_SIZES_CHR}" > "{out_bed}"',
            shell=True, check=True,
        )
        names = ["chrom", "pos", "pos1", "strand"] + perIS_col_base + chip_col + ["distance"]
        interdata = pd.read_csv(out_bed, sep="\t", names=names)
        interdata = interdata.groupby(perIS_col_base)[["distance"]].first().reset_index()
        interdata = interdata.rename(columns={"distance": f"distance_{histone}"})
        for col in perIS_col_base:
            interdata[col] = interdata[col].astype(str)
        merge_perIS = perIS.copy()
        for col in perIS_col_base:
            merge_perIS[col] = merge_perIS[col].astype(str)
        perIS = pd.merge(merge_perIS, interdata, on=perIS_col_base, how="left")
    return perIS


def compute_perIS_distance(perIS_original_cols_only):
    perIS, sorted_bed = add_gene_distance(perIS_original_cols_only)
    perIS = add_regulatory_distance(perIS, sorted_bed)
    perIS = add_histone_distance(perIS)
    return perIS


def main():
    os.makedirs(BED_DIR, exist_ok=True)
    os.makedirs(OUT_DIR, exist_ok=True)

    perIS = load_perIS()

    annotated = annotate(perIS)
    annotated.to_csv(INTERSECT_OUT, index=False)

    distance = compute_perIS_distance(perIS)
    distance.to_csv(DISTANCE_OUT, index=False)


if __name__ == "__main__":
    main()
