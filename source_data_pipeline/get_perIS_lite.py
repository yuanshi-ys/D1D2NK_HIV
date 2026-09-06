import os
import pandas as pd
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.dirname(SCRIPT_DIR)

ISS_LINKAGE_DIR = "/Users/yuanshi/Library/CloudStorage/Box-Box/Sequence Analysis/2026/D1D2 paper consistency check/ISS_linkage"
GSE_INDEX_PATH = os.path.join(UPLOAD_DIR, "input", "ref.index")
RNA_CSV_PATH = os.path.join(UPLOAD_DIR, "input", "RNA.csv")
OUT_TABLE = os.path.join(UPLOAD_DIR, "output", "perIS_lite.csv")

LINKAGE_COLUMNS = [
    "UMI", "read count", "barcode", "barcode_confidence", "2nd_barcode_freq",
    "provirus type", "site", "site_confidence", "2nd_site_freq", "mapping length",
]


def get_rc(seq):
    seq = seq.upper()
    matching = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return "".join(matching.get(b, "N") for b in seq)[::-1]


def load_gse_index():
    df = pd.read_csv(GSE_INDEX_PATH, sep="\t")
    df["Mouse"] = df["Mouse"].astype(int)
    df["Treatment"] = df["Condition"].map({"D1D2": "D1D2-NK", "GFP": "GFP-NK", "No": "No NK"})
    df["AnimalKey"] = df["Condition"] + "_mouse_" + df["Mouse"].astype(str)
    return df


def load_linkage(udp):
    path = os.path.join(ISS_LINKAGE_DIR, f"linkage_{udp}.txt")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=LINKAGE_COLUMNS)
    return pd.read_csv(path, sep="\t", names=LINKAGE_COLUMNS, dtype={"UMI": str})


def qc_filter(tdata, threshold=5):
    return tdata[tdata["read count"] > threshold]


def build_per_site(gse_index):
    fields = ["UDP", "AnimalKey", "Organ", "Treatment"]
    per_site_frames = []
    for _, entry in gse_index.iterrows():
        raw = load_linkage(entry["UDP"])
        if raw.empty:
            continue
        filtered = qc_filter(raw, threshold=5)
        filtered = filtered[filtered["provirus type"] == "hg38"]
        if filtered.empty:
            continue
        tdata = filtered.groupby(["barcode", "site"]).size().reset_index(name="UMI")
        for f in fields:
            tdata[f] = entry[f]
        per_site_frames.append(tdata)
    return pd.concat(per_site_frames, ignore_index=True) if per_site_frames else pd.DataFrame()


def assign_pair_to_dominant_animal(per_is):
    per_is = per_is[~per_is["barcode"].str.contains("N")].copy()
    per_is["barcode"] = per_is["barcode"].apply(get_rc)

    pair_animal_umi = per_is.groupby(["site", "barcode", "AnimalKey"])["UMI"].sum().reset_index()
    dominant = pair_animal_umi.loc[pair_animal_umi.groupby(["site", "barcode"])["UMI"].idxmax()]
    dominant_animal = dominant.set_index(["site", "barcode"])["AnimalKey"].to_dict()

    mask = per_is.apply(lambda r: dominant_animal.get((r["site"], r["barcode"])) == r["AnimalKey"], axis=1)
    return per_is[mask].reset_index(drop=True)


def add_proliferation(per_is):
    prolif_dict = per_is.groupby(["site", "AnimalKey"]).size().to_dict()
    per_is = per_is.copy()
    per_is["Proliferation"] = per_is.apply(lambda x: prolif_dict[(x["site"], x["AnimalKey"])] > 1, axis=1)
    per_is["Proliferation_tissue"] = per_is["UMI"] > 1
    return per_is


def add_rna_columns(per_is, rna_csv_path=RNA_CSV_PATH, gse_mouse_col="AnimalKey"):
    rna = pd.read_csv(rna_csv_path)

    rna["freq_organ"] = rna.groupby(["MouseID_D1D2_BIseq", "Organ"])["UMI"].transform(lambda s: s / s.sum())
    organ_lookup = rna.set_index(["MouseID_D1D2_BIseq", "Organ", "RNA_barcode"])["freq_organ"].to_dict()

    non_tb = rna[rna["Organ"] != "TB"].copy()
    non_tb["freq_organ"] = non_tb.groupby(["MouseID_D1D2_BIseq", "Organ"])["UMI"].transform(lambda s: s / s.sum())
    animal_lookup = {}
    for mouse, g in non_tb.groupby("MouseID_D1D2_BIseq"):
        pivot = g.pivot_table(index="RNA_barcode", columns="Organ", values="freq_organ", fill_value=0)
        animal_lookup[mouse] = pivot.mean(axis=1).to_dict()

    tb = rna[rna["Organ"] == "TB"].copy()
    tb["freq_tb"] = tb.groupby("MouseID_D1D2_BIseq")["UMI"].transform(lambda s: s / s.sum())
    tb_lookup = {}
    for mouse, g in tb.groupby("MouseID_D1D2_BIseq"):
        tb_lookup[mouse] = g.set_index("RNA_barcode")["freq_tb"].to_dict()

    per_is = per_is.copy()
    gse_mouse = per_is[gse_mouse_col]

    per_is["Transcription_Organ"] = [
        organ_lookup.get((m, o, bc), 0) if pd.notna(m) else np.nan
        for m, o, bc in zip(gse_mouse, per_is["Organ"], per_is["barcode"])
    ]
    per_is["Transcription_Animal"] = [
        animal_lookup.get(m, {}).get(bc, 0) if pd.notna(m) else np.nan
        for m, bc in zip(gse_mouse, per_is["barcode"])
    ]
    per_is["Transcription_TB"] = [
        tb_lookup.get(m, {}).get(bc, 0) if pd.notna(m) else np.nan
        for m, bc in zip(gse_mouse, per_is["barcode"])
    ]
    per_is["Viremia"] = per_is["Transcription_TB"] > 0
    return per_is


def build_perIS_lite():
    gse_index = load_gse_index()
    per_is_raw = build_per_site(gse_index)
    per_is = assign_pair_to_dominant_animal(per_is_raw)
    per_is = add_proliferation(per_is)
    per_is = add_rna_columns(per_is)

    final_cols = [
        "site", "AnimalKey", "Organ", "Treatment",
        "UMI", "barcode", "Proliferation", "Transcription_Organ", "Transcription_Animal",
        "Transcription_TB", "Viremia", "Proliferation_tissue",
    ]
    return per_is[final_cols].sort_values(["AnimalKey", "site"]).reset_index(drop=True)


def main():
    os.makedirs(os.path.dirname(OUT_TABLE), exist_ok=True)
    build_perIS_lite().to_csv(OUT_TABLE, index=False)


if __name__ == "__main__":
    main()
