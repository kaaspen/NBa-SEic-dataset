#!/usr/bin/env python
# pipeline.py
#
# Minimal script that reproduces the four output tables (Set2__Table1–4.csv)
# from ROSE-generated peak files.  Usage:
#
#   python pipeline.py --data_dir  /path/to/rose_bed_files   \
#                      --out_dir   results
#
# The directory given by --data_dir must contain, for every sample,
#   SE_<sample>_SE_mm10.bed   and   SE_<sample>_TE_mm10.bed
# obtained with the standard ROSE pipeline.

import argparse
import pathlib
import subprocess
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

# ────────────────────────────── CLI ──────────────────────────────
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate non-binary SE dataset tables from ROSE BED peaks")
    p.add_argument("--data_dir", required=True,
                   help="Folder with ROSE BED files (SE_*_SE_mm10.bed / TE_*.bed)")
    p.add_argument("--out_dir", required=True,
                   help="Folder to store Tables/ and Intermediate_tables/")
    return p.parse_args()


ARGS = parse_args()
DATA_DIR = pathlib.Path(ARGS.data_dir)
OUT_DIR = pathlib.Path(ARGS.out_dir)
TABLES = OUT_DIR / "Tables"
INTER = OUT_DIR / "Intermediate_tables"
TABLES.mkdir(parents=True, exist_ok=True)
INTER.mkdir(exist_ok=True)

# ────────────────────────────── I/O helpers ──────────────────────
def run_cmd(cmd: str):
    """Run shell command, raise if it fails"""
    subprocess.run(cmd, shell=True, check=True)


def bed_sort(in_bed: pathlib.Path, out_bed: pathlib.Path):
    run_cmd(f"sort -k1,1 -k2,2n {in_bed} > {out_bed}")


def bed_merge_sorted(sorted_bed: pathlib.Path, out_bed: pathlib.Path):
    run_cmd(
        f"bedtools merge -i {sorted_bed} "
        "-c 4,4,5,5 -o count_distinct,collapse,count_distinct,collapse "
        "-d 12500 > {out_bed}".format(out_bed=out_bed)
    )


def bed_intersect(a_bed: pathlib.Path, b_bed: pathlib.Path, out_bed: pathlib.Path):
    run_cmd(
        f"bedtools intersect -a {a_bed} -b {b_bed} "
        "-wo -f 0.5 -F 1 -e > {out_bed}"
    )


# ────────────────────────────── Workflow ────────────────────────
def main() -> None:
    samples = ["11_0077", "12_0118", "12_0449", "12_0450"]

    SE_table, TE_table, df_ste = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    # -------- 1.  Read individual peak files ----------
    for s in samples:
        se_bed = DATA_DIR / f"SE_{s}_SE_mm10.bed"
        te_bed = DATA_DIR / f"SE_{s}_TE_mm10.bed"

        df_se = pd.read_csv(se_bed, sep="\t")
        df_te = pd.read_csv(te_bed, sep="\t")

        df_se["region_length"] = df_se.se_end - df_se.se_start
        df_te["region_length"] = df_te.TE_end - df_te.TE_start

        # simple cleaning
        df_te = (
            df_te.dropna(subset=["cell_id", "TE_cas_value"])
            .astype({"TE_con_value": float})
            .query("TE_con_value < TE_cas_value")
            .drop_duplicates()
        )

        df_te["avg_rpm_diff"] = df_te.TE_cas_value - df_te.TE_con_value
        df_se["avg_rpm_diff"] = df_se.se_cas_value - df_se.se_con_value

        df_ste_curr = pd.concat(
            [
                df_se[
                    ["cell_id", "se_id", "se_chr", "se_start", "se_end",
                     "se_rank", "region_length", "avg_rpm_diff"]
                ].rename(
                    columns={
                        "se_id": "ste_id",
                        "se_rank": "ste_rank",
                        "se_chr": "ste_chr",
                        "se_start": "ste_start",
                        "se_end": "ste_end",
                    }
                ),
                df_te[
                    ["cell_id", "TE_id", "TE_chr", "TE_start", "TE_end",
                     "TE_rank", "region_length", "avg_rpm_diff"]
                ].rename(
                    columns={
                        "TE_id": "ste_id",
                        "TE_rank": "ste_rank",
                        "TE_chr": "ste_chr",
                        "TE_start": "ste_start",
                        "TE_end": "ste_end",
                    }
                ),
            ],
            ignore_index=True,
        )

        SE_table = pd.concat([SE_table, df_se], ignore_index=True)
        TE_table = pd.concat([TE_table, df_te], ignore_index=True)
        df_ste = pd.concat([df_ste, df_ste_curr], ignore_index=True)

    # -------- 2.  Consolidate SE loci -----------------
    concat_bed = INTER / "Set2__SE_concatinated.bed"
    SE_table[["se_chr", "se_start", "se_end", "cell_id", "se_id"]].to_csv(
        concat_bed, sep="\t", header=False, index=False
    )

    sorted_bed = INTER / "sorted.Set2__SE_concatinated.bed"
    merged_bed = INTER / "merged.sorted.Set2__SE_concatinated.bed"
    bed_sort(concat_bed, sorted_bed)
    bed_merge_sorted(sorted_bed, merged_bed)

    df_merged = pd.read_csv(
        merged_bed,
        sep="\t",
        header=None,
        names=[
            "se_chr", "se_start", "se_end",
            "cell_id_count", "cell_id_list", "se_id_count", "se_id_list",
        ],
    )

    for s in samples:
        df_merged[f"SE_{s}"] = df_merged.cell_id_list.str.contains(f"SE_{s}").astype(int)

    df_merged["se_locus_id"] = [f"se_region_{i+1}" for i in range(len(df_merged))]
    df_merged = df_merged[
        ["se_locus_id", "se_chr", "se_start", "se_end", "se_id_list"]
        + [f"SE_{s}" for s in samples]
    ].rename(columns={"se_chr": "chr", "se_start": "start", "se_end": "end"})
    df_merged.to_csv(TABLES / "Set2__Table1.csv", index=False)

    # -------- 3.  Intersect with union SE+TE ----------
    loci_bed = INTER / "Set2__consolidated_SE_loci.bed"
    df_merged[["chr", "start", "end", "se_locus_id"]].to_csv(
        loci_bed, sep="\t", header=False, index=False
    )

    ste_bed = INTER / "Set2__SE_TE_concatinated.bed"
    df_ste[
        ["ste_chr", "ste_start", "ste_end",
         "cell_id", "ste_id", "ste_rank", "region_length", "avg_rpm_diff"]
    ].to_csv(ste_bed, sep="\t", header=False, index=False)

    ste_sorted = INTER / "sorted.Set2__SE_TE_concatinated.bed"
    bed_sort(ste_bed, ste_sorted)

    intersect_bed = INTER / "Set2__STE_within_consolidated_SE_loci.bed"
    bed_intersect(ste_sorted, loci_bed, intersect_bed)

    df_ste_within = pd.read_csv(
        intersect_bed,
        sep="\t",
        header=None,
        names=[
            "ste_chr", "ste_start", "ste_end", "cell_id",
            "ste_id", "ste_rank", "region_length", "avg_rpm_diff",
            "chr", "start", "end", "se_locus_id", "overlap",
        ],
    )

    df_ste_within = df_ste_within.merge(
        df_ste_within.groupby(["se_locus_id", "cell_id"], as_index=False)
        .agg(locus_ste_overlap_total=("overlap", np.sum)),
        on=["se_locus_id", "cell_id"],
    )
    df_ste_within["ste_weight_within_locus"] = (
        df_ste_within.overlap / df_ste_within.locus_ste_overlap_total
    )

    cols = [
        "se_locus_id", "cell_id", "ste_id", "ste_chr", "ste_start", "ste_end",
        "ste_rank", "avg_rpm_diff", "overlap", "ste_weight_within_locus",
    ]
    df_ste_within[cols].to_csv(TABLES / "Set2__Table2.csv", index=False)

    # -------- 4.  Locus-level activity ----------------
    df_ste_within["is_SE"] = df_ste_within.ste_id.str.contains("SE").astype(int)
    df_ste_within["weighted_avg_rpm_diff"] = (
        df_ste_within.avg_rpm_diff * df_ste_within.ste_weight_within_locus
    )

    agg = df_ste_within.groupby(["se_locus_id", "cell_id"], as_index=False).agg(
        avg_rpm_diff__max=("avg_rpm_diff", "max"),
        avg_rpm_diff__weighted=("weighted_avg_rpm_diff", "sum"),
        max_rank=("ste_rank", "max"),
        min_rank=("ste_rank", "min"),
        active_SE=("is_SE", "max"),
        active_SE_count=("is_SE", "sum"),
        active_TE_count=("is_SE", lambda x: (1 - x).sum()),
    )
    df_locus_activity = df_merged[["se_locus_id"]].merge(agg, on="se_locus_id")
    df_locus_activity.to_csv(TABLES / "Set2__Table3.csv", index=False)

    # -------- 5.  Feature matrix ----------------------
    df_features = df_locus_activity.pivot_table(
        index="se_locus_id", columns="cell_id", values="avg_rpm_diff__max"
    ).sort_index()

    # median normalisation, imputation, log1p, z-score
    for col in df_features.columns:
        med = df_features[col].median()
        df_features[f"{col}_medianNormalized"] = df_features[col] / med

    imputed = np.nan_to_num(df_features.iloc[:, len(samples):].values, nan=1e-6)
    df_features[df_features.columns[len(samples):] + "_imputed"] = imputed

    X = np.log1p(imputed)
    scaler = StandardScaler().fit(X)
    df_features[df_features.columns[len(samples)*2:] + "_zscaled"] = scaler.transform(X)

    # add binary presence
    df_features = df_features.reset_index().merge(
        df_merged[["se_locus_id"] + [f"SE_{s}" for s in samples]],
        on="se_locus_id",
    )

    df_features.to_csv(TABLES / "Set2__Table4.csv", index=False)
    print("✓ Finished.  Tables written to", TABLES)


# ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    main()
