#!/usr/bin/env python3
"""
Script: 03_aggregate_by_population.py
Purpose: Merangkum colocalization_summary.csv menjadi tiga matrix agregat per populasi
         yang siap digunakan untuk analisis statistik di Phase 4.

Input:
  - results/colocalization_summary.csv

Output:
  - results/amr_abundance_matrix.csv     : jumlah hit per AMR class per sampel
  - results/mge_distribution_matrix.csv  : proporsi MGE type per negara
  - results/amr_mge_association_matrix.csv : cross-tabulation AMR class x MGE type x negara
"""

import pandas as pd
import os
import sys


def main():
    input_file = "results/colocalization_summary.csv"
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("[Phase 3C] Agregasi data per populasi")
    print("=" * 60)

    # --- Load data ---
    print(f"\n[1] Membaca {input_file}...")
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"[ERROR] File tidak ditemukan: {input_file}", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("[ERROR] colocalization_summary.csv kosong.", file=sys.stderr)
        sys.exit(1)

    print(f"    Total baris   : {len(df)}")
    print(f"    Total sampel  : {df['Sample_ID'].nunique()}")
    print(f"    Kolom         : {list(df.columns)}\n")

    # Pastikan kolom kritis ada
    required_cols = ["Sample_ID", "Country", "Region", "Drug Class", "MGE_Type"]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"[ERROR] Kolom berikut tidak ditemukan: {missing}", file=sys.stderr)
        sys.exit(1)

    # -------------------------------------------------------
    # OUTPUT 1: AMR Abundance Matrix (AMR class × sampel)
    # -------------------------------------------------------
    print("[2] Membuat AMR Abundance Matrix...")
    amr_abundance = (
        df.groupby(["Sample_ID", "Country", "Region", "Drug Class"])
        .size()
        .reset_index(name="Count")
    )
    amr_abundance_pivot = amr_abundance.pivot_table(
        index=["Sample_ID", "Country", "Region"],
        columns="Drug Class",
        values="Count",
        fill_value=0
    ).reset_index()
    amr_abundance_pivot.to_csv(f"{output_dir}/amr_abundance_matrix.csv", index=False)
    print(f"    -> {output_dir}/amr_abundance_matrix.csv  ({len(amr_abundance_pivot)} sampel x {len(amr_abundance_pivot.columns)} kolom)")

    # -------------------------------------------------------
    # OUTPUT 2: MGE Distribution Matrix (MGE type × negara)
    # -------------------------------------------------------
    print("[3] Membuat MGE Distribution Matrix...")
    mge_dist = (
        df.groupby(["Country", "Country_Name", "Region", "MGE_Type"])
        .size()
        .reset_index(name="Count")
    )
    # Hitung proporsi per negara
    mge_total = mge_dist.groupby("Country")["Count"].transform("sum")
    mge_dist["Proportion"] = mge_dist["Count"] / mge_total
    mge_dist.to_csv(f"{output_dir}/mge_distribution_matrix.csv", index=False)
    print(f"    -> {output_dir}/mge_distribution_matrix.csv  ({len(mge_dist)} baris)")

    # -------------------------------------------------------
    # OUTPUT 3: AMR × MGE Association Matrix per Negara
    # -------------------------------------------------------
    print("[4] Membuat AMR-MGE Association Matrix...")
    amr_mge = (
        df.groupby(["Country", "Region", "Drug Class", "MGE_Type"])
        .size()
        .reset_index(name="Count")
    )
    amr_mge.to_csv(f"{output_dir}/amr_mge_association_matrix.csv", index=False)
    print(f"    -> {output_dir}/amr_mge_association_matrix.csv  ({len(amr_mge)} baris)")

    # -------------------------------------------------------
    # SUMMARY PRINT
    # -------------------------------------------------------
    print("\n" + "=" * 60)
    print("RANGKUMAN HASIL:")
    print("=" * 60)

    print("\nTotal gen AMR per negara:")
    country_summary = df.groupby(["Country", "Country_Name", "Region"]).size().reset_index(name="Total_AMR_Hits")
    print(country_summary.to_string(index=False))

    print("\nDistribusi MGE Type (semua sampel):")
    mge_total_all = df["MGE_Type"].value_counts()
    mge_pct = (mge_total_all / mge_total_all.sum() * 100).round(1)
    for mge, count in mge_total_all.items():
        print(f"  {mge:20s}: {count:5d} ({mge_pct[mge]}%)")

    print("\n[SELESAI] Tiga matrix tersimpan di folder results/")
    print("          Siap untuk digunakan di 04_run_stats.R\n")


if __name__ == "__main__":
    main()
