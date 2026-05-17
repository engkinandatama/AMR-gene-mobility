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

import argparse
import pandas as pd
import os
import sys

def main():
    parser = argparse.ArgumentParser(description="Agregasi data colocalization per populasi")
    parser.add_argument("--input", required=True, help="Path ke colocalization_summary.csv")
    parser.add_argument("--output-dir", required=True, help="Folder untuk menyimpan hasil agregasi")
    parser.add_argument("--pop", default="combined", help="Nama populasi (untuk prefix file)")
    args = parser.parse_args()

    input_file = args.input
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print(f"[Phase 3C] Agregasi data populasi: {args.pop}")
    print("=" * 60)

    # --- Load data ---
    print(f"\n[1] Membaca {input_file}...")
    try:
        df = pd.read_csv(input_file)
    except FileNotFoundError:
        print(f"[ERROR] File tidak ditemukan: {input_file}", file=sys.stderr)
        sys.exit(1)

    if df.empty:
        print("[ERROR] File input kosong.", file=sys.stderr)
        sys.exit(1)

    # Filter data berdasarkan Region/Populasi jika bukan "combined" atau "all"
    if args.pop.lower() not in ["combined", "all"]:
        print(f"    [INFO] Memfilter data untuk populasi (Region): {args.pop}")
        df = df[df["Region"] == args.pop]
        if df.empty:
            print(f"    [WARN] Tidak ada data sampel yang ditemukan untuk Region '{args.pop}'.")

    # -------------------------------------------------------
    # OUTPUT 1: AMR Abundance Matrix
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
    amr_abundance_pivot.to_csv(f"{output_dir}/{args.pop}_amr_abundance.csv", index=False)

    # -------------------------------------------------------
    # OUTPUT 2: MGE Distribution Matrix
    # -------------------------------------------------------
    print("[3] Membuat MGE Distribution Matrix...")
    mge_dist = (
        df.groupby(["Country", "Region", "MGE_Type"])
        .size()
        .reset_index(name="Count")
    )
    mge_total = mge_dist.groupby("Country")["Count"].transform("sum")
    mge_dist["Proportion"] = mge_dist["Count"] / mge_total
    mge_dist.to_csv(f"{output_dir}/{args.pop}_mge_distribution.csv", index=False)

    # -------------------------------------------------------
    # OUTPUT 3: AMR × MGE Association Matrix
    # -------------------------------------------------------
    print("[4] Membuat AMR-MGE Association Matrix...")
    amr_mge = (
        df.groupby(["Country", "Region", "Drug Class", "MGE_Type"])
        .size()
        .reset_index(name="Count")
    )
    amr_mge.to_csv(f"{output_dir}/{args.pop}_combined.csv", index=False) # Ini file utama yang diminta Snakefile
    print(f"    -> Hasil tersimpan di {output_dir}")

    print("\n[SELESAI]\n")

if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
