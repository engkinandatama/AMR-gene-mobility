#!/usr/bin/env python3
"""
Script: 02_find_colocalization.py
Purpose: Mengintegrasikan hasil deteksi AMR (RGI) dengan deteksi MGE dari 3 tools berbeda
         (MOBsuite, IntegronFinder, ISEScan) untuk menentukan posisi setiap gen AMR
         apakah berada di Plasmid, Integron, Transposon/IS, atau Kromosom.
         Kolom country & region diekstrak dari sample_map.csv.

Input:
  - results/rgi/{sample}_rgi.txt
  - results/mobsuite/{sample}_plasmid.txt
  - results/integron/{sample}_integrons.tsv
  - results/isescan/{sample}_is.tsv
  - data/metadata/sample_map.csv

Output:
  - results/colocalization_summary.csv
"""

import pandas as pd
import argparse
import glob
import os
import sys

# ============================================================
# Helper functions
# ============================================================

def load_rgi(rgi_file):
    """Membaca output RGI, mengembalikan DataFrame. Return kosong jika gagal."""
    try:
        df = pd.read_csv(rgi_file, sep="\t")
        if df.empty:
            return pd.DataFrame()
        # Bersihkan nama contig: RGI kadang menambahkan suffix '_1' dll
        df["contig_clean"] = df["Contig"].astype(str).apply(
            lambda x: "_".join(x.split("_")[:-1]) if x.count("_") > 1 else x
        )
        return df
    except Exception as e:
        print(f"  [WARN] Gagal membaca RGI file {rgi_file}: {e}", file=sys.stderr)
        return pd.DataFrame()


def load_mobsuite(mob_file):
    """Membaca output MOBsuite, mengembalikan set contig_id yang teridentifikasi PLASMID."""
    try:
        df = pd.read_csv(mob_file, sep="\t")
        plasmid_df = df[df["molecule_type"].str.lower() == "plasmid"]
        plasmid_contigs = set(plasmid_df["contig_id"].astype(str))
        return plasmid_contigs, plasmid_df[["contig_id", "rep_type(s)", "relaxase_type(s)"]].rename(
            columns={"contig_id": "contig_clean", "rep_type(s)": "plasmid_rep_type", "relaxase_type(s)": "plasmid_relaxase"}
        )
    except Exception as e:
        print(f"  [WARN] Gagal membaca MOBsuite file {mob_file}: {e}", file=sys.stderr)
        return set(), pd.DataFrame()


def load_integronfinder(integron_file):
    """Membaca output IntegronFinder, mengembalikan set contig_id yang mengandung integron."""
    try:
        df = pd.read_csv(integron_file, sep="\t", comment="#")
        if df.empty:
            return set(), pd.DataFrame()
        # Kolom khas IntegronFinder v2: ID_replicon
        if "ID_replicon" not in df.columns:
            return set(), pd.DataFrame()
        
        # Semua baris di file *.integrons adalah bagian dari integron
        integron_contigs = set(df["ID_replicon"].astype(str))
        summary = df[["ID_replicon"]].drop_duplicates().rename(
            columns={"ID_replicon": "contig_clean"}
        )
        summary["integron_type"] = "Integron"
        return integron_contigs, summary
    except Exception as e:
        print(f"  [WARN] Gagal membaca IntegronFinder file {integron_file}: {e}", file=sys.stderr)
        return set(), pd.DataFrame()


def load_isescan(is_file):
    """Membaca output ISEScan, mengembalikan set contig_id yang mengandung IS/Transposon."""
    try:
        df = pd.read_csv(is_file, sep="\t")
        if df.empty:
            return set(), pd.DataFrame()
        # Rename seqID to seqid to handle case-insensitivity
        df.rename(columns=lambda x: x.lower(), inplace=True)
        # Kolom khas ISEScan: seqid (nama contig), family, cluster
        if "seqid" not in df.columns:
            return set(), pd.DataFrame()
        is_contigs = set(df["seqid"].astype(str))
        summary = df[["seqid", "family"]].drop_duplicates().rename(
            columns={"seqid": "contig_clean", "family": "is_family"}
        )
        return is_contigs, summary
    except Exception as e:
        print(f"  [WARN] Gagal membaca ISEScan file {is_file}: {e}", file=sys.stderr)
        return set(), pd.DataFrame()


def determine_mge_type(row):
    """
    Menentukan tipe MGE berdasarkan hirarki:
    Plasmid > Integron > IS/Transposon > Chromosomal
    """
    if row["On_Plasmid"]:
        return "Plasmid"
    elif row["On_Integron"]:
        return "Integron"
    elif row["On_IS"]:
        return "IS/Transposon"
    else:
        return "Chromosomal"


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Integrasi AMR + MGE untuk analisis co-localization dan HGT linkage"
    )
    parser.add_argument("--out", required=True, help="Path file CSV output")
    parser.add_argument("--sample-map", default="data/metadata/sample_map.csv",
                        help="Path ke sample_map.csv (country/region mapping)")
    args = parser.parse_args()

    # Load sample map
    print("  [INFO] Membaca sample_map.csv...")
    try:
        sample_map = pd.read_csv(args.sample_map)
        # Index: sample_id -> dict
        map_dict = sample_map.set_index("sample_id").to_dict("index")
        print(f"  [INFO] {len(map_dict)} sampel ditemukan di sample_map.")
    except Exception as e:
        print(f"  [ERROR] Gagal membaca sample_map.csv: {e}", file=sys.stderr)
        sys.exit(1)

    # Temukan semua file RGI yang ada
    rgi_files = glob.glob("results/rgi/*_rgi.txt")
    if not rgi_files:
        print("  [ERROR] Tidak ada file RGI ditemukan di results/rgi/.", file=sys.stderr)
        sys.exit(1)

    all_results = []

    for rgi_file in sorted(rgi_files):
        sample_id = os.path.basename(rgi_file).replace("_rgi.txt", "")
        print(f"\n  [SAMPLE] {sample_id}")

        # --- Metadata ---
        meta = map_dict.get(sample_id, {})
        country      = meta.get("country", "Unknown")
        country_name = meta.get("country_name", "Unknown")
        region       = meta.get("region", "Unknown")

        # --- Load semua tools ---
        rgi_df = load_rgi(rgi_file)
        if rgi_df.empty:
            print(f"    [SKIP] RGI kosong, melewati sampel ini.")
            continue

        mob_contigs, mob_meta = load_mobsuite(
            f"results/mobsuite/{sample_id}_plasmid.txt"
        )
        integron_contigs, _ = load_integronfinder(
            f"results/integron/{sample_id}_integrons.tsv"
        )
        is_contigs, _ = load_isescan(
            f"results/isescan/{sample_id}_is.tsv"
        )

        print(f"    AMR genes    : {len(rgi_df)}")
        print(f"    Plasmid ctg  : {len(mob_contigs)}")
        print(f"    Integron ctg : {len(integron_contigs)}")
        print(f"    IS/Tn ctg    : {len(is_contigs)}")

        # --- Anotasi tiap gen AMR ---
        rgi_df["On_Plasmid"]  = rgi_df["contig_clean"].isin(mob_contigs)
        rgi_df["On_Integron"] = rgi_df["contig_clean"].isin(integron_contigs)
        rgi_df["On_IS"]       = rgi_df["contig_clean"].isin(is_contigs)
        rgi_df["MGE_Type"]    = rgi_df.apply(determine_mge_type, axis=1)

        # --- Tambah metadata populasi ---
        rgi_df["Sample_ID"]    = sample_id
        rgi_df["Country"]      = country
        rgi_df["Country_Name"] = country_name
        rgi_df["Region"]       = region

        # Kolom yang akan disimpan (subset dari output RGI)
        keep_cols = [
            "Sample_ID", "Country", "Country_Name", "Region",
            "Contig", "contig_clean", "Best_Hit_ARO",
            "Drug Class", "Resistance Mechanism", "AMR Gene Family",
            "Best_Identities", "On_Plasmid", "On_Integron", "On_IS", "MGE_Type"
        ]
        # Ambil kolom yang tersedia saja (nama kolom RGI bisa berbeda antar versi)
        available = [c for c in keep_cols if c in rgi_df.columns]
        all_results.append(rgi_df[available])

        mobile_count = rgi_df["On_Plasmid"].sum() + rgi_df["On_Integron"].sum() + rgi_df["On_IS"].sum()
        print(f"    -> {mobile_count}/{len(rgi_df)} gen AMR terkait MGE")

    # --- Gabung dan simpan ---
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        final_df.to_csv(args.out, index=False)
        print(f"\n  [DONE] Colocalization summary tersimpan: {args.out}")
        print(f"         Total baris: {len(final_df)}")
        print(f"         Total sampel: {final_df['Sample_ID'].nunique()}")
    else:
        print("\n  [WARN] Tidak ada sampel yang berhasil diproses.")
        pd.DataFrame().to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
