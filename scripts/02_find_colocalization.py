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
    parser.add_argument("--base-dir", required=True, help="Base results directory (e.g. results/pilot_run)")
    parser.add_argument("--sample", required=True, help="Sample ID to process")
    parser.add_argument("--sample-map", default="data/metadata/sample_map.csv",
                        help="Path ke sample_map.csv (country/region mapping)")
    args = parser.parse_args()

    base_dir = args.base_dir.rstrip("/")
    sample_id = args.sample

    # Load sample map
    try:
        sample_map = pd.read_csv(args.sample_map)
        # Drop duplicate sample_id rows: to_dict("index") requires a unique index,
        # and the Snakefile itself already deduplicates via df_map["sample_id"].unique().
        sample_map = sample_map.drop_duplicates(subset="sample_id", keep="first")
        map_dict = sample_map.set_index("sample_id").to_dict("index")
    except Exception as e:
        print(f"  [ERROR] Gagal membaca sample_map.csv: {e}", file=sys.stderr)
        sys.exit(1)

    # --- Metadata ---
    meta = map_dict.get(sample_id, {})
    country      = meta.get("country", "Unknown")
    country_name = meta.get("country_name", "Unknown")
    region       = meta.get("region", "Unknown")
    age          = meta.get("age", "Unknown")
    gender       = meta.get("gender", "Unknown")
    study_name   = meta.get("study_name", "Unknown")

    # --- Path Input ---
    rgi_file = f"{base_dir}/analysis/rgi/{sample_id}.tsv"
    mob_file = f"{base_dir}/analysis/mobsuite/{sample_id}/contig_report.txt"
    integron_file = f"{base_dir}/analysis/integron/{sample_id}.tsv"
    is_file = f"{base_dir}/analysis/isescan/{sample_id}.tsv"

    print(f"\n  [PROCESS] Sample: {sample_id}")

    # --- Load RGI (Wajib ada) ---
    rgi_df = load_rgi(rgi_file)
    if rgi_df.empty:
        print(f"    [ERROR] RGI file tidak ditemukan atau kosong: {rgi_file}")
        sys.exit(1)

    # --- Load MGEs (Opsional, kalau tidak ada kita anggap kosong) ---
    mob_contigs, _ = load_mobsuite(mob_file) if os.path.exists(mob_file) else (set(), None)
    integron_contigs, _ = load_integronfinder(integron_file) if os.path.exists(integron_file) else (set(), None)
    is_contigs, _ = load_isescan(is_file) if os.path.exists(is_file) else (set(), None)

    # --- Anotasi tiap gen AMR ---
    rgi_df["On_Plasmid"]  = rgi_df["contig_clean"].isin(mob_contigs)
    rgi_df["On_Integron"] = rgi_df["contig_clean"].isin(integron_contigs)
    rgi_df["On_IS"]       = rgi_df["contig_clean"].isin(is_contigs)
    rgi_df["MGE_Type"]    = rgi_df.apply(determine_mge_type, axis=1)

    # --- Tambah metadata ---
    rgi_df["Sample_ID"]    = sample_id
    rgi_df["Country"]      = country
    rgi_df["Country_Name"] = country_name
    rgi_df["Region"]       = region
    rgi_df["Age"]          = age
    rgi_df["Gender"]       = gender
    rgi_df["Study"]        = study_name

    # Kolom yang akan disimpan
    keep_cols = [
        "Sample_ID", "Country", "Country_Name", "Region", "Age", "Gender", "Study",
        "Contig", "contig_clean", "Best_Hit_ARO",
        "Drug Class", "Resistance Mechanism", "AMR Gene Family",
        "Best_Identities", "On_Plasmid", "On_Integron", "On_IS", "MGE_Type"
    ]
    available = [c for c in keep_cols if c in rgi_df.columns]
    final_df = rgi_df[available]

    # --- Simpan ---
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    final_df.to_csv(args.out, index=False)
    print(f"  [DONE] Colocalization saved: {args.out} ({len(final_df)} genes)")


if __name__ == "__main__":
    main()
