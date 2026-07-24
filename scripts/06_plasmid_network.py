#!/usr/bin/env python3
"""
Script: 06_plasmid_network.py
Purpose: Plasmid Replicon Typing & "Pandemic Plasmid" Network Analysis

Mengambil output MOB-suite (contig_report.txt) dari semua sampel dan melakukan:
1. Membangun Replicon Frequency Matrix per region/country
2. Mendeteksi "Pandemic Plasmid Candidates" — replicon yang muncul di ≥2 region berbeda
3. Membangun bipartite network: Replicon ↔ Country
4. Menghitung centrality metrics (degree, betweenness) per replicon

Input:
  - results/{run_id}/analysis/mobsuite/{sample}/contig_report.txt (semua sampel)
  - data/metadata/sample_map.csv

Output:
  - results/{run_id}/analysis/plasmid_network/{pop}/replicon_frequency.csv
  - results/{run_id}/analysis/plasmid_network/{pop}/pandemic_candidates.csv
  - results/{run_id}/analysis/plasmid_network/{pop}/replicon_country_network.csv
  - results/{run_id}/analysis/plasmid_network/{pop}/network_metrics.csv
  - results/{run_id}/analysis/plasmid_network/{pop}/network_done.flag
"""

import argparse
import os
import sys
import glob
import pandas as pd
from collections import defaultdict


# ============================================================
# Helper functions
# ============================================================

def load_mobsuite_report(filepath, sample_id, country, country_name, region):
    """Membaca contig_report.txt dari MOB-suite dan mengembalikan DataFrame plasmid."""
    try:
        df = pd.read_csv(filepath, sep="\t")
        if df.empty:
            return pd.DataFrame()

        # Filter hanya kontig yang teridentifikasi sebagai plasmid
        plasmid_df = df[df["molecule_type"].str.lower() == "plasmid"].copy()
        if plasmid_df.empty:
            return pd.DataFrame()

        plasmid_df["sample_id"]    = sample_id
        plasmid_df["country"]      = country
        plasmid_df["country_name"] = country_name
        plasmid_df["region"]       = region

        return plasmid_df

    except Exception as e:
        print(f"  [WARN] Gagal membaca {filepath}: {e}", file=sys.stderr)
        return pd.DataFrame()


def parse_replicon_types(rep_type_str):
    """
    Memecah string rep_type(s) dari MOB-suite menjadi list.
    Contoh: 'IncF,IncI1' -> ['IncF', 'IncI1']
    """
    if pd.isna(rep_type_str) or str(rep_type_str).strip() in ["", "-", "NA"]:
        return ["Unknown"]
    return [r.strip() for r in str(rep_type_str).split(",") if r.strip()]


def build_replicon_frequency(all_plasmids_df):
    """Membangun Replicon Frequency Matrix: berapa kali setiap replicon muncul per country."""
    rows = []
    for _, row in all_plasmids_df.iterrows():
        replicons = parse_replicon_types(row.get("rep_type(s)", "Unknown"))
        for rep in replicons:
            rows.append({
                "replicon_type": rep,
                "sample_id":     row["sample_id"],
                "country":       row["country"],
                "country_name":  row.get("country_name", row["country"]),
                "region":        row["region"],
                "ptu":           row.get("primary_cluster_id", "Unknown"),
                "relaxase_type": row.get("relaxase_type(s)", "Unknown"),
            })

    if not rows:
        return pd.DataFrame()

    freq_df = pd.DataFrame(rows)
    return freq_df


def detect_pandemic_candidates(freq_df, min_regions=2):
    """
    Mendeteksi "Pandemic Plasmid Candidates" — replicon yang muncul di >= min_regions region berbeda.
    Ini adalah replicon yang paling mungkin menyebar lintas benua.
    """
    if freq_df.empty:
        return pd.DataFrame()

    # Hitung di berapa region berbeda setiap replicon muncul
    region_spread = (
        freq_df[freq_df["replicon_type"] != "Unknown"]
        .groupby("replicon_type")
        .agg(
            n_regions=("region", "nunique"),
            regions_found=("region", lambda x: "|".join(sorted(x.unique()))),
            n_countries=("country", "nunique"),
            countries_found=("country", lambda x: "|".join(sorted(x.unique()))),
            total_occurrences=("sample_id", "count"),
            ptu=("ptu", lambda x: "|".join(sorted(x.dropna().unique())))
        )
        .reset_index()
        .sort_values("n_regions", ascending=False)
    )

    pandemic_candidates = region_spread[region_spread["n_regions"] >= min_regions].copy()
    pandemic_candidates["pandemic_candidate"] = True

    return pandemic_candidates


def build_bipartite_network(freq_df):
    """
    Membangun edge list bipartite network: Replicon <-> Country.
    Weight = jumlah sampel yang memiliki replicon tersebut di negara tersebut.
    """
    if freq_df.empty:
        return pd.DataFrame()

    network_df = (
        freq_df[freq_df["replicon_type"] != "Unknown"]
        .groupby(["replicon_type", "country", "country_name", "region"])
        .agg(weight=("sample_id", "count"))
        .reset_index()
    )
    return network_df


def compute_network_metrics(network_df):
    """
    Menghitung metrik sederhana untuk setiap replicon:
    - degree: jumlah negara berbeda yang memiliki replicon ini
    - total_weight: total kemunculan di semua negara
    """
    if network_df.empty:
        return pd.DataFrame()

    metrics = (
        network_df.groupby("replicon_type")
        .agg(
            degree=("country", "nunique"),
            total_weight=("weight", "sum"),
            regions=("region", lambda x: "|".join(sorted(x.unique())))
        )
        .reset_index()
        .sort_values("degree", ascending=False)
    )
    return metrics


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Plasmid Replicon Typing & Pandemic Plasmid Network Analysis"
    )
    parser.add_argument("--base-dir",   required=True, help="Base results directory (e.g. results/pilot_run)")
    parser.add_argument("--output-dir", required=True, help="Output directory for this population")
    parser.add_argument("--pop",        required=True, help="Population/region name (e.g. East_Asia, combined)")
    parser.add_argument("--sample-map", default="data/metadata/pilot_map.csv",
                        help="Path ke sample_map.csv")
    parser.add_argument("--min-regions", type=int, default=2,
                        help="Minimum jumlah region untuk disebut pandemic candidate (default: 2)")
    parser.add_argument("--flag-file",  required=True, help="Path ke flag file output")
    args = parser.parse_args()

    base_dir   = args.base_dir.rstrip("/")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)

    print("=" * 60)
    print(f"[06_plasmid_network] Population: {args.pop}")
    print("=" * 60)

    # --- Load sample map ---
    try:
        sample_map = pd.read_csv(args.sample_map)
        # Drop duplicate sample_id rows so set_index(...).to_dict("index") stays valid.
        sample_map = sample_map.drop_duplicates(subset="sample_id", keep="first")
        map_dict   = sample_map.set_index("sample_id").to_dict("index")
    except Exception as e:
        print(f"[ERROR] Gagal membaca sample_map: {e}", file=sys.stderr)
        sys.exit(1)

    # Filter sampel berdasarkan populasi jika bukan "combined"
    if args.pop.lower() not in ["combined", "all"]:
        region_col = "Region" if "Region" in sample_map.columns else "region"
        sample_map = sample_map[sample_map[region_col] == args.pop]

    samples = sample_map["sample_id"].tolist()
    print(f"[INFO] Jumlah sampel yang diproses: {len(samples)}")

    # --- Kumpulkan semua data plasmid dari MOB-suite ---
    all_plasmids = []
    for sample_id in samples:
        meta         = map_dict.get(sample_id, {})
        country      = meta.get("country", "Unknown")
        country_name = meta.get("country_name", "Unknown")
        region_key   = "Region" if "Region" in meta else "region"
        region       = meta.get(region_key, meta.get("region", "Unknown"))

        mob_file = f"{base_dir}/analysis/mobsuite/{sample_id}/contig_report.txt"
        if not os.path.exists(mob_file):
            print(f"  [SKIP] File tidak ditemukan: {mob_file}")
            continue

        df = load_mobsuite_report(mob_file, sample_id, country, country_name, region)
        if not df.empty:
            all_plasmids.append(df)

    if not all_plasmids:
        print("[WARN] Tidak ada data plasmid yang berhasil dikumpulkan. Output akan kosong.")
        # Buat file kosong agar Snakemake tidak fail
        pd.DataFrame().to_csv(f"{output_dir}/replicon_frequency.csv", index=False)
        pd.DataFrame().to_csv(f"{output_dir}/pandemic_candidates.csv", index=False)
        pd.DataFrame().to_csv(f"{output_dir}/replicon_country_network.csv", index=False)
        pd.DataFrame().to_csv(f"{output_dir}/network_metrics.csv", index=False)
        open(args.flag_file, "w").close()
        return

    all_plasmids_df = pd.concat(all_plasmids, ignore_index=True)
    print(f"[INFO] Total baris plasmid: {len(all_plasmids_df)}")

    # --- 1. Replicon Frequency ---
    print("\n[1] Membangun Replicon Frequency Matrix...")
    freq_df = build_replicon_frequency(all_plasmids_df)
    freq_df.to_csv(f"{output_dir}/replicon_frequency.csv", index=False)
    print(f"    -> {len(freq_df)} baris tersimpan di replicon_frequency.csv")

    # --- 2. Pandemic Candidates ---
    print(f"\n[2] Mendeteksi Pandemic Plasmid Candidates (min {args.min_regions} region)...")
    pandemic_df = detect_pandemic_candidates(freq_df, min_regions=args.min_regions)
    pandemic_df.to_csv(f"{output_dir}/pandemic_candidates.csv", index=False)
    if not pandemic_df.empty:
        print(f"    -> {len(pandemic_df)} pandemic candidate(s) ditemukan:")
        for _, row in pandemic_df.head(10).iterrows():
            print(f"       * {row['replicon_type']} — {row['n_regions']} regions ({row['regions_found']})")
    else:
        print("    -> Tidak ada pandemic candidate yang memenuhi threshold.")

    # --- 3. Bipartite Network ---
    print("\n[3] Membangun Bipartite Network (Replicon ↔ Country)...")
    network_df = build_bipartite_network(freq_df)
    network_df.to_csv(f"{output_dir}/replicon_country_network.csv", index=False)
    print(f"    -> {len(network_df)} edges tersimpan di replicon_country_network.csv")

    # --- 4. Network Metrics ---
    print("\n[4] Menghitung Network Centrality Metrics...")
    metrics_df = compute_network_metrics(network_df)
    metrics_df.to_csv(f"{output_dir}/network_metrics.csv", index=False)
    print(f"    -> {len(metrics_df)} replicon type dianalisis")

    # --- Summary Stats ---
    print("\n" + "=" * 60)
    print("[SUMMARY]")
    if not freq_df.empty:
        top_replicons = (
            freq_df[freq_df["replicon_type"] != "Unknown"]
            .groupby("replicon_type")["sample_id"].count()
            .sort_values(ascending=False)
            .head(5)
        )
        print(f"  Top 5 Replicon Types (by occurrence):")
        for rep, count in top_replicons.items():
            print(f"    {rep}: {count} occurrences")

    print(f"\n  Pandemic candidates: {len(pandemic_df)}")
    print(f"  Output dir: {output_dir}")
    print("=" * 60)

    # --- Done Flag ---
    open(args.flag_file, "w").close()
    print(f"\n[DONE] Flag file dibuat: {args.flag_file}")


if __name__ == "__main__":
    main()
