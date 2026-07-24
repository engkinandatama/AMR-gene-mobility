#!/usr/bin/env python3
"""
Script: 07_integron_deep.py
Purpose: Deep Integron Analysis — Cassette Composition & Class Distribution

Mengambil output IntegronFinder dari semua sampel dan melakukan:
1. Klasifikasi integron: filter hanya tipe 'complete' untuk analisis utama
2. CALIN disimpan sebagai output terpisah (evidence for historical HGT)
3. Ekstraksi gene cassette composition per sampel/negara
4. Identifikasi cassette yang bersifat region-specific vs universal
5. Link ke AMR class dari data colocalization

Input:
  - results/{run_id}/analysis/integron/{sample}.tsv (semua sampel)
  - results/{run_id}/analysis/aggregated/{pop}_all_coloc.csv
  - data/metadata/sample_map.csv

Output:
  - results/{run_id}/analysis/integron_deep/{pop}/complete_integrons.csv   [analisis utama]
  - results/{run_id}/analysis/integron_deep/{pop}/calin_cassettes.csv      [historical HGT evidence]
  - results/{run_id}/analysis/integron_deep/{pop}/cassette_matrix.csv      [pivot: cassette x country]
  - results/{run_id}/analysis/integron_deep/{pop}/class_distribution.csv   [integron class per negara]
  - results/{run_id}/analysis/integron_deep/{pop}/region_specificity.csv   [universal vs specific cassettes]
  - results/{run_id}/analysis/integron_deep/{pop}/amr_integron_link.csv    [link cassette ke AMR class]
  - results/{run_id}/analysis/integron_deep/{pop}/integron_report_done.flag
"""

import argparse
import os
import sys
import pandas as pd
from collections import defaultdict


# ============================================================
# Helper functions
# ============================================================

def load_integron_file(filepath, sample_id, country, country_name, region):
    """
    Membaca file .integrons (TSV) dari IntegronFinder.
    Kolom utama: ID_replicon, element, pos_beg, pos_end, strand,
                 evalue, type_elt, annotation, model, type (integron class)
    """
    try:
        # IntegronFinder kadang menyertakan baris komentar (#)
        df = pd.read_csv(filepath, sep="\t", comment="#")
        if df.empty:
            return pd.DataFrame()

        # Kolom wajib ada
        if "ID_replicon" not in df.columns:
            print(f"  [WARN] Kolom ID_replicon tidak ditemukan di {filepath}", file=sys.stderr)
            return pd.DataFrame()

        df["sample_id"]    = sample_id
        df["country"]      = country
        df["country_name"] = country_name
        df["region"]       = region
        return df

    except Exception as e:
        print(f"  [WARN] Gagal membaca {filepath}: {e}", file=sys.stderr)
        return pd.DataFrame()


def classify_integron_elements(df):
    """
    IntegronFinder melaporkan SETIAP elemen dalam integron (bukan hanya integronnya).
    Kita perlu identifikasi struktur integron dari kolom 'type_elt' dan 'type':
    - type = integron class (complete, CALIN, In0)
    - type_elt = jenis elemen (intI = integrase, attC = recombination site, protein = gene cassette)
    """
    # Pastikan kolom lowercase untuk konsistensi
    if "type" in df.columns:
        df["integron_class"] = df["type"].str.lower().str.strip()
    else:
        df["integron_class"] = "unknown"

    if "type_elt" in df.columns:
        df["element_type"] = df["type_elt"].str.lower().str.strip()
    else:
        df["element_type"] = "unknown"

    return df


def extract_gene_cassettes(df):
    """
    Ekstraksi gene cassettes: hanya baris dengan type_elt = 'protein' atau 'attC'
    yang bukan integrase (annotation != 'intI*').
    Kolom 'annotation' berisi nama gen cassette.
    """
    if "annotation" not in df.columns:
        return pd.DataFrame()

    # Gene cassettes = elemen dengan annotation yang bukan integrase
    cassette_mask = (
        ~df["annotation"].str.lower().str.startswith("inti", na=True) &
        df["annotation"].notna() &
        (df["annotation"].str.strip() != "") &
        (df["annotation"].str.strip() != "NA")
    )
    cassettes = df[cassette_mask].copy()
    return cassettes


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Deep Integron Analysis — Cassette Composition & Class Distribution"
    )
    parser.add_argument("--base-dir",   required=True, help="Base results directory (e.g. results/pilot_run)")
    parser.add_argument("--output-dir", required=True, help="Output directory for this population")
    parser.add_argument("--pop",        required=True, help="Population/region name")
    parser.add_argument("--sample-map", default="data/metadata/pilot_map.csv")
    parser.add_argument("--coloc-csv",  required=True,
                        help="Path ke {pop}_all_coloc.csv untuk link ke AMR class")
    parser.add_argument("--flag-file",  required=True, help="Path ke flag file output")
    args = parser.parse_args()

    base_dir   = args.base_dir.rstrip("/")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(f"{output_dir}/figures", exist_ok=True)

    print("=" * 60)
    print(f"[07_integron_deep] Population: {args.pop}")
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

    # Filter sampel berdasarkan populasi
    if args.pop.lower() not in ["combined", "all"]:
        region_col = "Region" if "Region" in sample_map.columns else "region"
        sample_map = sample_map[sample_map[region_col] == args.pop]

    samples = sample_map["sample_id"].tolist()
    print(f"[INFO] Jumlah sampel yang diproses: {len(samples)}")

    # --- Kumpulkan semua data integron ---
    all_integrons = []
    for sample_id in samples:
        meta         = map_dict.get(sample_id, {})
        country      = meta.get("country", "Unknown")
        country_name = meta.get("country_name", "Unknown")
        region       = meta.get("Region", meta.get("region", "Unknown"))

        integron_file = f"{base_dir}/analysis/integron/{sample_id}.tsv"
        if not os.path.exists(integron_file):
            print(f"  [SKIP] File tidak ditemukan: {integron_file}")
            continue

        # Cek apakah file kosong (bisa terjadi jika tidak ada integron di sampel)
        if os.path.getsize(integron_file) == 0:
            print(f"  [SKIP] File kosong: {integron_file}")
            continue

        df = load_integron_file(integron_file, sample_id, country, country_name, region)
        if not df.empty:
            all_integrons.append(df)

    if not all_integrons:
        print("[WARN] Tidak ada data integron yang berhasil dikumpulkan.")
        # Buat file kosong agar Snakemake tidak fail
        for fname in ["complete_integrons.csv", "calin_cassettes.csv",
                      "cassette_matrix.csv", "class_distribution.csv",
                      "region_specificity.csv", "amr_integron_link.csv"]:
            pd.DataFrame().to_csv(f"{output_dir}/{fname}", index=False)
        open(args.flag_file, "w").close()
        return

    all_df = pd.concat(all_integrons, ignore_index=True)
    all_df = classify_integron_elements(all_df)
    print(f"[INFO] Total baris data integron: {len(all_df)}")

    # ============================================================
    # OUTPUT 1: Complete Integrons (Analisis Utama)
    # ============================================================
    print("\n[1] Memfilter COMPLETE integrons untuk analisis utama...")
    complete_df = all_df[all_df["integron_class"] == "complete"].copy()
    complete_df.to_csv(f"{output_dir}/complete_integrons.csv", index=False)
    print(f"    -> {len(complete_df)} baris complete integron tersimpan")

    # ============================================================
    # OUTPUT 2: CALIN Cassettes (Historical HGT Evidence)
    # ============================================================
    print("\n[2] Menyimpan CALIN cassettes (historical HGT evidence)...")
    calin_df = all_df[all_df["integron_class"] == "calin"].copy()
    calin_df.to_csv(f"{output_dir}/calin_cassettes.csv", index=False)
    print(f"    -> {len(calin_df)} baris CALIN cassette tersimpan (untuk Discussion section)")

    # ============================================================
    # OUTPUT 3: Integron Class Distribution per Country
    # ============================================================
    print("\n[3] Menghitung Integron Class Distribution per Country...")
    class_dist = (
        all_df
        .groupby(["country", "country_name", "region", "integron_class"])
        .agg(count=("ID_replicon", "count"),
             n_samples=("sample_id", "nunique"))
        .reset_index()
    )
    # Tambahkan proporsi per negara
    class_total = class_dist.groupby("country")["count"].transform("sum")
    class_dist["proportion"] = class_dist["count"] / class_total
    class_dist.to_csv(f"{output_dir}/class_distribution.csv", index=False)
    print(f"    -> {len(class_dist)} baris tersimpan di class_distribution.csv")

    # ============================================================
    # OUTPUT 4: Cassette Matrix (hanya dari complete integrons)
    # ============================================================
    print("\n[4] Membangun Gene Cassette Composition Matrix...")
    cassettes = extract_gene_cassettes(complete_df)
    if not cassettes.empty:
        # Hitung frekuensi setiap cassette per negara
        cassette_freq = (
            cassettes
            .groupby(["annotation", "country", "country_name", "region"])
            .agg(count=("sample_id", "count"),
                 n_samples=("sample_id", "nunique"))
            .reset_index()
        )
        cassette_freq.to_csv(f"{output_dir}/cassette_matrix.csv", index=False)
        print(f"    -> {len(cassette_freq)} baris cassette tersimpan")

        # ============================================================
        # OUTPUT 5: Region Specificity Analysis
        # ============================================================
        print("\n[5] Menganalisis Region Specificity cassettes...")
        specificity = (
            cassette_freq
            .groupby("annotation")
            .agg(
                n_regions=("region", "nunique"),
                regions_found=("region", lambda x: "|".join(sorted(x.unique()))),
                n_countries=("country", "nunique"),
                total_count=("count", "sum")
            )
            .reset_index()
        )
        # Klasifikasikan: universal = muncul di semua region, specific = hanya 1 region
        n_total_regions = all_df["region"].nunique()
        specificity["specificity"] = specificity["n_regions"].apply(
            lambda n: "universal" if n >= n_total_regions
                      else ("multi-regional" if n >= 2 else "region-specific")
        )
        specificity = specificity.sort_values("n_regions", ascending=False)
        specificity.to_csv(f"{output_dir}/region_specificity.csv", index=False)
        print(f"    -> {len(specificity)} unique cassettes dianalisis")
        if not specificity.empty:
            for label in ["universal", "multi-regional", "region-specific"]:
                n = len(specificity[specificity["specificity"] == label])
                print(f"       {label}: {n} cassettes")
    else:
        print("    [WARN] Tidak ada gene cassette yang ditemukan dari complete integrons.")
        pd.DataFrame().to_csv(f"{output_dir}/cassette_matrix.csv", index=False)
        pd.DataFrame().to_csv(f"{output_dir}/region_specificity.csv", index=False)

    # ============================================================
    # OUTPUT 6: AMR-Integron Link (via colocalization data)
    # ============================================================
    print("\n[6] Menghubungkan Integron dengan AMR class dari data colocalization...")
    try:
        coloc_df = pd.read_csv(args.coloc_csv)

        # Filter hanya AMR yang ada di integron (On_Integron == True)
        if "On_Integron" in coloc_df.columns:
            integron_amr = coloc_df[coloc_df["On_Integron"] == True].copy()
            integron_amr = integron_amr.filter(
                items=["Sample_ID", "Country", "Country_Name", "Region",
                       "Best_Hit_ARO", "Drug Class", "AMR Gene Family",
                       "Best_Identities", "Contig"]
            )
            integron_amr.to_csv(f"{output_dir}/amr_integron_link.csv", index=False)
            print(f"    -> {len(integron_amr)} AMR genes on integrons tersimpan")

            # Summary per drug class
            if not integron_amr.empty and "Drug Class" in integron_amr.columns:
                summary = (
                    integron_amr
                    .groupby(["Drug Class", "Region" if "Region" in integron_amr.columns else "Country"])
                    .size()
                    .reset_index(name="count")
                    .sort_values("count", ascending=False)
                )
                print(f"       Top AMR classes on integrons:")
                for _, row in summary.head(5).iterrows():
                    print(f"         {row['Drug Class']}: {row['count']}")
        else:
            print("    [WARN] Kolom On_Integron tidak ditemukan di colocalization data.")
            pd.DataFrame().to_csv(f"{output_dir}/amr_integron_link.csv", index=False)
    except Exception as e:
        print(f"    [WARN] Gagal membaca colocalization data: {e}", file=sys.stderr)
        pd.DataFrame().to_csv(f"{output_dir}/amr_integron_link.csv", index=False)

    # --- Done ---
    print("\n" + "=" * 60)
    print("[SELESAI] Deep Integron Analysis complete.")
    print(f"  Output tersimpan di: {output_dir}")
    print("=" * 60)

    open(args.flag_file, "w").close()
    print(f"\n[DONE] Flag file dibuat: {args.flag_file}")


if __name__ == "__main__":
    main()
