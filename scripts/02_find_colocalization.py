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

def first_token(x):
    """Return the first whitespace-delimited token of a contig id.
    MEGAHIT headers look like 'k99_10000 flag=0 multi=3.9 len=1560'. RGI (prodigal)
    and IntegronFinder keep only 'k99_10000', but MOB-suite stores the WHOLE header
    as contig_id. Normalising every tool to the first token is what makes the
    co-localization join line up (otherwise plasmid/IS hits never match)."""
    s = str(x)
    parts = s.split()
    return parts[0] if parts else s


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
        # MOB-suite stores the full FASTA header in contig_id; normalise to first token.
        plasmid_contigs = set(plasmid_df["contig_id"].astype(str).map(first_token))
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
        integron_contigs = set(df["ID_replicon"].astype(str).map(first_token))
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
        is_contigs = set(df["seqid"].astype(str).map(first_token))
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

    # Columns kept in the colocalization output (defined once, reused for the empty case).
    keep_cols = [
        "Sample_ID", "Country", "Country_Name", "Region", "Age", "Gender", "Study",
        "Contig", "contig_clean", "Best_Hit_ARO",
        "Drug Class", "Resistance Mechanism", "AMR Gene Family",
        "Best_Identities", "On_Plasmid", "On_Integron", "On_IS", "MGE_Type"
    ]

    # --- Load RGI ---
    rgi_df = load_rgi(rgi_file)
    if rgi_df.empty:
        # RGI ran but found no AMR genes (or the file is header-only). This is a valid
        # biological outcome, not an error: emit an empty colocalization so aggregation
        # for the whole population is not blocked by one AMR-free sample.
        print(f"    [INFO] Tidak ada gen AMR di {sample_id} (RGI kosong). Menulis colocalization kosong.")
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        pd.DataFrame(columns=keep_cols).to_csv(args.out, index=False)
        sys.exit(0)

    # --- Load MGEs (Opsional, kalau tidak ada kita anggap kosong) ---
    mob_contigs, _ = load_mobsuite(mob_file) if os.path.exists(mob_file) else (set(), None)
    integron_contigs, _ = load_integronfinder(integron_file) if os.path.exists(integron_file) else (set(), None)
    is_contigs, _ = load_isescan(is_file) if os.path.exists(is_file) else (set(), None)

    # --- Anotasi tiap gen AMR ---
    rgi_df["On_Plasmid"]  = rgi_df["contig_clean"].isin(mob_contigs)
    rgi_df["On_Integron"] = rgi_df["contig_clean"].isin(integron_contigs)
    rgi_df["On_IS"]       = rgi_df["contig_clean"].isin(is_contigs)
    rgi_df["MGE_Type"]    = rgi_df.apply(determine_mge_type, axis=1)

    # --- QC: contig-namespace match diagnostic ---
    # The whole co-localization result depends on AMR contig-ids matching the contig-ids
    # each MGE tool reports. If a tool finds contigs but NONE overlap the AMR namespace,
    # every gene silently falls back to "Chromosomal" — this surfaces that failure mode.
    amr_contigs = set(rgi_df["contig_clean"])
    n_amr = len(amr_contigs)
    qc_rows = []
    print(f"  [QC] Contig-namespace match for {sample_id} (n_amr_contigs={n_amr}):")
    for name, mge_set in [("Plasmid", mob_contigs), ("Integron", integron_contigs), ("IS", is_contigs)]:
        n_tool = len(mge_set)
        n_hit  = len(amr_contigs & mge_set)
        pct    = (100.0 * n_hit / n_amr) if n_amr else 0.0
        warn   = (n_tool > 0 and n_hit == 0)
        flag   = "  <-- WARNING: tool found contigs but none match AMR namespace" if warn else ""
        print(f"        {name}: {n_hit}/{n_amr} matched ({pct:.1f}%); tool_contigs={n_tool}{flag}")
        qc_rows.append({"sample_id": sample_id, "mge_tool": name,
                        "n_amr_contigs": n_amr, "n_matched": n_hit,
                        "pct_matched": round(pct, 2), "tool_contigs": n_tool,
                        "namespace_mismatch_warning": warn})
    try:
        qc_path = os.path.join(os.path.dirname(args.out), f"{sample_id}_matchqc.csv")
        pd.DataFrame(qc_rows).to_csv(qc_path, index=False)
    except Exception as e:
        print(f"  [WARN] Gagal menulis QC match file: {e}", file=sys.stderr)

    # --- Tambah metadata ---
    rgi_df["Sample_ID"]    = sample_id
    rgi_df["Country"]      = country
    rgi_df["Country_Name"] = country_name
    rgi_df["Region"]       = region
    rgi_df["Age"]          = age
    rgi_df["Gender"]       = gender
    rgi_df["Study"]        = study_name

    # Keep only the declared columns (keep_cols defined above).
    available = [c for c in keep_cols if c in rgi_df.columns]
    final_df = rgi_df[available]

    # --- Simpan ---
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    final_df.to_csv(args.out, index=False)
    print(f"  [DONE] Colocalization saved: {args.out} ({len(final_df)} genes)")


if __name__ == "__main__":
    main()
