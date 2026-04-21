#!/usr/bin/env python3
"""
Skrip: Mencari co-localization AMR dan MGE
Fungsi: Membaca output RGI dan MOB-suite, lalu menyebarkan gen AMR mana saja yang terletak di dalam Plasmid.
"""

import pandas as pd
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, help="File CSV output")
    args = parser.parse_args()
    
    rgi_files = glob.glob("results/rgi/*_rgi.txt")
    all_colocalized = []
    
    print("Memulai pemindaian persilangan antar gen...")
    
    for rgi_file in rgi_files:
        sample_id = os.path.basename(rgi_file).replace('_rgi.txt', '')
        mob_file = f"results/mobsuite/{sample_id}_plasmid.txt"
        
        if not os.path.exists(mob_file):
            continue
            
        try:
            rgi_df = pd.read_csv(rgi_file, sep="\t")
            mob_df = pd.read_csv(mob_file, sep="\t")
        except:
            continue
            
        if rgi_df.empty or mob_df.empty:
            continue
            
        # Bersihkan nama contig RGI
        rgi_df['Contig_Bersih'] = rgi_df['Contig'].apply(lambda x: "_".join(x.split('_')[:-1]) if '_' in str(x) else x)
        
        # Ambil daftar ID contig Plasmid
        plasmids_df = mob_df[mob_df['molecule_type'] == 'plasmid']
        plasmid_contigs = set(plasmids_df['contig_id'].astype(str))
        
        # Cek persilangan
        rgi_df['On_Plasmid'] = rgi_df['Contig_Bersih'].astype(str).isin(plasmid_contigs)
        rgi_df['Sample_ID'] = sample_id
        
        all_colocalized.append(rgi_df)
        
    if all_colocalized:
        final_df = pd.concat(all_colocalized, ignore_index=True)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        final_df.to_csv(args.out, index=False)
        print(f"Sukses! Data silang berhasil disimpan ke {args.out}")
    else:
        print("Tidak ada sampel yang berhasil disilangkan.")
        pd.DataFrame().to_csv(args.out, index=False)

if __name__ == "__main__":
    main()
