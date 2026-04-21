"""
Snakefile: AMR Gene Mobility Pipeline
Versi: 3.0 (V3 - Full Comprehensive Pipeline)

Cara penggunaan:
  HPC / Lokal (Phase 3 saja):
    snakemake phase3 --cores all --use-conda

  Statistik & Network (Phase 4 & 5, jalankan setelah Phase 3 selesai):
    snakemake phase4 --cores 4 --use-conda
    snakemake phase5 --cores 4 --use-conda

  Semua sekaligus:
    snakemake all --cores all --use-conda
"""

import pandas as pd
import os

# ============================================================
# KONFIGURASI: Baca sample list dari sample_map.csv
# ============================================================
SAMPLE_MAP_FILE = "data/metadata/sample_map.csv"

if os.path.exists(SAMPLE_MAP_FILE):
    sample_df = pd.read_csv(SAMPLE_MAP_FILE)
    SAMPLES = sample_df["sample_id"].tolist()
else:
    # Fallback: deteksi otomatis dari file .fasta yang sudah ada
    import glob
    SAMPLES = [os.path.basename(f).replace(".fasta", "")
               for f in glob.glob("data/contigs/*.fasta")]

print(f"[Snakemake] Sampel terdeteksi: {SAMPLES}")

# ============================================================
# TARGET RULES (Entrypoint)
# ============================================================

rule all:
    """Target utama: jalankan seluruh pipeline end-to-end."""
    input:
        expand("results/rgi/{sample}_rgi.txt",         sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv",      sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv",
        "results/figures/Fig2_MGE_distribution.pdf",
        "results/figures/Fig3_Network_combined_done.flag"

rule phase3:
    """Jalankan hanya deteksi AMR + MGE + integrasi data."""
    input:
        expand("results/rgi/{sample}_rgi.txt",           sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt",  sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv",sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv",        sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv"

rule phase4:
    """Jalankan statistik dan visualisasi saja (setelah Phase 3 selesai)."""
    input:
        "results/figures/Fig2_MGE_distribution.pdf"

rule phase5:
    """Jalankan network comparison (opsional)."""
    input:
        "results/figures/Fig3_Network_combined_done.flag"

# ============================================================
# PHASE 3A: Deteksi AMR dengan RGI + CARD Database
# ============================================================

rule run_rgi:
    """
    Deteksi gen resistensi antibiotik dari contigs menggunakan RGI + CARD database.
    Output: tabel hit AMR per contig.
    """
    input:
        fasta = "data/contigs/{sample}.fasta"
    output:
        rgi   = "results/rgi/{sample}_rgi.txt"
    conda:
        "envs/rgi.yaml"
    log:
        "logs/rgi_{sample}.log"
    threads: 4
    shell:
        """
        echo "[RGI] Memproses {wildcards.sample}..." | tee {log}
        rgi main \
            --input_sequence {input.fasta} \
            --output_file results/rgi/{wildcards.sample}_rgi \
            --input_type contig \
            --alignment_tool BLAST \
            --clean \
            --num_threads {threads} >> {log} 2>&1

        # RGI menambahkan '.txt' otomatis ke nama output
        if [ ! -f {output.rgi} ]; then
            mv results/rgi/{wildcards.sample}_rgi.txt {output.rgi}
        fi
        echo "[RGI] Selesai. $(grep -c '' {output.rgi}) baris." | tee -a {log}
        """

# ============================================================
# PHASE 3B: Deteksi Plasmid dengan MOB-suite
# ============================================================

rule run_mobsuite:
    """
    Identifikasi contig yang berasal dari plasmid dan tipe mobilisasinya.
    Output: laporan per contig dengan kolom molecule_type dan mobility.
    """
    input:
        fasta = "data/contigs/{sample}.fasta"
    output:
        report = "results/mobsuite/{sample}_plasmid.txt"
    conda:
        "envs/mobsuite.yaml"
    log:
        "logs/mobsuite_{sample}.log"
    threads: 4
    shell:
        """
        echo "[MOBsuite] Memproses {wildcards.sample}..." | tee {log}
        mob_recon \
            --infile {input.fasta} \
            --outdir results/mobsuite/{wildcards.sample}_dir \
            --num_threads {threads} >> {log} 2>&1

        # Pindahkan file laporan ke path output
        cp results/mobsuite/{wildcards.sample}_dir/mob_recon_report.txt {output.report}
        echo "[MOBsuite] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 3C: Deteksi Integron dengan IntegronFinder
# ============================================================

rule run_integronfinder:
    """
    Deteksi integron dan gene cassette pada contigs.
    Output: TSV dengan koordinat intI, attC, dan gene cassette.
    """
    input:
        fasta = "data/contigs/{sample}.fasta"
    output:
        tsv   = "results/integron/{sample}_integrons.tsv"
    conda:
        "envs/integronfinder.yaml"
    log:
        "logs/integron_{sample}.log"
    threads: 2
    shell:
        """
        echo "[IntegronFinder] Memproses {wildcards.sample}..." | tee {log}
        mkdir -p results/integron/{wildcards.sample}_dir

        integron_finder \
            {input.fasta} \
            --outdir results/integron/{wildcards.sample}_dir \
            --cpu {threads} \
            --local-max >> {log} 2>&1

        # Gabungkan semua file .integrons menjadi satu TSV
        find results/integron/{wildcards.sample}_dir -name "*.integrons" \
            | xargs cat > {output.tsv} 2>>{log} || touch {output.tsv}
        echo "[IntegronFinder] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 3D: Deteksi IS / Transposon dengan ISEScan
# ============================================================

rule run_isescan:
    """
    Deteksi Insertion Sequences (IS) dan transposon dari contigs.
    Output: TSV dengan koordinat IS per contig.
    """
    input:
        fasta = "data/contigs/{sample}.fasta"
    output:
        tsv   = "results/isescan/{sample}_is.tsv"
    conda:
        "envs/isescan.yaml"
    log:
        "logs/isescan_{sample}.log"
    threads: 2
    shell:
        """
        echo "[ISEScan] Memproses {wildcards.sample}..." | tee {log}
        mkdir -p results/isescan/{wildcards.sample}_dir

        isescan.py \
            --seqfile {input.fasta} \
            --output results/isescan/{wildcards.sample}_dir \
            --nthread {threads} >> {log} 2>&1

        # Ambil file summary IS
        find results/isescan/{wildcards.sample}_dir -name "*.tsv" \
            | head -1 | xargs cp -t results/isescan/ 2>>{log} \
            && mv results/isescan/*.tsv {output.tsv} 2>>{log} \
            || touch {output.tsv}
        echo "[ISEScan] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 3E: HGT Linkage — Co-localization Integration
# ============================================================

rule colocalization:
    """
    Mengintegrasikan hasil keempat tools di atas untuk setiap gen AMR.
    Menentukan apakah gen AMR berada di Plasmid, Integron, IS, atau Kromosom.
    Menambahkan metadata Country dan Region dari sample_map.csv.
    """
    input:
        rgi      = expand("results/rgi/{sample}_rgi.txt",            sample=SAMPLES),
        mob      = expand("results/mobsuite/{sample}_plasmid.txt",   sample=SAMPLES),
        integron = expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        isescan  = expand("results/isescan/{sample}_is.tsv",         sample=SAMPLES),
        map      = SAMPLE_MAP_FILE
    output:
        csv = "results/colocalization_summary.csv"
    conda:
        "envs/python.yaml"
    log:
        "logs/colocalization.log"
    shell:
        """
        echo "[Co-localization] Mengintegrasikan data MGE dan AMR..." | tee {log}
        python scripts/02_find_colocalization.py \
            --out {output.csv} \
            --sample-map {input.map} \
            2>&1 | tee -a {log}
        echo "[Co-localization] Selesai. Output: {output.csv}" | tee -a {log}
        """

# ============================================================
# PHASE 3F: Agregasi per Populasi
# ============================================================

rule aggregate_by_population:
    """
    Merangkum colocalization_summary.csv menjadi tiga matrix agregat
    yang siap digunakan untuk analisis statistik.
    """
    input:
        csv = "results/colocalization_summary.csv"
    output:
        abund = "results/amr_abundance_matrix.csv",
        mge   = "results/mge_distribution_matrix.csv",
        assoc = "results/amr_mge_association_matrix.csv"
    conda:
        "envs/python.yaml"
    log:
        "logs/aggregate.log"
    shell:
        """
        echo "[Aggregate] Merangkum data per populasi..." | tee {log}
        python scripts/03_aggregate_by_population.py 2>&1 | tee -a {log}
        """

# ============================================================
# PHASE 4: Statistik & Visualisasi
# ============================================================

rule run_statistics:
    """
    Analisis statistik komprehensif:
    - 4A: Alpha & beta diversity (Shannon, PERMANOVA)
    - 4B: Kruskal-Wallis per AMR class antar populasi
    - 4C: Chi-square distribusi MGE type
    - 4D: Fisher's exact AMR × MGE association (FDR-corrected)
    - 4E: Heatmap, stacked bar, supplementary tables
    """
    input:
        abund = "results/amr_abundance_matrix.csv",
        mge   = "results/mge_distribution_matrix.csv",
        assoc = "results/amr_mge_association_matrix.csv",
        coloc = "results/colocalization_summary.csv"
    output:
        fig2  = "results/figures/Fig2_MGE_distribution.pdf"
    conda:
        "envs/r_stats.yaml"
    log:
        "logs/statistics.log"
    shell:
        """
        echo "[Statistics] Menjalankan 04_run_stats.R..." | tee {log}
        Rscript scripts/04_run_stats.R 2>&1 | tee -a {log}
        echo "[Statistics] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 5: Network Comparison (Opsional)
# ============================================================

rule network_analysis:
    """
    Analisis jaringan bipartit AMR-MGE per populasi.
    Membandingkan struktur network antar DNK, CHN, IND.
    Menghasilkan Jaccard similarity dan hub genes.
    """
    input:
        coloc = "results/colocalization_summary.csv"
    output:
        flag  = "results/figures/Fig3_Network_combined_done.flag"
    conda:
        "envs/r_stats.yaml"
    log:
        "logs/network.log"
    shell:
        """
        echo "[Network] Menjalankan 05_network_analysis.R..." | tee {log}
        Rscript scripts/05_network_analysis.R 2>&1 | tee -a {log}
        touch {output.flag}
        echo "[Network] Selesai." | tee -a {log}
        """