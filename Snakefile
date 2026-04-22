# ============================================================
# Load konfigurasi dari config.yaml
# ============================================================
configfile: "config.yaml"

"""
Snakefile: AMR Gene Mobility Pipeline
Versi: 3.1 (V3 - Syntax Fixed)

Cara penggunaan:
  HPC / Lokal (Phase 3 saja):
    snakemake phase3 --cores all --use-conda

  Statistik & Network (Phase 4 & 5):
    snakemake phase4 --cores 4 --use-conda
    snakemake phase5 --cores 4 --use-conda

  Semua sekaligus:
    snakemake all --cores all --use-conda

  Dry-run (validasi pipeline):
    snakemake phase3 --dry-run --cores 4
"""

import pandas as pd
import os
import glob as glob_mod

SAMPLE_MAP_FILE = config["paths"]["sample_map"]

if os.path.exists(SAMPLE_MAP_FILE):
    sample_df = pd.read_csv(SAMPLE_MAP_FILE)
    SAMPLES = sample_df["sample_id"].tolist()
else:
    # Fallback: deteksi otomatis dari file .fasta yang sudah ada
    SAMPLES = [os.path.basename(f).replace(".fasta", "")
               for f in glob_mod.glob("data/contigs/*.fasta")]

print(f"[Snakemake] Sampel terdeteksi: {SAMPLES}")

# ============================================================
# TARGET RULES (Entrypoints)
# ============================================================

rule all:
    """Target utama: jalankan seluruh pipeline end-to-end."""
    input:
        expand("results/rgi/{sample}_rgi.txt",            sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt",   sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv",         sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv",
        "results/figures/Fig2_MGE_distribution.pdf",
        "results/figures/Fig3_Network_combined_done.flag"

rule phase3:
    """Jalankan hanya deteksi AMR + MGE + integrasi data."""
    input:
        expand("results/rgi/{sample}_rgi.txt",            sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt",   sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv",         sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv"

rule phase4:
    """Jalankan statistik dan visualisasi (setelah Phase 3 selesai)."""
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
        rgi = "results/rgi/{sample}_rgi.txt"
    params:
        input_type     = config["rgi"]["input_type"],
        alignment_tool = config["rgi"]["alignment_tool"]
    conda:
        "envs/rgi.yaml"
    log:
        "logs/rgi_{sample}.log"
    threads: config["resources"]["threads_rgi"]
    shell:
        """
        echo "[RGI] Memproses {wildcards.sample}..." | tee {log}

        # Auto-load CARD database dari conda environment
        CARD_JSON=$(ls ${{CONDA_PREFIX}}/lib/python*/site-packages/app/_data/card.json 2>/dev/null | head -1)
        if [ -n "$CARD_JSON" ]; then
            echo "[RGI] Loading CARD database: $CARD_JSON" | tee -a {log}
            rgi load --card_json "$CARD_JSON" >> {log} 2>&1
        else
            echo "[RGI] Downloading CARD database..." | tee -a {log}
            rgi database --download >> {log} 2>&1
            rgi load --card_json card.json >> {log} 2>&1
        fi

        rgi main \
            --input_sequence {input.fasta} \
            --output_file results/rgi/{wildcards.sample}_rgi \
            --input_type {params.input_type} \
            --alignment_tool {params.alignment_tool} \
            --clean \
            --num_threads {threads} >> {log} 2>&1

        if [ ! -f {output.rgi} ]; then
            touch {output.rgi}
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
    threads: config["resources"]["threads_mobsuite"]
    shell:
        """
        echo "[MOBsuite] Memproses {wildcards.sample}..." | tee {log}
        mob_recon \
            --infile {input.fasta} \
            --outdir results/mobsuite/{wildcards.sample}_dir \
            --force \
            --num_threads {threads} >> {log} 2>&1

        cp results/mobsuite/{wildcards.sample}_dir/contig_report.txt {output.report}
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
        tsv = "results/integron/{sample}_integrons.tsv"
    params:
        mode = config["integronfinder"]["mode"]
    conda:
        "envs/integronfinder.yaml"
    log:
        "logs/integron_{sample}.log"
    threads: config["resources"]["threads_integron"]
    shell:
        """
        echo "[IntegronFinder] Memproses {wildcards.sample}..." | tee {log}
        mkdir -p results/integron/{wildcards.sample}_dir

        integron_finder \
            {input.fasta} \
            --outdir results/integron/{wildcards.sample}_dir \
            --cpu {threads} \
            --{params.mode} >> {log} 2>&1

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
        tsv = "results/isescan/{sample}_is.tsv"
    conda:
        "envs/isescan.yaml"
    log:
        "logs/isescan_{sample}.log"
    threads: config["resources"]["threads_isescan"]
    shell:
        """
        echo "[ISEScan] Memproses {wildcards.sample}..." | tee {log}
        mkdir -p results/isescan/{wildcards.sample}_dir

        isescan.py \
            --seqfile {input.fasta} \
            --output results/isescan/{wildcards.sample}_dir \
            --nthread {threads} >> {log} 2>&1

        find results/isescan/{wildcards.sample}_dir -name "*.tsv" \
            | head -1 | xargs -I_F_ cp _F_ {output.tsv} 2>>{log} \
            || touch {output.tsv}
        echo "[ISEScan] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 3E: HGT Linkage — Co-localization Integration
# ============================================================

rule colocalization:
    """
    Mengintegrasikan hasil keempat tools untuk setiap gen AMR.
    Menentukan posisi: Plasmid / Integron / IS / Kromosom.
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
    Merangkum colocalization_summary.csv menjadi tiga matrix agregat.
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
    - 4A: Alpha & beta diversity (Shannon, PERMANOVA)   -> Sub-Q1
    - 4B: Kruskal-Wallis per AMR class antar populasi  -> Sub-Q1
    - 4C: Chi-square distribusi MGE type               -> Sub-Q2
    - 4D: Fisher's exact AMR x MGE (FDR-corrected)    -> Sub-Q3
    - 4E: Heatmap + stacked bar + supplementary tables
    """
    input:
        abund = "results/amr_abundance_matrix.csv",
        mge   = "results/mge_distribution_matrix.csv",
        assoc = "results/amr_mge_association_matrix.csv",
        coloc = "results/colocalization_summary.csv"
    output:
        fig2 = "results/figures/Fig2_MGE_distribution.pdf"
    conda:
        "envs/r_stats.yaml"
    log:
        "logs/statistics.log"
    threads: config["resources"]["threads_stats"]
    shell:
        """
        echo "[Statistics] Menjalankan 04_run_stats.R..." | tee {log}
        Rscript scripts/04_run_stats.R 2>&1 | tee -a {log}
        echo "[Statistics] Selesai." | tee -a {log}
        """

# ============================================================
# PHASE 5: Network Comparison (Opsional - novelty figure)
# ============================================================

rule network_analysis:
    """
    Analisis jaringan bipartit AMR-MGE per populasi.
    Membandingkan struktur network antar DNK, CHN, IND.
    Output: Jaccard similarity + hub genes + network figures.
    """
    input:
        coloc = "results/colocalization_summary.csv"
    output:
        flag = "results/figures/Fig3_Network_combined_done.flag"
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