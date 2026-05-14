configfile: "config.yaml"

"""
Snakefile: AMR Gene Mobility Pipeline
Versi: 4.0 (End-to-End)

Pipeline lengkap dari accession ID hingga analisis AMR-MGE.

Cara penggunaan:
  End-to-end (250 sampel di HPC):
    snakemake all --cores all --use-conda

  Dari contigs saja (backward compatible):
    snakemake phase3 --cores all --use-conda

  Statistik & Network:
    snakemake phase4 --cores 4 --use-conda
    snakemake phase5 --cores 4 --use-conda

  Dry-run:
    snakemake all --dry-run --cores 4
"""

import pandas as pd
import os
import glob as glob_mod

SAMPLE_MAP_FILE = config["paths"]["sample_map"]

if os.path.exists(SAMPLE_MAP_FILE):
    sample_df = pd.read_csv(SAMPLE_MAP_FILE)
    SAMPLES = sample_df["sample_id"].tolist()
    ACCESSION_MAP = dict(zip(sample_df["sample_id"], sample_df["accession"]))
else:
    SAMPLES = [os.path.basename(f).replace(".fasta", "")
               for f in glob_mod.glob("data/contigs/*.fasta")]
    ACCESSION_MAP = {}

print(f"[Snakemake] Sampel terdeteksi: {SAMPLES}")
print(f"[Snakemake] Total: {len(SAMPLES)} sampel")

def get_input_type(wildcards):
    """
    Detect input type for sample.
    Returns: 'contig' if FASTA exists, 'fastq' if FASTQ exists, 'accession' otherwise
    """
    fasta_path = f"data/contigs/{wildcards.sample}.fasta"
    fastq_path = f"data/raw_reads/{wildcards.sample}_1.fastq.gz"
    
    if os.path.exists(fasta_path):
        return "contig"
    elif os.path.exists(fastq_path):
        return "fastq"
    else:
        return "accession"

def get_assembly_input(wildcards):
    """
    Determine input for assembly rule based on available files.
    """
    input_type = get_input_type(wildcards)
    
    if input_type == "contig":
        return f"data/contigs/{wildcards.sample}.fasta"
    else:
        return {
            "r1": f"data/nonhost_reads/{wildcards.sample}_1.fastq.gz",
            "r2": f"data/nonhost_reads/{wildcards.sample}_2.fastq.gz"
        }

rule all:
    """Target utama: jalankan seluruh pipeline end-to-end.
    Termasuk cleanup intermediate FASTQ setelah assembly berhasil
    untuk menghemat kuota storage (500 GB scratch Mahameru).
    """
    input:
        expand("results/rgi/{sample}_rgi.txt", sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv", sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv",
        "results/figures/Fig2_MGE_distribution.pdf",
        "results/figures/Fig3_Network_combined_done.flag",
        expand("logs/cleanup_{sample}.done", sample=SAMPLES)

rule all_full:
    """Target untuk full pipeline dari accession (Phase 1-5)."""
    input:
        expand("data/contigs/{sample}.fasta", sample=SAMPLES),
        expand("results/rgi/{sample}_rgi.txt", sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv", sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv",
        "results/figures/Fig2_MGE_distribution.pdf",
        "results/figures/Fig3_Network_combined_done.flag"

rule phase3:
    """Jalankan hanya deteksi AMR + MGE + integrasi data."""
    input:
        expand("results/rgi/{sample}_rgi.txt", sample=SAMPLES),
        expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES),
        expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        expand("results/isescan/{sample}_is.tsv", sample=SAMPLES),
        "results/colocalization_summary.csv",
        "results/amr_abundance_matrix.csv"

rule phase4:
    """Jalankan statistik dan visualisasi."""
    input:
        "results/figures/Fig2_MGE_distribution.pdf"

rule phase5:
    """Jalankan network comparison."""
    input:
        "results/figures/Fig3_Network_combined_done.flag"

# PHASE 1: Download SRA Data

rule download_sra:
    """
    Download raw reads dari NCBI SRA.
    Menggunakan sra-tools dengan fallback mechanism.
    """
    input:
        sample_map = SAMPLE_MAP_FILE
    output:
        r1 = "data/raw_reads/{sample}_1.fastq.gz",
        r2 = "data/raw_reads/{sample}_2.fastq.gz"
    params:
        accession = lambda wildcards: ACCESSION_MAP.get(wildcards.sample, wildcards.sample)
    conda:
        "envs/sra-tools.yaml"
    log:
        "logs/download_{sample}.log"
    benchmark:
        "benchmarks/download_{sample}.tsv"
    threads: config["resources"]["threads_sra"]
    resources:
        mem_mb   = 4000,
        disk_mb  = 15000,
        partition = "short",
        runtime  = "01:00:00"
    shell:
        """
        echo "[Download] Mengunduh {wildcards.sample} ({params.accession})..." | tee {log}
        
        mkdir -p data/raw_reads tmp_download_{wildcards.sample}
        
        # Split accession by semicolon
        IFS=';' read -ra ADDR <<< "{params.accession}"
        
        for acc in "${{ADDR[@]}}"; do
            echo "[Download] Memproses run: $acc" >> {log}
            
            # Prefetch
            prefetch "$acc" --max-size 50G 2>>{log} || \
                echo "[WARN] prefetch $acc gagal, mencoba fasterq-dump langsung..." >>{log}
            
            # Download FASTQ
            fasterq-dump --split-files --threads {threads} "$acc" \
                --outdir tmp_download_{wildcards.sample} 2>>{log}
            
            # Gabungkan (append) ke file temporer
            cat tmp_download_{wildcards.sample}/"$acc"_1.fastq >> tmp_combined_{wildcards.sample}_1.fastq
            cat tmp_download_{wildcards.sample}/"$acc"_2.fastq >> tmp_combined_{wildcards.sample}_2.fastq
            
            # Bersihkan temporary per-run
            rm -rf "$acc" tmp_download_{wildcards.sample}/"$acc"*
        done
        
        # Kompresi file gabungan final
        echo "[Download] Kompresi file gabungan..." >> {log}
        pigz -p {threads} -c tmp_combined_{wildcards.sample}_1.fastq > {output.r1}
        pigz -p {threads} -c tmp_combined_{wildcards.sample}_2.fastq > {output.r2}
        
        # Bersihkan temporary gabungan
        rm tmp_combined_{wildcards.sample}_1.fastq tmp_combined_{wildcards.sample}_2.fastq
        rm -rf tmp_download_{wildcards.sample}
        
        echo "[Download] Selesai. Ukuran file gabungan:" >> {log}
        ls -lh {output.r1} {output.r2} >> {log}
        """

# PHASE 2A: Quality Control dengan fastp

rule qc_fastp:
    """
    Quality control dan adapter trimming menggunakan fastp.
    Menghasilkan cleaned reads dan QC report.
    """
    input:
        r1 = "data/raw_reads/{sample}_1.fastq.gz",
        r2 = "data/raw_reads/{sample}_2.fastq.gz"
    output:
        qc1 = "data/qc_reads/{sample}_1.fastq.gz",
        qc2 = "data/qc_reads/{sample}_2.fastq.gz",
        html = "logs/qc/fastp_{sample}.html",
        json = "logs/qc/fastp_{sample}.json"
    params:
        qualified_quality = config.get("fastp", {}).get("qualified_quality", 20),
        length_required = config.get("fastp", {}).get("length_required", 50)
    conda:
        "envs/fastp.yaml"
    log:
        "logs/fastp_{sample}.log"
    benchmark:
        "benchmarks/fastp_{sample}.tsv"
    threads: config["resources"]["threads_fastp"]
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "00:30:00"
    shell:
        """
        echo "[Fastp] Quality control untuk {wildcards.sample}..." | tee {log}
        
        mkdir -p data/qc_reads logs/qc
        
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.qc1} -O {output.qc2} \
            --thread {threads} \
            --html {output.html} \
            --json {output.json} \
            --qualified_quality_phred {params.qualified_quality} \
            --length_required {params.length_required} \
            --detect_adapter_for_pe \
            --correction \
            2>>{log}
        
        echo "[Fastp] Selesai. QC report: {output.html}" >>{log}
        """

# PHASE 2B: Host Removal dengan Bowtie2

rule download_hg38_index:
    """
    Download dan index hg38 reference genome untuk host removal.
    Hanya dijalankan sekali, output diproteksi.
    """
    output:
        index = protected("databases/hg38/hg38.1.bt2")
    params:
        url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    conda:
        "envs/bowtie2.yaml"
    log:
        "logs/download_hg38.log"
    benchmark:
        "benchmarks/download_hg38.tsv"
    threads: 4
    resources:
        mem_mb   = 16000,
        partition = "medium-small",
        runtime  = "04:00:00"
    shell:
        """
        echo "[HG38] Downloading human reference genome..." | tee {log}
        
        mkdir -p databases/hg38
        cd databases/hg38
        
        wget -O hg38.fa.gz {params.url} 2>>{log} || \
            curl -L -o hg38.fa.gz {params.url} 2>>{log}
        
        gunzip hg38.fa.gz 2>>{log}
        
        bowtie2-build --threads {threads} hg38.fa hg38 2>>{log}
        
        echo "[HG38] Index build complete." >>{log}
        """

rule host_removal:
    """
    Remove human reads menggunakan Bowtie2 alignment ke hg38.
    Output: non-host reads yang siap untuk assembly.
    """
    input:
        r1 = "data/qc_reads/{sample}_1.fastq.gz",
        r2 = "data/qc_reads/{sample}_2.fastq.gz",
        index = "databases/hg38/hg38.1.bt2"
    output:
        nh1 = "data/nonhost_reads/{sample}_1.fastq.gz",
        nh2 = "data/nonhost_reads/{sample}_2.fastq.gz",
        stats = "logs/host_removal/{sample}_bowtie2.stats"
    params:
        index_base = "databases/hg38/hg38"
    conda:
        "envs/bowtie2.yaml"
    log:
        "logs/host_removal_{sample}.log"
    benchmark:
        "benchmarks/host_removal_{sample}.tsv"
    threads: config["resources"]["threads_bowtie"]
    resources:
        mem_mb   = 16000,
        partition = "short",
        runtime  = "02:00:00"
    shell:
        """
        echo "[Host Removal] {wildcards.sample} - removing human reads..." | tee {log}
        
        mkdir -p data/nonhost_reads logs/host_removal
        
        bowtie2 \
            -x {params.index_base} \
            -1 {input.r1} -2 {input.r2} \
            --un-conc-gz data/nonhost_reads/{wildcards.sample}_%.fastq.gz \
            --very-sensitive \
            -p {threads} \
            2> {output.stats} 1>>{log}
        
        echo "[Host Removal] Selesai." >>{log}
        
        grep "overall alignment rate" {output.stats} >>{log} || true
        """

# PHASE 2C: De Novo Assembly dengan MEGAHIT

rule assembly_megahit:
    """
    Metagenome assembly menggunakan MEGAHIT.
    Menghasilkan contigs dari non-host reads.
    """
    input:
        r1 = "data/nonhost_reads/{sample}_1.fastq.gz",
        r2 = "data/nonhost_reads/{sample}_2.fastq.gz"
    output:
        contigs = "data/contigs/{sample}.fasta",
        log_file = "logs/assembly/{sample}_megahit.log"
    params:
        preset = config.get("megahit", {}).get("preset", "meta-large"),
        min_contig = config.get("megahit", {}).get("min_contig_len", 1000),
        outdir = "data/assembly/{sample}"
    conda:
        "envs/megahit.yaml"
    log:
        "logs/assembly_{sample}.log"
    benchmark:
        "benchmarks/assembly_{sample}.tsv"
    threads: config["resources"]["threads_megahit"]
    resources:
        mem_mb   = 32000,
        partition = "medium-small",
        runtime  = "06:00:00"
    shell:
        """
        echo "[MEGAHIT] Assembly untuk {wildcards.sample}..." | tee {log}
        
        mkdir -p data/contigs logs/assembly
        
        # Hapus output dir lama jika ada (mencegah error re-run)
        rm -rf {params.outdir}
        
        megahit \
            -1 {input.r1} -2 {input.r2} \
            -o {params.outdir} \
            --presets {params.preset} \
            --min-contig-len {params.min_contig} \
            --num-cpu-threads {threads} \
            --memory 0.9 \
            2>&1 | tee {output.log_file}
        
        mv {params.outdir}/final.contigs.fa {output.contigs}
        
        echo "[MEGAHIT] Selesai. Total contigs:" >>{log}
        grep -c "^>" {output.contigs} >>{log} || echo "0" >>{log}
        """

# PHASE 3A: Deteksi AMR dengan RGI + CARD Database

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
        input_type = config["rgi"]["input_type"],
        alignment_tool = config["rgi"]["alignment_tool"]
    conda:
        "envs/rgi.yaml"
    log:
        "logs/rgi_{sample}.log"
    benchmark:
        "benchmarks/rgi_{sample}.tsv"
    threads: config["resources"]["threads_rgi"]
    resources:
        mem_mb   = 16000,
        partition = "short",
        runtime  = "03:00:00"
    shell:
        """
        echo "[RGI] Memproses {wildcards.sample}..." | tee {log}
        
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

# PHASE 3B: Deteksi Plasmid dengan MOB-suite

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
    benchmark:
        "benchmarks/mobsuite_{sample}.tsv"
    threads: config["resources"]["threads_mobsuite"]
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "01:00:00"
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

# PHASE 3C: Deteksi Integron dengan IntegronFinder

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
    benchmark:
        "benchmarks/integron_{sample}.tsv"
    threads: config["resources"]["threads_integron"]
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "01:00:00"
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

# PHASE 3D: Deteksi IS / Transposon dengan ISEScan

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
    benchmark:
        "benchmarks/isescan_{sample}.tsv"
    threads: config["resources"]["threads_isescan"]
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "01:00:00"
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

# PHASE 3E: HGT Linkage — Co-localization Integration

rule colocalization:
    """
    Mengintegrasikan hasil keempat tools untuk setiap gen AMR.
    Menentukan posisi: Plasmid / Integron / IS / Kromosom.
    Menambahkan metadata Country dan Region dari sample_map.csv.
    """
    input:
        rgi = expand("results/rgi/{sample}_rgi.txt", sample=SAMPLES),
        mob = expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES),
        integron = expand("results/integron/{sample}_integrons.tsv", sample=SAMPLES),
        isescan = expand("results/isescan/{sample}_is.tsv", sample=SAMPLES),
        map = SAMPLE_MAP_FILE
    output:
        csv = "results/colocalization_summary.csv"
    conda:
        "envs/python.yaml"
    log:
        "logs/colocalization.log"
    benchmark:
        "benchmarks/colocalization.tsv"
    resources:
        mem_mb   = 4000,
        partition = "short",
        runtime  = "00:30:00"
    shell:
        """
        echo "[Co-localization] Mengintegrasikan data MGE dan AMR..." | tee {log}
        python scripts/02_find_colocalization.py \
            --out {output.csv} \
            --sample-map {input.map} \
            2>&1 | tee -a {log}
        echo "[Co-localization] Selesai. Output: {output.csv}" | tee -a {log}
        """

# PHASE 3F: Agregasi per Populasi

rule aggregate_by_population:
    """
    Merangkum colocalization_summary.csv menjadi tiga matrix agregat.
    """
    input:
        csv = "results/colocalization_summary.csv"
    output:
        abund = "results/amr_abundance_matrix.csv",
        mge = "results/mge_distribution_matrix.csv",
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

# PHASE 4: Statistik & Visualisasi

rule run_statistics:
    """
    Analisis statistik komprehensif:
    - 4A: Alpha & beta diversity (Shannon, PERMANOVA)
    - 4B: Kruskal-Wallis per AMR class antar populasi
    - 4C: Chi-square distribusi MGE type
    - 4D: Fisher's exact AMR x MGE (FDR-corrected)
    - 4E: Heatmap + stacked bar + supplementary tables
    """
    input:
        abund = "results/amr_abundance_matrix.csv",
        mge = "results/mge_distribution_matrix.csv",
        assoc = "results/amr_mge_association_matrix.csv",
        coloc = "results/colocalization_summary.csv"
    output:
        fig2 = "results/figures/Fig2_MGE_distribution.pdf"
    conda:
        "envs/r_stats.yaml"
    log:
        "logs/statistics.log"
    benchmark:
        "benchmarks/statistics.tsv"
    threads: config["resources"]["threads_stats"]
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "01:00:00"
    shell:
        """
        echo "[Statistics] Menjalankan 04_run_stats.R..." | tee {log}
        Rscript scripts/04_run_stats.R 2>&1 | tee -a {log}
        echo "[Statistics] Selesai." | tee -a {log}
        """

# PHASE 5: Network Comparison

rule network_analysis:
    """
    Analisis jaringan bipartit AMR-MGE per populasi.
    Membandingkan struktur network antar populasi.
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
    benchmark:
        "benchmarks/network.tsv"
    resources:
        mem_mb   = 8000,
        partition = "short",
        runtime  = "01:00:00"
    shell:
        """
        echo "[Network] Menjalankan 05_network_analysis.R..." | tee {log}
        Rscript scripts/05_network_analysis.R 2>&1 | tee -a {log}
        touch {output.flag}
        echo "[Network] Selesai." | tee -a {log}
        """

# Cleanup temporary files (optional)

rule cleanup_intermediates:
    """
    Hapus file intermediate (raw reads, QC reads, non-host reads)
    untuk menghemat storage setelah assembly sukses.
    """
    input:
        contigs = "data/contigs/{sample}.fasta"
    output:
        touch("logs/cleanup_{sample}.done")
    shell:
        """
        echo "[Cleanup] Removing intermediate files for {wildcards.sample}..."
        
        rm -f data/raw_reads/{wildcards.sample}_*.fastq.gz 2>/dev/null || true
        rm -f data/qc_reads/{wildcards.sample}_*.fastq.gz 2>/dev/null || true
        rm -f data/nonhost_reads/{wildcards.sample}_*.fastq.gz 2>/dev/null || true
        rm -rf data/assembly/{wildcards.sample} 2>/dev/null || true
        
        echo "[Cleanup] Done."
        """
