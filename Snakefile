configfile: "config.yaml"

import os
import pandas as pd

# --- CONFIGURATION & PATHS ---
# Bapak bisa ganti nama RUN_ID ini lewat command line: --config run_id=nama_baru
RUN_ID = config.get("run_id", "pilot_run")
OUT = f"results/{RUN_ID}"

# Otomatis buat folder utama agar SLURM tidak ngambek
# logs/slurm di root = untuk system log dari SLURM (sbatch --output/--error)
# {OUT}/logs/slurm  = untuk per-rule log dari pipeline
os.makedirs("logs/slurm", exist_ok=True)
os.makedirs(f"{OUT}/logs/slurm", exist_ok=True)
os.makedirs(f"{OUT}/benchmarks", exist_ok=True)
os.makedirs(f"{OUT}/tmp", exist_ok=True)

# Rule yang harus jalan di Login Node
localrules: all, download_sra, cleanup_intermediates, download_hg38_index

# --- METADATA LOADING ---
# Support dua cara: --config sample_map=... atau dari config.yaml
_sample_map = config.get("sample_map", config.get("paths", {}).get("sample_map", "data/metadata/pilot_map.csv"))
SAMPLES = pd.read_csv(_sample_map)["sample_id"].unique()

rule all:
    input:
        expand(f"{OUT}/analysis/statistics/{{population}}/summary_stats.txt", population=["East_Asia", "Europe"]),
        expand(f"{OUT}/analysis/networks/{{population}}/network_done.flag", population=["East_Asia", "Europe"])

# PHASE 1: SRA Download & HG38 Index
rule download_sra:
    """Mengunduh data FASTQ dari NCBI SRA"""
    output:
        r1 = f"{OUT}/data/raw_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/raw_reads/{{sample}}_2.fastq.gz"
    params:
        accession = lambda wildcards: pd.read_csv(config["sample_map"]).set_index("sample_id").loc[wildcards.sample, "accession"]
    conda:
        "envs/sra-tools.yaml"
    log:
        f"{OUT}/logs/download_{{sample}}.log"
    threads: 2
    shell:
        """
        echo "[Download] Mengunduh {wildcards.sample} ({params.accession})..." | tee "{log}"
        mkdir -p "{OUT}/data/raw_reads" "{OUT}/tmp/download_{wildcards.sample}"
        
        IFS=';' read -ra ADDR <<< "{params.accession}"
        for acc in "${{ADDR[@]}}"; do
            echo "[Download] Memproses run: $acc" >> "{log}"
            prefetch "$acc" --max-size 50G 2>>"{log}" || echo "[WARN] prefetch $acc gagal..." >>"{log}"
            fasterq-dump --split-files --threads {threads} "$acc" --outdir "{OUT}/tmp/download_{wildcards.sample}" 2>>"{log}"
            cat "{OUT}/tmp/download_{wildcards.sample}/${{acc}}_1.fastq" >> "{OUT}/tmp/combined_{wildcards.sample}_1.fastq"
            cat "{OUT}/tmp/download_{wildcards.sample}/${{acc}}_2.fastq" >> "{OUT}/tmp/combined_{wildcards.sample}_2.fastq"
            rm -rf "$acc" "{OUT}/tmp/download_{wildcards.sample}/${{acc}}*"
        done
        
        pigz -p {threads} -c "{OUT}/tmp/combined_{wildcards.sample}_1.fastq" > "{output.r1}"
        pigz -p {threads} -c "{OUT}/tmp/combined_{wildcards.sample}_2.fastq" > "{output.r2}"
        rm "{OUT}/tmp/combined_{wildcards.sample}_1.fastq" "{OUT}/tmp/combined_{wildcards.sample}_2.fastq"
        rm -rf "{OUT}/tmp/download_{wildcards.sample}"
        """

rule qc_fastp:
    """Quality Control menggunakan fastp"""
    input:
        r1 = f"{OUT}/data/raw_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/raw_reads/{{sample}}_2.fastq.gz"
    output:
        qc1 = f"{OUT}/data/qc_reads/{{sample}}_1.fastq.gz",
        qc2 = f"{OUT}/data/qc_reads/{{sample}}_2.fastq.gz",
        html = f"{OUT}/logs/qc/{{sample}}.html",
        json = f"{OUT}/logs/qc/{{sample}}.json"
    conda:
        "envs/fastp.yaml"
    log:
        f"{OUT}/logs/qc/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        partition = "short",
        runtime = 30
    shell:
        """
        mkdir -p "{OUT}/data/qc_reads" "{OUT}/logs/qc"
        fastp -i "{input.r1}" -I "{input.r2}" -o "{output.qc1}" -O "{output.qc2}" \
              --thread {threads} --html "{output.html}" --json "{output.json}" 2>>"{log}"
        """

rule download_hg38_index:
    """Download hg38 reference untuk host removal"""
    output:
        index = protected("databases/hg38/hg38.1.bt2")
    conda:
        "envs/bowtie2.yaml"
    log:
        f"{OUT}/logs/download_hg38.log"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "medium-small",
        runtime = 240
    shell:
        """
        echo "[HG38] Downloading..." | tee "{log}"
        mkdir -p databases/hg38
        wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > databases/hg38/hg38.fa 2>>"{log}"
        bowtie2-build --threads {threads} databases/hg38/hg38.fa databases/hg38/hg38 2>>"{log}"
        """

rule host_removal:
    """Menghapus kontaminasi DNA manusia"""
    input:
        r1 = f"{OUT}/data/qc_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/qc_reads/{{sample}}_2.fastq.gz",
        index = "databases/hg38/hg38.1.bt2"
    output:
        r1 = f"{OUT}/data/nonhost_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/nonhost_reads/{{sample}}_2.fastq.gz",
        stats = f"{OUT}/logs/host_removal/{{sample}}.stats"
    params:
        index_base = "databases/hg38/hg38"
    conda:
        "envs/bowtie2.yaml"
    log:
        f"{OUT}/logs/host_removal/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "medium-small",
        runtime = 120
    shell:
        """
        bowtie2 -p {threads} -x "{params.index_base}" -1 "{input.r1}" -2 "{input.r2}" \
                --un-conc-gz "{OUT}/data/nonhost_reads/{wildcards.sample}_%.fastq.gz" \
                --very-sensitive > /dev/null 2>"{output.stats}"
        mv "{OUT}/data/nonhost_reads/{wildcards.sample}_1.fastq.gz" "{output.r1}"
        mv "{OUT}/data/nonhost_reads/{wildcards.sample}_2.fastq.gz" "{output.r2}"
        """

rule assembly_megahit:
    """Perakitan metagenome de novo"""
    input:
        r1 = f"{OUT}/data/nonhost_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/nonhost_reads/{{sample}}_2.fastq.gz"
    output:
        contigs = f"{OUT}/data/contigs/{{sample}}.fa",
        log_file = f"{OUT}/logs/assembly/{{sample}}.megahit.log"
    params:
        outdir = f"{OUT}/tmp/assembly/{{sample}}",
        preset = "meta-sensitive",
        min_contig = 500
    conda:
        "envs/megahit.yaml"
    log:
        f"{OUT}/logs/assembly/{{sample}}.log"
    threads: 16
    resources:
        mem_mb = 32000,
        partition = "medium-small",
        runtime = 360
    shell:
        """
        rm -rf "{params.outdir}"
        megahit -1 "{input.r1}" -2 "{input.r2}" -o "{params.outdir}" \
                --out-prefix "{wildcards.sample}" -t {threads} 2>&1 | tee "{output.log_file}"
        mv "{params.outdir}/{wildcards.sample}.contigs.fa" "{output.contigs}"
        """

rule run_rgi:
    """Prediksi gen AMR dengan RGI/CARD"""
    input:
        fasta = f"{OUT}/data/contigs/{{sample}}.fa"
    output:
        rgi = f"{OUT}/analysis/rgi/{{sample}}.tsv"
    conda:
        "envs/rgi.yaml"
    log:
        f"{OUT}/logs/rgi/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 16000,
        partition = "short",
        runtime = 180
    shell:
        """
        rgi main --input_sequence "{input.fasta}" --output_file "{OUT}/analysis/rgi/{wildcards.sample}" \
                 --input_type contig --clean --num_threads {threads} >> "{log}" 2>&1
        mv "{OUT}/analysis/rgi/{wildcards.sample}.txt" "{output.rgi}"
        """

rule run_mobsuite:
    """Deteksi plasmid dengan MOB-suite"""
    input:
        fasta = f"{OUT}/data/contigs/{{sample}}.fa"
    output:
        report = f"{OUT}/analysis/mobsuite/{{sample}}/contig_report.txt"
    conda:
        "envs/mobsuite.yaml"
    log:
        f"{OUT}/logs/mobsuite/{{sample}}.log"
    threads: 8
    resources:
        mem_mb = 8000,
        partition = "short",
        runtime = 60
    shell:
        """
        mob_recon --infile "{input.fasta}" --outdir "{OUT}/analysis/mobsuite/{wildcards.sample}" \
                  --force --num_threads {threads} >> "{log}" 2>&1
        """

rule run_integronfinder:
    """Deteksi integron dengan IntegronFinder 2"""
    input:
        fasta = f"{OUT}/data/contigs/{{sample}}.fa"
    output:
        tsv = f"{OUT}/analysis/integron/{{sample}}.tsv"
    conda:
        "envs/integronfinder.yaml"
    log:
        f"{OUT}/logs/integron/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        partition = "short",
        runtime = 60
    shell:
        """
        integron_finder "{input.fasta}" --outdir "{OUT}/tmp/integron/{wildcards.sample}" --cpu {threads} --local-max >> "{log}" 2>&1
        find "{OUT}/tmp/integron/{wildcards.sample}" -name "*.integrons" | xargs cat > "{output.tsv}" 2>>"{log}" || touch "{output.tsv}"
        """

rule run_isescan:
    """Deteksi Insertion Sequences (IS)"""
    input:
        fasta = f"{OUT}/data/contigs/{{sample}}.fa"
    output:
        tsv = f"{OUT}/analysis/isescan/{{sample}}.tsv"
    conda:
        "envs/isescan.yaml"
    log:
        f"{OUT}/logs/isescan/{{sample}}.log"
    threads: 4
    resources:
        mem_mb = 8000,
        partition = "short",
        runtime = 60
    shell:
        """
        isescan.py --seqfile "{input.fasta}" --output "{OUT}/tmp/isescan/{wildcards.sample}" --nthread {threads} >> "{log}" 2>&1
        find "{OUT}/tmp/isescan/{wildcards.sample}" -name "*.tsv" | head -1 | xargs -I_F_ cp _F_ "{output.tsv}" 2>>"{log}" || touch "{output.tsv}"
        """

rule colocalization:
    """Integrasi AMR dan MGE (Co-localization)"""
    input:
        amr = f"{OUT}/analysis/rgi/{{sample}}.tsv",
        plasmid = f"{OUT}/analysis/mobsuite/{{sample}}/contig_report.txt",
        integron = f"{OUT}/analysis/integron/{{sample}}.tsv",
        is_elements = f"{OUT}/analysis/isescan/{{sample}}.tsv"
    output:
        csv = f"{OUT}/analysis/colocalization/{{sample}}_coloc.csv"
    conda:
        "envs/r_stats.yaml"
    log:
        f"{OUT}/logs/colocalization/{{sample}}.log"
    shell:
        """
        python scripts/02_find_colocalization.py --out "{output.csv}" --sample-map "{config[sample_map]}" 2>&1 | tee -a "{log}"
        """

rule aggregate_by_population:
    """Agregasi data berdasarkan Populasi"""
    input:
        csvs = expand(f"{OUT}/analysis/colocalization/{{sample}}_coloc.csv", sample=SAMPLES)
    output:
        csv = f"{OUT}/analysis/aggregated/{{population}}_combined.csv"
    conda:
        "envs/r_stats.yaml"
    log:
        f"{OUT}/logs/aggregate_{{population}}.log"
    shell:
        """
        python scripts/03_aggregate_by_population.py "{output.csv}" 2>&1 | tee -a "{log}"
        """

rule run_statistics:
    """Analisis Statistik dengan R"""
    input:
        aggregated = f"{OUT}/analysis/aggregated/{{population}}_combined.csv"
    output:
        stats = f"{OUT}/analysis/statistics/{{population}}/summary_stats.txt"
    conda:
        "envs/r_stats.yaml"
    log:
        f"{OUT}/logs/stats_{{population}}.log"
    shell:
        """
        Rscript scripts/04_run_stats.R --input "{input.aggregated}" --output_dir "{OUT}/analysis/statistics/{wildcards.population}" 2>&1 | tee -a "{log}"
        """

rule network_analysis:
    """Analisis Jaringan AMR-MGE"""
    input:
        aggregated = f"{OUT}/analysis/aggregated/{{population}}_combined.csv"
    output:
        flag = f"{OUT}/analysis/networks/{{population}}/network_done.flag"
    conda:
        "envs/r_stats.yaml"
    log:
        f"{OUT}/logs/network_{{population}}.log"
    shell:
        """
        Rscript scripts/05_network_analysis.R --input "{input.aggregated}" --output_dir "{OUT}/analysis/networks/{wildcards.population}" 2>&1 | tee -a "{log}"
        touch "{output.flag}"
        """

rule cleanup_intermediates:
    """Membersihkan file FASTQ setelah assembly selesai (opsional)"""
    input:
        contigs = f"{OUT}/data/contigs/{{sample}}.fa"
    shell:
        """
        rm -f "{OUT}/data/raw_reads/{wildcards.sample}"*.fastq.gz
        rm -f "{OUT}/data/qc_reads/{wildcards.sample}"*.fastq.gz
        rm -f "{OUT}/data/nonhost_reads/{wildcards.sample}"*.fastq.gz
        """
