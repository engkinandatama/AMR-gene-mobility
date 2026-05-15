configfile: "config.yaml"

import os
import pandas as pd

# --- CONFIGURATION & PATHS ---
# Bapak bisa ganti nama RUN_ID ini lewat command line: --config run_id=nama_baru
RUN_ID = config.get("run_id", "pilot_run")
OUT = f"results/{RUN_ID}"

# Otomatis buat folder utama agar SLURM tidak ngambek
# Kita satukan semua log ke dalam folder hasil agar root tetap bersih
os.makedirs(f"{OUT}/logs/slurm", exist_ok=True)
os.makedirs(f"{OUT}/benchmarks", exist_ok=True)
os.makedirs(f"{OUT}/tmp", exist_ok=True)

# Rule yang harus jalan di Login Node
localrules: all, download_sra, cleanup_intermediates, download_hg38_fasta

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
        r1 = temp(f"{OUT}/data/raw_reads/{{sample}}_1.fastq"),
        r2 = temp(f"{OUT}/data/raw_reads/{{sample}}_2.fastq")
    params:
        accession = lambda wildcards: pd.read_csv(config["sample_map"]).set_index("sample_id").loc[wildcards.sample, "accession"]
    conda:
        "envs/sra-tools.yaml"
    log:
        f"{OUT}/logs/download_{{sample}}.log"
    threads: 1
    shell:
        """
        exec > "{log}" 2>&1
        set -x
        
        echo "[Download] Mengunduh {wildcards.sample} ({params.accession})..."
        mkdir -p "{OUT}/data/raw_reads"
        tmp_run="{OUT}/tmp/download_{wildcards.sample}"
        mkdir -p "$tmp_run"
        
        IFS=';' read -ra ADDR <<< "{params.accession}"
        for acc in "${{ADDR[@]}}"; do
            echo "[Download] Memproses run: ${{acc}}"
            # Pakai -O saja yang paling standar untuk semua versi sra-tools
            prefetch "${{acc}}" -O "$tmp_run" --max-size 50G || echo "[WARN] prefetch ${{acc}} gagal, lanjut fasterq-dump..."
            
            # Jika prefetch sukses, file ada di $tmp_run/acc/acc.sra
            # Jika gagal, fasterq-dump akan coba download sendiri ke outdir
            target="$tmp_run/${{acc}}"
            [ ! -d "$target" ] && target="${{acc}}" 
            
            fasterq-dump --split-files --threads {threads} "$target" --outdir "$tmp_run" --temp "$tmp_run"
            
            # Deteksi: Apakah Paired-End (_1 dan _2) atau Single-End (.fastq saja)?
            if [ -f "$tmp_run/${{acc}}_1.fastq" ]; then
                cat "$tmp_run/${{acc}}_1.fastq" >> "$tmp_run/combined_1.fastq"
                cat "$tmp_run/${{acc}}_2.fastq" >> "$tmp_run/combined_2.fastq"
                rm -f "$tmp_run/${{acc}}_1.fastq" "$tmp_run/${{acc}}_2.fastq"
            elif [ -f "$tmp_run/${{acc}}.fastq" ]; then
                echo "[INFO] Sampel ${{acc}} terdeteksi Single-End."
                cat "$tmp_run/${{acc}}.fastq" >> "$tmp_run/combined_1.fastq"
                touch "$tmp_run/combined_2.fastq"
                rm -f "$tmp_run/${{acc}}.fastq"
            fi
            
            rm -rf "${{acc}}"
        done
        
        # Pindahkan file mentah ke folder tujuan
        if [ -f "$tmp_run/combined_1.fastq" ]; then
            mv "$tmp_run/combined_1.fastq" "{output.r1}"
        else
            echo "[ERROR] File combined_1.fastq tidak ditemukan!"
            exit 1
        fi
        
        if [ -f "$tmp_run/combined_2.fastq" ]; then
            mv "$tmp_run/combined_2.fastq" "{output.r2}"
        else
            touch "{output.r2}"
        fi
        
        # Bersihkan folder temp
        rm -rf "$tmp_run"
        """

rule compress_sra:
    """Mengompres file FASTQ di Cluster Node"""
    input:
        r1 = f"{OUT}/data/raw_reads/{{sample}}_1.fastq",
        r2 = f"{OUT}/data/raw_reads/{{sample}}_2.fastq"
    output:
        r1 = f"{OUT}/data/raw_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/raw_reads/{{sample}}_2.fastq.gz"
    conda:
        "envs/sra-tools.yaml"
    threads: 4
    resources:
        mem_mb=16000,
        partition="short",
        runtime=60
    shell:
        """
        pigz -p {threads} -c {input.r1} > {output.r1} && rm -f {input.r1}
        pigz -p {threads} -c {input.r2} > {output.r2} && rm -f {input.r2}
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

# ------------------------------------------------------------
# 1. HOST REMOVAL (hg38 Index)
# ------------------------------------------------------------

rule download_hg38_fasta:
    """Download hg38 reference genome FASTA (Localrule)"""
    output:
        fasta = "databases/hg38/hg38.fa"
    log:
        f"results/{RUN_ID}/logs/download_hg38_fasta.log"
    localrule: True
    shell:
        """
        echo "[HG38] Downloading FASTA..." | tee -a "{log}"
        mkdir -p databases/hg38
        if [ ! -f "{output.fasta}" ]; then
            wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | gunzip -c > "{output.fasta}" 2>> "{log}"
        fi
        """

rule build_hg38_index:
    """Build Bowtie2 index for hg38 (SLURM)"""
    input:
        fasta = "databases/hg38/hg38.fa"
    output:
        index = "databases/hg38/hg38.1.bt2"
    log:
        f"results/{RUN_ID}/logs/build_hg38_index.log"
    threads: 16
    resources:
        mem_mb = 32000,
        runtime = 360
    conda:
        "envs/bowtie2.yaml"
    shell:
        """
        echo "[HG38] Building Bowtie2 Index (this may take 1-2 hours)..." | tee -a "{log}"
        bowtie2-build --threads {threads} "{input.fasta}" databases/hg38/hg38 2>> "{log}"
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
        # Cek ukuran file R2 (jika < 100 byte, berarti dummy/empty)
        r2_size=$(stat -c%s "{input.r2}")
        
        if [ "$r2_size" -lt 100 ]; then
            echo "[Host Removal] Sampel {wildcards.sample} terdeteksi Single-End (SE). Menggunakan mode -U." | tee "{log}"
            bowtie2 -p {threads} -x "databases/hg38/hg38" -U "{input.r1}" \
                --un-gz "{output.r1}" \
                --very-sensitive > /dev/null 2>"{output.stats}"
            # Buat dummy R2 non-host agar rule berikutnya tidak error
            touch "dummy_r2.fastq"
            pigz -c "dummy_r2.fastq" > "{output.r2}"
            rm "dummy_r2.fastq"
        else
            echo "[Host Removal] Sampel {wildcards.sample} terdeteksi Paired-End (PE). Menggunakan mode -1 -2." | tee "{log}"
            bowtie2 -p {threads} -x "databases/hg38/hg38" -1 "{input.r1}" -2 "{input.r2}" \
                --un-conc-gz "{OUT}/data/nonhost_reads/{wildcards.sample}_%.fastq.gz" \
                --very-sensitive > /dev/null 2>"{output.stats}"
        fi
        """

rule assembly_megahit:
    """Metagenome Assembly menggunakan MEGAHIT"""
    input:
        r1 = f"{OUT}/data/nonhost_reads/{{sample}}_1.fastq.gz",
        r2 = f"{OUT}/data/nonhost_reads/{{sample}}_2.fastq.gz"
    output:
        fa = f"{OUT}/data/contigs/{{sample}}.fa",
        m_log = f"{OUT}/logs/assembly/{{sample}}.megahit.log"
    log:
        f"{OUT}/logs/assembly/{{sample}}.log"
    conda:
        "envs/megahit.yaml"
    threads: 16
    resources:
        mem_mb=64000,
        partition="medium-small",
        runtime=360
    shell:
        """
        exec > "{log}" 2>&1
        set -x
        
        mkdir -p "{OUT}/logs/assembly"
        rm -rf "{OUT}/tmp/assembly/{wildcards.sample}"
        
        megahit -1 "{input.r1}" -2 "{input.r2}" \
                -o "{OUT}/tmp/assembly/{wildcards.sample}" \
                --out-prefix "{wildcards.sample}" \
                --mem-flag 2 \
                -t {threads} 2>&1 | tee "{output.m_log}"
        
        if [ -f "{OUT}/tmp/assembly/{wildcards.sample}/{wildcards.sample}.contigs.fa" ]; then
            mv "{OUT}/tmp/assembly/{wildcards.sample}/{wildcards.sample}.contigs.fa" "{output.fa}"
            rm -rf "{OUT}/tmp/assembly/{wildcards.sample}"
        else
            echo "[ERROR] Megahit failed to produce contigs!"
            exit 1
        fi
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
