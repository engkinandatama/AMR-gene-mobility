# Snakefile: Menganalisis hasil contigs dari Galaxy
# Jalankan ini nanti dengan perintah: snakemake --cores all --use-conda

import glob
import os

# Otomatis mendeteksi file .fasta apa pun yang nanti Anda taruh di dalam folder data/contigs/
SAMPLES = [os.path.basename(f).replace('.fasta', '') for f in glob.glob("data/contigs/*.fasta")]

rule all:
    input:
        "results/colocalization_summary.csv"

rule run_rgi:
    """Deteksi gen resistensi dari bakteri (membutuhkan memori ringan, dijalankan per sampel)"""
    input:
        "data/contigs/{sample}.fasta"
    output:
        "results/rgi/{sample}_rgi.txt"
    conda:
        "envs/rgi.yaml"
    log:
        "logs/rgi_{sample}.log"
    shell:
        "rgi main -i {input} -o results/rgi/{wildcards.sample}_rgi -t contig -a BLAST --clean > {log} 2>&1"

rule run_mobsuite:
    """Deteksi MGE/Plasmid (membutuhkan memori ringan, dijalankan per sampel)"""
    input:
        "data/contigs/{sample}.fasta"
    output:
        "results/mobsuite/{sample}_plasmid.txt"
    conda:
        "envs/mobsuite.yaml"
    log:
        "logs/mobsuite_{sample}.log"
    shell:
        "mob_recon --infile {input} --outdir results/mobsuite/{wildcards.sample} > {log} 2>&1 && "
        "mv results/mobsuite/{wildcards.sample}/mob_recon_report.txt {output}"

rule data_integration:
    """Menyatukan data RGI dan MOBsuite memakan Skrip Python yang akan kita buat nanti"""
    input:
        rgi=expand("results/rgi/{sample}_rgi.txt", sample=SAMPLES),
        mob=expand("results/mobsuite/{sample}_plasmid.txt", sample=SAMPLES)
    output:
        "results/colocalization_summary.csv"
    conda:
        "envs/python.yaml"
    shell:
        "python scripts/02_find_colocalization.py --out {output}"