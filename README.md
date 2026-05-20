# AMR-MGE Co-localization Pipeline

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> **A reproducible end-to-end bioinformatics pipeline for comparative analysis of AMR gene–mobile genetic element associations across metagenomes from different study cohorts or geographic populations.**

---

## Overview

This repository provides a reproducible, end-to-end framework to analyze whether antimicrobial resistance (AMR) gene–mobile genetic element (MGE) associations and co-localization patterns differ across populations, treatment groups, or geographic regions. It resolves both the *what* (AMR gene diversity/abundance) and the *how* (mechanistic linkage via plasmids, integrons, and insertion sequences/transposons).

### Research Questions Addressed
This pipeline is designed to address key ecological and biological questions regarding resistome mobility, such as:
1. Do AMR gene diversity and abundance differ significantly across different study cohorts or populations?
2. Which MGE types (plasmid, integron, insertion sequence/transposon) predominantly carry AMR genes in each cohort?
3. Do AMR class–MGE type co-localization patterns differ across cohorts, suggesting distinct horizontal gene transfer (HGT) routes?

---

## Pipeline Architecture

```
                       INPUT: sample_map.csv
                                 │
 ┌───────────────────────────────┴───────────────────────────────┐
 │  snakemake assembly                                           │
 │  ==================                                           │
 │  1. Download Raw Reads (sra-tools: prefetch + fasterq-dump)   │
 │  2. Quality Control (fastp: adapter/quality filtering)        │
 │  3. Host Depletion (Bowtie2 vs hg38: remove human reads)      │
 │  4. Metagenome Assembly (MEGAHIT: contig generation)          │
 └───────────────────────────────┬───────────────────────────────┘
                                 ▼
                     results/{run_id}/data/contigs/
                                 │
 ┌───────────────────────────────┴───────────────────────────────┐
 │  snakemake annotation_aggregation                             │
 │  ================================                             │
 │  1. AMR Detection (RGI + CARD database)                       │
 │  2. MGE Profiling:                                            │
 │     - Plasmids (MOB-suite)                                    │
 │     - Integrons (IntegronFinder)                              │
 │     - Insertion Sequences & Transposons (ISEScan)             │
 │  3. Co-localization Integration (02_find_colocalization.py)   │
 │  4. Regional Matrix Aggregation (03_aggregate_by_population)  │
 └───────────────────────────────┬───────────────────────────────┘
                                 ▼
                     results/{run_id}/analysis/aggregated/
                                 │
 ┌───────────────────────────────┴───────────────────────────────┐
 │  snakemake statistics                                         │
 │  ====================                                         │
 │  Resistome comparisons, alpha/beta diversity,                 │
 │  association significance (Fisher's exact, PERMANOVA, etc.)   │
 └───────────────────────────────┬───────────────────────────────┘
                                 ▼
                     results/{run_id}/analysis/statistics/
                                 │
 ┌───────────────────────────────┴───────────────────────────────┐
 │  snakemake networks                                           │
 │  ==================                                           │
 │  Bipartite AMR class – MGE type co-occurrence networks        │
 │  and topology centrality analysis                             │
 └───────────────────────────────┬───────────────────────────────┘
                                 ▼
                     results/{run_id}/analysis/networks/
```

### Hybrid Execution Mode

The pipeline supports multiple entry points:

1. **End-to-End (From Raw Reads)**: Run the full workflow starting from NCBI SRA accession IDs to download, clean, assemble, annotate, and analyze the samples.
2. **From Pre-assembled Contigs (Annotation & Analysis)**: If you already have assembled contigs, bypass read downloading and assembly, and start directly with AMR/MGE detection and aggregation.
3. **Analysis-only**: Run statistics and network comparisons directly on pre-computed co-localization and abundance matrices.

---

## Repository Structure

```
AMR-gene-mobility/
├── config.yaml                # ← Central config: threads, tool parameters
├── Snakefile                  # ← Pipeline orchestrator (v4.1 - End-to-End & Modular)
├── profiles/
│   └── slurm/
│       └── config.yaml        # SLURM profile for HPC BRIN
├── data/
│   ├── raw_reads/             # Downloaded FASTQ files (temporary)
│   ├── qc_reads/              # QC'd reads (temporary)
│   ├── nonhost_reads/         # Non-human reads (temporary)
│   ├── contigs/               # Assembled contigs (permanent)
│   └── metadata/
│       └── sample_map.csv     # Sample ID → country/region/accession mapping
├── databases/
│   └── hg38/                  # Human reference genome index (downloaded once)
├── envs/                      # Conda environment definitions
│   ├── sra-tools.yaml         # SRA download
│   ├── fastp.yaml             # Quality control
│   ├── bowtie2.yaml           # Host removal
│   ├── megahit.yaml           # Metagenome assembly
│   ├── rgi.yaml               # AMR detection
│   ├── mobsuite.yaml          # Plasmid detection
│   ├── integronfinder.yaml    # Integron detection
│   ├── isescan.yaml           # IS/Transposon detection
│   ├── python.yaml            # Co-localization & aggregation
│   └── r_stats.yaml           # Statistics & network analysis
├── scripts/
│   ├── 01_fetch_metadata.R          # Metadata download & sample cohort filtering
│   ├── 01b_fetch_taxonomy.R         # Taxonomy profiling of cohort samples
│   ├── filter_contigs.py            # Contig filtering by minimum length
│   ├── 02_find_colocalization.py    # AMR-MGE co-localization mapping
│   ├── 03_aggregate_by_population.py # Abundance & co-localization matrix aggregation
│   ├── 04_run_stats.R               # Statistical testing and resistome profiling
│   ├── 05_network_analysis.R        # Bipartite network construction & topology metrics
│   └── analyze_benchmarks.py        # Utility: Compile run duration & memory benchmarks
└── results/                         # Generated outputs (grouped by run_id, e.g. pilot_run)
    └── pilot_run/
        ├── data/                    # Contigs
        ├── analysis/                # Subdivided into rgi, mobsuite, integron, isescan, aggregated, statistics, networks
        ├── logs/                    # Rule log files and SLURM logs
        └── benchmarks/              # Job execution time and memory profiles
```

---

## Prerequisites

### 1. Package Manager (Conda / Miniconda)

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Miniforge](https://github.com/conda-forge/miniforge) if not already available.

### 2. Snakemake Installation

To avoid package conflicts in your base environment, it is highly recommended to install Snakemake in a dedicated environment. You can use standard **Conda** or **Mamba** (if installed):

#### Option A: Using Conda (Standard)
```bash
# Create a dedicated environment and install Snakemake
conda create -c conda-forge -c bioconda -n snakemake snakemake -y

# Activate the environment
conda activate snakemake
```

#### Option B: Using Mamba (Alternative)
```bash
# Create environment using mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake -y

# Activate the environment
conda activate snakemake
```

> **Note**: All other software tools (R, Python, MEGAHIT, Bowtie2, RGI, MOB-suite, etc.) are **automatically downloaded and installed** inside isolated environments by Snakemake during the first run using the `--use-conda` flag. No manual setup is needed!

### 3. Hardware Requirements

| Component | Minimum | Recommended (HPC) |
|-----------|---------|-------------------|
| RAM | 16 GB | 32-64 GB per parallel job |
| CPU | 4 cores | 8+ cores per parallel job |
| Storage | 100 GB | ~500 GB (for 250 samples, thanks to automatic intermediate file cleanup) |
| Internet | Required | Required (for reference database & SRA download) |

---

## Usage

### Step 1 — Clone the repository

```bash
git clone https://github.com/engkinandatama/AMR-gene-mobility.git
cd AMR-gene-mobility
```

### Step 2 — Prepare sample metadata

Generate sample list using the R script:

```bash
Rscript scripts/01_fetch_metadata.R
```

This will create `data/metadata/sample_map.csv` with curated samples from `curatedMetagenomicData`.

Or manually create `data/metadata/sample_map.csv`:

```csv
sample_id,country,country_name,region,accession
sample_01,DNK,Denmark,Europe,ERR321618
sample_02,CHN,China,East_Asia,SRR9108951
sample_03,IND,India,South_Asia,ERR2017479
```

### Step 3 — Run the pipeline

#### Option A: End-to-End (From Accession IDs)

```bash
# Dry-run to validate
snakemake all --dry-run --cores 4

# Run on HPC with SLURM (uses profiles/slurm/config.yaml configuration)
snakemake all --profile profiles/slurm --use-conda

# Or run locally (not recommended for large batches due to assembly load)
snakemake all --cores all --use-conda
```

#### Option B: From Assembled Contigs (Backward Compatible)

Place your assembled `.fa` or `.fasta` files in `results/{run_id}/data/contigs/` (or specify via `paths.contigs_dir` in `config.yaml`):

```
results/pilot_run/data/contigs/
├── DNK_ERR321618.fa        # Denmark (Europe)
├── CHN_SRR9108951.fa       # China (East Asia)
└── IND_ERR2017479.fa       # India (South Asia)
```

Then run the remaining annotation and analysis phases:

```bash
# Run annotation and aggregation only
snakemake annotation_aggregation --cores all --use-conda

# Or run everything from annotation to final plots
snakemake all --cores all --use-conda
```

#### Option C: Functional Component Execution

```bash
# Run Metagenome Assembly (Download, QC, and MEGAHIT)
snakemake assembly --cores all --use-conda

# Run AMR & MGE Detection, Co-localization, and Aggregation
snakemake annotation_aggregation --cores all --use-conda

# Run Statistical Analyses & Figures
snakemake statistics --cores all --use-conda

# Run Bipartite Network Construction & Centrality Comparison
snakemake networks --cores all --use-conda
```

### Step 4 — Monitor progress

```bash
# Check SLURM queue
squeue -u $USER

# View logs
tail -f logs/download_*.log
tail -f logs/assembly_*.log
```

### Step 5 — Clean up intermediate files (optional)

```bash
# Remove raw reads and intermediate files to save storage
snakemake cleanup_intermediates --cores 1
```

---

## Configuration

All tunable parameters are centralized in `config.yaml`:

### Resource Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `resources.threads_rgi` | 8 | CPU threads for RGI |
| `resources.threads_mobsuite` | 4 | CPU threads for MOB-suite |
| `resources.threads_integron` | 4 | CPU threads for IntegronFinder |
| `resources.threads_isescan` | 4 | CPU threads for ISEScan |
| `resources.threads_stats` | 4 | CPU threads for R statistics |
| `resources.max_memory_mb` | 32000 | Maximum memory per job (MB) |

### Tool Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `fastp.qualified_quality` | 20 | Minimum Phred quality score |
| `fastp.length_required` | 50 | Minimum read length after trimming |
| `megahit.preset` | `meta-large` | Assembly preset for gut metagenome |
| `megahit.min_contig_len` | 1000 | Minimum contig length |
| `rgi.alignment_tool` | `BLAST` | `BLAST` (accurate) or `DIAMOND` (fast) |
| `integronfinder.mode` | `local-max` | IntegronFinder search mode |
| `stats.p_adjust_method` | `BH` | Multiple testing correction method |
| `stats.alpha` | `0.05` | Significance threshold |

### Sample Selection Criteria

Edit `config.yaml` to customize sample selection:

```yaml
sample_selection:
  n_per_country: 50          # Target samples per country
  min_reads: 10000000        # Minimum 10 million reads
  disease_status: "healthy"  # Only healthy subjects
  countries:                 # Target countries
    DNK: "Europe"
    CHN: "East_Asia"
    IND: "South_Asia"
    CMR: "Sub-Saharan_Africa"
    MDG: "Sub-Saharan_Africa"
```

For HPC with many cores, update thread values accordingly:
```yaml
resources:
  threads_rgi:      16
  threads_mobsuite: 8
  threads_integron: 8
  threads_isescan:  8
```

---

## Outputs

All output files are organized per run in the `results/{run_id}/` directory (default: `results/pilot_run/`).

### Primary Output Files

| File | Description |
|------|-------------|
| `results/{run_id}/data/contigs/{sample}.fa` | Assembled contigs filtered by length |
| `results/{run_id}/analysis/aggregated/{population}_all_coloc.csv` | Full co-localization table for the population (contains AMR-MGE association positions) |
| `results/{run_id}/analysis/aggregated/{population}_amr_abundance.csv` | AMR drug class abundance matrix per sample |
| `results/{run_id}/analysis/aggregated/{population}_mge_distribution.csv` | Proportional distribution of MGE types (Plasmid, Integron, IS) per country |
| `results/{run_id}/analysis/aggregated/{population}_combined.csv` | Cross-tabulated AMR class × MGE type × Country count matrix |

### Figures

Figures are generated under the respective population statistics and networks directories:

| File | Description |
|------|-------------|
| `results/{run_id}/analysis/statistics/{population}/figures/Fig1_Heatmap_AMR.pdf` | Heatmap of AMR drug class abundance across samples |
| `results/{run_id}/analysis/statistics/{population}/figures/Fig2_MGE_distribution.pdf` | Stacked bar plot of MGE type distribution per country |
| `results/{run_id}/analysis/statistics/{population}/figures/Fig_S1_alpha_diversity.pdf` | Boxplot comparing AMR gene richness across countries |
| `results/{run_id}/analysis/statistics/{population}/figures/Fig_S2_permdisp_plot.pdf` | PERMDISP group dispersion visualization |
| `results/{run_id}/analysis/statistics/{population}/figures/Fig_S3_procrustes_plot.pdf` | Procrustes overlay plot (Taxonomy vs. Resistome) |
| `results/{run_id}/analysis/statistics/{population}/figures/Fig_S4_mobility_index.pdf` | Distribution of AMR Mobility Index per country |
| `results/{run_id}/analysis/networks/{population}/figures/Fig3_Network_{country}.pdf` | Bipartite AMR class – MGE type co-occurrence network per country |

### Statistical & Network Tables

| File | Description |
|------|-------------|
| `results/{run_id}/analysis/statistics/{population}/tables/Table1_KruskalWallis_AMR_Class.csv` | Kruskal-Wallis and Dunn post-hoc test results comparing AMR class abundances (Sub-Q1) |
| `results/{run_id}/analysis/statistics/{population}/tables/Table2_Fisher_AMR_MGE.csv` | Fisher's Exact test results for association significance (Sub-Q3) |
| `results/{run_id}/analysis/statistics/{population}/tables/Table_PERMANOVA.csv` | PERMANOVA test comparing resistome composition across countries |
| `results/{run_id}/analysis/statistics/{population}/tables/Table_Network_Metrics.csv` | Topology comparison metrics across countries (Sub-Q3) |
| `results/{run_id}/analysis/statistics/{population}/tables/Supplementary_Statistics.xlsx` | Comprehensive multi-sheet workbook containing all tables and values |

### QC Reports

| File | Description |
|------|-------------|
| `logs/qc/fastp_{sample}.html` | FastP QC report per sample |
| `logs/host_removal/{sample}_bowtie2.stats` | Host removal statistics |
| `logs/assembly/{sample}_megahit.log` | Assembly log and statistics |

---

## HPC Deployment (SLURM)

### Resource Estimates (250 Samples)

| Resource | Total Required |
|----------|----------------|
| CPU-hours | ~5,000 |
| Storage | 3 TB |
| Wall time | 3-5 days (parallel) |
| RAM per job | 16-32 GB |

### SLURM Configuration

The pipeline includes a pre-configured SLURM profile template at `profiles/slurm/config.yaml`. 

Before running, adjust the settings to match your cluster's scheduler and partitions, and update the default resource configuration (e.g., default partition, memory limits, and billing account code if required by your HPC).

### Submit Jobs

```bash
# Load modules (if required by HPC)
module load miniconda3

# Activate conda
conda activate snakemake

# Run pipeline with SLURM scheduler integration
snakemake all --profile profiles/slurm --use-conda --jobs 20
```

### Storage Management

The pipeline automatically manages storage:

- **Temporary files**: Raw FASTQ, QC reads, non-host reads are deleted after assembly
- **Permanent files**: Contigs, results, logs are preserved
- **Peak storage**: ~500 GB active (not 1.8 TB)

To manually clean up:
```bash
snakemake cleanup_intermediates --cores 1
```

---

## Data Availability

- Metagenomic metadata sourced from public databases (e.g., [`curatedMetagenomicData`](https://bioconductor.org/packages/curatedMetagenomicData/)) or user-provided datasets.
- Raw sequencing data available from NCBI SRA under accession numbers defined in your metadata map.
- Target cohorts: Supports arbitrary geographic populations, treatment groups, or study cohorts defined in your metadata map.

---

## Tools & Databases Used

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| [sra-tools](https://github.com/ncbi/sra-tools) | ≥3.0 | SRA download | NCBI |
| [fastp](https://github.com/OpenGene/fastp) | ≥0.23 | Quality control | Chen et al., 2018 |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) | ≥2.5 | Host removal | Langmead & Salzberg, 2012 |
| [MEGAHIT](https://github.com/voutcn/megahit) | ≥1.2 | Metagenome assembly | Li et al., 2015 |
| [RGI](https://github.com/arpcard/rgi) | ≥6.0 | AMR gene detection | Alcock et al., 2023 |
| [CARD](https://card.mcmaster.ca/) | ≥3.2 | AMR database | Alcock et al., 2023 |
| [MOB-suite](https://github.com/phac-nml/mob-suite) | ≥3.1 | Plasmid detection | Robertson & Nash, 2018 |
| [IntegronFinder](https://github.com/gem-pasteur/Integron_Finder) | ≥2.0 | Integron detection | Cury et al., 2016 |
| [ISEScan](https://github.com/xiezhq/ISEScan) | ≥1.7 | IS/transposon detection | Xie & Tang, 2017 |
| [Snakemake](https://snakemake.readthedocs.io) | ≥7.0 | Workflow management | Mölder et al., 2021 |

---

## Troubleshooting

### Common Issues

**1. SRA download fails**
```bash
# Try alternative download method
prefetch --max-size 50G SRR123456
fasterq-dump --split-files SRR123456
```

**2. Assembly runs out of memory**
```bash
# Reduce memory fraction in config.yaml
megahit:
  preset: "meta-sensitive"  # Less memory-intensive preset
```

**3. Conda environment conflicts**
```bash
# Clean conda cache
conda clean --all

# Recreate environments
snakemake --conda-cleanup-pkgs cache --use-conda
```

**4. SLURM job timeout**
```yaml
# Increase walltime in profiles/slurm/config.yaml
assembly_megahit:
  time: "08:00:00"  # Increase from 4h to 8h
```

---

## Citation

If you use this pipeline, please cite:

> Nandatama, E. (2026). *AMR-MGE Co-localization Pipeline: Comparative analysis of resistome mobility in metagenomes*. GitHub. https://github.com/engkinandatama/AMR-gene-mobility

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Acknowledgments

We thank the **Badan Riset dan Inovasi Nasional (BRIN) HPC Team** for providing the computational services and resources (Mahameru HPC System) used to perform the metagenomic assembly and analysis in our study.

We also thank the developers and maintainers of the **curatedMetagenomicData** Bioconductor package for providing the curated human gut metagenomic metadata and taxonomy profiles that served as the sample selection basis for this work.

---

## Contact

**Engki Nandatama**  
GitHub: [@engkinandatama](https://github.com/engkinandatama)
