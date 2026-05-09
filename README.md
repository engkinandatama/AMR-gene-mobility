# AMR Gene Mobility in Human Gut Metagenomes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> **A reproducible end-to-end bioinformatics pipeline for comparative analysis of AMR gene–mobile genetic element associations across human gut metagenomes from Europe, Asia, and Africa.**

---

## Overview

This repository contains the full computational pipeline for the study:

**"Cross-Continental Patterns of AMR Gene-MGE Co-localization in Human Gut Metagenomes: A One Health Perspective"**

This pipeline investigates whether AMR gene–MGE associations differ across populations with different antibiotic usage profiles (Denmark, China, India, Cameroon, Madagascar), addressing both the *what* (AMR gene diversity) and the *how* (mechanistic linkage via plasmids, integrons, and transposons).

### Research Questions
1. Do AMR gene diversity and abundance differ significantly across populations from Europe, Asia, and Africa?
2. Which MGE types (plasmid, integron, insertion sequence/transposon) predominantly carry AMR genes in each population?
3. Do AMR class–MGE type associations differ across populations, indicating distinct horizontal gene transfer (HGT) routes?

---

## Pipeline Architecture

```
END-TO-END PIPELINE (HPC BRIN)
================================================================================

INPUT: Accession IDs (SRR/ERR) from sample_map.csv
         │
         ▼ [Phase 1] Download Raw Reads
    sra-tools (prefetch + fasterq-dump)
         │
         ▼ [Phase 2A] Quality Control
    fastp (adapter trimming, quality filtering)
         │
         ▼ [Phase 2B] Host Depletion
    Bowtie2 vs hg38 (remove human reads)
         │
         ▼ [Phase 2C] De Novo Assembly
    MEGAHIT (metagenome assembler)
         │
         ▼ data/contigs/{sample}.fasta
         │
         ├──────────────────────────────────────┐
         │                                      │
         ▼ [Phase 3A] AMR Detection             ▼ [Phase 3B-D] MGE Profiling
    RGI + CARD Database                  MOB-suite (Plasmid)
                                         IntegronFinder (Integron)
                                         ISEScan (IS elements)
         │                                      │
         └──────────────┬───────────────────────┘
                        │
                        ▼ [Phase 3E] Co-localization
                  HGT Linkage Analysis
                  (AMR gene → MGE type per contig)
                        │
                        ▼ [Phase 3F] Aggregation
                  Matrices: AMR abundance, MGE distribution
                        │
                        ├──────────────────────┐
                        │                      │
                        ▼ [Phase 4]            ▼ [Phase 5]
                  Statistics              Network Analysis
                  Kruskal-Wallis          Bipartite Networks
                  Fisher's Exact         Jaccard Similarity
                  PERMANOVA              Hub Identification
                        │                      │
                        └──────────┬───────────┘
                                   ▼
                           OUTPUT: Publication-ready
                           Tables, Figures, Networks
```

### Hybrid Execution Mode

The pipeline supports multiple input types:

1. **End-to-end (Full Pipeline)**: Start from accession IDs → download → QC → assembly → analysis
2. **Assembly-only (Phase 3+)**: Start from assembled contigs → AMR/MGE detection → analysis
3. **Analysis-only (Phase 4-5)**: Run statistics and network analysis on existing results

---

## Repository Structure

```
AMR-gene-mobility/
├── config.yaml                # ← Central config: threads, tool parameters
├── Snakefile                  # ← Pipeline orchestrator (v4.0 - End-to-End)
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
│   ├── sra-tools.yaml         # NEW: SRA download
│   ├── fastp.yaml             # NEW: Quality control
│   ├── bowtie2.yaml           # NEW: Host removal
│   ├── megahit.yaml           # NEW: Assembly
│   ├── rgi.yaml
│   ├── mobsuite.yaml
│   ├── integronfinder.yaml
│   ├── isescan.yaml
│   ├── python.yaml
│   └── r_stats.yaml
├── scripts/
│   ├── 01_fetch_metadata.R          # Phase 0: Sample selection from curatedMetagenomicData
│   ├── 02_find_colocalization.py    # Phase 3: AMR-MGE integration
│   ├── 03_aggregate_by_population.py # Phase 3: Summary matrices
│   ├── 04_run_stats.R               # Phase 4: Statistics & visualization
│   └── 05_network_analysis.R        # Phase 5: Network comparison
├── results/                   # Generated outputs
│   ├── rgi/
│   ├── mobsuite/
│   ├── integron/
│   ├── isescan/
│   ├── figures/
│   └── tables/
└── logs/                      # Per-rule log files
```

---

## Prerequisites

### 1. Conda / Miniconda

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if not already available.

### 2. Snakemake

```bash
conda install -n base -c conda-forge -c bioconda snakemake -y
```

> All tools are **automatically installed** by Snakemake via `--use-conda`. No manual installation required.

### 3. Hardware Requirements

| Component | Minimum | Recommended (HPC) |
|-----------|---------|-------------------|
| RAM | 16 GB | 32-64 GB per job |
| CPU | 4 cores | 8+ cores per job |
| Storage | 100 GB | 3 TB (for 250 samples) |
| Internet | Required | Required (for SRA download) |

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
DNK_ERR321618,DNK,Denmark,Europe,ERR321618
CHN_SRR9108951,CHN,China,East_Asia,SRR9108951
IND_ERR2017479,IND,India,South_Asia,ERR2017479
```

### Step 3 — Run the pipeline

#### Option A: End-to-End (From Accession IDs)

```bash
# Dry-run to validate
snakemake all_full --dry-run --cores 4

# Run on HPC with SLURM
snakemake all_full --profile profiles/slurm --use-conda

# Or run locally (not recommended for large batches)
snakemake all_full --cores all --use-conda
```

#### Option B: From Assembled Contigs (Backward Compatible)

Place your assembled `.fasta` files in `data/contigs/`:

```
data/contigs/
├── DNK_ERR321618.fasta     # Denmark (Europe)
├── CHN_SRR9108951.fasta    # China (East Asia)
└── IND_ERR2017479.fasta    # India (South Asia)
```

Then run:

```bash
# Phase 3 only (AMR + MGE detection)
snakemake phase3 --cores all --use-conda

# Or run everything
snakemake all --cores all --use-conda
```

#### Option C: Phase-by-Phase

```bash
# Phase 1-2: Download, QC, Assembly
snakemake all_full --cores all --use-conda

# Phase 3: AMR + MGE Detection
snakemake phase3 --cores all --use-conda

# Phase 4: Statistics
snakemake phase4 --cores 4 --use-conda

# Phase 5: Network Analysis
snakemake phase5 --cores 4 --use-conda
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

### Primary Output Files

| File | Description |
|------|-------------|
| `data/contigs/{sample}.fasta` | Assembled contigs (permanent) |
| `results/colocalization_summary.csv` | Per-gene AMR–MGE co-localization table |
| `results/amr_abundance_matrix.csv` | AMR class abundance matrix per sample |
| `results/mge_distribution_matrix.csv` | MGE type proportions per country |
| `results/amr_mge_association_matrix.csv` | AMR class × MGE type cross-tabulation |

### Figures

| File | Description |
|------|-------------|
| `results/figures/Fig1_Heatmap_AMR.pdf` | Heatmap: AMR abundance across samples |
| `results/figures/Fig2_MGE_distribution.pdf` | Stacked bar: MGE distribution per population |
| `results/figures/Fig3_Network_*.pdf` | Bipartite AMR–MGE network per population |

### Statistical Tables

| File | Description |
|------|-------------|
| `results/tables/Table1_KruskalWallis_AMR_Class.csv` | Kruskal-Wallis results (Sub-Q1) |
| `results/tables/Table2_Fisher_AMR_MGE.csv` | Fisher's exact results (Sub-Q3) |
| `results/tables/Table_PERMANOVA.csv` | PERMANOVA results |
| `results/tables/Table_Network_Metrics.csv` | Network comparison metrics |
| `results/tables/Supplementary_Statistics.xlsx` | All statistical outputs |

### QC Reports

| File | Description |
|------|-------------|
| `logs/qc/fastp_{sample}.html` | FastP QC report per sample |
| `logs/host_removal/{sample}_bowtie2.stats` | Host removal statistics |
| `logs/assembly/{sample}_megahit.log` | Assembly log and statistics |

---

## HPC BRIN Deployment

### Resource Estimates (250 Samples)

| Resource | Total Required |
|----------|----------------|
| CPU-hours | ~5,000 |
| Storage | 3 TB |
| Wall time | 3-5 days (parallel) |
| RAM per job | 16-32 GB |

### SLURM Configuration

The pipeline includes a pre-configured SLURM profile at `profiles/slurm/config.yaml`. 

Before running, update the account name:

```yaml
__default__:
  account: "your_project_code"  # Replace with your BRIN project code
```

### Submit Jobs

```bash
# Load modules (if required by HPC)
module load miniconda3

# Activate conda
conda activate snakemake

# Run pipeline
snakemake all_full --profile profiles/slurm --use-conda --jobs 100
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

- Metagenomic metadata sourced from [`curatedMetagenomicData`](https://bioconductor.org/packages/curatedMetagenomicData/)
- Raw sequencing data available from NCBI SRA under accession numbers listed in `data/metadata/sample_map.csv`
- Target cohorts: Denmark (Europe), China (East Asia), India (South Asia), Cameroon (Sub-Saharan Africa), Madagascar (Sub-Saharan Africa)

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

> Nandatama, E. (2026). *AMR Gene Mobility in Human Gut Metagenomes: Cross-Population Comparative Pipeline*. GitHub. https://github.com/engkinandatama/AMR-gene-mobility

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Engki Nandatama**  
GitHub: [@engkinandatama](https://github.com/engkinandatama)
