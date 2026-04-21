# AMR Gene Mobility in Human Gut Metagenomes

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

> **A reproducible bioinformatics pipeline for comparative analysis of AMR gene–mobile genetic element associations across human gut metagenomes from Europe and Asia.**

---

## Overview

This repository contains the full computational pipeline for the study:

**"Cross-Population Comparative Analysis of AMR Gene Mobility in Human Gut Metagenomes: Association with Mobile Genetic Elements across European and Asian Cohorts"**

This pipeline investigates whether AMR gene–MGE associations differ across populations with different antibiotic usage profiles (Denmark, China, India), addressing both the *what* (AMR gene diversity) and the *how* (mechanistic linkage via plasmids, integrons, and transposons).

### Research Questions
1. Do AMR gene diversity and abundance differ significantly between European and Asian gut metagenomes?
2. Which MGE types (plasmid, integron, insertion sequence/transposon) predominantly carry AMR genes in each population?
3. Do AMR class–MGE type associations differ across populations, indicating distinct horizontal gene transfer (HGT) routes?

---

## Pipeline Architecture

```
Galaxy (cloud)                    HPC / Local
──────────────────────            ─────────────────────────────────────────────
Raw reads (FASTQ)          →      Assembled contigs (.fasta)
  └─ QC: fastp                       │
  └─ Assembly: MEGAHIT               ▼
                                   Snakemake (--use-conda)
                                   ├── RGI          → AMR gene detection (CARD)
                                   ├── MOB-suite    → Plasmid identification
                                   ├── IntegronFinder → Integron detection
                                   └── ISEScan      → IS / Transposon detection
                                         │
                                         ▼
                                   co-localization_summary.csv
                                         │
                                   ├── Aggregation by population
                                   ├── Statistical analysis (R)
                                   │    ├── Shannon / PERMANOVA (Sub-Q1)
                                   │    ├── Kruskal-Wallis (Sub-Q1)
                                   │    ├── Chi-square (Sub-Q2)
                                   │    └── Fisher's exact + FDR (Sub-Q3)
                                   └── Network comparison (igraph/ggraph)
```

---

## Repository Structure

```
AMR-gene-mobility/
├── config.yaml                # ← Central config: threads, tool parameters
├── Snakefile                  # ← Pipeline orchestrator
├── data/
│   ├── contigs/               # Place downloaded .fasta files here (gitignored)
│   └── metadata/
│       └── sample_map.csv     # Sample ID → country/region mapping
├── envs/                      # Conda environment definitions
│   ├── rgi.yaml
│   ├── mobsuite.yaml
│   ├── integronfinder.yaml
│   ├── isescan.yaml
│   ├── python.yaml
│   └── r_stats.yaml
├── scripts/
│   ├── 01_fetch_metadata.R          # Phase 1: Sample selection
│   ├── 02_find_colocalization.py    # Phase 3: AMR-MGE integration
│   ├── 03_aggregate_by_population.py # Phase 3: Summary matrices
│   ├── 04_run_stats.R               # Phase 4: Statistics & visualization
│   └── 05_network_analysis.R        # Phase 5: Network comparison
├── results/                   # Generated outputs (gitignored)
│   ├── rgi/
│   ├── mobsuite/
│   ├── integron/
│   ├── isescan/
│   ├── figures/
│   └── tables/
└── logs/                      # Per-rule log files (gitignored)
```

---

## Prerequisites

### 1. Conda / Miniconda

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) if not already available.

### 2. Snakemake

```bash
conda install -n base -c conda-forge -c bioconda snakemake -y
```

> All other tools (RGI, MOB-suite, IntegronFinder, ISEScan, R packages) are **automatically installed** by Snakemake via `--use-conda`. No manual installation required.

---

## Usage

### Step 1 — Clone the repository

```bash
git clone https://github.com/engkinandatama/AMR-gene-mobility.git
cd AMR-gene-mobility
```

### Step 2 — Prepare your assembled contigs

Place your assembled `.fasta` files in `data/contigs/` following this naming convention:

```
data/contigs/
├── DNK_ERR321618.fasta     # Denmark (Europe)
├── CHN_SRR9108951.fasta    # China (East Asia)
└── IND_ERR2017479.fasta    # India (South Asia)
```

> **Naming convention:** `{COUNTRY_CODE}_{SRA_ACCESSION}.fasta`  
> Country codes follow ISO 3166-1 alpha-3 (DNK, CHN, IND, etc.)

For large batches, download contigs directly from UseGalaxy to HPC:
```bash
wget -O data/contigs/DNK_ERR321618.fasta "YOUR_GALAXY_DOWNLOAD_LINK"
```

### Step 3 — Validate the pipeline (dry-run)

```bash
snakemake phase3 --dry-run --cores 4
```

### Step 4 — Run Phase 3 (AMR + MGE Detection, HPC recommended)

```bash
# In tmux to prevent session loss
tmux new -s amr
snakemake phase3 --cores all --use-conda
```

### Step 5 — Run Phase 4 & 5 (Statistics & Network, can run locally)

```bash
snakemake phase4 --cores 4 --use-conda
snakemake phase5 --cores 4 --use-conda  # Optional but recommended
```

### Run everything at once

```bash
snakemake all --cores all --use-conda
```

---

## Configuration

All tunable parameters are centralized in `config.yaml`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `resources.threads_rgi` | 4 | CPU threads for RGI |
| `resources.threads_mobsuite` | 4 | CPU threads for MOB-suite |
| `resources.threads_integron` | 2 | CPU threads for IntegronFinder |
| `resources.threads_isescan` | 2 | CPU threads for ISEScan |
| `rgi.alignment_tool` | `BLAST` | `BLAST` (accurate) or `DIAMOND` (fast) |
| `rgi.input_type` | `contig` | Input type for RGI |
| `integronfinder.mode` | `local-max` | IntegronFinder search mode |
| `stats.p_adjust_method` | `BH` | Multiple testing correction method |
| `stats.alpha` | `0.05` | Significance threshold |

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

| File | Description |
|------|-------------|
| `results/colocalization_summary.csv` | Per-gene AMR–MGE co-localization table with country/region metadata |
| `results/amr_abundance_matrix.csv` | AMR class abundance matrix per sample |
| `results/mge_distribution_matrix.csv` | MGE type proportions per country |
| `results/amr_mge_association_matrix.csv` | AMR class × MGE type cross-tabulation |
| `results/figures/Fig1_Heatmap_AMR.pdf` | Heatmap: AMR abundance across samples |
| `results/figures/Fig2_MGE_distribution.pdf` | Stacked bar: MGE distribution per population |
| `results/figures/Fig3_Network_*.pdf` | Bipartite AMR–MGE network per population |
| `results/tables/Table1_KruskalWallis_AMR_Class.csv` | Kruskal-Wallis results (Sub-Q1) |
| `results/tables/Table2_Fisher_AMR_MGE.csv` | Fisher's exact results (Sub-Q3) |
| `results/tables/Supplementary_Statistics.xlsx` | All statistical outputs |

---

## Data Availability

- Metagenomic metadata sourced from [`curatedMetagenomicData`](https://bioconductor.org/packages/curatedMetagenomicData/)
- Raw sequencing data available from NCBI SRA under accession numbers listed in `data/metadata/sample_map.csv`
- Pilot cohorts: Danish gut microbiome (SRP065114), Chinese gut metagenome (SRA045646), Indian gut cohort

---

## Tools & Databases Used

| Tool | Version | Purpose | Reference |
|------|---------|---------|-----------|
| [RGI](https://github.com/arpcard/rgi) | ≥6.0 | AMR gene detection | Alcock et al., 2023 |
| [CARD](https://card.mcmaster.ca/) | ≥3.2 | AMR database | Alcock et al., 2023 |
| [MOB-suite](https://github.com/phac-nml/mob-suite) | ≥3.1 | Plasmid detection | Robertson & Nash, 2018 |
| [IntegronFinder](https://github.com/gem-pasteur/Integron_Finder) | ≥2.0 | Integron detection | Cury et al., 2016 |
| [ISEScan](https://github.com/xiezhq/ISEScan) | ≥1.7 | IS/transposon detection | Xie & Tang, 2017 |
| [MEGAHIT](https://github.com/voutcn/megahit) | ≥1.2 | Metagenome assembly (Galaxy) | Li et al., 2015 |
| [Snakemake](https://snakemake.readthedocs.io) | ≥7.0 | Workflow management | Mölder et al., 2021 |

---

## Citation

If you use this pipeline, please cite:

> Nandatama, E. (2025). *AMR Gene Mobility in Human Gut Metagenomes: Cross-Population Comparative Pipeline*. GitHub. https://github.com/engkinandatama/AMR-gene-mobility

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## Contact

**Engki Nandatama**  
GitHub: [@engkinandatama](https://github.com/engkinandatama)
