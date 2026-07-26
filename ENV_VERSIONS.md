# Software and database versions

Exact versions of every tool, database and language runtime used to produce the
`full_run` results. Captured directly from the Conda environments that processed
the dataset, so this is the ground truth for the manuscript Methods section and
for reproducing the analysis.

## Tools

| Tool | Version | Role |
|------|---------|------|
| SRA Toolkit (prefetch/fasterq-dump) | 3.0.7 | Read download |
| fastp | 1.3.3 | Read QC / trimming |
| Bowtie2 | 2.5.5 | Host (human) read removal |
| MEGAHIT | 1.2.9 | Metagenome assembly |
| RGI | 6.0.5 | AMR gene detection |
| MOB-suite (mob_recon) | 3.1.9 | Plasmid detection / replicon typing |
| IntegronFinder | 2.0.6 | Integron detection |
| ISEScan | 1.7.3 | Insertion sequence / transposon detection |

## Databases

| Database | Version | Notes |
|----------|---------|-------|
| CARD (via RGI) | 3.2.7 | Antibiotic resistance ontology used by RGI |
| Human reference | hg38 | UCSC goldenPath `hg38.fa.gz`, indexed with Bowtie2 |

## Language runtimes and key libraries

Python environment (co-localization / aggregation / deep analysis scripts):

- Python 3.11.15
- pandas 3.0.3

R environment (statistics / networks — `envs/r_stats.yaml`):

- R 4.3.2
- vegan 2.7.1
- lme4 1.1.37 (binomial GLMM for AMR mobility)
- igraph 2.1.4
- ggraph 2.2.1
- dunn.test 1.3.7
- ComplexHeatmap 2.18.0
- plus standard tidyverse packages (dplyr, tidyr, ggplot2, tibble), writexl, RColorBrewer, scales, readr

## Pinning the Conda environments (apply AFTER the run finishes)

The `envs/*.yaml` files are intentionally left partly unpinned during the run:
editing them changes each environment's hash, which would make Snakemake rebuild
the environment and could process the remaining samples with a different tool
version than the ones already completed. Once `snakemake ... all` reports
`Nothing to be done`, pin the versions below so a future clean rebuild is exactly
reproducible:

```yaml
# envs/sra-tools.yaml
  - sra-tools=3.0.7
# envs/fastp.yaml
  - fastp=1.3.3
# envs/bowtie2.yaml
  - bowtie2=2.5.5
# envs/megahit.yaml
  - megahit=1.2.9
# envs/rgi.yaml
  - rgi=6.0.5
# envs/mobsuite.yaml
  - mob_suite=3.1.9
# envs/integronfinder.yaml
  - integron_finder=2.0.6
# envs/isescan.yaml
  - isescan=1.7.3
# envs/python.yaml
  - python=3.11.15
  - pandas=3.0.3
# envs/r_stats.yaml
  - r-base=4.3.2
  - r-lme4=1.1.37
  - r-vegan=2.7.1
  # (other R packages already resolve consistently under r-base=4.3.2)
```

Note: the CARD database version (3.2.7) is not set by the Conda YAML — it is the
database RGI loaded at build time. Record it explicitly in the Methods section.
