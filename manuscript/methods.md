# Methods

> Draft for the AMR–MGE co-localization manuscript. Values in `[square brackets]`
> depend on the completed run and must be filled once the pipeline finishes.
> Every parameter below reflects what the pipeline **actually executed**, verified
> against run logs and tool outputs rather than the configuration file.

## Study design and sample selection

Publicly available human gut metagenomes were obtained through the
`curatedMetagenomicData` R package [version], which provides uniformly processed
metadata and taxonomic profiles for published shotgun metagenomic studies. We
selected samples from five countries spanning four geographic regions: Denmark
(Europe), China (East Asia), India (South Asia), and Cameroon and Madagascar
(sub-Saharan Africa).

Samples were retained when they satisfied all of the following criteria: stool
body site; host disease status recorded as healthy and study condition recorded
as control; no current antibiotic use; sequencing on Illumina HiSeq or NextSeq
platforms; at least 10 million reads; and an available NCBI accession. To
preserve statistical independence, only one sample was retained per subject, the
sample with the highest read count being kept where a subject contributed
several. Fifty samples per country were then drawn by stratified random sampling
with a fixed random seed (42).

This yielded 248 unique samples for analysis (Denmark 50, China 50, India 48,
Cameroon 50, Madagascar 50), corresponding to 50 samples each for Europe, East
Asia and South Asia, and 100 for sub-Saharan Africa. Two accessions in the Indian
subset resolved to identical sample identifiers and were collapsed to a single
record.

## Read retrieval, quality control and host depletion

Raw sequencing reads were retrieved from the NCBI Sequence Read Archive using the
SRA Toolkit v3.0.7 (`prefetch` followed by `fasterq-dump`). Samples distributed
across multiple sequencing runs had their runs concatenated into a single
paired-end read set.

Reads were quality-filtered and adapter-trimmed with fastp v1.3.3 using default
parameters, which apply a base quality threshold of Q15, a minimum read length of
15 bp, and automatic adapter detection for paired-end libraries. Filtering was
minimally invasive on this dataset, consistent with the curated origin of the
data: in a representative sample, 68,471,370 of 68,472,136 reads (>99.99%) were
retained.

Human host reads were removed by aligning the filtered reads against the human
reference genome hg38 (UCSC `goldenPath` assembly) with Bowtie2 v2.5.5 in
`--very-sensitive` mode. Read pairs failing to align concordantly to the host
genome were retained as the non-host fraction.

## Metagenomic assembly

Non-host reads were assembled per sample with MEGAHIT v1.2.9 using default
parameters, corresponding to a multiple k-mer strategy over k = 21, 29, 39, 59,
79, 99, 119 and 141. Contigs shorter than 1,000 bp were discarded prior to
downstream annotation. Assemblies yielded a median contig length of
approximately 1.8 kb [update with final assembly statistics].

## Antimicrobial resistance gene detection

Antimicrobial resistance genes (ARGs) were predicted from assembled contigs with
the Resistance Gene Identifier (RGI) v6.0.5 against the Comprehensive Antibiotic
Resistance Database (CARD) v3.2.7, run in contig mode with BLAST as the alignment
tool. RGI's default reporting criteria were used, which retain Perfect and Strict
hits and exclude Loose hits.

Because RGI's Strict criterion applies curated bitscore cutoffs rather than a
sequence-identity threshold, it admits distant homologs that are unlikely to
confer resistance. In this dataset the unfiltered Strict set had a median
identity to the CARD reference of 40.1%, and these low-identity predictions were
predominantly chromosomal (95.6%), consistent with detection of conserved
housekeeping homologs of resistance gene families rather than acquired
resistance determinants.

We therefore restricted the primary analysis to high-confidence ARGs, defined as
hits with at least 95% identity to the CARD reference sequence and at least 95%
coverage of the reference length. Raw RGI output was retained so that the
permissive (unfiltered) resistome can be reported as a sensitivity analysis
[see Supplementary Table S1].

## Mobile genetic element profiling

Three classes of mobile genetic element (MGE) were profiled independently on the
same contigs:

- **Plasmids.** MOB-suite v3.1.9 (`mob_recon`) was used to reconstruct and type
  plasmid sequences; contigs classified as `plasmid` were treated as
  plasmid-derived, and replicon and relaxase types were recorded.
- **Integrons.** IntegronFinder v2.0.6 was run in `--local-max` mode. Contigs
  carrying complete integrons, CALIN elements (clusters of attC sites lacking a
  neighbouring integrase) or In0 elements were recorded, together with the
  annotated gene cassettes.
- **Insertion sequences and transposons.** ISEScan v1.7.3 was run with default
  parameters to identify contigs carrying insertion sequence elements.

## AMR–MGE co-localization

Co-localization was defined at the contig level: an ARG was considered associated
with an MGE when the contig carrying it was independently identified as
plasmid-derived, integron-bearing or IS-bearing. Because the three MGE tools and
RGI report contig identifiers in different formats — MOB-suite retains the full
assembler FASTA header whereas RGI and IntegronFinder report only the leading
identifier — all contig identifiers were normalized to the first
whitespace-delimited token before matching.

Each ARG was assigned to a single, mutually exclusive category using the
hierarchy Plasmid > Integron > IS/Transposon > Chromosomal, so that an ARG on a
contig carrying multiple element types is attributed to the highest-ranked class.
ARGs on contigs with no detected MGE were classified as chromosomal, which is a
conservative assignment (see Limitations).

To guard against silent identifier-namespace mismatches between tools, the
pipeline emits a per-sample diagnostic recording, for each MGE tool, the number
and proportion of ARG-bearing contigs matched, flagging samples in which a tool
reported contigs but none intersected the ARG contig set.

A per-sample **Resistome Mobility Index** was defined as the proportion of
high-confidence ARGs assigned to any non-chromosomal category.

## Statistical analysis

All statistical analyses were performed in R v4.3.2. Unless stated otherwise,
p-values were adjusted for multiple testing with the Benjamini–Hochberg
procedure and evaluated at α = 0.05, and permutation tests used 999 permutations.

**Resistome diversity.** ARG abundance matrices were constructed per sample at
the drug-class level. Alpha diversity was summarized by the Shannon index and
observed richness (vegan v2.7.1). Beta diversity was computed as Bray–Curtis
dissimilarity and tested against country of origin by PERMANOVA (`adonis2`).
Because PERMANOVA is sensitive to heterogeneous within-group dispersion, the
homogeneity assumption was tested explicitly with PERMDISP (`betadisper`
followed by `permutest`), and pairwise PERMANOVA was used for post-hoc
comparisons between countries. A multi-factor PERMANOVA additionally including
participant age and sex was fitted to assess whether geographic differences
persisted after adjustment for these covariates.

**Resistome versus community composition.** To distinguish resistome variation
that tracks bacterial community structure from variation that does not,
species-level relative abundance profiles (MetaPhlAn, obtained through
`curatedMetagenomicData`) were compared with ARG profiles using a Mantel test
(Spearman correlation on Bray–Curtis distances) and Procrustes analysis on the
corresponding principal coordinate ordinations. A dual-Mantel decoupling analysis
was then performed, correlating taxonomic distance separately with the mobile and
the non-mobile (chromosomal) fractions of the resistome, on the rationale that
vertically inherited resistance should track taxonomy more closely than
horizontally mobilized resistance.

**Differential abundance.** Differences in ARG abundance between countries were
tested per drug class by Kruskal–Wallis test. The sample-by-drug-class matrix was
completed with structural zeros before testing, so that samples lacking a given
drug class contribute zero counts rather than being dropped as missing.

**Resistome mobility.** The Resistome Mobility Index was compared across
countries by Kruskal–Wallis test with Dunn's post-hoc test. As the primary
inferential test of whether ARG mobility differs geographically, a
generalized linear mixed model was fitted with lme4 v1.1.37, modelling the
gene-level binary outcome (mobile versus chromosomal) as a function of country
with a random intercept for sample: `mobile ~ Country + (1 | Sample_ID)`. The
random effect accounts for the clustering of multiple ARGs within a sample, which
a pooled gene-level contingency test would treat as independent observations. The
overall effect of country was assessed by likelihood-ratio test against the
corresponding intercept-only model. Pooled per-drug-class Fisher exact tests and
a chi-square test of MGE-type distribution are additionally reported as
descriptive, exploratory summaries only, and are not interpreted inferentially.

## Network analysis

Bipartite co-occurrence networks linking ARGs to MGE types were constructed per
country from mobile ARGs, with edges weighted by the number of co-occurrences
(igraph v2.1.4). Degree and betweenness centrality were computed to identify
putative hub ARGs, and network density and topology were compared between
countries. Structural similarity between country networks was quantified by the
Jaccard index of their edge sets. Networks were visualized with ggraph v2.2.1.

Plasmid replicon content was additionally analysed across regions to identify
replicon types detected in two or more geographic regions, and integron gene
cassettes were classified as region-specific, multi-regional or universal
according to the number of regions in which they were detected.

## Reproducibility and code availability

The complete workflow was implemented in Snakemake v7.32.4 with per-rule Conda
environments and executed on an HPC cluster under SLURM. The workflow definition,
environment specifications, analysis scripts and exact software versions are
available at [repository URL / DOI]. Exact versions of all tools and databases
are listed in `ENV_VERSIONS.md`.

## Limitations to state in the Discussion

- **Short-read fragmentation of repetitive elements.** All libraries were Illumina
  short-read (224 HiSeq, 26 NextSeq). Insertion sequences (~0.7–2.5 kb) exceed the
  read and typical insert length, so de Bruijn graph assembly breaks contigs at
  repeat boundaries. ARGs physically flanked by IS elements may therefore be
  assembled onto contigs lacking the element and conservatively classified as
  chromosomal. Reported MGE co-localization — particularly for IS/transposons —
  is a lower bound. Because the bias is systematic across all samples, comparative
  conclusions remain valid, but absolute mobilization fractions are
  under-estimated. Long-read or hybrid assembly would be required to recover the
  full ARG–MGE linkage.
- **Contig-level co-localization** establishes physical linkage on the same
  assembled fragment but does not by itself demonstrate that an element mobilizes
  the ARG.
- **Uneven geographic representation**, with sub-Saharan Africa represented by two
  countries and the remaining regions by one each.
- **Cross-sectional, publicly archived data**, precluding inference about
  temporal dynamics or causal drivers such as antibiotic consumption.
