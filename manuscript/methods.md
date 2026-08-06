# Methods

> Draft for the AMR–MGE co-localization manuscript. Values in `[square brackets]`
> depend on the completed run and must be filled once the pipeline finishes.
> Every parameter below reflects what the pipeline **actually executed**, verified
> against run logs and tool outputs rather than the configuration file.
>
> **Citation check required.** References were compiled from domain knowledge and
> are believed correct, but DOIs, volume/page numbers and years must be verified
> against the primary sources before submission. Do not submit without checking.

## Study design and sample selection

Publicly available human gut metagenomes were obtained through the
`curatedMetagenomicData` R package [version] [1], which provides uniformly
processed sequencing metadata and taxonomic profiles across published shotgun
metagenomic studies. Using a single curated resource rather than assembling
studies independently reduces heterogeneity in metadata definitions and
taxonomic profiling, a recognized source of batch effect in cross-study
metagenomic comparisons [1].

Samples were drawn from five countries spanning four geographic regions:
Denmark (Europe), China (East Asia), India (South Asia), and Cameroon and
Madagascar (sub-Saharan Africa). This design follows earlier evidence that the
human gut resistome varies systematically with country and regional antibiotic
use practices [2,3].

Samples were retained when they satisfied all of the following criteria: stool
body site; host disease status recorded as healthy and study condition recorded
as control; no current antibiotic use; sequencing on Illumina HiSeq or NextSeq
platforms; a minimum of 10 million reads; and an available NCBI accession.
Restricting to healthy, antibiotic-free control participants removes the two
strongest known confounders of resistome composition — disease state and recent
antibiotic exposure [2] — so that residual variation is attributable primarily
to geography. Restricting to a single sequencing chemistry avoids
platform-driven differences in error profile and read length that would
otherwise propagate into assembly quality.

To preserve statistical independence, only one sample was retained per subject,
the sample with the highest read count being kept where a subject contributed
several; treating repeated samples from one individual as independent
observations would constitute pseudoreplication [4]. Fifty samples per country
were then drawn by stratified random sampling with a fixed random seed (42) to
make the selection reproducible.

This yielded 248 unique samples for analysis (Denmark 50, China 50, India 48,
Cameroon 50, Madagascar 50), corresponding to 50 samples each for Europe, East
Asia and South Asia, and 100 for sub-Saharan Africa. Two accessions in the
Indian subset resolved to identical sample identifiers and were collapsed to a
single record.

## Read retrieval, quality control and host depletion

Raw sequencing reads were retrieved from the NCBI Sequence Read Archive [5]
using the SRA Toolkit v3.0.7 (`prefetch` followed by `fasterq-dump`). Samples
distributed across multiple sequencing runs had their runs concatenated into a
single paired-end read set prior to processing.

Reads were quality-filtered and adapter-trimmed with fastp v1.3.3 [6] using
default parameters, which apply a base quality threshold of Q15, a minimum read
length of 15 bp, and automatic adapter detection for paired-end libraries.
Filtering was minimally invasive on this dataset, consistent with its curated
origin: in a representative sample, 68,471,370 of 68,472,136 reads (>99.99%)
were retained. Because assembly graph construction is itself robust to residual
low-quality bases, and because the selection criteria already restricted the
dataset to deeply sequenced libraries, more aggressive trimming was not required.

Human host reads were removed by aligning filtered reads against the human
reference genome GRCh38/hg38 [7] with Bowtie2 v2.5.5 [8] in `--very-sensitive`
mode; read pairs failing to align concordantly to the host genome were retained
as the non-host fraction. The most sensitive preset was chosen deliberately, as
residual host sequence inflates assembly size and can generate spurious
annotation on human repetitive elements.

## Metagenomic assembly

Non-host reads were assembled per sample with MEGAHIT v1.2.9 [9] using default
parameters, corresponding to a multiple k-mer strategy over k = 21, 29, 39, 59,
79, 99, 119 and 141. MEGAHIT was selected for its succinct de Bruijn graph
implementation, which makes per-sample assembly of deeply sequenced gut
metagenomes tractable at the memory available on the compute nodes used [9].

Contigs shorter than 1,000 bp were discarded prior to annotation. This threshold
is a deliberate trade-off: MGE detection depends on flanking sequence context —
integron identification requires an integrase and adjacent *attC* sites, and
plasmid classification depends on replicon and relaxase content — so short
contigs cannot support reliable element calling, and retaining them would
increase false-negative co-localization while adding annotation noise.
Assemblies yielded a median contig length of approximately 1.8 kb
[update with final assembly statistics].

## Antimicrobial resistance gene detection

Antimicrobial resistance genes (ARGs) were predicted from assembled contigs with
the Resistance Gene Identifier (RGI) v6.0.5 against the Comprehensive Antibiotic
Resistance Database (CARD) v3.2.7 [10], run in contig mode with BLAST as the
alignment tool. CARD was chosen because its ontology is expert-curated and its
models carry gene-family-specific bitscore cutoffs, rather than relying on a
single global similarity threshold [10]. RGI's default reporting criteria were
used, retaining Perfect and Strict hits and excluding Loose hits.

Because the Strict criterion applies curated bitscore cutoffs rather than a
sequence-identity threshold, it admits distant homologs that are unlikely to
confer a resistance phenotype. In this dataset the unfiltered Strict set had a
median identity to the CARD reference of 40.1%, and low-identity predictions
were predominantly chromosomal (95.6%), a pattern consistent with detection of
conserved housekeeping homologs of resistance gene families rather than acquired
resistance determinants. This distinction between a sequence that resembles a
resistance gene and one that plausibly functions as a transferable resistance
determinant is well recognized as a central problem in resistome analysis
[11,12].

We therefore restricted the primary analysis to high-confidence ARGs, defined as
hits with at least 95% identity to the CARD reference sequence and at least 95%
coverage of the reference length. Comparable identity and coverage thresholds
are standard practice when the analytical goal is acquired, potentially mobile
resistance genes rather than the full space of resistance-related homologs
[12]. Raw RGI output was retained so that the permissive (unfiltered) resistome
can be reported as a sensitivity analysis [see Supplementary Table S1].

## Mobile genetic element profiling

Three classes of mobile genetic element (MGE) were profiled independently on the
same contigs, chosen because plasmids, integrons and insertion sequences
represent the principal, mechanistically distinct vehicles of resistance gene
mobilization in bacteria [13]:

- **Plasmids.** MOB-suite v3.1.9 (`mob_recon`) [14,15] was used to reconstruct
  and type plasmid sequences; contigs classified as `plasmid` were treated as
  plasmid-derived, and replicon and relaxase types were recorded. MOB-suite was
  selected because it types plasmids by replicon and relaxase content and
  assigns cluster identity by whole-sequence similarity, which is more robust on
  fragmented metagenomic assemblies than replicon detection alone [15].
- **Integrons.** IntegronFinder v2.0.6 [16,17] was run in `--local-max` mode,
  which improves *attC* site detection sensitivity in fragmented sequence at the
  cost of runtime. Contigs carrying complete integrons, CALIN elements (clusters
  of *attC* sites lacking a neighbouring integrase) or In0 elements were
  recorded, together with annotated gene cassettes. CALIN elements were retained
  separately as evidence of historical integron activity on contigs where the
  integrase was not co-assembled.
- **Insertion sequences and transposons.** ISEScan v1.7.3 [18] was run with
  default parameters to identify contigs carrying insertion sequence elements.

## AMR–MGE co-localization

Co-localization was defined at the contig level: an ARG was considered
associated with an MGE when the contig carrying it was independently identified
as plasmid-derived, integron-bearing or IS-bearing. Physical linkage on a single
assembled fragment is the standard assembly-based proxy for mobilization
potential, since it places the resistance gene within the genetic context of the
element [13].

Because the MGE tools and RGI report contig identifiers in different formats —
MOB-suite retains the full assembler FASTA header, whereas RGI and
IntegronFinder report only the leading identifier — all contig identifiers were
normalized to the first whitespace-delimited token before matching. To guard
against silent identifier-namespace mismatches, the pipeline emits a per-sample
diagnostic recording, for each MGE tool, the number and proportion of ARG-bearing
contigs matched, and flags samples in which a tool reported contigs but none
intersected the ARG contig set.

Each ARG was assigned to a single, mutually exclusive category using the
hierarchy Plasmid > Integron > IS/Transposon > Chromosomal, so that an ARG on a
contig carrying several element types is attributed to the highest-ranked class.
The ordering reflects decreasing capacity for autonomous inter-cellular
transfer: plasmids can be conjugatively transferred as replicons, integrons
capture and express cassettes but depend on a vector for cell-to-cell movement,
and insertion sequences mobilize sequence intracellularly [13]. ARGs on contigs
with no detected MGE were classified as chromosomal, a conservative assignment
given the assembly limitations discussed below.

A per-sample **Resistome Mobility Index** was defined as the proportion of
high-confidence ARGs assigned to any non-chromosomal category, providing a
single sample-level summary of resistome mobilization that is comparable across
samples with differing total ARG loads.

## Statistical analysis

All statistical analyses were performed in R v4.3.2 [19]. Unless stated
otherwise, p-values were adjusted for multiple testing by the
Benjamini–Hochberg false discovery rate procedure [20] and evaluated at
α = 0.05, and permutation tests used 999 permutations.

**Resistome diversity.** ARG abundance matrices were constructed per sample at
the drug-class level. Alpha diversity was summarized by the Shannon index and
observed richness using vegan [21]. Beta diversity was computed as Bray–Curtis
dissimilarity [22], which is appropriate for abundance data in which shared
absences should not contribute to similarity, and tested against country of
origin by PERMANOVA (`adonis2`) [23]. Because PERMANOVA can return significant
results when groups differ in within-group dispersion rather than in location,
the homogeneity assumption was tested explicitly with PERMDISP
(`betadisper` followed by `permutest`) [24], and pairwise PERMANOVA was used for
post-hoc comparisons between countries. A multi-factor PERMANOVA additionally
including participant age and sex was fitted to assess whether geographic
differences persisted after adjustment for these covariates.

**Resistome versus community composition.** To distinguish resistome variation
that simply tracks bacterial community structure from variation that does not,
species-level relative abundance profiles (MetaPhlAn [25], obtained through
`curatedMetagenomicData`) were compared with ARG profiles by Mantel test
(Spearman correlation on Bray–Curtis distances) [26] and by Procrustes analysis
on the corresponding principal coordinate ordinations, the latter providing
greater power and a more interpretable measure of ordination concordance than
the Mantel test alone [27]. A dual-Mantel decoupling analysis then correlated
taxonomic distance separately with the mobile and the non-mobile (chromosomal)
fractions of the resistome, on the rationale that vertically inherited
resistance should track taxonomy more closely than horizontally mobilized
resistance; divergence between the two correlations is therefore interpretable
as a signature of horizontal transfer.

**Differential abundance.** Differences in ARG abundance between countries were
tested per drug class by Kruskal–Wallis test [28], a rank-based test chosen
because ARG count distributions are strongly right-skewed and zero-inflated. The
sample-by-drug-class matrix was completed with structural zeros before testing,
so that samples lacking a given drug class contribute zero counts rather than
being silently dropped as missing observations.

**Resistome mobility.** The Resistome Mobility Index was compared across
countries by Kruskal–Wallis test [28] with Dunn's post-hoc test [29]. As the
primary inferential test of whether ARG mobility differs geographically, a
generalized linear mixed model was fitted with lme4 v1.1.37 [30], modelling the
gene-level binary outcome (mobile versus chromosomal) as a function of country
with a random intercept for sample: `mobile ~ Country + (1 | Sample_ID)`. The
random effect accounts for clustering of multiple ARGs within a sample; pooling
gene-level observations across samples and treating them as independent would
constitute pseudoreplication and produce anticonservative p-values [4,31]. The
overall effect of country was assessed by likelihood-ratio test against the
corresponding intercept-only model. Pooled per-drug-class Fisher exact tests and
a chi-square test of MGE-type distribution are additionally reported as
descriptive, exploratory summaries only, and are explicitly not interpreted
inferentially, for the same reason.

## Network analysis

Bipartite co-occurrence networks linking ARGs to MGE types were constructed per
country from mobile ARGs, with edges weighted by the number of co-occurrences,
using igraph [32]. Degree and betweenness centrality were computed to identify
putative hub ARGs — genes associated with multiple element classes, and
therefore with the broadest apparent mobilization potential — and network
density and topology were compared between countries. Structural similarity
between country networks was quantified by the Jaccard index of their edge sets,
which measures how far the same ARG–element associations recur across
populations. Networks were visualized with ggraph [33].

Plasmid replicon content was additionally analysed across regions to identify
replicon types detected in two or more geographic regions, and integron gene
cassettes were classified as region-specific, multi-regional or universal
according to the number of regions in which they were detected. Both analyses
address whether particular mobilization vehicles are geographically widespread
or locally restricted [3,34].

## Reproducibility and code availability

The complete workflow was implemented in Snakemake v7.32.4 [35] with per-rule
Conda environments and executed on an HPC cluster under SLURM. Pinning each rule
to an isolated environment, and recording the resolved versions of every tool and
database, allows the analysis to be re-executed with identical software. The
workflow definition, environment specifications, analysis scripts and exact
software versions are available at [repository URL / DOI]; exact versions of all
tools and databases are additionally listed in `ENV_VERSIONS.md`.

## Limitations to state in the Discussion

- **Short-read fragmentation of repetitive elements.** All libraries were
  Illumina short-read (224 HiSeq, 26 NextSeq). Insertion sequences (~0.7–2.5 kb)
  exceed the read and typical insert length, so de Bruijn graph assembly breaks
  contigs at repeat boundaries; the difficulty of reconstructing plasmids and
  mobile elements from short-read data alone is well documented [36,37]. ARGs
  physically flanked by IS elements may therefore be assembled onto contigs
  lacking the element and conservatively classified as chromosomal. Reported MGE
  co-localization — particularly for IS/transposons — is a lower bound. Because
  the bias is systematic across all samples, comparative conclusions remain
  valid, but absolute mobilization fractions are under-estimated. Long-read or
  hybrid assembly would be required to recover full ARG–MGE linkage [37].
- **Contig-level co-localization** establishes physical linkage on the same
  assembled fragment but does not by itself demonstrate that the element
  mobilizes the ARG; functional or conjugation assays would be required.
- **Uneven geographic representation**, with sub-Saharan Africa represented by
  two countries and the remaining regions by one each.
- **Cross-sectional, publicly archived data**, precluding inference about
  temporal dynamics or causal drivers such as national antibiotic consumption,
  which could only be addressed by linkage to consumption data [2].

---

**References:** see the shared, continuously numbered reference list in
`references.md`. BibTeX entries are in `references.bib`.
