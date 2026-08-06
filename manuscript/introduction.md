# Introduction

> Draft for the AMR–MGE co-localization manuscript. Citation numbers refer to the
> shared reference list in `references.md`. See the citation-verification note at
> the top of that file before submission.

Antimicrobial resistance (AMR) is among the most consequential threats to global
health. An estimated 4.95 million deaths were associated with bacterial AMR in
2019, of which 1.27 million were directly attributable to resistant infections,
with the heaviest burden falling on sub-Saharan Africa and South Asia [38]. The
human gut harbours one of the densest microbial communities known and functions
as a substantial reservoir of antimicrobial resistance genes (ARGs), collectively
termed the gut resistome [39,44]. Because commensal gut bacteria live in close
physical association with transient and opportunistic pathogens, this reservoir
is not inert: the gut is an ecological setting in which resistance determinants
can be exchanged between community members and clinically relevant organisms
[40,44].

What converts such a reservoir into a public health risk is not the mere presence
of resistance genes but their capacity to move. The clinical relevance of an ARG
depends heavily on its genetic context — whether it is confined to a
chromosomal backbone or embedded within machinery capable of transferring it to
another cell [11]. Mobile genetic elements (MGEs) supply that machinery through
mechanistically distinct routes: conjugative plasmids transfer as autonomous
replicons between cells; integrons capture, arrange and express gene cassettes
via site-specific recombination; and insertion sequences mobilize flanking DNA,
forming composite transposons and altering expression of neighbouring genes
[13]. These elements interact, and their interplay shapes how resistance
disseminates across bacterial lineages [34]. Empirically, ARGs situated on
mobile elements spread across taxonomically distant hosts far more readily than
those that are not [46], making mobilization a central determinant of
dissemination potential.

Resistome composition is known to vary systematically with geography. Country of
residence is a strong correlate of gut resistome content, tracking national
antibiotic use practices [2], and global surveys of sewage and human populations
have documented pronounced regional structure in resistance gene profiles
[3,43]. Resistomes in low-income settings are frequently distinct from and more
diverse than those in high-income settings, reflecting differences in sanitation,
antibiotic availability and environmental exposure [41]. These patterns are
unfolding against a backdrop of rising and geographically converging antibiotic
consumption [42], and resistome composition further clusters by ecological
context rather than by taxonomy alone [47].

Despite this, population-scale resistome surveys have overwhelmingly characterized
gene content — which ARGs are present and at what abundance — while the genetic
context of those genes has rarely been assessed at comparable scale. This is
partly methodological: most large surveys quantify ARGs directly from unassembled
reads, an approach that is sensitive and computationally efficient but cannot
resolve the genomic neighbourhood of a gene, and therefore cannot determine
whether a given ARG resides on a plasmid, within an integron cassette array,
adjacent to an insertion sequence, or in the chromosomal backbone [12]. The
consequence is an asymmetry in what is known: the geography of *which* resistance
genes circulate is comparatively well described, whereas the geography of *how*
they are mobilized remains largely uncharacterized.

This gap matters because abundance and mobility are not interchangeable. Two
populations may carry statistically indistinguishable ARG loads yet differ
substantially in the transmissibility of that resistance if one population's
genes are predominantly chromosomal while the other's are plasmid-borne. Risk
assessment frameworks for resistance genes accordingly weight mobility, alongside
abundance and pathogen host range, as a primary criterion [11,45]. Resolving
mobilization routes at population scale is therefore necessary both to interpret
resistome surveillance in risk terms and to identify which vehicles — rather than
which genes — drive dissemination in a given region.

Here we address that gap using assembly-based co-localization across 248 gut
metagenomes from healthy, antibiotic-free participants sampled in five countries
spanning four geographic regions. Each metagenome was assembled independently,
high-confidence ARGs were called from the resulting contigs, and plasmids,
integrons and insertion sequences were profiled on the same contigs using three
orthogonal tools, allowing every resistance gene to be assigned to the mobile
element with which it is physically linked. We ask three questions: whether ARG
diversity and abundance differ across populations; which MGE classes
predominantly carry ARGs in each population; and whether ARG-class–MGE-type
association patterns differ geographically in a manner consistent with distinct
routes of horizontal transfer. To do so we summarize per-sample mobilization with
a Resistome Mobility Index, test whether resistome structure is decoupled from
the underlying bacterial community composition — a signature expected of
horizontally acquired rather than vertically inherited resistance — and compare
ARG–MGE association networks between populations to distinguish conserved
mobilization hubs from region-specific ones.
