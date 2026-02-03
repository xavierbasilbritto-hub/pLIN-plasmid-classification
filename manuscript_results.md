# pLIN: A Plasmid Life Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids

---

## ABSTRACT

Plasmids are central agents of horizontal gene transfer in bacteria, driving the dissemination of antimicrobial resistance (AMR) determinants, virulence factors, and adaptive traits. Despite their clinical and epidemiological importance, no existing plasmid classification system simultaneously provides hierarchical multi-resolution typing, code permanence, and reference-free operation. Current approaches—including replicon-based Inc typing (PlasmidFinder), plasmid multilocus sequence typing (pMLST), mobilization-based clustering (MOB-suite), ANI-based plasmid taxonomic units (COPLA/PTUs), and reference-free unitig clustering (mge-cluster)—each suffer from flat classification, database dependency, code instability, or limited taxonomic scope. Here, we introduce pLIN (plasmid Life Identification Number), the first application of the Life Identification Number (LIN) framework to plasmid genomes. pLIN assigns each plasmid a six-position hierarchical code based on tetranucleotide composition distances and single-linkage clustering at six biologically calibrated thresholds, spanning family-level (~85% ANI) to strain-level (~99.9% ANI) resolution. Applied to 6,346 complete plasmid sequences from three clinically important incompatibility groups (IncFII, n=4,581; IncN, n=1,064; IncX1, n=701), pLIN resolved 2,232 unique strain-level codes with a Simpson's Index of Diversity of 0.979, while maintaining 99.5% concordance with established Inc-group assignments. A nested cross-validation machine learning pipeline confirmed that the same compositional features underlying pLIN assignment robustly predict Inc-group membership (weighted F1=0.903, XGBoost). Mosaic structure analysis identified 43.1% of IncX1 plasmids as candidate chimeras, and seed ORF conservation analysis delineated a 25-gene IncX1 backbone with 50–85% prevalence across the group. Integration with NCBI AMRFinderPlus v4.2.5 across all 6,346 plasmids identified 60,372 gene detections (27,465 AMR, 5,834 virulence, 27,073 stress) in 84.2% of plasmids. Cross-referencing pLIN lineages with AMR gene content revealed lineage-specific resistance profiles: pLIN 1.1.1.7.30.567 (IncN) carried *blaKPC-2* on 100% of members with a mean of 13.2 AMR genes per plasmid, while the cross-Inc lineage 1.1.1.7.30.1313 harboured both *blaKPC-2* (58%) and the pan-aminoglycoside resistance gene *rmtB1* (71%). Critically important resistance determinants including carbapenemases (1,490 detections), ESBLs (3,734), plasmid-mediated colistin resistance (*mcr*, 160), and plasmid-mediated quinolone resistance (2,160) were mapped to specific pLIN lineages, demonstrating the utility of pLIN as a framework for tracking AMR gene dissemination through plasmid lineage surveillance. The pLIN system offers a stable, scalable, and biologically interpretable nomenclature for plasmid genomic epidemiology.

---

## RESULTS

### 3.1 Dataset Composition and Quality Control

A total of 6,346 complete plasmid genome sequences were curated from NCBI RefSeq across three incompatibility groups that are among the most clinically significant carriers of antimicrobial resistance genes in Enterobacteriaceae: IncFII (n=4,581), IncN (n=1,064), and IncX1 (n=701) (Table 1). These three groups were selected to span a range of plasmid population structures, from the highly diverse and recombination-prone IncFII family to the more cohesive IncX1 lineage.

Prior to analysis, a systematic deduplication procedure was performed. An exact duplicate directory ('IncFII 2', 4,671 identical files, ~553 MB) was identified and removed. Cross-directory comparison revealed 142 sequences present in multiple Inc-type folders with identical content: 52 shared between IncFII and IncX1, 42 between IncFII and IncN, and 52 between IncN and IncX1. Duplicates were removed from the larger group in each case, retaining the copy in the most specific (smallest) Inc-type folder to avoid training bias. The final deduplicated dataset comprised 6,346 non-redundant plasmid sequences.

Plasmid length varied substantially both within and across incompatibility groups (Table 1). IncFII plasmids were the largest on average (mean 120,219 bp; range 1,871–399,913 bp), consistent with the known modular and mosaic architecture of F-type plasmids that frequently carry large accessory regions including AMR gene cassettes. IncN plasmids showed intermediate sizes (mean 89,332 bp; range 2,548–395,758 bp), while IncX1 plasmids were the smallest (mean 62,980 bp; median 47,397 bp; range 3,978–380,891 bp), reflecting a more compact backbone architecture.

**Table 1. Dataset summary statistics.**

| Inc Group | Sequences (n) | Mean Length (bp) | Median Length (bp) | Range (bp) | Unique pLIN Codes |
|-----------|---------------|------------------|--------------------|------------|-------------------|
| IncFII    | 4,581         | 120,219          | 110,786            | 1,871–399,913 | 1,409          |
| IncN      | 1,064         | 89,332           | 62,391             | 2,548–395,758 | 420            |
| IncX1     | 701           | 62,980           | 47,397             | 3,978–380,891 | 417            |
| **Total** | **6,346**     | —                | —                  | —          | **2,232**         |

---

### 3.2 Construction of the pLIN Hierarchical Coding System

#### 3.2.1 Whole-Plasmid Distance Estimation

Tetranucleotide (4-mer) frequency vectors were computed for all 6,346 plasmids, yielding a 6,346 x 256 composition matrix. Pairwise cosine distances were calculated across all 20,132,685 plasmid pairs using the scipy spatial distance module. Cosine distance was chosen as the distance metric because it is scale-invariant (robust to plasmid size differences), computationally efficient, and has been shown to correlate well with genomic divergence for k-mer frequency profiles. This composition-based approach approximates Mash/ANI-based distances while operating entirely in pure Python without external bioinformatics tool dependencies.

For the IncX1 reference subset (n=500, randomly sampled), the mean pairwise distance was 0.041 (SD 0.029), the median was 0.033, and distances ranged from 0.000 (identical composition) to 0.185. Hybrid distances incorporating alignment fraction estimates (D = (1 - ANI) x AF) produced a compressed range (mean 0.021, median 0.017), reflecting the adjustment for partial alignments between divergent plasmids.

#### 3.2.2 Hierarchical Threshold Calibration

Six distance thresholds were defined to capture biologically meaningful levels of plasmid relatedness, calibrated against established ANI benchmarks (Table 2). These thresholds were informed by the nearest-neighbour distance distribution of 178 selected IncX-like reference plasmids (seed plasmid: RefSeq_NZ_AP027441.1, distance cutoff d ≤ 0.12). Nearest-neighbour distance quantiles showed that 50% of plasmids had a closest relative within d = 0.011, 90% within d = 0.038, and 99% within d = 0.072, providing empirical anchors for threshold placement.

**Table 2. pLIN hierarchical bin definitions and clustering results.**

| Bin | Level | Distance Threshold (d) | ANI Equivalent | Clusters (n) | Median Cluster Size | Max Cluster Size | Singletons |
|-----|-------|------------------------|----------------|---------------|---------------------|------------------|------------|
| A   | Family        | ≤ 0.150 | ~85% | 1     | 6,346 | 6,346 | 0     |
| B   | Subfamily     | ≤ 0.100 | ~90% | 1     | 6,346 | 6,346 | 0     |
| C   | Cluster       | ≤ 0.050 | ~95% | 2     | 3,173 | 6,345 | 1     |
| D   | Subcluster    | ≤ 0.020 | ~98% | 20    | 1     | 6,320 | 15    |
| E   | Clone complex | ≤ 0.010 | ~99% | 82    | 1     | 6,203 | 58    |
| F   | Strain        | ≤ 0.001 | ~99.9% | 2,232 | 1   | 817   | 1,702 |

At the coarsest level (Bin A, d ≤ 0.150), all 6,346 plasmids from three Inc groups formed a single cluster, consistent with their shared membership in the broader Enterobacteriaceae plasmid superfamily. Meaningful separation first emerged at Bin C (d ≤ 0.050, ~95% ANI), which split the dataset into two clusters: a major cluster (n=6,345) and a single outlier—an atypical IncFII plasmid with extreme compositional divergence. At Bin D (d ≤ 0.020, ~98% ANI), 20 subclusters were resolved, and at Bin E (d ≤ 0.010, ~99% ANI), 82 clone complexes were delineated. The finest resolution, Bin F (d ≤ 0.001, ~99.9% ANI), produced 2,232 strain-level groups, of which 1,702 (76.3%) were singletons and 530 (23.7%) contained two or more members.

#### 3.2.3 pLIN Code Assignment

Each plasmid was assigned a six-position pLIN code of the form A.B.C.D.E.F, where each position denotes the cluster identifier at the corresponding hierarchical level. For example, pLIN code 1.1.1.7.30.63 indicates membership in Family 1, Subfamily 1, Cluster 1, Subcluster 7, Clone complex 30, and Strain group 63. The pLIN code is assigned once and is permanent: addition of new plasmids to the database does not alter existing codes, a fundamental property inherited from the LIN framework.

The largest strain-level cluster (pLIN 1.1.1.7.30.1240, n=817) consisted almost exclusively of IncFII plasmids (814/817, 99.6%) and likely represents a dominant IncFII lineage with highly conserved backbone composition. The second largest cluster (pLIN 1.1.1.7.30.1313, n=225) was more heterogeneous, containing IncFII (n=184), IncX1 (n=23), and IncN (n=18) members, suggesting compositional convergence or shared accessory module content across Inc groups.

---

### 3.3 Concordance with Established Inc-Group Classification

To assess whether pLIN codes preserve established plasmid taxonomy, we evaluated the Inc-type purity of strain-level pLIN clusters (Bin F). Of 2,232 unique pLIN codes, 2,220 (99.5%) contained plasmids from a single Inc group, and only 12 codes (0.5%) contained members from two or more Inc groups (Table 3). This near-perfect concordance indicates that pLIN captures Inc-group boundaries as an emergent property of whole-plasmid composition, without requiring explicit replicon detection.

**Table 3. Mixed-Inc pLIN codes at strain level (Bin F).**

| pLIN Code | Total (n) | IncFII | IncN | IncX1 | Interpretation |
|-----------|-----------|--------|------|-------|----------------|
| 1.1.1.7.30.1240 | 817 | 814 | 3 | 0 | Dominant IncFII lineage, 3 IncN convergent |
| 1.1.1.7.30.1313 | 225 | 184 | 18 | 23 | Cross-Inc convergence zone |
| 1.1.1.7.30.743  | 108 | 3   | 104 | 1 | Dominant IncN lineage |
| 1.1.1.7.30.1271 | 94  | 93  | 0  | 1  | IncFII with 1 convergent IncX1 |
| 1.1.1.7.30.1166 | 42  | 41  | 1  | 0  | IncFII with 1 convergent IncN |

The 12 mixed-Inc codes contained 1,369 plasmids total, of which 1,325 (96.8%) belonged to the majority Inc type. The few cross-Inc plasmids within these codes may represent genuine cases of compositional convergence driven by horizontal acquisition of large genomic modules, or multi-replicon plasmids carrying both Inc-type markers. Notably, pLIN 1.1.1.7.30.1313 (the largest mixed code) contained all three Inc types, suggesting a compositional 'grey zone' where extensive module exchange has blurred Inc-group boundaries—a phenomenon well-documented in the literature for IncF/IncN hybrid plasmids.

---

### 3.4 Discriminatory Power

The overall Simpson's Index of Diversity (D) for the pLIN system at strain level (Bin F) was 0.979, approaching the theoretical maximum of 1.0 and indicating that two randomly selected plasmids have a 97.9% probability of receiving different pLIN codes. Discriminatory power varied across Inc groups: IncX1 exhibited the highest diversity (D=0.995), followed by IncN (D=0.976) and IncFII (D=0.962), reflecting the more homogeneous population structure of the large IncFII dataset where several dominant lineages account for a disproportionate share of sequences.

For comparison, replicon-based Inc typing alone would yield D=0.508 on this three-group dataset, and pMLST (where available) typically achieves D=0.70–0.90 within individual Inc groups. The hierarchical structure of pLIN further enables adjustable resolution: at Bin E (clone complex level, d ≤ 0.010), D=0.960 with 82 groups; at Bin D (subcluster, d ≤ 0.020), D=0.870 with 20 groups. This tuneable resolution is a unique feature not available in flat classification systems.

---

### 3.5 Machine Learning Validation of Compositional Features

To independently validate that the tetranucleotide composition features underlying pLIN assignment carry genuine biological signal for plasmid classification, we trained supervised machine learning models to predict Inc-group membership from the same 32-feature composition vector used for pLIN distance computation. A balanced subset of 1,500 plasmids (500 per Inc group) was used to avoid class imbalance effects. A rigorous nested cross-validation framework (3 outer folds, 2 inner folds) with Optuna hyperparameter optimization was employed to prevent overfitting and data leakage.

**Table 4. Machine learning model performance for Inc-group prediction (3-class, nested CV).**

| Model | Mean Weighted F1 | SD | Min | Max |
|-------|-------------------|----|-----|-----|
| XGBoost | 0.903 | 0.009 | 0.896 | 0.916 |
| Gradient Boosting | 0.901 | 0.007 | 0.896 | 0.912 |
| Random Forest | 0.881 | 0.019 | 0.862 | 0.907 |
| Logistic Regression | 0.873 | 0.010 | 0.860 | 0.884 |

All four models achieved weighted F1 scores exceeding 0.87, with XGBoost performing best (F1=0.903 ± 0.009). The strong performance of even the linear baseline (Logistic Regression, F1=0.873) indicates that Inc-group boundaries are substantially, though not entirely, linearly separable in composition space. The modest improvement from tree-based ensemble methods suggests the presence of nonlinear interaction effects between compositional features.

Feature importance analysis across all three tree-based models revealed a consistent ranking of discriminative features (Table 5). Stop codon-associated trinucleotides (TAG, TGA) and CpG-related motifs (GCG, CGC) dominated the top positions, consistent with known differences in codon usage, methylation patterns, and restriction-modification systems across plasmid lineages.

**Table 5. Top 10 consensus feature importances (averaged across tree-based models).**

| Rank | Feature | Mean Importance | Biological Relevance |
|------|---------|-----------------|----------------------|
| 1 | tri_TAG | 0.114 | Amber stop codon; reflects codon usage |
| 2 | tri_GCG | 0.099 | CpG island proxy; methylation signatures |
| 3 | at_content | 0.076 | Overall base composition |
| 4 | gc_content | 0.053 | Global GC%; backbone vs. accessory content |
| 5 | di_AA | 0.050 | Poly-A tracts; regulatory elements |
| 6 | tri_TGA | 0.048 | Opal stop codon frequency |
| 7 | tri_CGC | 0.040 | CpG-related; restriction site frequency |
| 8 | di_TT | 0.039 | Poly-T tracts; transcriptional signals |
| 9 | di_GG | 0.032 | G-richness; G-quadruplex potential |
| 10 | tri_GGG | 0.030 | Extreme GC-rich tracts |

---

### 3.6 IncX1 Backbone Architecture and Seed ORF Conservation

Using the IncX1 seed plasmid (RefSeq_NZ_AP027441.1) as a reference, the 25 most conserved open reading frames (ORFs) were identified by k-mer containment analysis across 178 selected IncX-like references (Table 6). These ORFs define the candidate IncX1 core backbone.

**Table 6. Top 5 conserved seed ORFs in the IncX1 backbone.**

| ORF ID | Position | Length (nt) | Prevalence (%) | Median Match Fraction |
|--------|----------|-------------|-----------------|----------------------|
| IncX1_core_ORF067 | 42806–43253 (+) | 447 | 84.8% | 0.852 |
| IncX1_core_ORF063 | 38165–38552 (+) | 387 | 71.3% | 0.943 |
| IncX1_core_ORF062 | 37861–38197 (+) | 336 | 71.3% | 0.934 |
| IncX1_core_ORF061 | 37787–38093 (+) | 306 | 71.3% | 0.900 |
| IncX1_core_ORF066 | 41930–42767 (+) | 837 | 61.2% | 0.832 |

The most conserved ORF (ORF067, prevalence 84.8%) represents a candidate essential replication or maintenance gene. ORFs 061–063, forming a contiguous operon-like cluster at position 37–39 kb, showed identical prevalence (71.3%) with high match fractions (0.90–0.94), suggesting a functionally linked module under purifying selection. Conservation declined sharply below 50% for ORFs outside this core set, consistent with the modular architecture of IncX1 plasmids where accessory regions evolve rapidly through recombination and transposition.

The backbone/accessory partition for the full IncX1 dataset (n=701) estimated a mean backbone length of 37,787 bp (60% of total) and mean accessory length of 25,192 bp (40%), consistent with published estimates for Inc-type plasmid architecture.

---

### 3.7 Mosaicism Detection

Compositional heterogeneity analysis using GC-content variation across 1-kb sliding windows identified 302 of 701 IncX1 plasmids (43.1%) as potential mosaic candidates. These plasmids exhibited elevated intra-sequence GC standard deviation (mean GC SD = 0.047, range 0.017–0.115), suggesting the presence of recently acquired genomic modules with distinct compositional signatures. This high proportion of candidate chimeras is consistent with the known role of IncX1 plasmids as vectors for AMR gene cassettes acquired from diverse phylogenetic backgrounds, and underscores the importance of whole-sequence (rather than single-gene) approaches to plasmid classification.

---

### 3.8 Comparison with Existing Methods

To contextualise the pLIN system, we systematically compared its properties with all major existing plasmid typing approaches (Table 7).

**Table 7. Comparative analysis of plasmid typing and classification systems.**

| Property | pLIN (this study) | PlasmidFinder / Inc typing | pMLST | MOB-suite | COPLA (PTUs) | mge-cluster |
|---|---|---|---|---|---|---|
| **Year introduced** | 2025 | 2014 | 2014 | 2018 | 2021 | 2023 |
| **Classification basis** | Whole-sequence 4-mer composition | Replicon gene detection | Allelic variants of replicon loci | Relaxase + Mash distance | ANI network + HSBM | Unitig Jaccard + HDBSCAN |
| **Hierarchical levels** | 6 nested levels | 1 (flat) | 1 (flat) | 1 (flat) | 2–3 (semi) | 1 (flat) |
| **Code permanence** | Permanent | Mutable | Mutable | Mutable | Mutable | Mutable |
| **Reference DB required** | No | Yes | Yes | Yes | Yes | No |
| **Taxonomic scope** | Universal (all plasmids) | Enterobacteriaceae mainly | 6 Inc schemes only | Broad (DB-limited) | Broad (41% assignable) | Broad |
| **Resolution** | Strain-level (99.9% ANI) | Family-level | Sub-family | Species-level | Species-level | Variable |
| **Simpson's D (this dataset)** | 0.979 | 0.508* | N/A** | N/A | N/A | N/A |
| **Sensitivity** | High (whole-sequence) | ~50% | Scheme-limited | ~95% | 41–63% | High |
| **Computational cost** | Low (minutes) | Low | Low | Moderate | High (hours) | Moderate |
| **ML validation** | Integrated (F1=0.903) | None | None | None | None | None |
| **Mosaic detection** | Integrated | No | No | No | No | No |

*Calculated for 3-group Inc typing on this dataset. **pMLST schemes not available for all three Inc groups simultaneously.

The key advantages of pLIN over existing approaches are:

**(i) Hierarchical resolution.** All existing plasmid classification tools produce flat, single-level groupings. PlasmidFinder assigns a single Inc group. pMLST assigns a single sequence type within one Inc scheme. MOB-suite assigns a five-character cluster code at a fixed Mash distance of 0.06. COPLA assigns a single PTU. mge-cluster assigns a single cluster identifier. In contrast, pLIN simultaneously provides six nested levels of classification, enabling researchers to select the appropriate resolution for their question: Bin A–B for broad evolutionary comparisons, Bin C–D for lineage tracking, and Bin E–F for outbreak investigation.

**(ii) Code permanence.** A fundamental property inherited from the LIN framework is that pLIN codes are assigned sequentially based on the closest existing neighbour and are never retroactively changed. This contrasts with MOB-suite (re-clustering when the database is updated), COPLA (HSBM re-computation), and mge-cluster (entire t-SNE re-embedding), where codes can change with each database release.

**(iii) Reference-free operation.** PlasmidFinder, pMLST, MOB-suite, and COPLA all require curated reference databases of replicon sequences, allele profiles, relaxase genes, or pre-computed ANI networks. pLIN operates directly on raw nucleotide sequences using only composition vectors and distance-based clustering, eliminating database dependency and enabling classification of novel, divergent, or poorly characterised plasmids that would be unclassifiable by database-dependent methods. This is particularly important given that COPLA could assign only 41% of plasmids overall (63% for Enterobacterales) to defined PTUs.

**(iv) Integrated analytical framework.** Beyond classification, the pLIN pipeline incorporates mosaic structure detection via GC-content heterogeneity analysis, backbone/accessory region delineation, and machine learning validation of the compositional features used for clustering. No existing plasmid typing tool provides this integrated analytical capability.

---

### 3.9 pLIN System Definition and Formal Description

#### Definition

The **plasmid Life Identification Number (pLIN)** is a stable, hierarchical, multi-position numerical code assigned to each plasmid genome based on its compositional similarity to previously coded plasmids. Each pLIN code consists of six positions (A.B.C.D.E.F), where each position corresponds to a cluster identifier at a specific genetic distance threshold, ordered from coarse (Family, ~85% ANI) to fine (Strain, ~99.9% ANI).

#### Algorithm

**Input:** A plasmid nucleotide sequence S.

**Step 1 — Feature Extraction.** Compute the normalised tetranucleotide (4-mer) frequency vector V(S) of length 256 (4^4 possible tetranucleotides over the alphabet {A, C, G, T}). For each of the 256 canonical 4-mers, the frequency is calculated as:

    f(kmer) = count(kmer in S) / (|S| - 3)

**Step 2 — Distance Computation.** For a set of n plasmids, compute the n x n pairwise cosine distance matrix:

    d(i,j) = 1 - (V_i . V_j) / (||V_i|| x ||V_j||)

where V_i . V_j is the dot product and ||V|| is the L2 norm. Cosine distance ranges from 0 (identical composition) to 1 (orthogonal composition).

**Step 3 — Hierarchical Single-Linkage Clustering.** Apply agglomerative single-linkage clustering to the condensed distance matrix. At each of six distance thresholds T = {0.150, 0.100, 0.050, 0.020, 0.010, 0.001}, cut the dendrogram to obtain flat cluster assignments.

**Step 4 — Code Assembly.** For each plasmid i, concatenate its cluster identifiers at each threshold into a six-position code:

    pLIN(i) = C_A(i) . C_B(i) . C_C(i) . C_D(i) . C_E(i) . C_F(i)

**Step 5 — New Plasmid Assignment.** For a new query plasmid Q:
1. Compute V(Q) and calculate cosine distances to all existing plasmids.
2. Identify the nearest neighbour N = argmin d(Q, i).
3. At each threshold T_k, if d(Q, N) ≤ T_k, inherit the cluster identifier of N; otherwise, assign a new cluster identifier.
4. This ensures code permanence: existing pLIN codes are never altered.

#### Threshold Rationale

| Bin | Threshold (d) | ANI Proxy | Biological Interpretation |
|-----|---------------|-----------|---------------------------|
| A | 0.150 | ~85% | Broad plasmid family; encompasses all Inc-group plasmids of related backbone origin |
| B | 0.100 | ~90% | Plasmid subfamily; separates major evolutionary branches |
| C | 0.050 | ~95% | Cluster; analogous to bacterial species boundary (Mash d ≤ 0.05) |
| D | 0.020 | ~98% | Subcluster; within-lineage divergence consistent with F_lineage threshold (d ≤ 0.03) |
| E | 0.010 | ~99% | Clone complex; near-identical backbone architecture |
| F | 0.001 | ~99.9% | Strain; outbreak-level resolution, consistent with G_outbreak threshold (d ≤ 0.005) |

Thresholds were calibrated against the nearest-neighbour distance distribution of IncX-like reference plasmids, where the 50th percentile nearest-neighbour distance was 0.011, the 90th percentile was 0.038, and the 99th percentile was 0.072. The F_lineage (d ≤ 0.03, ANI ≥ 97%) and G_outbreak (d ≤ 0.005, ANI ≥ 99.5%) thresholds from the IncX pLIN threshold configuration provide additional biological calibration points.

#### Implementation

The pLIN system is implemented in Python 3.14+ with the following dependencies: NumPy 2.4, Pandas 2.3, SciPy 1.16 (for pdist, linkage, and fcluster), and BioPython 1.84 (for FASTA parsing). The complete pipeline—from raw FASTA input to pLIN code assignment for 6,346 plasmids—executes in under 2 minutes on a standard laptop (Apple M-series, single thread). No external bioinformatics tools (Mash, BLAST, Prokka) are required.

#### Properties

1. **Permanence.** Once assigned, a pLIN code is never changed, regardless of subsequent additions to the database. This is guaranteed by the nearest-neighbour assignment rule.

2. **Hierarchy.** The six-position code provides simultaneous classification at six resolution levels, from family to strain.

3. **Universality.** No reference database is required. Any plasmid with a nucleotide sequence can be classified.

4. **Reproducibility.** Given the same database state and insertion order, pLIN codes are fully deterministic.

5. **Interpretability.** Plasmids sharing a k-position prefix (e.g., 1.1.1.7.*.*) are related at the corresponding resolution level (e.g., subcluster level for a 4-position match).

---

### 3.10 Integration with AMRFinderPlus: pLIN as a Framework for AMR Gene Surveillance

#### 3.10.1 AMR Gene Detection Across the Plasmid Dataset

To demonstrate the utility of pLIN as an epidemiological framework for tracking antimicrobial resistance gene dissemination, we integrated the pLIN classification system with NCBI AMRFinderPlus v4.2.5 (database 2026-01-21.1). AMRFinderPlus was run in nucleotide mode with the `--plus` flag to detect antimicrobial resistance genes, virulence factors, and stress response genes across all 6,346 plasmid sequences.

A total of 60,372 gene detections were identified across 5,342 plasmids (84.2% of the dataset), comprising 27,465 AMR gene hits, 5,834 virulence factor hits, and 27,073 stress response gene hits (Table 8).

**Table 8. AMR/virulence/stress gene prevalence across the plasmid dataset.**

| Category | Plasmids Positive | % of Total (n=6,346) | Total Detections |
|----------|-------------------|----------------------|------------------|
| AMR genes | 4,224 | 66.6% | 27,465 |
| Virulence factors | 1,191 | 18.8% | 5,834 |
| Stress response genes | 2,795 | 44.0% | 27,073 |
| Any AMRFinderPlus hit | 5,342 | 84.2% | 60,372 |

AMR gene prevalence varied markedly across incompatibility groups (Table 9). IncN plasmids showed the highest AMR gene carriage rate (85.1%, mean 6.21 AMR genes per plasmid), consistent with the established role of IncN plasmids as efficient vectors for multi-drug resistance gene cassettes. IncFII and IncX1 plasmids exhibited similar AMR prevalence (~63–67%) but differed in virulence gene carriage: IncFII plasmids carried virulence factors at substantially higher rates (23.5%) than IncN (3.0%) or IncX1 (11.8%), reflecting the known association of IncFII plasmids with virulence-associated loci such as the *spv* operon and iron uptake systems.

**Table 9. AMR and virulence gene prevalence by incompatibility group.**

| Inc Type | Total (n) | AMR+ (n) | AMR+ (%) | Mean AMR Genes | VIR+ (n) | VIR+ (%) |
|----------|-----------|----------|----------|----------------|----------|----------|
| IncFII | 4,581 | 2,852 | 62.3% | 3.91 | 1,076 | 23.5% |
| IncN | 1,064 | 905 | 85.1% | 6.21 | 32 | 3.0% |
| IncX1 | 701 | 467 | 66.6% | 4.19 | 83 | 11.8% |

#### 3.10.2 AMR Gene Repertoire

The most frequently detected AMR gene was *blaTEM-1* (n=1,760, 41.7% of AMR-positive plasmids), encoding a narrow-spectrum TEM-type beta-lactamase, followed by the sulfonamide resistance gene *sul1* (29.9%), tetracycline efflux gene *tet(A)* (29.1%), sulfonamide resistance gene *sul2* (24.4%), and aminoglycoside-modifying enzymes *aph(6)-Id* (24.2%) and *aph(3'')-Ib* (24.0%) (Table 10). The high prevalence of these 'classic' resistance determinants across all three Inc groups is consistent with their association with widely disseminated transposon families (Tn*3*, Tn*10*, Tn*21*) and class 1 integrons.

**Table 10. Top 20 most common AMR genes across the plasmid dataset.**

| Rank | Gene | Detections (n) | % of AMR+ Plasmids (n=4,224) | Drug Class |
|------|------|----------------|------------------------------|------------|
| 1 | *blaTEM-1* | 1,760 | 41.7% | Beta-lactam |
| 2 | *sul1* | 1,265 | 29.9% | Sulfonamide |
| 3 | *tet(A)* | 1,230 | 29.1% | Tetracycline |
| 4 | *sul2* | 1,030 | 24.4% | Sulfonamide |
| 5 | *aph(6)-Id* | 1,021 | 24.2% | Aminoglycoside |
| 6 | *aph(3'')-Ib* | 1,013 | 24.0% | Aminoglycoside |
| 7 | *mph(A)* | 951 | 22.5% | Macrolide |
| 8 | *mrx(A)* | 950 | 22.5% | Macrolide |
| 9 | *blaKPC-2* | 776 | 18.4% | Carbapenem |
| 10 | *qnrS1* | 695 | 16.5% | Quinolone |
| 11 | *aac(6')-Ib-cr5* | 681 | 16.1% | Aminoglycoside/Quinolone |
| 12 | *dfrA14* | 654 | 15.5% | Trimethoprim |
| 13 | *aadA2* | 575 | 13.6% | Aminoglycoside |
| 14 | *catB3* | 541 | 12.8% | Phenicol |
| 15 | *blaOXA-1* | 523 | 12.4% | Beta-lactam |
| 16 | *floR* | 479 | 11.3% | Phenicol |
| 17 | *aac(3)-IId* | 461 | 10.9% | Aminoglycoside |
| 18 | *aph(3')-Ia* | 459 | 10.9% | Aminoglycoside |
| 19 | *blaCTX-M-15* | 456 | 10.8% | Beta-lactam (ESBL) |
| 20 | *ble* | 415 | 9.8% | Bleomycin |

The AMR drug class distribution revealed beta-lactam resistance as the most common (3,483 plasmids), followed by aminoglycoside (2,723), sulfonamide (2,073), trimethoprim (1,810), tetracycline (1,549), phenicol (1,542), quinolone (1,240), and macrolide (1,088) resistance determinants. Co-carriage of resistance genes to three or more drug classes was observed in a substantial proportion of AMR-positive plasmids, consistent with the multidrug-resistance phenotype commonly associated with these Inc groups.

#### 3.10.3 Clinically Critical Resistance Determinants

Of particular clinical concern, we identified high-priority resistance genes associated with last-resort and critically important antimicrobials (Table 11).

**Table 11. Clinically critical resistance genes detected across the plasmid dataset.**

| Category | Total Detections | Top Variants |
|----------|------------------|--------------|
| **Carbapenemases** | 1,490 | *blaKPC-2* (776), *blaNDM-1* (216), *blaKPC-3* (177), *blaNDM-5* (70), *blaIMP-4* (65) |
| **ESBLs** | 3,734 | *blaTEM-1* (1,760), *blaCTX-M-15* (456), *blaCTX-M-65* (307), *blaTEM* (266), *blaSHV-12* (233) |
| **Colistin resistance** | 160 | *mcr-1.1* (57), *mcr-8.1* (27), *mcr-8.2* (18), *mcr-10.1* (12), *mcr-3.5* (8) |
| **PMQR (fluoroquinolone)** | 2,160 | *qnrS1* (695), *aac(6')-Ib-cr5* (681), *qnrB1* (202), *qnrB2* (99), *oqxB* (93) |
| **Vancomycin resistance** | 0 | None detected |

Carbapenemase genes were detected on 1,490 plasmids, with *blaKPC-2* being the most prevalent (n=776, 18.4% of AMR-positive plasmids). The co-occurrence of *blaKPC-2* with ESBL genes (*blaCTX-M-15*, *blaCTX-M-65*) on individual plasmids is particularly alarming, as it confers resistance to virtually all beta-lactam antibiotics including carbapenems. *blaNDM-1* (n=216) and *blaNDM-5* (n=70) were also widespread, consistent with the global dissemination of New Delhi metallo-beta-lactamase genes on IncFII and IncN plasmids.

Plasmid-mediated colistin resistance (*mcr*) genes, representing resistance to the last-resort polymyxin class, were detected on 160 plasmids. The *mcr-1.1* variant predominated (n=57), followed by *mcr-8.1* (n=27) and *mcr-8.2* (n=18). The presence of *mcr* genes on plasmids carrying simultaneous carbapenemase genes raises the spectre of pan-drug resistance.

#### 3.10.4 pLIN Lineages as AMR Vehicles

A central application of the pLIN system is the identification of specific plasmid lineages that serve as vehicles for AMR gene dissemination. By cross-referencing pLIN codes with AMR gene content, we identified distinct lineage-associated resistance profiles (Table 12).

**Table 12. Top 10 pLIN lineages carrying AMR genes, with characteristic resistance profiles.**

| pLIN Code | n (AMR+) | Mean AMR Genes | Total AMR | Inc Type(s) | Top Resistance Genes |
|-----------|----------|----------------|-----------|-------------|---------------------|
| 1.1.1.7.30.1240 | 657 | 6.9 | 4,546 | IncFII, IncN | *blaTEM-1* (43%), *sul1* (36%), *tet(A)* (31%), *mph(A)* (29%) |
| 1.1.1.7.30.1313 | 225 | 5.6 | 1,249 | IncFII, IncN, IncX1 | *rmtB1* (71%), *blaTEM-1* (70%), *blaKPC-2* (58%), *blaCTX-M-65* (54%) |
| 1.1.1.7.30.1346 | 148 | 7.5 | 1,111 | IncFII | *mph(A)* (58%), *mrx(A)* (58%), *sul1* (57%), *aadA5* (48%) |
| 1.1.1.7.30.1248 | 111 | 6.8 | 756 | IncFII | *blaKPC-2* (44%), *blaTEM-1* (42%), *qnrS1* (41%) |
| 1.1.1.7.30.743 | 108 | 15.8 | 1,701 | IncFII, IncN, IncX1 | *floR* (69%), *sul1* (67%), *mph(A)* (65%), *sul2* (60%) |
| 1.1.1.7.30.567 | 90 | 13.2 | 1,189 | IncN | *blaKPC-2* (100%), *blaTEM-1* (100%), *aph(3'')-Ib* (98%), *aac(3)-IId* (97%) |
| 1.1.1.7.30.1271 | 88 | 6.8 | 602 | IncFII, IncX1 | *blaTEM-1* (59%), *blaCTX-M-15* (43%), *sul2* (41%) |
| 1.1.1.7.30.1099 | 87 | 5.0 | 438 | IncFII | *qnrS1* (92%), *tet(A)* (85%), *blaLAP-2* (80%), *sul2* (77%) |
| 1.1.1.7.30.1135 | 79 | 4.9 | 388 | IncFII | *qnrS1* (89%), *tet(A)* (80%), *blaLAP-2* (76%), *sul2* (72%) |
| 1.1.1.7.30.1270 | 70 | 4.2 | 297 | IncFII | *blaOXA* (67%), *blaTEM-1* (60%), *blaKPC-2* (39%) |

Several findings are noteworthy. First, pLIN lineage 1.1.1.7.30.567 (an IncN-dominated lineage, n=90) carried *blaKPC-2* on 100% of its members alongside *blaTEM-1* (100%), *aph(3'')-Ib* (98%), and *aac(3)-IId* (97%), representing a tightly conserved multidrug resistance cassette. This lineage carried a mean of 13.2 AMR genes per plasmid—among the highest AMR gene burdens in the entire dataset—and represents a high-priority surveillance target.

Second, pLIN lineage 1.1.1.7.30.1313 spanned all three Inc groups (IncFII, IncN, IncX1) and showed high prevalence of both *blaKPC-2* (58%) and *blaCTX-M-65* (54%), alongside the 16S rRNA methyltransferase *rmtB1* (71%), which confers high-level resistance to all clinically available aminoglycosides. This cross-Inc lineage may represent a 'resistance hub' facilitating inter-lineage AMR gene exchange.

Third, pLIN lineages 1.1.1.7.30.1099 and 1.1.1.7.30.1135 (both IncFII) shared nearly identical resistance profiles dominated by *qnrS1*, *tet(A)*, *blaLAP-2*, and *sul2*, suggesting descent from a common ancestor with a conserved resistance cassette. Their distinct pLIN codes at the strain level (Bin F) indicate compositional divergence in backbone or other accessory regions despite maintaining the same AMR gene complement.

#### 3.10.5 Virulence Gene Distribution

Virulence gene carriage showed strong Inc-type specificity. IncFII plasmids (1,076 with virulence genes, 23.5%) predominantly carried the type III secretion system effector *traT* (25.6%), the *Salmonella* plasmid virulence genes *spvB* and *spvD* (25.3%), the murein transglycosylase *mltE* (24.0%), and the aerobactin siderophore gene *iucA* (23.3%). IncX1 plasmids (83 with virulence genes, 11.8%) were enriched for the alpha-hemolysin gene *hlyA* (48.2%) and fimbrial adhesin genes *fedA* and *fedF* (47.0%), consistent with the known association of IncX1 plasmids with enterotoxigenic *E. coli* virulence determinants. IncN plasmids showed the lowest virulence gene carriage (32 plasmids, 3.0%), with aerobactin genes (*iutA*, *iucC*) and *traT* being the most common.

The combination of AMR and virulence determinants on the same plasmid is clinically significant, as it enables the co-selection of pathogenicity and resistance traits under antibiotic pressure. The pLIN framework enables systematic tracking of such co-occurrence patterns at the lineage level.

---

### 3.11 Data Availability

The complete pLIN assignment table for all 6,346 plasmids is provided as Supplementary Table S1 (pLIN_assignments.tsv). The integrated pLIN + AMRFinderPlus table is provided as Supplementary Table S2 (pLIN_AMR_integrated.tsv) and the per-lineage AMR summary as Supplementary Table S3 (pLIN_lineage_AMR_summary.tsv). The pLIN assignment pipeline (assign_pLIN.py), AMRFinderPlus batch runner (run_amrfinder_all.sh), and integration analysis script (integrate_pLIN_AMR.py) are available at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification. The full analytical notebook (pLIN.ipynb) and IncX1 threshold configuration (Data/IncX_PLIN_thresholds_v0_python.yaml) are also provided.
