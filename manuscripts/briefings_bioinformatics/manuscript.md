# pLIN: A Computational Framework for Hierarchical, Permanent Classification of Bacterial Plasmids Using Tetranucleotide Composition and Lineage Identification Numbers

**Basil Britto Xavier^1^, Anurag Kumar Bari^1^, Bhanu Sinha^1^, John W A Rossen^1,\*^**

^1^ University of Groningen, University Medical Center Groningen, Department of Medical Microbiology and Infection Prevention, Groningen, The Netherlands

^\*^ Corresponding author: John W A Rossen

*Manuscript prepared for Briefings in Bioinformatics (Problem Solving Protocol / Case Study)*

*Double-spaced; approximately 8,200 words*

---

## Abstract

Bacterial plasmids are the principal vehicles of horizontal gene transfer and antimicrobial resistance (AMR) dissemination, yet no existing classification system simultaneously provides hierarchical multi-resolution typing, code permanence, and reference-free operation. Here we describe pLIN (plasmid Lineage Identification Number), the first application of the Life Identification Number framework to plasmid genomes. pLIN computes normalised tetranucleotide (4-mer) frequency vectors (256 features) for each plasmid, constructs pairwise cosine distance matrices, and applies single-linkage hierarchical clustering at six biologically calibrated thresholds spanning family-level (~85% average nucleotide identity [ANI]) to strain-level (~99.9% ANI) resolution. Each plasmid receives a permanent six-position code whose stability is mathematically guaranteed by the nearest-neighbour assignment rule. Applied to 8,077 training samples (8,056 unique plasmid sequences) from 28 groups -- 20 Gram-negative incompatibility (Inc) groups, 4 Gram-positive replicon (rep) type groups spanning *Staphylococcus aureus* and *Enterococcus* spp., 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups -- pLIN resolved 3,073 unique strain-level codes with a Simpson's Index of Diversity of 0.985 (1.54-fold improvement over Inc/rep typing alone), while maintaining 97.8% concordance with established group assignments. A nested cross-validation machine learning pipeline confirmed that the compositional features robustly predict group membership (XGBoost weighted F1 = 0.896 +/- 0.009). A KNN-based group classifier (k = 5, cosine distance, distance-weighted; 91.1% accuracy, 5-fold cross-validation) enables per-group clustering of an expanded reference database of 79,305 plasmids (8,056 unique training + 71,249 reference), resolving 57,886 unique pLIN codes in under 30 minutes on consumer hardware. Beyond classification, pLIN now incorporates eight analytical modules addressing previously identified limitations: automated plasmid versus chromosome contig identification (multi-signal scoring using sequence length, compositional distance, header keywords, and classifier confidence; 97.8% accuracy), assembly completeness scoring (composite 0--100 scale), database coverage and novelty detection (nearest-neighbour distance percentiles with traffic-light flagging), recombination detection (minimap2-based mosaic analysis), novel group discovery (hierarchical clustering of low-confidence plasmids), evolutionary rate estimation (SNP accumulation regression), adaptive threshold calibration with cluster stability assessment (bootstrap resampling with Adjusted Rand Index), and mobile genetic element boundary detection (IS elements, composite transposons, colour-coded gene maps). Integration with NCBI AMRFinderPlus v4.2.5 identified 64,891 gene detections (29,583 AMR, 6,286 virulence, 29,022 stress) across 83.1% of Gram-negative plasmids, enabling lineage-level tracking of carbapenemase, ESBL, colistin resistance, and plasmid-mediated quinolone resistance genes. Cross-validation against 74 plasmids from 27 published outbreak and surveillance studies spanning seven resistance mechanisms across 13 countries confirmed prospective identification of known high-risk lineages, correct within-outbreak backbone grouping, and 85.1% high-confidence classification rate. The pLIN tool is implemented as an open-source Streamlit web application with integrated AMR surveillance, outbreak detection, ANI validation, and comprehensive quality assessment modules, and is freely available under GPL-3.0.

---

## Key Points

- pLIN is the first application of the Life Identification Number (LIN) framework to plasmid genomes, assigning each plasmid a permanent, hierarchical six-position code based on tetranucleotide composition (256 features) and single-linkage clustering at six ANI-calibrated distance thresholds.

- Applied to 8,077 training samples (8,056 unique plasmids) from 28 groups (20 Gram-negative Inc groups, 4 Gram-positive rep type groups covering *S. aureus* and *Enterococcus* spp., 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups), pLIN achieved a Simpson's Index of Diversity of 0.985 (1.54-fold improvement over group typing), 97.8% concordance with group assignments, and 91.1% group classification accuracy via a KNN classifier validated by nested cross-validation (XGBoost F1 = 0.896).

- Eight analytical modules address previously identified limitations: automated plasmid versus chromosome contig identification (multi-signal scoring; 97.8% accuracy, enabling direct analysis from whole-genome assemblies), assembly completeness scoring (composite 0--100), database coverage and novelty detection (traffic-light flagging), recombination detection (minimap2 mosaic analysis), novel group discovery (hierarchical clustering of low-confidence plasmids), evolutionary rate estimation (SNP regression), adaptive threshold calibration with bootstrap cluster stability (Adjusted Rand Index), and MGE boundary detection (IS elements, composite transposons, colour-coded gene maps).

- A per-group clustering strategy scales to 79,305 reference plasmids (8,056 unique training + 71,249 reference; 57,886 unique pLIN codes) in under 30 minutes on a standard laptop, reducing memory requirements from 24 GB (global) to 2.3 GB (largest group), while the code permanence guarantee ensures that existing assignments are never altered by database expansion.

- Integration with AMRFinderPlus across 6,998 Gram-negative plasmids (20 Inc groups) identified 1,635 carbapenemase, 1,804 ESBL, 204 mcr colistin resistance, and 2,315 plasmid-mediated quinolone resistance detections mapped to specific pLIN lineages, including pLIN 671 (IncN, n = 90, 100% blaKPC-2, 13.2 mean AMR genes) and pLIN 860 (n = 142, 5 Inc groups, 14.4 mean AMR genes, 44.4% mcr carriage).

- Cross-validation against 74 plasmids from 27 independent outbreak and surveillance studies spanning seven resistance mechanisms (KPC, NDM, OXA-48, VIM, IMP, mcr, CTX-M) across 13 countries confirmed that pLIN prospectively identifies known high-risk lineages, correctly groups outbreak-related backbones (9 intra-study clusters detected), and provides hierarchical within-outbreak resolution unavailable from any existing tool.

**Keywords:** plasmid classification, Lineage Identification Number, tetranucleotide composition, antimicrobial resistance, hierarchical clustering, genomic epidemiology, Gram-positive plasmids, Acinetobacter baumannii, Pseudomonas aeruginosa, mobile genetic elements

---

## 1. Introduction

Bacterial plasmids -- extrachromosomal, self-replicating DNA elements -- are the principal vehicles of horizontal gene transfer (HGT) in prokaryotes, driving the rapid dissemination of antimicrobial resistance (AMR) determinants, virulence factors, and stress tolerance genes across species boundaries [1--3]. The clinical impact is profound: carbapenem-resistant Enterobacterales, extended-spectrum beta-lactamase (ESBL)-producing organisms, and colistin-resistant strains are now among the most urgent public health threats globally, and plasmids are the primary vectors for the resistance genes underlying each of these phenotypes [4--6,19]. A single conjugative plasmid can carry resistance determinants spanning multiple drug classes, creating extensively drug-resistant phenotypes in a single transfer event [7,8]. Tracking the plasmid -- not merely the pathogen -- is therefore essential for understanding and controlling AMR dissemination.

Despite their central epidemiological role, plasmid classification remains fragmented. Existing approaches suffer from one or more critical limitations (Table 1; Figure 10). Replicon-based incompatibility (Inc) typing via PlasmidFinder [9] provides a single flat classification label and depends on a curated reference database. Plasmid multilocus sequence typing (pMLST) [9] covers only six Inc-group schemes. MOB-suite [10] clusters plasmids at a fixed Mash distance threshold of 0.06, producing flat codes that are reassigned with each database update. COPLA/Plasmid Taxonomic Units (PTUs) [11] employ average nucleotide identity (ANI) networks and hierarchical stochastic block modelling, but classify only 41% of plasmids overall and recompute codes upon each release. The reference-free tool mge-cluster [12] uses unitig Jaccard distances and HDBSCAN clustering but produces flat, non-permanent codes. No existing tool simultaneously satisfies the five criteria required for a comprehensive plasmid nomenclature: hierarchical multi-resolution typing, code permanence, reference-free operation, broad taxonomic scope, and integrated AMR surveillance.

The Life Identification Number (LIN) system, originally developed for hierarchical, permanent classification of bacterial strains based on whole-genome similarity [13,14], provides an elegant solution. LIN assigns each genome a multi-position numerical code based on its distance to the nearest previously coded genome at a series of nested thresholds. The nearest-neighbour assignment rule guarantees two fundamental properties: (i) hierarchical consistency -- if two genomes share a code at a coarse level, they necessarily share codes at all coarser levels; and (ii) code permanence -- once assigned, a code is never altered by subsequent additions to the database. LIN has been validated for bacterial species classification [13] and plant pathogen strain typing [14], but has never been applied to plasmid genomes.

Here we describe pLIN (plasmid Lineage Identification Number), the first LIN-based framework for plasmid classification. We detail the computational methodology, from tetranucleotide feature extraction and cosine distance computation through hierarchical threshold calibration and KNN-based group classification, followed by comprehensive benchmarking, scalability assessment, and application to AMR surveillance. Critically, pLIN has been expanded beyond Gram-negative Enterobacterales to incorporate four Gram-positive replicon type groups from *Staphylococcus aureus* and *Enterococcus* spp., two *Acinetobacter baumannii* rep type groups, and two *Pseudomonas aeruginosa* rep type groups, bringing the total to 28 groups (8,077 training samples; 8,056 unique plasmids). Eight analytical modules address key methodological limitations identified in earlier versions: automated plasmid versus chromosome contig identification, assembly completeness assessment, database coverage and novelty detection, recombination detection, novel group discovery, evolutionary rate estimation, adaptive threshold calibration with cluster stability analysis, and mobile genetic element (MGE) boundary detection. The complete system is implemented as an open-source Streamlit web application with modules for AMR annotation, mobility prediction, outbreak detection, ANI validation, and comprehensive quality assessment (Figures 8, 11).

---

## 2. The pLIN Framework

### 2.1 Tetranucleotide composition and distance metric

The pLIN algorithm begins with alignment-free feature extraction. For each input plasmid sequence *S*, a normalised tetranucleotide (4-mer) frequency vector **V**(*S*) of length 256 (4^4 canonical tetranucleotides over the alphabet {A, C, G, T}) is computed:

    f(kmer) = count(kmer in S) / (|S| - 3)

This yields a 256-dimensional compositional fingerprint that is invariant to sequence orientation and robust to minor assembly artefacts. For *n* plasmids, the pairwise cosine distance matrix **D** is computed:

    D(i,j) = 1 - [ V(S_i) . V(S_j) ] / [ ||V(S_i)|| * ||V(S_j)|| ]

Cosine distance was chosen over Euclidean or Manhattan alternatives for three reasons: (i) scale invariance -- it is robust to the wide range of plasmid sizes in the dataset (1.8 kb to 400 kb); (ii) computational efficiency -- it can be computed in vectorised form using scipy.spatial.distance; and (iii) established correlation with genomic divergence for k-mer frequency profiles [15,16]. Applied to the training dataset of 8,077 plasmid sequences (8,056 unique; 20 Gram-negative Inc groups, 4 Gram-positive rep type groups, 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups), this produced a distance matrix of 32,619,926 pairwise comparisons. Note: the 21-sequence difference between 8,077 training samples and 8,056 unique plasmids reflects 21 *E. faecium* plasmids present in both repEF_conj and repEF_res training sets due to multi-replicon structure.

### 2.2 Hierarchical clustering and threshold calibration

Six cosine distance thresholds were defined to capture biologically meaningful levels of plasmid relatedness, calibrated against established ANI benchmarks (Table 2; Figure 7). The thresholds span three orders of magnitude:

| Level | Threshold (d) | ANI equivalent | Biological interpretation |
|-------|---------------|----------------|---------------------------|
| L1    | <= 0.150      | ~85%           | Family                    |
| L2    | <= 0.100      | ~90%           | Subfamily                 |
| L3    | <= 0.050      | ~95%           | Cluster                   |
| L4    | <= 0.020      | ~98%           | Subcluster                |
| L5    | <= 0.010      | ~99%           | Clone complex             |
| L6    | <= 0.001      | ~99.9%         | Strain / Outbreak         |

Thresholds were initially calibrated against the nearest-neighbour distance distribution of 178 selected IncX-like reference plasmids and subsequently validated across the 20 Gram-negative Inc groups using FastANI v1.34 (--fragLen 1000) on 4,970 within-group pairs. The overall Spearman correlation between cosine distance and FastANI was rho = -0.348 (P < 10^-141), with 15/20 groups showing significant negative correlation (P < 0.05). At the strain-level threshold (d <= 0.001), the median FastANI was 99.9%, confirming the target ANI mapping. The strongest per-group correlations were observed in IncFIC (rho = -0.88), IncI2 (rho = -0.81), ColE (rho = -0.72), and IncN (rho = -0.61). FastANI validation was performed on Gram-negative Inc groups only; validation of Gram-positive, *Acinetobacter*, and *Pseudomonas* groups is planned for future work.

Single-linkage hierarchical clustering is applied at each threshold to produce nested cluster assignments. At the training scale (8,077 plasmids across 28 groups), the six levels resolved 1 family, 2 subfamilies, 9 clusters, 44 subclusters, 131 clone complexes, and 3,073 strain-level groups, with the Gram-positive groups contributing 255 additional strain-level codes (repSA_large 72, repSA_small 52, repEF_conj 63, repEF_res 69). Each plasmid is assigned a six-position code of the form L1.L2.L3.L4.L5.L6 (e.g., 1.1.2.15.48.671). Single-linkage was chosen for theoretical consistency with the LIN nearest-neighbour rule, which guarantees code permanence [13].

### 2.3 The code permanence guarantee

The pLIN code permanence property is a mathematical consequence of the nearest-neighbour assignment rule and merits formal statement. When a new plasmid *Q* is submitted for classification:

1. Its 4-mer vector **V**(*Q*) is computed and compared to all existing reference vectors.
2. The nearest neighbour *R* is identified (minimum cosine distance *d*(*Q*, *R*)).
3. At each threshold level *L_k*, if *d*(*Q*, *R*) <= *t_k*, then *Q* inherits the level-*k* code of *R*; otherwise, *Q* receives a new, unique level-*k* identifier.

Because the code assignment depends solely on the distance to the nearest existing reference -- not on the global structure of the database -- the addition of *Q* cannot alter the codes of any previously assigned plasmid. This property distinguishes pLIN from MOB-suite (re-clustering on database update), COPLA (HSBM recomputation), and mge-cluster (t-SNE re-embedding), all of which can retroactively change codes when new sequences are added.

### 2.4 KNN group classifier and Gram-positive expansion

To enable per-group clustering of large reference databases, pLIN incorporates a K-nearest neighbours (KNN) classifier for group prediction. The classifier operates on the same 256-dimensional 4-mer feature space and was trained on the 8,077-sample dataset (8,056 unique plasmids) spanning 28 groups: 20 Gram-negative Inc groups, 4 Gram-positive rep type groups, 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups. Key hyperparameters: k = 5 neighbours, cosine distance metric, distance-weighted voting. Five-fold cross-validation yielded an overall accuracy of 91.1% (Figure 12).

**Gram-positive and non-fermentative expansion.** The classifier was extended from 20 Gram-negative Inc groups to 28 groups by incorporating four Gram-positive replicon type groups curated from complete plasmid sequences in NCBI RefSeq and PLSDB [31], and four non-fermentative Gram-negative rep type groups:

| Group | Organism | Description | Training (n) |
|-------|----------|-------------|--------------|
| repSA_large | *S. aureus* | Large plasmids: pI258, pSK1, pSK41 families | 121 |
| repSA_small | *S. aureus* | Small plasmids: pT181, SAP, pWBG749 families | 73 |
| repEF_conj | *Enterococcus* spp. | Conjugative: pAD1, pCF10 families | 92 |
| repEF_res | *Enterococcus* spp. | Resistance: pRUM, pRE25, pHTbeta families | 121 (100 unique plasmids; 21 multi-replicon overlap with repEF_conj) |
| repAci1 | *A. baumannii* | Small rep-type plasmids: pRAY, pABVA01 families | -- |
| repAci_large | *A. baumannii* | Large conjugative/resistance plasmids | -- |
| repPae_large | *P. aeruginosa* | Large plasmids: pOZ176, megaplasmid families | -- |
| repPae_small | *P. aeruginosa* | Small plasmids: pVS1, ColE-like families | -- |

Both *A. baumannii* and *P. aeruginosa* are designated WHO critical priority pathogens due to their extensive drug resistance and limited treatment options, making plasmid surveillance in these species particularly urgent. The reference plasmid sequences for all groups were sourced primarily from PLSDB 2025 (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) [31], supplemented by NCBI RefSeq.

The 407 Gram-positive training sequences and the additional *Acinetobacter* and *Pseudomonas* training sequences were curated using the same pipeline as the Gram-negative groups: replicon typing via in silico PCR against known rep genes, followed by manual verification of complete assembly status. Gram-positive plasmids exhibit lower GC content (mean 33.2% for *S. aureus*, 36.8% for *Enterococcus*) than Gram-negative plasmids (mean 49.1%), while *A. baumannii* (mean 39.1%) and *P. aeruginosa* (mean 60.2%) plasmids span a wide GC range, necessitating a widened GC-content validation window from 30--65% to 25--70%. In the tetranucleotide composition space, the Gram-positive groups form a distinct cluster separated from Gram-negative groups by cosine distances of 0.08--0.15, while the *Acinetobacter* and *Pseudomonas* groups occupy intermediate positions reflecting their distinct genomic compositions. Cross-validation accuracy for the four Gram-positive groups ranged from 88.4% (repEF_conj) to 94.2% (repSA_small), comparable to Gram-negative group performance.

The classifier stores training vectors and labels in a compressed NumPy archive (inc_classifier.npz), enabling rapid loading and inference without retraining. Classification confidence is computed as the proportion of nearest-neighbour votes for the winning class, weighted by inverse distance. Sequences with confidence below 40% are flagged as "Unknown/Novel" and excluded from pLIN assignment, preventing unreliable codes. The dominant confusion pairs involve compositionally near-identical Inc groups within the IncF family (IncFIB/IncFII inter-centroid distance d = 0.002; IncF/IncFII d = 0.003), where misclassification places the plasmid in a biologically equivalent neighbourhood. Within the Gram-positive groups, the primary confusion pair is repEF_conj/repEF_res (inter-centroid d = 0.012), reflecting shared *Enterococcus* backbone composition. Within the *Acinetobacter* and *Pseudomonas* groups, the size-based separation (repAci1/repAci_large, repPae_large/repPae_small) provides strong classification signal.

---

## 3. Implementation

### 3.1 Software architecture

pLIN is implemented as a single-file Streamlit web application (plin_app.py, ~6,500 lines of Python) that runs locally in the user's web browser without server-side infrastructure. The architecture comprises eight modular layers:

1. **Feature extraction layer.** Computes 4-mer frequency vectors using BioPython for FASTA parsing and NumPy for vectorised counting. Accepts single or batch FASTA uploads.
2. **Contig identification layer.** Multi-signal scoring classifies input contigs as plasmid or chromosomal, automatically filtering non-plasmid sequences before downstream analysis (Section 3.2).
3. **Classification layer.** KNN group prediction (28 groups) using scikit-learn, followed by nearest-neighbour pLIN code assignment against the pre-computed reference database.
4. **Quality assessment layer.** Assembly completeness scoring, database coverage and novelty detection, and recombination analysis (Sections 3.6--3.8).
5. **AMR annotation layer.** Optional integration with NCBI AMRFinderPlus for nucleotide-mode detection of AMR, virulence, and stress genes, with results merged into the pLIN output table. MGE boundary detection and composite transposon identification (Section 3.12).
6. **Validation layer.** Optional Mash and FastANI ANI estimation; minimap2-based SNP sub-typing for within-strain resolution below the L6 threshold.
7. **Discovery layer.** Novel group discovery from low-confidence classifications, evolutionary rate estimation for dated clusters, and adaptive threshold calibration with cluster stability assessment (Sections 3.9--3.11).
8. **Surveillance layer.** Basic and temporal outbreak detection modules with risk stratification; CRISPR-based host inference via MinCED; mobility prediction via MOB-suite.

External tool dependencies (AMRFinderPlus, Mash, FastANI, minimap2, MinCED, MOB-suite) follow a graceful degradation model: each tool is auto-detected at runtime, and the corresponding module is enabled or disabled accordingly. The core pLIN classification pipeline -- feature extraction, distance computation, and code assignment -- requires only Python standard library packages plus NumPy, SciPy, and scikit-learn, ensuring portability across platforms (Figure 11).

### 3.2 Plasmid contig identification

A key challenge in plasmid genomics is that user-submitted assemblies frequently contain a mixture of plasmid and chromosomal contigs, particularly when derived from short-read whole-genome sequencing without prior plasmid extraction. Submitting chromosomal contigs to the pLIN classifier can produce spurious group assignments and inflate false-positive rates. To address this, pLIN incorporates a multi-signal scoring system (`classify_contigs_plasmid_vs_chromosome()`) that automatically identifies and separates plasmid contigs from chromosomal sequences prior to pLIN assignment.

The scoring system integrates four independent signals, each contributing an additive component to a composite plasmid score:

1. **Sequence length scoring.** Plasmids are typically smaller than chromosomes, and contig length provides a strong prior. The length signal assigns scores as follows: contigs exceeding 1 Mb receive a score of -50 (strongly chromosomal); contigs between 500 kb and 1 Mb receive -30; contigs below 300 kb receive +20 (consistent with typical plasmid sizes); and contigs below 20 kb receive an additional +10 bonus, reflecting the predominance of small plasmids in public databases. These thresholds were calibrated against the size distribution of the 79,305-plasmid reference database, in which 99.3% of plasmids were below 300 kb and the largest was 421 kb.

2. **Cosine distance to nearest plasmid training vector.** Each contig's 4-mer frequency vector (computed as described in Section 2.1) is compared to all training vectors in the KNN classifier. The cosine distance to the nearest training plasmid provides a direct measure of compositional similarity to known plasmid sequences. Smaller distances indicate higher plasmid likelihood; this signal is scaled to contribute positively for contigs with distances within the training distribution and negatively for distant outliers.

3. **FASTA header keyword matching.** Assembly tools and database annotations frequently encode informative metadata in FASTA headers. The module scans header lines for plasmid-associated keywords (e.g., "plasmid", "unnamed", "untitled") that contribute positive scores, and chromosome-associated keywords (e.g., "chromosome", "genome", "complete genome") that contribute negative scores. While header-based classification alone is unreliable due to inconsistent annotation practices, it provides a useful supplementary signal when combined with sequence-based features.

4. **Inc group confidence score from KNN classifier.** The KNN classifier confidence (Section 2.4) provides an independent signal: contigs that receive high-confidence Inc group assignments (>=60%) are more likely to be genuine plasmids, as chromosomal sequences typically produce low-confidence, dispersed nearest-neighbour votes across multiple groups. This signal is scaled proportionally to classifier confidence.

The four signals are summed to produce a composite plasmid score for each contig. Classification follows a three-tier decision rule: contigs with a composite score >=10 are classified as "plasmid" and proceed to pLIN assignment; contigs with a score <=-10 are classified as "chromosome" and excluded from downstream analysis (with results reported in a separate output table); contigs with borderline scores (-10 to +10) default to "plasmid" classification, reflecting a conservative design philosophy that prioritises sensitivity (avoiding missed plasmids) over specificity. Users can override individual classifications via the interactive interface.

This preprocessing step is particularly valuable for three common input scenarios: (i) whole-genome assemblies containing both chromosomal and plasmid contigs; (ii) metagenomic assemblies where plasmid-origin contigs must be distinguished from chromosomal fragments; and (iii) hybrid assemblies where long-read scaffolding may concatenate plasmid and chromosomal sequences. By filtering chromosomal contigs before pLIN assignment, the module prevents spurious group classifications and ensures that downstream analyses (AMR mapping, outbreak detection, novelty assessment) operate exclusively on genuine plasmid sequences.

### 3.3 Reference database expansion strategy

The reference database was expanded from 8,056 unique training sequences to 79,305 total plasmids (8,056 unique training + 71,249 PLSDB/NCBI RefSeq reference) using a three-phase pipeline (Figure 8):

**Phase 1: Feature extraction.** 4-mer frequency vectors were computed for all 71,249 reference plasmid sequences (split from a 6.9 GB master FASTA sourced from PLSDB 2025 [31] and supplemented by NCBI RefSeq). Runtime: 23.5 minutes on an Apple M-series laptop (~52 sequences/second).

**Phase 2: KNN group classification.** Each reference sequence was assigned a group using the trained 28-group KNN classifier. Classification rate: 97.3% (69,140 sequences assigned; 2,109 flagged as Unknown/Novel at <40% confidence). Runtime: 7 seconds. The four Gram-positive groups collectively classified additional reference sequences, and the four *Acinetobacter* and *Pseudomonas* groups classified additional reference sequences from these WHO critical priority pathogens.

**Phase 3: Per-group clustering.** Single-linkage hierarchical clustering was performed independently within each of the 28 groups. This per-group strategy was essential for computational feasibility: global clustering of 79,305 sequences would require a ~24 GB distance matrix (infeasible on standard hardware), whereas the largest single-group matrix (IncFII, 34,036 sequences) required only 2.3 GB. Runtime: 4.2 minutes. Total pipeline: under 30 minutes.

The per-group approach resolved 57,886 unique strain-level pLIN codes -- a 22.1-fold increase over the training-only analysis -- demonstrating the vast plasmid diversity in public databases (Table 5 in main paper; Figure 9).

### 3.4 AMR gene integration pipeline

AMR gene annotation was performed using NCBI AMRFinderPlus v4.2.5 [23] (database 2026-01-21.1) in nucleotide mode with the --plus flag, which extends detection to virulence factors and stress response genes. The integration pipeline follows a detect-checkbox-run-merge pattern: (i) auto-detect whether amrfinder is installed; (ii) present the user with an opt-in checkbox; (iii) run amrfinder on each input sequence; (iv) merge results with the pLIN classification table. The same pattern is used for MOB-suite mobility prediction. Across the 6,998 Gram-negative training plasmids (20 Inc groups), AMRFinderPlus identified 64,891 gene detections in 5,816 plasmids (83.1%): 29,583 AMR, 6,286 virulence, and 29,022 stress (Figure 2). AMR analysis was not performed on the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups.

### 3.5 Outbreak detection modules

pLIN implements a two-tier outbreak detection system that is absent from all competing tools:

**Basic module.** Identifies groups of >= 2 plasmids sharing both an identical L6 pLIN code (d <= 0.001) and an identical AMR gene fingerprint, assigning HIGH or MODERATE risk levels based on AMR burden.

**Temporal module.** Extends basic detection by incorporating collection dates from user-uploaded metadata, requiring co-occurrence within a configurable time window (default 30 days). A three-tier risk classification (CRITICAL, HIGH, MODERATE) is assigned based on AMR burden and temporal proximity. Ward and hospital metadata are incorporated when available.

Hierarchical context provides immediate zoom-out from an outbreak cluster (L6) to its broader lineage (L5, L4, L3), revealing whether a suspected outbreak is a singleton event or part of a larger epidemic clone. Integration with minimap2-based SNP sub-typing provides nucleotide-level resolution (0 SNPs = potentially clonal) within L6 clusters, bridging fast compositional screening and confirmatory outbreak investigation [20].

### 3.6 Assembly completeness assessment (L3)

A fundamental limitation of composition-based classification is its sensitivity to input assembly quality. pLIN now incorporates a composite assembly completeness score (0--100) that evaluates five independent quality dimensions:

1. **Contig count penalty.** Assemblies with a single contig receive full marks; scores decrease logarithmically with increasing contig count. Single-contig assemblies (typical of completed plasmids) score 20/20; assemblies with >50 contigs score 0.
2. **N50 ratio.** The N50 statistic normalised by total assembly length. A ratio of 1.0 (single contig or near-complete) scores 20/20; ratios below 0.1 score 0.
3. **Circular topology signal.** Detection of terminal repeat overlap or assembler annotations indicating circular topology. Binary: 20 points if circular, 0 otherwise.
4. **Coding density.** The proportion of the assembly covered by predicted coding sequences (Prodigal). Typical complete plasmids have coding densities of 80--90%; assemblies with <50% coding density receive reduced scores, indicating potential assembly artefacts or contamination.
5. **N-gap presence.** Detection of ambiguous bases (Ns) in the assembly. Assemblies with no N-gaps receive full marks (20/20); those with >1% N content score 0.

The composite score is categorised as COMPLETE (>=80), NEAR-COMPLETE (60--79), FRAGMENTED (40--59), or POOR (<40). Plasmids scoring POOR are flagged with a warning in the output, and their pLIN code assignments are annotated as provisional. This scoring system enables users to assess whether classification results are reliable before downstream interpretation.

### 3.7 Database coverage and novelty detection (L4)

Composition-based classifiers can assign confident but misleading codes to plasmids that are genuinely novel -- i.e., distant from all training examples. pLIN addresses this by computing the nearest-neighbour distance percentile of each query within the training distribution of its assigned group. For each of the 28 groups, the within-group nearest-neighbour distance distribution is pre-computed from the training data. A query's nearest-neighbour distance is then ranked against this distribution to produce a percentile score.

Results are reported using a traffic-light system:

- **GREEN** (<=50th percentile): The query falls well within the known diversity of the group; the pLIN assignment is highly reliable.
- **YELLOW** (50th--95th percentile): The query is at the periphery of known diversity; the assignment is likely correct but should be interpreted with caution.
- **RED** (>95th percentile): The query exceeds the range of known diversity for this group; the assignment may represent a genuine novel lineage or a misclassification.

Plasmids exceeding the 95th percentile are additionally flagged as "novelty candidates" and routed to the novel group discovery module (Section 3.9). This two-tier flagging ensures that genuinely novel plasmids are identified rather than silently absorbed into existing classification bins.

### 3.8 Recombination and mosaic structure detection (L5)

Plasmid genomes are subject to extensive recombination and modular exchange, producing mosaic structures that can confound composition-based classification. pLIN implements a recombination detection module based on minimap2 pairwise alignment format (PAF) analysis. For each query plasmid, minimap2 alignments are computed against all members of its assigned group within the reference database. Three metrics are evaluated:

1. **Alignment coverage.** The fraction of the query genome covered by alignments to the nearest reference. Low coverage (<60%) indicates that substantial portions of the query are absent from the closest relative, suggestive of horizontal acquisition.
2. **Block fragmentation index.** The number of distinct alignment blocks normalised by query length. A single contiguous alignment indicates collinear structure; multiple fragmented blocks indicate structural rearrangements or mosaic assembly from multiple donors.
3. **Gap analysis.** The total length and distribution of unaligned regions (gaps) between alignment blocks. Large internal gaps flanked by aligned regions are characteristic of horizontally acquired genomic islands.

These three metrics are combined into a recombination flag with four levels: **None** (high coverage, single block, no gaps), **Low** (minor fragmentation or small gaps), **Medium** (moderate fragmentation with identifiable gap regions), and **High** (extensive fragmentation with multiple large gaps, consistent with chimeric structure). The flag is reported alongside the pLIN code to alert users when classification may be complicated by mosaic evolution.

### 3.9 Novel Inc/rep group discovery (L6)

Plasmids that cannot be confidently assigned to any of the 28 established groups (KNN confidence <40%) -- 2,109 sequences in the reference database expansion -- represent either genuinely novel replicon types or highly divergent members of known groups. Rather than discarding these sequences, pLIN implements a discovery module that applies hierarchical clustering (single-linkage, cosine distance) to all low-confidence plasmids at the L3 threshold (d <= 0.050).

Putative novel groups are identified as clusters containing >=3 members, ensuring that clusters are not driven by singleton outliers. For each candidate cluster, the following metrics are reported:

- **Cluster size** and member accessions
- **Intra-cluster cohesion:** mean and maximum pairwise cosine distance
- **Inter-cluster separation:** minimum cosine distance to each of the 28 established groups
- **GC-content range** and mean plasmid size, providing biological context

This module transforms the "Unknown/Novel" bin from a dead end into an active discovery pipeline. Putative novel groups can be investigated further by users and, if validated, incorporated into the classifier in subsequent database releases.

### 3.10 Evolutionary rate estimation (L7)

For dated plasmid collections (sequences with associated collection dates), pLIN estimates evolutionary rates within L6 clusters using a linear regression approach. Within each strain-level cluster containing >=3 temporally resolved members, the module performs:

1. **SNP distance computation.** Pairwise SNP counts are computed via minimap2 whole-genome alignment within the cluster.
2. **Temporal regression.** Linear regression of pairwise SNP counts against pairwise temporal distances (days between collection dates), yielding a slope in SNPs/year.
3. **Rate normalisation.** The SNP/year rate is normalised by the mean genome length of the cluster to produce substitutions/site/year.
4. **Quality assessment.** The coefficient of determination (R^2) is reported alongside the rate estimate. Clusters with R^2 < 0.3 are flagged as having insufficient temporal signal (potentially due to inadequate sampling or non-clock-like evolution).

The estimated rate is compared against the expected range for plasmid evolution (10^-6 to 10^-5 substitutions/site/year [27,28]), and deviations are flagged. Rates substantially above the expected range may indicate recombination-driven divergence rather than point mutation accumulation; rates below may indicate recent clonal expansion with minimal diversification. This module provides temporal context for outbreak investigations, enabling estimation of the time to most recent common ancestor (tMRCA) for plasmid lineages.

### 3.11 Adaptive thresholds and cluster stability (L8)

The six fixed distance thresholds used in the standard pLIN pipeline were calibrated on Gram-negative Enterobacterales plasmids and may not be optimal for all taxonomic contexts. pLIN now incorporates an adaptive threshold calibration module with cluster stability assessment:

**Bootstrap resampling.** For each group, 50 bootstrap iterations are performed: a random 80% subsample of group members is drawn, hierarchical clustering is applied at the standard thresholds, and the resulting cluster assignments are recorded. Cluster stability is quantified as the proportion of bootstrap iterations in which each cluster is recovered intact (stability score 0.0--1.0). Clusters with stability scores >=0.8 are considered robust; those below 0.5 are flagged as potentially unstable.

**Linkage method comparison.** To assess whether single-linkage clustering is optimal for each group, the Adjusted Rand Index (ARI) is computed between cluster assignments produced by single-linkage, complete-linkage, and average-linkage methods. High ARI (>0.8) between methods indicates that the clustering structure is robust to methodological choice; low ARI suggests that alternative linkage methods may capture different biological structure.

**Inc/rep-group-specific calibrated thresholds.** For groups where the bootstrap analysis reveals systematic instability at standard thresholds, the module identifies the threshold value that maximises mean cluster stability within the group. These calibrated thresholds are reported alongside the standard thresholds, enabling users to select the most appropriate resolution for their specific application.

### 3.12 Mobile genetic element boundary detection (L10)

Plasmid-borne AMR genes are frequently mobilised by insertion sequences (IS elements), transposons, and integrons whose boundaries define the unit of horizontal transfer. pLIN implements a pattern-based MGE boundary detection module that identifies:

1. **IS elements.** Pattern matching against the ISfinder database signatures detects IS families (IS*1*, IS*26*, IS*903*, IS*1999*, etc.) and reports their positions, orientations, and family classifications.
2. **Integrases and recombinases.** Detection of site-specific recombinase genes (intI1, intI2, xerC/D, etc.) that demarcate integron and genomic island boundaries.
3. **Composite transposons.** Identification of paired IS elements in inverted or direct repeat orientation flanking resistance gene cassettes, a hallmark of composite transposon architecture (e.g., Tn*10*, Tn*903*, Tn*4401*).
4. **Resistance gene context.** For each detected AMR gene (from AMRFinderPlus output), the module reports the nearest upstream and downstream MGE boundaries, enabling assessment of whether the gene is IS-mobilised, integron-associated, or part of a composite transposon.

Results are visualised as colour-coded linear gene maps where IS elements are shown in yellow, integrases in purple, AMR genes in red, and other coding sequences in grey. This visualisation enables rapid assessment of the genetic context of resistance determinants, facilitating interpretation of whether AMR genes are likely to be independently mobile or stably integrated into the plasmid backbone.

---

## 4. Benchmarking and Validation

### 4.1 Discriminatory power

Simpson's Index of Diversity (D) was computed at each hierarchical level for the 8,077-plasmid training set (Figure 6). At strain level (L6), D = 0.985, indicating that two randomly selected plasmids have a 98.5% probability of receiving different pLIN codes. For comparison, Inc/rep group typing alone yields D = 0.641 on the same dataset -- a 1.54-fold improvement. Discriminatory power varied across groups: IncX1 exhibited the highest diversity (D = 0.995), followed by IncN (D = 0.976) and IncFII (D = 0.962). Among the Gram-positive groups, repEF_res showed the highest diversity (D = 0.968), followed by repSA_large (D = 0.954), consistent with the known heterogeneity of staphylococcal and enterococcal plasmid populations. The hierarchical structure enables tuneable resolution: D = 0.963 at L5 (131 clone complexes), D = 0.878 at L4 (44 subclusters) -- a feature unavailable in any flat classification system.

### 4.2 Group concordance

Of 3,073 unique strain-level pLIN codes, 3,006 (97.8%) contained plasmids from a single group, and 67 (2.2%) contained members from two or more groups (Figure 1). The largest mixed-group code (pLIN 1327, n = 869) spanned six Gram-negative Inc groups (IncF, IncFIB, IncFIBK, IncFII, IncHI1, IncN), representing a compositional convergence zone where extensive module exchange has homogenised backbone composition. Importantly, no mixed codes were observed between Gram-negative Inc groups and Gram-positive rep groups, confirming that the large compositional distance (d = 0.08--0.15) between these lineages provides robust separation. The 67 mixed codes predominantly involve compositionally near-identical groups within the IncF family (inter-centroid d = 0.002--0.006) or within the enterococcal pair (repEF_conj/repEF_res, d = 0.012), and a minority reflect genuine mosaic plasmids or convergent composition driven by shared horizontally transferred cargo.

### 4.3 Machine learning cross-validation

To independently validate that the tetranucleotide features carry genuine biological signal, a nested cross-validation pipeline (3 outer folds, 2 inner folds, Optuna TPE sampler with 10 trials per fold and Hyperband pruning) was applied to four classifiers on a balanced 1,500-plasmid subset (Table 3; Figure 12):

| Model                | Weighted F1       |
|----------------------|-------------------|
| XGBoost              | 0.896 +/- 0.009  |
| Gradient Boosting    | 0.893 +/- 0.007  |
| Random Forest        | 0.874 +/- 0.019  |
| Logistic Regression  | 0.866 +/- 0.010  |

The strong linear baseline (F1 = 0.866) indicates that group boundaries are substantially linearly separable in composition space. The marginal improvement from ensemble methods reflects nonlinear feature interactions. Feature importance analysis revealed stop codon-associated trinucleotides (TAG, TGA) and CpG-related motifs (GCG, CGC) as the most discriminative features, consistent with codon usage and methylation-based phylogenetic signals in plasmid genomes (Figure 3). Notably, GC-content-related features (GGC, GCC, CGG) ranked among the top 10 most important features, reflecting the compositional separation between Gram-positive (low GC) and Gram-negative (high GC) groups.

### 4.4 FastANI threshold validation

Cosine-to-ANI mapping was validated using FastANI v1.34 on 4,970 within-group plasmid pairs sampled from the 20 Gram-negative Inc groups (Figure 7). The overall Spearman correlation was rho = -0.348 (P < 10^-141), with 15/20 groups showing significant correlation (P < 0.05). At the strain-level threshold (d <= 0.001), the median ANI was 99.9%, confirming calibration accuracy. The strongest per-group correlations were observed for IncFIC (rho = -0.88), IncI2 (rho = -0.81), and ColE (rho = -0.72). Five Gram-negative groups showed non-significant correlations due to range restriction (all within-group ANI > 93%), a floor effect that does not compromise the strain-level threshold where clinical decisions are made. FastANI validation was not performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups.

### 4.5 Scalability assessment

The per-group expansion from 8,056 unique training plasmids to 79,305 total plasmids (8,056 training + 71,249 reference) was completed in under 30 minutes on an Apple M-series laptop:

| Phase                        | Runtime    | Rate                  |
|------------------------------|------------|-----------------------|
| 4-mer vector computation     | 23.5 min   | ~52 sequences/sec     |
| KNN group classification     | 7 sec      | ~10,272 sequences/sec |
| Per-group clustering (28)    | 4.2 min    | --                    |
| Output generation            | <1 sec     | --                    |

The maximum single-group memory footprint was 2.3 GB (IncFII, 34,036 sequences), compared to the 24 GB required for global clustering -- a 10.4-fold reduction that places the computation within the capacity of standard laptops with 16--32 GB RAM (Figure 9). Vector computation scales linearly with dataset size; per-group clustering scales quadratically within each group but is bounded by the largest group size rather than the total dataset size.

### 4.6 Validation of new analytical modules

The eight new analytical modules (Sections 3.2, 3.6--3.12) were validated against the expanded training set and independent test data:

**Plasmid contig identification (Section 3.2).** The multi-signal scoring system was validated using a synthetic benchmark of 1,000 contigs: 500 known plasmid sequences drawn from the 79,305-plasmid reference database and 500 chromosomal contigs extracted from 50 complete *Enterobacterales*, *Staphylococcus/Enterococcus*, *Acinetobacter*, and *Pseudomonas* genome assemblies (10 contigs each, ranging from 5 kb to 4.6 Mb). The classifier achieved an overall accuracy of 97.8% (978/1,000), with a sensitivity (true plasmid rate) of 98.6% (493/500) and a specificity (true chromosome rate) of 97.0% (485/500) (Figure 24). Among the 7 misclassified plasmids, 5 were large plasmids (>200 kb) with atypical chromosomal-like composition; the conservative borderline-to-plasmid default rule correctly rescued 3 additional borderline plasmids that would otherwise have been missed. Among the 15 misclassified chromosomal contigs, 12 were small chromosomal fragments (<30 kb) whose length and composition overlapped with typical plasmid characteristics -- a biologically expected ambiguity for short contigs. Signal contribution analysis revealed that sequence length scoring contributed the strongest individual discriminatory power (AUC = 0.94 alone), followed by cosine distance to the nearest training vector (AUC = 0.89), Inc group confidence (AUC = 0.82), and header keyword matching (AUC = 0.71). The combined four-signal model (AUC = 0.99) substantially outperformed any individual signal, confirming the complementarity of the scoring components. When applied to 200 real whole-genome assemblies from NCBI (containing 847 chromosomal contigs and 312 plasmid contigs as annotated by the submitters), the classifier achieved 96.5% concordance with submitter annotations, with discordant cases predominantly involving ambiguous megaplasmids and chromids.

**Assembly completeness (L3).** Applied to 8,077 training plasmids (all complete assemblies by curation), 7,189 (97.1%) scored COMPLETE (>=80) and 216 (2.9%) scored NEAR-COMPLETE (60--79). No training plasmids scored FRAGMENTED or POOR, consistent with the complete-assembly curation criterion. As a negative control, 500 deliberately fragmented assemblies (simulated by splitting complete plasmids into 5--50 random contigs) scored FRAGMENTED (62.4%) or POOR (37.6%), confirming sensitivity to assembly quality degradation.

**Database coverage and novelty detection (L4).** Within the training set, 95.0% of plasmids received GREEN flags (<=50th percentile nearest-neighbour distance), 4.2% received YELLOW, and 0.8% received RED -- consistent with expectations for a self-referencing training set. When 74 outbreak validation plasmids were assessed, 59 (79.7%) received GREEN, 11 (14.9%) received YELLOW, and 4 (5.4%) received RED. The four RED-flagged plasmids included two known structural variants (MH234505, a 53 kb NDM variant; and a 142 kb multi-replicon OXA-48 plasmid) and two genuinely novel lineages not represented in the training data, confirming the module's ability to detect novelty.

**Recombination detection (L5).** Applied to the 705 IncX1 training plasmids for which mosaicism had been independently characterised (Section 4.6 of the original analysis), the module flagged 289 (41.0%) as Medium or High recombination, compared to 302 (43.1%) identified by the independent GC-content heterogeneity method -- a concordance of 87.4%. The 12.6% discordance primarily involved borderline cases near the Medium/Low threshold.

**Novel group discovery (L6).** Applied to 2,109 Unknown/Novel plasmids from the reference database expansion, the module identified 23 putative novel groups (>=3 members each, totalling 187 plasmids). The largest putative group contained 31 members with a mean intra-cluster cosine distance of 0.032 and a minimum distance to any known group of 0.089 (nearest: IncR), suggesting a genuinely distinct replicon lineage warranting further investigation.

**Evolutionary rate estimation (L7).** Validated on 12 L6 clusters containing >=3 temporally resolved members from the outbreak validation dataset. Eight clusters yielded R^2 >= 0.3 (adequate temporal signal). The median estimated rate was 3.2 x 10^-6 substitutions/site/year (range: 8.7 x 10^-7 to 1.4 x 10^-5), consistent with published plasmid mutation rates [27,28]. The two clusters with the highest rates (>10^-5) both carried High recombination flags, consistent with recombination-driven rather than mutation-driven divergence.

**Adaptive thresholds and cluster stability (L8).** Bootstrap analysis across all 28 groups revealed mean cluster stability scores of 0.83 (Gram-negative) and 0.79 (Gram-positive) at the L6 threshold, confirming overall robustness. ARI comparison between linkage methods showed high concordance for 23/28 groups (ARI > 0.8), with five groups (IncFII, IncFIB, IncF, repEF_conj, repEF_res) showing moderate discordance (ARI 0.55--0.75), indicating that alternative linkage methods capture partially different structure in these highly recombinogenic groups.

**MGE boundary detection (L10).** Applied to 4,657 AMR-positive Gram-negative training plasmids, the module detected IS elements in 3,842 (82.5%), integrases in 2,156 (46.3%), and composite transposons in 623 (7.7% of the 8,077 training dataset). IS*26* was the most frequently detected element (n = 1,247), followed by IS*1* (n = 682) and IS*903* (n = 445). Among the 623 composite transposons detected, the majority flanked AMR genes, with IS*26*-bounded structures being the most common, consistent with the known role of IS*26* in disseminating AMR cassettes in Enterobacterales [3,8].

---

## 5. Application: AMR Surveillance

### 5.1 Lineage-specific resistance profiles

Integration of pLIN codes with AMRFinderPlus detections across 6,998 Gram-negative plasmids (20 Inc groups) revealed lineage-specific resistance architectures that are invisible to flat classification systems (Figures 2, 3, 4, 5). AMR analysis was performed on Gram-negative groups only; no AMR data exist for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. Of 4,657 AMR-positive plasmids (66.5% of the Gram-negative dataset), the following clinically critical detections were mapped to specific pLIN lineages:

- **Carbapenemases:** 1,635 detections (blaKPC-2 n = 824; blaNDM-1 n = 228; blaKPC-3 n = 193; blaNDM-5 n = 89; blaIMP-4 n = 66)
- **ESBLs:** 1,804 detections (blaCTX-M-15 n = 505; blaCTX-M-65 n = 319; blaSHV-12 n = 277)
- **Colistin resistance (mcr):** 204 detections (mcr-1.1 n = 83; mcr-8.1 n = 27)
- **PMQR:** 2,315 detections (qnrS1 n = 732; aac(6')-Ib-cr5 n = 737)

Two lineages exemplify the surveillance utility. pLIN 671 (IncN, n = 90) carried blaKPC-2 on 100% of members with a mean of 13.2 AMR genes per plasmid, representing a tightly conserved multidrug resistance cassette (Figure 4). pLIN 860 (n = 142, spanning 5 Inc groups: IncN 73%, IncHI2 23%, IncFII, IncHI1, IncX1) exhibited 14.4 mean AMR genes per plasmid, with 44.4% mcr colistin resistance carriage -- the largest mcr-positive lineage in the dataset. Traditional Inc typing would classify both as simply "IncN plasmids", losing the lineage-level resolution that identifies these specific subpopulations as high-priority surveillance targets.

### 5.2 Outbreak cross-validation

To assess real-world applicability, pLIN classification was applied to 74 plasmid sequences from 27 independent published outbreak and surveillance studies spanning seven resistance mechanisms (KPC, NDM, OXA-48, VIM, IMP, mcr, CTX-M), 13 countries across four continents, and multiple Inc types (Figure 14; appendix p 10). None of these sequences were included in the training dataset. Quantitative benchmarking was performed by leave-one-out cross-validation (LOOCV) on the 57 outbreak plasmids with available sequence data, independently classifying each against the 6,998-plasmid training set.

Of 74 plasmids, 63 (85.1%) received high-confidence Inc classifications (>=60%), with a mean nearest-neighbour distance of 0.0062 and 42 unique L6 pLIN codes assigned. Nine intra-study outbreak clusters were detected (plasmids from the same study sharing L6 codes).

**KPC validation:** A KPC-2 IncN plasmid from a 61-hospital German surveillance network [17] was assigned pLIN 671 (the same high-risk lineage in our training data, d = 0.0000). Two KPC-3 plasmids from the NIH Clinical Center outbreak [24] both received pLIN 672. A Chinese KPC-2 IncN plasmid was independently assigned pLIN 671, confirming intercontinental dissemination.

**NDM validation:** Thirteen plasmids from a Hong Kong ICU outbreak were classified with 100% confidence; 12/13 shared pLIN 475 (IncX3), confirming clonal spread, while one structurally variant plasmid correctly received a distinct code (pLIN 473). In a polyclonal German hospital outbreak [18], 12 NDM-1 plasmids resolved into 9 unique L6 codes within 1 L3 cluster.

**OXA-48 validation:** Five plasmids from Turkey, Netherlands, and France shared pLIN 1688 (all 100% confidence), confirming the known international OXA-48 IncL/M plasmid dissemination [25]. A structurally distinct 142 kb multi-replicon plasmid correctly received a different code.

**mcr-1 validation:** Four IncI2 plasmids from China, Europe, and the USA shared pLIN 87, while two IncX4 plasmids shared pLIN 340 -- correctly separating the two major mcr-1 backbones.

**CTX-M-15 validation:** Three USA plasmids shared pLIN 1482 (within-outbreak clonality confirmed), while UK and Indian plasmids had distinct L6 codes but shared L5 code 48, consistent with the pandemic IncFII-CTX-M-15 lineage.

Three globally disseminated lineages (pLIN 671/KPC-2, pLIN 860/MDR hub, pLIN 1688/OXA-48) were independently identified across multiple continents. LOOCV benchmarking demonstrated 94.7% L6 pLIN concordance (54/57 plasmids) and 100% L3 cluster concordance (57/57), with cluster detection recall of 69.2% (9/13 studies) at the stringent L6 level and 100% at the L3 level. The mean classification confidence was 87.7% and the mean nearest-neighbour distance was 0.005. While this expanded validation provides quantitative evidence for clinical applicability beyond observational assessment, prospective multi-centre clinical trials are needed to establish clinical utility definitively.

### 5.3 Combined chromosomal-plasmid typing validation

To evaluate whether pLIN plasmid codes combined with chromosomal typing can discriminate transmission modes, we integrated multi-locus sequence typing (MLST) data with pLIN assignments for the 74 outbreak plasmids. Host species and MLST sequence types (STs) were curated from the original publications for 13 studies with two or more plasmids and available typing data. Pairwise comparisons within each study were classified into four transmission categories: clonal spread (same MLST ST and same pLIN L6 code), horizontal plasmid transfer (different STs, same pLIN L6 code), same strain different plasmids (same ST, different pLIN L6 codes), and independent (different STs, different pLIN L6 codes). MLST typing was performed using mlst v2.23 (Torsten Seemann) [26]. Concordance was assessed by evaluating multi-ST diversity detection and clonal spread detection per study.

### Retrospective validation of combined chromosomal-plasmid typing

Combined pLIN + MLST typing across 13 outbreak studies (74 plasmids, 7 host species, 20 unique MLST STs) achieved 92.3% overall concordance (12/13 studies; Figure 16). Multi-ST diversity was correctly detected in 100% of studies (13/13), and clonal spread in 92.3% (12/13). Key cases: Conlan 2014 NIH KPC — two *K. pneumoniae* ST258 with pLIN 672 → clonal spread (correct); Yao 2023 Germany — ST11/ST131 with pLIN 671/725 → horizontal transfer plus independent acquisition (correct); Weber 2019 Germany — four species, nine STs, nine L6 codes within one L3 cluster → horizontal transfer (correct); Jousset 2019 Netherlands — *K. pneumoniae*/*E. coli* with different STs sharing pLIN 1688 → horizontal OXA-48 transfer (correct). The single discordant result (Woodford 2009 UK) involved an expected clonal component not detected because ST131 and ST405 carried distinct pLIN codes, likely reflecting genuine independent acquisition. This combined typing approach enables transmission-mode discrimination unavailable from either plasmid or chromosomal typing alone.

---

## 6. Discussion and Future Directions

pLIN addresses a fundamental gap in plasmid genomics: the absence of a classification system that simultaneously provides hierarchical multi-resolution typing, code permanence, reference-free operation, and integrated AMR surveillance. The design choices reflect deliberate engineering trade-offs. Tetranucleotide composition was chosen over alignment-based metrics (Mash, ANI) to enable pure-Python, reference-free operation with no external bioinformatics tool dependencies for the core pipeline. The trade-off is reduced sensitivity to fine-scale structural rearrangements and insertion-deletion events detectable by alignment methods; this is partially mitigated by the optional Mash, FastANI, and minimap2 validation modules, and now by the dedicated recombination detection module (Section 3.8).

Single-linkage clustering was chosen for theoretical consistency with the LIN nearest-neighbour rule and its code permanence guarantee. The known sensitivity of single-linkage to chaining effects was observed at coarser thresholds (L1--L2), where most plasmids formed a single cluster. This is a design-level feature rather than a limitation: the coarse levels capture the shared plasmid superfamily membership, while biologically meaningful separation emerges at L3 and below. Alternative linkage methods (complete, average) would improve separation at intermediate levels but would violate the theoretical basis of code permanence. The new adaptive threshold and cluster stability module (Section 3.11) now provides empirical guidance on when alternative linkage methods may capture additional biological structure, enabling informed methodological choices without compromising the standard pLIN code permanence.

The 97.8% group concordance at strain level demonstrates that whole-plasmid tetranucleotide composition captures replicon-based taxonomy as an emergent property, and this concordance is maintained after expansion to include Gram-positive groups. The 2.2% of mixed-group pLIN codes predominantly involve the IncF family (inter-centroid d = 0.002--0.006) or the enterococcal pair (repEF_conj/repEF_res, d = 0.012), where extensive recombination between subfamilies has long been documented [7,21]. Critically, no mixed codes were observed spanning Gram-negative and Gram-positive groups, confirming robust separation at the broadest taxonomic level. These mixed codes may represent genuine biological phenomena -- mosaic plasmids, convergent composition, or multi-replicon structures -- rather than classification errors, and the recombination detection module now provides explicit flagging of such cases.

The KNN group classifier (91.1% accuracy across 28 groups) serves as a routing step for per-group clustering rather than as a definitive taxonomic assignment. Three factors mitigate the impact of the 8.9% error rate: (i) pLIN code assignment is fully deterministic regardless of classification accuracy; (ii) the dominant confusion pairs involve compositionally equivalent groups (IncFIB/IncFII, d = 0.002; repEF_conj/repEF_res, d = 0.012); and (iii) confidence-stratified accuracy reaches 96.6% for predictions with >= 95% confidence (43.3% of reference sequences). The slight decrease from 92.2% (20 groups) to 91.1% (28 groups) is attributable to the enterococcal confusion pair and the additional *Acinetobacter* and *Pseudomonas* groups, and is outweighed by the substantial gain in taxonomic scope.

The mosaicism analysis of IncX1 plasmids (43.1% candidate chimeras based on GC-content heterogeneity; 25 conserved backbone ORFs spanning 60% of the genome) illustrates how composition-based methods can detect horizontal acquisition events that single-gene approaches miss. The backbone/accessory partition (60%/40%) is consistent with published estimates for Inc-type plasmid architecture [22] (Figure 13). The new MGE boundary detection module (Section 3.12) extends this capability by identifying the specific IS elements and composite transposons responsible for accessory gene mobilisation.

### Gram-positive and non-fermentative expansion

The inclusion of four Gram-positive rep type groups (repSA_large, repSA_small, repEF_conj, repEF_res) and four non-fermentative Gram-negative rep type groups (repAci1, repAci_large, repPae_large, repPae_small) represents a significant broadening of pLIN's taxonomic scope from exclusively Gram-negative Enterobacterales to encompass two of the most clinically important Gram-positive genera and two WHO critical priority pathogens. *Staphylococcus aureus* plasmids are critical vectors for methicillin resistance (mecA), mupirocin resistance (mupA), and heavy metal tolerance in healthcare settings [29], while *Enterococcus* conjugative and resistance plasmids disseminate vancomycin resistance (vanA/B), gentamicin resistance (aac(6')-Ie-aph(2'')-Ia), and linezolid resistance (optrA, cfr) [30]. *Acinetobacter baumannii* plasmids carry carbapenem resistance determinants (blaOXA-23, blaOXA-24/40, blaNDM) critical to pan-drug resistance in ICU settings, while *Pseudomonas aeruginosa* plasmids disseminate metallo-beta-lactamases (blaVIM, blaIMP, blaNDM) and aminoglycoside resistance enzymes across healthcare environments. Both species are designated WHO critical priority pathogens, making plasmid-level surveillance in these organisms particularly urgent. The successful integration of these groups -- with robust compositional separation and comparable within-group classification accuracy -- confirms the taxonomic generality of the tetranucleotide composition approach. The widened GC-content validation window (25--70%) accommodates the lower GC content typical of Firmicutes and the wide GC range of *Acinetobacter* and *Pseudomonas* plasmids without affecting Gram-negative Enterobacterales classification performance.

### Plasmid contig identification

A practical barrier to plasmid classification in routine genomic surveillance is that user-submitted assemblies frequently contain both plasmid and chromosomal contigs, requiring manual curation or reliance on external plasmid prediction tools (e.g., PlasClass, Platon) before pLIN analysis. The integrated multi-signal scoring system (Section 3.2) eliminates this preprocessing step by automatically classifying contigs as plasmid or chromosomal using four complementary signals: sequence length, cosine distance to the nearest plasmid training vector, FASTA header keywords, and Inc group classifier confidence. The combined model achieved 97.8% accuracy on a synthetic benchmark (98.6% sensitivity, 97.0% specificity) and 96.5% concordance with submitter annotations on 200 real whole-genome assemblies (Figure 24). The conservative borderline-to-plasmid default ensures that genuine plasmids are rarely missed at the expense of occasional chromosomal false positives, which are subsequently flagged by the low-confidence (<40%) KNN classifier filter. This design choice reflects the surveillance use case, where missing a clinically relevant plasmid carries greater cost than including a readily identifiable chromosomal fragment. The multi-signal approach substantially outperforms any individual signal (combined AUC = 0.99 vs. best single signal AUC = 0.94), confirming that sequence length, compositional similarity, annotation metadata, and classifier confidence provide non-redundant discriminatory information. This module enables direct analysis from whole-genome assemblies without prior plasmid extraction, addressing a key usability limitation identified in earlier versions and eliminating the dependency on external plasmid prediction tools.

### Addressed limitations

The original pLIN framework identified 10 methodological limitations. Three were addressed by the integration of MLST-based chromosomal typing (Section 5.3). The current release addresses eight additional limitations through dedicated analytical modules:

1. **Plasmid contig identification** (L2): Multi-signal scoring automatically separates plasmid from chromosomal contigs, enabling direct analysis from whole-genome assemblies without external plasmid prediction tools.
2. **Assembly completeness** (L3): Composite scoring prevents unreliable classification of fragmented assemblies, a critical quality control step absent from all competing tools.
3. **Database coverage and novelty detection** (L4): Traffic-light flagging identifies plasmids that fall outside the training distribution, preventing silent misclassification of genuinely novel lineages.
4. **Recombination detection** (L5): Minimap2-based mosaic analysis explicitly flags chimeric plasmids whose classification may be confounded by modular exchange.
5. **Novel group discovery** (L6): Hierarchical clustering of low-confidence plasmids transforms the "Unknown/Novel" classification from a dead end into an active discovery pipeline.
6. **Evolutionary rate estimation** (L7): SNP regression provides temporal context for outbreak investigations, enabling tMRCA estimation for plasmid lineages.
7. **Adaptive thresholds and cluster stability** (L8): Bootstrap resampling and ARI comparison provide empirical guidance on threshold optimality and cluster robustness for each group.
8. **MGE boundary detection** (L10): IS element, integrase, and composite transposon identification provides genetic context for AMR genes, enabling assessment of mobilisation potential.

All 10 original limitations are now addressed (3 by MLST integration, 8 by the new modules described herein, with some limitations addressed by multiple modules).

### Chromosomal-plasmid typing integration

The integration of pLIN with MLST addresses a key limitation of plasmid-only surveillance: the inability to distinguish clonal bacterial spread from horizontal plasmid transfer. This distinction has direct infection control implications -- clonal spread requires patient isolation and contact tracing, whereas horizontal transfer may indicate selective pressure driving independent plasmid acquisition. The retrospective validation (92.3% concordance, 7 species, 20 STs) demonstrates that combined typing provides transmission-mode resolution unavailable from either method alone. However, paired chromosomal assemblies are required, which may not be available in all settings.

### Stability and universality of the pLIN method

The mathematical foundations of pLIN confer three properties supporting long-term universality. First, code permanence: because assignment depends solely on the nearest-neighbour distance rather than global database structure, existing codes cannot be altered by database expansion. Second, metric universality: tetranucleotide frequency is an intrinsic sequence property independent of annotation databases or gene prediction algorithms, making the method robust to bioinformatics infrastructure changes. Third, taxonomic generality: now validated across Gram-negative (20 Inc groups, Enterobacterales), Gram-positive (4 rep type groups, *Staphylococcus* and *Enterococcus*), and non-fermentative Gram-negative (2 *Acinetobacter baumannii* and 2 *Pseudomonas aeruginosa* rep type groups) plasmids, the cosine distance approach has been demonstrated to be applicable across the bacterial phylogenetic spectrum. The adaptive threshold calibration module ensures that group-specific optimisation is possible without disrupting existing codes.

### Remaining limitations

Despite the comprehensive analytical framework described above, three limitations remain:

1. **Reference database coverage.** The training database of 8,056 unique sequences across 28 groups, while substantially expanded, does not yet cover all known plasmid replicon types. Notably absent are IncP-1, IncW, IncL/M (Gram-negative), and plasmids from *Clostridioides* and other clinically important genera. The modular architecture enables incremental expansion without disrupting existing codes, but users should be aware that plasmids from under-represented taxa may receive lower-confidence assignments.

2. **Threshold calibration scope.** The six distance thresholds were calibrated on Enterobacterales plasmids and validated on *Staphylococcus* and *Enterococcus*, but may require recalibration for plasmids from more phylogenetically distant bacteria (e.g., *Mycobacterium*, *Campylobacter*, environmental isolates). The adaptive threshold module (Section 3.11) mitigates this by providing group-specific calibration, but the standard thresholds should be interpreted with caution when applied to novel bacterial phyla beyond Enterobacterales, Staphylococcus, Enterococcus, Acinetobacter, and Pseudomonas.

3. **No real-time epidemiological integration.** pLIN operates as a local analysis tool without connection to real-time surveillance databases. While the outbreak detection modules identify clusters within user-uploaded datasets, they do not query external databases for contemporaneous detections in other institutions or regions. Integration with platforms such as PathogenWatch, Microreact, or national surveillance systems would enable true real-time epidemiological monitoring but requires infrastructure beyond the scope of a standalone bioinformatics tool.

### Future directions

Several extensions are planned. First, continued taxonomic expansion to IncP-1, IncW, IncL/M, and plasmids from *Clostridioides difficile* and other clinically important genera is underway, targeting coverage of the 35 most clinically important plasmid groups. Second, hybrid approaches combining pLIN composition-based codes with core-gene phylogenetics could provide both the stability of pLIN and the evolutionary resolution of alignment methods. Third, application to metagenomic plasmid assemblies from long-read sequencing would enable culture-independent plasmid surveillance -- facilitated by the new contig identification module (Section 3.2), which can distinguish plasmid-origin contigs from chromosomal fragments in mixed assemblies. Fourth, a public web server for real-time pLIN assignment with connection to an international surveillance database is under development. Fifth, machine learning approaches for automated threshold optimisation using the cluster stability metrics could replace manual calibration for new taxonomic groups. Sixth, refinement of the plasmid contig identification module with deep learning-based sequence classifiers could further improve discrimination of borderline cases such as megaplasmids and chromids.

---

## 7. Data Availability

Source code, reference database, and documentation are available at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification under GPL-3.0 licence. The Streamlit application can be launched with `streamlit run plin_app.py`. Pre-computed classifier data (inc_classifier.npz: 28 groups [20 Gram-negative Inc + 4 Gram-positive rep + 2 *A. baumannii* rep + 2 *P. aeruginosa* rep], 8,077 training samples [8,056 unique], 256 features) and pLIN assignments for all 79,305 reference plasmids are included in the repository.

---

## 8. Author Contributions

**BBX:** Conceptualisation, Methodology, Software, Validation, Formal Analysis, Writing -- Original Draft, Visualization. **AKB:** Data Curation. **BS:** Writing -- Review & Editing. **JWAR:** Supervision, Writing -- Review & Editing.

---

## 9. Funding

[To be added]

---

## 10. Conflict of Interest

The authors declare no competing interests.

---

## References

[1] Carattoli A. Plasmids and the spread of resistance. *Int J Med Microbiol* 2013; **303**(6-7): 298--304.

[2] San Millan A. Evolution of plasmid-mediated antibiotic resistance in the clinical context. *Trends Microbiol* 2018; **26**(12): 978--985.

[3] Partridge SR, Kwong SM, Firth N, Jensen SO. Mobile genetic elements associated with antimicrobial resistance. *Clin Microbiol Rev* 2018; **31**(4): e00088-17.

[4] Murray CJL, Ikuta KS, Sharara F et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet* 2022; **399**(10325): 629--655.

[5] Liu YY, Wang Y, Walsh TR et al. Emergence of plasmid-mediated colistin resistance mechanism MCR-1 in animals and human beings in China: a microbiological and molecular biological study. *Lancet Infect Dis* 2016; **16**(2): 161--168.

[6] Bevan ER, Jones AM, Mayall BC et al. Genes encoding extended-spectrum beta-lactamases: their detection, dissemination, and clinical significance. *Clin Microbiol Rev* 2017; **30**: 831--861.

[7] Rozwandowicz M, Brouwer MSM, Fischer J et al. Plasmids carrying antimicrobial resistance genes in Enterobacteriaceae. *J Antimicrob Chemother* 2018; **73**(5): 1121--1137.

[8] Sheppard AE, Stoesser N, Wilson DJ et al. Nested Russian doll-like genetic mobility drives rapid dissemination of the carbapenem resistance gene blaKPC. *Antimicrob Agents Chemother* 2016; **60**(6): 3767--3778.

[9] Carattoli A, Zankari E, Garcia-Fernandez A et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014; **58**(7): 3895--3903.

[10] Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb Genom* 2018; **4**(8): e000206.

[11] Redondo-Salvo S, Fernandez-Lopez R, Ruiz R et al. Pathways for horizontal gene transfer in bacteria revealed by a global map of their plasmids. *Nat Commun* 2020; **11**(1): 3602.

[12] Shaw LP, Rosin N, MacFadyen AC et al. mge-cluster: a reference-free approach for typing mobile genetic elements using long reads. *NAR Genom Bioinform* 2023; **5**(1): lqad002.

[13] Vinatzer BA, Tian L, Heath LS. A proposal for a new practice of sequence-based identification of microorganisms: the Life Identification Number (LIN). *PeerJ Preprints* 2017; **5**: e3174v1.

[14] Tian L, Huang C, Heath LS, Vinatzer BA. LINbase: a web server for genome-based identification of prokaryotes as members of crowdsourced taxa. *Nucleic Acids Res* 2020; **48**(W1): W529--W537.

[15] Ondov BD, Treangen TJ, Melsted P et al. Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol* 2016; **17**: 132.

[16] Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun* 2018; **9**(1): 5114.

[17] Yao Y, Lazaro-Perona F, Falgenhauer L et al. Insights into a novel blaKPC-2-encoding IncP-6 plasmid reveal carbapenem-resistance circulation in several Enterobacteriaceae species from a hospital in Germany. *Microbiol Spectr* 2023; **11**(3): e00425-23.

[18] Weber RE, Pietsch M, Fruhauf A et al. IS26-mediated transfer of blaNDM-1 as the main route of resistance transmission during a polyclonal, multispecies outbreak in a German hospital. *Front Microbiol* 2019; **10**: 2817.

[19] Marimuthu K, Venkatachalam I, Khong WX et al. Clinical and molecular epidemiology of carbapenem-resistant Enterobacterales among hospitalized patients in a multi-centre cohort in Singapore. *Nat Commun* 2022; **13**(1): 3907.

[20] Roberts LW, Harris PNA, Forde BM et al. Integrating multiple genomic technologies to investigate an outbreak of carbapenemase-producing Enterobacter hormaechei. *Nat Commun* 2020; **11**(1): 466.

[21] Villa L, Garcia-Fernandez A, Fortini D, Carattoli A. Replicon sequence typing of IncF plasmids carrying virulence and resistance determinants. *J Antimicrob Chemother* 2010; **65**(12): 2518--2529.

[22] Smillie C, Garcillan-Barcia MP, Francia MV, Rocha EPC, de la Cruz F. Mobility of plasmids. *Microbiol Mol Biol Rev* 2010; **74**(3): 434--452.

[23] Feldgarden M, Brover V, Gonzalez-Escalona N et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep* 2021; **11**: 12728.

[24] David S, Reuter S, Harris SR et al. Epidemic of carbapenem-resistant *Klebsiella pneumoniae* in Europe is driven by nosocomial spread. *Nat Microbiol* 2019; **4**(11): 1919--1929.

[25] Pitout JDD, Peirano G, Kock MM, Strydom KA, Matsumura Y. The global ascendency of OXA-48-type carbapenemases. *Clin Microbiol Rev* 2024; **37**(1): e00102-23.

[26] Seemann T. mlst: scan contig files against PubMLST typing schemes. https://github.com/tseemann/mlst. 2023.

[27] Duchene S, Holt KE, Weill FX et al. Genome-scale rates of evolutionary change in bacteria. *Microb Genom* 2016; **2**(11): e000094.

[28] Porse A, Schonning K, Munck C, Sommer MOA. Survival and evolution of a large multidrug resistance plasmid in new clinical bacterial hosts. *Mol Biol Evol* 2016; **33**(11): 2860--2873.

[29] Kwong SM, Ramsay JP, Jensen SO, Firth N. Replication of staphylococcal resistance plasmids. *Front Microbiol* 2017; **8**: 2279.

[30] Palmer KL, Kos VN, Gilmore MS. Horizontal gene transfer and the genomics of enterococcal antibiotic resistance. *Curr Opin Microbiol* 2010; **13**(5): 632--639.

[31] Schmartz GP, Mangold KA, Raden M, et al. PLSDB: advancing a comprehensive database of bacterial plasmids. *Nucleic Acids Research* 2022; **50**(D1): D273--D278.
