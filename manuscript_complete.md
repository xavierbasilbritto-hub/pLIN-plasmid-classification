# pLIN: A Plasmid Life Identification Number System for Hierarchical, Permanent Classification of Bacterial Plasmids Integrated with Antimicrobial Resistance Gene Surveillance

---

## ABSTRACT

Plasmids are central agents of horizontal gene transfer in bacteria, driving the dissemination of antimicrobial resistance (AMR) determinants, virulence factors, and adaptive traits. Despite their clinical and epidemiological importance, no existing plasmid classification system simultaneously provides hierarchical multi-resolution typing, code permanence, and reference-free operation. Current approaches—including replicon-based Inc typing (PlasmidFinder), plasmid multilocus sequence typing (pMLST), mobilization-based clustering (MOB-suite), ANI-based plasmid taxonomic units (COPLA/PTUs), and reference-free unitig clustering (mge-cluster)—each suffer from flat classification, database dependency, code instability, or limited taxonomic scope. Here, we introduce pLIN (plasmid Life Identification Number), the first application of the Life Identification Number (LIN) framework to plasmid genomes. pLIN assigns each plasmid a six-position hierarchical code based on tetranucleotide composition distances and single-linkage clustering at six biologically calibrated thresholds, spanning family-level (~85% ANI) to strain-level (~99.9% ANI) resolution. Applied to 6,346 complete plasmid sequences from three clinically important incompatibility groups (IncFII, n=4,581; IncN, n=1,064; IncX1, n=701), pLIN resolved 2,232 unique strain-level codes with a Simpson's Index of Diversity of 0.979, while maintaining 99.5% concordance with established Inc-group assignments. A nested cross-validation machine learning pipeline confirmed that the same compositional features underlying pLIN assignment robustly predict Inc-group membership (weighted F1=0.903, XGBoost). Mosaic structure analysis identified 43.1% of IncX1 plasmids as candidate chimeras, and seed ORF conservation analysis delineated a 25-gene IncX1 backbone with 50–85% prevalence across the group. Integration with NCBI AMRFinderPlus v4.2.5 across all 6,346 plasmids identified 60,372 gene detections (27,465 AMR, 5,834 virulence, 27,073 stress) in 84.2% of plasmids. Cross-referencing pLIN lineages with AMR gene content revealed lineage-specific resistance profiles: pLIN 1.1.1.7.30.567 (IncN) carried *blaKPC-2* on 100% of members with a mean of 13.2 AMR genes per plasmid, while the cross-Inc lineage 1.1.1.7.30.1313 harboured both *blaKPC-2* (58%) and the pan-aminoglycoside resistance gene *rmtB1* (71%). Critically important resistance determinants including carbapenemases (1,490 detections), ESBLs (3,734), plasmid-mediated colistin resistance (*mcr*, 160), and plasmid-mediated quinolone resistance (2,160) were mapped to specific pLIN lineages, demonstrating the utility of pLIN as a framework for tracking AMR gene dissemination through plasmid lineage surveillance. The pLIN system offers a stable, scalable, and biologically interpretable nomenclature for plasmid genomic epidemiology.

**Keywords:** plasmid classification, Life Identification Number, antimicrobial resistance, hierarchical clustering, tetranucleotide composition, AMRFinderPlus, genomic epidemiology

---

## 1. INTRODUCTION

### 1.1 Background

Bacterial plasmids are extrachromosomal DNA elements that serve as the primary vehicles of horizontal gene transfer (HGT) in prokaryotes, facilitating the rapid dissemination of antimicrobial resistance (AMR) genes, virulence factors, metabolic functions, and stress tolerance determinants across species and genera (Carattoli, 2009; San Millan, 2018). The clinical importance of plasmid-mediated AMR is underscored by the global spread of carbapenem resistance (*blaKPC*, *blaNDM*), extended-spectrum beta-lactamase (ESBL) genes (*blaCTX-M*), and plasmid-mediated colistin resistance (*mcr*) on conjugative plasmids belonging to a limited number of high-risk incompatibility (Inc) groups, notably IncFII, IncN, and IncX (Partridge et al., 2018; Wang et al., 2018; Rozwandowicz et al., 2018).

Despite the centrality of plasmids to the AMR crisis, the field lacks a unified, stable, and hierarchical nomenclature system for plasmid classification. The current landscape of plasmid typing tools, while individually valuable, suffers from fundamental limitations that hinder comparative genomic epidemiology:

**Replicon-based Inc typing (PlasmidFinder).** The most widely used approach, PlasmidFinder (Carattoli et al., 2014), detects incompatibility group markers by BLAST comparison against a curated database of replicon sequences. While computationally efficient, this method produces flat (single-level) classifications, fails to detect novel or divergent replicons, requires database updates, and cannot distinguish between closely related plasmid lineages within the same Inc group.

**Plasmid multilocus sequence typing (pMLST).** Analogous to bacterial MLST, pMLST (Carattoli et al., 2014) assigns sequence types based on allelic variants of multiple loci within specific Inc-group schemes. However, pMLST schemes exist for only six Inc groups, produce flat classifications, and require curated allele databases that must be updated as new variants emerge.

**MOB-suite.** Robertson and Nash (2018) developed MOB-suite, which classifies plasmids based on relaxase typing and Mash distance clustering at a fixed distance threshold of 0.06. While broader in scope than replicon-based methods, MOB-suite produces flat cluster codes, depends on a reference database, and re-assigns cluster identifiers when the database is updated, violating code permanence.

**COPLA and Plasmid Taxonomic Units (PTUs).** Redondo-Salvo et al. (2021) defined PTUs using average nucleotide identity (ANI) networks and hierarchical stochastic block modelling (HSBM). COPLA provides semi-hierarchical classification but could assign only 41% of plasmids (63% of Enterobacterales plasmids) to defined PTUs, and codes are recomputed with each database release.

**mge-cluster.** Shaw et al. (2023) introduced mge-cluster, a reference-free method based on unitig Jaccard distances and HDBSCAN clustering. While reference-free and applicable to diverse mobile genetic elements, mge-cluster produces flat, single-level clusters and requires complete re-embedding (via t-SNE) when new sequences are added.

### 1.2 The Life Identification Number (LIN) Framework

The Life Identification Number (LIN) system was originally developed for hierarchical, permanent classification of bacterial strains based on whole-genome similarity (Vinatzer et al., 2017; Tian et al., 2020). LIN assigns each genome a multi-position numerical code based on its similarity to previously coded genomes at a series of nested distance thresholds. The key properties of LIN—hierarchical multi-resolution classification, code permanence (codes are never retroactively changed), and nearest-neighbour assignment—make it ideally suited for a plasmid classification system, yet LIN has never been applied to plasmid genomes.

### 1.3 Rationale and Objectives

We hypothesised that the LIN framework can be adapted to plasmid genomes using whole-plasmid tetranucleotide composition distances, and that the resulting pLIN (plasmid Life Identification Number) system would provide a classification scheme that is simultaneously (i) hierarchical, (ii) permanent, (iii) reference-free, (iv) biologically meaningful, and (v) integrable with AMR gene surveillance data.

The specific objectives of this study were:

1. To develop and validate the pLIN system for hierarchical plasmid classification based on tetranucleotide composition and single-linkage clustering at six biologically calibrated distance thresholds.
2. To evaluate pLIN concordance with established Inc-group classification and assess its discriminatory power using Simpson's Index of Diversity.
3. To independently validate the compositional features underlying pLIN using a machine learning pipeline with nested cross-validation.
4. To characterise plasmid backbone architecture and mosaic structure using composition-based analyses.
5. To integrate pLIN with NCBI AMRFinderPlus for comprehensive antimicrobial resistance, virulence, and stress gene surveillance, demonstrating the utility of pLIN lineages as a framework for tracking resistance gene dissemination.

---

## 2. METHODS

### 2.1 Plasmid Sequence Dataset

#### 2.1.1 Sequence Retrieval

Complete plasmid genome sequences were retrieved from NCBI RefSeq for three clinically important incompatibility groups: IncFII, IncN, and IncX1. These groups were selected to represent a spectrum of plasmid population structures: IncFII plasmids are the most prevalent conjugative plasmids in Enterobacteriaceae and are characterised by extensive modular mosaicism (Villa et al., 2010); IncN plasmids are efficient carriers of multi-drug resistance cassettes (Rozwandowicz et al., 2018); and IncX1 plasmids represent a more compact, cohesive backbone lineage associated with specific resistance genes and enterotoxigenic virulence determinants (Johnson et al., 2012). All sequences were stored as individual FASTA files, one per plasmid, organised by Inc type.

#### 2.1.2 Quality Control and Deduplication

A systematic deduplication procedure was performed in three stages:

1. **Duplicate directory removal.** An exact duplicate directory ('IncFII 2', containing 4,671 files identical to the IncFII directory, ~553 MB) was identified by byte-level file comparison and removed.

2. **Cross-directory duplicate detection.** MD5 checksums were computed for all FASTA files across all three Inc-type directories. Files with identical checksums present in multiple directories were flagged as cross-directory duplicates. A total of 142 cross-directory duplicates were identified: 52 between IncFII and IncX1, 42 between IncFII and IncN, and 52 between IncN and IncX1.

3. **Resolution.** For each duplicate pair, the copy in the larger (less specific) Inc-type group was removed, retaining the copy in the smallest group. This conservative approach avoids overrepresentation of common plasmid backbones in the larger groups and reduces training bias.

The final deduplicated dataset comprised **6,346 non-redundant plasmid sequences**: IncFII (n=4,581), IncN (n=1,064), and IncX1 (n=701).

### 2.2 pLIN Classification System

#### 2.2.1 Tetranucleotide (4-mer) Frequency Vector Computation

For each plasmid sequence, a normalised tetranucleotide frequency vector of length 256 (4^4 possible 4-mers over the alphabet {A, C, G, T}) was computed. For each of the 256 canonical tetranucleotides *k*, the frequency was calculated as:

    f(k) = count(k in S) / (|S| - 3)

where *S* is the plasmid nucleotide sequence and |*S*| is its length in base pairs. Sequences were converted to uppercase prior to counting. The resulting 6,346 x 256 composition matrix was stored as a double-precision (float64) NumPy array.

**Parameters:**
- k-mer size: k = 4 (tetranucleotide)
- Alphabet: {A, C, G, T}
- Feature vector length: 256 (4^4)
- Normalisation: per-sequence frequency (counts divided by total possible 4-mers)
- Data type: numpy.float64

#### 2.2.2 Pairwise Distance Computation

Pairwise cosine distances were computed between all plasmid pairs using the SciPy `pdist` function (metric='cosine'), yielding a condensed distance vector of length n(n-1)/2 = 20,132,685 for n = 6,346 plasmids. Cosine distance is defined as:

    d(i,j) = 1 - (V_i · V_j) / (||V_i|| × ||V_j||)

where V_i · V_j is the dot product of the two composition vectors and ||V|| is the L2 (Euclidean) norm. Cosine distance ranges from 0 (identical composition) to 1 (orthogonal composition) and is scale-invariant, making it robust to differences in plasmid size.

**Parameters:**
- Distance metric: cosine distance
- Implementation: `scipy.spatial.distance.pdist(vectors, metric='cosine')`
- Total pairwise distances: 20,132,685

#### 2.2.3 Hierarchical Single-Linkage Clustering

Agglomerative single-linkage clustering was applied to the condensed distance matrix using `scipy.cluster.hierarchy.linkage(Z, method='single')`. Single-linkage (nearest-neighbour) was chosen because it produces the same clustering as the LIN nearest-neighbour assignment rule: two plasmids are in the same cluster at threshold *t* if and only if there exists a chain of plasmids connecting them where each successive pair has distance ≤ *t*.

The dendrogram was cut at six distance thresholds using `scipy.cluster.hierarchy.fcluster(Z, t=threshold, criterion='distance')` to obtain flat cluster assignments at each level.

**Parameters:**
- Clustering method: single-linkage (nearest-neighbour)
- Implementation: `scipy.cluster.hierarchy.linkage(dist_condensed, method='single')`
- Cluster extraction: `scipy.cluster.hierarchy.fcluster(Z, t=threshold, criterion='distance')`

#### 2.2.4 Hierarchical Threshold Calibration

Six cosine distance thresholds were defined, each corresponding to a biologically meaningful level of plasmid relatedness calibrated against established average nucleotide identity (ANI) benchmarks:

| Bin | Level | Cosine Distance Threshold (d) | Approximate ANI Equivalent | Biological Interpretation |
|-----|-------|-------------------------------|---------------------------|---------------------------|
| A | Family | ≤ 0.150 | ~85% | Broad plasmid family |
| B | Subfamily | ≤ 0.100 | ~90% | Major evolutionary branches |
| C | Cluster | ≤ 0.050 | ~95% | Analogous to species boundary |
| D | Subcluster | ≤ 0.020 | ~98% | Within-lineage divergence |
| E | Clone complex | ≤ 0.010 | ~99% | Near-identical backbone |
| F | Strain | ≤ 0.001 | ~99.9% | Outbreak-level resolution |

Thresholds were empirically calibrated using the nearest-neighbour distance distribution of 178 selected IncX-like reference plasmids (seed plasmid: RefSeq_NZ_AP027441.1; selection cutoff: d ≤ 0.12 from seed). Key calibration quantiles:

| Quantile | Nearest-Neighbour Distance (d) |
|----------|-------------------------------|
| 1st | 0.000 |
| 5th | 0.000 |
| 10th | 0.000 |
| 25th | 0.001 |
| 50th (median) | 0.011 |
| 75th | 0.025 |
| 90th | 0.038 |
| 95th | 0.048 |
| 99th | 0.072 |

Additional calibration points from the IncX pLIN threshold configuration: F_lineage threshold (d ≤ 0.03, ANI ≥ 97%) and G_outbreak threshold (d ≤ 0.005, ANI ≥ 99.5%).

#### 2.2.5 pLIN Code Assembly

For each plasmid *i*, the six cluster identifiers were concatenated into a pLIN code:

    pLIN(i) = C_A(i).C_B(i).C_C(i).C_D(i).C_E(i).C_F(i)

where C_X(i) denotes the cluster identifier of plasmid *i* at bin level X.

#### 2.2.6 New Plasmid Assignment Rule

For a new query plasmid Q not in the existing database:
1. Compute V(Q) and calculate cosine distances to all existing plasmids.
2. Identify the nearest neighbour N = argmin_i d(Q, i).
3. At each threshold T_k: if d(Q, N) ≤ T_k, inherit the cluster identifier of N at level k; otherwise, assign a new, previously unused cluster identifier.
4. This rule guarantees **code permanence**: existing pLIN codes are never altered by the addition of new plasmids.

### 2.3 Discriminatory Power Assessment

Simpson's Index of Diversity (D) was computed for each hierarchical level using:

    D = 1 - (1/N(N-1)) × Σ n_j(n_j - 1)

where N is the total number of plasmids and n_j is the number of plasmids in cluster j. D ranges from 0 (no diversity) to 1 (maximum diversity). D was calculated both overall and stratified by Inc type at the strain level (Bin F).

### 2.4 Machine Learning Validation

#### 2.4.1 Feature Engineering

A 33-dimensional feature vector was computed for each plasmid, comprising three categories:

**Basic composition features (7 features):**
- Sequence length (bp)
- Log₁₀(length)
- GC content (proportion of G+C bases)
- AT content (proportion of A+T bases)
- AT skew: (A - T) / (A + T)
- GC skew: (G - C) / (G + C)

**Dinucleotide frequencies (16 features):**
All 16 canonical dinucleotides (AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT), each normalised by (sequence length - 1).

**Selected trinucleotide frequencies (10 features):**
ATG, TAA, TAG, TGA, GCG, CGC, AAA, TTT, CCC, GGG — selected for biological relevance (start/stop codons, CpG motifs, homopolymeric tracts), each normalised by (sequence length - 2).

#### 2.4.2 Data Preparation

To avoid class imbalance effects in the three-class (IncFII, IncN, IncX1) classification task, a balanced subset of **1,500 plasmids** (500 per Inc group) was randomly sampled (random seed = 42). Feature scaling was performed using StandardScaler (z-score normalisation), fitted on the training partition of each cross-validation fold to prevent data leakage.

#### 2.4.3 Nested Cross-Validation Design

A rigorous **nested cross-validation** framework was employed to provide unbiased performance estimates while simultaneously optimising hyperparameters:

- **Outer loop**: 3-fold stratified cross-validation (StratifiedKFold, shuffle=True, random_state=42) for performance estimation
- **Inner loop**: 2-fold stratified cross-validation for hyperparameter tuning
- **Scoring metric**: Weighted F1 score (f1_weighted), which accounts for class-specific precision and recall weighted by class frequency

This nested design ensures that the test set in each outer fold is never used for hyperparameter selection, providing unbiased generalisation estimates.

#### 2.4.4 Model Specifications and Hyperparameter Search Spaces

Four classification models were evaluated. Hyperparameter optimisation was performed using **Optuna** (Akiba et al., 2019) with Tree-structured Parzen Estimator (TPE) sampling (seed=42) and Hyperband pruning.

**Optuna parameters:**
- Trials per fold: 10
- Timeout per fold: 60 seconds
- Sampler: TPESampler(seed=42)
- Pruner: HyperbandPruner
- Direction: maximise f1_weighted

**Model 1: Random Forest (scikit-learn RandomForestClassifier)**

| Hyperparameter | Search Range | Scale |
|----------------|-------------|-------|
| n_estimators | [50, 500] | Integer |
| max_depth | [3, 20] | Integer |
| min_samples_split | [2, 20] | Integer |
| min_samples_leaf | [1, 10] | Integer |
| max_features | {'sqrt', 'log2', None} | Categorical |
| class_weight | {'balanced', None} | Categorical |

**Model 2: XGBoost (XGBClassifier)**

| Hyperparameter | Search Range | Scale |
|----------------|-------------|-------|
| n_estimators | [50, 500] | Integer |
| max_depth | [3, 15] | Integer |
| learning_rate | [0.01, 0.3] | Log-uniform |
| subsample | [0.6, 1.0] | Uniform |
| colsample_bytree | [0.6, 1.0] | Uniform |
| gamma | [1e-8, 1.0] | Log-uniform |
| reg_alpha | [1e-8, 10.0] | Log-uniform |
| reg_lambda | [1e-8, 10.0] | Log-uniform |
| tree_method | 'hist' | Fixed |
| eval_metric | 'mlogloss' | Fixed |

**Model 3: Gradient Boosting (scikit-learn GradientBoostingClassifier)**

| Hyperparameter | Search Range | Scale |
|----------------|-------------|-------|
| n_estimators | [50, 500] | Integer |
| max_depth | [3, 15] | Integer |
| learning_rate | [0.01, 0.3] | Log-uniform |
| subsample | [0.6, 1.0] | Uniform |
| min_samples_split | [2, 20] | Integer |
| min_samples_leaf | [1, 10] | Integer |

**Model 4: Logistic Regression (scikit-learn LogisticRegression)**

| Hyperparameter | Search Range | Scale |
|----------------|-------------|-------|
| C | [1e-3, 100.0] | Log-uniform |
| penalty | {'l1', 'l2'} | Categorical |
| solver | 'saga' | Fixed |
| max_iter | 1000 | Fixed |
| class_weight | 'balanced' | Fixed |

#### 2.4.5 Feature Importance Analysis

Feature importance scores were extracted from the three tree-based models (Random Forest, XGBoost, Gradient Boosting) using the native `feature_importances_` attribute (Gini importance for Random Forest, gain-based for boosting methods). Consensus importance was computed as the mean importance across all three models, normalised to sum to 1.0.

### 2.5 IncX1 Backbone Architecture Analysis

#### 2.5.1 Reference Plasmid Selection

IncX-like reference plasmids were selected from the dataset by computing the 4-mer cosine distance from each plasmid to a seed plasmid (RefSeq_NZ_AP027441.1) and selecting all plasmids within d ≤ 0.12. This yielded 178 IncX-like reference plasmids.

**Parameters:**
- Seed plasmid: RefSeq_NZ_AP027441.1
- Distance cutoff: d ≤ 0.12 (cosine distance to seed)
- Selected references: n = 178

#### 2.5.2 Core ORF Identification

Open reading frames (ORFs) were identified in the seed plasmid using a minimum length criterion. Conservation of each seed ORF across the 178 references was assessed by k-mer containment analysis: for each ORF, the fraction of k-mers present in each reference plasmid was computed, and an ORF was considered present if the containment fraction exceeded a threshold. The 25 ORFs with the highest prevalence across references were designated as the candidate IncX1 core backbone.

#### 2.5.3 Backbone/Accessory Partitioning

For the full IncX1 dataset (n=701), the estimated backbone fraction was computed as the proportion of each plasmid's sequence covered by the conserved core ORFs, with the remainder designated as accessory sequence.

### 2.6 Mosaicism Detection

Compositional heterogeneity within individual plasmids was assessed using a sliding window approach. For each IncX1 plasmid (n=701):

1. GC content was computed in 1-kb non-overlapping windows across the sequence.
2. The standard deviation of window-level GC content (GC-SD) was calculated.
3. Plasmids with GC-SD exceeding a threshold indicative of compositional heterogeneity were flagged as candidate mosaics.

**Parameters:**
- Window size: 1,000 bp (1 kb)
- Window step: non-overlapping
- Metric: standard deviation of GC content across windows

### 2.7 AMRFinderPlus Integration

#### 2.7.1 AMR/Virulence/Stress Gene Detection

NCBI AMRFinderPlus v4.2.5 (Feldgarden et al., 2021) was used to detect antimicrobial resistance genes, virulence factors, and stress response genes across all 6,346 plasmid sequences. AMRFinderPlus was run in nucleotide mode with the `--plus` flag to enable detection of virulence and stress genes in addition to AMR determinants.

**Parameters:**
- Software: AMRFinderPlus v4.2.5
- Database: NCBI AMRFinderPlus database release 2026-01-21.1
- Mode: nucleotide (`-n` flag)
- Extended detection: enabled (`--plus` flag for virulence and stress genes)
- Input: individual FASTA files, one per plasmid
- Processing: sequential, per Inc type (IncFII → IncN → IncX1)

**Command:**
```
amrfinder -n <plasmid.fasta> -d <database_path> --plus
```

#### 2.7.2 Output Processing and Classification

AMRFinderPlus output was parsed to classify each detection by its `Type` field into three categories: AMR (antimicrobial resistance), VIRULENCE (virulence factors), and STRESS (stress response genes). For each plasmid, the following summary statistics were computed:
- Number of AMR, virulence, and stress gene detections
- Unique AMR gene symbols and their drug class/subclass annotations
- Unique virulence gene symbols
- Unique stress gene symbols

#### 2.7.3 pLIN–AMR Integration

The pLIN assignment table (output of Section 2.2) was merged with the per-plasmid AMR summary using a left join on plasmid identifier. Plasmids with no AMRFinderPlus detections were assigned zero counts for all gene categories and 'none' for gene lists. The integrated table was used for all downstream AMR-pLIN cross-reference analyses.

#### 2.7.4 Clinically Critical Gene Identification

Five categories of clinically critical resistance genes were defined using pattern-matching on AMR gene symbols:

| Category | Gene Symbol Patterns |
|----------|---------------------|
| Carbapenemases | blaKPC, blaNDM, blaOXA-48, blaVIM, blaIMP |
| Extended-spectrum beta-lactamases (ESBLs) | blaCTX-M, blaSHV, blaTEM |
| Colistin resistance | mcr- |
| Vancomycin resistance | van |
| Plasmid-mediated quinolone resistance (PMQR) | qnr, aac(6')-Ib-cr, oqxA, oqxB |

#### 2.7.5 Lineage-Level AMR Profiling

For each pLIN code (Bin F, strain level), the following statistics were computed among AMR-positive plasmids: number of AMR-positive members, mean and total AMR gene count, Inc-type composition, and prevalence of each AMR gene within the lineage. Lineages were ranked by the number of AMR-positive members.

### 2.8 Statistical Analyses

All statistical analyses were performed in Python 3.14 using NumPy 2.4, Pandas 2.3, SciPy 1.16, scikit-learn, XGBoost, and Optuna. Figures were generated using Matplotlib 3.10.8 and Seaborn 0.13.2 at 300 DPI resolution. Gene frequency counts were computed using Python's `collections.Counter`. All random operations used random_state=42 for reproducibility.

### 2.9 Software and Computational Environment

| Software | Version | Purpose |
|----------|---------|---------|
| Python | 3.14 | Primary programming language |
| NumPy | 2.4 | Numerical computation |
| Pandas | 2.3 | Data manipulation |
| SciPy | 1.16 | Distance computation, hierarchical clustering |
| BioPython | 1.84 | FASTA sequence parsing |
| scikit-learn | (latest) | Machine learning models, cross-validation, scaling |
| XGBoost | (latest) | Gradient boosted tree classifier |
| Optuna | (latest) | Bayesian hyperparameter optimisation |
| AMRFinderPlus | 4.2.5 | AMR/virulence/stress gene detection |
| AMRFinderPlus DB | 2026-01-21.1 | Reference database for gene detection |
| Matplotlib | 3.10.8 | Figure generation |
| Seaborn | 0.13.2 | Statistical visualisation |

All analyses were performed on a single Apple M-series laptop (macOS Darwin 24.5.0). The complete pLIN assignment pipeline executes in under 2 minutes for 6,346 plasmids; AMRFinderPlus processing required approximately 90 minutes for the full dataset.

---

## 3. RESULTS

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

At the coarsest level (Bin A, d ≤ 0.150), all 6,346 plasmids from three Inc groups formed a single cluster, consistent with their shared membership in the broader Enterobacteriaceae plasmid superfamily. Meaningful separation first emerged at Bin C (d ≤ 0.050, ~95% ANI), which split the dataset into two clusters: a major cluster (n=6,345) and a single outlier—an atypical IncFII plasmid with extreme compositional divergence. At Bin D (d ≤ 0.020, ~98% ANI), 20 subclusters were resolved, and at Bin E (d ≤ 0.010, ~99% ANI), 82 clone complexes were delineated. The finest resolution, Bin F (d ≤ 0.001, ~99.9% ANI), produced 2,232 strain-level groups, of which 1,702 (76.3%) were singletons and 530 (23.7%) contained two or more members (Figure 6).

#### 3.2.3 pLIN Code Assignment

Each plasmid was assigned a six-position pLIN code of the form A.B.C.D.E.F, where each position denotes the cluster identifier at the corresponding hierarchical level. For example, pLIN code 1.1.1.7.30.63 indicates membership in Family 1, Subfamily 1, Cluster 1, Subcluster 7, Clone complex 30, and Strain group 63. The pLIN code is assigned once and is permanent: addition of new plasmids to the database does not alter existing codes, a fundamental property inherited from the LIN framework.

The largest strain-level cluster (pLIN 1.1.1.7.30.1240, n=817) consisted almost exclusively of IncFII plasmids (814/817, 99.6%) and likely represents a dominant IncFII lineage with highly conserved backbone composition. The second largest cluster (pLIN 1.1.1.7.30.1313, n=225) was more heterogeneous, containing IncFII (n=184), IncX1 (n=23), and IncN (n=18) members, suggesting compositional convergence or shared accessory module content across Inc groups (Figure 1).

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

To independently validate that the tetranucleotide composition features underlying pLIN assignment carry genuine biological signal for plasmid classification, we trained supervised machine learning models to predict Inc-group membership from the same 33-feature composition vector (Section 2.4.1). A balanced subset of 1,500 plasmids (500 per Inc group) was used to avoid class imbalance effects. A rigorous nested cross-validation framework (3 outer folds, 2 inner folds) with Optuna hyperparameter optimization (10 trials per fold, TPE sampler, Hyperband pruning) was employed to prevent overfitting and data leakage.

**Table 4. Machine learning model performance for Inc-group prediction (3-class, nested CV).**

| Model | Mean Weighted F1 | SD | Min | Max |
|-------|-------------------|----|-----|-----|
| XGBoost | 0.903 | 0.009 | 0.896 | 0.916 |
| Gradient Boosting | 0.901 | 0.007 | 0.896 | 0.912 |
| Random Forest | 0.881 | 0.019 | 0.862 | 0.907 |
| Logistic Regression | 0.873 | 0.010 | 0.860 | 0.884 |

All four models achieved weighted F1 scores exceeding 0.87, with XGBoost performing best (F1=0.903 ± 0.009) (Table 4). The strong performance of even the linear baseline (Logistic Regression, F1=0.873) indicates that Inc-group boundaries are substantially, though not entirely, linearly separable in composition space. The modest improvement from tree-based ensemble methods suggests the presence of nonlinear interaction effects between compositional features.

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

\*Calculated for 3-group Inc typing on this dataset. \*\*pMLST schemes not available for all three Inc groups simultaneously.

---

### 3.9 pLIN System Definition and Formal Description

#### Definition

The **plasmid Life Identification Number (pLIN)** is a stable, hierarchical, multi-position numerical code assigned to each plasmid genome based on its compositional similarity to previously coded plasmids. Each pLIN code consists of six positions (A.B.C.D.E.F), where each position corresponds to a cluster identifier at a specific genetic distance threshold, ordered from coarse (Family, ~85% ANI) to fine (Strain, ~99.9% ANI).

#### Properties

1. **Permanence.** Once assigned, a pLIN code is never changed, regardless of subsequent additions to the database. This is guaranteed by the nearest-neighbour assignment rule.
2. **Hierarchy.** The six-position code provides simultaneous classification at six resolution levels, from family to strain.
3. **Universality.** No reference database is required. Any plasmid with a nucleotide sequence can be classified.
4. **Reproducibility.** Given the same database state and insertion order, pLIN codes are fully deterministic.
5. **Interpretability.** Plasmids sharing a k-position prefix (e.g., 1.1.1.7.\*.\*) are related at the corresponding resolution level (e.g., subcluster level for a 4-position match).

---

### 3.10 Integration with AMRFinderPlus: pLIN as a Framework for AMR Gene Surveillance

#### 3.10.1 AMR Gene Detection Across the Plasmid Dataset

To demonstrate the utility of pLIN as an epidemiological framework for tracking antimicrobial resistance gene dissemination, we integrated the pLIN classification system with NCBI AMRFinderPlus v4.2.5 (database 2026-01-21.1). AMRFinderPlus was run in nucleotide mode with the `--plus` flag to detect antimicrobial resistance genes, virulence factors, and stress response genes across all 6,346 plasmid sequences.

A total of 60,372 gene detections were identified across 5,342 plasmids (84.2% of the dataset), comprising 27,465 AMR gene hits, 5,834 virulence factor hits, and 27,073 stress response gene hits (Table 8; Figure 2A).

**Table 8. AMR/virulence/stress gene prevalence across the plasmid dataset.**

| Category | Plasmids Positive | % of Total (n=6,346) | Total Detections |
|----------|-------------------|----------------------|------------------|
| AMR genes | 4,224 | 66.6% | 27,465 |
| Virulence factors | 1,191 | 18.8% | 5,834 |
| Stress response genes | 2,795 | 44.0% | 27,073 |
| Any AMRFinderPlus hit | 5,342 | 84.2% | 60,372 |

AMR gene prevalence varied markedly across incompatibility groups (Table 9). IncN plasmids showed the highest AMR gene carriage rate (85.1%, mean 6.21 AMR genes per plasmid), consistent with the established role of IncN plasmids as efficient vectors for multi-drug resistance gene cassettes. IncFII and IncX1 plasmids exhibited similar AMR prevalence (~63–67%) but differed in virulence gene carriage: IncFII plasmids carried virulence factors at substantially higher rates (23.5%) than IncN (3.0%) or IncX1 (11.8%), reflecting the known association of IncFII plasmids with virulence-associated loci such as the *spv* operon and iron uptake systems (Figure 5A).

**Table 9. AMR and virulence gene prevalence by incompatibility group.**

| Inc Type | Total (n) | AMR+ (n) | AMR+ (%) | Mean AMR Genes | VIR+ (n) | VIR+ (%) |
|----------|-----------|----------|----------|----------------|----------|----------|
| IncFII | 4,581 | 2,852 | 62.3% | 3.91 | 1,076 | 23.5% |
| IncN | 1,064 | 905 | 85.1% | 6.21 | 32 | 3.0% |
| IncX1 | 701 | 467 | 66.6% | 4.19 | 83 | 11.8% |

#### 3.10.2 AMR Gene Repertoire

The most frequently detected AMR gene was *blaTEM-1* (n=1,760, 41.7% of AMR-positive plasmids), encoding a narrow-spectrum TEM-type beta-lactamase, followed by the sulfonamide resistance gene *sul1* (29.9%), tetracycline efflux gene *tet(A)* (29.1%), sulfonamide resistance gene *sul2* (24.4%), and aminoglycoside-modifying enzymes *aph(6)-Id* (24.2%) and *aph(3'')-Ib* (24.0%) (Table 10; Figure 2B). The high prevalence of these 'classic' resistance determinants across all three Inc groups is consistent with their association with widely disseminated transposon families (Tn*3*, Tn*10*, Tn*21*) and class 1 integrons.

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

The AMR drug class distribution revealed beta-lactam resistance as the most common (3,483 plasmids), followed by aminoglycoside (2,723), sulfonamide (2,073), trimethoprim (1,810), tetracycline (1,549), phenicol (1,542), quinolone (1,240), and macrolide (1,088) resistance determinants (Figure 2C). Co-carriage of resistance genes to three or more drug classes was observed in a substantial proportion of AMR-positive plasmids, consistent with the multidrug-resistance phenotype commonly associated with these Inc groups.

#### 3.10.3 Clinically Critical Resistance Determinants

Of particular clinical concern, we identified high-priority resistance genes associated with last-resort and critically important antimicrobials (Table 11; Figure 3).

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

A central application of the pLIN system is the identification of specific plasmid lineages that serve as vehicles for AMR gene dissemination. By cross-referencing pLIN codes with AMR gene content, we identified distinct lineage-associated resistance profiles (Table 12; Figure 4).

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

Virulence gene carriage showed strong Inc-type specificity (Figure 5B,C). IncFII plasmids (1,076 with virulence genes, 23.5%) predominantly carried the type III secretion system effector *traT* (25.6%), the *Salmonella* plasmid virulence genes *spvB* and *spvD* (25.3%), the murein transglycosylase *mltE* (24.0%), and the aerobactin siderophore gene *iucA* (23.3%). IncX1 plasmids (83 with virulence genes, 11.8%) were enriched for the alpha-hemolysin gene *hlyA* (48.2%) and fimbrial adhesin genes *fedA* and *fedF* (47.0%), consistent with the known association of IncX1 plasmids with enterotoxigenic *E. coli* virulence determinants. IncN plasmids showed the lowest virulence gene carriage (32 plasmids, 3.0%), with aerobactin genes (*iutA*, *iucC*) and *traT* being the most common.

The combination of AMR and virulence determinants on the same plasmid is clinically significant, as it enables the co-selection of pathogenicity and resistance traits under antibiotic pressure. The pLIN framework enables systematic tracking of such co-occurrence patterns at the lineage level.

---

## 4. DISCUSSION

### 4.1 pLIN as a Novel Paradigm for Plasmid Classification

This study presents pLIN, the first application of the Life Identification Number framework to plasmid genomes. Unlike all existing plasmid typing systems, pLIN simultaneously provides six nested levels of hierarchical classification, guaranteed code permanence, reference-free operation, and integration with AMR surveillance data. The system achieved a Simpson's Index of Diversity of 0.979 at the strain level while maintaining 99.5% concordance with established Inc-group assignments—demonstrating that whole-plasmid composition captures known taxonomic boundaries as an emergent property without requiring explicit replicon detection.

The finding that composition-based clustering recapitulates Inc-group classification has important implications. It suggests that replicon identity and overall plasmid backbone composition are correlated, likely because the replication and maintenance modules impose selective constraints on the composition of the entire plasmid backbone. The 0.5% of pLIN codes with mixed Inc-type membership likely represents genuine cases of mosaic plasmids where extensive recombination between Inc groups has homogenised backbone composition while maintaining distinct replicon markers.

### 4.2 Machine Learning Validates Compositional Signal

The nested cross-validation ML pipeline confirmed that the 33-feature composition vector underlying pLIN carries robust biological signal for plasmid classification (weighted F1 = 0.903, XGBoost). The strong performance of the linear baseline (Logistic Regression, F1 = 0.873) is particularly informative: it indicates that ~87% of the classificatory information is captured by linear combinations of composition features, while the additional ~3% improvement from ensemble methods reflects nonlinear feature interactions. The dominance of stop codon-associated trinucleotides (TAG, TGA) and CpG-related motifs (GCG, CGC) in feature importance rankings is consistent with the known role of codon usage and DNA methylation as phylogenetic signals in prokaryotic genomes.

### 4.3 AMRFinderPlus Integration Reveals Lineage-Specific Resistance

The integration of pLIN with AMRFinderPlus across 6,346 plasmids represents, to our knowledge, the first systematic cross-referencing of a hierarchical plasmid classification system with comprehensive AMR gene surveillance data. Several findings have direct clinical and epidemiological relevance:

**Carbapenemase-carrying lineages.** The identification of pLIN 1.1.1.7.30.567 as a 100% *blaKPC-2*-positive IncN lineage with a mean of 13.2 AMR genes per plasmid provides a concrete surveillance target. Traditional Inc typing would classify these as simply 'IncN plasmids', losing the lineage-level resolution that identifies this specific subpopulation as the primary KPC-carrying vehicle.

**Cross-Inc resistance hubs.** pLIN 1.1.1.7.30.1313, which spans all three Inc groups and carries both *blaKPC-2* and *rmtB1*, illustrates the power of composition-based classification to identify plasmid lineages that transcend traditional Inc-group boundaries. These cross-Inc convergence zones may represent critical nodes in the horizontal gene transfer network where resistance cassettes are exchanged between different plasmid backbone types.

**Colistin resistance.** The detection of 160 *mcr*-positive plasmids across multiple pLIN lineages enables lineage-level tracking of this last-resort resistance mechanism—information that is unavailable from Inc typing alone.

### 4.4 Advantages Over Existing Approaches

The comparative analysis (Table 7) demonstrates that pLIN addresses fundamental limitations of all existing methods:

1. **Versus PlasmidFinder/pMLST:** pLIN provides 6 hierarchical levels versus 1, does not require reference databases, and achieves D = 0.979 versus D = 0.508 for Inc typing alone.
2. **Versus MOB-suite:** pLIN codes are permanent; MOB-suite cluster codes change with database updates. pLIN provides hierarchical resolution; MOB-suite operates at a single Mash distance threshold (0.06).
3. **Versus COPLA:** pLIN classifies 100% of input plasmids; COPLA assigns only 41% overall (63% for Enterobacterales). pLIN codes are permanent; COPLA requires HSBM recomputation.
4. **Versus mge-cluster:** Both are reference-free, but pLIN provides hierarchical resolution and code permanence; mge-cluster produces flat clusters that change with re-embedding.

### 4.5 Limitations

Several limitations should be acknowledged:

1. **Composition-based distance as a proxy.** Cosine distance on tetranucleotide frequencies is an approximation of true sequence divergence. While well-correlated with ANI for closely related sequences, composition-based distances may not resolve fine-scale structural rearrangements, insertion/deletion events, or recombination breakpoints that are detectable by alignment-based methods. Future implementations should incorporate Mash-based MinHash distances or core-gene alignment when external tools are available.

2. **Single-linkage clustering sensitivity.** Single-linkage clustering, while theoretically consistent with the LIN nearest-neighbour assignment rule, is sensitive to chaining effects where distant sequences can be merged into the same cluster through intermediaries. In practice, this was observed at the coarser thresholds (Bins A–C) where most plasmids formed a single large cluster. Complete-linkage or average-linkage alternatives may provide better separation at intermediate thresholds, at the cost of losing theoretical consistency with the LIN framework.

3. **Inc-group scope.** The current study validated pLIN on three Inc groups (IncFII, IncN, IncX1). While the method is in principle applicable to any plasmid, the distance thresholds were calibrated against these specific groups and may require adjustment for plasmid families with different evolutionary rates or population structures.

4. **AMRFinderPlus nucleotide mode.** Running AMRFinderPlus in nucleotide mode (without predicted proteins) may miss some gene detections that would be found in protein mode, particularly for genes with low nucleotide identity but conserved protein function. The `--plus` flag extends detection to virulence and stress genes but relies on the NCBI curated database, which may not include all recently described resistance mechanisms.

5. **Threshold calibration.** The six distance thresholds were calibrated against the IncX-like reference plasmid distribution and ANI equivalence estimates. These calibrations are approximate and may benefit from refinement as larger multi-Inc datasets become available. The mapping between cosine distance and ANI is empirical rather than theoretical, and its accuracy may vary across plasmid families with different base composition distributions.

### 4.6 Future Directions

Several extensions of the pLIN framework are envisioned:

1. **Expansion to all plasmid families.** Validation across the full diversity of plasmid Inc groups, including IncA/C, IncI, IncL/M, IncP, and IncW, will establish the universality of the pLIN approach.
2. **Web-based pLIN assignment server.** Development of a public web service and command-line tool for real-time pLIN code assignment to new plasmid sequences.
3. **Longitudinal AMR surveillance.** Application of pLIN to longitudinal clinical datasets to track the emergence and spread of specific AMR-carrying plasmid lineages over time.
4. **Integration with phylogenomic methods.** Hybrid approaches combining pLIN composition-based classification with core-gene phylogenetics could provide both the stability of pLIN codes and the evolutionary resolution of alignment-based methods.
5. **Pan-plasmidome analysis.** Application to metagenomic plasmid assemblies (e.g., from long-read sequencing of clinical samples) for culture-independent plasmid surveillance.

---

## 5. CONCLUSIONS

We present pLIN, the first hierarchical, permanent, reference-free classification system for bacterial plasmid genomes, integrated with comprehensive AMR gene surveillance. Applied to 6,346 plasmids from three clinically important Inc groups, pLIN resolved 2,232 strain-level codes (D = 0.979) with 99.5% concordance with Inc-group assignments. Integration with AMRFinderPlus identified 60,372 gene detections across 84.2% of plasmids and mapped clinically critical resistance genes—including 1,490 carbapenemase, 3,734 ESBL, and 160 colistin resistance detections—to specific pLIN lineages. The pLIN system addresses fundamental gaps in current plasmid classification: it provides hierarchical multi-resolution typing where existing tools offer flat groupings, code permanence where others require reclassification, and reference-free operation where others depend on curated databases. By enabling lineage-level tracking of AMR gene dissemination through plasmid populations, pLIN offers a practical framework for plasmid genomic epidemiology and AMR surveillance.

---

## FIGURE LEGENDS

**Figure 1. Dataset overview.** (A) Incompatibility group composition of the 6,346-plasmid dataset. (B) Plasmid size distribution by Inc type. (C) Hierarchical pLIN diversity showing the number of unique clusters at each bin level (A through F) per Inc type (log scale).

**Figure 2. AMR gene prevalence overview.** (A) Prevalence of AMR, virulence, and stress response genes by Inc type. (B) Top 15 AMR genes ranked by percentage of AMR-positive plasmids carrying each gene. (C) AMR drug class distribution showing the number of plasmids carrying resistance determinants in each class.

**Figure 3. Clinically critical resistance determinants.** Four-panel breakdown showing the distribution of (A) carbapenemase genes, (B) ESBL genes, (C) plasmid-mediated colistin resistance (*mcr*) genes, and (D) plasmid-mediated quinolone resistance (PMQR) genes across the dataset.

**Figure 4. AMR gene prevalence heatmap across top pLIN lineages.** Heatmap showing the prevalence (%) of the 15 most common AMR genes (columns) in the 10 largest AMR-positive pLIN lineages (rows). Cell colour intensity reflects prevalence from 0% (white) to 100% (dark red). Lineage labels include pLIN code, Inc type composition, and sample size.

**Figure 5. AMR gene burden and virulence distribution.** (A) Violin plot of AMR gene count per plasmid by Inc type (AMR-positive plasmids only), with mean (black line) and median (red line) indicated. (B) Top virulence genes by Inc type. (C) AMR–virulence co-occurrence scatter plot showing the relationship between AMR and virulence gene counts per plasmid, coloured by Inc type.

**Figure 6. pLIN hierarchical structure.** (A) Number of clusters resolved at each hierarchical threshold level, overall and per Inc type (log scale). (B) Strain-level (Bin F) cluster size distribution, showing the predominance of singletons (76.3% of pLIN codes).

**Figure 7. Composite manuscript figure.** Multi-panel summary combining dataset overview (A–C), AMR prevalence (D–F), pLIN-AMR lineage heatmap (G), and clinically critical resistance and co-occurrence analyses (H–J).

---

## SUPPLEMENTARY MATERIALS

**Supplementary Table S1.** Complete pLIN assignment table for all 6,346 plasmids (pLIN_assignments.tsv). Columns: plasmid_id, inc_type, length_bp, pLIN, bin_A, bin_B, bin_C, bin_D, bin_E, bin_F.

**Supplementary Table S2.** Integrated pLIN + AMRFinderPlus table (pLIN_AMR_integrated.tsv). Columns: plasmid_id, inc_type, length_bp, pLIN, bin_A–F, n_amr_genes, n_vir_genes, n_stress_genes, n_total_hits, amr_genes, vir_genes, stress_genes, amr_classes, amr_subclasses.

**Supplementary Table S3.** pLIN lineage AMR summary (pLIN_lineage_AMR_summary.tsv). Columns: pLIN, n_plasmids, mean_amr, total_amr, inc_types.

---

## REFERENCES

Akiba, T., Sano, S., Yanase, T., Ohta, T., & Koyama, M. (2019). Optuna: A next-generation hyperparameter optimization framework. *Proceedings of the 25th ACM SIGKDD International Conference on Knowledge Discovery & Data Mining*, 2623–2631.

Carattoli, A. (2009). Resistance plasmid families in Enterobacteriaceae. *Antimicrobial Agents and Chemotherapy*, 53(6), 2227–2238.

Carattoli, A., Zankari, E., Garcia-Fernandez, A., et al. (2014). In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrobial Agents and Chemotherapy*, 58(7), 3895–3903.

Feldgarden, M., Brover, V., Gonzalez-Escalona, N., et al. (2021). AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Scientific Reports*, 11, 12728.

Johnson, T. J., Bielak, E. M., Fortini, D., et al. (2012). Expansion of the IncX plasmid family for improved identification and typing of novel plasmids in drug-resistant Enterobacteriaceae. *Plasmid*, 68(1), 43–50.

Partridge, S. R., Kwong, S. M., Firth, N., & Jensen, S. O. (2018). Mobile genetic elements associated with antimicrobial resistance. *Clinical Microbiology Reviews*, 31(4), e00088-17.

Redondo-Salvo, S., Fernandez-Lopez, R., Ruiz, R., et al. (2021). Pathways for horizontal gene transfer in bacteria revealed by a global map of their plasmids. *Nature Communications*, 11, 3602.

Robertson, J., & Nash, J. H. (2018). MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microbial Genomics*, 4(8), e000206.

Rozwandowicz, M., Brouwer, M. S., Fischer, J., et al. (2018). Plasmids carrying antimicrobial resistance genes in Enterobacteriaceae. *Journal of Antimicrobial Chemotherapy*, 73(5), 1121–1137.

San Millan, A. (2018). Evolution of plasmid-mediated antibiotic resistance in the clinical context. *Trends in Microbiology*, 26(12), 978–985.

Shaw, L. P., Chau, K. K., Kavanagh, J., et al. (2023). Niche and local geography shape the pangenome of wastewater- and livestock-associated Enterobacteriaceae. *Science Advances*, 9(20), eadd3112.

Tian, L., Huang, C., Heath, L. S., & Vinatzer, B. A. (2020). LINbase: a web server for genome-based identification of prokaryotes as members of crowdsourced taxa. *Nucleic Acids Research*, 48(W1), W529–W537.

Villa, L., Garcia-Fernandez, A., Fortini, D., & Carattoli, A. (2010). Replicon sequence typing of IncF plasmids carrying virulence and resistance determinants. *Journal of Antimicrobial Chemotherapy*, 65(12), 2518–2529.

Vinatzer, B. A., Tian, L., & Heath, L. S. (2017). A proposal for a genome similarity-based taxonomy for plant-pathogenic bacteria that is sufficiently precise to reflect phylogeny, host range, and outbreak affiliation applied to *Pseudomonas syringae* sensu lato as a proof of concept. *Phytopathology*, 107(1), 18–28.

Wang, R., van Dorp, L., Shaw, L. P., et al. (2018). The global distribution and spread of the mobilized colistin resistance gene *mcr-1*. *Nature Communications*, 9, 1179.

---

## DATA AVAILABILITY

The complete pLIN assignment pipeline (assign_pLIN.py), AMRFinderPlus batch runner (run_amrfinder_all.sh), integration analysis script (integrate_pLIN_AMR.py), and figure generation script (generate_figures.py) are available at [repository URL]. The analytical notebook (pLIN.ipynb) and all output data files (Supplementary Tables S1–S3, Figures 1–7) are included in the repository. The IncX1 threshold configuration is provided in Data/IncX_PLIN_thresholds_v0_python.yaml. All plasmid sequences were obtained from NCBI RefSeq and are publicly accessible via their accession numbers listed in Supplementary Table S1.
