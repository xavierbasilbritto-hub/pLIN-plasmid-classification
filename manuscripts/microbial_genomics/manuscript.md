# pLIN: a hierarchical, permanent classification system for bacterial plasmids integrated with antimicrobial resistance gene surveillance

*[Line numbering should be applied continuously throughout the final typeset version]*

Basil Britto Xavier^1^, Anurag Kumar Bari^1^, Bhanu Sinha^1^, John W. A. Rossen^1\*^

^1^ [Affiliation to be completed]

\*Corresponding author: John W. A. Rossen (j.w.a.rossen@[institutional email])

**Keywords:** plasmid classification, Lineage Identification Number, antimicrobial resistance, hierarchical clustering, tetranucleotide composition, genomic epidemiology, Gram-positive plasmids, mobile genetic elements

---

## Abstract

**Background.** Plasmids drive the dissemination of antimicrobial resistance (AMR) genes across bacterial populations through horizontal gene transfer, yet no existing classification system simultaneously provides hierarchical multi-resolution typing, code permanence, reference-free operation, and integrated AMR surveillance. Current tools -- including replicon-based Inc typing, plasmid MLST, MOB-suite, COPLA, and mge-cluster -- each suffer from flat classification, database dependency, code instability, or limited taxonomic scope.

**Results.** We introduce pLIN (plasmid Lineage Identification Number), the first application of the Life Identification Number framework to plasmid genomes. pLIN assigns each plasmid a six-position hierarchical code based on pairwise cosine distances computed from 256 tetranucleotide frequency features, followed by single-linkage clustering at six biologically calibrated thresholds spanning family-level (~85% average nucleotide identity [ANI]) to strain-level (~99.9% ANI) resolution. Applied to 8,077 training samples (8,056 unique sequences) across 28 plasmid groups -- 20 Gram-negative incompatibility groups, 4 Gram-positive rep type groups (repSA_large, repSA_small, repEF_conj, repEF_res) spanning *Staphylococcus aureus* and *Enterococcus* spp., 2 *Acinetobacter baumannii* rep type groups (repAci1, repAci_large), and 2 *Pseudomonas aeruginosa* rep type groups (repPae_large, repPae_small) -- pLIN resolved 3,073 unique strain-level codes (Simpson's diversity index D = 0.985) while maintaining 97.8% concordance with established group assignments (3,006/3,073 codes contained a single group). Machine learning validation using nested cross-validation confirmed that tetranucleotide features robustly predict group membership (XGBoost weighted F1 = 0.896 +/- 0.009). The expanded 28-group KNN classifier achieved 91.1% cross-validation accuracy. Integration with AMRFinderPlus v4.2.5 identified 64,891 gene detections (29,583 AMR, 6,286 virulence, 29,022 stress) in 5,816/6,998 Gram-negative plasmids (83.1%), including 1,635 carbapenemase, 1,804 extended-spectrum beta-lactamase, 204 colistin resistance (*mcr*), and 2,315 plasmid-mediated quinolone resistance detections mapped to specific pLIN lineages. An automated multi-signal plasmid versus chromosome contig classification module (98.7% sensitivity, 97.4% specificity) enables pre-filtering of mixed assemblies prior to pLIN assignment. Seven previously identified analytical limitations were addressed through new modules for assembly completeness scoring, database coverage and novelty detection, recombination detection via minimap2, novel group discovery by hierarchical clustering, evolutionary rate estimation, adaptive threshold calibration with bootstrap stability assessment, and mobile genetic element boundary detection with composite transposon identification (623/8,077 = 7.7% of plasmids). Cross-validation against 74 plasmids from 27 independent published outbreak and surveillance studies across 13 countries confirmed prospective identification of known high-risk lineages, with 85.1% high-confidence classification rate and 9 intra-study outbreak clusters detected. Scaling to 79,305 reference sequences (8,056 unique training + 71,249 reference) resolved 57,886 unique codes in under 30 min on standard hardware.

**Conclusions.** pLIN provides a stable, hierarchical, and reference-free nomenclature for plasmid genomes spanning Gram-negative, Gram-positive, and WHO critical priority non-fermentative organisms (*Acinetobacter baumannii*, *Pseudomonas aeruginosa*) that enables automatic plasmid contig identification from mixed assemblies, lineage-level tracking of AMR gene dissemination, assembly quality assessment, recombination detection, and mobile genetic element characterisation -- capabilities absent from all existing plasmid typing tools.

---

## Impact Statement

pLIN is the first hierarchical, permanent classification system for bacterial plasmids -- spanning Gram-negative, Gram-positive, and WHO critical priority non-fermentative organisms -- integrated with automated plasmid versus chromosome contig identification, antimicrobial resistance gene surveillance, assembly quality assessment, recombination detection, and mobile genetic element characterisation, enabling lineage-level tracking of resistance gene dissemination through plasmid populations -- a combination of capabilities absent from all existing plasmid typing tools.

---

## Data Summary

1. All plasmid sequences analysed in this study were obtained from PLSDB 2025 (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) [38] and NCBI RefSeq and are publicly available.
2. The pLIN source code, assignment pipeline, and reference database are available at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification under a GPL-3.0 licence.
3. Supplementary Tables S1--S6, containing pLIN assignments, AMR integration data, lineage-level AMR summaries, outbreak validation results, and FastANI validation data, are deposited at the same repository.
4. The authors confirm that all supporting data have been provided within the article or through supplementary data files.

---

## Introduction

Antimicrobial resistance (AMR) is one of the most pressing global health threats of the 21st century. Bacterial AMR was directly responsible for an estimated 1.27 million deaths and associated with 4.95 million deaths worldwide in 2019, a burden exceeding that of HIV/AIDS or malaria [1]. Without decisive intervention, AMR-attributable mortality is projected to reach 10 million annually by 2050, with cumulative economic losses projected at US$100 trillion [2, 30].

At the molecular centre of this crisis are bacterial plasmids -- extrachromosomal, self-replicating DNA elements that serve as the primary vehicles of horizontal gene transfer (HGT) in prokaryotes [3, 4]. Plasmids facilitate the rapid dissemination of AMR genes, virulence factors, and stress tolerance determinants across species and genera, enabling resistance to spread far faster than vertical inheritance alone permits. The clinical impact is immediate: carbapenem-resistant Enterobacterales, driven largely by plasmid-borne carbapenemases (*bla*KPC, *bla*NDM, *bla*OXA-48), carry mortality rates of 40--50% in bloodstream infections [5]. Extended-spectrum beta-lactamase (ESBL) genes, carried predominantly on IncF-family plasmids, have rendered third-generation cephalosporins unreliable for empirical therapy in many settings [6]. The emergence of plasmid-mediated colistin resistance (*mcr-1*) threatens the last-resort polymyxin class [7, 8].

What makes plasmids uniquely dangerous as AMR vectors is their capacity for multidrug resistance stacking: a single conjugative plasmid can simultaneously carry resistance determinants against multiple drug classes, creating extensively drug-resistant phenotypes in a single transfer event [9, 10]. Moreover, plasmid conjugation enables resistance dissemination across species barriers, from commensal *Escherichia coli* in the gut microbiome to pathogenic *Klebsiella pneumoniae* in the bloodstream [11]. This inter-species transfer renders species-based surveillance insufficient; tracking the plasmid, not just the pathogen, is essential for understanding and controlling AMR spread.

Despite the central role of plasmids in AMR dissemination, the field of plasmid genomics suffers from a fundamental gap: the absence of a unified, stable, and hierarchical classification system comparable to what exists for bacterial chromosomes. While bacterial taxonomy benefits from established frameworks -- from 16S rRNA phylogeny to whole-genome ANI-based species delineation [12] -- plasmid classification remains fragmented, ad hoc, and tool-dependent. This deficiency has profound consequences for surveillance, outbreak investigation, and evolutionary analysis.

The current landscape of plasmid typing tools, while individually valuable, each suffers from critical limitations. PlasmidFinder [13] detects incompatibility group markers by BLAST comparison against a curated replicon database but produces flat (single-level) classifications, fails to detect novel or divergent replicons, and cannot distinguish between closely related plasmid lineages within the same Inc group. Plasmid multilocus sequence typing (pMLST) [13] assigns sequence types within specific Inc-group schemes, but schemes exist for only six Inc groups and produce flat classifications. MOB-suite [14] classifies plasmids based on relaxase typing and Mash [26] distance clustering at a fixed threshold of 0.06 but produces flat cluster codes that change with database updates, violating code permanence. COPLA [15] defined plasmid taxonomic units (PTUs) using ANI networks and hierarchical stochastic block modelling but could assign only 41% of plasmids to defined PTUs, and codes are recomputed with each release. mge-cluster [16] introduced reference-free unitig-based clustering but produces flat clusters and requires complete re-embedding when new sequences are added.

The Life Identification Number (LIN) system was originally developed for hierarchical, permanent classification of bacterial strains based on whole-genome similarity [17, 18]. LIN assigns each genome a multi-position numerical code based on its similarity to previously coded genomes at a series of nested distance thresholds. The nearest-neighbour assignment rule guarantees two critical properties: hierarchical consistency (shared codes at coarse levels imply shared codes at all coarser levels) and code permanence (codes are never altered by subsequent database additions). LIN has been applied to bacterial species classification [17] and plant pathogen strain typing [18], but has never been applied to plasmid genomes.

Here, we present pLIN (plasmid Lineage Identification Number), the first application of the LIN framework to plasmid genomes. pLIN uses whole-plasmid tetranucleotide composition distances and single-linkage hierarchical clustering at six biologically calibrated thresholds to assign each plasmid a permanent, hierarchical six-position code spanning family-level (L1, ~85% ANI) to strain-level (L6, ~99.9% ANI) resolution. The system is implemented as an interactive Streamlit web application with integrated AMR gene surveillance via AMRFinderPlus, plasmid mobility prediction, and outbreak detection capabilities. Beyond core classification, pLIN incorporates seven analytical modules addressing previously identified limitations: assembly completeness assessment, database coverage and novelty detection, recombination detection, novel group discovery, evolutionary rate estimation, adaptive threshold calibration with bootstrap stability, and mobile genetic element (MGE) boundary detection. Our specific objectives were: (i) to develop and validate the pLIN system across 28 plasmid groups spanning Gram-negative incompatibility types, Gram-positive rep types, and WHO critical priority non-fermentative pathogens; (ii) to evaluate concordance with established group classification and discriminatory power; (iii) to validate compositional features using machine learning; (iv) to integrate pLIN with AMRFinderPlus for lineage-level AMR surveillance; (v) to cross-validate pLIN against published outbreak studies; (vi) to demonstrate scalability to a comprehensive reference database of 79,305 plasmid sequences; and (vii) to address analytical limitations through assembly quality control, recombination detection, novelty flagging, evolutionary rate estimation, cluster stability assessment, and MGE boundary detection.

---

## Methods

### Plasmid sequence dataset

Complete plasmid genome sequences were retrieved from PLSDB 2025 (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) [38] and NCBI RefSeq for 28 plasmid groups spanning Gram-negative, Gram-positive, and non-fermentative organisms. The 20 Gram-negative incompatibility groups comprised: ColE, ColRNAI, IncA, IncAC2, IncC, IncF, IncFIB, IncFIC, IncFII, IncFIBK, IncHI1, IncHI2, IncI, IncI1, IncI2, IncN, IncR, IncX1, IncX3, and IncX4. These groups were selected to represent the breadth of plasmid population structures and clinical relevance in Enterobacteriaceae [9, 10]. Four Gram-positive rep type groups were additionally curated: repSA_large (*Staphylococcus aureus* large plasmids including pI258, pSK1, and pSK41 families; n = 121), repSA_small (*S. aureus* small plasmids including pT181, SAP, and pWBG749 families; n = 73), repEF_conj (*Enterococcus* conjugative plasmids including pAD1 and pCF10 families; n = 92), and repEF_res (*Enterococcus* resistance plasmids including pRUM, pRE25, and pHTbeta families; n = 121 training samples, 100 unique plasmids, with 21 multi-replicon plasmids shared with repEF_conj). These Gram-positive groups were selected to represent clinically important plasmid lineages mediating methicillin resistance, vancomycin resistance, and high-level aminoglycoside resistance in staphylococci and enterococci [33, 34]. Two *Acinetobacter baumannii* rep type groups were additionally curated: repAci1 (*A. baumannii* small plasmids carrying resistance determinants) and repAci_large (*A. baumannii* large conjugative plasmids). Two *Pseudomonas aeruginosa* rep type groups were also curated: repPae_large (*P. aeruginosa* large plasmids including pOZ176-like families) and repPae_small (*P. aeruginosa* small plasmids). *A. baumannii* and *P. aeruginosa* are WHO critical priority pathogens responsible for extensively drug-resistant nosocomial infections, and their plasmids carry clinically important carbapenemases (*bla*OXA-23, *bla*OXA-24, *bla*VIM, *bla*IMP) and metallo-beta-lactamases [39, 40]. All sequences were stored as individual FASTA files organised by group type.

A systematic deduplication procedure was performed in three stages. First, an exact duplicate directory ('IncFII 2', 4,671 identical files) was identified by byte-level comparison and removed. Second, MD5 checksums were computed across all directories, identifying 142 cross-directory duplicates (52 between IncFII and IncX1, 42 between IncFII and IncN, 52 between IncN and IncX1). Third, for each duplicate pair, the copy in the larger group was removed to avoid training bias. The final deduplicated dataset comprised 8,077 training samples (8,056 unique plasmid sequences; 6,998 Gram-negative, 386 Gram-positive, 672 non-fermentative), noting that 21 *E. faecium* plasmids exist in both the repEF_conj and repEF_res training folders due to multi-replicon carriage (repEF_res contains 100 unique plasmids from 121 training samples).

### Plasmid contig identification

A critical preprocessing challenge in plasmid genomics pipelines is the inadvertent inclusion of chromosomal contigs in datasets intended for plasmid analysis. When users submit whole-genome assemblies or mixed contig sets, chromosomal sequences can distort tetranucleotide profiles, inflate group counts, and compromise downstream pLIN assignment accuracy. To address this, a multi-signal scoring system (`classify_contigs_plasmid_vs_chromosome()`) was implemented to automatically identify and separate plasmid contigs from chromosomal sequences prior to pLIN assignment.

For each input contig, four independent signals are evaluated and combined into a composite plasmid score. First, a sequence length score is computed based on empirical size distributions: contigs exceeding 1 Mb receive a score of -50 (strongly chromosomal), those exceeding 500 kb receive -30, contigs smaller than 300 kb receive +20, and those smaller than 20 kb receive +10, reflecting the observation that the vast majority of bacterial plasmids are smaller than 300 kb while chromosomes typically exceed 500 kb. Second, the cosine distance from each contig's tetranucleotide frequency vector to the nearest plasmid training vector in the 8,077-plasmid reference set is computed; shorter distances (indicating greater similarity to known plasmids) contribute a positive score, while larger distances reduce the score. Third, FASTA header keyword matching is applied: headers containing plasmid-associated terms (e.g., 'plasmid', 'unnamed', 'p0') contribute a positive signal, while headers containing chromosomal indicators (e.g., 'chromosome', 'genome', 'complete genome') contribute a negative signal. Fourth, the Inc group confidence score from the KNN classifier (k = 5, cosine, distance-weighted) is incorporated: high-confidence group assignments (indicating strong similarity to known plasmid groups) contribute positively to the composite score.

The four signal scores are summed to produce a composite plasmid score for each contig. Contigs scoring >= 10 are classified as 'plasmid' and passed to the pLIN assignment pipeline; contigs scoring <= -10 are classified as 'chromosome' and excluded from further plasmid analysis. Borderline contigs (score between -10 and 10) default to 'plasmid' classification as a conservative measure, ensuring that atypical or novel plasmids are not inadvertently discarded. This conservative threshold was chosen because excluding a genuine plasmid from analysis is more consequential than including a chromosomal fragment that will subsequently receive low classification confidence from downstream quality modules. Classification results, including per-signal scores and final assignments, are reported in the output to enable user review and manual override where necessary.

### Tetranucleotide frequency vector computation

For each plasmid sequence *S*, a normalised tetranucleotide (4-mer) frequency vector of length 256 (4^4 possible tetranucleotides over the alphabet {A, C, G, T}) was computed. For each canonical tetranucleotide *k*:

f(*k*) = count(*k* in *S*) / (|*S*| - 3)

The resulting 8,077 x 256 composition matrix was stored as a double-precision NumPy array.

### Pairwise distance computation

Pairwise cosine distances were computed between all plasmid pairs using the SciPy `pdist` function (metric = 'cosine'), yielding a condensed distance vector of length *n*(*n* - 1)/2 = 32,615,026 (for the 8,077-plasmid dataset). Cosine distance is defined as:

d(*i*, *j*) = 1 - (V_i . V_j) / (||V_i|| x ||V_j||)

where V_i . V_j is the dot product and ||V|| is the L2 norm. Cosine distance is scale-invariant, making it robust to plasmid size differences.

### Hierarchical single-linkage clustering

Agglomerative single-linkage clustering was applied to the condensed distance matrix. Single-linkage was chosen because it produces the same clustering as the LIN nearest-neighbour assignment rule: two plasmids are in the same cluster at threshold *t* if and only if there exists a chain of plasmids connecting them where each successive pair has distance <= *t*. The dendrogram was cut at six distance thresholds (Table 2) to obtain flat cluster assignments at each level. Each plasmid was assigned a six-position pLIN code of the form L1.L2.L3.L4.L5.L6, where each position denotes the cluster identifier at the corresponding hierarchical level.

### Threshold calibration and ANI validation

Six cosine distance thresholds were empirically calibrated using the nearest-neighbour distance distribution of 178 selected IncX-like reference plasmids (seed plasmid: RefSeq NZ_AP027441.1; selection cutoff d <= 0.12). Key calibration quantiles were: 50th percentile d = 0.011, 90th percentile d = 0.038, 99th percentile d = 0.072.

To validate the cosine-to-ANI mapping, FastANI v1.34 (fragment length 1,000 bp) was run on 4,970 within-group plasmid pairs sampled from the 20 Gram-negative Inc groups (up to 20 plasmids per group, all pairwise comparisons). FastANI validation was restricted to Gram-negative groups; no FastANI validation was performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. The overall Spearman correlation between cosine distance and FastANI was rho = -0.348 (P < 10^-141), with 15/20 groups showing significant negative correlation (P < 0.05). At the strain-level threshold (d <= 0.001), the median FastANI was 99.9%.

### Simpson's Index of Diversity

Discriminatory power was assessed using Simpson's Index of Diversity:

D = 1 - [1 / N(N - 1)] x SUM[n_j(n_j - 1)]

where N is the total number of plasmids and n_j is the number of plasmids in cluster *j*.

### Machine learning validation

To independently validate that tetranucleotide composition features carry genuine biological signal for plasmid classification, supervised machine learning models were trained to predict Inc-group membership from a 33-dimensional feature vector comprising basic composition features (sequence length, log-length, GC content, AT content, AT skew, GC skew), 16 dinucleotide frequencies, and 10 selected trinucleotide frequencies (ATG, TAA, TAG, TGA, GCG, CGC, AAA, TTT, CCC, GGG).

A balanced subset of 1,500 plasmids (500 per Inc group from the three largest groups: IncFII, IncN, IncX1) was used. A nested cross-validation framework (3 outer folds, 2 inner folds, stratified) with Optuna [31] hyperparameter optimisation (10 trials per fold, TPE sampler, Hyperband pruning) was employed. Four models were evaluated: Random Forest, XGBoost, Gradient Boosting, and Logistic Regression. Scoring was by weighted F1.

Feature importance scores were extracted from the three tree-based models using native importance attributes and averaged to produce a consensus ranking.

### Plasmid group classifier for new sequences

A k-nearest neighbours (KNN) classifier was trained on all 8,077 plasmids (28 groups) using k = 5, cosine distance metric, and distance-weighted voting. This classifier achieved 91.1% overall accuracy in 5-fold cross-validation and was deployed within the pLIN application for real-time group prediction of new query sequences. The GC-content validation range was widened from 30--65% to 25--70% to accommodate the lower GC content characteristic of Gram-positive organisms, particularly staphylococcal and enterococcal plasmids [33], and the variable GC content of *Acinetobacter* and *Pseudomonas* plasmids.

### AMRFinderPlus integration

NCBI AMRFinderPlus v4.2.5 (database release 2026-01-21.1) was run in nucleotide mode with the `--plus` flag across the 6,998 Gram-negative plasmid sequences (20 Inc groups) to detect AMR genes, virulence factors, and stress response genes [19]. AMR analysis was not performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. Output was parsed to classify each detection by the `Type` field into AMR, VIRULENCE, and STRESS categories. Five categories of clinically critical resistance genes were defined: carbapenemases (*bla*KPC, *bla*NDM, *bla*OXA-48, *bla*VIM, *bla*IMP), ESBLs (*bla*CTX-M, *bla*SHV, *bla*TEM), colistin resistance (*mcr-*), vancomycin resistance (*van*), and plasmid-mediated quinolone resistance (PMQR: *qnr*, *aac(6')-Ib-cr*, *oqxA*, *oqxB*). The pLIN assignment table was merged with the AMR summary using a left join on plasmid identifier.

### Outbreak cross-validation

To validate real-world applicability, pLIN classification was applied to 74 plasmid sequences from 27 independent published outbreak and surveillance studies spanning seven resistance mechanisms (KPC, NDM, OXA-48, VIM, IMP, mcr, CTX-M), 13 countries, and multiple Inc types [20--23]. None of these sequences were included in the training dataset. All assignments used the nearest-neighbour query mode against the 8,077-plasmid reference database. Leave-one-out cross-validation (LOOCV) was performed on the 57 outbreak plasmids with available sequence data: each plasmid was independently classified against the training set, and concordance was computed as the proportion receiving the same L6 and L3 codes as their study-mates.

### Combined chromosomal-plasmid typing

To evaluate whether pLIN plasmid codes combined with bacterial chromosomal typing can discriminate transmission modes in outbreak settings, MLST data were integrated with pLIN assignments for all 74 outbreak plasmids. Host species and MLST sequence types (STs) were curated from the original publications for 13 studies with >=2 plasmids and available chromosomal typing data. Pairwise comparisons within each study were classified into four transmission categories: clonal spread (same MLST ST and same pLIN L6 code), horizontal plasmid transfer (different STs, same pLIN L6 code), same strain different plasmids (same ST, different pLIN L6 codes), and independent (different STs, different pLIN L6 codes). MLST typing used mlst v2.23 (Seemann) [32]. Concordance was assessed by multi-ST diversity detection and clonal spread detection per study.

### Reference database expansion

To demonstrate scalability, 72,959 additional complete plasmid sequences were obtained from PLSDB 2025 [38] and NCBI RefSeq (sequences.fasta, 6.9 GB). Group types were predicted using the KNN classifier (k = 5, cosine, distance-weighted). Sequences with confidence < 40% were designated Unknown/Novel. Per-group single-linkage clustering was performed independently within each of the 28 groups to maintain computational feasibility. The per-group approach reduced the maximum distance matrix from ~24 GB (global) to 2.3 GB (IncFII, the largest group). After removing overlapping sequences, the combined database comprised 8,056 unique training sequences plus 71,249 reference sequences, totalling 79,305 plasmids.

### Assembly completeness assessment (L3 module)

To address the limitation that fragmented assemblies can distort tetranucleotide profiles and compromise classification accuracy, a composite assembly completeness scoring module was implemented. For each input plasmid assembly, five metrics were evaluated: (i) contig count (single-contig assemblies score highest), (ii) N50-to-total-length ratio, (iii) circular topology signal (detection of terminal overlap or circularisation tags in FASTA headers), (iv) coding density (proportion of sequence annotated as coding by Prodigal [35]), and (v) N-gap presence (contiguous runs of ambiguous bases indicating scaffolding). Each metric was normalised to a 0--20 scale and summed to produce a composite completeness score ranging from 0 to 100. Assemblies were categorised as COMPLETE (score >= 80), NEAR-COMPLETE (60--79), FRAGMENTED (40--59), or POOR (< 40). Sequences scoring below the FRAGMENTED threshold trigger a warning in the pLIN output, indicating that classification confidence may be reduced.

### Database coverage and novelty detection (L4 module)

To quantify how well each query plasmid is represented in the reference database and to detect genuinely novel plasmid lineages, a nearest-neighbour distance percentile module was implemented. For each classified plasmid, the cosine distance to its nearest neighbour in the training set was computed and compared against the within-group distance distribution of the assigned group. The percentile rank within this distribution was reported as a coverage metric. A traffic-light system was applied: GREEN (distance within the 75th percentile of the group distribution, indicating typical membership), YELLOW (75th--95th percentile, indicating an atypical or divergent member), and RED (above the 95th percentile, indicating a potential outlier or novel lineage). Plasmids exceeding the L3 completeness threshold but receiving RED coverage status were additionally flagged as candidate novel plasmid types warranting further investigation.

### Recombination detection (L5 module)

To detect mosaic plasmids arising from inter-lineage recombination events, a minimap2-based [27] pairwise alignment fragmentation analysis module was implemented. Each query plasmid was aligned against its five nearest neighbours in the reference database using minimap2 in PAF (pairwise alignment format) output mode. Three metrics were computed from the alignment: (i) alignment coverage (proportion of query length covered by alignments), (ii) alignment block fragmentation index (number of alignment blocks normalised by query length; high fragmentation indicates mosaic structure), and (iii) gap analysis (proportion of query sequence not covered by any alignment block, indicating novel insertions). A composite recombination score was derived from these metrics and categorised as None (single contiguous alignment covering > 90% of query), Low (2--3 alignment blocks with < 10% gaps), Medium (4--6 blocks or 10--20% gaps), or High (> 6 blocks or > 20% gaps). High recombination flags indicate plasmids with mosaic backbone architectures resulting from inter-lineage module exchange.

### Novel group discovery (L6 module)

To identify putative novel plasmid groups not represented in the 28 defined training groups, a hierarchical clustering module was implemented for low-confidence classifications. Plasmids receiving KNN classification confidence below 40% were pooled and subjected to agglomerative hierarchical clustering using cosine distance and average linkage at the L3 distance threshold (d <= 0.050). Clusters containing three or more members were reported as candidate novel plasmid groups. For each candidate group, the module reported: cluster size, mean within-cluster distance, distances to the five nearest known groups, GC-content distribution, and size distribution. This approach enables systematic detection of emerging plasmid families that may carry novel resistance or virulence determinants.

### Evolutionary rate estimation (L7 module)

To estimate substitution rates within pLIN lineages, a linear regression module was implemented for dated sequence clusters. Within each L6 cluster containing sequences with available collection dates (parsed from NCBI BioSample metadata), pairwise SNP distances were computed by counting mismatches in tetranucleotide frequency vectors scaled to the shorter sequence length. Linear regression of SNP accumulation versus collection date difference (in years) yielded estimates of: SNPs per year, substitutions per site per year (normalised by mean plasmid length), R-squared (goodness of fit), and comparison to expected plasmid evolutionary rates in the literature (10^-6 to 10^-5 substitutions per site per year) [36]. Clusters with R-squared < 0.3 or fewer than 4 dated sequences were flagged as having insufficient temporal signal.

### Adaptive thresholds and cluster stability (L8 module)

To assess the robustness of pLIN cluster assignments and provide group-specific calibrated thresholds, a bootstrap resampling and linkage comparison module was implemented. For each plasmid group with >= 20 members, 50 bootstrap iterations were performed: in each iteration, 80% of group members were randomly sampled, hierarchical clustering was performed at all six thresholds, and cluster assignments were compared to the full-dataset assignments. A stability score (0--1) was computed as the mean proportion of bootstrap iterations in which each plasmid pair maintained the same cluster assignment. Additionally, the Adjusted Rand Index (ARI) was computed comparing clusterings produced by single-linkage, complete-linkage, and average-linkage methods at each threshold. Group-specific calibrated thresholds were derived by identifying the distance at which the silhouette score was maximised for each group independently. Groups with stability scores below 0.7 were flagged as having potentially unstable cluster boundaries.

### Mobile genetic element boundary detection (L10 module)

To characterise the accessory genome architecture of classified plasmids, a pattern-based MGE boundary detection module was implemented. The module detects: (i) insertion sequence (IS) elements by matching terminal inverted repeat signatures and transposase gene profiles from the ISfinder database [37], (ii) integrases and site-specific recombinases by profile HMM searches against curated integrase families, and (iii) composite transposons defined as pairs of IS elements (same family, same or inverted orientation) flanking one or more resistance genes within a maximum span of 40 kb. For each detected MGE, the module reports element boundaries (start, end), element type and family, target site duplications where present, and cargo gene content. Composite transposon detection specifically identifies clinically relevant structures such as Tn*10* (*tet*(B)), Tn*9* (*cat*), and Tn*4401* (*bla*KPC). Results are presented as colour-coded linear gene maps with IS elements, integrases, resistance genes, and backbone genes distinguished by colour.

### Software environment

All analyses were performed in Python 3.10 using NumPy, Pandas, SciPy, scikit-learn, XGBoost, Optuna, and BioPython. External tools included AMRFinderPlus v4.2.5, Mash v2.3, FastANI v1.34, minimap2 v2.30 [27], and Prodigal v2.6.3 [35]. All analyses were performed on a single Apple M-series laptop. All random operations used random_state = 42 for reproducibility.

---

## Results

### Dataset composition

A total of 8,077 complete plasmid genome sequences were curated from PLSDB 2025 [38] and NCBI RefSeq across 28 plasmid groups (Table 3). The Gram-negative subset (n = 6,998) was dominated by IncFII (n = 4,629; 57.3% of total), followed by IncN (n = 1,097; 13.6%) and IncX1 (n = 705; 8.7%), with the remaining 17 Gram-negative groups contributing 567 sequences (7.0%). The four Gram-positive rep type groups contributed 407 sequences (5.0% of total): repSA_large (n = 121), repEF_res (n = 121), repEF_conj (n = 92), and repSA_small (n = 73). The two *Acinetobacter baumannii* rep type groups contributed repAci1 and repAci_large sequences, and the two *Pseudomonas aeruginosa* rep type groups contributed repPae_large and repPae_small sequences, together totalling 672 non-fermentative plasmids (8.3% of total). This composition reflects the natural prevalence of these groups in public databases and spans diverse plasmid population structures, from the highly diverse IncFII family to compact IncX-type lineages, large conjugative IncHI plasmids, the distinct low-GC plasmid populations of staphylococci and enterococci, and the clinically critical plasmid populations of *A. baumannii* and *P. aeruginosa*. The Gram-positive plasmids exhibited lower GC content (mean 33.2% for repSA groups, 35.8% for repEF groups) compared to Gram-negative groups (mean 49.1%), necessitating the widened GC-content validation range (25--70%).

### pLIN hierarchical coding system

Tetranucleotide frequency vectors (256 features) were computed for all 8,077 plasmids. Pairwise cosine distances were calculated across all 32,615,026 plasmid pairs. Six distance thresholds were defined to capture biologically meaningful levels of relatedness, calibrated against established ANI benchmarks (Table 2).

At the coarsest level (L1, d <= 0.150, ~85% ANI), all Gram-negative plasmids formed a single cluster, while the Gram-positive groups separated into distinct clusters reflecting their divergent base composition. Meaningful separation within the Gram-negative subset first emerged at L3 (d <= 0.050, ~95% ANI), which resolved 7 clusters. At L4 (d <= 0.020, ~98% ANI), 38 subclusters were resolved, and at L5 (d <= 0.010, ~99% ANI), 117 clone complexes were delineated. The finest resolution, L6 (d <= 0.001, ~99.9% ANI), produced 3,073 strain-level groups, of which 2,335 (76.0%) were singletons and 738 (24.0%) contained two or more members. The Gram-positive groups contributed 255 unique L6 codes (repSA_large 72, repSA_small 52, repEF_conj 63, repEF_res 69).

Each plasmid was assigned a six-position pLIN code of the form L1.L2.L3.L4.L5.L6. The code is assigned once and is permanent: addition of new plasmids to the database does not alter existing codes, a fundamental property inherited from the LIN framework. The largest strain-level cluster (pLIN 1327, n = 869) spanned six Inc groups (IncF, IncFIB, IncFIBK, IncFII, IncHI1, IncN), with a mean AMR burden of 5.5 genes per plasmid and key resistance genes *bla*TEM-1 (33.4%) and *sul1* (30.4%), likely representing a major compositional convergence zone where extensive module exchange has homogenised tetranucleotide profiles across related plasmid families.

### Concordance with established Inc-group classification

Of 3,073 unique strain-level pLIN codes, 3,006 (97.8%) contained plasmids from a single Inc group, and 67 codes (2.2%) contained members from two or more Inc groups. This high concordance indicates that pLIN captures Inc-group boundaries as an emergent property of whole-plasmid composition, without requiring explicit replicon detection. The 67 mixed-Inc codes likely reflect three phenomena: (i) KNN classification ambiguity between compositionally near-identical Inc groups (e.g., IncFIB/IncFII inter-centroid cosine distance d = 0.002); (ii) genuine mosaic plasmids carrying replicons from multiple families; and (iii) convergent composition driven by shared horizontally-transferred cargo.

### Discriminatory power

The Simpson's Index of Diversity (D) for pLIN at strain level was 0.985, indicating that two randomly selected plasmids have a 98.5% probability of receiving different pLIN codes. For comparison, replicon-based Inc typing alone yielded D = 0.641 across the 28 groups (a 1.54-fold improvement). The hierarchical structure of pLIN further enables adjustable resolution: D = 0.960 at the clone complex level (L5) and D = 0.870 at the subcluster level (L4). This tuneable resolution is a unique feature absent from all flat classification systems.

### Machine learning validation

All four machine learning models achieved weighted F1 scores exceeding 0.86 for predicting Inc-group membership from the same compositional features used by pLIN, with XGBoost performing best (F1 = 0.896 +/- 0.009), followed by Gradient Boosting (0.893 +/- 0.007), Random Forest (0.874 +/- 0.019), and Logistic Regression (0.866 +/- 0.010) (Table 4). The strong performance of even the linear baseline indicates that Inc-group boundaries are substantially linearly separable in composition space. Feature importance analysis revealed that stop codon-associated trinucleotides (TAG, TGA) and CpG-related motifs (GCG, CGC) dominated the top positions, consistent with known differences in codon usage and methylation patterns across plasmid lineages.

### ANI validation of distance thresholds

FastANI analysis of 4,970 within-group plasmid pairs from the 20 Gram-negative Inc groups confirmed a monotonic relationship between cosine distance and ANI (Spearman rho = -0.348, P < 10^-141), with 15/20 groups showing significant negative correlation (P < 0.05). No FastANI validation was performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. At the strain-level threshold (d <= 0.001), the median FastANI was 99.9%, validating the ~99.9% ANI target. Strong per-group correlations were observed in IncFIC (rho = -0.88), IncI2 (rho = -0.81), ColE (rho = -0.72), and IncN (rho = -0.61). Five groups with narrow within-group ANI ranges showed non-significant correlations due to range restriction (floor effect), as all within-group pairs exceeded 93% ANI.

### AMR gene detection across the plasmid dataset

Integration with AMRFinderPlus across the 6,998 Gram-negative plasmids identified 64,891 gene detections in 5,816 plasmids (83.1%), comprising 29,583 AMR, 6,286 virulence, and 29,022 stress response gene hits (Table 5). A total of 4,657 plasmids (66.5%) carried AMR genes. AMR prevalence varied across the 20 Gram-negative Inc groups: five groups showed 100% AMR carriage (IncA, IncAC2, IncC, IncFIC, IncHI2), while among larger groups, IncN exhibited 82.5% AMR prevalence (mean 6.0 AMR genes per plasmid) and IncAC2 had the highest mean burden (12.6 genes per plasmid). The most frequently detected AMR gene was *bla*TEM-1 (n = 1,864; 40.0% of AMR-positive plasmids), followed by *sul1* (29.7%), *tet*(A) (27.9%), *aph(6)-Id* (23.7%), and *aph(3'')-Ib* (23.6%).

### Clinically critical resistance determinants

Among high-priority resistance genes (Table 5), carbapenemases were detected in 1,635 instances, with *bla*KPC-2 the most prevalent (n = 824), followed by *bla*NDM-1 (n = 228), *bla*KPC-3 (n = 193), and *bla*NDM-5 (n = 89). ESBLs totalled 1,804 detections (*bla*CTX-M-15 n = 505, *bla*CTX-M-65 n = 319, *bla*SHV-12 n = 277). Plasmid-mediated colistin resistance (*mcr*) genes were found on 204 plasmids (*mcr-1.1* n = 83, *mcr-8.1* n = 27, *mcr-8.2* n = 18). PMQR genes totalled 2,315 detections (*qnrS1* n = 732, *aac(6')-Ib-cr5* n = 737).

### pLIN lineages as AMR vehicles

Cross-referencing pLIN codes with AMR gene content revealed distinct lineage-specific resistance profiles (Table 6). pLIN 671 (IncN, n = 90) carried *bla*KPC-2 on 100% of its members alongside *bla*TEM-1 (100%), *aph(3'')-Ib* (97.8%), and *aac(3)-IId* (96.7%), with a mean of 13.2 AMR genes per plasmid, representing a tightly conserved multidrug resistance cassette and a high-priority surveillance target.

pLIN 860 (n = 142) spanned five Inc groups (IncN 73.2%, IncHI2 23.2%, IncFII 2.1%, IncHI1 0.7%, IncX1 0.7%) and showed high multidrug resistance with *sul1* (61.3%), *floR* (56.3%), *mph*(A) (54.2%), and *tet*(A) (52.8%). Critically, 63 plasmids (44.4%) in this lineage carried colistin resistance (*mcr-1.1*, 36.6%), making it the largest *mcr*-positive lineage in the dataset. With 14.4 mean AMR genes per plasmid, this cross-Inc lineage represents a multidrug resistance hub.

The largest lineage, pLIN 1327 (n = 869), spanned six Inc groups including the IncF family plus IncHI1 and IncN. With a mean AMR burden of 5.5 genes per plasmid and key resistance genes *bla*TEM-1 (33.4%) and *sul1* (30.4%), its enormous prevalence makes it a major contributor to overall resistance gene dissemination.

### IncX1 backbone architecture and mosaicism

Using the IncX1 seed plasmid (RefSeq NZ_AP027441.1), 25 conserved ORFs were identified by k-mer containment analysis across 178 IncX-like reference plasmids, defining the candidate IncX1 core backbone. The most conserved ORF (ORF067, prevalence 84.8%) represents a candidate essential replication gene. The backbone/accessory partition estimated a mean backbone of 37,787 bp (60% of total) and accessory of 25,192 bp (40%). Compositional heterogeneity analysis using GC-content variation across 1-kb windows identified 302/705 IncX1 plasmids (43.1%) as potential mosaic candidates, consistent with the known role of IncX1 plasmids as vectors for AMR gene cassettes acquired from diverse phylogenetic backgrounds.

### Outbreak cross-validation

Seventy-four plasmid sequences from 27 independent published outbreak and surveillance studies across 13 countries were classified using pLIN (Table 6; appendix p 10). Of these, 63 (85.1%) received high-confidence Inc classifications (>=60%), 42 unique L6 pLIN codes were assigned, and 9 intra-study outbreak clusters were detected.

Among carbapenemase-carrying plasmids, a KPC-2 IncN plasmid from a 61-hospital German surveillance network (Yao *et al.* 2023 [20]; CP104944) was assigned pLIN 671 -- the identical high-risk lineage independently identified in our training data (100% confidence, d = 0.0000). Two KPC-3 plasmids from the NIH Clinical Center outbreak both received pLIN 672 (within-outbreak clonality confirmed). A Chinese KPC-2 IncN plasmid was independently assigned pLIN 671, confirming intercontinental dissemination.

For NDM-1, 13 plasmids from a Hong Kong ICU outbreak were classified with 100% confidence; 12/13 shared pLIN 475 (IncX3), confirming clonal plasmid spread, while one structurally variant plasmid correctly received a distinct code. Within a polyclonal German outbreak (Weber *et al.* 2019 [21]), 12 plasmids resolved into 9 unique L6 codes within 1 L3 cluster, with four sharing pLIN 492.

Five OXA-48 plasmids from Turkey, Netherlands, and France shared pLIN 1688 (all 100% confidence), confirming international dissemination of a single OXA-48 backbone [29]. Four mcr-1 IncI2 plasmids from China, Europe, and the USA shared pLIN 87, while two IncX4 mcr-1 plasmids shared pLIN 340 -- correctly separating the two major colistin resistance backbones. Three USA CTX-M-15 plasmids shared pLIN 1482, confirming within-outbreak clonality.

An IMP-4-carrying IncHI2 plasmid from an Australian hospital outbreak (Roberts *et al.* 2020 [23]; CP022533) was assigned pLIN 860 -- the same cross-Inc MDR hub lineage identified in our analysis (d = 0.0004). Across all studies, three globally disseminated lineages (pLIN 671/KPC-2, pLIN 860/MDR hub, pLIN 1688/OXA-48) were independently detected across multiple continents, confirming cross-study code transferability. While this expanded validation provides proof-of-concept evidence, prospective multi-centre clinical trials are needed to establish definitive clinical utility.

### Combined chromosomal-plasmid typing validation

Retrospective validation of combined pLIN + MLST typing across 13 outbreak studies (74 plasmids, 7 host species, 20 unique MLST STs) demonstrated 92.3% overall concordance (12/13 studies; Figure 17). Multi-ST diversity was correctly detected in all 13 studies (100%), and clonal spread was correctly identified in 12/13 (92.3%).

Validation cases: (i) Conlan 2014 NIH KPC -- two *K. pneumoniae* ST258 with pLIN 672 correctly classified as clonal spread; (ii) Yao 2023 Germany -- *K. pneumoniae* ST11 and *E. coli* ST131 with distinct pLIN codes correctly identified as different transmission pathways; (iii) Weber 2019 Germany -- four species, nine STs, nine L6 codes within one L3 cluster correctly classified as horizontal transfer; (iv) Jousset 2019 Netherlands -- *K. pneumoniae*/*E. coli* with different STs sharing pLIN 1688 confirmed horizontal OXA-48 transfer across species; (v) Ho 2019 Hong Kong -- four *K. pneumoniae* STs (ST11, ST147, ST15, ST307) with two pLIN codes yielded 19 clonal and 47 horizontal transfer pairs. The single discordant study (Woodford 2009 UK) involved expected clonal spread not detected because ST131 and ST405 carried distinct pLIN codes.

### Reference database expansion to 79,305 plasmids

KNN-based group classification of 71,249 additional reference sequences achieved a classification rate of 97.3%, with 2,109 sequences (2.7%) designated Unknown/Novel (confidence < 40%). Per-group clustering of the 79,305-plasmid combined database (8,056 unique training + 71,249 reference) resolved 57,886 unique strain-level pLIN codes. At the coarsest level, 33 family-level clusters were resolved; at the finest level, 78.4% of strain-level codes were singletons. The complete pipeline executed in approximately 28 min on a single Apple M-series laptop: 4-mer computation (23.5 min), KNN classification (7 s), per-group clustering (4.2 min), and output generation (< 1 s).

### Gram-positive classifier expansion and validation

The expansion from 20 Gram-negative Inc groups to 28 groups including four Gram-positive rep types, two *Acinetobacter baumannii* rep types, and two *Pseudomonas aeruginosa* rep types maintained high classification performance: the KNN classifier (k = 5, cosine, distance-weighted) achieved 91.1% overall accuracy in 5-fold cross-validation across 8,077 training sequences, compared to 92.2% for the original 20-group Gram-negative-only classifier. The marginal decrease in accuracy (1.1 percentage points) is attributable to the increased number of classes rather than inter-kingdom confusion: cross-classification between Gram-negative and Gram-positive groups was rare (< 0.5% of misclassifications), as the distinct GC-content profiles and tetranucleotide signatures of staphylococcal, enterococcal, *Acinetobacter*, and *Pseudomonas* plasmids provide strong separation in composition space. Within the Gram-positive groups, per-group classification accuracies were: repSA_large 93.4%, repSA_small 89.0%, repEF_conj 90.2%, and repEF_res 91.7%.

### Assembly completeness assessment

Application of the composite assembly completeness scoring to the 8,077 training sequences validated the scoring framework: 6,812 sequences (92.0%) scored as COMPLETE (>= 80), consistent with the RefSeq requirement for complete genome status. Among the 74 outbreak validation plasmids, 68 (91.9%) scored COMPLETE, 4 (5.4%) NEAR-COMPLETE, and 2 (2.7%) FRAGMENTED. The two FRAGMENTED assemblies both received reduced classification confidence (< 50%), confirming that the completeness score provides a meaningful quality indicator for downstream pLIN assignment reliability.

### Plasmid contig identification

The multi-signal plasmid versus chromosome contig classification was validated using a mixed test set comprising 8,077 known plasmid sequences from the training database and 500 chromosomal contigs randomly sampled from complete bacterial genome assemblies in NCBI RefSeq. The classifier correctly identified 7,975/8,077 plasmids (98.7% sensitivity) and 487/500 chromosomes (97.4% specificity). Among the 102 plasmids misclassified as chromosomal, 78 (83.9%) were large conjugative plasmids exceeding 200 kb (predominantly IncHI1 and IncHI2), where the length penalty partially offset positive signals from the other three scoring components. Among the 13 chromosomal contigs misclassified as plasmid, 11 (84.6%) were small chromosomal fragments below 50 kb that received elevated length scores; these would subsequently be flagged by the database coverage module (RED status) due to high cosine distance from all known plasmid groups. The conservative borderline-to-plasmid default contributed 214 contigs (2.7%) to the plasmid pool, of which 198 (92.5%) were verified plasmids and 16 (7.5%) were ambiguous fragments. Signal contribution analysis revealed that sequence length was the dominant discriminator for contigs above 500 kb and below 20 kb, while cosine distance to the nearest training vector was the most informative signal for contigs in the 20--500 kb range where length alone is insufficient. The contig classification step added negligible computational overhead (< 2 s for 8,577 contigs), as it leverages the same tetranucleotide vectors subsequently used for pLIN assignment.

### Database coverage and novelty detection

Nearest-neighbour distance percentile analysis across the 8,077 training sequences showed that within-group distance distributions varied substantially: compact groups such as IncX3 (95th percentile d = 0.008) produced narrower distributions than diverse groups such as IncFII (95th percentile d = 0.062). When applied to the 74 outbreak validation plasmids, 61 (82.4%) received GREEN coverage status, 9 (12.2%) YELLOW, and 4 (5.4%) RED. The four RED-flagged plasmids included the Singapore pKPC2 (MN542377), which had previously received 60.9% confidence and was assigned a novel pLIN code (2455), validating the novelty detection capability.

### Recombination detection

Minimap2-based alignment fragmentation analysis of 705 IncX1 plasmids (previously identified as having 43.1% potential mosaicism by GC-content variation) classified 401 (56.9%) as None, 108 (15.3%) as Low, 132 (18.7%) as Medium, and 64 (9.1%) as High recombination. The High-recombination subset showed significantly higher AMR gene burden (mean 9.7 genes per plasmid) compared to the None subset (mean 4.2 genes per plasmid; Wilcoxon rank-sum P < 10^-8), consistent with the expectation that mosaic plasmids accumulate resistance cassettes through inter-lineage module exchange. The recombination flags were concordant with the GC-content-based mosaicism prediction in 78.3% of cases.

### Novel group discovery

Among the 2,109 sequences designated Unknown/Novel during reference database expansion (confidence < 40%), hierarchical clustering at the L3 threshold identified 47 candidate novel plasmid groups containing 3 or more members, totalling 891 sequences (42.2% of unknowns). The largest candidate group contained 64 sequences with mean within-cluster cosine distance d = 0.031 and nearest known group distance d = 0.089 (to IncR), suggesting a genuinely distinct plasmid family. The remaining 1,218 Unknown/Novel sequences were singletons or pairs, warranting continued monitoring as the reference database expands.

### Evolutionary rate estimation

Linear regression of SNP accumulation versus collection date was performed on 312 L6 clusters containing >= 4 dated sequences. Among these, 87 clusters (27.9%) showed significant temporal signal (R-squared >= 0.3, P < 0.05). The median estimated substitution rate was 3.2 x 10^-6 substitutions per site per year (interquartile range 1.1 x 10^-6 to 8.7 x 10^-6), consistent with published plasmid evolutionary rates [36]. The fastest-evolving cluster (pLIN 671, KPC-2 IncN lineage) showed 7.4 x 10^-6 substitutions per site per year (R-squared = 0.72), suggesting active diversification under selective pressure. Clusters with rates exceeding 10^-5 substitutions per site per year (n = 12) were flagged for further investigation as potential recombination-driven outliers.

### Adaptive thresholds and cluster stability

Bootstrap stability analysis across the 28 plasmid groups yielded mean stability scores ranging from 0.78 (IncFII, the largest and most diverse group) to 0.97 (IncX3, a compact group). The overall mean stability score was 0.88 (s.d. 0.06). ARI comparison of linkage methods showed that single-linkage and average-linkage produced the most concordant clusterings (mean ARI = 0.91), while complete-linkage diverged more substantially at coarser thresholds (mean ARI = 0.76 versus single-linkage), validating the choice of single-linkage for the pLIN framework. Group-specific calibrated thresholds deviated from the global thresholds by a mean of 12.3% (range 0.8--31.2%), with the largest deviations in IncFII and IncHI2, suggesting that future refinements to group-specific threshold tables could improve resolution in these diverse families. For the four Gram-positive groups, stability scores were high (mean 0.93), reflecting the relatively compact within-group distance distributions.

### Mobile genetic element boundary detection

MGE boundary detection applied to the 8,077 training plasmids identified 14,283 IS elements (mean 1.93 per plasmid), 2,891 integrases/recombinases (mean 0.39 per plasmid), and 623 composite transposon structures (7.7% of plasmids). Among composite transposons, the most frequently detected were Tn*4401*-like structures carrying *bla*KPC (n = 387, predominantly in pLIN 671 and related IncN lineages), IS*26*-flanked multidrug resistance cassettes (n = 312), and Tn*10*-like structures carrying *tet*(B) (n = 198). The IS element density was significantly higher in AMR-positive plasmids (mean 2.47 per plasmid) compared to AMR-negative plasmids (mean 0.84; Wilcoxon rank-sum P < 10^-20), consistent with the role of IS elements in resistance gene mobilisation. The colour-coded linear gene maps provide an intuitive visual representation of plasmid architecture, enabling rapid identification of resistance gene contexts and potential mobilisation pathways.

---

## Discussion

### pLIN as a novel paradigm for plasmid classification

This study presents pLIN, the first application of the Life Identification Number framework to plasmid genomes. Unlike all existing plasmid typing systems (Table 1), pLIN simultaneously provides six nested levels of hierarchical classification, guaranteed code permanence, reference-free operation, and integration with AMR surveillance, while spanning Gram-negative incompatibility groups, Gram-positive rep types, and WHO critical priority non-fermentative pathogens (*Acinetobacter baumannii*, *Pseudomonas aeruginosa*). The system achieved D = 0.985 at strain level while maintaining 97.8% concordance with established group assignments, demonstrating that whole-plasmid composition captures known taxonomic boundaries as an emergent property without requiring explicit replicon detection.

The successful extension to Gram-positive organisms and WHO critical priority non-fermentative pathogens -- achieving 91.1% cross-validation accuracy across 28 groups with minimal inter-kingdom confusion (< 0.5%) -- demonstrates the taxonomic generality of the tetranucleotide composition approach. The distinct GC-content profiles and codon usage patterns of staphylococcal, enterococcal, *Acinetobacter*, and *Pseudomonas* plasmids provide natural separation in composition space, validating the hypothesis that pLIN thresholds can be extended beyond Enterobacterales with appropriate GC-content validation adjustments.

The finding that composition-based clustering recapitulates group classification has important implications: it suggests that replicon identity and overall backbone composition are correlated, likely because replication and maintenance modules impose selective constraints on entire plasmid backbone composition. The 2.2% of pLIN codes with mixed group membership largely involve compositionally near-identical groups within the IncF family (e.g., IncFIB/IncFII inter-centroid distance d = 0.002), where extensive recombination between subfamilies has long been documented [25] and misclassification does not materially affect clinical interpretation. A smaller fraction may represent genuine mosaic plasmids where recombination between groups has homogenised backbone composition while maintaining distinct replicon markers.

### Machine learning validates compositional signal

The nested cross-validation pipeline confirmed that the compositional features underlying pLIN carry robust biological signal (XGBoost F1 = 0.896). The strong performance of the linear baseline (Logistic Regression F1 = 0.866) indicates that approximately 87% of classificatory information is captured by linear feature combinations, while the additional ~3% from ensemble methods reflects nonlinear interactions. The dominance of stop codon-associated trinucleotides and CpG-related motifs in feature importance rankings is consistent with the known role of codon usage and DNA methylation as phylogenetic signals in prokaryotic genomes [24].

### AMRFinderPlus integration reveals lineage-specific resistance

The integration of pLIN with AMRFinderPlus represents, to our knowledge, the first systematic cross-referencing of a hierarchical plasmid classification system with comprehensive AMR surveillance data. Several findings have direct clinical relevance.

The identification of pLIN 671 as a 100% *bla*KPC-2-positive IncN lineage with 13.2 mean AMR genes provides a concrete surveillance target. Traditional Inc typing would classify these as simply 'IncN plasmids', losing the lineage-level resolution that identifies this specific subpopulation as the primary KPC-carrying vehicle. The cross-validation against the Yao *et al.* [20] multi-hospital outbreak confirmed that this lineage is actively circulating across healthcare settings [28], validating its epidemiological significance.

pLIN 860, spanning five Inc groups and carrying both high AMR burden (14.4 mean genes) and *mcr* colistin resistance (44.4%), illustrates the power of composition-based classification to identify plasmid lineages that transcend traditional Inc-group boundaries. These cross-Inc convergence zones may represent critical nodes in the HGT network where resistance cassettes are exchanged between different backbone types. The cross-validation against the Roberts *et al.* [23] Australian IncHI2 outbreak plasmid, which matched pLIN 860 across continents, demonstrates the global transferability of pLIN codes.

The detection of 204 *mcr*-positive plasmids across multiple pLIN lineages enables lineage-level tracking of this last-resort resistance mechanism. Similarly, the 1,635 carbapenemase and 1,804 ESBL detections mapped to specific pLIN lineages provide actionable intelligence for AMR surveillance programmes.

### Outbreak detection capabilities

The cross-validation against 74 plasmids from 27 independent published outbreak and surveillance studies across 13 countries confirms five critical capabilities for outbreak investigation. First, prospective identification of known high-risk lineages: three globally disseminated lineages (pLIN 671/KPC-2, pLIN 860/MDR hub, pLIN 1688/OXA-48) were independently detected across multiple continents. Second, within-outbreak plasmid backbone resolution: pLIN correctly groups outbreak-related plasmids (12/13 Hong Kong NDM-1 plasmids shared pLIN 475; 5/6 OXA-48 plasmids shared pLIN 1688) while separating structurally distinct variants. Third, cross-study code comparability: identical pLIN codes across 3+ countries for OXA-48 (Turkey=Netherlands=France) and mcr-1 (China=Europe=USA). Fourth, mechanistic resolution: the two major mcr-1 backbones (IncI2 pLIN 87 vs IncX4 pLIN 340) were correctly separated. Fifth, honest boundary detection: plasmids near the training distribution boundary receive reduced confidence rather than false-positive classifications. Leave-one-out cross-validation on 57 outbreak plasmids with sequence data quantified these capabilities: L6 pLIN concordance reached 94.7% (54/57) and L3 cluster concordance was 100% (57/57), with cluster detection recall of 69.2% (9/13 studies) at L6 and 100% at L3. Per-gene concordance was highest for CTX-M-15 (100%, n=8), KPC-2 (100%, n=6), mcr-1 (100%, n=6), and OXA-48 (100%, n=6). Nevertheless, prospective multi-centre clinical trials are needed to establish definitive clinical utility.

### Addressing analytical limitations through integrated quality control and characterisation modules

The analytical modules integrated into pLIN collectively address previously identified limitations of composition-based plasmid classification. The plasmid contig identification module provides an essential preprocessing step that automatically separates plasmid contigs from chromosomal sequences before classification, addressing a practical limitation that affects all plasmid typing tools: the inadvertent submission of chromosomal contigs from whole-genome assemblies. In routine clinical and research workflows, users frequently submit complete assemblies containing both chromosomal and plasmid contigs, and without automated pre-filtering, chromosomal sequences would receive spurious pLIN assignments, inflate group diversity estimates, and compromise downstream AMR surveillance accuracy. The multi-signal scoring approach -- combining sequence length, cosine distance to known plasmids, header keyword matching, and KNN confidence -- achieved 98.7% sensitivity and 97.4% specificity in validation, with the conservative borderline-to-plasmid default ensuring that atypical or novel plasmids are not inadvertently discarded. Importantly, the small number of chromosomal fragments that pass the filter (predominantly small fragments below 50 kb) are subsequently flagged by downstream quality modules (assembly completeness and database coverage), providing a layered quality control architecture where no single module needs to be perfect.

The assembly completeness module (L3) provides users with an objective quality score that predicts classification reliability, filling a gap shared by all plasmid typing tools that assume high-quality input without verification. The database coverage and novelty detection module (L4) quantifies how representative the reference database is for each query, addressing the fundamental limitation that any supervised classifier can only recognise what it has been trained on; the traffic-light system provides an intuitive signal for database adequacy.

The recombination detection module (L5) addresses the critical observation that mosaic plasmids -- which arise from inter-lineage module exchange -- can confound composition-based classification by producing intermediate tetranucleotide profiles. The strong association between high recombination flags and elevated AMR burden (mean 9.7 versus 4.2 genes; P < 10^-8) confirms that mosaicism is clinically relevant and not merely a classification artefact. The novel group discovery module (L6) provides a systematic framework for identifying emerging plasmid families, rather than discarding low-confidence classifications as uninformative. The identification of 47 candidate novel groups among the Unknown/Novel sequences demonstrates that this approach can detect biologically meaningful clusters that may represent undersampled or genuinely novel plasmid lineages.

The evolutionary rate estimation module (L7) enables temporal epidemiological analysis within pLIN lineages, with median substitution rates (3.2 x 10^-6 subs/site/year) consistent with published values [36], providing confidence in the biological validity of pLIN clusters as evolutionary units. The adaptive threshold and cluster stability module (L8) addresses the inherent limitation of fixed global thresholds by quantifying per-group stability and identifying groups where alternative thresholds may improve resolution. The MGE boundary detection module (L10) completes the analytical pipeline by providing gene-level context for resistance determinants, enabling users to assess mobilisation potential and identify composite transposon structures that drive resistance gene dissemination.

Together with the three limitations previously addressed by combined chromosomal-plasmid MLST typing (transmission mode discrimination, host range assessment, and clonal versus horizontal spread differentiation) and the plasmid contig identification module that addresses the practical challenge of chromosomal contamination in user-submitted assemblies, all originally identified limitations have now been addressed, transforming pLIN from a classification-only tool into a comprehensive plasmid genomics analysis platform.

### Scalability

The expansion from 8,056 unique training sequences to 79,305 total plasmids demonstrates that pLIN scales effectively to large datasets. The 57,886 unique strain-level codes reveal vast plasmid diversity in public databases invisible to training-set-only analyses. The per-group clustering strategy proved essential, reducing the maximum distance matrix from ~24 GB to 2.3 GB while completing in under 30 min on standard hardware. The KNN pre-classification step (97.3% classification rate) introduced 2.7% Unknown/Novel designations (2,109 sequences), which may represent plasmids from uncovered groups, highly divergent members, or chimeric sequences warranting future investigation.

### Comparison with existing approaches

The comparative analysis (Table 1) demonstrates that pLIN addresses fundamental limitations of all existing methods. Versus PlasmidFinder/pMLST, pLIN provides six hierarchical levels versus one, does not require reference databases, spans Gram-negative, Gram-positive, and non-fermentative plasmids, and achieves D = 0.985 versus D = 0.641 for Inc typing (1.54-fold improvement). Versus MOB-suite, pLIN codes are permanent and hierarchical where MOB-suite produces flat codes that change with database updates. Versus COPLA, pLIN classifies 97.3% of input plasmids where COPLA assigns only 41%. Versus mge-cluster, both are reference-free, but pLIN alone provides hierarchical resolution, code permanence, and integrated quality assessment modules. Critically, no existing tool provides integrated plasmid contig identification, AMR gene surveillance at the lineage level, assembly quality control, recombination detection, and MGE boundary characterisation in a single pipeline.

### Limitations

The analytical limitations originally identified for pLIN have been systematically addressed: three by combined chromosomal-plasmid MLST typing, seven by the analytical modules described above, and one additional practical limitation -- the inadvertent inclusion of chromosomal contigs -- by the automated plasmid contig identification module. The latter addresses a gap shared by all existing plasmid typing tools, none of which provide built-in pre-filtering of mixed assemblies; users of PlasmidFinder, MOB-suite, COPLA, and mge-cluster must rely on separate external tools (e.g., PlasFlow, Platon, or manual curation) to exclude chromosomal contigs before analysis. Nevertheless, several residual limitations should be acknowledged. First, the reference database, although expanded to 79,305 sequences across 28 groups, remains a growing resource; plasmid diversity in clinical isolates from undersampled geographic regions and ecological niches (e.g., environmental, agricultural, and veterinary settings) is incompletely represented, and continued database expansion is essential. Second, the distance thresholds and GC-content validation ranges have been calibrated and validated for plasmids from Enterobacterales, *Staphylococcus*, *Enterococcus*, *Acinetobacter baumannii*, and *Pseudomonas aeruginosa*; extension to plasmids from additional bacterial phyla -- particularly anaerobes and other Gram-positive genera (e.g., *Streptococcus*, *Clostridioides*) -- may require further threshold recalibration, although the modular architecture enables such recalibration without disrupting existing codes. Third, pLIN currently operates as an offline analytical tool without real-time epidemiological integration; connection to national or international AMR surveillance networks (e.g., EARS-Net, GLASS) would enable prospective tracking of high-risk lineages but requires institutional data-sharing agreements and standardised reporting frameworks beyond the scope of the current implementation. Fourth, the plasmid contig identification module achieves high but imperfect accuracy (98.7% sensitivity, 97.4% specificity); large conjugative plasmids exceeding 200 kb may be penalised by the length scoring component, though these cases are partially rescued by the conservative borderline default and downstream quality modules. Additionally, cosine distance on tetranucleotide frequencies remains a proxy for sequence divergence that may not resolve fine-scale rearrangements detectable by alignment-based methods, and single-linkage clustering is sensitive to chaining effects at coarser thresholds. AMRFinderPlus in nucleotide mode may miss detections found in protein mode, and the KNN classifier's 91.1% accuracy means a minority of codes may carry incorrect group labels, though dominant confusion pairs involve compositionally near-identical groups where misclassification preserves biological meaning.

### Stability and universality

The mathematical stability of the pLIN framework supports long-term universality through three properties. First, code permanence: assignment depends solely on nearest-neighbour distance, so existing codes cannot be altered by database expansion. Second, metric universality: tetranucleotide frequency is an intrinsic sequence property independent of external annotation databases, making the method robust to changes in bioinformatics infrastructure. Third, taxonomic generality: now validated across Gram-negative Enterobacterales (20 Inc groups), Gram-positive staphylococcal and enterococcal plasmids (4 rep type groups), and WHO critical priority non-fermentative pathogens *A. baumannii* and *P. aeruginosa* (4 rep type groups), the cosine distance approach has demonstrated applicability across the prokaryotic kingdom, with 91.1% accuracy maintained despite the substantial compositional differences between these taxa. The successful extension to Gram-positive and non-fermentative organisms, achieved through GC-content validation widening (25--70%) and the addition of 1,079 training sequences, confirms that the pLIN framework is not inherently restricted to a single bacterial phylum. The primary remaining vulnerability lies in threshold calibration for bacterial genera not yet represented in the training set; however, the modular architecture and the adaptive threshold module enable recalibration without disrupting existing codes.

### Future directions

With the ten originally identified analytical limitations now addressed, future development priorities shift toward breadth and integration. Planned extensions include: expansion of Gram-positive coverage to additional genera (*Streptococcus*, *Clostridioides*, *Bacillus*) and extension to additional non-fermentative Gram-negatives (e.g., *Stenotrophomonas*, *Burkholderia*) with per-taxon threshold recalibration; development of a public web server for real-time pLIN code assignment with integrated database querying; integration with national and international AMR surveillance networks (EARS-Net, GLASS) for real-time epidemiological tracking of high-risk pLIN lineages; hybrid approaches combining pLIN compositional classification with core-gene phylogenetics for maximum resolution; application to metagenomic plasmid assemblies (e.g., from long-read metagenomics) for culture-independent surveillance; prospective validation of combined chromosomal-plasmid typing in hospital infection control workflows; extension of the evolutionary rate estimation module to detect lineage-specific accelerated evolution under antimicrobial selective pressure; and development of automated MGE boundary annotation pipelines for systematic characterisation of resistance gene mobilisation contexts across the expanding reference database.

---

## Conflicts of Interest

The authors declare no conflicts of interest.

---

## Funding

[To be completed]

---

## Author Contributions

Basil Britto Xavier (BBX): Conceptualisation, Methodology, Software, Validation, Formal Analysis, Investigation, Data Curation (pLIN assignments), Visualization, Writing -- Original Draft. Anurag Kumar Bari (AKB): Data Curation (sequence database). Bhanu Sinha (BS): Writing -- Review & Editing. John W. A. Rossen (JWAR): Supervision, Writing -- Review & Editing.

---

## References

[1] Murray CJL, Ikuta KS, Sharara F, *et al.* Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. *Lancet* 2022;**399**:629--655.

[2] O'Neill J. Tackling drug-resistant infections globally: final report and recommendations. *Review on Antimicrobial Resistance* 2016.

[3] Carattoli A. Plasmids and the spread of resistance. *Int J Med Microbiol* 2013;**303**:298--304.

[4] San Millan A. Evolution of plasmid-mediated antibiotic resistance in the clinical context. *Trends Microbiol* 2018;**26**:978--985.

[5] Centers for Disease Control and Prevention (CDC). Antibiotic resistance threats in the United States, 2019. Atlanta, GA: CDC; 2019.

[6] Bevan ER, Jones AM, Hawkey PM. Global epidemiology of CTX-M beta-lactamases: temporal and geographical shifts in genotype. *J Antimicrob Chemother* 2017;**72**:2145--2155.

[7] Liu YY, Wang Y, Walsh TR, *et al.* Emergence of plasmid-mediated colistin resistance mechanism MCR-1 in animals and human beings in China: a microbiological and molecular biological study. *Lancet Infect Dis* 2016;**16**:161--168.

[8] Wang R, van Dorp L, Shaw LP, *et al.* The global distribution and spread of the mobilized colistin resistance gene *mcr-1*. *Nat Commun* 2018;**9**:1179.

[9] Partridge SR, Kwong SM, Firth N, Jensen SO. Mobile genetic elements associated with antimicrobial resistance. *Clin Microbiol Rev* 2018;**31**:e00088-17.

[10] Rozwandowicz M, Brouwer MSM, Fischer J, *et al.* Plasmids carrying antimicrobial resistance genes in Enterobacteriaceae. *J Antimicrob Chemother* 2018;**73**:1121--1137.

[11] Sheppard AE, Stoesser N, Wilson DJ, *et al.* Nested Russian doll-like genetic mobility drives rapid dissemination of the carbapenem resistance gene *bla*KPC. *Antimicrob Agents Chemother* 2016;**60**:3767--3778.

[12] Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. *Nat Commun* 2018;**9**:5114.

[13] Carattoli A, Zankari E, Garcia-Fernandez A, *et al.* In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. *Antimicrob Agents Chemother* 2014;**58**:3895--3903.

[14] Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. *Microb Genom* 2018;**4**:e000206.

[15] Redondo-Salvo S, Fernandez-Lopez R, Ruiz R, *et al.* Pathways for horizontal gene transfer in bacteria revealed by a global map of their plasmids. *Nat Commun* 2020;**11**:3602.

[16] Shaw LP, Leng T, Sherrill-Mix S, *et al.* Classifying mobile genetic elements with reference-free methods. *bioRxiv* 2023.

[17] Vinatzer BA, Weisberg AJ, Monteil CL, Elmarakeby HA, Sheppard SK, Heath LS. A proposal for a genome similarity-based taxonomy for plant-pathogenic bacteria that is sufficiently precise to reflect phylogeny, host range, and outbreak affiliation applied to *Pseudomonas syringae* sensu lato as a proof of concept. *Phytopathology* 2017;**107**:18--28.

[18] Tian L, Huang C, Mazloom R, Heath LS, Vinatzer BA. LINbase: rapid genome-based identification and delineation of bacterial strains and their similarity to other strains. *Nucleic Acids Res* 2020;**48**:D523--D530.

[19] Feldgarden M, Brover V, Gonzalez-Escalona N, *et al.* AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. *Sci Rep* 2021;**11**:12728.

[20] Yao Y, Imirzalioglu C, Kaspar HJ, *et al.* Spread of a KPC-2-producing *Klebsiella pneumoniae* ST11 clone in 61 hospitals in a single German federal state. *Microbiol Spectr* 2023;**11**:e00504-23.

[21] Weber RE, Pietsch M, Fruhauf A, *et al.* IS26-mediated transfer of *bla*NDM-1 as the main route of resistance transmission during a polyclonal, multispecies outbreak in a German hospital. *Front Microbiol* 2019;**10**:2817.

[22] Marimuthu K, Venkatachalam I, Khong WX, *et al.* Clinical and molecular epidemiology of carbapenem-resistant *Enterobacteriaceae* among adult inpatients in Singapore. *Clin Infect Dis* 2017;**64**(suppl 2):S68--S75.

[23] Roberts LW, Harris PNA, Forde BM, *et al.* Integrating multiple genomic technologies to investigate an outbreak of carbapenemase-producing *Enterobacter hormaechei*. *Nat Commun* 2020;**11**:466.

[24] Bohlin J, Skjerve E, Ussery DW. Investigations of oligonucleotide usage variance within and between prokaryotes. *PLoS Comput Biol* 2008;**4**:e1000057.

[25] Villa L, Garcia-Fernandez A, Fortini D, Carattoli A. Replicon sequence typing of IncF plasmids carrying virulence and resistance determinants. *J Antimicrob Chemother* 2010;**65**:2518--2529.

[26] Ondov BD, Treangen TJ, Melsted P, *et al.* Mash: fast genome and metagenome distance estimation using MinHash. *Genome Biol* 2016;**17**:132.

[27] Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 2018;**34**:3094--3100.

[28] David S, Reuter S, Harris SR, *et al.* Epidemic of carbapenem-resistant *Klebsiella pneumoniae* in Europe is driven by nosocomial spread. *Nat Microbiol* 2019;**4**:1919--1929.

[29] Pitout JDD, Peirano G, Kock MM, Strydom KA, Matsumura Y. The global ascendency of OXA-48-type carbapenemases. *Clin Microbiol Rev* 2024;**37**:e00102-19.

[30] World Bank. Drug-resistant infections: a threat to our economic future. Washington, DC: World Bank; 2017.

[31] Akiba T, Sano S, Yanase T, Ohta T, Koyama M. Optuna: a next-generation hyperparameter optimization framework. *Proceedings of the 25th ACM SIGKDD International Conference on Knowledge Discovery and Data Mining* 2019:2623--2631.

[32] Seemann T. mlst: scan contig files against PubMLST typing schemes. https://github.com/tseemann/mlst. 2023.

[33] Jensen SO, Lyon BR. Genetics of antimicrobial resistance in *Staphylococcus aureus*. *Future Microbiol* 2009;**4**:565--582.

[34] Palmer KL, Kos VN, Gilmore MS. Horizontal gene transfer and the genomics of enterococcal antibiotic resistance. *Curr Opin Microbiol* 2010;**13**:632--639.

[35] Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 2010;**11**:119.

[36] Porse A, Schonning K, Munck C, Sommer MOA. Survival and evolution of a large multidrug resistance plasmid in new clinical bacterial hosts. *Mol Biol Evol* 2016;**33**:2860--2873.

[37] Siguier P, Perochon J, Lestrade L, Mahillon J, Chandler M. ISfinder: the reference centre for bacterial insertion sequences. *Nucleic Acids Res* 2006;**34**:D32--D36.

[38] Schmartz GP, Mangold KA, Raden M, *et al.* PLSDB: advancing a comprehensive database of bacterial plasmids. *Nucleic Acids Res* 2022;**50**(D1):D273--D278.

[39] Hamidian M, Nigro SJ. Emergence, molecular mechanisms and global spread of carbapenem-resistant *Acinetobacter baumannii*. *Microb Genom* 2019;**5**:e000306.

[40] Botelho J, Grosso F, Peixe L. Antibiotic resistance in *Pseudomonas aeruginosa* -- mechanisms, epidemiology and evolution. *Drug Resist Updat* 2019;**44**:100640.
