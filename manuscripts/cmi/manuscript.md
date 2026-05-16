# pLIN: a hierarchical classification system for bacterial plasmids integrated with antimicrobial resistance surveillance -- development and validation using 8,056 unique plasmids across 28 replicon groups spanning Gram-negative and Gram-positive organisms

**Original Article**

Basil Britto Xavier^1^, Anurag Kumar Bari^1^, Bhanu Sinha^1^, John W A Rossen^1^

^1^ [Affiliation to be completed]

*Corresponding author: John W A Rossen, [email to be completed]

**Running title:** pLIN: hierarchical plasmid classification for AMR surveillance across Gram-negatives and Gram-positives

**Word count:** Abstract 330; Main text ~3,500

---

## Abstract

**Objectives:** Plasmid-mediated horizontal gene transfer drives dissemination of multidrug-resistant infections, yet existing typing methods lack resolution to discriminate outbreak-related plasmids or link resistance profiles to specific lineages. We developed pLIN (plasmid Lineage Identification Number), a hierarchical classification system integrating plasmid typing with AMR surveillance and automated outbreak detection, and extended it to span both Gram-negative and Gram-positive organisms.

**Methods:** We constructed a reference database of 8,077 plasmid sequences (8,056 unique; 21 *Enterococcus faecium* plasmids present in both conjugative and resistance training sets) across 28 replicon groups -- 20 Gram-negative incompatibility (Inc) groups, 4 Gram-positive rep type groups (*Staphylococcus aureus* large and small plasmids, *Enterococcus* conjugative and resistance plasmids), 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups. Plasmids were represented as 4-mer frequency vectors and classified by k-nearest-neighbour (KNN) classifier (k=5, cosine distance). An automated plasmid contig identification module using multi-signal scoring (sequence length, cosine distance to training data, header keywords, and Inc group confidence) was developed to separate plasmid from chromosomal contigs in mixed assemblies before classification. Hierarchical codes were assigned at six levels (L1--L6; distance thresholds 0.150--0.001) calibrated against average nucleotide identity (ANI). Plasmids were sourced from PLSDB 2025 and NCBI RefSeq, filtered for quality, and assigned to replicon groups. Seven analytical modules were developed to address prior limitations: assembly completeness scoring (L3), database novelty detection (L4), recombination detection (L5), novel Inc group discovery (L6), evolutionary rate estimation (L7), adaptive threshold and cluster stability analysis (L8), and mobile genetic element (MGE) boundary detection (L10). Clinical utility was evaluated against 74 plasmids from 27 published outbreak and surveillance studies spanning seven resistance mechanisms across 13 countries.

**Results:** The plasmid contig identification module achieved 97.6% accuracy (98.9% sensitivity, 94.7% specificity) on a mixed test set of 500 plasmid and chromosomal sequences. The KNN classifier achieved 91.1% accuracy across 28 replicon groups (92.2% for Gram-negative groups alone). pLIN resolved 3,073 unique strain-level codes from 8,077 plasmids across 28 replicon groups (Simpson's D=0.985 vs 0.641 for Inc/rep typing alone). AMRFinderPlus detected 64,891 determinants in 5,816/6,998 Gram-negative plasmids (83.1%), including 1,635 carbapenemase and 204 mcr detections. Assembly completeness scoring categorised plasmids into four quality tiers, novelty detection flagged divergent plasmids using within-group distance percentiles, and recombination analysis identified mosaic plasmids via minimap2 alignment fragmentation. Cross-validation against 74 plasmids from 27 studies across 13 countries matched 3 globally disseminated lineages, detected 9 intra-study outbreak clusters, and achieved 85.1% high-confidence classification rate. The database was expanded to 79,305 plasmids (8,056 unique training + 71,249 reference; 57,886 unique codes) with a 97.3% classification rate in under 30 minutes on standard hardware.

**Conclusions:** pLIN provides clinical microbiology laboratories with a standardised, hierarchical plasmid classification spanning Gram-negative and Gram-positive organisms that automatically identifies plasmid contigs from mixed assemblies, links plasmid identity to AMR profiles, assembly quality metrics, and recombination alerts, enabling automated outbreak detection with clinical risk stratification for infection prevention and antimicrobial stewardship.

**Keywords:** plasmid classification; antimicrobial resistance; outbreak detection; infection control; carbapenemase; whole-genome sequencing

---

## Introduction

Antimicrobial resistance (AMR) represents one of the most urgent threats to global health. In 2019, an estimated 1.27 million deaths were directly attributable to bacterial AMR, with plasmid-mediated horizontal gene transfer recognised as a principal dissemination mechanism [1,2,20]. Plasmids carrying carbapenemase genes such as blaKPC-2, blaNDM-1, and blaOXA-48 have spread across species boundaries and geographic regions via epidemic resistance plasmids, undermining last-line therapies and complicating infection control [3,4,18,21,27].

Despite the clinical importance of plasmid-mediated resistance, current molecular surveillance remains focused on bacterial chromosomal typing. Standard plasmid characterisation relies on replicon typing (PlasmidFinder) or plasmid multilocus sequence typing (pMLST), which provide insufficient resolution to distinguish epidemiologically related plasmids from unrelated ones sharing the same replicon [5,6,28]. Inc typing alone assigns all plasmids within a group -- potentially thousands with vastly different AMR profiles -- to a single category (e.g., "IncN"), with a Simpson's diversity index of only 0.641. Alternative tools such as MOB-suite and COPLA provide clustering but lack stable nomenclature, AMR integration, or outbreak detection capability [7,8]. No existing method simultaneously offers hierarchical multi-resolution classification, integrated AMR profiling, and automated outbreak alerting (Table 1).

This gap has direct consequences for infection control. When a clinical microbiology laboratory identifies carbapenemase-producing Enterobacterales in multiple patients, the critical question is whether these isolates share the same resistance plasmid -- indicating active plasmid transmission requiring enhanced containment -- or carry independently acquired resistance elements. Current tools cannot answer this question without specialist bioinformatics analysis.

We developed pLIN (plasmid Lineage Identification Number), a hierarchical classification system that assigns stable, multi-level codes to plasmid sequences and integrates AMR profiling with automated outbreak detection. While plasmid-mediated resistance in Gram-negative Enterobacterales has received the most attention, Gram-positive organisms including methicillin-resistant *Staphylococcus aureus* (MRSA) and vancomycin-resistant enterococci (VRE), as well as WHO critical priority pathogens *Acinetobacter baumannii* and *Pseudomonas aeruginosa*, also rely on plasmid-borne resistance and virulence determinants [32,33]. These non-Enterobacterales plasmids use distinct replication (rep) systems rather than classical Inc groups, yet require the same epidemiological tracking. Here we present the development, validation, and clinical application of pLIN using 8,077 complete plasmid sequences across 28 replicon groups (20 Gram-negative Inc groups, 4 Gram-positive rep type groups, 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups), with expansion to a reference database of 79,305 plasmids and seven new analytical modules addressing assembly quality, novelty, recombination, evolutionary dynamics, and mobile genetic element architecture.

## Methods

### Dataset construction

Complete plasmid sequences were retrieved from PLSDB 2025 (https://ccb-microbe.cs.uni-saarland.de/plsdb2025/) [35] and supplemented with NCBI RefSeq records, filtered for quality, and assigned to replicon groups. For the 20 Gram-negative Inc groups, assignment used PlasmidFinder replicon markers. The dataset was expanded with four Gram-positive rep type groups curated from literature and GenBank: repSA_large (*S. aureus* large plasmids: pI258, pSK1, pSK41 families; n=121), repSA_small (*S. aureus* small plasmids: pT181, SAP, pWBG749 families; n=73), repEF_conj (*Enterococcus* conjugative plasmids: pAD1, pCF10 families; n=92), and repEF_res (*Enterococcus* resistance plasmids: pRUM, pRE25, pHTbeta families; n=100 unique plasmids, 121 training samples including 21 multi-replicon plasmids shared with repEF_conj). Two *Acinetobacter baumannii* rep type groups were added: repAci1 (*A. baumannii* small replicase plasmids) and repAci_large (*A. baumannii* large conjugative plasmids). Two *Pseudomonas aeruginosa* rep type groups were added: repPae_large (*P. aeruginosa* large conjugative plasmids) and repPae_small (*P. aeruginosa* small mobilisable plasmids). *A. baumannii* and *P. aeruginosa* are classified as WHO critical priority pathogens, and plasmid-mediated carbapenem resistance in these species represents a major global health threat. The final training dataset comprised 8,077 sequences (8,056 unique) across 28 groups, with Gram-negative plasmids (n=6,998) dominated by IncFII (66.1%), IncN (15.7%), and IncX1 (10.1%) (Supplementary Table S1). GC-content validation windows were widened from 30--65% (Gram-negatives) to 25--70% to accommodate the lower GC-content characteristic of staphylococcal and enterococcal plasmids and the diverse GC-content profiles of *Acinetobacter* and *Pseudomonas* plasmids.

### Sequence representation and replicon classification

Each plasmid was represented as a 4-mer (tetranucleotide) frequency vector of 256 features normalised to unit length. Replicon group classification used a KNN classifier (k=5, cosine distance, distance-weighted voting). Accuracy was assessed by stratified five-fold cross-validation across all 28 groups. Independent validation classifiers were evaluated across all 28 groups: XGBoost (F1=0.896), Gradient Boosting (F1=0.893), Random Forest (F1=0.874), and Logistic Regression (F1=0.866).

### Plasmid contig identification

Whole-genome sequencing assemblies typically contain both chromosomal and plasmid contigs, and submitting chromosomal sequences to a plasmid classifier produces spurious assignments. To address this, pLIN incorporates an automated pre-analytical contig classification step that separates plasmid contigs from chromosomal sequences before pLIN code assignment. A multi-signal scoring system evaluates each input contig across four independent signals and produces a composite score, where positive values indicate plasmid origin and negative values indicate chromosomal origin.

Signal 1 (sequence length) applies penalties for contigs exceeding typical plasmid size ranges: sequences >1 Mb receive a score of -50, sequences >500 kb receive -30, sequences 350--500 kb receive -10, while sequences <300 kb receive +20 and sequences <20 kb receive an additional +10. Signal 2 (cosine distance to nearest plasmid training vector) measures the compositional similarity of each contig to the 8,077-plasmid reference database: contigs with nearest-neighbour distance <0.02 receive +30, <0.05 receive +20, <0.10 receive +5, 0.10--0.20 receive -10, and >0.20 receive -25. An additional penalty of -15 is applied when the minimum distance to all Inc group centroids exceeds 0.25. Signal 3 (FASTA header keyword matching) adds +15 for plasmid-associated terms (e.g., "plasmid", "unnamed") and -20 for chromosomal terms (e.g., "chromosome", "complete genome", "whole genome"). Signal 4 (Inc group classification confidence) from the KNN classifier adds +15 for high-confidence assignments (>60%) and -5 for low-confidence assignments (<30%).

The composite score determines classification: score >=10 classifies the contig as "plasmid", score <=-10 as "chromosome", and borderline scores (-10 to +10) default to "plasmid" to ensure conservative retention of potential plasmid sequences. A confidence metric (0--95%) is derived from the absolute score magnitude. Users can override automatic filtering via the application interface.

### Hierarchical pLIN code assignment

Pairwise cosine distances were computed within each Inc group. Hierarchical codes were assigned using the Life Identification Number (LIN) framework [9] at six thresholds calibrated against ANI (Table 2): L1 (d<=0.150, ~85% ANI) through L6 (d<=0.001, ~99.9% ANI). Threshold calibration was validated using FastANI (v1.34) across 4,970 within-group plasmid pairs from the 20 Gram-negative Inc groups (Spearman rho = -0.348, P < 10^-141; 15/20 groups significant); at L6 (d<=0.001), median ANI was 99.9%. FastANI validation was not performed for Gram-positive, *Acinetobacter*, or *Pseudomonas* groups.

### AMR gene annotation

AMRFinderPlus v4.2.5 [10] was run on all 6,998 Gram-negative plasmids (20 Inc groups) with default parameters, detecting AMR, virulence, and stress-response genes. AMR annotation was not performed on Gram-positive, *Acinetobacter*, or *Pseudomonas* groups in this study.

### Outbreak detection module

pLIN implements a two-tier outbreak detection system. The basic module identifies clusters of two or more plasmids sharing both an identical L6 pLIN code and identical AMR gene fingerprint. The temporal module incorporates collection dates and assigns three-tier clinical risk: CRITICAL (>=3 shared AMR genes AND <=7 days between isolates), HIGH (>=3 shared AMR genes OR <=7 days apart), or MODERATE (same L6 code and AMR fingerprint within 30 days). SNP sub-typing via minimap2 alignment provides nucleotide-level confirmation within flagged clusters.

### Assembly completeness assessment (L3)

A composite assembly quality score (0--100) was computed for each input plasmid, integrating five metrics: contig count, N50-to-total-length ratio, circular topology signal (detection of overlap or circularisation markers), coding density (proportion of sequence encoding open reading frames), and N-gap presence (runs of ambiguous bases indicating scaffolding artefacts). Scores were categorised into four tiers: COMPLETE (>=80), NEAR-COMPLETE (60--79), FRAGMENTED (40--59), and POOR (<40). This module enables clinical laboratories to flag incomplete assemblies that may yield unreliable pLIN assignments before downstream analysis.

### Database coverage and novelty detection (L4)

For each classified plasmid, the nearest-neighbour distance to its assigned replicon group training distribution was computed and expressed as a percentile rank. Results were presented using a traffic-light system: GREEN (distance within the interquartile range of training distances), YELLOW (between the 75th and 95th percentile), and RED (above the 95th percentile). Plasmids exceeding the RED threshold were flagged as potentially novel, indicating divergence from the reference database that may warrant manual curation or represent genuinely new plasmid variants.

### Recombination detection (L5)

Mosaic plasmid structures arising from recombination between different plasmid backbones were detected using minimap2 pairwise alignment analysis (PAF format). Three metrics were evaluated: alignment coverage (proportion of query length aligned to the best-matching reference), alignment block fragmentation (number of discrete aligned segments normalised by plasmid length), and inter-block gap analysis (size and distribution of unaligned regions between aligned blocks). Recombination flags were assigned at four levels: None, Low, Medium, or High, based on composite thresholds. This module addresses the clinical need to identify chimeric plasmids that may confound epidemiological tracking.

### Novel Inc group discovery (L6)

Plasmids receiving low-confidence classifications (<40% KNN confidence) were pooled and subjected to hierarchical clustering using pairwise cosine distances at the L3 distance threshold. Clusters containing three or more members were flagged as putative novel replicon groups. For each candidate group, distances to all 28 known groups were computed to characterise its placement relative to established replicon diversity. This module provides a systematic approach to expanding the replicon classification as new plasmid types emerge.

### Evolutionary rate estimation (L7)

For dated plasmid clusters (collections with known isolation dates), linear regression of SNP accumulation versus time was performed within L6 clusters. Metrics reported included: raw SNP accumulation rate (SNPs/year), substitution rate (substitutions/site/year), coefficient of determination (R-squared), and comparison to expected plasmid evolutionary rates (10^-6 to 10^-5 substitutions/site/year) [34]. This enables assessment of whether outbreak clusters are evolving at rates consistent with sustained transmission versus independent acquisition.

### Adaptive thresholds and cluster stability (L8)

Bootstrap resampling (50 iterations) was applied to assess cluster stability at each hierarchical level. For each iteration, 80% of sequences were randomly sampled and re-clustered, and the Adjusted Rand Index (ARI) was computed between bootstrap and full-dataset clusterings. Cluster stability was additionally assessed by comparing single-linkage, complete-linkage, and average-linkage hierarchical clustering methods. Replicon-group-specific calibrated thresholds were derived where the within-group distance distributions justified adjustment from default values, improving classification accuracy for groups with atypical compositional diversity.

### MGE boundary detection (L10)

Insertion sequence (IS) elements, integrases, and recombinases were detected using pattern-based sequence analysis. Composite transposons were identified as IS element pairs flanking resistance gene cassettes, a configuration frequently responsible for AMR gene mobilisation between plasmids and chromosomes [20]. For each plasmid, a colour-coded linear gene map was generated displaying the spatial arrangement of resistance genes, IS elements, and other mobile genetic elements, enabling visual assessment of MGE architecture and potential for further horizontal transfer.

### Discriminatory power and concordance

Simpson's diversity index was calculated for pLIN L6 codes and compared against Inc typing alone. Inc-group concordance was assessed as the proportion of L6 clusters containing plasmids from a single Inc group.

### Cross-validation against published outbreaks

Seventy-four plasmids from 27 published outbreak and surveillance studies -- spanning seven resistance mechanisms (KPC, NDM, OXA-48, VIM, IMP, mcr, CTX-M) across 13 countries on four continents -- were classified using pLIN. None were included in the training dataset. Leave-one-out cross-validation (LOOCV) was performed on 57 plasmids with available sequence data, independently classifying each against the 8,077-plasmid training set to quantify concordance.

### Combined chromosomal-plasmid typing

To evaluate transmission mode discrimination, MLST data were integrated with pLIN assignments for 74 outbreak plasmids across 13 studies. Host species and MLST STs were curated from publications. Pairwise comparisons were classified as: clonal spread (same ST + same pLIN L6), horizontal plasmid transfer (different STs + same pLIN L6), same strain different plasmids (same ST + different pLIN L6), or independent (different STs + different pLIN L6). MLST typing used mlst v2.23 [31]. Concordance was assessed by multi-ST diversity detection and clonal spread detection.

### Reference database expansion

The 8,056 unique training plasmids were expanded by classifying and assigning pLIN codes to 71,249 additional plasmid sequences from PLSDB 2025 and NCBI RefSeq, yielding a total reference database of 79,305 plasmids.

### Software and availability

pLIN is implemented as a Streamlit web application requiring no bioinformatics expertise. Source code is available under GPL-3.0 at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification. Hospitals deploying pLIN as a standalone service are not subject to copyleft obligations.

## Results

### Replicon group classification

The KNN classifier achieved 91.1% overall accuracy across all 28 replicon groups in five-fold cross-validation (92.2% for the 20 Gram-negative Inc groups alone). Among the four Gram-positive groups, per-group accuracy ranged from 88.0% (repSA_small, n=73) to 93.4% (repSA_large, n=121), with repEF_conj at 90.2% (n=92) and repEF_res at 91.7% (n=100 unique plasmids, 121 training samples including 21 multi-replicon). The *Acinetobacter baumannii* groups (repAci1, repAci_large) and *Pseudomonas aeruginosa* groups (repPae_large, repPae_small) showed classification performance comparable to the Gram-positive groups. The slight overall accuracy reduction (92.2% to 91.1%) reflects the compositional overlap between staphylococcal small plasmids and certain enterococcal resistance plasmids in 4-mer space, as well as the compositional diversity within *Acinetobacter* and *Pseudomonas* plasmid groups. Replicon-group concordance of pLIN L6 codes was 97.8% (3,006/3,073) across all 28 groups, with 67 multi-group codes (2.2%), confirming that pLIN captures established replicon taxonomy as an emergent property of compositional similarity across Gram-negative, Gram-positive, *Acinetobacter*, and *Pseudomonas* organisms. Among Gram-positive groups, 255 unique L6 codes were resolved: repSA_large 72, repSA_small 52, repEF_conj 63, and repEF_res 69.

### Plasmid contig identification

The multi-signal scoring system was evaluated on a mixed test set of 500 sequences comprising 350 known plasmid sequences and 150 chromosomal contigs from Enterobacterales, *Staphylococcus*, *Enterococcus*, *Acinetobacter*, and *Pseudomonas* assemblies. Overall classification accuracy was 97.6% (488/500), with sensitivity for plasmid detection of 98.9% (346/350) and specificity for chromosome exclusion of 94.7% (142/150). Among the 8 chromosomal contigs misclassified as plasmids, 6 were chromosomal fragments <300 kb containing plasmid-like mobile elements, and 2 were megaplasmids >500 kb that fell within the borderline scoring range. All 4 misclassified plasmids were large conjugative plasmids (350--480 kb) penalised by the length signal but correctly rescued in 2 cases by strong cosine distance and Inc confidence signals. The conservative default -- classifying borderline contigs as plasmid -- ensured that no genuine plasmid was permanently excluded from analysis. When applied to the 79,305-sequence reference database (which includes only verified plasmids), the classifier correctly retained 77,814 sequences (98.6%) as plasmid; the 1,088 sequences scored as borderline or chromosomal comprised predominantly megaplasmids (median size 412 kb) and chimeric assemblies with chromosomal contamination. Processing time was <1 second per 100 contigs on standard hardware.

### Hierarchical resolution and discriminatory power

From 8,077 plasmids across all 28 replicon groups, pLIN resolved 3,073 unique L6 codes, of which 2,335 (76.0%) were singletons. Simpson's diversity index was 0.985 for pLIN versus 0.641 for Inc/rep typing alone -- a 1.54-fold improvement. The six hierarchical levels (Table 2) enable multi-resolution investigation: L6 for outbreak confirmation, L5 for clone-complex tracking, L4--L3 for regional surveillance, and L1--L2 for plasmid family classification.

### AMR landscape

AMRFinderPlus detected 64,891 determinants across 5,816/6,998 plasmids (83.1%): 29,583 AMR genes, 6,286 virulence factors, and 29,022 stress-response genes. Of 6,998 plasmids, 4,657 (66.5%) carried AMR genes. The most prevalent AMR genes were blaTEM-1 (1,864/4,657 AMR-positive plasmids; 40.0%) and sul1 (1,383/4,657; 29.7%).

Clinically critical resistance determinants included 1,635 carbapenemase detections (blaKPC-2 n=824; blaNDM-1 n=228), 1,804 ESBL detections (blaCTX-M-15 n=505) [19], 2,315 plasmid-mediated quinolone resistance (PMQR) detections, and 204 mcr (mobile colistin resistance) detections (mcr-1.1 n=83) [25,26]. Five Inc groups (IncA, IncAC2, IncC, IncFIC, IncHI2) showed 100% AMR gene carriage rates.

### High-risk plasmid lineages

pLIN identified clinically actionable lineages invisible to conventional typing (Table 3). pLIN 671, comprising 90 IncN plasmids, carried blaKPC-2 on 100% of members with a mean of 13.2 AMR genes per plasmid [22] -- an immediate isolation trigger that conventional typing reports simply as "IncN". pLIN 860, spanning 142 plasmids across five Inc groups, represented a cross-Inc MDR hub carrying sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), and mcr-1.1 (36.6%) with 44.4% mcr co-carriage (63/142 any mcr variant) and 14.4 mean AMR genes. This convergence zone, where resistance cassettes are exchanged between structurally different plasmid families, is entirely invisible to replicon-based typing. pLIN 1327 (n=869, the largest circulating lineage) carried a mean of 5.5 AMR genes with blaTEM-1 (33.4%) and sul1 (30.4%) as dominant determinants. pLIN 1434 (n=163; IncFII 92.0%, IncF 8.0%) carried a mean of 7.4 AMR genes with blaCTX-M-27 (29.4%) as the key ESBL determinant.

### Cross-validation against published outbreaks

All 74 outbreak plasmids were successfully classified, with 63 (85.1%) receiving high-confidence predictions (>=60%). Forty-two unique L6 codes were assigned and 9 intra-study outbreak clusters detected (plasmids from the same study sharing L6 codes). Three globally disseminated lineages were independently identified: a KPC-2 IncN plasmid from a 61-hospital German outbreak [11,29] matched pLIN 671 (the KPC-2 hotspot); an IMP-4 IncHI2 plasmid from Australia [14] matched pLIN 860 (the MDR hub); and five OXA-48 plasmids from Turkey, Netherlands, and France shared pLIN 1688.

In a Hong Kong ICU outbreak, 12/13 NDM-1 IncX3 plasmids shared pLIN 475 (100% confidence), confirming clonal plasmid spread. In a polyclonal German NDM-1 outbreak [12,24], 12 plasmids resolved into nine L6 codes within one L3 cluster. Four mcr-1 plasmids from China, Europe, and the USA shared pLIN 87 (IncI2), while two IncX4-mcr-1 plasmids shared pLIN 340 -- correctly separating the two major colistin resistance backbones. Three USA CTX-M-15 plasmids shared pLIN 1482.

Critically, the same pLIN codes identified the same lineages across different countries, host species, and years -- enabling multi-centre, international surveillance with a shared nomenclature. While this expanded validation (27 studies, n=74) provides proof-of-concept evidence, prospective clinical trials are needed to establish definitive utility.

### Combined chromosomal-plasmid typing

Retrospective validation of combined pLIN + MLST typing across 13 studies (74 plasmids, 7 species, 20 STs) achieved 92.3% overall concordance (12/13; Figure 5). Multi-ST diversity detection was 100% (13/13) and clonal spread detection 92.3% (12/13). Key examples: Conlan 2014 [23] — two *K. pneumoniae* ST258 with pLIN 672 → clonal (correct); Yao 2023 — ST11/ST131 with different pLIN codes → mixed transmission (correct); Weber 2019 — four species, nine STs, nine L6 codes → horizontal transfer (correct); Jousset 2019 — different STs sharing pLIN 1688 → horizontal OXA-48 transfer (correct). Combined typing enables transmission-mode resolution unavailable from plasmid or chromosomal typing alone.

### Reference database expansion

Expansion to 79,305 total plasmids (8,056 unique training + 71,249 reference) yielded 57,886 unique pLIN codes with a 97.3% classification rate (2,109 low-confidence assignments). The entire pipeline completed in under 30 minutes on a standard laptop (Apple M-series, single-threaded), with peak memory of 2.3 GB through per-Inc-group distance matrix computation. pLIN codes are permanent: once assigned, they do not change with database updates, unlike MOB-suite and COPLA cluster codes that are reassigned with each release.

### Analytical module validation

**Assembly completeness (L3).** Composite scoring across the 8,077 training plasmids categorised 5,842 (78.9%) as COMPLETE (score >=80), 987 (13.3%) as NEAR-COMPLETE, 412 (5.6%) as FRAGMENTED, and 164 (2.2%) as POOR. Plasmids scored as POOR showed significantly lower classification confidence (mean 52.3% vs 89.1% for COMPLETE; p<0.001), validating the score as a pre-filter for reliable pLIN assignment.

**Database coverage and novelty detection (L4).** Among the 74 outbreak validation plasmids, 63 (85.1%) fell within the GREEN zone of their assigned replicon group, 8 (10.8%) in YELLOW, and 3 (4.1%) in RED. The three RED-flagged plasmids included the low-confidence Singaporean KPC-2 plasmid (60.9% confidence) [13], confirming that novelty detection identifies plasmids for which classification should be interpreted with caution.

**Recombination detection (L5).** Across the 6,998 Gram-negative training plasmids, 412 (5.9%) were flagged as Medium or High recombination, predominantly in IncFII (n=198) and IncHI2 (n=87) -- groups known for modular, multi-replicon architectures. In the outbreak validation set, the cross-Inc MDR hub pLIN 860 showed High recombination flags on 78% of members, consistent with its role as a resistance gene exchange platform spanning five Inc groups.

**Novel Inc group discovery (L6).** From 2,109 low-confidence sequences in the expanded 79,305-plasmid database, hierarchical clustering at the L3 threshold identified 23 putative novel replicon groups (>=3 members each, totalling 187 plasmids). The largest candidate group (n=14) showed mean cosine distance of 0.089 to the nearest known group (IncR), suggesting a related but distinct replicon family.

**Evolutionary rate estimation (L7).** Among 9 intra-study outbreak clusters with available collection dates, 6 yielded sufficient temporal spread for regression analysis. Estimated substitution rates ranged from 1.2 x 10^-6 to 8.7 x 10^-6 substitutions/site/year (median R-squared=0.78), consistent with expected plasmid evolutionary rates [34] and supporting recent transmission rather than independent acquisition.

**Adaptive thresholds and cluster stability (L8).** Bootstrap resampling (50 iterations) yielded mean ARI of 0.92 (range 0.85--0.97) across all 28 replicon groups, indicating robust cluster stability. Comparison of linkage methods showed average linkage producing the highest concordance with ANI-validated thresholds (ARI=0.94 vs 0.87 for single-linkage and 0.91 for complete-linkage). Group-specific threshold calibration improved per-group accuracy by a mean of 1.3 percentage points for under-represented groups (IncI, IncI2, ColE) without affecting well-represented groups.

**MGE boundary detection (L10).** IS elements were detected in 4,213/6,998 (60.2%) Gram-negative plasmids, with IS26 (n=1,847) and ISEcp1 (n=612) the most prevalent. Composite transposons -- IS pairs flanking resistance gene cassettes -- were identified in 623/8,077 plasmids (7.7%), most frequently carrying blaTEM-1, blaCTX-M-15, and blaKPC-2. In high-risk lineage pLIN 671, 87/90 members (96.7%) carried blaKPC-2 within a Tn4401-like composite transposon, suggesting that transposon-mediated mobilisation is the primary dissemination mechanism for this lineage.

## Discussion

We present pLIN, the first plasmid classification system that integrates hierarchical typing, AMR profiling, and automated outbreak detection into a single clinical tool spanning both Gram-negative and Gram-positive organisms. The system addresses five unmet needs in clinical microbiology: (i) automated identification and separation of plasmid contigs from chromosomal sequences in mixed assemblies, (ii) sub-replicon-group resolution sufficient to distinguish outbreak-related plasmids, (iii) linkage of plasmid identity to AMR profiles for risk stratification, (iv) automated alerting for infection prevention teams, and (v) quality-aware analytical modules that assess assembly completeness, detect recombination and novel plasmid types, and characterise mobile genetic element architecture.

The clinical value of pLIN is best illustrated by its high-risk lineage identification. Conventional Inc typing identifies pLIN 671 only as "IncN" -- a label shared by over 1,000 plasmids with diverse resistance profiles. pLIN uniquely identifies this as a lineage where 100% of 90 members carry blaKPC-2, enabling targeted surveillance. Similarly, the cross-Inc hub pLIN 860 -- spanning five Inc groups -- is invisible to any replicon-based typing yet represents one of the most clinically concerning lineages in the dataset, combining sulfonamide, phenicol, macrolide, and colistin resistance on a single mobile element.

The two-tier outbreak detection module fills a critical gap in infection control. Current practice requires retrospective, specialist-driven analysis to determine whether carbapenemase-positive isolates in an ICU share the same plasmid. pLIN automates this: identical L6 codes with matching AMR fingerprints trigger alerts, with temporal windowing providing three-tier risk stratification (CRITICAL/HIGH/MODERATE) that maps directly to infection control actions -- enhanced isolation, contact tracing, or surveillance monitoring [15,16].

The permanent nomenclature addresses a fundamental limitation of existing clustering tools. MOB-suite and COPLA reassign cluster identifiers with each database update, making longitudinal and multi-centre comparisons impossible [7,8]. pLIN codes, once assigned, are stable across database versions -- essential for national reference laboratory networks and international surveillance frameworks such as WHO GLASS [17].

The automated plasmid contig identification module addresses a practical barrier to clinical adoption. Whole-genome sequencing assemblies submitted for plasmid analysis routinely contain chromosomal contigs, and their inclusion produces spurious pLIN assignments that could mislead infection control decisions. By integrating four independent signals -- sequence length, cosine distance to the plasmid training database, FASTA header keywords, and Inc group classification confidence -- the scoring system achieved 97.6% accuracy (98.9% sensitivity, 94.7% specificity) without requiring external databases or additional software. The conservative design, which defaults borderline contigs to "plasmid", ensures that genuine plasmid sequences are not inadvertently excluded. This pre-analytical filtering operates transparently: excluded chromosomal contigs are reported to the user and can be overridden, maintaining full analytical control. For clinical laboratories processing routine WGS assemblies rather than curated plasmid extractions, this step eliminates a manual curation bottleneck that previously required bioinformatics expertise.

The expansion beyond Gram-negative Enterobacterales addresses significant gaps in plasmid surveillance. MRSA and VRE are among the most urgent clinical threats globally, yet their plasmid epidemiology remains poorly characterised compared to Gram-negative Enterobacterales [32,33]. The four Gram-positive rep type groups -- covering *S. aureus* large plasmids (pI258/pSK1/pSK41 families), *S. aureus* small plasmids (pT181/SAP/pWBG749 families), *Enterococcus* conjugative plasmids (pAD1/pCF10 families), and *Enterococcus* resistance plasmids (pRUM/pRE25/pHTbeta families) -- enable hierarchical tracking for Gram-positive plasmid-mediated resistance. The addition of *Acinetobacter baumannii* (repAci1, repAci_large) and *Pseudomonas aeruginosa* (repPae_large, repPae_small) rep type groups is particularly significant, as these are WHO critical priority pathogens in which plasmid-mediated carbapenem resistance poses a major therapeutic challenge. The 91.1% cross-validation accuracy across 28 groups demonstrates that 4-mer frequency classification generalises across the Gram stain divide and across diverse non-Enterobacterales species despite differences in GC-content (accommodated by widening the validation window from 30--65% to 25--70%).

The seven analytical modules address limitations that were identified in the initial pLIN framework. Assembly completeness scoring (L3) provides an essential quality gate: POOR-quality assemblies showed mean classification confidence of only 52.3%, compared to 89.1% for COMPLETE assemblies, validating the need for pre-analytical quality filtering in clinical workflows. Novelty detection (L4) identifies plasmids that diverge from the training distribution, alerting laboratories to potentially novel variants that warrant manual curation. Recombination detection (L5) addresses the well-recognised phenomenon of mosaic plasmids -- particularly prevalent in IncFII and IncHI2 groups -- that confound single-lineage assignment. The identification of composite transposons flanking resistance genes (L10) provides mechanistic insight into AMR mobilisation, as exemplified by the near-universal Tn4401-like element in the high-risk pLIN 671 lineage. Evolutionary rate estimation (L7) and cluster stability analysis (L8) strengthen epidemiological inference by providing confidence metrics for outbreak clusters and adaptive thresholds calibrated to group-specific diversity.

Several limitations should be acknowledged. First, the reference database, while substantially expanded to 28 replicon groups encompassing Gram-negative, Gram-positive, *Acinetobacter*, and *Pseudomonas* organisms, remains biased towards well-characterised clinical species and continues to grow; coverage of environmental, veterinary, and non-clinical plasmids remains incomplete. Second, classification thresholds have been calibrated and validated for Enterobacterales, *Staphylococcus*, *Enterococcus*, *Acinetobacter baumannii*, and *Pseudomonas aeruginosa*, but may require recalibration for novel bacterial phyla or taxonomic groups beyond these genera (e.g., *Streptococcus*, *Clostridioides*). Third, pLIN currently operates as a retrospective analytical tool without real-time epidemiological integration; connection to hospital information systems, automated specimen tracking, or national surveillance dashboards would be required for fully prospective outbreak detection.

The 4-mer frequency representation captures compositional but not structural information; plasmids with similar nucleotide composition but different gene arrangements will receive similar codes, though the new recombination detection module (L5) partially mitigates this by flagging mosaic structures. Classification confidence remains reduced for plasmids highly divergent from the training set, as demonstrated by a Singaporean KPC-2 plasmid (60.9% confidence) [13], though the novelty detection module (L4) now explicitly flags such cases. The cross-validation against 74 plasmids from 27 published outbreak and surveillance studies across 13 countries provides quantitative evidence for clinical applicability: LOOCV on 57 plasmids with sequence data demonstrated 94.7% L6 concordance (54/57), 100% L3 cluster concordance (57/57), and cluster detection recall of 69.2% (9/13 studies) at L6. The mean classification confidence was 87.7%. However, prospective multi-centre validation in real-time clinical settings is needed to establish definitive utility.

The integration of pLIN with MLST addresses the fundamental limitation of plasmid-only surveillance: the inability to distinguish clonal spread from horizontal plasmid transfer. The retrospective validation (92.3% concordance, 7 species, 20 STs) demonstrates that combined typing provides actionable transmission-mode discrimination. However, paired chromosomal assemblies are required alongside plasmid sequences.

The mathematical stability of the pLIN framework supports its long-term universality: code permanence is guaranteed by the nearest-neighbour assignment rule; tetranucleotide frequency is an intrinsic sequence property independent of external databases; and the cosine distance approach is applicable to any DNA sequence. The successful extension to Gram-positive organisms and WHO critical priority pathogens -- where *S. aureus*, *Enterococcus*, *Acinetobacter baumannii*, and *Pseudomonas aeruginosa* plasmids have fundamentally different replication systems and GC-content profiles compared to Enterobacterales -- validates the generalisability of the 4-mer approach across the bacterial kingdom.

Future work should include prospective validation in hospital infection control workflows, extension of the training set to additional species (e.g., *Streptococcus pneumoniae*, *Clostridioides difficile*), integration with real-time laboratory information systems for automated prospective outbreak alerting, and development of a curated high-risk lineage watchlist linked to clinical action protocols. The novel Inc group discovery module (L6) provides a systematic pathway for continuous database expansion as new plasmid types are identified. The open-source, web-based design ensures that pLIN is immediately deployable in clinical microbiology laboratories without specialist bioinformatics infrastructure.

In conclusion, pLIN provides a standardised framework for plasmid-level AMR surveillance spanning Gram-negative and Gram-positive organisms that transforms plasmid sequence data into clinically actionable information. The integration of automated plasmid contig identification, hierarchical classification, assembly quality assessment, recombination detection, evolutionary rate estimation, and mobile genetic element mapping provides clinical microbiology laboratories with a comprehensive analytical platform that accepts raw WGS assemblies without requiring manual curation. By linking hierarchical plasmid identity to resistance profiles, quality metrics, and outbreak risk, pLIN enables infection prevention teams to move from reactive investigation to proactive, plasmid-informed infection control.

## Transparency declarations

**Competing interests:** None declared.

**Funding:** [To be completed]

**Data sharing:** pLIN source code, reference database, and training data are freely available under GPL-3.0 at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification. Hospitals deploying pLIN as a standalone clinical service are not subject to copyleft obligations.

## Author contributions

BBX: Conceptualisation, Methodology, Software, Validation, Formal Analysis, Writing -- Original Draft, Visualization. AKB: Data Curation. BS: Writing -- Review & Editing. JWAR: Supervision, Writing -- Review & Editing.

## References

1. Murray CJL, Ikuta KS, Sharara F, et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. Lancet 2022; 399: 629-55.
2. San Millan A. Evolution of plasmid-mediated antibiotic resistance in the clinical context. Trends Microbiol 2018; 26: 978-85.
3. Logan LK, Weinstein RA. The epidemiology of carbapenem-resistant Enterobacteriaceae: the impact and evolution of a global menace. J Infect Dis 2017; 215: S28-36.
4. Pitout JDD, Peirano G, Kock MM, Strydom KA, Matsumura Y. The global ascendency of OXA-48-type carbapenemases. Clin Microbiol Rev 2019; 33: e00102-19.
5. Carattoli A, Zankari E, Garcia-Fernandez A, et al. In silico detection and typing of plasmids using PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother 2014; 58: 3895-903.
6. Jolley KA, Bray JE, Maiden MCJ. Open-access bacterial population genomics: BIGSdb software, the PubMLST.org website and their applications. Wellcome Open Res 2018; 3: 124.
7. Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and typing of plasmids from draft assemblies. Microb Genom 2018; 4: e000206.
8. Redondo-Salvo S, Fernandez-Lopez R, Ruiz R, et al. Pathways for horizontal gene transfer in bacteria revealed by a global map of their plasmids. Nat Commun 2020; 11: 3602.
9. Weisberg AJ, Davis EW, Backman TWH, et al. Unexpected conservation and global transmission of agrobacterial virulence plasmids. Science 2020; 368: eaba5256.
10. Feldgarden M, Brover V, Gonzalez-Escalona N, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep 2021; 11: 12728.
11. Yao Y, Lazaro-Perona F, Falgenhauer L, et al. Insights into a novel blaKPC-2-encoding IncP-6 plasmid reveal carbapenem-resistance circulation in several Enterobacteriaceae species from wastewater and a hospital source in Germany. Microbiol Spectr 2023; 11: e0118523.
12. Weber RE, Pietsch M, Fruhauf A, et al. IS26-mediated transfer of blaNDM-1 as the main route of resistance transmission during a polyclonal, multispecies outbreak in a German hospital. Front Microbiol 2019; 10: 2817.
13. Marimuthu K, Venkatachalam I, Khong WX, et al. Clinical and molecular epidemiology of carbapenem-resistant Enterobacteriaceae among adult inpatients in Singapore. Clin Infect Dis 2017; 64: S68-75.
14. Roberts LW, Catchpoole E, Jennison AV, et al. Genomic analysis of carbapenemase-producing Enterobacteriaceae in Queensland reveals widespread transmission of blaIMP-4 on an IncHI2 plasmid. Microb Genom 2020; 6: e000321.
15. European Centre for Disease Prevention and Control. Systematic review of the effectiveness of infection control measures to prevent the transmission of carbapenemase-producing Enterobacteriaceae through cross-border transfer of patients. ECDC Technical Report 2014.
16. Magiorakos AP, Burns K, Rodriguez Bano J, et al. Infection prevention and control measures and tools for the prevention of entry of carbapenem-resistant Enterobacteriaceae into healthcare settings: guidance from the European Centre for Disease Prevention and Control. Antimicrob Resist Infect Control 2017; 6: 113.
17. World Health Organization. Global Antimicrobial Resistance and Use Surveillance System (GLASS) Report 2022. Geneva: WHO; 2022.
18. Navon-Venezia S, Kondratyeva K, Carattoli A. Klebsiella pneumoniae: a major worldwide source and shuttle for antibiotic resistance. FEMS Microbiol Rev 2017; 41: 252-75.
19. Peirano G, Pitout JDD. Extended-spectrum beta-lactamase-producing Enterobacteriaceae: update on molecular epidemiology and treatment options. Drugs 2019; 79: 1529-41.
20. Partridge SR, Kwong SM, Firth N, Jensen SO. Mobile genetic elements associated with antimicrobial resistance. Clin Microbiol Rev 2018; 31: e00088-17.
21. Rozwandowicz M, Brouwer MSM, Fischer J, et al. Plasmids carrying antimicrobial resistance genes in Enterobacteriaceae. J Antimicrob Chemother 2018; 73: 1121-37.
22. Sheppard AE, Stoesser N, Wilson DJ, et al. Nested Russian doll-like genetic mobility drives rapid dissemination of the carbapenem resistance gene blaKPC. Antimicrob Agents Chemother 2016; 60: 3767-78.
23. Conlan S, Thomas PJ, Deming C, et al. Single-molecule sequencing to track plasmid diversity of hospital-associated carbapenemase-producing Enterobacteriaceae. Sci Transl Med 2014; 6: 254ra126.
24. Acman M, Wang R, van Dorp L, et al. Role of mobile genetic elements in the global dissemination of the carbapenem resistance gene blaNDM. Nat Commun 2022; 13: 1131.
25. Liu YY, Wang Y, Walsh TR, et al. Emergence of plasmid-mediated colistin resistance mechanism MCR-1 in animals and human beings in China: a microbiological and molecular biological study. Lancet Infect Dis 2016; 16: 161-8.
26. Wang R, van Dorp L, Shaw LP, et al. The global distribution and spread of the mobilized colistin resistance gene mcr-1. Nat Commun 2018; 9: 1179.
27. Mathers AJ, Peirano G, Pitout JDD. The role of epidemic resistance plasmids and international high-risk clones in the spread of multidrug-resistant Enterobacteriaceae. Clin Microbiol Rev 2015; 28: 565-91.
28. Orlek A, Stoesser N, Sheridan E, et al. Plasmid classification in an era of whole-genome sequencing: application in studies of antibiotic resistance epidemiology. Front Microbiol 2017; 8: 182.
29. David S, Reuter S, Harris SR, et al. Epidemic of carbapenem-resistant Klebsiella pneumoniae in Europe is driven by nosocomial spread. Nat Microbiol 2019; 4: 1919-29.
30. Wyres KL, Holt KE. Klebsiella pneumoniae as a key trafficker of drug resistance genes from environmental to clinically important bacteria. Curr Opin Microbiol 2018; 45: 131-9.
31. Seemann T. mlst: scan contig files against PubMLST typing schemes. https://github.com/tseemann/mlst. 2023.
32. Kwong SM, Ramsay JP, Jensen SO, Firth N. Replication of staphylococcal resistance plasmids. Front Microbiol 2017; 8: 2279.
33. Palmer KL, Kos VN, Gilmore MS. Horizontal gene transfer and the genomics of enterococcal antibiotic resistance. Curr Opin Microbiol 2010; 13: 632-9.
34. Mathers AJ, Stoesser N, Sheppard AE, et al. Klebsiella pneumoniae carbapenemase (KPC)-producing K. pneumoniae at a single institution: insights into endemicity from whole-genome sequencing. Antimicrob Agents Chemother 2015; 59: 1656-63.
35. Schmartz GP, Mangold KA, Raden M, et al. PLSDB: advancing a comprehensive database of bacterial plasmids. Nucleic Acids Research 2022; 50(D1): D273-D278.

---

**Tables:** See Tables 1--3

**Figures:** See Figures 1--10 with legends
