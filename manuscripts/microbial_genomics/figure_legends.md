# Figure Legends

## Figure 1. Dataset overview and pLIN classification summary.

**(a)** Distribution of 8,077 complete plasmid sequences across 28 plasmid groups (20 Gram-negative Inc groups, 4 Gram-positive rep type groups, 2 *Acinetobacter baumannii* rep type groups, and 2 *Pseudomonas aeruginosa* rep type groups). IncFII dominates the dataset (n = 4,629; 57.3%), followed by IncN (n = 1,097; 13.6%) and IncX1 (n = 705; 8.7%). The four Gram-positive groups (repSA_large n = 121, repEF_res n = 121, repEF_conj n = 92, repSA_small n = 73) and four non-fermentative groups (repAci1, repAci_large, repPae_large, repPae_small) are highlighted with distinct colouring. **(b)** Number of unique strain-level (L6) pLIN codes per group, demonstrating that pLIN captures substantial within-group diversity. **(c)** Proportion of singleton versus multi-member pLIN codes at strain level: 2,335/3,073 codes (76.0%) are singletons. The high singleton rate reflects the extensive genomic diversity of plasmid populations. **(d)** GC-content distribution comparison between Gram-negative (mean 49.1%) and Gram-positive (mean 33.2--35.8%) plasmid groups, illustrating the widened GC validation range (25--70%). See Table 3 for full dataset statistics.

## Figure 2. AMR gene prevalence across the plasmid dataset.

**(a)** Overall detection summary: 64,891 total gene detections (29,583 AMR, 6,286 virulence, 29,022 stress) were identified in 5,816/6,998 Gram-negative plasmids (83.1%) using AMRFinderPlus v4.2.5. **(b)** Top 20 most frequently detected AMR genes, with *bla*TEM-1 (n = 1,864; 40.0%) being the most prevalent, followed by *sul1* (29.7%) and *tet*(A) (27.9%). **(c)** AMR drug class distribution showing beta-lactam resistance as the most common category (3,483 plasmids), followed by aminoglycoside (2,723), sulfonamide (2,073), trimethoprim (1,810), and tetracycline (1,549) resistance determinants. See Table 5 for detailed AMR statistics.

## Figure 3. Clinically critical resistance determinants mapped to plasmid groups.

Distribution of high-priority resistance genes across the 20 Gram-negative Inc groups: carbapenemases (n = 1,635; *bla*KPC-2 n = 824, *bla*NDM-1 n = 228, *bla*KPC-3 n = 193), extended-spectrum beta-lactamases (n = 1,804; *bla*CTX-M-15 n = 505, *bla*CTX-M-65 n = 319, *bla*SHV-12 n = 277), colistin resistance *mcr* genes (n = 204; *mcr-1.1* n = 83), and plasmid-mediated quinolone resistance (n = 2,315; *qnrS1* n = 732, *aac(6')-Ib-cr5* n = 737). AMR analysis was performed only on the 6,998 Gram-negative plasmids; no AMR data are available for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. Stacked bars indicate relative contribution of each Gram-negative group to each resistance category. See Table 5 Panel C.

## Figure 4. pLIN lineage-AMR heatmap showing resistance gene profiles of the top pLIN lineages.

Heatmap showing the prevalence of the 15 most common AMR genes (columns) across the top 10 pLIN lineages ranked by AMR-positive plasmid count (rows). Colour intensity represents the proportion of lineage members carrying each gene (0--100%). pLIN 671 (IncN, n = 90) shows 100% carriage of *bla*KPC-2 and *bla*TEM-1 with 13.2 mean AMR genes. pLIN 860 (5 Inc groups, n = 142) shows the highest mean AMR burden (14.4 genes per plasmid) with a distinctive multi-drug resistance signature spanning sulfonamides (*sul1* 61.3%), macrolides (*mph(A)*/*mrx(A)* 54.2%), aminoglycosides (*aph(6)-Id* 52.1%), tetracyclines (*tet(A)* 52.8%), and beta-lactams (*blaTEM-1* 38.0%); this lineage also carries 44.4% *mcr* colistin resistance (not shown among the top 15 genes). White cells indicate absence of the gene in that lineage. See Table 6 Panel A for numerical values.

## Figure 5. AMR burden distribution and virulence gene specificity by plasmid group.

**(a)** Boxplot of AMR gene counts per plasmid stratified by plasmid group, showing IncAC2 with the highest median burden (12.6 mean genes) and ColRNAI with the lowest (0.63 mean genes). Five groups (IncA, IncAC2, IncC, IncFIC, IncHI2) showed 100% AMR carriage. AMR data are shown only for the 20 Gram-negative Inc groups (6,998 plasmids); no AMR analysis was performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. **(b)** Virulence gene prevalence by group, with IncF plasmids showing the highest rate (62.7%) driven by *traT* and *spv* operon genes. **(c)** Group specificity of virulence determinants: IncFII plasmids enriched for *traT* and *spv* genes; IncX1 enriched for *hlyA* and fimbrial adhesin genes; IncN showing minimal virulence carriage (2.9%).

## Figure 6. pLIN hierarchical clustering structure.

Dendrogram illustrating the six-level hierarchical structure of pLIN codes for all 8,077 plasmids. Horizontal dashed lines indicate the six distance thresholds: L1 (d <= 0.150, ~85% ANI), L2 (d <= 0.100, ~90% ANI), L3 (d <= 0.050, ~95% ANI), L4 (d <= 0.020, ~98% ANI), L5 (d <= 0.010, ~99% ANI), and L6 (d <= 0.001, ~99.9% ANI). Gram-positive groups form distinct subtrees separated from Gram-negative clusters at L2--L3, reflecting their divergent composition. At L1--L2, most plasmids form a single cluster. Meaningful separation emerges at L3 and increases progressively to strain-level groups at L6. See Table 2 for threshold definitions.

## Figure 7. Composite manuscript figure showing the integrated pLIN analytical pipeline.

Overview of the pLIN system combining: **(a)** plasmid versus chromosome contig identification via multi-signal scoring (sequence length, cosine distance to nearest training vector, header keyword matching, KNN confidence; score >= 10 plasmid, <= -10 chromosome), **(b)** tetranucleotide frequency computation and cosine distance matrix, **(c)** single-linkage hierarchical clustering at six thresholds, **(d)** pLIN code assignment, **(e)** group concordance (97.8%), **(f)** AMRFinderPlus integration, **(g)** lineage-level AMR profiling, **(h)** assembly completeness assessment, **(i)** recombination detection, and **(j)** MGE boundary detection. Arrows indicate data flow from FASTA input through contig identification, classification, and quality assessment to epidemiological output.

## Figure 8. pLIN analysis pipeline overview.

Flowchart depicting the expanded pLIN workflow: (1) FASTA input, (2) plasmid versus chromosome contig identification via multi-signal scoring (contigs classified as plasmid [score >= 10], chromosome [score <= -10], or borderline [defaults to plasmid]; chromosomal contigs excluded from further analysis), (3) assembly completeness assessment (composite score 0--100), (4) tetranucleotide frequency vector computation (256 features), (5) pairwise cosine distance computation, (6) single-linkage hierarchical clustering, (7) six-threshold dendrogram cutting, (8) pLIN code assembly, (9) group classification (KNN, k = 5, 91.1% accuracy, 28 groups), (10) database coverage and novelty detection (traffic-light system), (11) AMRFinderPlus gene detection, (12) recombination detection via minimap2, (13) MGE boundary detection and composite transposon identification, (14) pLIN-AMR integration, (15) lineage-level AMR profiling, (16) evolutionary rate estimation, (17) adaptive threshold and cluster stability assessment, (18) novel group discovery for low-confidence plasmids, (19) outbreak detection (L6 + AMR fingerprint + temporal), and (20) interactive Streamlit web interface output.

## Figure 9. Threshold calibration: cosine distance versus ANI relationship.

**(a)** Nearest-neighbour distance distribution of 178 IncX-like reference plasmids used for threshold calibration, with key quantiles marked (50th percentile d = 0.011, 90th d = 0.038, 99th d = 0.072). **(b)** Scatter plot of cosine distance versus FastANI for 4,970 within-group plasmid pairs from the 20 Gram-negative Inc groups, showing the monotonic negative correlation (Spearman rho = -0.348, P < 10^-141; 15/20 groups significant at P < 0.05). The six pLIN thresholds are overlaid as vertical lines with their approximate ANI equivalents.

## Figure 10. Comparison of plasmid classification approaches.

Radar chart or comparative matrix visualising the properties of six plasmid typing systems (pLIN, PlasmidFinder, pMLST, MOB-suite, COPLA, mge-cluster) across key criteria: hierarchical levels, code permanence, reference-free operation, taxonomic scope (updated to include Gram-positive coverage), resolution, discriminatory power, AMR integration, quality assessment modules, and computational cost. pLIN is the only system that satisfies all criteria simultaneously. See Table 1 for detailed comparison.

## Figure 11. Gram-positive classifier expansion: composition space and validation.

**(a)** PCA projection of 8,077 plasmid tetranucleotide frequency vectors coloured by group, showing clear separation between Gram-negative Inc groups, Gram-positive rep type groups, and non-fermentative rep type groups along PC1 (driven primarily by GC content). **(b)** Confusion matrix for the 28-group KNN classifier (91.1% accuracy), with the 4 Gram-positive groups and 4 non-fermentative groups (*A. baumannii* and *P. aeruginosa*) highlighted. Inter-kingdom misclassification rate < 0.5%. **(c)** Per-group classification accuracy: repSA_large 93.4%, repEF_res 91.7%, repEF_conj 90.2%, repSA_small 89.0%. **(d)** GC-content distribution across all 28 groups illustrating the expanded validation range (25--70%).

## Figure 12. Assembly completeness and database coverage quality modules.

**(a)** Distribution of composite assembly completeness scores across 8,077 training sequences: 92.0% COMPLETE (>= 80), with breakdown by scoring component (contig count, N50 ratio, circularity, coding density, N-gaps). **(b)** Nearest-neighbour distance percentile distributions for all 28 groups, with traffic-light thresholds (GREEN/YELLOW/RED) overlaid. **(c)** Relationship between completeness score and classification confidence for the 74 outbreak validation plasmids, confirming that lower completeness scores predict reduced classification reliability.

## Figure 13. Plasmid versus chromosome contig classification validation.

**(a)** Receiver operating characteristic (ROC) curve for the multi-signal contig classifier evaluated on a mixed test set of 8,077 known plasmids and 500 chromosomal contigs, showing 98.7% sensitivity and 97.4% specificity. The composite score threshold (score >= 10 plasmid, score <= -10 chromosome) is marked on the curve. **(b)** Contribution of each scoring signal stratified by contig size range: sequence length dominates classification for contigs > 500 kb and < 20 kb, while cosine distance to the nearest plasmid training vector is the most informative signal in the 20--500 kb range. **(c)** Confusion matrix showing classification outcomes: 7,975 true plasmid, 487 true chromosome, 102 false chromosome (predominantly large conjugative plasmids > 200 kb from IncHI1/IncHI2), and 13 false plasmid (predominantly small chromosomal fragments < 50 kb). **(d)** Distribution of composite plasmid scores for true plasmids (green) and true chromosomes (red), with the decision thresholds (score = 10 and score = -10) and borderline zone shaded. The 214 borderline contigs (2.7%) are shown in the inset, with 92.5% verified as true plasmids, confirming the validity of the conservative borderline-to-plasmid default.

## Figure 14. Reference database expansion: scaling pLIN to 79,305 plasmids.

**(a)** Group distribution in the expanded 79,305-plasmid database, dominated by IncFII (34,036; 43.1%) and IncX1 (26,154; 33.1%). **(b)** Fold increase in cluster counts from training-only (n = 8,056 unique) to full database (n = 79,305; 8,056 unique training + 71,249 reference) at each hierarchical level, reaching 57,886 unique L6 codes. **(c)** Runtime breakdown: 23.5 min for 4-mer computation, 7 s for KNN classification, 4.2 min for per-group clustering. Total < 30 min on standard hardware. Classification rate: 97.3% (2,109 Unknown/Novel, 2.7%).

## Figure 15. Outbreak cross-validation against four published studies.

**(a)** Study 1 (Yao *et al.* 2023): KPC-2 IncN plasmid CP104944 assigned pLIN 671 (d = 0.0000), matching the known high-risk lineage (n = 90, 100% *bla*KPC-2). **(b)** Study 2 (Weber *et al.* 2019): 12 NDM-1 outbreak plasmids resolved into 9 L6 codes but 1 L3 cluster; 4 plasmids share pLIN 492. **(c)** Study 3 (Marimuthu *et al.* 2022): Singapore pKPC2 plasmid MN542377 assigned novel pLIN 2455 (IncN, 60.9% confidence, d = 0.0058). **(d)** Study 4 (Roberts *et al.* 2020): Australian IncHI2 outbreak plasmid CP022533 assigned pLIN 860 (d = 0.0004), matching the cross-Inc MDR hub. See Table 6 Panel B for detailed results.

## Figure 16. Resolution comparison across plasmid classification methods.

Alluvial flow diagram showing how 74 outbreak plasmids from 27 studies across 13 countries are resolved by six classification methods, ordered by increasing resolution. Ribbons coloured by resistance mechanism. PlasmidFinder: 11 flat Inc labels (100% coverage); pMLST: 10 sequence types (61% coverage); MOB-suite: 8 cluster IDs (61% coverage, grey 'unclassified' bands); COPLA/PTU: 6 PTU groups (57% coverage); mge-cluster: 22 clusters (100% coverage, no permanent codes, dashed borders); pLIN: 42 unique hierarchical codes (100% coverage). Callout annotations: Ho 2019 NDM (12/13 -> pLIN 475), OXA-48 cross-country (5/6 -> pLIN 1688), mcr-1 separation (pLIN 87 vs 340).

## Figure 17. Combined chromosomal-plasmid typing: transmission mode discrimination.

**(a)** Bipartite network showing MLST sequence type (ST) nodes (left, coloured by host species) connected to pLIN L6 code nodes (right, green squares). Edges represent co-occurrence within outbreak studies. Seven host species represented: *E. coli*, *K. pneumoniae*, *E. cloacae*, *C. freundii*, *E. hormaechei*, *S. enterica*, *A. baumannii*. **(b)** Stacked horizontal bar chart showing pairwise transmission mode distribution (clonal spread, horizontal transfer, independent) across 13 studies. **(c)** Summary validation metrics: 13 studies, 7 species, 20 STs, 92.3% overall concordance (12/13), 100% multi-ST detection, 92.3% clonal detection.

## Figure 18. Recombination detection and mosaic plasmid characterisation.

**(a)** Distribution of recombination flags (None/Low/Medium/High) across 705 IncX1 plasmids: 56.9% None, 15.3% Low, 18.7% Medium, 9.1% High. **(b)** Comparison of AMR gene burden between recombination categories: High-recombination plasmids carry significantly more AMR genes (mean 9.7) than None (mean 4.2; P < 10^-8). **(c)** Example minimap2 PAF alignment visualisation for a High-recombination mosaic plasmid, showing fragmented alignment blocks (coloured by reference match) interspersed with novel insertions (grey). **(d)** Concordance between minimap2-based recombination flags and GC-content-based mosaicism predictions (78.3% agreement).

## Figure 19. Novel group discovery among Unknown/Novel sequences.

**(a)** Hierarchical clustering dendrogram of 2,109 Unknown/Novel sequences (KNN confidence < 40%) at the L3 threshold, with 47 candidate novel groups (>= 3 members) highlighted. **(b)** Size distribution of candidate novel groups (range 3--64 members, median 8). **(c)** Nearest known group distances for each candidate, with the largest candidate group (n = 64) showing nearest distance d = 0.089 to IncR. **(d)** GC-content and plasmid size distributions for the three largest candidate novel groups compared to established groups.

## Figure 20. Evolutionary rate estimation within pLIN lineages.

**(a)** Distribution of estimated substitution rates across 87 L6 clusters with significant temporal signal (R-squared >= 0.3): median 3.2 x 10^-6 subs/site/year (IQR 1.1--8.7 x 10^-6), consistent with published plasmid rates. Grey shaded region indicates the expected range (10^-6 to 10^-5 subs/site/year). **(b)** Exemplar regression plot for pLIN 671 (KPC-2 IncN lineage): 7.4 x 10^-6 subs/site/year, R-squared = 0.72, suggesting active diversification. **(c)** Twelve outlier clusters with rates exceeding 10^-5 subs/site/year, flagged as potential recombination-driven outliers.

## Figure 21. Adaptive thresholds and cluster stability analysis.

**(a)** Bootstrap stability scores (50 iterations) across all 28 plasmid groups, ranging from 0.78 (IncFII) to 0.97 (IncX3); overall mean 0.88. Gram-positive groups show high stability (mean 0.93). **(b)** Adjusted Rand Index comparing single-linkage, complete-linkage, and average-linkage clustering at each of the six thresholds: single versus average linkage mean ARI = 0.91; single versus complete linkage mean ARI = 0.76. **(c)** Group-specific calibrated threshold deviations from global thresholds (mean 12.3%, range 0.8--31.2%), with the largest deviations in IncFII and IncHI2.

## Figure 22. Mobile genetic element boundary detection and composite transposon identification.

**(a)** IS element density across 8,077 plasmids: mean 1.93 IS elements per plasmid, with significantly higher density in AMR-positive plasmids (2.47 vs 0.84; P < 10^-20). **(b)** Frequency distribution of composite transposon types: Tn*4401*-like/*bla*KPC (n = 387), IS*26*-flanked MDR cassettes (n = 312), Tn*10*-like/*tet*(B) (n = 198). **(c)** Example colour-coded linear gene map of a pLIN 671 member showing Tn*4401* composite transposon carrying *bla*KPC-2 flanked by IS*Kpn7* elements, with backbone genes (blue), IS elements (orange), resistance genes (red), and integrases (green) distinguished by colour. **(d)** Integrase and recombinase distribution across plasmid groups (2,891 total detections; mean 0.39 per plasmid).

## Supplementary Figure S1. Cosine distance to ANI validation across all 20 Gram-negative Inc groups.

Per-Inc-group scatter plots of cosine distance versus FastANI (n = 4,970 pairs total from 20 Gram-negative Inc groups; up to 20 plasmids per group, all pairwise comparisons; Spearman rho = -0.348, P < 10^-141; 15/20 groups significant at P < 0.05). No FastANI validation was performed for the Gram-positive, *Acinetobacter*, or *Pseudomonas* groups. Strongest correlations observed in IncFIC (rho = -0.88), IncI2 (rho = -0.81), ColE (rho = -0.72), and IncN (rho = -0.61). Five groups (IncX3, IncX4, IncFIBK, IncFII, IncX1) show non-significant correlations due to range restriction (all within-group ANI values exceed 93%). At d <= 0.001, median FastANI = 99.9% across all groups.

## Supplementary Figure S2. Gram-positive plasmid group characterisation.

**(a)** Plasmid size distributions for the four Gram-positive rep type groups: repSA_large (median 28.4 kb), repSA_small (median 4.6 kb), repEF_conj (median 62.1 kb), and repEF_res (median 45.3 kb). **(b)** Within-group pairwise cosine distance distributions compared to the 20 Gram-negative groups. **(c)** Representative pLIN dendrogram subtrees for each Gram-positive group, illustrating within-group diversity and hierarchical structure. Note: repEF_res contains 100 unique plasmids (121 training samples including 21 multi-replicon plasmids shared with repEF_conj). The Gram-positive groups contributed 255 unique L6 codes: repSA_large 72, repSA_small 52, repEF_conj 63, repEF_res 69. **(d)** No AMR analysis was performed for the Gram-positive groups; AMR data are available only for the 20 Gram-negative Inc groups (6,998 plasmids).
