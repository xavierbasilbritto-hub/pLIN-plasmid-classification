# Appendix p 8: Machine learning validation — detailed methodology and results

## Feature engineering

A 33-dimensional feature vector was computed for each plasmid comprising:
- **Basic composition (7 features):** sequence length, log10(length), GC content, AT content, AT skew [(A-T)/(A+T)], GC skew [(G-C)/(G+C)]
- **Dinucleotide frequencies (16 features):** all 16 canonical dinucleotides, normalised by (length - 1)
- **Trinucleotide frequencies (10 features):** ATG, TAA, TAG, TGA, GCG, CGC, AAA, TTT, CCC, GGG — selected for biological relevance (start/stop codons, CpG motifs, homopolymeric tracts)

## Nested cross-validation design

| Parameter | Value |
|-----------|-------|
| Outer loop | 3-fold stratified CV (StratifiedKFold, shuffle=True, seed=42) |
| Inner loop | 2-fold stratified CV |
| Scoring metric | Weighted F1 (f1_weighted) |
| Balanced subset | 1,500 plasmids (500 IncFII, 500 IncN, 500 IncX1) |
| Feature scaling | StandardScaler (z-score), fitted per fold |

## Hyperparameter search spaces (Optuna TPE sampler, Hyperband pruning, 10 trials/fold)

| Model | Hyperparameter | Search range |
|-------|---------------|-------------|
| XGBoost | n_estimators | [50, 500] |
| | max_depth | [3, 15] |
| | learning_rate | [0.01, 0.3] (log-uniform) |
| | subsample | [0.6, 1.0] |
| | colsample_bytree | [0.6, 1.0] |
| | gamma | [1e-8, 1.0] (log-uniform) |
| | reg_alpha, reg_lambda | [1e-8, 10.0] (log-uniform) |
| Gradient Boosting | n_estimators | [50, 500] |
| | max_depth | [3, 15] |
| | learning_rate | [0.01, 0.3] (log-uniform) |
| | subsample | [0.6, 1.0] |
| Random Forest | n_estimators | [50, 500] |
| | max_depth | [3, 20] |
| | min_samples_split | [2, 20] |
| | min_samples_leaf | [1, 10] |
| Logistic Regression | C | [1e-3, 100] (log-uniform) |
| | penalty | {l1, l2} |
| | solver | saga (fixed) |

## Top 10 consensus feature importances

| Rank | Feature | Mean importance | Biological relevance |
|------|---------|----------------|---------------------|
| 1 | tri_TAG | 0.114 | Amber stop codon frequency; reflects codon usage bias |
| 2 | tri_GCG | 0.099 | CpG island proxy; DNA methylation signatures |
| 3 | at_content | 0.076 | Overall AT composition; distinguishes backbone vs accessory |
| 4 | gc_content | 0.053 | GC%; correlates with host range and evolutionary origin |
| 5 | di_AA | 0.050 | Poly-A tracts; regulatory element density |
| 6 | tri_TGA | 0.048 | Opal stop codon; alternative stop codon usage |
| 7 | tri_CGC | 0.040 | CpG-related; restriction-modification system signatures |
| 8 | di_TT | 0.039 | Poly-T tracts; rho-independent terminator signals |
| 9 | di_GG | 0.032 | G-richness; G-quadruplex potential |
| 10 | tri_GGG | 0.030 | Extreme GC-rich tracts |

**Clinical interpretation:** The dominance of stop codon and CpG-related features confirms that plasmid classification by 4-mer composition captures genuine phylogenetic signals -- not statistical noise. The strong performance of the linear baseline (Logistic Regression F1=0.866) indicates that ~87% of the classificatory information is linearly separable, with ensemble methods capturing an additional ~3% from nonlinear feature interactions.
