# Table 4. Machine learning model performance for Inc-group prediction from tetranucleotide composition features (3-class balanced subset, nested cross-validation). KNN classifier: 91.1% accuracy across 28 groups, 8,077 training samples.

| Model | Mean weighted F1 | SD | Min | Max |
|-------|------------------|----|-----|-----|
| XGBoost | 0.896 | 0.009 | 0.889 | 0.909 |
| Gradient Boosting | 0.893 | 0.007 | 0.888 | 0.904 |
| Random Forest | 0.874 | 0.019 | 0.855 | 0.900 |
| Logistic Regression | 0.866 | 0.010 | 0.853 | 0.877 |

A balanced subset of 1,500 plasmids (500 per group: IncFII, IncN, IncX1) was used. Nested cross-validation: 3 outer folds (stratified), 2 inner folds for hyperparameter tuning. Optimisation: Optuna with TPE sampler (seed = 42) and Hyperband pruning (10 trials per fold). Scoring: weighted F1. All models were trained on a 33-dimensional feature vector comprising basic composition metrics (sequence length, log-length, GC content, AT content, AT skew, GC skew), 16 dinucleotide frequencies, and 10 selected trinucleotide frequencies (ATG, TAA, TAG, TGA, GCG, CGC, AAA, TTT, CCC, GGG).

## Top 10 consensus feature importances (averaged across tree-based models)

| Rank | Feature | Mean importance | Biological relevance |
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
