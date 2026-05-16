# Table 4: Machine learning classifier performance for Inc/Rep group prediction

All models were trained on 256-dimensional tetranucleotide (4-mer) frequency vectors computed from 8,077 plasmid sequences across 28 Inc/Rep groups (20 Gram-negative Inc groups, 4 Gram-positive rep type groups, 2 *Acinetobacter* rep type groups, and 2 *Pseudomonas* rep type groups). Performance was assessed using nested cross-validation (outer: 5-fold stratified; inner: 5-fold for hyperparameter tuning where applicable).

| Model | Weighted F1 (mean +/- SD) | Key hyperparameters |
|-------|:---:|---|
| KNN (primary classifier) | 0.911* | k=5, cosine distance, distance-weighted |
| XGBoost | 0.896 +/- 0.011 | max_depth=6, n_estimators=300, learning_rate=0.1 |
| Gradient Boosting | 0.893 +/- 0.009 | max_depth=5, n_estimators=200, learning_rate=0.1 |
| Random Forest | 0.874 +/- 0.021 | n_estimators=500, max_features=sqrt |
| Logistic Regression | 0.866 +/- 0.012 | C=1.0, multi_class=multinomial, solver=lbfgs |

\* KNN accuracy reported as overall accuracy from 5-fold stratified cross-validation (91.1%). The weighted F1 score is comparable. The dominant confusion pairs involve near-identical Inc groups: IncFIB and IncFII (median inter-group cosine distance d=0.002), reflecting genuine biological overlap of IncF-family replicons rather than methodological failure. Among Gram-positive groups, repSA_large and repEF_res showed occasional cross-classification (d=0.008).

**Rationale for KNN as primary classifier:** KNN was selected because (1) it naturally integrates with the nearest-neighbour pLIN assignment framework, (2) it provides interpretable distance-based confidence scores, and (3) its performance is comparable to ensemble methods while requiring no hyperparameter tuning beyond k and the distance metric.
