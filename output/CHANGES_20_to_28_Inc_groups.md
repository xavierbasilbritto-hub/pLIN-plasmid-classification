# Changes: "20 Inc groups" → "20 Gram-negative Inc groups"

Use Find & Replace in your manuscript editor to apply these changes.
Each entry shows: FILE → LINE NUMBER → FIND → REPLACE

---

## 1. Appendix p01 — Title (Line 1) — ALL 4 JOURNALS + output/

**FIND:**
```
Complete dataset breakdown for all 20 Inc groups
```
**REPLACE WITH:**
```
Complete dataset breakdown for all 20 Gram-negative Inc groups
```

Files:
- `manuscripts/lancet_microbe/supplementary/Appendix_p01_Dataset_breakdown_20_Inc_groups.md` → L1
- `manuscripts/cmi/supplementary/Appendix_p01_Dataset_breakdown_20_Inc_groups.md` → L1
- `manuscripts/microbial_genomics/supplementary/Appendix_p01_Dataset_breakdown_20_Inc_groups.md` → L1
- `manuscripts/briefings_bioinformatics/supplementary/Appendix_p01_Dataset_breakdown_20_Inc_groups.md` → L1
- `output/lancet_microbe_v2_tables_appendix/Appendix_p01_Dataset_breakdown_20_Inc_groups.md` → L1

---

## 2. Appendix p09 — Section heading (Line 26) — ALL 4 JOURNALS + output/

**FIND:**
```
Multi-Inc validation (all 20 Inc groups, n=4,970 pairs)
```
**REPLACE WITH:**
```
Multi-Inc validation (all 20 Gram-negative Inc groups, n=4,970 pairs)
```

Files:
- `manuscripts/lancet_microbe/supplementary/Appendix_p09_Hierarchical_threshold_calibration.md` → L26
- `manuscripts/cmi/supplementary/Appendix_p09_Hierarchical_threshold_calibration.md` → L26
- `manuscripts/microbial_genomics/supplementary/Appendix_p09_Hierarchical_threshold_calibration.md` → L26
- `manuscripts/briefings_bioinformatics/supplementary/Appendix_p09_Hierarchical_threshold_calibration.md` → L26
- `output/lancet_microbe_v2_tables_appendix/Appendix_p09_Hierarchical_threshold_calibration.md` → L26

---

## 3. Appendix p09 — Body text (Line 28) — ALL 4 JOURNALS + output/

**FIND:**
```
sampled from all 20 Inc groups (up to 20 randomly
```
**REPLACE WITH:**
```
sampled from all 20 Gram-negative Inc groups (up to 20 randomly
```

Files: (same 5 files as #2 above, Line 28)

---

## 4. Appendix p12 — Runtime line (Line 45 or 47) — ALL 4 JOURNALS + output/

**FIND:**
```
Phase 3 (per-group clustering, 20 groups)
```
**REPLACE WITH:**
```
Phase 3 (per-group clustering, 20 Gram-negative groups)
```

Files:
- `manuscripts/lancet_microbe/supplementary/Appendix_p12_Reference_database_expansion_statistics.md` → L45
- `manuscripts/cmi/supplementary/Appendix_p12_Reference_database_expansion_statistics.md` → L45
- `manuscripts/microbial_genomics/supplementary/Appendix_p12_Reference_database_expansion_statistics.md` → L45
- `manuscripts/briefings_bioinformatics/supplementary/Appendix_p12_Reference_database_expansion_statistics.md` → L47
- `output/lancet_microbe_v2_tables_appendix/Appendix_p12_Reference_database_expansion_statistics.md` → L45

---

## 5. Table5 Reference Database (output/ only) — Line 7

**FIND:**
```
| Inc groups | 20 |
```
**REPLACE WITH:**
```
| Inc/Rep groups | 28 |
```

File:
- `output/lancet_microbe_v2_tables_appendix/Table5_Reference_database_78902_plasmids.md` → L7

---

## 6. Lancet Table2 — Intro paragraph (Line 3)

**FIND:**
```
validated against FastANI across 4,970 plasmid pairs from all 20 Inc groups.
```
**REPLACE WITH:**
```
validated against FastANI across 4,970 plasmid pairs from all 20 Gram-negative Inc groups.
```

File:
- `manuscripts/lancet_microbe/tables/Table2_thresholds.md` → L3

---

## 7. Lancet Table2 — FastANI validation note (Line 14)

**FIND:**
```
for 15 of 20 Inc groups
```
**REPLACE WITH:**
```
for 15 of 20 Gram-negative Inc groups
```

File:
- `manuscripts/lancet_microbe/tables/Table2_thresholds.md` → L14

---

## 8. Lancet Table3 — Title (Line 1)

**FIND:**
```
6,998 complete plasmid sequences across 20 Inc groups
```
**REPLACE WITH:**
```
6,998 complete plasmid sequences across 20 Gram-negative Inc groups
```

File:
- `manuscripts/lancet_microbe/tables/Table3_dataset.md` → L1

---

## 9. MG Table1 — Footnote (Line 21)

**FIND:**
```
pMLST schemes not available for all 20 Inc groups simultaneously.
```
**REPLACE WITH:**
```
pMLST schemes not available for all 20 Gram-negative Inc groups simultaneously.
```

File:
- `manuscripts/microbial_genomics/tables/Table1_comparison.md` → L21

---

## 10. BB Table1 — Footnote (Line 21)

**FIND:**
```
pMLST schemes not available for all 20 Inc groups simultaneously.
```
**REPLACE WITH:**
```
pMLST schemes not available for all 20 Gram-negative Inc groups simultaneously.
```

File:
- `manuscripts/briefings_bioinformatics/tables/Table1_comparison.md` → L21

---

## Summary

**In every case, the change is simply inserting the word "Gram-negative" before "Inc groups" (or "groups").**

Quick Find & Replace across all files:
1. `all 20 Inc groups` → `all 20 Gram-negative Inc groups`
2. `20 groups)` → `20 Gram-negative groups)`
3. `of 20 Inc groups` → `of 20 Gram-negative Inc groups`
4. `across 20 Inc groups` → `across 20 Gram-negative Inc groups`
5. `| Inc groups | 20 |` → `| Inc/Rep groups | 28 |`

**Total: 23 locations across 23 files**
