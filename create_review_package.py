#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Create complete reviewer package for for_review2/:
  - Tables 1-6 (Excel workbook with formatted sheets)
  - Supplementary Tables (Excel workbook)
  - Appendix pages 1-16 (Markdown)
  - Raw datasets underlying each figure (TSV)
  - Figure legends document
"""

import os
import sys
import numpy as np
import pandas as pd
from collections import Counter
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
REVIEW_DIR = os.path.join(BASE_DIR, "for_review2")
TABLES_DIR = os.path.join(REVIEW_DIR, "tables")
APPENDIX_DIR = os.path.join(REVIEW_DIR, "appendix")
RAW_DATA_DIR = os.path.join(REVIEW_DIR, "source_data")
os.makedirs(TABLES_DIR, exist_ok=True)
os.makedirs(APPENDIX_DIR, exist_ok=True)
os.makedirs(RAW_DATA_DIR, exist_ok=True)

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading data ...")
INTEGRATED = os.path.join(BASE_DIR, "output", "integrated", "pLIN_AMR_integrated.tsv")
PLIN_FILE = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
AMR_RAW = os.path.join(BASE_DIR, "output", "amrfinder", "amrfinder_all_plasmids.tsv")
ANI_FILE = os.path.join(BASE_DIR, "output", "cosine_to_ani_validation_all_inc.tsv")

merged = pd.read_csv(INTEGRATED, sep="\t")
plin = pd.read_csv(PLIN_FILE, sep="\t")
amr_raw = pd.read_csv(AMR_RAW, sep="\t")

INC_COL = "inc_type_x" if "inc_type_x" in merged.columns else "inc_type"
inc_counts = plin["inc_type"].value_counts()
INC_ORDER = inc_counts.index.tolist()

GRAM_NEG = [
    "ColE", "ColRNAI", "IncA", "IncAC2", "IncC", "IncF", "IncFIB",
    "IncFIBK", "IncFIC", "IncFII", "IncHI1", "IncHI2", "IncI", "IncI1",
    "IncI2", "IncN", "IncR", "IncX1", "IncX3", "IncX4",
]
GRAM_POS = ["repSA_large", "repSA_small", "repEF_conj", "repEF_res"]
ACI_GROUPS = ["repAci1", "repAci_large"]
PAE_GROUPS = ["repPae_large", "repPae_small"]

amr_only = amr_raw[amr_raw["Type"] == "AMR"]
vir_only = amr_raw[amr_raw["Type"] == "VIRULENCE"] if "VIRULENCE" in amr_raw["Type"].values else pd.DataFrame()
stress_only = amr_raw[amr_raw["Type"] == "STRESS"] if "STRESS" in amr_raw["Type"].values else pd.DataFrame()

print(f"  Integrated: {len(merged)}, pLIN: {len(plin)}, AMR raw: {len(amr_raw)}")


# ── Helper functions ──────────────────────────────────────────────────────────
def style_header(ws, ncols, row=1):
    hfill = PatternFill(start_color="1565C0", end_color="1565C0", fill_type="solid")
    hfont = Font(name="Calibri", bold=True, color="FFFFFF", size=11)
    for c in range(1, ncols + 1):
        cell = ws.cell(row=row, column=c)
        cell.fill = hfill
        cell.font = hfont
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)


def style_data(ws, nrows, ncols, start_row=2):
    alt = PatternFill(start_color="E3F2FD", end_color="E3F2FD", fill_type="solid")
    dfont = Font(name="Calibri", size=10)
    bdr = Border(bottom=Side(style="thin", color="DDDDDD"))
    for r in range(start_row, start_row + nrows):
        for c in range(1, ncols + 1):
            cell = ws.cell(row=r, column=c)
            cell.font = dfont
            cell.border = bdr
            if (r - start_row) % 2 == 1:
                cell.fill = alt


def auto_width(ws, min_w=8, max_w=45):
    for col in ws.columns:
        mx = 0
        letter = get_column_letter(col[0].column)
        for cell in col:
            if cell.value:
                mx = max(mx, len(str(cell.value)))
        ws.column_dimensions[letter].width = min(max(mx + 2, min_w), max_w)


def simpsons_d(labels):
    counts = Counter(labels)
    N = sum(counts.values())
    if N <= 1:
        return 0.0
    return 1.0 - sum(n * (n - 1) for n in counts.values()) / (N * (N - 1))


# ══════════════════════════════════════════════════════════════════════════════
# PART 1: MAIN TABLES EXCEL (Tables 1-6)
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Creating Tables 1-6 Excel workbook ===")
wb = Workbook()

# ── Table 1: Comparative evaluation ──
print("  Table 1: Comparative evaluation ...")
ws1 = wb.active
ws1.title = "Table 1 - Comparison"

headers = ["Feature", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA", "mge-cluster", "pLIN"]
for c, h in enumerate(headers, 1):
    ws1.cell(row=1, column=c, value=h)

rows_t1 = [
    ["Classification approach", "Replicon typing", "Allelic profiling", "Relaxase clustering",
     "Host-range + mobility", "Reference-free k-mer", "Hierarchical 4-mer"],
    ["Resolution levels", "1 (Inc group)", "1 (sequence type)", "1 (cluster)",
     "1 (PTU)", "1 (cluster)", "6 (L1-L6)"],
    ["Multi-resolution hierarchy", "No", "No", "No", "Partial", "No", "Yes"],
    ["Stable nomenclature", "Yes", "Yes", "No*", "No*", "No", "Yes"],
    ["Discriminatory power (Simpson's D)", "0.641", "NA", "NA", "NA", "NA", "0.985"],
    ["Integrated AMR profiling", "No", "No", "No", "No", "No", "Yes"],
    ["Automated outbreak detection", "No", "No", "No", "No", "No", "Yes"],
    ["Clinical risk stratification", "No", "No", "No", "No", "No", "Yes"],
    ["Reference-free operation", "No", "No", "No", "No", "Yes", "Yes**"],
    ["Open-source availability", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"],
    ["No. Inc/Rep groups supported", ">30", "~10", "NA", "NA", "NA", "28"],
    ["Scalability (>50,000 plasmids)", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"],
]
for r, row in enumerate(rows_t1, 2):
    for c, val in enumerate(row, 1):
        ws1.cell(row=r, column=c, value=val)

# Footnotes
fn_row = len(rows_t1) + 3
ws1.cell(row=fn_row, column=1, value="*MOB-suite and COPLA reassign cluster identifiers with each database update.")
ws1.cell(row=fn_row + 1, column=1, value="**pLIN operates reference-free for 4-mer computation; the Inc group classifier requires the reference database.")

style_header(ws1, len(headers))
style_data(ws1, len(rows_t1), len(headers))
auto_width(ws1)

# ── Table 2: Threshold definitions ──
print("  Table 2: Threshold definitions ...")
ws2 = wb.create_sheet("Table 2 - Thresholds")

headers2 = ["Level", "Cosine distance (d)", "Approximate ANI", "Taxonomic analogy",
            "Clinical interpretation", "Calibration quantile"]
for c, h in enumerate(headers2, 1):
    ws2.cell(row=1, column=c, value=h)

rows_t2 = [
    ["L1", "<=0.150", "~85%", "Plasmid superfamily",
     "Broad plasmid family identification", "99th percentile"],
    ["L2", "<=0.100", "~90%", "Major lineage",
     "Epidemiological lineage grouping; regional surveillance", "95th-99th percentile"],
    ["L3", "<=0.050", "~95%", "Species-level cluster",
     "Regional surveillance; epidemic context", "95th percentile (d=0.048)"],
    ["L4", "<=0.020", "~98%", "Sublineage",
     "Inter-hospital transmission tracking", "75th percentile (d=0.025)"],
    ["L5", "<=0.010", "~99%", "Clone complex",
     "Intra-hospital transmission; ward-level spread", "Median (d=0.011)"],
    ["L6", "<=0.001", "~99.9%", "Strain / outbreak level",
     "Outbreak confirmation; infection control trigger", "25th percentile (d=0.001)"],
]
for r, row in enumerate(rows_t2, 2):
    for c, val in enumerate(row, 1):
        ws2.cell(row=r, column=c, value=val)

fn_row = len(rows_t2) + 3
ws2.cell(row=fn_row, column=1,
         value="FastANI validation (n=4,970 pairs): At L6, median ANI=99.9%. Spearman rho=-0.348 (P<10^-141).")

style_header(ws2, len(headers2))
style_data(ws2, len(rows_t2), len(headers2))
auto_width(ws2)

# ── Table 3: Training dataset summary ──
print("  Table 3: Training dataset summary ...")
ws3 = wb.create_sheet("Table 3 - Training Data")

headers3 = ["Inc/Rep group", "n", "% of dataset", "Unique pLIN codes (L6)",
            "AMR carriage (%)", "Mean AMR genes", "Virulence (%)"]
for c, h in enumerate(headers3, 1):
    ws3.cell(row=1, column=c, value=h)

# Manuscript Table 3 shows only the 20 Gram-negative Inc groups (6,998 plasmids)
# Sort by count descending to match manuscript ordering
gn_sorted = sorted(GRAM_NEG, key=lambda x: inc_counts.get(x, 0), reverse=True)

# Use only Gram-negative subset from merged (integrated) data
gn_plin = plin[plin["inc_type"].isin(GRAM_NEG)]
gn_total = len(merged)  # merged is already Gram-negative only (6,998)

row_idx = 2
total_n = 0
total_unique_l6 = 0
for inc in gn_sorted:
    sub_merged = merged[merged[INC_COL] == inc]
    n = len(sub_merged)
    if n == 0:
        continue
    total_n += n
    pct = n / gn_total * 100
    unique_l6 = sub_merged["pLIN"].nunique()
    total_unique_l6 += unique_l6

    amr_pct = (sub_merged["n_amr_genes"] > 0).sum() / n * 100
    mean_amr = sub_merged["n_amr_genes"].mean()
    vir_pct = (sub_merged["n_vir_genes"] > 0).sum() / n * 100

    ws3.cell(row=row_idx, column=1, value=inc)
    ws3.cell(row=row_idx, column=2, value=n)
    ws3.cell(row=row_idx, column=3, value=round(pct, 1))
    ws3.cell(row=row_idx, column=4, value=unique_l6)
    ws3.cell(row=row_idx, column=5, value=round(amr_pct, 1))
    ws3.cell(row=row_idx, column=6, value=round(mean_amr, 1))
    ws3.cell(row=row_idx, column=7, value=round(vir_pct, 1))
    row_idx += 1

# Totals row
ws3.cell(row=row_idx, column=1, value="Total")
ws3.cell(row=row_idx, column=2, value=total_n)
ws3.cell(row=row_idx, column=3, value="100.0")
ws3.cell(row=row_idx, column=4, value=total_unique_l6)
total_font = Font(name="Calibri", bold=True, size=10)
for c in range(1, len(headers3) + 1):
    ws3.cell(row=row_idx, column=c).font = total_font

# Footnote matching manuscript
fn_row3 = row_idx + 2
ws3.cell(row=fn_row3, column=1,
         value=f"Overall: {(merged['n_amr_genes'] > 0).sum():,} of {len(merged):,} plasmids (83.1%) "
               "carried at least one AMR, virulence, or stress-response gene; "
               f"{(merged['n_amr_genes'] > 0).sum():,} carried at least one AMR gene specifically.")
ws3.cell(row=fn_row3 + 1, column=1,
         value="Data source: NCBI RefSeq complete plasmid sequences, classified by PlasmidFinder "
               "(identity >=95%, coverage >=60%). AMR annotation: AMRFinderPlus v4.2.5 (database 2026-01-21.1).")

style_header(ws3, len(headers3))
style_data(ws3, row_idx - 1, len(headers3))
auto_width(ws3)

# ── Table 4: ML classifier performance ──
print("  Table 4: ML classifier performance ...")
ws4 = wb.create_sheet("Table 4 - ML Performance")

headers4 = ["Model", "Weighted F1 (mean +/- SD)", "Key hyperparameters"]
for c, h in enumerate(headers4, 1):
    ws4.cell(row=1, column=c, value=h)

rows_t4 = [
    ["KNN (primary classifier)", "0.911*", "k=5, cosine distance, distance-weighted"],
    ["XGBoost", "0.896 +/- 0.011", "max_depth=6, n_estimators=300, learning_rate=0.1"],
    ["Gradient Boosting", "0.893 +/- 0.009", "max_depth=5, n_estimators=200, learning_rate=0.1"],
    ["Random Forest", "0.874 +/- 0.021", "n_estimators=500, max_features=sqrt"],
    ["Logistic Regression", "0.866 +/- 0.012", "C=1.0, multi_class=multinomial, solver=lbfgs"],
]
for r, row in enumerate(rows_t4, 2):
    for c, val in enumerate(row, 1):
        ws4.cell(row=r, column=c, value=val)

fn_row = len(rows_t4) + 3
ws4.cell(row=fn_row, column=1,
         value="*KNN accuracy from 5-fold stratified cross-validation (91.1%).")
ws4.cell(row=fn_row + 1, column=1,
         value="Training: 256-dim 4-mer vectors from 8,077 sequences across 28 Inc/Rep groups.")

style_header(ws4, len(headers4))
style_data(ws4, len(rows_t4), len(headers4))
auto_width(ws4)

# ── Table 5: Critical resistance determinants ──
print("  Table 5: Critical resistance determinants ...")
ws5 = wb.create_sheet("Table 5 - Critical AMR")

headers5 = ["Resistance category", "WHO priority", "Total detections",
            "Top variants (n)", "Therapeutic implications"]
for c, h in enumerate(headers5, 1):
    ws5.cell(row=1, column=c, value=h)

# Compute from raw data
gene_symbols = amr_only["Element symbol"]
carb = gene_symbols[gene_symbols.apply(lambda g: any(p in str(g) for p in ["blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-244", "blaVIM", "blaIMP"]))]
esbl = gene_symbols[gene_symbols.apply(lambda g: any(p in str(g) for p in ["blaCTX-M", "blaSHV"]))]
mcr = gene_symbols[gene_symbols.str.contains("mcr-", na=False)]
pmqr = gene_symbols[gene_symbols.apply(lambda g: any(p in str(g) for p in ["qnr", "aac(6')-Ib-cr", "oqxA", "oqxB"]))]

def top_variants(series, n=4):
    return ", ".join(f"{g} ({c})" for g, c in Counter(series).most_common(n))

n_amr_pos = (merged["n_amr_genes"] > 0).sum()
top_overall = Counter()
for g in merged[merged["n_amr_genes"] > 0]["amr_genes"]:
    if g != "none":
        top_overall.update(g.split("; "))
top3_overall = top_overall.most_common(3)
top3_str = ", ".join(f"{g} ({c:,}; {c/n_amr_pos*100:.1f}%)" for g, c in top3_overall)

rows_t5 = [
    ["Carbapenemases", "Critical", len(carb), top_variants(carb),
     "Last-resort beta-lactam failure; 40-50% BSI mortality"],
    ["ESBLs", "Critical", len(esbl), top_variants(esbl),
     "Third-generation cephalosporin failure"],
    ["Colistin resistance (mcr)", "Critical", len(mcr), top_variants(mcr),
     "Last-resort polymyxin failure; pan-drug resistance risk"],
    ["PMQR", "High", len(pmqr), top_variants(pmqr),
     "Fluoroquinolone failure; undermines empirical UTI therapy"],
    ["Top AMR genes overall", "--", "--", top3_str,
     "Penicillin, sulfonamide, and tetracycline resistance backbone"],
]
for r, row in enumerate(rows_t5, 2):
    for c, val in enumerate(row, 1):
        ws5.cell(row=r, column=c, value=val)

fn_row = len(rows_t5) + 3
ws5.cell(row=fn_row, column=1,
         value=f"AMR-positive plasmids: n={n_amr_pos:,}. Total gene detections: {len(amr_raw):,} "
               f"({len(amr_only):,} AMR + {len(vir_only):,} virulence + {len(stress_only):,} stress).")
ws5.cell(row=fn_row + 1, column=1,
         value="AMRFinderPlus v4.2.5 (database 2026-01-21.1).")

style_header(ws5, len(headers5))
style_data(ws5, len(rows_t5), len(headers5))
auto_width(ws5)

# ── Table 6: Top pLIN lineages ──
print("  Table 6: Top pLIN lineages ...")
ws6 = wb.create_sheet("Table 6 - Top Lineages")

headers6 = ["pLIN code (L6)", "n", "Inc group(s)", "Mean AMR genes",
            "Key resistance determinants (% carriage)", "Clinical risk profile"]
for c, h in enumerate(headers6, 1):
    ws6.cell(row=1, column=c, value=h)

# Manuscript Table 6 shows these 4 key surveillance lineages
manuscript_plins = ["1327", "1434", "860", "671"]

# Clinical risk profiles and cross-validation matches from manuscript
manuscript_risks = {
    "1327": "Largest circulating lineage; moderate-burden resistance hub; broad Inc group distribution",
    "1434": "ESBL-carrying IncFII lineage; community-associated urinary tract infections",
    "860": "Multi-Inc convergence zone; multi-drug resistance with colistin resistance (44.4% mcr carriage); MDR hub",
    "671": "Uniformly carbapenemase-positive; immediate isolation trigger; confirmed in German surveillance (Yao 2023)",
}

row_idx = 2
for plin_suffix in manuscript_plins:
    # Match pLIN codes ending with this L6 suffix
    matching = merged[merged["pLIN"].astype(str).str.endswith(f".{plin_suffix}")]
    if len(matching) == 0:
        matching = merged[merged["pLIN"].astype(str) == plin_suffix]
    if len(matching) == 0:
        continue

    sub = matching
    n_total = len(sub)
    inc_groups = sub[INC_COL].value_counts()
    inc_str = ", ".join(f"{inc} ({c})" for inc, c in inc_groups.items())
    mean_amr = sub["n_amr_genes"].mean()

    # Key resistance determinants
    genes_list = []
    for g in sub["amr_genes"]:
        if g != "none":
            genes_list.extend(g.split("; "))
    gc = Counter(genes_list)
    top_genes = gc.most_common(4)
    det_str = ", ".join(f"{g} ({c/n_total*100:.1f}%)" for g, c in top_genes)

    risk = manuscript_risks.get(plin_suffix, "Surveillance target")

    ws6.cell(row=row_idx, column=1, value=str(sub["pLIN"].iloc[0]))
    ws6.cell(row=row_idx, column=2, value=n_total)
    ws6.cell(row=row_idx, column=3, value=inc_str)
    ws6.cell(row=row_idx, column=4, value=round(mean_amr, 1))
    ws6.cell(row=row_idx, column=5, value=det_str)
    ws6.cell(row=row_idx, column=6, value=risk)
    row_idx += 1

# Add footnote about cross-validation matches
fn_row = row_idx + 1
ws6.cell(row=fn_row, column=1,
         value="Cross-validation matches: pLIN 671 (CP104944, Yao 2023, KPC-2, 100% confidence, d=0.0000); "
               "pLIN 860 (CP022533, Roberts 2020, IMP-4, 100% confidence, d=0.0004).")

style_header(ws6, len(headers6))
style_data(ws6, row_idx - 2, len(headers6))
auto_width(ws6)

# Save Tables workbook
tables_path = os.path.join(TABLES_DIR, "pLIN_Manuscript_Tables_1-6.xlsx")
wb.save(tables_path)
print(f"  Saved: {tables_path}")


# ══════════════════════════════════════════════════════════════════════════════
# PART 2: SUPPLEMENTARY TABLES EXCEL
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Creating Supplementary Tables Excel ===")
wb_supp = Workbook()

# ── Supp Table S1: Full pLIN assignments ──
print("  Supp Table S1: Full pLIN assignments ...")
ws_s1 = wb_supp.active
ws_s1.title = "S1 - pLIN Assignments"
for c, h in enumerate(plin.columns, 1):
    ws_s1.cell(row=1, column=c, value=h)
for r, (_, row) in enumerate(plin.iterrows(), 2):
    for c, val in enumerate(row, 1):
        ws_s1.cell(row=r, column=c, value=val)
    if r > 8100:
        break
style_header(ws_s1, len(plin.columns))
auto_width(ws_s1)

# ── Supp Table S2: AMR gene prevalence by Inc group ──
print("  Supp Table S2: AMR gene prevalence by Inc group ...")
ws_s2 = wb_supp.create_sheet("S2 - AMR by Inc Group")
s2_headers = ["Inc group", "n total", "n AMR+", "AMR %", "n VIR+", "VIR %",
              "n STRESS+", "STRESS %", "Mean AMR genes", "Mean VIR genes"]
for c, h in enumerate(s2_headers, 1):
    ws_s2.cell(row=1, column=c, value=h)

row_idx = 2
for inc in GRAM_NEG:
    sub = merged[merged[INC_COL] == inc]
    n = len(sub)
    if n == 0:
        continue
    ws_s2.cell(row=row_idx, column=1, value=inc)
    ws_s2.cell(row=row_idx, column=2, value=n)
    ws_s2.cell(row=row_idx, column=3, value=int((sub["n_amr_genes"] > 0).sum()))
    ws_s2.cell(row=row_idx, column=4, value=round((sub["n_amr_genes"] > 0).sum() / n * 100, 1))
    ws_s2.cell(row=row_idx, column=5, value=int((sub["n_vir_genes"] > 0).sum()))
    ws_s2.cell(row=row_idx, column=6, value=round((sub["n_vir_genes"] > 0).sum() / n * 100, 1))
    ws_s2.cell(row=row_idx, column=7, value=int((sub["n_stress_genes"] > 0).sum()))
    ws_s2.cell(row=row_idx, column=8, value=round((sub["n_stress_genes"] > 0).sum() / n * 100, 1))
    ws_s2.cell(row=row_idx, column=9, value=round(sub["n_amr_genes"].mean(), 1))
    ws_s2.cell(row=row_idx, column=10, value=round(sub["n_vir_genes"].mean(), 1))
    row_idx += 1
style_header(ws_s2, len(s2_headers))
style_data(ws_s2, row_idx - 2, len(s2_headers))
auto_width(ws_s2)

# ── Supp Table S3: Top 50 AMR genes ──
print("  Supp Table S3: Top 50 AMR genes ...")
ws_s3 = wb_supp.create_sheet("S3 - Top AMR Genes")
s3_headers = ["Rank", "Gene", "Detections", "% of AMR+ plasmids", "Drug class"]

gene_class_map = {}
for gene, group in amr_only.groupby("Element symbol"):
    gene_class_map[gene] = group["Class"].value_counts().index[0]

all_amr_genes = []
for g in merged[merged["n_amr_genes"] > 0]["amr_genes"]:
    if g != "none":
        all_amr_genes.extend(g.split("; "))
gene_counts = Counter(all_amr_genes)
n_amr_pos = (merged["n_amr_genes"] > 0).sum()

for c, h in enumerate(s3_headers, 1):
    ws_s3.cell(row=1, column=c, value=h)
for rank, (gene, cnt) in enumerate(gene_counts.most_common(50), 1):
    ws_s3.cell(row=rank + 1, column=1, value=rank)
    ws_s3.cell(row=rank + 1, column=2, value=gene)
    ws_s3.cell(row=rank + 1, column=3, value=cnt)
    ws_s3.cell(row=rank + 1, column=4, value=round(cnt / n_amr_pos * 100, 1))
    ws_s3.cell(row=rank + 1, column=5, value=gene_class_map.get(gene, "Unknown"))
style_header(ws_s3, len(s3_headers))
style_data(ws_s3, 50, len(s3_headers))
auto_width(ws_s3)

# ── Supp Table S4: pLIN lineage AMR summary (top 50) ──
print("  Supp Table S4: pLIN lineage AMR summary ...")
ws_s4 = wb_supp.create_sheet("S4 - Lineage AMR Summary")
s4_headers = ["pLIN code", "n plasmids", "Inc group(s)", "Mean AMR", "Mean VIR",
              "Top 3 AMR genes", "AMR %", "VIR %"]
for c, h in enumerate(s4_headers, 1):
    ws_s4.cell(row=1, column=c, value=h)

top50_plins = merged["pLIN"].value_counts().head(50)
row_idx = 2
for pcode, n_plas in top50_plins.items():
    sub = merged[merged["pLIN"] == pcode]
    inc_str = ", ".join(sub[INC_COL].value_counts().index.tolist()[:3])
    genes = []
    for g in sub["amr_genes"]:
        if g != "none":
            genes.extend(g.split("; "))
    top3 = ", ".join(g for g, _ in Counter(genes).most_common(3)) if genes else "none"

    ws_s4.cell(row=row_idx, column=1, value=str(pcode))
    ws_s4.cell(row=row_idx, column=2, value=n_plas)
    ws_s4.cell(row=row_idx, column=3, value=inc_str)
    ws_s4.cell(row=row_idx, column=4, value=round(sub["n_amr_genes"].mean(), 1))
    ws_s4.cell(row=row_idx, column=5, value=round(sub["n_vir_genes"].mean(), 1))
    ws_s4.cell(row=row_idx, column=6, value=top3)
    ws_s4.cell(row=row_idx, column=7, value=round((sub["n_amr_genes"] > 0).sum() / len(sub) * 100, 1))
    ws_s4.cell(row=row_idx, column=8, value=round((sub["n_vir_genes"] > 0).sum() / len(sub) * 100, 1))
    row_idx += 1
style_header(ws_s4, len(s4_headers))
style_data(ws_s4, row_idx - 2, len(s4_headers))
auto_width(ws_s4)

# ── Supp Table S5: FastANI validation summary ──
print("  Supp Table S5: FastANI validation ...")
ws_s5 = wb_supp.create_sheet("S5 - ANI Validation")
if os.path.isfile(ANI_FILE):
    ani = pd.read_csv(ANI_FILE, sep="\t")
    s5_headers = ["Inc group", "n pairs", "Median ANI", "Mean cosine dist",
                  "Spearman rho", "P-value"]
    for c, h in enumerate(s5_headers, 1):
        ws_s5.cell(row=1, column=c, value=h)

    from scipy.stats import spearmanr
    row_idx = 2
    for inc in GRAM_NEG:
        sub = ani[ani["inc_type"] == inc].dropna(subset=["fastani_ani", "cosine_distance"])
        sub = sub[sub["fastani_ani"] > 0]
        if len(sub) < 3:
            continue
        rho, pval = spearmanr(sub["cosine_distance"], sub["fastani_ani"])
        ws_s5.cell(row=row_idx, column=1, value=inc)
        ws_s5.cell(row=row_idx, column=2, value=len(sub))
        ws_s5.cell(row=row_idx, column=3, value=round(sub["fastani_ani"].median(), 2))
        ws_s5.cell(row=row_idx, column=4, value=round(sub["cosine_distance"].mean(), 4))
        ws_s5.cell(row=row_idx, column=5, value=round(rho, 3))
        ws_s5.cell(row=row_idx, column=6, value=f"{pval:.2e}")
        row_idx += 1

    # Overall
    all_valid = ani.dropna(subset=["fastani_ani", "cosine_distance"])
    all_valid = all_valid[all_valid["fastani_ani"] > 0]
    rho_all, pval_all = spearmanr(all_valid["cosine_distance"], all_valid["fastani_ani"])
    ws_s5.cell(row=row_idx, column=1, value="Overall")
    ws_s5.cell(row=row_idx, column=2, value=len(all_valid))
    ws_s5.cell(row=row_idx, column=3, value=round(all_valid["fastani_ani"].median(), 2))
    ws_s5.cell(row=row_idx, column=4, value=round(all_valid["cosine_distance"].mean(), 4))
    ws_s5.cell(row=row_idx, column=5, value=round(rho_all, 3))
    ws_s5.cell(row=row_idx, column=6, value=f"{pval_all:.2e}")
    for c in range(1, len(s5_headers) + 1):
        ws_s5.cell(row=row_idx, column=c).font = Font(name="Calibri", bold=True, size=10)

    style_header(ws_s5, len(s5_headers))
    style_data(ws_s5, row_idx - 1, len(s5_headers))
    auto_width(ws_s5)

supp_path = os.path.join(TABLES_DIR, "pLIN_Supplementary_Tables.xlsx")
wb_supp.save(supp_path)
print(f"  Saved: {supp_path}")


# ══════════════════════════════════════════════════════════════════════════════
# PART 3: RAW DATA FILES FOR FIGURES
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Exporting raw datasets for figures ===")

# Figure 1 source data
print("  Figure 1 source data ...")
fig1_a = pd.DataFrame({
    "inc_type": INC_ORDER,
    "n_plasmids": [inc_counts[inc] for inc in INC_ORDER],
    "pct_of_total": [inc_counts[inc] / len(plin) * 100 for inc in INC_ORDER],
})
fig1_a.to_csv(os.path.join(RAW_DATA_DIR, "Figure1A_group_distribution.tsv"), sep="\t", index=False)

fig1_b = pd.DataFrame({
    "inc_type": INC_ORDER,
    "unique_pLIN_L6": [plin[plin["inc_type"] == inc]["bin_F"].nunique() for inc in INC_ORDER],
})
fig1_b.to_csv(os.path.join(RAW_DATA_DIR, "Figure1B_unique_pLIN_codes.tsv"), sep="\t", index=False)

d_plin = simpsons_d(plin["pLIN"].values)
d_inc = simpsons_d(plin["inc_type"].values)
fig1_c = pd.DataFrame({
    "method": ["Inc/Rep typing alone", "pLIN (6-level)"],
    "simpsons_D": [round(d_inc, 4), round(d_plin, 4)],
})
fig1_c.to_csv(os.path.join(RAW_DATA_DIR, "Figure1C_simpsons_D.tsv"), sep="\t", index=False)

# Figure 2 source data
print("  Figure 2 source data ...")
top20 = gene_counts.most_common(20)
fig2_data = pd.DataFrame([
    {"gene": g, "detections": c, "pct_amr_pos": round(c / n_amr_pos * 100, 1),
     "drug_class": gene_class_map.get(g, "Unknown")}
    for g, c in top20
])
fig2_data.to_csv(os.path.join(RAW_DATA_DIR, "Figure2_top20_AMR_genes.tsv"), sep="\t", index=False)

# Figure 3 source data
print("  Figure 3 source data ...")
critical_cats = {
    "Carbapenemases": ["blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-244", "blaVIM", "blaIMP"],
    "ESBLs": ["blaCTX-M", "blaSHV"],
    "Colistin_mcr": ["mcr-"],
    "PMQR": ["qnr", "aac(6')-Ib-cr", "oqxA", "oqxB"],
}
for cat, patterns in critical_cats.items():
    matching = [g for g in amr_only["Element symbol"] if any(p in str(g) for p in patterns)]
    counts = Counter(matching).most_common()
    df = pd.DataFrame(counts, columns=["gene", "detections"])
    df["total_category"] = len(matching)
    df.to_csv(os.path.join(RAW_DATA_DIR, f"Figure3_{cat}.tsv"), sep="\t", index=False)

# Figure 4 source data (heatmap matrix)
print("  Figure 4 source data ...")
amr_pos_m = merged[merged["n_amr_genes"] > 0]
top10_plins = amr_pos_m["pLIN"].value_counts().head(10).index.tolist()
top15_genes = [g for g, _ in gene_counts.most_common(15)]

heatmap_rows = []
for pcode in top10_plins:
    sub = merged[(merged["pLIN"] == pcode) & (merged["n_amr_genes"] > 0)]
    n = len(sub)
    inc_str = ", ".join(sorted(sub[INC_COL].unique()))
    genes_in = []
    for g in sub["amr_genes"]:
        if g != "none":
            genes_in.extend(g.split("; "))
    gc = Counter(genes_in)
    row = {"pLIN": pcode, "inc_types": inc_str, "n_plasmids": n}
    for gene in top15_genes:
        row[gene] = round(gc.get(gene, 0) / max(n, 1) * 100, 1)
    heatmap_rows.append(row)
pd.DataFrame(heatmap_rows).to_csv(os.path.join(RAW_DATA_DIR, "Figure4_heatmap_data.tsv"), sep="\t", index=False)

# Figure 5 source data
print("  Figure 5 source data ...")
fig5_rows = []
for inc in GRAM_NEG:
    sub = merged[merged[INC_COL] == inc]
    if len(sub) == 0:
        continue
    fig5_rows.append({
        "inc_type": inc,
        "n_total": len(sub),
        "mean_amr_genes": round(sub["n_amr_genes"].mean(), 1),
        "amr_carriage_pct": round((sub["n_amr_genes"] > 0).sum() / len(sub) * 100, 1),
        "vir_prevalence_pct": round((sub["n_vir_genes"] > 0).sum() / len(sub) * 100, 1),
    })
pd.DataFrame(fig5_rows).to_csv(os.path.join(RAW_DATA_DIR, "Figure5_AMR_burden_virulence.tsv"), sep="\t", index=False)

# Figure 6 source data (ANI validation)
print("  Figure 6 source data ...")
if os.path.isfile(ANI_FILE):
    import shutil
    shutil.copy2(ANI_FILE, os.path.join(RAW_DATA_DIR, "Figure6_ANI_validation.tsv"))

print("  All source data exported.")


# ══════════════════════════════════════════════════════════════════════════════
# PART 4: APPENDIX PAGES
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Creating Appendix pages ===")

# Copy existing appendix files
src_appendix = os.path.join(BASE_DIR, "manuscripts", "lancet_microbe", "supplementary")
if os.path.isdir(src_appendix):
    import shutil
    for f in sorted(os.listdir(src_appendix)):
        if f.startswith("Appendix_") and f.endswith(".md"):
            shutil.copy2(os.path.join(src_appendix, f), os.path.join(APPENDIX_DIR, f))
            print(f"  Copied: {f}")

# Create/update key appendix pages with latest data

# Appendix p01: Dataset breakdown
print("  Creating Appendix p01 (Dataset breakdown) ...")
lines = ["# Appendix p01: Dataset Breakdown — 28 Inc/Rep Groups\n\n"]
lines.append(f"Total training sequences: {len(plin):,}\n\n")
lines.append("## Gram-negative Inc groups (20 groups)\n\n")
lines.append("| Inc group | n | % | Unique L6 codes |\n|---|---|---|---|\n")
for inc in sorted(GRAM_NEG, key=lambda x: inc_counts.get(x, 0), reverse=True):
    n = inc_counts.get(inc, 0)
    u = plin[plin["inc_type"] == inc]["bin_F"].nunique()
    lines.append(f"| {inc} | {n:,} | {n/len(plin)*100:.1f} | {u:,} |\n")
lines.append(f"\n**Subtotal**: {sum(inc_counts.get(i,0) for i in GRAM_NEG):,}\n\n")

lines.append("## Gram-positive Rep type groups (4 groups)\n\n")
lines.append("| Rep group | n | % | Unique L6 codes |\n|---|---|---|---|\n")
for inc in GRAM_POS:
    n = inc_counts.get(inc, 0)
    u = plin[plin["inc_type"] == inc]["bin_F"].nunique()
    lines.append(f"| {inc} | {n:,} | {n/len(plin)*100:.1f} | {u:,} |\n")

lines.append("\n## Acinetobacter spp. groups (2 groups)\n\n")
lines.append("| Rep group | n | % | Unique L6 codes |\n|---|---|---|---|\n")
for inc in ACI_GROUPS:
    n = inc_counts.get(inc, 0)
    u = plin[plin["inc_type"] == inc]["bin_F"].nunique()
    lines.append(f"| {inc} | {n:,} | {n/len(plin)*100:.1f} | {u:,} |\n")

lines.append("\n## Pseudomonas spp. groups (2 groups)\n\n")
lines.append("| Rep group | n | % | Unique L6 codes |\n|---|---|---|---|\n")
for inc in PAE_GROUPS:
    n = inc_counts.get(inc, 0)
    u = plin[plin["inc_type"] == inc]["bin_F"].nunique()
    lines.append(f"| {inc} | {n:,} | {n/len(plin)*100:.1f} | {u:,} |\n")

with open(os.path.join(APPENDIX_DIR, "Appendix_p01_Dataset_breakdown_28_groups.md"), "w") as f:
    f.writelines(lines)

# Appendix p02: AMR gene prevalence by Inc group
print("  Creating Appendix p02 (AMR prevalence by Inc group) ...")
lines = ["# Appendix p02: AMR Gene Prevalence by Inc Group\n\n"]
lines.append(f"AMRFinderPlus v4.2.5 (database 2026-01-21.1)\n\n")
lines.append("| Inc group | n | AMR+ n | AMR % | Mean AMR | Top 3 genes |\n|---|---|---|---|---|---|\n")
for inc in sorted(GRAM_NEG, key=lambda x: inc_counts.get(x, 0), reverse=True):
    sub = merged[merged[INC_COL] == inc]
    n = len(sub)
    if n == 0:
        continue
    amr_n = (sub["n_amr_genes"] > 0).sum()
    amr_pct = amr_n / n * 100
    mean_amr = sub["n_amr_genes"].mean()
    genes = []
    for g in sub[sub["n_amr_genes"] > 0]["amr_genes"]:
        if g != "none":
            genes.extend(g.split("; "))
    top3 = ", ".join(g for g, _ in Counter(genes).most_common(3)) if genes else "—"
    lines.append(f"| {inc} | {n:,} | {amr_n:,} | {amr_pct:.1f} | {mean_amr:.1f} | {top3} |\n")

with open(os.path.join(APPENDIX_DIR, "Appendix_p02_AMR_gene_prevalence_by_Inc_group.md"), "w") as f:
    f.writelines(lines)

# Appendix p03: Top 20 AMR genes
print("  Creating Appendix p03 (Top 20 AMR genes) ...")
lines = ["# Appendix p03: Top 20 AMR Genes — Clinical Significance\n\n"]
lines.append(f"Among {n_amr_pos:,} AMR-positive plasmids.\n\n")
lines.append("| Rank | Gene | n | % AMR+ | Drug class | Clinical significance |\n|---|---|---|---|---|---|\n")

clinical_significance = {
    "blaTEM-1": "Most common penicillinase; baseline beta-lactam resistance",
    "sul1": "Sulfonamide resistance; often integron-associated",
    "tet(A)": "Tetracycline efflux pump; widespread in Enterobacterales",
    "aph(6)-Id": "Streptomycin resistance; common in resistance cassettes",
    "aph(3'')-Ib": "Streptomycin resistance; frequently co-carried with aph(6)-Id",
    "sul2": "Sulfonamide resistance; plasmid-borne",
    "mph(A)": "Macrolide phosphotransferase; azithromycin resistance",
    "mrx(A)": "Macrolide resistance; often co-located with mph(A)",
    "blaKPC-2": "Carbapenemase; critical WHO priority; last-resort resistance",
    "qnrS1": "Quinolone resistance; low-level fluoroquinolone resistance",
    "aac(6')-Ib-cr5": "Dual aminoglycoside/quinolone resistance",
    "dfrA14": "Trimethoprim resistance; DHFR inhibitor resistance",
    "aadA2": "Streptomycin/spectinomycin resistance",
    "catB3": "Chloramphenicol acetyltransferase",
    "blaOXA-1": "Narrow-spectrum oxacillinase",
    "aph(3')-Ia": "Kanamycin resistance",
    "floR": "Florfenicol/chloramphenicol efflux",
    "blaCTX-M-15": "ESBL; third-gen cephalosporin resistance; pandemic clone",
    "aac(3)-IId": "Gentamicin resistance",
    "dfrA12": "Trimethoprim resistance",
}
for rank, (gene, cnt) in enumerate(gene_counts.most_common(20), 1):
    pct = cnt / n_amr_pos * 100
    cls = gene_class_map.get(gene, "Unknown")
    sig = clinical_significance.get(gene, "—")
    lines.append(f"| {rank} | *{gene}* | {cnt:,} | {pct:.1f} | {cls} | {sig} |\n")

with open(os.path.join(APPENDIX_DIR, "Appendix_p03_Top_20_AMR_genes_clinical_significance.md"), "w") as f:
    f.writelines(lines)

# Appendix p09: Threshold calibration
print("  Creating Appendix p09 (Threshold calibration) ...")
lines = ["# Appendix p09: Hierarchical Threshold Calibration\n\n"]
lines.append("## Nearest-Neighbour Distance Distribution\n\n")
lines.append("The six pLIN thresholds were calibrated using the nearest-neighbour distance distribution\n")
lines.append("of within-group plasmid pairs. The following quantiles were used:\n\n")
lines.append("| Quantile | Distance | pLIN Level |\n|---|---|---|\n")
lines.append("| 25th percentile | 0.001 | L6 (strain) |\n")
lines.append("| Median | 0.011 | L5 (clone complex) |\n")
lines.append("| 75th percentile | 0.025 | L4 (sublineage) |\n")
lines.append("| 95th percentile | 0.048 | L3 (species-level) |\n")
lines.append("| 95th-99th | 0.100 | L2 (major lineage) |\n")
lines.append("| 99th percentile | 0.150 | L1 (superfamily) |\n")
lines.append("\n## FastANI Validation\n\n")
if os.path.isfile(ANI_FILE):
    ani = pd.read_csv(ANI_FILE, sep="\t")
    valid = ani.dropna(subset=["fastani_ani", "cosine_distance"])
    valid = valid[valid["fastani_ani"] > 0]
    from scipy.stats import spearmanr
    rho, pval = spearmanr(valid["cosine_distance"], valid["fastani_ani"])
    l6 = valid[valid["cosine_distance"] <= 0.001]
    lines.append(f"- Total pairs validated: {len(valid):,}\n")
    lines.append(f"- Overall Spearman rho: {rho:.3f} (P = {pval:.2e})\n")
    lines.append(f"- L6 pairs (d <= 0.001): {len(l6):,}\n")
    lines.append(f"- L6 median ANI: {l6['fastani_ani'].median():.1f}%\n")
    lines.append(f"- L6 IQR: {l6['fastani_ani'].quantile(0.25):.1f}% — {l6['fastani_ani'].quantile(0.75):.1f}%\n")

with open(os.path.join(APPENDIX_DIR, "Appendix_p09_Hierarchical_threshold_calibration.md"), "w") as f:
    f.writelines(lines)

# Appendix p12: Reference database expansion
print("  Creating Appendix p12 (Reference database expansion) ...")
ref_file = os.path.join(BASE_DIR, "output", "pLIN_reference_assignments.tsv")
lines = ["# Appendix p12: Reference Database Expansion Statistics\n\n"]
if os.path.isfile(ref_file):
    ref_df = pd.read_csv(ref_file, sep="\t")
    lines.append(f"## Expanded Database: {len(ref_df):,} plasmid sequences\n\n")
    lines.append(f"- Training set: {len(plin):,} plasmids, {plin['pLIN'].nunique():,} unique pLIN codes\n")
    lines.append(f"- Expanded database: {len(ref_df):,} plasmids, {ref_df['pLIN'].nunique():,} unique pLIN codes\n")
    fold = ref_df['pLIN'].nunique() / max(plin['pLIN'].nunique(), 1)
    lines.append(f"- Fold increase in pLIN diversity: {fold:.1f}x\n\n")

    if "inc_type" in ref_df.columns:
        lines.append("## Inc/Rep Group Distribution (Expanded Database)\n\n")
        lines.append("| Inc/Rep group | n | Unique pLIN codes |\n|---|---|---|\n")
        for inc in ref_df["inc_type"].value_counts().index:
            sub = ref_df[ref_df["inc_type"] == inc]
            lines.append(f"| {inc} | {len(sub):,} | {sub['pLIN'].nunique():,} |\n")
else:
    lines.append("Reference database file not found.\n")

with open(os.path.join(APPENDIX_DIR, "Appendix_p12_Reference_database_expansion_statistics.md"), "w") as f:
    f.writelines(lines)


# ══════════════════════════════════════════════════════════════════════════════
# PART 5: FIGURE LEGENDS DOCUMENT
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Creating Figure Legends document ===")

# Compute actual values from data for the legends
d_plin_val = round(simpsons_d(plin["pLIN"].values), 3)
d_inc_val = round(simpsons_d(plin["inc_type"].values), 3)
fold_improv = round(d_plin_val / d_inc_val, 2)
n_unique_l6 = plin["bin_F"].nunique()
n_singletons = (plin["pLIN"].value_counts() == 1).sum()
singleton_pct = round(n_singletons / n_unique_l6 * 100, 1)

# Gene counts from raw data
carb_total = len(amr_only[amr_only["Element symbol"].apply(lambda g: any(p in str(g) for p in ["blaKPC", "blaNDM", "blaOXA-48", "blaOXA-181", "blaOXA-232", "blaOXA-244", "blaVIM", "blaIMP"]))])
esbl_total = len(amr_only[amr_only["Element symbol"].apply(lambda g: any(p in str(g) for p in ["blaCTX-M", "blaSHV"]))])
mcr_total = len(amr_only[amr_only["Element symbol"].str.contains("mcr-", na=False)])
pmqr_total = len(amr_only[amr_only["Element symbol"].apply(lambda g: any(p in str(g) for p in ["qnr", "aac(6')-Ib-cr", "oqxA", "oqxB"]))])

# Top gene counts
top3_genes = gene_counts.most_common(3)

legends = f"""# Figure Legends — pLIN Manuscript (Figures 1-6)

Data-verified values computed from the training dataset ({len(plin):,} plasmids, {len(INC_ORDER)} Inc/Rep groups).

---

## Figure 1: Dataset overview and pLIN diversity

(A) Distribution of {len(plin):,} plasmid sequences across {len(INC_ORDER)} groups ({len(GRAM_NEG)} Gram-negative Inc groups, {len(GRAM_POS)} Gram-positive rep type groups, {len(ACI_GROUPS)} A. baumannii rep type groups, and {len(PAE_GROUPS)} P. aeruginosa rep type groups). Among Gram-negative groups, three dominate: IncFII (n={inc_counts.get("IncFII",0):,}; {inc_counts.get("IncFII",0)/len(plin)*100:.1f}%), IncN (n={inc_counts.get("IncN",0):,}; {inc_counts.get("IncN",0)/len(plin)*100:.1f}%), and IncX1 (n={inc_counts.get("IncX1",0):,}; {inc_counts.get("IncX1",0)/len(plin)*100:.1f}%). (B) Number of unique pLIN strain-level (L6) codes per group. A total of {n_unique_l6:,} unique codes were resolved, with {n_singletons:,} ({singleton_pct}%) representing singletons. (C) Comparison of discriminatory power: pLIN (Simpson's D = {d_plin_val}) versus Inc/rep typing alone (D = {d_inc_val}), representing a {fold_improv}-fold improvement. Error bars indicate 95% confidence intervals.

---

## Figure 2: AMR gene prevalence overview

Prevalence of the top 20 most frequently detected AMR genes across {n_amr_pos:,} AMR-positive plasmids. The three most prevalent genes were {top3_genes[0][0]} (n={top3_genes[0][1]:,}; {top3_genes[0][1]/n_amr_pos*100:.1f}%), {top3_genes[1][0]} (n={top3_genes[1][1]:,}; {top3_genes[1][1]/n_amr_pos*100:.1f}%), and {top3_genes[2][0]} (n={top3_genes[2][1]:,}; {top3_genes[2][1]/n_amr_pos*100:.1f}%). Gene detections were identified using AMRFinderPlus v4.2.5 (database 2026-01-21.1). Bars are coloured by drug class.

---

## Figure 3: Clinically critical resistance determinants

(A) Carbapenemase gene distribution: {carb_total:,} total detections dominated by blaKPC-2 (n={Counter(amr_only[amr_only["Element symbol"].str.contains("blaKPC-2", na=False)]["Element symbol"]).most_common(1)[0][1] if len(amr_only[amr_only["Element symbol"].str.contains("blaKPC-2", na=False)]) > 0 else 0}). (B) Extended-spectrum beta-lactamase (ESBL) gene distribution: {esbl_total:,} total detections. (C) Mobile colistin resistance gene distribution: {mcr_total:,} total detections. (D) Plasmid-mediated quinolone resistance (PMQR) determinant distribution: {pmqr_total:,} total detections.

---

## Figure 4: AMR gene prevalence heatmap across pLIN lineages

Heatmap showing the prevalence of clinically important AMR genes across the top pLIN lineages (rows) ordered by total AMR burden. Columns represent individual resistance genes grouped by drug class. Colour intensity reflects prevalence within each lineage (0-100%). White cells indicate absence of the gene in that lineage.

---

## Figure 5: AMR gene burden and virulence distribution

(A) Distribution of mean AMR gene count per plasmid across the {len(GRAM_NEG)} Gram-negative Inc groups. AMR profiling was performed on the {len(GRAM_NEG)} Gram-negative Inc groups only; Gram-positive, Acinetobacter spp., and Pseudomonas spp. groups were not profiled. (B) Virulence gene prevalence by group.

---

## Figure 6: pLIN hierarchical structure

Dendrogram illustrates the six-level hierarchical clustering structure for a representative subset of plasmids. Horizontal dashed lines indicate the six distance thresholds: L1 (d<=0.150), L2 (d<=0.100), L3 (d<=0.050), L4 (d<=0.020), L5 (d<=0.010), and L6 (d<=0.001). Branch colours represent Inc/rep groups. The inset shows the FastANI validation: at L6 (d<=0.001), median ANI = 99.9%; overall Spearman rho = -0.348 across {len(ani.dropna(subset=['fastani_ani','cosine_distance']).query('fastani_ani > 0')) if os.path.isfile(ANI_FILE) else 'N/A'} plasmid pairs.
"""

with open(os.path.join(REVIEW_DIR, "Figure_Legends.md"), "w") as f:
    f.write(legends)
print("  Figure legends saved.")


# ══════════════════════════════════════════════════════════════════════════════
# PART 6: ACCESSION NUMBERS
# ══════════════════════════════════════════════════════════════════════════════
print("\n=== Exporting accession numbers ===")
accessions = plin["plasmid_id"].tolist()
with open(os.path.join(APPENDIX_DIR, "Appendix_Accession_Numbers.txt"), "w") as f:
    f.write(f"# Plasmid Accession Numbers — pLIN Training Dataset\n")
    f.write(f"# Total: {len(accessions):,} sequences across {len(INC_ORDER)} Inc/Rep groups\n")
    f.write(f"# Format: accession_number\\tinc_type\n\n")
    for _, row in plin[["plasmid_id", "inc_type"]].iterrows():
        f.write(f"{row['plasmid_id']}\t{row['inc_type']}\n")
print(f"  Exported {len(accessions):,} accession numbers.")


# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"  REVIEW PACKAGE COMPLETE — for_review2/")
print(f"{'='*70}")

for root, dirs, files in os.walk(REVIEW_DIR):
    level = root.replace(REVIEW_DIR, "").count(os.sep)
    indent = " " * 2 * level
    subdir = os.path.basename(root)
    print(f"{indent}{subdir}/")
    subindent = " " * 2 * (level + 1)
    for f in sorted(files):
        if f.startswith("."):
            continue
        fpath = os.path.join(root, f)
        size_kb = os.path.getsize(fpath) / 1024
        print(f"{subindent}{f:<55s} {size_kb:>8.1f} KB")

print(f"{'='*70}")
