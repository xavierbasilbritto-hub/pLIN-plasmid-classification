#!/usr/bin/env python3
"""Generate Excel tables for all four journal manuscripts."""

import os
try:
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter
except ImportError:
    print("openpyxl not found, installing...")
    import subprocess
    subprocess.check_call(["pip", "install", "openpyxl"])
    import openpyxl
    from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
    from openpyxl.utils import get_column_letter

BASE = "/Users/basilxavier/Desktop/PLASMID_TOOL/manuscripts"

# Style definitions
HEADER_FONT = Font(name="Arial", bold=True, size=10, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="2F5496", end_color="2F5496", fill_type="solid")
TITLE_FONT = Font(name="Arial", bold=True, size=12)
BODY_FONT = Font(name="Arial", size=10)
NOTE_FONT = Font(name="Arial", size=9, italic=True)
THIN_BORDER = Border(
    left=Side(style="thin", color="D9D9D9"),
    right=Side(style="thin", color="D9D9D9"),
    top=Side(style="thin", color="D9D9D9"),
    bottom=Side(style="thin", color="D9D9D9"),
)
ALT_FILL = PatternFill(start_color="D6E4F0", end_color="D6E4F0", fill_type="solid")
WRAP = Alignment(wrap_text=True, vertical="top")
CENTER = Alignment(horizontal="center", vertical="top", wrap_text=True)


def style_sheet(ws, title, headers, data, notes=None, col_widths=None):
    """Add a styled table to a worksheet."""
    # Title row
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(headers))
    cell = ws.cell(row=1, column=1, value=title)
    cell.font = TITLE_FONT
    cell.alignment = Alignment(wrap_text=True)

    # Header row
    for j, h in enumerate(headers, 1):
        cell = ws.cell(row=3, column=j, value=h)
        cell.font = HEADER_FONT
        cell.fill = HEADER_FILL
        cell.alignment = CENTER
        cell.border = THIN_BORDER

    # Data rows
    for i, row in enumerate(data, 4):
        for j, val in enumerate(row, 1):
            cell = ws.cell(row=i, column=j, value=val)
            cell.font = BODY_FONT
            cell.alignment = WRAP
            cell.border = THIN_BORDER
            if (i - 4) % 2 == 1:
                cell.fill = ALT_FILL

    # Notes
    if notes:
        note_row = 4 + len(data) + 1
        ws.merge_cells(start_row=note_row, start_column=1, end_row=note_row, end_column=len(headers))
        cell = ws.cell(row=note_row, column=1, value=notes)
        cell.font = NOTE_FONT
        cell.alignment = Alignment(wrap_text=True)

    # Column widths
    if col_widths:
        for j, w in enumerate(col_widths, 1):
            ws.column_dimensions[get_column_letter(j)].width = w
    else:
        for j in range(1, len(headers) + 1):
            ws.column_dimensions[get_column_letter(j)].width = 18


# ==============================================================================
# LANCET MICROBE
# ==============================================================================
def create_lancet_microbe():
    wb = openpyxl.Workbook()

    # Table 1: Comparison
    ws = wb.active
    ws.title = "Table 1 - Comparison"
    style_sheet(ws, "Table 1: Comparative evaluation of plasmid classification systems",
        ["Feature", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA", "mge-cluster", "pLIN"],
        [
            ["Classification approach", "Replicon typing", "Allelic profiling", "Relaxase clustering", "Host-range + mobility", "Reference-free k-mer", "Hierarchical 4-mer"],
            ["Resolution levels", "1 (Inc group)", "1 (sequence type)", "1 (cluster)", "1 (PTU)", "1 (cluster)", "6 (L1-L6)"],
            ["Multi-resolution hierarchy", "No", "No", "No", "Partial", "No", "Yes"],
            ["Stable nomenclature", "Yes", "Yes", "No*", "No*", "No", "Yes"],
            ["Discriminatory power (Simpson's D)", "0.641", "NA", "NA", "NA", "NA", "0.985"],
            ["Integrated AMR profiling", "No", "No", "No", "No", "No", "Yes"],
            ["Automated outbreak detection", "No", "No", "No", "No", "No", "Yes"],
            ["Clinical risk stratification", "No", "No", "No", "No", "No", "Yes"],
            ["Reference-free operation", "No", "No", "No", "No", "Yes", "Yes**"],
            ["Open-source availability", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"],
            ["No. Inc groups supported", ">30", "~10", "NA", "NA", "NA", "28 (20 Gram-neg + 4 Gram-pos + 2 Acinetobacter + 2 Pseudomonas)"],
            ["Scalability (>50,000 plasmids)", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"],
        ],
        notes="* MOB-suite and COPLA reassign cluster identifiers with each database update. ** pLIN operates reference-free for 4-mer computation; the Inc classifier requires the reference database.",
        col_widths=[30, 16, 16, 16, 16, 16, 16])

    # Table 2: Thresholds
    ws2 = wb.create_sheet("Table 2 - Thresholds")
    style_sheet(ws2, "Table 2: pLIN hierarchical threshold definitions",
        ["Level", "Cosine distance (d)", "Approx. ANI", "Taxonomic analogy", "Clinical interpretation", "Calibration quantile"],
        [
            ["L1", "≤0.150", "~85%", "Plasmid superfamily", "Broad plasmid family identification", "99th percentile"],
            ["L2", "≤0.100", "~90%", "Major lineage", "Epidemiological lineage grouping", "95th-99th percentile"],
            ["L3", "≤0.050", "~95%", "Species-level cluster", "Regional surveillance comparisons", "95th percentile (d=0.048)"],
            ["L4", "≤0.020", "~98%", "Sublineage", "Inter-hospital transmission tracking", "75th percentile (d=0.025)"],
            ["L5", "≤0.010", "~99%", "Clone group", "Intra-hospital transmission investigation", "Median (d=0.011)"],
            ["L6", "≤0.001", "~99.9%", "Lineage / outbreak", "Outbreak confirmation; infection control trigger", "25th percentile (d=0.001)"],
        ],
        notes="FastANI validation (n=4,970 pairs, 20 Gram-negative Inc groups): At L6 (d≤0.001), median ANI = 99.9%. Spearman ρ = -0.348 (P < 10⁻¹⁴¹). Significant correlations for 15/20 Inc groups (P < 0.05).",
        col_widths=[10, 20, 14, 22, 35, 28])

    # Table 3: Dataset
    ws3 = wb.create_sheet("Table 3 - Dataset")
    style_sheet(ws3, "Table 3: Training dataset — 8,077 plasmids across 28 Inc/Rep groups",
        ["Inc/Rep group", "n", "% of dataset", "Unique pLIN codes (L6)", "AMR carriage (%)", "Mean AMR genes", "Virulence (%)"],
        [
            ["IncFII", 4629, 57.3, 1421, 62.3, 4.2, 23.2],
            ["IncN", 1097, 13.6, 431, 82.5, 6.0, 3.8],
            ["IncX1", 705, 8.7, 420, 71.4, 5.1, 8.6],
            ["repAci_large", 203, 2.5, 112, "—", "—", "—"],
            ["repPae_large", 199, 2.5, 108, "—", "—", "—"],
            ["repPae_small", 150, 1.9, 82, "—", "—", "—"],
            ["repSA_large", 121, 1.5, 78, "—", "—", "—"],
            ["repEF_res", "100*", 1.2, 69, "—", "—", "—"],
            ["repAci1", 120, 1.5, 68, "—", "—", "—"],
            ["IncFIB", 97, 1.2, 42, 68.0, 4.8, 15.5],
            ["repEF_conj", 92, 1.1, 52, "—", "—", "—"],
            ["ColRNAI", 91, 1.1, 38, 45.1, 1.8, 2.2],
            ["IncF", 75, 0.9, 31, 72.0, 5.3, 18.7],
            ["repSA_small", 73, 0.9, 41, "—", "—", "—"],
            ["IncX3", 56, 0.7, 24, 87.5, 3.9, 1.8],
            ["IncHI2", 36, 0.4, 18, 100.0, 10.8, 5.6],
            ["IncI1", 27, 0.3, 14, 85.2, 4.7, 3.7],
            ["IncI2", 25, 0.3, 13, 80.0, 3.2, 4.0],
            ["IncX4", 24, 0.3, 12, 66.7, 2.5, 0.0],
            ["IncR", 21, 0.3, 11, 85.7, 7.1, 4.8],
            ["ColE", 19, 0.2, 9, 52.6, 1.4, 0.0],
            ["IncC", 16, 0.2, 8, 100.0, 9.4, 6.3],
            ["IncHI1", 16, 0.2, 8, 93.8, 8.7, 6.3],
            ["IncFIC", 14, 0.2, 7, 100.0, 6.9, 21.4],
            ["IncAC2", 14, 0.2, 7, 100.0, 12.6, 7.1],
            ["IncA", 14, 0.2, 7, 100.0, 8.9, 7.1],
            ["IncI", 11, 0.1, 6, 81.8, 4.5, 0.0],
            ["IncFIBK", 11, 0.1, 6, 63.6, 3.1, 9.1],
            ["Total", "8,056†", 100.0, 3073, "83.1*", "—", "—"],
        ],
        notes="* AMR run on 20 GN groups only: 5,816/6,998 (83.1%) carried ≥1 AMR/virulence/stress gene; 4,657 (66.5%) carried ≥1 AMR gene. GP/Aci/Pae AMR not yet performed (marked —). † 8,077 training samples includes 21 multi-replicon duplicates; 8,056 unique plasmids. * repEF_res: 100 unique + 21 shared with repEF_conj = 121 training samples.",
        col_widths=[14, 8, 12, 22, 16, 16, 14])

    # Table 4: ML Performance
    ws4 = wb.create_sheet("Table 4 - ML Performance")
    style_sheet(ws4, "Table 4: Machine learning classifier performance for Inc group prediction",
        ["Model", "Weighted F1 (mean ± SD)", "Key hyperparameters"],
        [
            ["KNN (primary classifier)", "0.911*", "k=5, cosine distance, distance-weighted"],
            ["XGBoost", "0.896 ± 0.011", "max_depth=6, n_estimators=300, learning_rate=0.1"],
            ["Gradient Boosting", "0.893 ± 0.009", "max_depth=5, n_estimators=200, learning_rate=0.1"],
            ["Random Forest", "0.874 ± 0.021", "n_estimators=500, max_features=sqrt"],
            ["Logistic Regression", "0.866 ± 0.012", "C=1.0, multi_class=multinomial, solver=lbfgs"],
        ],
        notes="* KNN overall accuracy from 5-fold stratified CV (91.1%) across 28 Inc/Rep groups (8,077 samples). Dominant confusion: IncFIB↔IncFII (d=0.002).",
        col_widths=[25, 25, 50])

    # Table 5: Critical AMR
    ws5 = wb.create_sheet("Table 5 - Critical AMR")
    style_sheet(ws5, "Table 5: Clinically critical resistance determinants (n=8,077 plasmids)",
        ["Resistance category", "WHO priority", "Total detections", "Top variants (n)", "Therapeutic implications"],
        [
            ["Carbapenemases", "Critical", 1635, "blaKPC-2 (824), blaNDM-1 (228), blaKPC-3 (193), blaNDM-5 (89)", "Last-resort beta-lactam failure; 40-50% BSI mortality"],
            ["ESBLs", "Critical", 1804, "blaCTX-M-15 (505), blaCTX-M-65 (319), blaSHV-12 (277)", "3rd-gen cephalosporin failure"],
            ["Colistin resistance (mcr)", "Critical", 204, "mcr-1.1 (83), mcr-8.1 (27), mcr-8.2 (18)", "Last-resort polymyxin failure; pan-drug resistance risk"],
            ["PMQR", "High", 2315, "qnrS1 (732), aac(6')-Ib-cr5 (737), qnrB1 (212)", "Fluoroquinolone failure"],
            ["Top AMR genes overall", "—", "—", "blaTEM-1 (1,864; 40.0%), sul1 (1,383; 29.7%), tet(A) (1,297; 27.9%)", "Penicillin, sulfonamide, tetracycline resistance backbone"],
        ],
        notes="AMRFinderPlus v4.2.5 (database 2026-01-21.1). Total: 64,891 detections (29,583 AMR + 6,286 virulence + 29,022 stress) in 5,816/6,998 Gram-negative plasmids (83.1%). 100% AMR groups: IncA, IncAC2, IncC, IncFIC, IncHI2. Highest mean burden: IncAC2 (12.6 genes/plasmid).",
        col_widths=[24, 14, 16, 50, 40])

    # Table 6: Lineages
    ws6 = wb.create_sheet("Table 6 - Lineages")
    style_sheet(ws6, "Table 6: Top pLIN lineages — clinical surveillance targets",
        ["pLIN code (L6)", "n", "Inc group(s)", "Mean AMR genes", "Key resistance determinants (% carriage)", "Clinical risk profile"],
        [
            ["1327", 869, "6 Inc groups: IncFII (814), IncFIB (38), IncFIBK (7), IncN (5), IncHI1 (4), IncF (1)", 5.5, "blaTEM-1 (33.4%), sul1 (30.4%), mrx(A) (25.3%)", "Largest lineage; moderate-burden hub"],
            ["1434", 163, "IncFII (92.0%), IncF (8.0%)", 7.4, "mph(A) (56.4%), mrx(A) (56.4%), sul1 (56.4%), blaCTX-M-27 (29.4%)", "ESBL-carrying IncFII; community UTI"],
            ["860", 142, "5 Inc groups: IncN (73.2%), IncHI2 (23.2%)", 14.4, "sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), mcr-1.1 (36.6%)", "Cross-Inc MDR hub; colistin + carbapenem resistance"],
            ["671", 90, "IncN", 13.2, "blaKPC-2 (100%), blaTEM-1 (100%), aph(3'')-Ib (98%)", "Uniformly KPC-2+; immediate isolation trigger"],
        ],
        notes="Cross-validation: pLIN 671 matched CP104944 (Yao 2023, Germany, d=0.0000); pLIN 860 matched CP022533 (Roberts 2020, Australia, d=0.0004).",
        col_widths=[16, 8, 30, 16, 45, 40])

    add_appendix_sheets(wb)
    path = os.path.join(BASE, "lancet_microbe", "tables", "Lancet_Microbe_Tables.xlsx")
    wb.save(path)
    print(f"  Saved: {path}")


# ==============================================================================
# BRIEFINGS IN BIOINFORMATICS
# ==============================================================================
def create_briefings():
    wb = openpyxl.Workbook()

    # Table 1: Comparison
    ws = wb.active
    ws.title = "Table 1 - Comparison"
    style_sheet(ws, "Table 1: Comparative evaluation of plasmid classification systems",
        ["Property", "pLIN", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA/PTUs", "mge-cluster"],
        [
            ["Year introduced", "2025", "2014", "2014", "2018", "2021", "2023"],
            ["Classification basis", "Whole-sequence 4-mer (256 features)", "Replicon gene detection", "Allelic variants", "Relaxase + Mash distance", "ANI network + HSBM", "Unitig Jaccard + HDBSCAN"],
            ["Hierarchical levels", "6 nested (L1-L6)", "1 (flat)", "1 (flat)", "1 (flat)", "2-3 (semi)", "1 (flat)"],
            ["Code permanence", "Yes (NN rule)", "Stable (DB-dependent)", "Stable (DB-dependent)", "No (re-clustering)", "No (recomputation)", "No (re-embedding)"],
            ["Reference-free", "Yes (core)", "No", "No", "No", "No", "Yes"],
            ["Taxonomic scope", "28 Inc/Rep groups; universal", "~30 replicons", "6 Inc schemes", "Broad (DB-limited)", "41% assignable", "Broad"],
            ["Resolution", "Lineage (~99.9% ANI)", "Family", "Sub-family (ST)", "Species (Mash 0.06)", "Species", "Variable"],
            ["Simpson's D (n=8,077)", "0.985", "0.641", "N/A", "N/A", "N/A", "N/A"],
            ["Computational cost", "Low (minutes)", "Low (seconds)", "Low (seconds)", "Moderate", "High (hours)", "Moderate"],
            ["ML validation", "XGBoost F1=0.896", "None", "None", "None", "None", "None"],
            ["AMR integration", "Yes (AMRFinderPlus)", "No", "No", "No", "No", "No"],
            ["Outbreak detection", "Yes (basic + temporal)", "No", "No", "No", "No", "No"],
        ],
        col_widths=[22, 28, 22, 22, 22, 22, 22])

    # Table 2: Thresholds
    ws2 = wb.create_sheet("Table 2 - Thresholds")
    style_sheet(ws2, "Table 2: pLIN threshold definitions and clustering results",
        ["Level", "Bin", "Distance (d)", "ANI equiv.", "Interpretation", "Clusters (training)", "Singletons", "Clusters (79,305)"],
        [
            ["L1", "A", "≤0.150", "~85%", "Family", 1, 0, 253],
            ["L2", "B", "≤0.100", "~90%", "Subfamily", 1, 0, 662],
            ["L3", "C", "≤0.050", "~95%", "Cluster", 7, 5, 5975],
            ["L4", "D", "≤0.020", "~98%", "Subcluster", 38, 27, 23651],
            ["L5", "E", "≤0.010", "~99%", "Clone group", 117, 82, 36100],
            ["L6", "F", "≤0.001", "~99.9%", "Lineage/outbreak", 3073, 2335, 57886],
        ],
        notes="FastANI validation: n=4,970 pairs (20 GN groups), Spearman ρ=-0.348, P<10⁻¹⁴¹. At d≤0.001: median ANI=99.9%. Strongest: IncFIC (ρ=-0.88), IncI2 (-0.81), ColE (-0.72).",
        col_widths=[10, 8, 14, 12, 18, 20, 14, 20])

    # Table 3: Benchmarking
    ws3 = wb.create_sheet("Table 3 - Benchmarking")
    style_sheet(ws3, "Table 3: Comprehensive benchmarking results",
        ["Metric", "Value", "Details"],
        [
            ["Simpson's D (pLIN L6)", "0.985", "3,073 unique codes from 8,077 plasmids (28 groups)"],
            ["Simpson's D (Inc/Rep typing)", "0.641", "28 Inc/Rep groups"],
            ["Improvement", "1.54-fold", "pLIN L6 vs Inc/Rep typing"],
            ["Inc concordance", "97.8%", "3,006/3,073 single-Inc L6 codes; 67 mixed-Inc"],
            ["Singletons", "76.0%", "2,335/3,073 L6 codes"],
            ["KNN accuracy", "91.1%", "k=5, cosine, distance-weighted, 5-fold CV, 28 groups"],
            ["XGBoost F1", "0.896 ± 0.011", "Nested CV: 3 outer, 2 inner folds"],
            ["Gradient Boosting F1", "0.893 ± 0.009", "Optuna TPE, Hyperband pruning"],
            ["Random Forest F1", "0.874 ± 0.021", "10 trials/fold"],
            ["Logistic Regression F1", "0.866 ± 0.012", "C=1.0, multinomial"],
            ["Classification rate (79,305)", "97.3%", "2,109 below 40% confidence threshold (2.7%)"],
            ["Total runtime", "<30 min", "Apple M-series, single-threaded"],
            ["Phase 1 (4-mer)", "23.5 min", "~52 seq/sec, 71,249 sequences"],
            ["Phase 2 (KNN)", "7 sec", "~10,272 seq/sec"],
            ["Phase 3 (clustering)", "4.2 min", "28 groups"],
            ["Max memory (IncFII)", "2.3 GB", "34,036 sequences"],
            ["Global matrix (avoided)", "~24 GB", "10.4× reduction"],
            ["Unique codes (79,305)", "57,886", "19.8× increase over training"],
            ["AMR detections", "64,891", "29,583 AMR + 6,286 virulence + 29,022 stress"],
            ["Plasmids with hits", "5,816/8,077 (83.1%)", "4,657 (66.5%) with AMR specifically"],
            ["Carbapenemases", "1,635", "blaKPC-2 (824), blaNDM-1 (228)"],
            ["ESBLs", "1,804", "blaCTX-M-15 (505), blaCTX-M-65 (319)"],
            ["mcr (colistin)", "204", "mcr-1.1 (83), mcr-8.1 (27)"],
            ["PMQR", "2,315", "qnrS1 (732), aac(6')-Ib-cr5 (737)"],
            ["Outbreak plasmids tested", "74", "26 studies, 13 countries, 7 resistance mechanisms"],
            ["High-risk lineage matches", "3", "pLIN 671 (KPC-2), pLIN 860 (MDR hub), pLIN 1688 (OXA-48)"],
        ],
        col_widths=[30, 22, 50])

    add_appendix_sheets(wb)
    path = os.path.join(BASE, "briefings_bioinformatics", "tables", "Briefings_Bioinformatics_Tables.xlsx")
    wb.save(path)
    print(f"  Saved: {path}")


# ==============================================================================
# MICROBIAL GENOMICS
# ==============================================================================
def create_microbial_genomics():
    wb = openpyxl.Workbook()

    # Table 1: Comparison
    ws = wb.active
    ws.title = "Table 1 - Comparison"
    style_sheet(ws, "Table 1: Comparative evaluation of plasmid typing systems",
        ["Property", "pLIN", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA", "mge-cluster"],
        [
            ["Year", "2025", "2014", "2014", "2018", "2020", "2023"],
            ["Basis", "4-mer composition + cosine", "Replicon BLAST", "Allelic variants", "Relaxase + Mash", "ANI + HSBM", "Unitig + HDBSCAN"],
            ["Hierarchical levels", "6 (L1-L6)", "1", "1", "1", "2-3", "1"],
            ["Code permanence", "Yes", "Yes (DB-dep.)", "Yes (DB-dep.)", "No", "No", "No"],
            ["Reference-free", "Yes", "No", "No", "No", "No", "Yes"],
            ["Scope", "28 Inc/Rep groups", "~30 replicons", "6 Inc schemes", "Broad", "41% assignable", "Broad"],
            ["Resolution", "Lineage (~99.9% ANI)", "Family", "Sub-family", "Species", "Species", "Variable"],
            ["Simpson's D", "0.985", "0.641", "N/A", "N/A", "N/A", "N/A"],
            ["Classification rate", "97.3%", "~50%", "Limited", "~95%", "41-63%", "High"],
            ["Computational cost", "<30 min (79,305)", "Seconds", "Seconds", "Minutes", "Hours", "Minutes"],
            ["AMR integration", "Yes", "No", "No", "No", "No", "No"],
            ["Outbreak detection", "Yes", "No", "No", "No", "No", "No"],
        ],
        col_widths=[22, 25, 20, 20, 20, 20, 20])

    # Table 2: Thresholds
    ws2 = wb.create_sheet("Table 2 - Thresholds")
    style_sheet(ws2, "Table 2: pLIN hierarchical thresholds and clustering results (n=8,077)",
        ["Level", "Designation", "Distance (d)", "ANI equiv.", "Clusters", "Max cluster", "Singletons"],
        [
            ["L1", "Family", "≤0.150", "~85%", 1, 8077, 0],
            ["L2", "Subfamily", "≤0.100", "~90%", 1, 8077, 0],
            ["L3", "Cluster", "≤0.050", "~95%", 7, 6990, 5],
            ["L4", "Subcluster", "≤0.020", "~98%", 38, 6947, 27],
            ["L5", "Clone group", "≤0.010", "~99%", 117, 6737, 82],
            ["L6", "Lineage", "≤0.001", "~99.9%", 3073, 869, 2335],
        ],
        notes="FastANI (v1.34, n=4,970 pairs, 20 GN groups): Spearman ρ=-0.348 (P<10⁻¹⁴¹). At d≤0.001: median ANI=99.9%.",
        col_widths=[10, 16, 14, 12, 12, 14, 14])

    # Table 3: Dataset
    ws3 = wb.create_sheet("Table 3 - Dataset")
    style_sheet(ws3, "Table 3: Dataset composition — 8,077 plasmids across 28 Inc/Rep groups",
        ["Inc group", "Sequences (n)", "% of dataset", "Unique pLIN codes (L6)"],
        [
            ["IncFII", 4629, 66.1, 1421],
            ["IncN", 1097, 15.7, 431],
            ["IncX1", 705, 10.1, 420],
            ["IncFIB", 97, 1.4, 42],
            ["ColRNAI", 91, 1.3, 38],
            ["IncF", 75, 1.1, 31],
            ["IncX3", 56, 0.8, 24],
            ["IncHI2", 36, 0.5, 18],
            ["IncI1", 27, 0.4, 14],
            ["IncI2", 25, 0.4, 13],
            ["IncX4", 24, 0.3, 12],
            ["IncR", 21, 0.3, 11],
            ["ColE", 19, 0.3, 9],
            ["IncC", 16, 0.2, 8],
            ["IncHI1", 16, 0.2, 8],
            ["IncFIC", 14, 0.2, 7],
            ["IncAC2", 14, 0.2, 7],
            ["IncA", 14, 0.2, 7],
            ["IncI", 11, 0.2, 6],
            ["IncFIBK", 11, 0.2, 6],
            ["Total", 8077, 100.0, 3073],
        ],
        notes="Simpson's D at L6 = 0.985 (28 groups); Inc/Rep typing alone D = 0.641. Note: GN groups shown; 28 total groups include 4 GP, 2 Aci, 2 Pae.",
        col_widths=[14, 16, 14, 24])

    # Table 4: ML
    ws4 = wb.create_sheet("Table 4 - ML Validation")
    style_sheet(ws4, "Table 4: Machine learning model performance (nested CV, balanced 3-class subset)",
        ["Model", "Mean weighted F1", "SD", "Min", "Max"],
        [
            ["XGBoost", 0.896, 0.011, 0.885, 0.907],
            ["Gradient Boosting", 0.893, 0.009, 0.884, 0.902],
            ["Random Forest", 0.874, 0.021, 0.853, 0.895],
            ["Logistic Regression", 0.866, 0.012, 0.854, 0.878],
        ],
        notes="1,500 plasmids (500 each: IncFII, IncN, IncX1). 3 outer folds, 2 inner folds. Optuna TPE (10 trials/fold, Hyperband). KNN primary classifier: 91.1% accuracy (5-fold CV).",
        col_widths=[22, 20, 10, 10, 10])

    # Table 5: AMR
    ws5 = wb.create_sheet("Table 5 - AMR")
    style_sheet(ws5, "Table 5: AMR gene detection (AMRFinderPlus v4.2.5, n=8,077)",
        ["Category", "Total detections", "Top variants (n)"],
        [
            ["Overall AMR genes", 29583, "blaTEM-1 (1,864; 40.0%), sul1 (1,383; 29.7%), tet(A) (1,297; 27.9%)"],
            ["Virulence factors", 6286, "traT, spvB, iucA (IncFII-enriched)"],
            ["Stress response", 29022, "Various stress tolerance genes"],
            ["Total detections", 64891, "In 5,816/6,998 plasmids (83.1%); 4,657 (66.5%) with AMR"],
            ["Carbapenemases", 1635, "blaKPC-2 (824), blaNDM-1 (228), blaKPC-3 (193), blaNDM-5 (89)"],
            ["ESBLs", 1804, "blaCTX-M-15 (505), blaCTX-M-65 (319), blaSHV-12 (277)"],
            ["Colistin (mcr)", 204, "mcr-1.1 (83), mcr-8.1 (27), mcr-8.2 (18)"],
            ["PMQR", 2315, "qnrS1 (732), aac(6')-Ib-cr5 (737), qnrB1 (212)"],
        ],
        notes="100% AMR carriage: IncA, IncAC2, IncC, IncFIC, IncHI2. Highest mean burden: IncAC2 (12.6 genes/plasmid).",
        col_widths=[22, 18, 65])

    # Table 6: Lineages + Cross-validation
    ws6 = wb.create_sheet("Table 6 - Lineages")
    style_sheet(ws6, "Table 6: High-risk pLIN lineages and outbreak cross-validation",
        ["pLIN (L6)", "Members", "AMR+", "Mean AMR", "Inc group(s)", "Key resistance genes (% carriage)"],
        [
            ["1327", 869, 700, 5.5, "6 Inc groups: IncFII (814), IncFIB (38), +4", "blaTEM-1 (33.4%), sul1 (30.4%), mrx(A) (25.3%)"],
            ["860", 142, 142, 14.4, "IncN (73%), IncHI2 (23%), +3", "sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), mcr-1.1 (36.6%)"],
            ["671", 90, 90, 13.2, "IncN", "blaKPC-2 (100%), blaTEM-1 (100%), aph(3'')-Ib (98%)"],
            ["1434", 163, 161, 7.4, "IncFII (92.0%), IncF (8.0%)", "mph(A) (56.4%), mrx(A) (56.4%), sul1 (56.4%), blaCTX-M-27 (29.4%)"],
        ],
        notes="Cross-validation: 74 plasmids from 26 studies across 13 countries; 85.1% high-confidence; 3 globally disseminated lineages matched (pLIN 671/KPC-2, 860/MDR hub, 1688/OXA-48); 9 intra-study clusters detected.",
        col_widths=[14, 12, 10, 12, 28, 50])

    add_appendix_sheets(wb)
    path = os.path.join(BASE, "microbial_genomics", "tables", "Microbial_Genomics_Tables.xlsx")
    wb.save(path)
    print(f"  Saved: {path}")


# ==============================================================================
# CMI
# ==============================================================================
def create_cmi():
    wb = openpyxl.Workbook()

    # Table 1: Comparison
    ws = wb.active
    ws.title = "Table 1 - Comparison"
    style_sheet(ws, "Table 1: Comparison of plasmid typing approaches for clinical infection control",
        ["Feature", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA", "pLIN"],
        [
            ["Typing basis", "Replicon markers", "Allelic profiles", "Relaxase + replicon", "Network clustering", "4-mer composition (LIN)"],
            ["Resolution levels", "1 (Inc group)", "1 (sequence type)", "1 (MOB cluster)", "1 (PTU)", "6 (L1-L6)"],
            ["Simpson's D", "0.641", "N/A", "N/A", "N/A", "0.985"],
            ["Stable nomenclature", "Yes", "Yes", "No", "No", "Yes"],
            ["AMR profiling", "No", "No", "No", "No", "Yes (AMRFinderPlus)"],
            ["Lineage-level AMR profiles", "No", "No", "No", "No", "Yes"],
            ["Outbreak detection", "No", "No", "No", "No", "Yes (2-tier)"],
            ["Risk stratification", "No", "No", "No", "No", "Yes (CRITICAL/HIGH/MOD)"],
            ["SNP sub-typing", "No", "No", "No", "No", "Yes"],
            ["Multi-resolution", "No", "No", "No", "Partial", "Yes (family to strain)"],
            ["Reference-free", "No", "No", "No", "No", "Yes"],
            ["Standard hardware (<30 min)", "Yes", "Yes", "Yes", "Yes", "Yes"],
        ],
        col_widths=[28, 18, 18, 18, 18, 28])

    # Table 2: Thresholds
    ws2 = wb.create_sheet("Table 2 - Thresholds")
    style_sheet(ws2, "Table 2: pLIN thresholds, ANI equivalence, and clinical interpretation",
        ["Level", "Distance (d)", "Approx. ANI", "Median ANI (validated)", "Clinical interpretation", "Infection control use case"],
        [
            ["L1", "≤0.150", "~85%", "97.8%", "Plasmid superfamily", "Broad family identification"],
            ["L2", "≤0.100", "~90%", "97.8%", "Major lineage", "National/international surveillance"],
            ["L3", "≤0.050", "~95%", "97.8%", "Species-level cluster", "Regional surveillance; healthcare networks"],
            ["L4", "≤0.020", "~98%", "97.9%", "Sublineage", "Inter-hospital transmission tracking"],
            ["L5", "≤0.010", "~99%", "98.2%", "Clone group", "Intra-hospital; ward-level spread"],
            ["L6", "≤0.001", "~99.9%", "99.9%", "Lineage/outbreak", "Outbreak confirmation; direct transmission"],
        ],
        notes="Validated with FastANI v1.34 on 4,970 within-group plasmid pairs across 20 Gram-negative Inc groups. Spearman ρ=-0.348 (P<10⁻¹⁴¹).",
        col_widths=[10, 14, 14, 20, 22, 35])

    # Table 3: AMR lineages
    ws3 = wb.create_sheet("Table 3 - AMR Lineages")
    style_sheet(ws3, "Table 3: High-risk pLIN lineages — clinical surveillance targets",
        ["pLIN (L6)", "n", "Inc type(s)", "Mean AMR", "Key resistance determinants", "Clinical risk", "Suggested action"],
        [
            ["671", 90, "IncN", 13.2, "blaKPC-2 (100%), blaTEM-1 (100%), aph(3'')-Ib (98%)", "Carbapenem + aminoglycoside resistance", "Immediate isolation; carbapenem restriction"],
            ["860", 142, "5 Inc groups", 14.4, "sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), mcr-1.1 (36.6%)", "Pan-aminoglycoside + carbapenem + colistin", "Enhanced isolation; limited therapeutic options"],
            ["1327", 869, "6 Inc groups", 5.5, "blaTEM-1 (33.4%), sul1 (30.4%)", "Largest lineage; moderate-burden hub", "Surveillance monitoring; trend analysis"],
            ["1434", 163, "IncFII (92%), IncF (8%)", 7.4, "mph(A) (56.4%), mrx(A) (56.4%), sul1 (56.4%), blaCTX-M-27 (29.4%)", "ESBL-carrying IncFII; community UTI", "Community surveillance"],
        ],
        notes="Cross-validation: pLIN 671 matched CP104944 (Yao 2023, 61 hospitals, Germany); pLIN 860 matched CP022533 (Roberts 2020, Australia).",
        col_widths=[12, 8, 16, 12, 42, 32, 35])

    add_appendix_sheets(wb)
    path = os.path.join(BASE, "cmi", "tables", "CMI_Tables.xlsx")
    wb.save(path)
    print(f"  Saved: {path}")


def add_appendix_sheets(wb):
    """Add all 16 appendix/supplementary tables as sheets to the workbook."""

    # --- Appendix p1: Dataset breakdown ---
    ws = wb.create_sheet("App p1 - Dataset")
    style_sheet(ws, "Appendix p 1: Complete dataset breakdown for all 28 Inc/Rep groups",
        ["Inc/Rep group", "n", "% of dataset", "Unique pLIN codes (L6)"],
        [
            ["IncFII", "4,629", "57.3%", 1421],
            ["IncN", "1,097", "13.6%", 431],
            ["IncX1", 705, "8.7%", 420],
            ["repAci_large", 203, "2.5%", 112],
            ["repPae_large", 199, "2.5%", 108],
            ["repPae_small", 150, "1.9%", 82],
            ["repSA_large", 121, "1.5%", 78],
            ["repEF_res", "100*", "1.2%", 69],
            ["repAci1", 120, "1.5%", 68],
            ["IncFIB", 97, "1.2%", 42],
            ["repEF_conj", 92, "1.1%", 52],
            ["ColRNAI", 91, "1.1%", 38],
            ["IncF", 75, "0.9%", 31],
            ["repSA_small", 73, "0.9%", 41],
            ["IncX3", 56, "0.7%", 24],
            ["IncHI2", 36, "0.4%", 18],
            ["IncI1", 27, "0.3%", 14],
            ["IncI2", 25, "0.3%", 13],
            ["IncX4", 24, "0.3%", 12],
            ["IncR", 21, "0.3%", 11],
            ["ColE", 19, "0.2%", 9],
            ["IncC", 16, "0.2%", 8],
            ["IncHI1", 16, "0.2%", 8],
            ["IncFIC", 14, "0.2%", 7],
            ["IncAC2", 14, "0.2%", 7],
            ["IncA", 14, "0.2%", 7],
            ["IncI", 11, "0.1%", 6],
            ["IncFIBK", 11, "0.1%", 6],
            ["Total", "8,077*", "100%", 3073],
        ],
        notes="* 8,077 total training samples include 21 multi-replicon E. faecium plasmids counted in both repEF_conj and repEF_res (8,056 unique plasmids). repEF_res*: 121 training FASTAs, 100 unique after removing 21 multi-replicon overlap with repEF_conj.",
        col_widths=[16, 10, 14, 24])

    # --- Appendix p2: AMR prevalence by Inc group ---
    ws = wb.create_sheet("App p2 - AMR by Inc")
    style_sheet(ws, "Appendix p 2: AMR gene prevalence by Inc group (20 Gram-negative groups)",
        ["Inc group", "n", "AMR+ (n)", "AMR+ (%)", "Mean AMR genes", "VIR+ (n)", "VIR+ (%)", "Stress+ (n)", "Stress+ (%)"],
        [
            ["IncFII", "4,629", "2,882", "62.3%", 3.91, "1,076", "23.2%", "1,882", "40.7%"],
            ["IncN", "1,097", 936, "85.3%", 6.23, 32, "2.9%", 698, "63.6%"],
            ["IncX1", 705, 470, "66.7%", 4.20, 83, "11.8%", 298, "42.3%"],
            ["IncFIB", 97, 58, "59.8%", 3.12, 14, "14.4%", 42, "43.3%"],
            ["ColRNAI", 91, 42, "46.2%", 2.18, 6, "6.6%", 28, "30.8%"],
            ["IncF", 75, 51, "68.0%", 4.35, 11, "14.7%", 36, "48.0%"],
            ["IncX3", 56, 39, "69.6%", 5.87, 4, "7.1%", 21, "37.5%"],
            ["IncHI2", 36, 28, "77.8%", 7.42, 2, "5.6%", 24, "66.7%"],
            ["IncI1", 27, 19, "70.4%", 4.61, 3, "11.1%", 14, "51.9%"],
            ["IncI2", 25, 18, "72.0%", 5.33, 2, "8.0%", 12, "48.0%"],
            ["IncX4", 24, 16, "66.7%", 3.44, 1, "4.2%", 8, "33.3%"],
            ["IncR", 21, 15, "71.4%", 4.80, 1, "4.8%", 10, "47.6%"],
            ["ColE", 19, 8, "42.1%", 1.75, 2, "10.5%", 5, "26.3%"],
            ["IncC", 16, 12, "75.0%", 6.33, 1, "6.3%", 9, "56.3%"],
            ["IncHI1", 16, 11, "68.8%", 5.09, 2, "12.5%", 8, "50.0%"],
            ["IncFIC", 14, 9, "64.3%", 3.56, 3, "21.4%", 6, "42.9%"],
            ["IncAC2", 14, 10, "71.4%", 4.70, 1, "7.1%", 7, "50.0%"],
            ["IncA", 14, 10, "71.4%", 5.10, 1, "7.1%", 6, "42.9%"],
            ["IncI", 11, 8, "72.7%", 4.50, 0, "0.0%", 5, "45.5%"],
            ["IncFIBK", 11, 7, "63.6%", 2.86, 1, "9.1%", 4, "36.4%"],
            ["Total (GN)", "6,998", "4,657", "66.5%", 4.43, "1,265", "18.1%", "3,005", "42.9%"],
        ],
        notes="AMRFinderPlus analysis was performed on 20 Gram-negative Inc groups (6,998 plasmids) only. IncN: highest per-plasmid AMR risk (85.3%, 6.23 mean). IncHI2: 2nd highest burden (7.42 mean). IncFII: highest virulence (23.2%).",
        col_widths=[12, 10, 10, 10, 14, 10, 10, 12, 12])

    # --- Appendix p3: Top 20 AMR genes ---
    ws = wb.create_sheet("App p3 - Top 20 AMR")
    style_sheet(ws, "Appendix p 3: Top 20 AMR genes — clinical significance",
        ["Rank", "Gene", "Detections", "% of AMR+", "Drug class", "Clinical significance"],
        [
            [1, "blaTEM-1", "1,864", "40.0%", "Beta-lactam", "Narrow-spectrum beta-lactamase; ampicillin/penicillin resistance"],
            [2, "sul1", "1,403", "30.1%", "Sulfonamide", "Class 1 integron-associated; MDR cassette marker"],
            [3, "tet(A)", "1,365", "29.3%", "Tetracycline", "Efflux pump; first-line therapy failure"],
            [4, "sul2", "1,147", "24.6%", "Sulfonamide", "Broad-host-range; cotrimoxazole failure"],
            [5, "aph(6)-Id", "1,136", "24.4%", "Aminoglycoside", "Streptomycin resistance; linked with aph(3'')-Ib"],
            [6, "aph(3'')-Ib", "1,127", "24.2%", "Aminoglycoside", "Streptomycin resistance; Tn5393-like transposons"],
            [7, "mph(A)", "1,058", "22.7%", "Macrolide", "Macrolide phosphotransferase; azithromycin failure"],
            [8, "mrx(A)", "1,057", "22.7%", "Macrolide", "Co-located with mph(A)"],
            [9, "blaKPC-2", 824, "17.7%", "Carbapenem", "WHO critical; 40-50% BSI mortality"],
            [10, "qnrS1", 732, "15.7%", "Quinolone", "PMQR; empirical UTI/BSI therapy failure"],
            [11, "aac(6')-Ib-cr5", 737, "15.8%", "Aminoglycoside/Quinolone", "Bifunctional; dual aminoglycoside + ciprofloxacin resistance"],
            [12, "dfrA14", 689, "14.8%", "Trimethoprim", "DHF reductase; cotrimoxazole failure"],
            [13, "aadA2", 604, "13.0%", "Aminoglycoside", "Spectinomycin/streptomycin; integron-borne"],
            [14, "catB3", 569, "12.2%", "Phenicol", "Chloramphenicol acetyltransferase"],
            [15, "blaOXA-1", 549, "11.8%", "Beta-lactam", "Oxacillinase; co-carried with ESBLs"],
            [16, "floR", 503, "10.8%", "Phenicol", "Florfenicol/chloramphenicol efflux"],
            [17, "aac(3)-IId", 485, "10.4%", "Aminoglycoside", "Gentamicin resistance"],
            [18, "aph(3')-Ia", 482, "10.4%", "Aminoglycoside", "Kanamycin/neomycin resistance"],
            [19, "blaCTX-M-15", 505, "10.8%", "Beta-lactam (ESBL)", "WHO critical; 3rd-gen cephalosporin hydrolysis"],
            [20, "ble", 436, "9.4%", "Bleomycin", "Bleomycin-binding protein; co-selected on Tn"],
        ],
        col_widths=[8, 18, 12, 12, 22, 50])

    # --- Appendix p4: Top 10 high-risk lineages ---
    ws = wb.create_sheet("App p4 - Top10 Lineages")
    style_sheet(ws, "Appendix p 4: Full top-10 high-risk pLIN lineages — detailed clinical profiles",
        ["Rank", "pLIN code", "n (AMR+)", "Mean AMR", "Total det.", "Inc type(s)", "Top resistance genes", "Clinical risk assessment"],
        [
            [1, "1.1.2.15.48.1327", 869, 5.5, "4,780", "IncF, IncFIB, IncFIBK, IncFII, IncHI1, IncN", "blaTEM-1 (33.4%), sul1 (30.4%), tet(A) (28%), mph(A) (25%)", "Largest lineage; moderate AMR; 6 Inc groups; community-associated"],
            [2, "1.1.2.15.48.860", 142, 14.4, "2,045", "IncFII, IncHI1, IncHI2, IncN, IncX1", "sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), mcr-1.1 (36.6%)", "CRITICAL: Cross-Inc MDR hub; colistin + multi-drug resistance"],
            [3, "1.1.2.15.48.1434", 163, 7.4, "1,206", "IncFII (92.0%), IncF (8.0%)", "mph(A) (56.4%), mrx(A) (56.4%), sul1 (56.4%), blaCTX-M-27 (29.4%)", "ESBL lineage; IncFII-dominant; CTX-M-27 marker"],
            [4, "1.1.2.15.48.1465", 226, 8.2, "1,853", "IncN, IncFII", "rmtB1 (70%), blaKPC-2 (58%), blaCTX-M-65 (54%)", "CRITICAL: Pan-aminoglycoside + carbapenem; KPC-2 hub"],
            [5, "1.1.2.15.48.1335", 116, 6.1, 708, "IncFII", "blaTEM-1 (45%), sul1 (38%), tet(A) (35%)", "IncFII lineage; moderate AMR burden"],
            [6, "1.1.2.15.48.1361", 95, 5.8, 551, "IncFII", "blaTEM-1 (42%), mph(A) (36%), sul1 (33%)", "IncFII lineage; macrolide + sulfonamide resistance"],
            [7, "1.1.2.15.48.671", 90, 13.2, "1,188", "IncN", "blaKPC-2 (100%), blaTEM-1 (100%), aph(3'')-Ib (98%), aac(3)-IId (97%)", "CRITICAL: 100% carbapenemase; immediate isolation trigger"],
            [8, "1.1.2.15.48.1213", 88, 5.4, 475, "IncFII", "blaTEM-1 (40%), sul1 (34%), tet(A) (31%)", "IncFII lineage; moderate MDR burden"],
            [9, "1.1.2.15.48.1247", 79, 5.1, 403, "IncFII", "blaTEM-1 (38%), sul1 (32%), mph(A) (28%)", "IncFII lineage; community-associated"],
            [10, "1.1.2.15.48.1360", 73, 5.6, 409, "IncFII", "blaTEM-1 (44%), sul1 (37%), tet(A) (33%)", "IncFII lineage; moderate AMR; surveillance target"],
        ],
        col_widths=[8, 22, 12, 12, 12, 32, 50, 50])

    # --- Appendix p5: Critical resistance determinants ---
    ws = wb.create_sheet("App p5 - Carbapenems+mcr")
    style_sheet(ws, "Appendix p 5: Clinically critical resistance determinants — detailed breakdown",
        ["Category", "Variant", "Detections", "% of AMR+", "Primary Inc types", "Clinical impact"],
        [
            ["Carbapenemase", "blaKPC-2", 824, "17.7%", "IncN, IncFII", "Class A serine; hydrolyses all beta-lactams; inhibited by avibactam"],
            ["Carbapenemase", "blaNDM-1", 228, "4.9%", "IncFII, IncX3", "Metallo-BL; NOT inhibited by avibactam"],
            ["Carbapenemase", "blaKPC-3", 193, "4.1%", "IncN, IncFII", "Enhanced carbapenem hydrolysis"],
            ["Carbapenemase", "blaNDM-5", 89, "1.9%", "IncFII, IncX3", "Increased carbapenem hydrolysis; global spread on IncX3"],
            ["Carbapenemase", "blaIMP-4", 66, "1.4%", "IncHI2, IncN", "Metallo-BL; endemic Asia-Pacific"],
            ["Carbapenemase", "blaOXA-48-like", 42, "0.9%", "IncFII, IncR", "Low-level hydrolysis; often combined with ESBLs"],
            ["Carbapenemase", "blaVIM variants", 38, "0.8%", "IncN, IncI", "Metallo-BL; endemic Mediterranean"],
            ["Colistin", "mcr-1.1", 83, "—", "IncI2, IncX4, IncHI2", "Globally disseminated; first detected China 2015"],
            ["Colistin", "mcr-8.1", 27, "—", "IncFII", "Predominantly Asian; K. pneumoniae associated"],
            ["Colistin", "mcr-8.2", 19, "—", "IncFII", "Variant of mcr-8; similar geographic range"],
            ["Colistin", "mcr-10.1", 13, "—", "IncFII", "Emerging variant; limited epidemiological data"],
            ["Colistin", "mcr-3.5", 9, "—", "IncFII, IncHI2", "Asian and European isolates"],
            ["Colistin", "Other mcr", 53, "—", "Various", "Multiple emerging variants under surveillance"],
        ],
        notes="KPC: susceptible to ceftazidime-avibactam. NDM/VIM/IMP: resistant to all BL/BLI except cefiderocol + aztreonam-avibactam. Co-carriage KPC+NDM detected on 12 plasmids. mcr + carbapenemase co-carriage in 18 instances.",
        col_widths=[16, 18, 12, 12, 22, 55])

    # --- Appendix p6: Drug class co-carriage ---
    ws = wb.create_sheet("App p6 - Co-carriage")
    style_sheet(ws, "Appendix p 6: AMR drug class co-carriage analysis",
        ["Drug class", "Plasmids positive", "% of AMR+ (n=4,657)"],
        [
            ["Beta-lactam", "3,483", "74.8%"],
            ["Aminoglycoside", "2,723", "58.5%"],
            ["Sulfonamide", "2,073", "44.5%"],
            ["Trimethoprim", "1,810", "38.9%"],
            ["Tetracycline", "1,549", "33.3%"],
            ["Phenicol", "1,542", "33.1%"],
            ["Quinolone", "1,240", "26.6%"],
            ["Macrolide", "1,088", "23.4%"],
            ["", "", ""],
            ["Multi-drug resistance", "", ""],
            ["2+ drug classes", "3,891", "83.5% of AMR+"],
            ["3+ drug classes", "3,012", "64.7% of AMR+"],
            ["5+ drug classes", "1,847", "39.7% of AMR+"],
            ["7+ drug classes", 589, "12.6% of AMR+"],
        ],
        notes="Nearly two-thirds of AMR+ plasmids carry resistance to 3+ drug classes. 39.7% with 5+ classes would confer XDR in a single conjugation event.",
        col_widths=[24, 20, 24])

    # --- Appendix p7: Virulence distribution ---
    ws = wb.create_sheet("App p7 - Virulence")
    style_sheet(ws, "Appendix p 7: Virulence gene distribution by Inc type",
        ["Inc type", "VIR+ (%)", "Top virulence genes", "Clinical implication"],
        [
            ["IncFII", "23.2%", "traT (25.6%), spvB/spvD (25.3%), mltE (24.0%), iucA (23.3%)", "Salmonella virulence (spv), iron acquisition; enhanced invasive risk"],
            ["IncX1", "11.8%", "hlyA (48.2%), fedA/fedF (47.0%)", "Alpha-haemolysin, fimbrial adhesins (ETEC); diarrhoeal severity"],
            ["IncFIB", "14.4%", "traT (31.2%), iutA (22.1%)", "Serum resistance, aerobactin iron uptake; bacteraemia persistence"],
            ["IncN", "2.9%", "iutA (18.8%), iucC (15.6%), traT (12.5%)", "Low virulence but high AMR; primarily resistance vehicles"],
        ],
        notes="AMR-virulence co-carriage on IncFII (23.2% VIR + 62.3% AMR) creates co-selection risk under antibiotic pressure.",
        col_widths=[12, 12, 50, 50])

    # --- Appendix p8: ML validation detailed ---
    ws = wb.create_sheet("App p8 - ML Detailed")
    # Part A: nested CV design
    style_sheet(ws, "Appendix p 8: Machine learning validation — detailed methodology",
        ["Parameter", "Value"],
        [
            ["Outer loop", "3-fold stratified CV (StratifiedKFold, shuffle=True, seed=42)"],
            ["Inner loop", "2-fold stratified CV"],
            ["Scoring metric", "Weighted F1 (f1_weighted)"],
            ["Balanced subset", "1,500 plasmids (500 IncFII, 500 IncN, 500 IncX1)"],
            ["Feature scaling", "StandardScaler (z-score), fitted per fold"],
            ["Optimisation", "Optuna TPE sampler, Hyperband pruning, 10 trials/fold"],
            ["", ""],
            ["Feature vector (33-D):", ""],
            ["Basic composition", "7 features: length, log10(length), GC, AT, AT skew, GC skew"],
            ["Dinucleotide freq.", "16 canonical dinucleotides, normalised"],
            ["Trinucleotide freq.", "10 selected: ATG, TAA, TAG, TGA, GCG, CGC, AAA, TTT, CCC, GGG"],
        ],
        col_widths=[22, 65])
    # Part B: Top 10 features (start below existing data)
    feat_start = 4 + 11 + 3  # after data rows + gap
    ws.merge_cells(start_row=feat_start, start_column=1, end_row=feat_start, end_column=2)
    c = ws.cell(row=feat_start, column=1, value="Top 10 consensus feature importances")
    c.font = Font(name="Arial", bold=True, size=11)
    feat_headers = ["Rank", "Feature", "Mean importance", "Biological relevance"]
    for j, h in enumerate(feat_headers, 1):
        c = ws.cell(row=feat_start + 1, column=j, value=h)
        c.font = HEADER_FONT; c.fill = HEADER_FILL; c.alignment = CENTER; c.border = THIN_BORDER
    feat_data = [
        [1, "tri_TAG", 0.114, "Amber stop codon frequency; codon usage bias"],
        [2, "tri_GCG", 0.099, "CpG island proxy; DNA methylation signatures"],
        [3, "at_content", 0.076, "Overall AT composition; backbone vs accessory"],
        [4, "gc_content", 0.053, "GC%; correlates with host range"],
        [5, "di_AA", 0.050, "Poly-A tracts; regulatory element density"],
        [6, "tri_TGA", 0.048, "Opal stop codon; alternative stop codon usage"],
        [7, "tri_CGC", 0.040, "CpG-related; restriction-modification signatures"],
        [8, "di_TT", 0.039, "Poly-T tracts; rho-independent terminator signals"],
        [9, "di_GG", 0.032, "G-richness; G-quadruplex potential"],
        [10, "tri_GGG", 0.030, "Extreme GC-rich tracts"],
    ]
    for i, row in enumerate(feat_data):
        for j, val in enumerate(row, 1):
            c = ws.cell(row=feat_start + 2 + i, column=j, value=val)
            c.font = BODY_FONT; c.alignment = WRAP; c.border = THIN_BORDER
            if i % 2 == 1: c.fill = ALT_FILL
    for j, w in enumerate([8, 16, 18, 50], 1):
        ws.column_dimensions[get_column_letter(j)].width = max(
            ws.column_dimensions[get_column_letter(j)].width or 0, w)

    # --- Appendix p9: Threshold calibration ---
    ws = wb.create_sheet("App p9 - Thresholds")
    style_sheet(ws, "Appendix p 9: Hierarchical threshold calibration (IncX n=178 + FastANI n=4,970)",
        ["Quantile / Level", "NN distance (d)", "Calibration role / ANI equiv."],
        [
            ["Initial calibration (IncX-like, n=178)", "", ""],
            ["25th percentile", "0.001", "L6 threshold (lineage/outbreak)"],
            ["50th (median)", "0.011", "L5 threshold (clone group)"],
            ["75th percentile", "0.025", "L4 threshold (subcluster)"],
            ["95th percentile", "0.048", "L3 threshold (~95% ANI)"],
            ["99th percentile", "0.072", "L2 threshold range"],
            ["", "", ""],
            ["FastANI validation (20 Gram-negative Inc groups)", "", ""],
            ["L6 (d<=0.001)", "n=726 pairs", "Median ANI 99.9% (5th %ile: 97.4%)"],
            ["L5 (d<=0.010)", "n=2,696 pairs", "Median ANI 98.2% (5th %ile: 90.4%)"],
            ["L4 (d<=0.020)", "n=3,566 pairs", "Median ANI 97.9% (5th %ile: 89.5%)"],
            ["L3 (d<=0.050)", "n=4,554 pairs", "Median ANI 97.8% (5th %ile: 87.9%)"],
            ["L2 (d<=0.100)", "n=4,944 pairs", "Median ANI 97.8% (5th %ile: 85.3%)"],
            ["L1 (d<=0.150)", "n=4,966 pairs", "Median ANI 97.8% (5th %ile: 85.2%)"],
        ],
        notes="Overall: Spearman rho=-0.348 (P<10^-141). 15/20 Gram-negative groups significant (P<0.001). Strongest: IncFIC (rho=-0.883), IncI2 (-0.810), ColE (-0.722). FastANI validation was performed on 20 Gram-negative Inc groups only.",
        col_widths=[32, 18, 45])

    # --- Appendix p9 (cont): Per-Inc group correlations ---
    ws = wb.create_sheet("App p9b - Per-Inc ANI")
    style_sheet(ws, "Appendix p 9 (cont): Per-Inc-group Spearman correlation (cosine vs FastANI)",
        ["Inc group", "n pairs", "Spearman rho", "P-value", "Significance", "Mean ANI"],
        [
            ["IncFIC", 182, -0.883, "<0.001", "***", "98.5%"],
            ["IncI2", 308, -0.810, "<0.001", "***", "98.3%"],
            ["ColE", 254, -0.722, "<0.001", "***", "93.2%"],
            ["IncI", 90, -0.612, "<0.001", "***", "96.1%"],
            ["IncN", 230, -0.612, "<0.001", "***", "95.2%"],
            ["IncHI1", 75, -0.575, "<0.001", "***", "94.0%"],
            ["IncAC2", 182, -0.536, "<0.001", "***", "99.1%"],
            ["IncF", 343, -0.514, "<0.001", "***", "93.3%"],
            ["ColRNAI", 202, -0.502, "<0.001", "***", "96.4%"],
            ["IncFIB", 238, -0.475, "<0.001", "***", "94.5%"],
            ["IncI1", 306, -0.297, "<0.001", "***", "97.8%"],
            ["IncA", 182, -0.288, "<0.001", "***", "94.5%"],
            ["IncHI2", 380, -0.251, "<0.001", "***", "98.3%"],
            ["IncR", 378, -0.169, "<0.001", "***", "95.8%"],
            ["IncC", 240, -0.148, "0.022", "*", "98.6%"],
            ["IncX1", 274, -0.076, "0.207", "n.s.", "93.5%"],
            ["IncFII", 240, -0.070, "0.281", "n.s.", "89.4%"],
            ["IncFIBK", 106, -0.061, "0.531", "n.s.", "98.2%"],
            ["IncX4", 380, -0.016, "0.754", "n.s.", "98.7%"],
            ["IncX3", 380, 0.038, "0.467", "n.s.", "99.6%"],
        ],
        notes="5 non-significant groups (IncX1, IncFII, IncFIBK, IncX4, IncX3) have narrow ANI ranges (>89%), producing floor effects.",
        col_widths=[12, 10, 14, 12, 14, 12])

    # --- Appendix p10a: IncX1 backbone/mosaicism ---
    ws = wb.create_sheet("App p10a - IncX1 Backbone")
    style_sheet(ws, "Appendix p 10a: IncX1 backbone architecture and mosaicism (n=701)",
        ["ORF ID", "Position", "Length (nt)", "Prevalence (%)", "Median match", "Putative function"],
        [
            ["IncX1_core_ORF067", "42806-43253 (+)", 447, "84.8%", 0.852, "Candidate essential replication gene"],
            ["IncX1_core_ORF063", "38165-38552 (+)", 387, "71.3%", 0.943, "Operon cluster; maintenance function"],
            ["IncX1_core_ORF062", "37861-38197 (+)", 336, "71.3%", 0.934, "Operon cluster; linked to ORF063"],
            ["IncX1_core_ORF061", "37787-38093 (+)", 306, "71.3%", 0.900, "Operon cluster; functionally linked"],
            ["IncX1_core_ORF066", "41930-42767 (+)", 837, "61.2%", 0.832, "Large conserved ORF; potential transfer"],
            ["", "", "", "", "", ""],
            ["Backbone stats:", "", "", "", "", ""],
            ["Mean backbone length", "37,787 bp", "60%", "", "", "of total plasmid"],
            ["Mean accessory length", "25,192 bp", "40%", "", "", "of total plasmid"],
            ["Backbone range", "48-82%", "", "", "", "across individual plasmids"],
            ["Mosaic candidates", "302/701", "43.1%", "", "", "GC-content heterogeneity"],
            ["Mosaic mean GC-SD", "0.047", "", "", "", "range 0.017-0.115"],
            ["Non-mosaic GC-SD", "0.021", "", "", "", "range 0.008-0.032"],
        ],
        notes="43.1% mosaic IncX1 plasmids confirms active acquisition of foreign DNA modules including AMR cassettes.",
        col_widths=[22, 22, 14, 16, 16, 40])

    # --- Appendix p10b: Outbreak cross-validation (expanded) ---
    ws = wb.create_sheet("App p10b - Outbreak Valid")
    style_sheet(ws, "Appendix p 10b: Cross-validation against 27 published outbreak studies (74 plasmids, 13 countries)",
        ["Gene class", "n", "High conf.", "Unique pLIN", "Key studies", "Key finding"],
        [
            ["blaKPC-2/3", "13", "12 (92%)", "10", "Yao 2023, Conlan 2014, Li 2018", "pLIN 671 hotspot confirmed (Germany + China)"],
            ["blaNDM-1/5", "36", "30 (83%)", "19", "Ho 2019, Weber 2019, Rojas 2017", "12/13 HK outbreak share pLIN 475; 4 German share 492"],
            ["blaOXA-48", "6", "5 (83%)", "2", "Potron 2013, Jousset 2019", "5/6 share pLIN 1688 (Turkey=NL=France)"],
            ["blaVIM-1", "3", "1 (33%)", "2", "Arcari 2020", "2/3 share outbreak backbone pLIN 490"],
            ["blaIMP-4/6", "2", "2 (100%)", "2", "Roberts 2020, Tada 2015", "pLIN 860 MDR hub match (Australia)"],
            ["mcr-1", "6", "6 (100%)", "2", "Liu 2016, Zheng 2017, Hasman 2015", "IncI2 (pLIN 87) vs IncX4 (pLIN 340) separated"],
            ["blaCTX-M-15", "8", "7 (88%)", "6", "Sheppard 2016, Woodford 2009", "USA outbreak: 3 share pLIN 1482"],
            ["", "", "", "", "", ""],
            ["TOTAL", "74", "63 (85.1%)", "42", "26 studies, 13 countries", "9 intra-study clusters; 3 global lineages matched"],
        ],
        notes="Expanded from original pilot (n=17, 4 studies) to 74 plasmids from 26 studies across 13 countries. Mean NN distance=0.0062. Three globally disseminated lineages: pLIN 671/KPC-2, pLIN 860/MDR hub, pLIN 1688/OXA-48.",
        col_widths=[16, 8, 14, 14, 34, 50])

    # --- Appendix p11: Outbreak detection algorithm ---
    ws = wb.create_sheet("App p11 - Outbreak Algo")
    style_sheet(ws, "Appendix p 11: Outbreak detection algorithm specification",
        ["Component", "Parameter", "Value / Rule"],
        [
            ["Basic module", "Grouping criterion", "Same L6 pLIN code + same AMR fingerprint"],
            ["Basic module", "Minimum cluster size", ">=2 plasmids"],
            ["Basic module", "Risk: HIGH", ">=3 shared AMR genes"],
            ["Basic module", "Risk: MODERATE", "<3 shared AMR genes"],
            ["", "", ""],
            ["Temporal module", "Time window", "Configurable (default 30 days)"],
            ["Temporal module", "Risk: CRITICAL", ">=3 AMR genes AND <=7 days apart"],
            ["Temporal module", "Risk: HIGH", ">=3 AMR genes OR <=7 days apart"],
            ["Temporal module", "Risk: MODERATE", "<3 AMR genes AND >7 days apart"],
            ["", "", ""],
            ["SNP sub-typing", "Method", "minimap2 -cx asm5 --cs"],
            ["SNP sub-typing", "0 SNPs", "Potentially clonal; strongest transmission evidence"],
            ["SNP sub-typing", "1-5 SNPs", "Highly related; recent divergence"],
            ["SNP sub-typing", "6-20 SNPs", "Related but divergent; indirect transmission"],
            ["SNP sub-typing", ">20 SNPs", "Distinct within strain cluster; independent events"],
            ["", "", ""],
            ["Design rationale", "Dual criteria", "L6 identity + AMR fingerprint minimises false positives"],
        ],
        notes="Plasmids sharing L6 code but different AMR fingerprints = independent acquisitions on same backbone, not clonal spread.",
        col_widths=[20, 22, 55])

    # --- Appendix p12: Reference database expansion ---
    ws = wb.create_sheet("App p12 - Ref Expansion")
    style_sheet(ws, "Appendix p 12: Reference database expansion — 79,305 plasmids",
        ["Inc group", "Training", "Reference", "Total", "L6 pLIN codes", "Dist. matrix size"],
        [
            ["IncFII", "4,629", "29,407", "34,036", "14,403", "2.3 GB (chunked)"],
            ["IncX1", 705, "25,449", "26,154", "14,060", "1.4 GB (chunked)"],
            ["IncN", "1,097", "4,206", "5,303", "1,912", "56 MB"],
            ["ColRNAI", 91, "1,836", "1,927", 618, "7.4 MB"],
            ["IncX3", 56, "1,403", "1,459", 323, "4.3 MB"],
            ["IncFIB", 97, "1,022", "1,119", 279, "2.5 MB"],
            ["IncI1", 27, 920, 947, 82, "1.8 MB"],
            ["IncI2", 25, 851, 876, 229, "1.5 MB"],
            ["IncHI1", 16, 697, 713, 115, "1.0 MB"],
            ["IncAC2", 14, 693, 707, 262, "1.0 MB"],
            ["IncX4", 24, 669, 693, 195, "0.96 MB"],
            ["IncHI2", 36, 607, 643, 152, "0.82 MB"],
            ["IncF", 75, 412, 487, 120, "0.47 MB"],
            ["IncA", 14, 331, 345, 188, "0.24 MB"],
            ["IncC", 16, 324, 340, 107, "0.23 MB"],
            ["ColE", 19, 307, 326, 74, "0.21 MB"],
            ["IncI", 11, 305, 316, 60, "0.20 MB"],
            ["IncR", 21, 70, 91, 67, "0.02 MB"],
            ["IncFIC", 14, 28, 42, 22, "<0.01 MB"],
            ["IncFIBK", 11, 14, 25, 8, "<0.01 MB"],
        ],
        notes="Classification: 77,196 (97.3%) classified with >=40% confidence, 2,109 (2.7%) low-confidence. Runtime: 28 min total (Apple M-series, single-threaded).",
        col_widths=[12, 12, 12, 10, 16, 20])

    # --- Appendix p13: Clinical use case scenarios ---
    ws = wb.create_sheet("App p13 - Use Cases")
    style_sheet(ws, "Appendix p 13: Clinical use case scenarios",
        ["Scenario", "Setting", "Steps", "Outcome", "Without pLIN"],
        [
            ["1. ICU carbapenemase outbreak", "Hospital ICU; 5 CRE isolates in 2 weeks",
             "Upload FASTA -> pLIN assigns codes -> all 5 get pLIN 671 -> AMRFinderPlus: blaKPC-2 on all -> temporal: CRITICAL -> SNP: 0-2 SNPs -> clonal",
             "Enhanced contact precautions, environmental screening, stewardship review",
             "Lab reports 'IncN plasmid with KPC' x5; cannot distinguish clonal spread from independent acquisitions"],
            ["2. Regional colistin surveillance", "12 hospitals over 6 months; MCR+ isolates",
             "Assign pLIN -> 3 hospitals share L6 codes -> 2 others share L4 prefix -> rest unrelated",
             "Regional spread pattern identified; clonal expansion vs independent emergence distinguished",
             "All reported as 'MCR-positive'; spread pattern invisible"],
            ["3. Longitudinal surveillance", "Hospital A (Year 1) and Hospital B (Year 3, different country)",
             "pLIN 4.11.83.485.1175.3594 assigned Year 1 -> same code detected Year 3 -> international spread",
             "Cross-border lineage tracking with stable nomenclature",
             "MOB-suite/COPLA codes change with each DB update; longitudinal comparison impossible"],
        ],
        col_widths=[24, 28, 55, 40, 45])

    # --- Appendix p14: Software environment ---
    ws = wb.create_sheet("App p14 - Software")
    style_sheet(ws, "Appendix p 14: Software and computational environment",
        ["Software", "Version", "Purpose"],
        [
            ["Python", "3.10.19", "Primary programming language"],
            ["NumPy", "2.2.6", "Numerical computation"],
            ["Pandas", "2.3.3", "Data manipulation"],
            ["SciPy", "1.13.1", "Distance computation, hierarchical clustering"],
            ["BioPython", "1.86", "FASTA sequence parsing"],
            ["scikit-learn", "1.7.2", "KNN classifier, ML validation, cross-validation"],
            ["Streamlit", "1.54.0", "Interactive web application framework"],
            ["AMRFinderPlus", "4.2.5", "AMR/virulence/stress gene detection"],
            ["AMRFinderPlus DB", "2026-01-21.1", "Reference database for gene detection"],
            ["minimap2", "2.30", "SNP sub-typing within L6 clusters"],
            ["BLAST+", "2.16.0", "CRISPR spacer-plasmid matching"],
            ["MinCED", "0.4.2", "CRISPR spacer extraction"],
            ["Prodigal", "2.6.3", "Gene/ORF prediction"],
            ["FastANI", "1.34", "True ANI computation"],
            ["Matplotlib", "3.10.8", "Figure generation"],
            ["Seaborn", "0.13.2", "Statistical visualisation"],
            ["openpyxl", "3.1.5", "Excel workbook generation"],
            ["DendroPy", "5.0.8", "Phylogenetic tree handling"],
            ["PyTorch", "2.10.0", "Nucleotide Transformer embeddings (optional)"],
        ],
        notes="Hardware: Apple M-series laptop (macOS Darwin 24.5.0). No HPC, GPU, or cloud required.",
        col_widths=[20, 16, 50])

    # --- Appendix p15: Data availability ---
    ws = wb.create_sheet("App p15 - Data Avail")
    style_sheet(ws, "Appendix p 15: Data availability and supplementary tables",
        ["Table", "Description", "Key columns"],
        [
            ["S1", "Complete pLIN assignment (8,077 training plasmids)", "plasmid_id, inc_type, length_bp, pLIN, bin_A through bin_F"],
            ["S2", "Integrated pLIN + AMRFinderPlus (8,077 plasmids)", "plasmid_id, inc_type, pLIN, n_amr_genes, n_vir_genes, n_stress_genes, amr_genes, amr_classes"],
            ["S3", "pLIN lineage AMR summary (3,073 lineages)", "pLIN, n_plasmids, mean_amr, total_amr, inc_types, top_amr_genes"],
            ["S4", "Full reference database pLIN (79,305 plasmids)", "plasmid_id, inc_type, inc_confidence, pLIN, source (training/reference), is_novel"],
            ["S5", "Critical resistance gene detections", "plasmid_id, pLIN, gene_symbol, gene_category, drug_class"],
        ],
        notes="All data at https://github.com/xavierbasilbritto-hub/pLIN-plasmid-classification (GPL-3.0). All sequences from NCBI RefSeq.",
        col_widths=[10, 45, 60])

    return wb


if __name__ == "__main__":
    print("Generating Excel tables for all four journals...")
    print("(Including main tables + all 16 appendix/supplementary sheets)")
    print()
    print("1. Lancet Microbe (6 main tables + 16 appendix sheets)")
    create_lancet_microbe()
    print()
    print("2. Briefings in Bioinformatics (3 main tables + 16 appendix sheets)")
    create_briefings()
    print()
    print("3. Microbial Genomics (6 main tables + 16 appendix sheets)")
    create_microbial_genomics()
    print()
    print("4. CMI (3 main tables + 16 appendix sheets)")
    create_cmi()
    print()
    print("Done! All Excel files created with appendix/supplementary tables.")
