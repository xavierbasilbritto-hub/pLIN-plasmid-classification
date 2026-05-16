#!/usr/bin/env python3
"""
Create Excel workbook with all manuscript tables and supplementary data.
"""

import os
import pandas as pd
import numpy as np
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(BASE_DIR, "output", "manuscripts")
os.makedirs(OUT_DIR, exist_ok=True)


def style_header(ws, ncols, row=1):
    """Apply header styling."""
    header_fill = PatternFill(start_color="1565C0", end_color="1565C0", fill_type="solid")
    header_font = Font(name="Calibri", bold=True, color="FFFFFF", size=11)
    for col in range(1, ncols + 1):
        cell = ws.cell(row=row, column=col)
        cell.fill = header_fill
        cell.font = header_font
        cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)


def style_data(ws, nrows, ncols, start_row=2):
    """Apply data styling with alternating rows."""
    alt_fill = PatternFill(start_color="E3F2FD", end_color="E3F2FD", fill_type="solid")
    data_font = Font(name="Calibri", size=10)
    thin_border = Border(
        bottom=Side(style="thin", color="DDDDDD")
    )
    for row in range(start_row, start_row + nrows):
        for col in range(1, ncols + 1):
            cell = ws.cell(row=row, column=col)
            cell.font = data_font
            cell.border = thin_border
            if (row - start_row) % 2 == 1:
                cell.fill = alt_fill


def auto_width(ws):
    """Auto-adjust column widths."""
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            if cell.value:
                max_len = max(max_len, len(str(cell.value)))
        ws.column_dimensions[col_letter].width = min(max_len + 4, 50)


def create_table1_comparison(wb):
    """Table 1: Comparison of plasmid classification methods."""
    ws = wb.create_sheet("Table 1 - Method Comparison")

    headers = ["Feature", "PlasmidFinder", "pMLST", "MOB-suite", "COPLA", "mge-cluster", "pLIN"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["Classification type", "Flat (Inc group)", "Flat (ST)", "Flat (cluster)", "Semi-hierarchical (PTU)", "Flat (cluster)", "Hierarchical (6 levels)"],
        ["Resolution levels", "1", "1", "1", "2", "1", "6"],
        ["Code permanence", "Yes (stable DB)", "Yes (stable DB)", "No (reassigned)", "No (recomputed)", "No (re-embedded)", "Yes (guaranteed)"],
        ["Reference-free", "No", "No", "Partial", "No", "Yes", "Yes"],
        ["Contig identification", "No", "No", "No", "No", "No", "Yes (multi-signal scoring)"],
        ["Discriminatory power (D)", "0.641", "N/A", "N/A", "N/A", "N/A", "0.985"],
        ["AMR integration", "No", "No", "No", "No", "No", "Yes (AMRFinderPlus)"],
        ["Outbreak detection", "No", "No", "No", "No", "No", "Yes (2-tier)"],
        ["Risk stratification", "No", "No", "No", "No", "No", "Yes (3-tier)"],
        ["Coverage rate", "Varies", "61%", "61%", "41%", "100%", "97.3%"],
        ["Inc groups covered", "All known", "6 schemes", "All", "N/A", "N/A", "28 (20 Gram-neg + 4 Gram-pos + 2 Acinetobacter + 2 Pseudomonas)"],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_table2_thresholds(wb):
    """Table 2: pLIN hierarchical distance thresholds."""
    ws = wb.create_sheet("Table 2 - Thresholds")

    headers = ["Level", "Threshold (d)", "ANI Equivalent", "Biological Interpretation", "Clusters (n=8,077)", "Clusters (n=79,305)"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["L1", "<=0.150", "~85%", "Family / Superfamily", 1, 33],
        ["L2", "<=0.100", "~90%", "Subfamily / Major lineage", 1, 86],
        ["L3", "<=0.050", "~95%", "Cluster / Species-level", 7, 447],
        ["L4", "<=0.020", "~98%", "Subcluster / Sublineage", 38, 2772],
        ["L5", "<=0.010", "~99%", "Clone complex", 117, 5540],
        ["L6", "<=0.001", "~99.9%", "Strain / Outbreak level", 3073, 57886],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_table3_dataset(wb):
    """Table 3: Dataset composition by Inc group."""
    ws = wb.create_sheet("Table 3 - Dataset Composition")

    headers = ["Inc Group", "Count", "Percentage (%)", "Unique L6 Codes", "Singletons", "Simpson's D"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["IncFII", 4629, 66.1, 1220, 902, 0.962],
        ["IncN", 1097, 15.7, 580, 432, 0.976],
        ["IncX1", 705, 10.1, 402, 328, 0.995],
        ["IncFIB", 92, 1.3, 62, 49, 0.988],
        ["IncHI2", 80, 1.1, 30, 18, 0.945],
        ["IncX3", 67, 1.0, 27, 17, 0.961],
        ["IncFIBK", 60, 0.9, 34, 24, 0.978],
        ["ColRNAI", 47, 0.7, 32, 27, 0.988],
        ["IncI1", 37, 0.5, 27, 22, 0.990],
        ["IncC", 30, 0.4, 10, 6, 0.889],
        ["IncF", 29, 0.4, 16, 12, 0.966],
        ["IncR", 28, 0.4, 10, 5, 0.878],
        ["ColE", 24, 0.3, 18, 15, 0.993],
        ["IncI", 22, 0.3, 12, 8, 0.961],
        ["IncFIC", 17, 0.2, 7, 4, 0.882],
        ["IncA", 14, 0.2, 4, 2, 0.769],
        ["IncX4", 10, 0.1, 8, 7, 0.978],
        ["IncI2", 6, 0.1, 5, 4, 0.933],
        ["IncAC2", 3, 0.04, 2, 1, 0.667],
        ["IncHI1", 1, 0.01, 1, 1, 0.000],
        ["repSA_large", 121, 1.6, 78, 62, 0.985],
        ["repSA_small", 73, 1.0, 41, 33, 0.978],
        ["repEF_conj", 92, 1.1, 52, 41, 0.983],
        ["repEF_res", "100*", 1.2, 69, 55, 0.982],
        ["repAci1", 120, 1.5, 68, 54, 0.981],
        ["repAci_large", 203, 2.5, 112, 89, 0.984],
        ["repPae_large", 199, 2.5, 108, 85, 0.982],
        ["repPae_small", 150, 1.9, 82, 65, 0.982],
        ["TOTAL", "8,077*", 100.0, 3073, 2335, 0.985],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_table4_ml(wb):
    """Table 4: Machine learning cross-validation results."""
    ws = wb.create_sheet("Table 4 - ML Validation")

    headers = ["Model", "Weighted F1", "SD", "Accuracy (%)", "Notes"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["KNN (k=5, cosine)", "0.911", "N/A", "91.1", "Primary classifier; 5-fold CV; 28 groups"],
        ["XGBoost", "0.896", "0.011", "89.6", "Nested CV (5 outer, 5 inner)"],
        ["Gradient Boosting", "0.893", "0.009", "89.3", "Nested CV"],
        ["Random Forest", "0.874", "0.021", "87.4", "Nested CV"],
        ["Logistic Regression", "0.866", "0.012", "86.6", "Linear baseline"],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_table5_amr(wb):
    """Table 5: AMR gene detection summary."""
    ws = wb.create_sheet("Table 5 - AMR Summary")

    headers = ["Category", "Gene/Class", "Detections (n)", "Plasmids (%)", "Notes"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["Overall", "Total detections", 64891, "83.1% (5,816/6,998)", "AMR + Virulence + Stress; 20 GN groups only"],
        ["Overall", "AMR genes", 29583, "66.5% (4,657/6,998)", ""],
        ["Overall", "Virulence", 6286, "", ""],
        ["Overall", "Stress response", 29022, "", ""],
        ["Carbapenemase", "blaKPC-2", 824, "", "Most prevalent carbapenemase"],
        ["Carbapenemase", "blaNDM-1", 228, "", ""],
        ["Carbapenemase", "blaKPC-3", 193, "", ""],
        ["Carbapenemase", "blaNDM-5", 89, "", ""],
        ["Carbapenemase", "blaIMP-4", 66, "", ""],
        ["Carbapenemase", "Total", 1635, "", ""],
        ["ESBL", "blaCTX-M-15", 505, "", "Most prevalent ESBL"],
        ["ESBL", "blaCTX-M-65", 319, "", ""],
        ["ESBL", "blaSHV-12", 277, "", ""],
        ["ESBL", "Total", 1804, "", ""],
        ["Colistin (mcr)", "mcr-1.1", 83, "", ""],
        ["Colistin (mcr)", "mcr-8.1", 27, "", ""],
        ["Colistin (mcr)", "mcr-8.2", 18, "", ""],
        ["Colistin (mcr)", "Total", 204, "", ""],
        ["PMQR", "qnrS1", 732, "", ""],
        ["PMQR", "aac(6')-Ib-cr5", 737, "", ""],
        ["PMQR", "Total", 2315, "", ""],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_table6_lineages(wb):
    """Table 6: High-risk pLIN lineages."""
    ws = wb.create_sheet("Table 6 - High-Risk Lineages")

    headers = ["pLIN Code", "Inc Group(s)", "Members (n)", "Mean AMR Genes",
               "Key Resistance", "mcr (%)", "Clinical Significance"]
    for i, h in enumerate(headers, 1):
        ws.cell(row=1, column=i, value=h)

    data = [
        ["pLIN 1327", "6 Inc groups", 869, 5.5, "blaTEM-1 (32.8%), sul1 (30.3%)", "N/A", "Largest convergence zone"],
        ["pLIN 1434", "IncFII (92.0%), IncF (8.0%)", 163, 7.4, "blaCTX-M-27 (29.4%)", "0%", "ESBL lineage"],
        ["pLIN 860", "5 Inc groups", 142, 14.4, "Multi-class", "44.4%", "Cross-Inc MDR hub"],
        ["pLIN 671", "IncN", 90, 13.2, "blaKPC-2 (100%)", "0%", "KPC-2 hotspot lineage"],
    ]
    for r, row_data in enumerate(data, 2):
        for c, val in enumerate(row_data, 1):
            ws.cell(row=r, column=c, value=val)

    style_header(ws, len(headers))
    style_data(ws, len(data), len(headers))
    auto_width(ws)


def create_outbreak_validation_sheet(wb):
    """Supplementary: Outbreak cross-validation results."""
    ws = wb.create_sheet("S1 - Outbreak Validation")

    # Load the actual data
    tsv_path = os.path.join(BASE_DIR, "output", "outbreak_validation_combined_results.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, len(df), len(headers))
        auto_width(ws)
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def create_host_plasmid_validation_sheet(wb):
    """Supplementary: Combined chromosomal-plasmid validation."""
    ws = wb.create_sheet("S2 - Host-Plasmid Validation")

    tsv_path = os.path.join(BASE_DIR, "output", "retrospective_host_plasmid_validation.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, len(df), len(headers))
        auto_width(ws)
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def create_study_results_sheet(wb):
    """Supplementary: Per-study validation results."""
    ws = wb.create_sheet("S3 - Study Results")

    tsv_path = os.path.join(BASE_DIR, "output", "retrospective_validation_study_results.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, len(df), len(headers))
        auto_width(ws)
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def create_fastani_sheet(wb):
    """Supplementary: FastANI validation data."""
    ws = wb.create_sheet("S4 - FastANI Validation")

    tsv_path = os.path.join(BASE_DIR, "output", "cosine_to_ani_per_inc_summary.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, len(df), len(headers))
        auto_width(ws)
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def create_loocv_sheet(wb):
    """Supplementary: LOOCV benchmark results."""
    ws = wb.create_sheet("S5 - LOOCV Benchmark")

    tsv_path = os.path.join(BASE_DIR, "output", "outbreak_benchmark_results.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, len(df), len(headers))
        auto_width(ws)
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def create_plin_assignments_sheet(wb):
    """Supplementary: pLIN assignments for all training plasmids."""
    ws = wb.create_sheet("S6 - pLIN Assignments")

    tsv_path = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
    if os.path.exists(tsv_path):
        df = pd.read_csv(tsv_path, sep="\t")
        # Limit to first 5000 rows for Excel manageability
        if len(df) > 5000:
            df = df.head(5000)
            note = True
        else:
            note = False
        headers = list(df.columns)
        for i, h in enumerate(headers, 1):
            ws.cell(row=1, column=i, value=h)
        for r, (_, row) in enumerate(df.iterrows(), 2):
            for c, val in enumerate(row.values, 1):
                ws.cell(row=r, column=c, value=str(val) if pd.notna(val) else "")
        style_header(ws, len(headers))
        style_data(ws, min(len(df), 5000), len(headers))
        auto_width(ws)
        if note:
            ws.cell(row=5003, column=1,
                    value="Note: Only first 5,000 of 8,077 rows shown. Full data in pLIN_assignments.tsv")
    else:
        ws.cell(row=1, column=1, value="File not found: " + tsv_path)


def main():
    print("Creating manuscript tables Excel workbook...")

    wb = Workbook()
    # Remove default sheet
    wb.remove(wb.active)

    # Create all sheets
    create_table1_comparison(wb)
    create_table2_thresholds(wb)
    create_table3_dataset(wb)
    create_table4_ml(wb)
    create_table5_amr(wb)
    create_table6_lineages(wb)
    create_outbreak_validation_sheet(wb)
    create_host_plasmid_validation_sheet(wb)
    create_study_results_sheet(wb)
    create_fastani_sheet(wb)
    create_loocv_sheet(wb)
    create_plin_assignments_sheet(wb)

    out_path = os.path.join(OUT_DIR, "pLIN_manuscript_tables.xlsx")
    wb.save(out_path)
    print(f"Saved: {out_path}")
    print(f"Sheets: {wb.sheetnames}")


if __name__ == "__main__":
    main()
