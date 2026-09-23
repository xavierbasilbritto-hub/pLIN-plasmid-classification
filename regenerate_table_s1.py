#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
"""
Regenerate the two pLIN-code-dependent sheets in the npjAMR submission's
Supplementary Data Tables workbook ("pLIN Assignments" and "Lineage AMR
Summary") from the corrected, deterministic pLIN_assignments.tsv and the
freshly re-integrated AMR data.

The other three sheets (AMR by Replicon Group, Top AMR Genes, ANI
Validation) are keyed by Inc group / gene name, not by pLIN code, so they
are unaffected by the code-renumbering fix and are carried over unchanged.
"""

import os
import shutil
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter

BASE_DIR = "/Users/basilxavier/Desktop/PLASMID_TOOL"
SUBMISSION_DIR = os.path.join(
    BASE_DIR, "output", "manuscripts", "Revisions_Journal_Infection", "npjAMR_Submission_Full")
XLSX_PATH = os.path.join(SUBMISSION_DIR, "tables", "pLIN_Supplementary_Data_Tables.xlsx")

PLIN_ASSIGNMENTS_TSV = os.path.join(BASE_DIR, "output", "pLIN_assignments.tsv")
AMR_INTEGRATED_TSV = os.path.join(BASE_DIR, "output", "integrated", "pLIN_AMR_integrated.tsv")
LINEAGE_SUMMARY_TSV = os.path.join(BASE_DIR, "output", "integrated", "pLIN_lineage_AMR_summary.tsv")

HEADER_FONT = Font(name="Calibri", size=11, bold=True, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="1565C0", end_color="1565C0", fill_type="solid")
DATA_FONT = Font(name="Calibri", size=10)


def write_sheet(ws, df, col_widths=None):
    ws.delete_rows(1, ws.max_row)
    headers = list(df.columns)
    for j, h in enumerate(headers, 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font = HEADER_FONT
        c.fill = HEADER_FILL
        c.alignment = Alignment(horizontal="center")
    for i, row in enumerate(df.itertuples(index=False), 2):
        for j, val in enumerate(row, 1):
            c = ws.cell(row=i, column=j, value=val)
            c.font = DATA_FONT
    widths = col_widths or [max(12, len(h) + 2) for h in headers]
    for j, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "A2"


def build_plin_assignments_df():
    df = pd.read_csv(PLIN_ASSIGNMENTS_TSV, sep="\t")
    # Column order: keep original 10 columns first (for backward compatibility
    # with anything that indexes by position), then the two new columns
    # added for reviewer items 4 and 6.
    ordered = ["plasmid_id", "inc_type", "length_bp", "pLIN",
               "bin_A", "bin_B", "bin_C", "bin_D", "bin_E", "bin_F",
               "resolved_inc_type"]
    return df[ordered]


def build_lineage_summary_df():
    """Rebuild the Lineage AMR Summary sheet in the same shape as the
    original submitted sheet (pLIN code, n plasmids, Inc group(s), Mean AMR,
    Mean VIR, Top 3 AMR genes, AMR %, VIR %), sourced from the freshly
    regenerated AMR integration.
    """
    amr = pd.read_csv(AMR_INTEGRATED_TSV, sep="\t")

    rows = []
    for code, g in amr.groupby("pLIN"):
        n = len(g)
        inc_groups = ", ".join(sorted(g["inc_type_x"].unique()))
        mean_amr = g["n_amr_genes"].mean()
        mean_vir = g["n_vir_genes"].mean()
        amr_pct = (g["n_amr_genes"] > 0).mean() * 100
        vir_pct = (g["n_vir_genes"] > 0).mean() * 100

        gene_counts = {}
        for genes in g["amr_genes"].dropna():
            if genes == "none":
                continue
            for gene in genes.split(";"):
                gene = gene.strip()
                if gene:
                    gene_counts[gene] = gene_counts.get(gene, 0) + 1
        top3 = ", ".join(f"{g_}" for g_, _ in sorted(gene_counts.items(), key=lambda x: -x[1])[:3])

        rows.append({
            "pLIN code": code,
            "n plasmids": n,
            "Inc group(s)": inc_groups,
            "Mean AMR": round(mean_amr, 1),
            "Mean VIR": round(mean_vir, 1),
            "Top 3 AMR genes": top3,
            "AMR %": round(amr_pct, 1),
            "VIR %": round(vir_pct, 1),
        })

    out = pd.DataFrame(rows).sort_values("n plasmids", ascending=False).reset_index(drop=True)
    return out


def main():
    backup_path = XLSX_PATH + ".pre_pipeline_fix_backup.xlsx"
    if not os.path.exists(backup_path):
        shutil.copy(XLSX_PATH, backup_path)
        print(f"Backed up original workbook to: {backup_path}")

    wb = openpyxl.load_workbook(XLSX_PATH)

    plin_df = build_plin_assignments_df()
    ws1 = wb["pLIN Assignments (n=8,057)"]
    # Rename sheet to reflect the corrected, exact row count
    new_title = f"pLIN Assignments (n={len(plin_df):,})"
    ws1.title = new_title
    write_sheet(ws1, plin_df, col_widths=[29, 14, 12, 20, 8, 8, 8, 8, 8, 8, 18])
    print(f"Rewrote '{new_title}': {len(plin_df)} rows")

    lineage_df = build_lineage_summary_df()
    ws2 = wb["Lineage AMR Summary"]
    write_sheet(ws2, lineage_df, col_widths=[20, 12, 30, 10, 10, 45, 8, 8])
    print(f"Rewrote 'Lineage AMR Summary': {len(lineage_df)} rows")

    wb.save(XLSX_PATH)
    print(f"\nSaved: {XLSX_PATH}")


if __name__ == "__main__":
    main()
