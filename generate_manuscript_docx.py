#!/usr/bin/env python3
"""
Generate DOCX manuscripts and Excel supplementary files for all 4 journals.

Outputs:
  manuscripts/<journal>/
    <Journal>_Manuscript.docx          -- main manuscript + figure legends
    <Journal>_Supplementary_Tables.xlsx -- all supplementary tables + appendices
"""

import os
import re
import glob
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter

BASE = "/Users/basilxavier/Desktop/PLASMID_TOOL/manuscripts"

JOURNALS = {
    "lancet_microbe": {
        "name": "Lancet Microbe",
        "short": "Lancet_Microbe",
        "font": "Times New Roman",
        "font_size": 12,
        "line_spacing": 2.0,
    },
    "cmi": {
        "name": "Clinical Microbiology and Infection",
        "short": "CMI",
        "font": "Times New Roman",
        "font_size": 12,
        "line_spacing": 2.0,
    },
    "microbial_genomics": {
        "name": "Microbial Genomics",
        "short": "Microbial_Genomics",
        "font": "Times New Roman",
        "font_size": 12,
        "line_spacing": 2.0,
    },
    "briefings_bioinformatics": {
        "name": "Briefings in Bioinformatics",
        "short": "Briefings_Bioinformatics",
        "font": "Times New Roman",
        "font_size": 12,
        "line_spacing": 2.0,
    },
}

# ============================================================================
# DOCX HELPERS
# ============================================================================

def set_paragraph_format(para, font_name, font_size, line_spacing,
                          bold=False, italic=False, alignment=None,
                          space_after=Pt(6), space_before=Pt(0),
                          color=None, superscript_nums=False):
    """Apply formatting to a paragraph."""
    pf = para.paragraph_format
    pf.space_after = space_after
    pf.space_before = space_before
    pf.line_spacing = line_spacing
    if alignment is not None:
        pf.alignment = alignment
    for run in para.runs:
        run.font.name = font_name
        run.font.size = Pt(font_size)
        run.font.bold = bold
        run.font.italic = italic
        if color:
            run.font.color.rgb = color


def add_formatted_text(para, text, font_name, font_size, bold=False,
                       italic=False, superscript=False, color=None):
    """Add a run of formatted text to a paragraph."""
    run = para.add_run(text)
    run.font.name = font_name
    run.font.size = Pt(font_size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.superscript = superscript
    if color:
        run.font.color.rgb = color
    return run


def parse_inline_formatting(text):
    """Parse markdown inline formatting to list of (text, bold, italic, superscript) tuples."""
    segments = []
    # Handle bold+italic (***text*** or ___text___)
    # Handle bold (**text** or __text__)
    # Handle italic (*text* or _text_)
    # Handle superscript (^text^)

    pattern = re.compile(
        r'(\*\*\*(.+?)\*\*\*)'   # bold+italic
        r'|(\*\*(.+?)\*\*)'       # bold
        r'|(\*(.+?)\*)'           # italic
        r'|(\^(\d+(?:,\d+)*(?:--?\d+)?)\^)'  # superscript for references
    )

    last_end = 0
    for m in pattern.finditer(text):
        # Add text before match
        if m.start() > last_end:
            segments.append((text[last_end:m.start()], False, False, False))

        if m.group(2):  # bold+italic
            segments.append((m.group(2), True, True, False))
        elif m.group(4):  # bold
            segments.append((m.group(4), True, False, False))
        elif m.group(6):  # italic
            segments.append((m.group(6), False, True, False))
        elif m.group(8):  # superscript
            segments.append((m.group(8), False, False, True))

        last_end = m.end()

    # Add remaining text
    if last_end < len(text):
        segments.append((text[last_end:], False, False, False))

    if not segments:
        segments = [(text, False, False, False)]

    return segments


def add_markdown_paragraph(doc, text, font_name, font_size, line_spacing,
                           bold=False, italic=False, alignment=None,
                           space_after=Pt(6), space_before=Pt(0)):
    """Add a paragraph with markdown inline formatting parsed."""
    para = doc.add_paragraph()
    pf = para.paragraph_format
    pf.space_after = space_after
    pf.space_before = space_before
    pf.line_spacing = line_spacing
    if alignment:
        pf.alignment = alignment

    segments = parse_inline_formatting(text)
    for seg_text, seg_bold, seg_italic, seg_super in segments:
        add_formatted_text(para, seg_text, font_name, font_size,
                          bold=bold or seg_bold, italic=italic or seg_italic,
                          superscript=seg_super)
    return para


def parse_markdown_table(lines):
    """Parse a markdown table into headers and rows."""
    headers = []
    rows = []
    for line in lines:
        line = line.strip()
        if not line.startswith('|'):
            continue
        cells = [c.strip() for c in line.split('|')[1:-1]]
        # Skip separator rows (contain only dashes, colons, spaces)
        if all(re.match(r'^[-:]+$', c) for c in cells if c):
            continue
        if not headers:
            headers = cells
        else:
            rows.append(cells)
    return headers, rows


def add_docx_table(doc, headers, rows, font_name, font_size):
    """Add a formatted table to the document."""
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.style = 'Table Grid'
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    # Header row
    hdr_cells = table.rows[0].cells
    for i, h in enumerate(headers):
        hdr_cells[i].text = h
        for para in hdr_cells[i].paragraphs:
            para.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in para.runs:
                run.font.name = font_name
                run.font.size = Pt(font_size - 1)
                run.font.bold = True
        # Shade header
        shading = hdr_cells[i]._element.get_or_add_tcPr()
        shading_elm = shading.makeelement(qn('w:shd'), {
            qn('w:val'): 'clear',
            qn('w:color'): 'auto',
            qn('w:fill'): '2F5496'
        })
        shading.append(shading_elm)
        for para in hdr_cells[i].paragraphs:
            for run in para.runs:
                run.font.color.rgb = RGBColor(255, 255, 255)

    # Data rows
    for i, row in enumerate(rows):
        cells = table.rows[i + 1].cells
        for j, val in enumerate(row):
            if j < len(cells):
                cells[j].text = val
                for para in cells[j].paragraphs:
                    for run in para.runs:
                        run.font.name = font_name
                        run.font.size = Pt(font_size - 1)

    return table


def md_to_docx(md_path, doc, journal_cfg, is_supplement=False):
    """Convert a markdown file to DOCX content appended to doc."""
    font_name = journal_cfg["font"]
    font_size = journal_cfg["font_size"]
    line_sp = journal_cfg["line_spacing"]

    with open(md_path, 'r', encoding='utf-8') as f:
        content = f.read()

    lines = content.split('\n')
    i = 0
    table_buffer = []
    in_table = False

    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Handle tables
        if stripped.startswith('|'):
            if not in_table:
                in_table = True
                table_buffer = []
            table_buffer.append(stripped)
            i += 1
            continue
        elif in_table:
            # End of table
            headers, rows = parse_markdown_table(table_buffer)
            if headers and rows:
                add_docx_table(doc, headers, rows, font_name, font_size)
                doc.add_paragraph()  # spacing after table
            in_table = False
            table_buffer = []

        # Skip empty lines
        if not stripped:
            i += 1
            continue

        # Horizontal rule
        if stripped in ('---', '***', '___'):
            # Add a thin horizontal line via paragraph border
            para = doc.add_paragraph()
            pf = para.paragraph_format
            pf.space_before = Pt(6)
            pf.space_after = Pt(6)
            # Add bottom border
            pPr = para._element.get_or_add_pPr()
            pBdr = pPr.makeelement(qn('w:pBdr'), {})
            bottom = pBdr.makeelement(qn('w:bottom'), {
                qn('w:val'): 'single',
                qn('w:sz'): '6',
                qn('w:space'): '1',
                qn('w:color'): 'auto'
            })
            pBdr.append(bottom)
            pPr.append(pBdr)
            i += 1
            continue

        # Headings
        if stripped.startswith('#'):
            level = len(stripped) - len(stripped.lstrip('#'))
            heading_text = stripped.lstrip('#').strip()
            # Remove markdown formatting from heading
            heading_text = re.sub(r'\*\*(.+?)\*\*', r'\1', heading_text)
            heading_text = re.sub(r'\*(.+?)\*', r'\1', heading_text)

            if level == 1:
                para = doc.add_heading(heading_text, level=1)
            elif level == 2:
                para = doc.add_heading(heading_text, level=2)
            elif level == 3:
                para = doc.add_heading(heading_text, level=3)
            else:
                para = doc.add_heading(heading_text, level=min(level, 4))

            # Style the heading
            for run in para.runs:
                run.font.name = font_name
                run.font.size = Pt(font_size + (4 - level) * 2)

            i += 1
            continue

        # Regular paragraph (may span multiple lines until blank line)
        para_text = stripped
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].strip().startswith('#') \
              and not lines[i].strip().startswith('|') and not lines[i].strip() in ('---', '***', '___'):
            para_text += ' ' + lines[i].strip()
            i += 1

        add_markdown_paragraph(doc, para_text, font_name, font_size, line_sp)

    # Flush any remaining table
    if in_table and table_buffer:
        headers, rows = parse_markdown_table(table_buffer)
        if headers and rows:
            add_docx_table(doc, headers, rows, font_name, font_size)


def generate_manuscript_docx(journal_key, journal_cfg):
    """Generate a DOCX file for one journal's manuscript + figure legends."""
    journal_dir = os.path.join(BASE, journal_key)
    ms_path = os.path.join(journal_dir, "manuscript.md")
    fl_path = os.path.join(journal_dir, "figure_legends.md")

    if not os.path.exists(ms_path):
        print(f"  WARNING: {ms_path} not found, skipping")
        return None

    doc = Document()

    # Set default font
    style = doc.styles['Normal']
    font = style.font
    font.name = journal_cfg["font"]
    font.size = Pt(journal_cfg["font_size"])

    # Set narrow margins (1 inch = 914400 EMU)
    for section in doc.sections:
        section.top_margin = Cm(2.54)
        section.bottom_margin = Cm(2.54)
        section.left_margin = Cm(2.54)
        section.right_margin = Cm(2.54)

    # Convert manuscript
    print(f"  Converting manuscript.md...")
    md_to_docx(ms_path, doc, journal_cfg)

    # Page break before figure legends
    if os.path.exists(fl_path):
        doc.add_page_break()
        print(f"  Converting figure_legends.md...")
        md_to_docx(fl_path, doc, journal_cfg)

    # Add main tables inline
    table_dir = os.path.join(journal_dir, "tables")
    table_files = sorted(glob.glob(os.path.join(table_dir, "*.md")))
    if table_files:
        doc.add_page_break()
        heading = doc.add_heading("Tables", level=1)
        for run in heading.runs:
            run.font.name = journal_cfg["font"]

        for tf in table_files:
            print(f"  Converting table: {os.path.basename(tf)}")
            md_to_docx(tf, doc, journal_cfg)
            doc.add_paragraph()  # spacing between tables

    # Save
    out_path = os.path.join(journal_dir, f"{journal_cfg['short']}_Manuscript.docx")
    doc.save(out_path)
    print(f"  Saved: {out_path}")
    return out_path


# ============================================================================
# EXCEL SUPPLEMENTARY HELPERS
# ============================================================================

HEADER_FONT_XL = Font(name="Arial", bold=True, size=10, color="FFFFFF")
HEADER_FILL = PatternFill(start_color="2F5496", end_color="2F5496", fill_type="solid")
TITLE_FONT_XL = Font(name="Arial", bold=True, size=12)
BODY_FONT_XL = Font(name="Arial", size=10)
NOTE_FONT_XL = Font(name="Arial", size=9, italic=True)
THIN_BORDER = Border(
    left=Side(style="thin", color="D9D9D9"),
    right=Side(style="thin", color="D9D9D9"),
    top=Side(style="thin", color="D9D9D9"),
    bottom=Side(style="thin", color="D9D9D9"),
)
ALT_FILL = PatternFill(start_color="D6E4F0", end_color="D6E4F0", fill_type="solid")
WRAP_XL = Alignment(wrap_text=True, vertical="top")
CENTER_XL = Alignment(horizontal="center", vertical="top", wrap_text=True)


def sanitize_sheet_name(name):
    """Create a valid Excel sheet name (max 31 chars, no special chars)."""
    # Remove path and extension
    name = os.path.splitext(os.path.basename(name))[0]
    # Shorten common prefixes
    name = name.replace("Appendix_", "App_")
    name = name.replace("Table", "Tbl")
    name = name.replace("_detailed", "")
    name = name.replace("_and_", "_")
    name = name.replace("_of_", "_")
    name = name.replace("_the_", "_")
    # Remove invalid chars
    for ch in ['\\', '/', '*', '?', ':', '[', ']']:
        name = name.replace(ch, '')
    return name[:31]


def parse_md_sections(md_path):
    """Parse a markdown file into title, tables, and text blocks."""
    with open(md_path, 'r', encoding='utf-8') as f:
        content = f.read()

    lines = content.split('\n')
    title = ""
    tables = []
    text_blocks = []
    current_table_lines = []
    in_table = False
    current_text = []

    for line in lines:
        stripped = line.strip()

        # Title (first H1 heading)
        if stripped.startswith('# ') and not title:
            title = stripped.lstrip('#').strip()
            continue

        # Sub-headings
        if stripped.startswith('#'):
            if current_text:
                text_blocks.append('\n'.join(current_text))
                current_text = []
            heading = stripped.lstrip('#').strip()
            # Clean markdown formatting
            heading = re.sub(r'\*\*(.+?)\*\*', r'\1', heading)
            heading = re.sub(r'\*(.+?)\*', r'\1', heading)
            text_blocks.append(f"[HEADING] {heading}")
            continue

        # Table rows
        if stripped.startswith('|'):
            if not in_table:
                if current_text:
                    text_blocks.append('\n'.join(current_text))
                    current_text = []
                in_table = True
                current_table_lines = []
            current_table_lines.append(stripped)
            continue
        elif in_table:
            headers, rows = parse_markdown_table(current_table_lines)
            if headers:
                tables.append((headers, rows))
            in_table = False
            current_table_lines = []

        # Regular text
        if stripped:
            # Clean markdown formatting for Excel
            cleaned = re.sub(r'\*\*(.+?)\*\*', r'\1', stripped)
            cleaned = re.sub(r'\*(.+?)\*', r'\1', cleaned)
            cleaned = re.sub(r'\^(\d+)\^', r'[\1]', cleaned)
            current_text.append(cleaned)
        elif current_text:
            text_blocks.append('\n'.join(current_text))
            current_text = []

    # Flush
    if in_table and current_table_lines:
        headers, rows = parse_markdown_table(current_table_lines)
        if headers:
            tables.append((headers, rows))
    if current_text:
        text_blocks.append('\n'.join(current_text))

    return title, tables, text_blocks


def add_md_to_excel_sheet(ws, md_path):
    """Convert a markdown file into an Excel worksheet."""
    title, tables, text_blocks = parse_md_sections(md_path)

    row = 1

    # Title
    if title:
        ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=8)
        cell = ws.cell(row=row, column=1, value=title)
        cell.font = TITLE_FONT_XL
        cell.alignment = Alignment(wrap_text=True)
        row += 2

    # Write tables first (primary content)
    for headers, data_rows in tables:
        ncols = len(headers)
        # Header row
        for j, h in enumerate(headers, 1):
            cell = ws.cell(row=row, column=j, value=h)
            cell.font = HEADER_FONT_XL
            cell.fill = HEADER_FILL
            cell.alignment = CENTER_XL
            cell.border = THIN_BORDER
        row += 1

        # Data rows
        for di, dr in enumerate(data_rows):
            for j, val in enumerate(dr, 1):
                if j <= ncols:
                    # Try to convert numeric values
                    try:
                        if ',' in val and val.replace(',', '').replace('.', '').replace('-', '').isdigit():
                            val_clean = val.replace(',', '')
                            cell_val = float(val_clean) if '.' in val_clean else int(val_clean)
                        else:
                            cell_val = val
                    except (ValueError, AttributeError):
                        cell_val = val

                    cell = ws.cell(row=row, column=j, value=cell_val)
                    cell.font = BODY_FONT_XL
                    cell.alignment = WRAP_XL
                    cell.border = THIN_BORDER
                    if di % 2 == 1:
                        cell.fill = ALT_FILL
            row += 1

        row += 1  # gap between tables

    # Write text blocks as notes
    if text_blocks:
        row += 1
        for block in text_blocks:
            if block.startswith('[HEADING]'):
                heading_text = block.replace('[HEADING] ', '')
                ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=8)
                cell = ws.cell(row=row, column=1, value=heading_text)
                cell.font = Font(name="Arial", bold=True, size=11)
                cell.alignment = Alignment(wrap_text=True)
            else:
                ws.merge_cells(start_row=row, start_column=1, end_row=row, end_column=8)
                cell = ws.cell(row=row, column=1, value=block)
                cell.font = NOTE_FONT_XL
                cell.alignment = Alignment(wrap_text=True)
            row += 1

    # Set column widths
    for col in range(1, 9):
        ws.column_dimensions[get_column_letter(col)].width = 20


def generate_supplementary_excel(journal_key, journal_cfg):
    """Generate Excel workbook with all supplementary tables + appendices."""
    journal_dir = os.path.join(BASE, journal_key)
    supp_dir = os.path.join(journal_dir, "supplementary")
    table_dir = os.path.join(journal_dir, "tables")

    wb = openpyxl.Workbook()
    # Remove default sheet
    wb.remove(wb.active)

    # 1. Add main tables
    table_files = sorted(glob.glob(os.path.join(table_dir, "*.md")))
    for tf in table_files:
        sheet_name = sanitize_sheet_name(tf)
        print(f"  Adding sheet: {sheet_name}")
        ws = wb.create_sheet(title=sheet_name)
        add_md_to_excel_sheet(ws, tf)

    # 2. Add supplementary tables (Tbl1-Tbl5, Panel, Novel Findings)
    supp_tables = sorted(glob.glob(os.path.join(supp_dir, "Table*.md")))
    supp_tables += sorted(glob.glob(os.path.join(supp_dir, "Panel*.md")))
    novel = os.path.join(supp_dir, "pLIN_Novel_Findings.md")
    if os.path.exists(novel):
        supp_tables.append(novel)

    for sf in supp_tables:
        sheet_name = sanitize_sheet_name(sf)
        if sheet_name in [ws.title for ws in wb.worksheets]:
            sheet_name = "S_" + sheet_name[:29]
        print(f"  Adding sheet: {sheet_name}")
        ws = wb.create_sheet(title=sheet_name)
        add_md_to_excel_sheet(ws, sf)

    # 3. Add appendix pages
    appendix_files = sorted(glob.glob(os.path.join(supp_dir, "Appendix_*.md")))
    for af in appendix_files:
        sheet_name = sanitize_sheet_name(af)
        if sheet_name in [ws.title for ws in wb.worksheets]:
            sheet_name = sheet_name[:28] + "_2"
        print(f"  Adding sheet: {sheet_name}")
        ws = wb.create_sheet(title=sheet_name)
        add_md_to_excel_sheet(ws, af)

    # Save
    out_path = os.path.join(journal_dir, f"{journal_cfg['short']}_Supplementary_Tables.xlsx")
    wb.save(out_path)
    print(f"  Saved: {out_path}")
    return out_path


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("pLIN Manuscript DOCX & Supplementary Excel Generator")
    print("=" * 70)

    docx_files = []
    xlsx_files = []

    for journal_key, journal_cfg in JOURNALS.items():
        print(f"\n{'='*60}")
        print(f"  {journal_cfg['name']}")
        print(f"{'='*60}")

        # Generate DOCX
        print(f"\n--- Generating DOCX manuscript ---")
        docx_path = generate_manuscript_docx(journal_key, journal_cfg)
        if docx_path:
            docx_files.append(docx_path)

        # Generate Excel supplementary
        print(f"\n--- Generating Supplementary Excel ---")
        xlsx_path = generate_supplementary_excel(journal_key, journal_cfg)
        if xlsx_path:
            xlsx_files.append(xlsx_path)

    # Summary
    print(f"\n{'='*70}")
    print("GENERATION COMPLETE")
    print(f"{'='*70}")
    print(f"\nDOCX manuscripts ({len(docx_files)}):")
    for f in docx_files:
        size = os.path.getsize(f) / 1024
        print(f"  {f}  ({size:.1f} KB)")

    print(f"\nExcel supplementary ({len(xlsx_files)}):")
    for f in xlsx_files:
        size = os.path.getsize(f) / 1024
        print(f"  {f}  ({size:.1f} KB)")


if __name__ == "__main__":
    main()
