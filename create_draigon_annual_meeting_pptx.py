#!/usr/bin/env python3
"""
Generate a comprehensive PPTX slide deck for the DRAIGON Annual Consortium Meeting.
Output: output/pLIN_DRAIGON_Annual_Meeting.pptx

Target audience: DRAIGON consortium partners (clinical, bioinformatics, HTA, industry)
Duration: ~15-20 minute presentation (15 slides)
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

# ── Colours (DRAIGON / EU palette) ──────────────────────────────────────────
EU_BLUE     = RGBColor(0x00, 0x3F, 0x9A)
EU_YELLOW   = RGBColor(0xFF, 0xCC, 0x00)
DARK_BLUE   = RGBColor(0x0D, 0x47, 0xA1)
MED_BLUE    = RGBColor(0x1E, 0x88, 0xE5)
LIGHT_BLUE  = RGBColor(0xBB, 0xDE, 0xFB)
ACCENT_BLUE = RGBColor(0xE3, 0xF2, 0xFD)
NAVY        = RGBColor(0x0A, 0x2A, 0x6E)
WHITE       = RGBColor(0xFF, 0xFF, 0xFF)
BLACK       = RGBColor(0x00, 0x00, 0x00)
DARK_GRAY   = RGBColor(0x33, 0x33, 0x33)
MED_GRAY    = RGBColor(0x75, 0x75, 0x75)
LIGHT_GRAY  = RGBColor(0xF5, 0xF5, 0xF5)
GREEN       = RGBColor(0x2E, 0x7D, 0x32)
LIGHT_GREEN = RGBColor(0xC8, 0xE6, 0xC9)
ORANGE      = RGBColor(0xEF, 0x6C, 0x00)
RED         = RGBColor(0xC6, 0x28, 0x28)
TEAL        = RGBColor(0x00, 0x79, 0x6B)
PURPLE      = RGBColor(0x6A, 0x1B, 0x9A)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)


# ── Helpers ──────────────────────────────────────────────────────────────────

def add_bg(slide, color=WHITE):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_shape(slide, left, top, width, height, fill_color,
              shape_type=MSO_SHAPE.ROUNDED_RECTANGLE,
              border_color=None, border_width=Pt(1)):
    shape = slide.shapes.add_shape(shape_type, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if border_color:
        shape.line.color.rgb = border_color
        shape.line.width = border_width
    else:
        shape.line.fill.background()
    return shape


def add_text(slide, left, top, width, height, text, font_size=14, bold=False,
             color=DARK_GRAY, alignment=PP_ALIGN.LEFT, font_name="Calibri",
             italic=False):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    tf.auto_size = None
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.italic = italic
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    return txBox


def add_paragraphs(slide, left, top, width, height, items,
                   font_size=12, color=DARK_GRAY, font_name="Calibri",
                   bullet_char="▸", bullet_color=DARK_BLUE,
                   space_after=Pt(4)):
    """Add multi-line paragraphs with bullet character."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        run_b = p.add_run()
        run_b.text = f"{bullet_char}  "
        run_b.font.size = Pt(font_size)
        run_b.font.color.rgb = bullet_color
        run_b.font.bold = True
        run_b.font.name = font_name
        run_t = p.add_run()
        run_t.text = item
        run_t.font.size = Pt(font_size)
        run_t.font.color.rgb = color
        run_t.font.name = font_name
        p.space_after = space_after
        p.space_before = Pt(0)
    return txBox


def add_header_bar(slide, slide_num, total, title, subtitle=None):
    """Header bar with slide counter."""
    add_shape(slide, Inches(0), Inches(0), prs.slide_width, Inches(1.0),
              DARK_BLUE, shape_type=MSO_SHAPE.RECTANGLE)
    # Yellow accent line at bottom of header
    add_shape(slide, Inches(0), Inches(1.0), prs.slide_width, Inches(0.05),
              EU_YELLOW, shape_type=MSO_SHAPE.RECTANGLE)
    # Slide counter (top right)
    add_text(slide, Inches(11.5), Inches(0.15), Inches(1.5), Inches(0.3),
             f"Slide {slide_num} of {total}", font_size=10,
             color=LIGHT_BLUE, alignment=PP_ALIGN.RIGHT)
    # Title
    add_text(slide, Inches(0.5), Inches(0.15), Inches(11), Inches(0.5),
             title, font_size=24, bold=True, color=WHITE)
    if subtitle:
        add_text(slide, Inches(0.5), Inches(0.6), Inches(11), Inches(0.35),
                 subtitle, font_size=13, color=LIGHT_BLUE, italic=True)


def add_footer(slide):
    add_shape(slide, Inches(0), Inches(7.15), prs.slide_width, Inches(0.35),
              LIGHT_GRAY, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, Inches(0.4), Inches(7.18), Inches(6), Inches(0.3),
             "DRAIGON Annual Meeting | pLIN Work Package Update",
             font_size=9, color=MED_GRAY)
    add_text(slide, Inches(7), Inches(7.18), Inches(6), Inches(0.3),
             "UMCG | Horizon Europe GA No. 101137383",
             font_size=9, color=MED_GRAY, alignment=PP_ALIGN.RIGHT)


TOTAL_SLIDES = 15

# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 1: Title
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, NAVY)

# EU funding bar at top
add_shape(slide, Inches(0), Inches(0), prs.slide_width, Inches(0.55),
          EU_BLUE, shape_type=MSO_SHAPE.RECTANGLE)
# EU stars representation (simple)
add_text(slide, Inches(0.5), Inches(0.13), Inches(12), Inches(0.35),
         "★  Funded by the European Union  |  Horizon Europe  |  Grant Agreement No. 101137383  ★",
         font_size=12, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)

# DRAIGON badge
add_shape(slide, Inches(4.3), Inches(1.0), Inches(4.7), Inches(0.6),
          EU_YELLOW)
add_text(slide, Inches(4.3), Inches(1.05), Inches(4.7), Inches(0.5),
         "DRAIGON ANNUAL CONSORTIUM MEETING",
         font_size=16, bold=True, color=NAVY, alignment=PP_ALIGN.CENTER)

# Main title
add_text(slide, Inches(0.5), Inches(2.0), Inches(12.3), Inches(1.0),
         "pLIN: Plasmid Lineage Identification Number",
         font_size=42, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

# Subtitle
add_text(slide, Inches(0.5), Inches(3.1), Inches(12.3), Inches(0.7),
         "A Hierarchical Classification & AMR Surveillance Platform",
         font_size=22, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER, italic=True)

# Year 1 / Year 2 progress badge
add_shape(slide, Inches(4.5), Inches(4.1), Inches(4.3), Inches(0.5),
          DARK_BLUE, border_color=MED_BLUE)
add_text(slide, Inches(4.5), Inches(4.15), Inches(4.3), Inches(0.4),
         "Annual Progress Report | Year 2",
         font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

# Authors block
add_text(slide, Inches(0.5), Inches(5.0), Inches(12.3), Inches(0.5),
         "Basil Britto Xavier, Anurag Kumar Bari, Bhanu Sinha, John W.A. Rossen",
         font_size=18, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(5.5), Inches(12.3), Inches(0.4),
         "on behalf of the DRAIGON Consortium",
         font_size=14, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER, italic=True)

# Affiliation
add_text(slide, Inches(0.5), Inches(6.1), Inches(12.3), Inches(0.4),
         "Department of Medical Microbiology and Infection Prevention",
         font_size=12, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(6.45), Inches(12.3), Inches(0.4),
         "AGE Research Group | UMCG, University of Groningen, The Netherlands",
         font_size=12, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)

# Date
add_text(slide, Inches(0.5), Inches(6.95), Inches(12.3), Inches(0.4),
         "Annual Meeting 2026",
         font_size=11, color=MED_GRAY, alignment=PP_ALIGN.CENTER, italic=True)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: Outline / Agenda
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 2, TOTAL_SLIDES, "Presentation Outline",
               "What we'll cover today")
add_footer(slide)

agenda_items = [
    ("1.", "The Challenge", "Why current plasmid typing fails AMR surveillance",       DARK_BLUE),
    ("2.", "DRAIGON Alignment",   "Mapping pLIN to consortium objectives",                  TEAL),
    ("3.", "Methods Overview",    "Tetranucleotide composition + hierarchical clustering",  MED_BLUE),
    ("4.", "Reference Database",  "8,077 → 79,305 plasmids across 28 Inc/Rep groups",   PURPLE),
    ("5.", "Classifier Performance", "91.1% accuracy across Gram-neg, Gram-pos, ESKAPE",    GREEN),
    ("6.", "AMR Landscape",       "64,891 resistance gene detections mapped to lineages",   ORANGE),
    ("7.", "Outbreak Validation", "27 published studies, 13 countries, 7 mechanisms",      RED),
    ("8.", "Seven Analytical Modules", "QC, novelty, recombination, MGE, evolution",        TEAL),
    ("9.", "Clinical Deployment", "Streamlit GUI, cross-platform installers",               DARK_BLUE),
    ("10.", "Year 2 Deliverables & Next Steps", "Status, milestones, partner integration",  PURPLE),
]

for i, (num, title, desc, color) in enumerate(agenda_items):
    col = i // 5
    row = i % 5
    x = Inches(0.6) + col * Inches(6.3)
    y = Inches(1.5) + row * Inches(1.05)

    # Number circle
    add_shape(slide, x, y, Inches(0.7), Inches(0.7), color, shape_type=MSO_SHAPE.OVAL)
    add_text(slide, x, y + Inches(0.12), Inches(0.7), Inches(0.5),
             num, font_size=18, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

    # Title
    add_text(slide, x + Inches(0.85), y + Inches(0.05), Inches(5.2), Inches(0.4),
             title, font_size=16, bold=True, color=color)
    # Description
    add_text(slide, x + Inches(0.85), y + Inches(0.45), Inches(5.2), Inches(0.4),
             desc, font_size=11, color=MED_GRAY, italic=True)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: The Challenge
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 3, TOTAL_SLIDES,
               "The Challenge: Why Plasmid Typing Matters for AMR",
               "The gap that pLIN — and DRAIGON — must close")
add_footer(slide)

# Top: AMR burden statistic
add_shape(slide, Inches(0.5), Inches(1.4), Inches(12.3), Inches(1.1),
          RED, border_color=RED)
add_text(slide, Inches(0.7), Inches(1.5), Inches(12), Inches(0.4),
         "4.95 million deaths associated with bacterial AMR in 2019",
         font_size=22, bold=True, color=WHITE)
add_text(slide, Inches(0.7), Inches(1.95), Inches(12), Inches(0.5),
         "Plasmid-mediated horizontal gene transfer is the principal driver of resistance dissemination "
         "(Murray et al., Lancet 2022)",
         font_size=13, color=LIGHT_GRAY, italic=True)

# Three problems
problems = [
    ("Limited Resolution",
     "PlasmidFinder, pMLST, MOB-suite produce flat labels. 'IncN' encompasses >1,000 plasmids "
     "with vastly different resistance profiles — unable to track individual lineages.",
     ORANGE),
    ("Code Instability",
     "MOB-suite reassigns identifiers with each database update. COPLA recomputes PTUs. "
     "No tool guarantees longitudinal comparability across DRAIGON sites and time points.",
     RED),
    ("No AMR Integration",
     "No existing tool combines hierarchical classification, permanent nomenclature, "
     "and AMR gene profiling — required for clinically actionable surveillance.",
     PURPLE),
]

for i, (title, desc, color) in enumerate(problems):
    x = Inches(0.5) + i * Inches(4.27)
    y = Inches(2.8)

    add_shape(slide, x, y, Inches(4.0), Inches(2.5),
              LIGHT_GRAY, border_color=color, border_width=Pt(2))
    add_shape(slide, x, y, Inches(4.0), Inches(0.06), color,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x + Inches(0.2), y + Inches(0.15), Inches(3.6), Inches(0.5),
             title, font_size=18, bold=True, color=color)
    add_text(slide, x + Inches(0.2), y + Inches(0.75), Inches(3.6), Inches(1.7),
             desc, font_size=12, color=DARK_GRAY)

# Bottom: Clinical consequence
add_shape(slide, Inches(0.5), Inches(5.5), Inches(12.3), Inches(1.4),
          ACCENT_BLUE, border_color=DARK_BLUE)
add_text(slide, Inches(0.8), Inches(5.6), Inches(11.7), Inches(0.4),
         "Clinical consequence",
         font_size=15, bold=True, color=DARK_BLUE)
add_text(slide, Inches(0.8), Inches(6.0), Inches(11.7), Inches(0.85),
         "When a CPE outbreak occurs, infection-control teams need to know whether the resistance gene is "
         "spreading on a single plasmid lineage (clonal spread) or on multiple independent backbones "
         "(selective pressure). Without resolution, the right intervention cannot be chosen.",
         font_size=12, color=DARK_GRAY)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: DRAIGON Alignment
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 4, TOTAL_SLIDES,
               "DRAIGON Alignment: pLIN as a Core Deliverable",
               "Direct mapping to consortium objectives")
add_footer(slide)

# DRAIGON mission box
add_shape(slide, Inches(0.5), Inches(1.4), Inches(12.3), Inches(0.85),
          NAVY)
add_text(slide, Inches(0.7), Inches(1.5), Inches(12), Inches(0.35),
         "DRAIGON Mission",
         font_size=14, bold=True, color=EU_YELLOW)
add_text(slide, Inches(0.7), Inches(1.85), Inches(12), Inches(0.4),
         "Diagnose MDR infections using AI-powered genomic antibiotic susceptibility "
         "prediction from long-read sequencing data",
         font_size=12, color=WHITE, italic=True)

# 5 objectives → pLIN deliverables
mappings = [
    ("O1: Rapid AI Diagnostics", DARK_BLUE,
     "KNN classifier (91.1%) + 7 analytical modules → plasmid characterisation in <30 min"),
    ("O2: AMR Surveillance", TEAL,
     "AMRFinderPlus integration: 64,891 gene hits, lineage-linked resistance profiles"),
    ("O3: Outbreak Detection", ORANGE,
     "L6 pLIN codes (~99.9% ANI) + MLST integration to distinguish clonal vs HGT spread"),
    ("O4: Antibiotic Stewardship", GREEN,
     "3-tier risk stratification (Critical/High/Moderate) for clinical decision support"),
    ("O5: Cross-Site Deployability", PURPLE,
     "Streamlit GUI, cross-platform installers, no specialist bioinformatics required"),
]

for i, (obj, color, deliverable) in enumerate(mappings):
    y = Inches(2.5) + i * Inches(0.78)

    # Objective tag (left)
    add_shape(slide, Inches(0.5), y, Inches(0.06), Inches(0.65), color,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_shape(slide, Inches(0.56), y, Inches(3.4), Inches(0.65),
              LIGHT_GRAY, border_color=RGBColor(0xCC, 0xCC, 0xCC))
    add_text(slide, Inches(0.7), y + Inches(0.12), Inches(3.2), Inches(0.4),
             obj, font_size=13, bold=True, color=color)

    # Arrow
    add_shape(slide, Inches(4.05), y + Inches(0.18), Inches(0.45), Inches(0.3),
              color, shape_type=MSO_SHAPE.RIGHT_ARROW)

    # Deliverable (right)
    add_shape(slide, Inches(4.6), y, Inches(8.3), Inches(0.65),
              ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(0.5))
    txBox = slide.shapes.add_textbox(
        Inches(4.75), y + Inches(0.05), Inches(8.1), Inches(0.6))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    run_a = p.add_run()
    run_a.text = "✓  "
    run_a.font.size = Pt(13)
    run_a.font.color.rgb = GREEN
    run_a.font.bold = True
    run_a.font.name = "Calibri"
    run_b = p.add_run()
    run_b.text = deliverable
    run_b.font.size = Pt(12)
    run_b.font.color.rgb = DARK_GRAY
    run_b.font.name = "Calibri"

# Bottom note
add_text(slide, Inches(0.5), Inches(6.55), Inches(12.3), Inches(0.5),
         "pLIN is the genomic surveillance layer of DRAIGON's diagnostic stack — "
         "linking pathogen ID to AMR cargo to outbreak intelligence in a single open-source platform.",
         font_size=12, color=DARK_BLUE, italic=True, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: Methods Overview
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 5, TOTAL_SLIDES, "Methods: How pLIN Works",
               "Tetranucleotide composition + hierarchical single-linkage clustering")
add_footer(slide)

# Pipeline boxes (horizontal flow)
steps = [
    ("Input", "FASTA\nplasmid sequence\n(any length)", DARK_BLUE),
    ("4-mer Vectors", "256 features\nper plasmid\n(L2-normalised)", MED_BLUE),
    ("Cosine Distances", "Pairwise distances\nbetween all\nplasmids", TEAL),
    ("Single-Linkage\nClustering", "6 thresholds\ncalibrated to ANI\n(L1–L6)", GREEN),
    ("pLIN Code\nAssignment", "Permanent\nA.B.C.D.E.F\nidentifier", ORANGE),
    ("AMR / Reports", "AMRFinderPlus\n+ outbreak\ndetection", RED),
]

step_w = Inches(1.85)
step_h = Inches(1.6)
gap_x = Inches(0.18)
start_x = Inches(0.4)

for i, (title, desc, color) in enumerate(steps):
    x = start_x + i * (step_w + gap_x)
    y = Inches(1.6)

    add_shape(slide, x, y, step_w, step_h, color)
    add_text(slide, x + Inches(0.1), y + Inches(0.15), step_w - Inches(0.2), Inches(0.5),
             title, font_size=13, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.1), y + Inches(0.7), step_w - Inches(0.2), Inches(0.85),
             desc, font_size=10, color=WHITE, alignment=PP_ALIGN.CENTER)

    # Arrow between boxes
    if i < len(steps) - 1:
        arrow_x = x + step_w + Inches(0.01)
        add_shape(slide, arrow_x, y + Inches(0.65), Inches(0.16), Inches(0.3),
                  MED_GRAY, shape_type=MSO_SHAPE.RIGHT_ARROW)

# Hierarchical thresholds table
add_text(slide, Inches(0.5), Inches(3.5), Inches(6), Inches(0.4),
         "Six Hierarchical Levels Calibrated Against ANI",
         font_size=15, bold=True, color=DARK_BLUE)

threshold_data = [
    ("Level", "Distance",  "ANI",      "Resolution",          DARK_BLUE),
    ("L1",   "d ≤ 0.150", "~85%",  "Plasmid superfamily",   MED_BLUE),
    ("L2",   "d ≤ 0.100", "~90%",  "Major lineage",         MED_BLUE),
    ("L3",   "d ≤ 0.050", "~95%",  "Species-level cluster", TEAL),
    ("L4",   "d ≤ 0.020", "~98%",  "Sublineage",            TEAL),
    ("L5",   "d ≤ 0.010", "~99%",  "Clone complex",         GREEN),
    ("L6",   "d ≤ 0.001", "~99.9%", "Strain / outbreak",     RED),
]

table_x = Inches(0.5)
table_y = Inches(3.95)
col_widths = [Inches(0.7), Inches(1.4), Inches(1.0), Inches(2.5)]
row_h = Inches(0.36)

for r, row in enumerate(threshold_data):
    bg = DARK_BLUE if r == 0 else (LIGHT_GRAY if r % 2 == 0 else WHITE)
    txt = WHITE if r == 0 else DARK_GRAY
    bold = (r == 0)

    cum_x = table_x
    for c in range(4):
        add_shape(slide, cum_x, table_y + r * row_h, col_widths[c], row_h,
                  bg, shape_type=MSO_SHAPE.RECTANGLE,
                  border_color=RGBColor(0xCC, 0xCC, 0xCC), border_width=Pt(0.5))
        add_text(slide, cum_x + Inches(0.05), table_y + r * row_h + Inches(0.05),
                 col_widths[c] - Inches(0.1), row_h - Inches(0.05),
                 row[c], font_size=11, bold=bold, color=txt,
                 alignment=PP_ALIGN.CENTER if c < 3 else PP_ALIGN.LEFT)
    # Apply per-row accent (level letter colour)
    if r > 0:
        # color the level cell
        add_shape(slide, table_x, table_y + r * row_h, Inches(0.07), row_h,
                  row[4], shape_type=MSO_SHAPE.RECTANGLE)

# Right side: example pLIN code
add_text(slide, Inches(7.5), Inches(3.5), Inches(5.5), Inches(0.4),
         "Example pLIN Code",
         font_size=15, bold=True, color=DARK_BLUE)

add_shape(slide, Inches(7.5), Inches(3.95), Inches(5.3), Inches(2.5),
          NAVY)
add_text(slide, Inches(7.7), Inches(4.05), Inches(5), Inches(0.4),
         "1.1.2.15.48.671",
         font_size=28, bold=True, color=EU_YELLOW, font_name="Consolas",
         alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(7.7), Inches(4.55), Inches(5), Inches(0.3),
         "(IncN → KPC-2 carrier, n=90)",
         font_size=11, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER, italic=True)

example_lines = [
    ("L1 = 1", "Family"),
    ("L2 = 1", "Subfamily"),
    ("L3 = 2", "Cluster"),
    ("L4 = 15", "Subcluster"),
    ("L5 = 48", "Clone complex"),
    ("L6 = 671", "Strain / outbreak"),
]
for j, (code, level) in enumerate(example_lines):
    y_e = Inches(4.95) + j * Inches(0.23)
    add_text(slide, Inches(8.0), y_e, Inches(1.2), Inches(0.22),
             code, font_size=11, bold=True, color=EU_YELLOW, font_name="Consolas")
    add_text(slide, Inches(9.5), y_e, Inches(3), Inches(0.22),
             level, font_size=11, color=LIGHT_BLUE)

# Key principle
add_text(slide, Inches(0.5), Inches(6.45), Inches(12.3), Inches(0.45),
         "Key principle: codes are PERMANENT — never reassigned with database growth (LIN framework)",
         font_size=13, bold=True, color=DARK_BLUE,
         alignment=PP_ALIGN.CENTER, italic=True)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 6: Reference Database
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 6, TOTAL_SLIDES,
               "Reference Database: 79,305 Plasmids Across 28 Groups",
               "Most comprehensive plasmid classification reference to date")
add_footer(slide)

# Top: scale comparison
add_shape(slide, Inches(0.5), Inches(1.4), Inches(6), Inches(2.0),
          LIGHT_GREEN, border_color=GREEN, border_width=Pt(1))
add_text(slide, Inches(0.7), Inches(1.5), Inches(5.6), Inches(0.4),
         "Database Growth (Year 1 → Year 2)",
         font_size=15, bold=True, color=GREEN)

growth_lines = [
    ("Training set", "8,077 curated plasmids"),
    ("Reference set", "+71,249 PLSDB / NCBI RefSeq"),
    ("Combined", "79,305 unique plasmids"),
    ("pLIN codes generated", "57,886 unique L6 codes"),
    ("Classification rate", "97.3%"),
]
for j, (label, val) in enumerate(growth_lines):
    y_g = Inches(1.95) + j * Inches(0.28)
    add_text(slide, Inches(0.8), y_g, Inches(2.7), Inches(0.25),
             label, font_size=12, color=DARK_GRAY)
    add_text(slide, Inches(3.5), y_g, Inches(2.8), Inches(0.25),
             val, font_size=12, bold=True, color=DARK_BLUE)

# Right: 28 groups breakdown
add_shape(slide, Inches(6.8), Inches(1.4), Inches(6.0), Inches(2.0),
          ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(1))
add_text(slide, Inches(7.0), Inches(1.5), Inches(5.6), Inches(0.4),
         "28 Inc/Rep Groups Covered",
         font_size=15, bold=True, color=DARK_BLUE)

group_lines = [
    ("Gram-negative Inc",  "20 groups (IncF, IncN, IncX, IncHI, IncC, ColE…)"),
    ("Gram-positive rep",  "4 groups (S. aureus, Enterococcus)"),
    ("ESKAPE coverage",    "+ A. baumannii (2), P. aeruginosa (2)"),
    ("WHO priority",       "Carbapenem-resistant Acinetobacter & Pseudomonas"),
    ("New in Year 2",      "Gram-positive + non-fermenter expansion"),
]
for j, (label, val) in enumerate(group_lines):
    y_g = Inches(1.95) + j * Inches(0.28)
    add_text(slide, Inches(7.1), y_g, Inches(2.0), Inches(0.25),
             label, font_size=11, color=DARK_GRAY)
    add_text(slide, Inches(9.1), y_g, Inches(3.6), Inches(0.25),
             val, font_size=11, bold=True, color=DARK_BLUE)

# Bottom: per-group breakdown table (top 10 + summary)
add_text(slide, Inches(0.5), Inches(3.6), Inches(12.3), Inches(0.4),
         "Top Groups by Sequence Count",
         font_size=15, bold=True, color=DARK_BLUE)

top_groups = [
    ("IncFII", "4,629", "1,421"),
    ("IncN", "1,097", "431"),
    ("IncX1", "705", "420"),
    ("repEF (Enterococcus)", "192", "132"),
    ("repSA (S. aureus)", "194", "124"),
    ("IncFIB", "97", "42"),
    ("ColRNAI", "91", "38"),
    ("IncF", "75", "31"),
    ("IncX3", "56", "24"),
    ("IncHI2", "36", "18"),
    ("repAci (Acinetobacter)", "~ 280", "173"),
    ("repPae (Pseudomonas)", "~ 320", "200"),
]

t_x = Inches(0.5)
t_y = Inches(4.05)
c_w = [Inches(3.5), Inches(1.5), Inches(1.5)]
r_h = Inches(0.3)

# Header row
header = ["Inc/Rep Group", "Sequences (n)", "L6 codes (n)"]
for i, h in enumerate(header):
    cum_x = t_x + sum(c_w[:i], Inches(0))
    add_shape(slide, cum_x, t_y, c_w[i], r_h, DARK_BLUE,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, cum_x + Inches(0.05), t_y + Inches(0.04),
             c_w[i] - Inches(0.1), r_h - Inches(0.04),
             h, font_size=11, bold=True, color=WHITE)

# Data rows in two columns (6 each)
for col in range(2):
    sub = top_groups[col*6:(col+1)*6]
    for j, row in enumerate(sub):
        x_off = Inches(0) if col == 0 else Inches(6.6)
        y = t_y + r_h + j * r_h
        for i, val in enumerate(row):
            cum_x = t_x + x_off + sum(c_w[:i], Inches(0))
            bg = LIGHT_GRAY if j % 2 == 0 else WHITE
            add_shape(slide, cum_x, y, c_w[i], r_h, bg,
                      shape_type=MSO_SHAPE.RECTANGLE,
                      border_color=RGBColor(0xCC, 0xCC, 0xCC), border_width=Pt(0.3))
            align = PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER
            add_text(slide, cum_x + Inches(0.05), y + Inches(0.04),
                     c_w[i] - Inches(0.1), r_h - Inches(0.04),
                     val, font_size=10, color=DARK_GRAY, alignment=align)

# Headers for second column
for i, h in enumerate(header):
    cum_x = t_x + Inches(6.6) + sum(c_w[:i], Inches(0))
    add_shape(slide, cum_x, t_y, c_w[i], r_h, DARK_BLUE,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, cum_x + Inches(0.05), t_y + Inches(0.04),
             c_w[i] - Inches(0.1), r_h - Inches(0.04),
             h, font_size=11, bold=True, color=WHITE)

# Note
add_text(slide, Inches(0.5), Inches(6.7), Inches(12.3), Inches(0.3),
         "Source: PLSDB 2025 + NCBI RefSeq | GC-content window widened (25–70%) for Gram-positive coverage",
         font_size=10, color=MED_GRAY, italic=True, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 7: Classifier Performance
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 7, TOTAL_SLIDES,
               "Classifier Performance: 91.1% Accuracy Across 28 Groups",
               "Five-fold cross-validation + independent ML benchmark")
add_footer(slide)

# Big stat
add_shape(slide, Inches(0.5), Inches(1.4), Inches(4), Inches(2.5),
          DARK_BLUE)
add_text(slide, Inches(0.5), Inches(1.55), Inches(4), Inches(0.4),
         "Headline accuracy",
         font_size=14, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(1.95), Inches(4), Inches(1.2),
         "91.1%",
         font_size=72, bold=True, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(3.15), Inches(4), Inches(0.5),
         "KNN k=5, cosine metric\ndistance-weighted\n5-fold stratified CV",
         font_size=12, color=WHITE, alignment=PP_ALIGN.CENTER)

# ML benchmark table
add_text(slide, Inches(4.8), Inches(1.4), Inches(8), Inches(0.4),
         "Independent ML Benchmark (Inc-group prediction)",
         font_size=14, bold=True, color=DARK_BLUE)

ml_data = [
    ["Model",            "Weighted F1", "SD",    "Min",   "Max"],
    ["XGBoost",          "0.896",       "0.011", "0.885", "0.907"],
    ["Gradient Boosting","0.893",       "0.009", "0.884", "0.902"],
    ["Random Forest",    "0.874",       "0.021", "—","—"],
    ["Logistic Reg.",    "0.866",       "0.012", "—","—"],
]

ml_x = Inches(4.8)
ml_y = Inches(1.85)
ml_cw = [Inches(2.0), Inches(1.5), Inches(1.0), Inches(1.0), Inches(1.0)]
ml_rh = Inches(0.35)

for r, row in enumerate(ml_data):
    bg = DARK_BLUE if r == 0 else (LIGHT_GRAY if r % 2 == 0 else WHITE)
    txt = WHITE if r == 0 else DARK_GRAY
    bold_v = (r == 0)
    for c, val in enumerate(row):
        cum_x = ml_x + sum(ml_cw[:c], Inches(0))
        add_shape(slide, cum_x, ml_y + r * ml_rh, ml_cw[c], ml_rh,
                  bg, shape_type=MSO_SHAPE.RECTANGLE,
                  border_color=RGBColor(0xCC, 0xCC, 0xCC), border_width=Pt(0.3))
        align = PP_ALIGN.LEFT if c == 0 else PP_ALIGN.CENTER
        add_text(slide, cum_x + Inches(0.05), ml_y + r * ml_rh + Inches(0.05),
                 ml_cw[c] - Inches(0.1), ml_rh - Inches(0.05),
                 val, font_size=11, bold=bold_v, color=txt, alignment=align)

# Per-group highlights
add_text(slide, Inches(0.5), Inches(4.2), Inches(12.3), Inches(0.4),
         "Per-Class Accuracy: New Gram-Positive & ESKAPE Groups",
         font_size=15, bold=True, color=DARK_BLUE)

groups_perf = [
    ("repSA_large", "S. aureus large", "93.4%", GREEN),
    ("repSA_small", "S. aureus small", "95.9%", GREEN),
    ("repEF_conj",  "Enterococcus conj.", "90.2%", GREEN),
    ("repEF_res",   "Enterococcus resist.", "88.4%", ORANGE),
    ("repAci",      "Acinetobacter", "Pending FastANI", MED_GRAY),
    ("repPae",      "Pseudomonas", "Pending FastANI", MED_GRAY),
]

for i, (code, name, acc, color) in enumerate(groups_perf):
    col = i // 3
    row = i % 3
    x = Inches(0.5) + col * Inches(6.3)
    y = Inches(4.7) + row * Inches(0.55)

    add_shape(slide, x, y, Inches(6.0), Inches(0.45),
              LIGHT_GRAY, border_color=RGBColor(0xCC, 0xCC, 0xCC))
    add_shape(slide, x, y, Inches(0.06), Inches(0.45), color,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x + Inches(0.2), y + Inches(0.08),
             Inches(1.8), Inches(0.3),
             code, font_size=12, bold=True, color=DARK_GRAY)
    add_text(slide, x + Inches(2.0), y + Inches(0.08),
             Inches(2.5), Inches(0.3),
             name, font_size=11, color=MED_GRAY, italic=True)
    add_text(slide, x + Inches(4.5), y + Inches(0.08),
             Inches(1.5), Inches(0.3),
             acc, font_size=12, bold=True, color=color, alignment=PP_ALIGN.RIGHT)

# Footer note
add_text(slide, Inches(0.5), Inches(6.55), Inches(12.3), Inches(0.5),
         "Modest reduction from 92.2% (Year 1, 20 groups) to 91.1% (Year 2, 28 groups) reflects added "
         "complexity of low-GC Gram-positive groups — expected and acceptable.",
         font_size=11, color=MED_GRAY, italic=True, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 8: AMR Landscape
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 8, TOTAL_SLIDES,
               "AMR Landscape: 64,891 Resistance Gene Detections",
               "AMRFinderPlus integration mapped to plasmid lineages")
add_footer(slide)

# Stats grid - top row
amr_stats = [
    ("64,891", "Total gene\nhits",          DARK_BLUE),
    ("83.1%",  "Plasmids carrying\nAMR genes", TEAL),
    ("1,635",  "Carbapenemase\ngenes",        RED),
    ("1,804",  "ESBL\ngenes",                 ORANGE),
    ("204",    "Mobile colistin\n(mcr)",       PURPLE),
    ("2,315",  "PMQR\n(quinolone)",           GREEN),
]

for i, (num, label, color) in enumerate(amr_stats):
    x = Inches(0.5) + i * Inches(2.13)
    y = Inches(1.4)
    add_shape(slide, x, y, Inches(2.0), Inches(1.5),
              LIGHT_GRAY, border_color=color, border_width=Pt(1.5))
    add_shape(slide, x, y, Inches(2.0), Inches(0.06),
              color, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x, y + Inches(0.2), Inches(2.0), Inches(0.6),
             num, font_size=24, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.1), y + Inches(0.85),
             Inches(1.8), Inches(0.6),
             label, font_size=11, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Top resistance genes
add_text(slide, Inches(0.5), Inches(3.1), Inches(6), Inches(0.4),
         "Top Carbapenemase Detections",
         font_size=14, bold=True, color=DARK_BLUE)

carba = [
    ("blaKPC-2",  "824", 824),
    ("blaNDM-1",  "228", 228),
    ("blaKPC-3",  "193", 193),
    ("blaNDM-5",  "89",   89),
]
max_v = max(v for _, _, v in carba)
for i, (gene, val, v) in enumerate(carba):
    y = Inches(3.55) + i * Inches(0.35)
    add_text(slide, Inches(0.5), y, Inches(1.5), Inches(0.3),
             gene, font_size=11, bold=True, italic=True, color=DARK_GRAY)
    bar_w = Inches(3.5 * v / max_v)
    add_shape(slide, Inches(2.0), y + Inches(0.05), bar_w, Inches(0.2),
              RED, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, Inches(5.6), y, Inches(1), Inches(0.3),
             val, font_size=11, bold=True, color=RED)

# ESBL
add_text(slide, Inches(7.0), Inches(3.1), Inches(6), Inches(0.4),
         "Top ESBL & Colistin Detections",
         font_size=14, bold=True, color=DARK_BLUE)

esbl = [
    ("blaCTX-M-15", "505", 505, ORANGE),
    ("blaCTX-M-65", "319", 319, ORANGE),
    ("blaSHV-12",   "277", 277, ORANGE),
    ("mcr-1.1",     "83",   83, PURPLE),
]
max_v2 = max(v for _, _, v, _ in esbl)
for i, (gene, val, v, c) in enumerate(esbl):
    y = Inches(3.55) + i * Inches(0.35)
    add_text(slide, Inches(7.0), y, Inches(1.7), Inches(0.3),
             gene, font_size=11, bold=True, italic=True, color=DARK_GRAY)
    bar_w = Inches(2.7 * v / max_v2)
    add_shape(slide, Inches(8.7), y + Inches(0.05), bar_w, Inches(0.2),
              c, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, Inches(11.7), y, Inches(1), Inches(0.3),
             val, font_size=11, bold=True, color=c)

# Bottom - lineage-linked AMR
add_shape(slide, Inches(0.5), Inches(5.2), Inches(12.3), Inches(1.7),
          ACCENT_BLUE, border_color=DARK_BLUE)
add_text(slide, Inches(0.7), Inches(5.3), Inches(12), Inches(0.35),
         "Why this matters: Lineage-Linked AMR Surveillance",
         font_size=15, bold=True, color=DARK_BLUE)

linkage_items = [
    "First system to map every resistance gene detection to a permanent, hierarchical plasmid code",
    "Enables prospective tracking of high-risk lineages across DRAIGON sites and time points",
    "Distinguishes 'KPC-2 in IncN' (one specific lineage) from 'KPC-2 in 50 different IncN backbones'",
    "Direct input for DRAIGON’s antibiotic stewardship objective — right drug, right dose, right time",
]
add_paragraphs(slide, Inches(0.8), Inches(5.7), Inches(11.7), Inches(1.2),
               linkage_items, font_size=11, color=DARK_GRAY,
               bullet_color=DARK_BLUE, space_after=Pt(2))


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 9: Outbreak Validation
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 9, TOTAL_SLIDES,
               "Outbreak Validation: 27 Studies, 13 Countries",
               "Cross-validation against published surveillance data")
add_footer(slide)

# Top stats row
out_stats = [
    ("74",  "Outbreak\nplasmids",     DARK_BLUE),
    ("27",  "Published\nstudies",     TEAL),
    ("13",  "Countries\n(4 continents)", ORANGE),
    ("7",   "Resistance\nmechanisms", RED),
    ("85.1%", "High-confidence\nclassification", GREEN),
    ("9",   "Intra-study\nclusters detected", PURPLE),
]
for i, (num, label, color) in enumerate(out_stats):
    x = Inches(0.5) + i * Inches(2.13)
    y = Inches(1.4)
    add_shape(slide, x, y, Inches(2.0), Inches(1.4),
              color)
    add_text(slide, x, y + Inches(0.15), Inches(2.0), Inches(0.6),
             num, font_size=24, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.1), y + Inches(0.8),
             Inches(1.8), Inches(0.55),
             label, font_size=11, color=WHITE, alignment=PP_ALIGN.CENTER)

# Two high-risk lineages
add_text(slide, Inches(0.5), Inches(3.1), Inches(12.3), Inches(0.4),
         "Example: High-Risk Lineages Detected by pLIN",
         font_size=15, bold=True, color=DARK_BLUE)

lineages = [
    ("pLIN 671 (IncN)",
     RED,
     [("Members", "n = 90"),
      ("blaKPC-2 carriage", "100%"),
      ("Mean AMR genes", "13.2"),
      ("Risk tier", "CRITICAL")],
     "Dominant KPC-2 lineage spanning multiple countries — invisible under conventional 'IncN' typing."),
    ("pLIN 860 (multi-Inc)",
     ORANGE,
     [("Members", "n = 142, 5 Inc groups"),
      ("mcr carriage", "44.4%"),
      ("Mean AMR genes", "14.4"),
      ("Risk tier", "CRITICAL")],
     "Cross-Inc lineage with mobile colistin resistance — demonstrates MDR-stacking on a single backbone."),
]

for idx, (name, color, stats, note) in enumerate(lineages):
    x_b = Inches(0.5) + idx * Inches(6.4)
    y_b = Inches(3.6)
    add_shape(slide, x_b, y_b, Inches(6.0), Inches(2.7),
              LIGHT_GRAY, border_color=color, border_width=Pt(1.5))
    add_shape(slide, x_b, y_b, Inches(6.0), Inches(0.06),
              color, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x_b + Inches(0.2), y_b + Inches(0.15),
             Inches(5.6), Inches(0.4),
             name, font_size=18, bold=True, color=color)
    for j, (k, v) in enumerate(stats):
        y_s = y_b + Inches(0.65) + j * Inches(0.32)
        add_text(slide, x_b + Inches(0.3), y_s, Inches(2.5), Inches(0.28),
                 k + ":", font_size=11, bold=True, color=DARK_GRAY)
        add_text(slide, x_b + Inches(2.8), y_s, Inches(3), Inches(0.28),
                 v, font_size=11, bold=True, color=color)
    add_text(slide, x_b + Inches(0.3), y_b + Inches(2.0),
             Inches(5.5), Inches(0.65),
             note, font_size=10, color=MED_GRAY, italic=True)

# Bottom: combined chromosome-plasmid typing
add_shape(slide, Inches(0.5), Inches(6.45), Inches(12.3), Inches(0.55),
          NAVY)
add_text(slide, Inches(0.7), Inches(6.5), Inches(12), Inches(0.45),
         "Combined MLST + pLIN typing distinguishes clonal spread vs. horizontal plasmid transfer "
         "— directly answers infection-control questions in real outbreaks.",
         font_size=12, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 10: Seven Analytical Modules
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 10, TOTAL_SLIDES,
               "Seven Analytical Modules: Beyond Classification",
               "Year 2 expansion addressing surveillance limitations")
add_footer(slide)

modules = [
    ("L3", "Assembly\nCompleteness",
     "Composite 0–100 score\n5 metrics (N50, contigs, gaps...)",
     "94.2% correct categorisation", DARK_BLUE),
    ("L4", "Database\nCoverage",
     "Nearest-neighbour distance\npercentile, traffic-light alert",
     "GREEN/YELLOW/RED tiers", TEAL),
    ("L5", "Recombination\nDetection",
     "minimap2 PAF analysis,\nblock fragmentation index",
     "87.3% sens / 94.1% spec", ORANGE),
    ("L6", "Novel Group\nDiscovery",
     "Cluster low-confidence\nplasmids at L3 threshold",
     "18 putative novel groups", PURPLE),
    ("L7", "Evolutionary\nRate Estimation",
     "SNP × time regression\nwithin L6 clusters",
     "1.8–8.4 × 10⁻⁶ subs/site/yr", GREEN),
    ("L8", "Cluster\nStability",
     "50-iter bootstrap +\nlinkage method comparison",
     "ARI > 0.85 in 20/28 grps", RED),
    ("L10", "MGE Boundary\nDetection",
     "IS elements, integrases,\ncomposite transposon scan",
     "94.7% IS sensitivity", NAVY),
]

# 4 + 3 layout
for i, (limit, name, desc, perf, color) in enumerate(modules):
    col = i % 4
    row = i // 4
    x = Inches(0.5) + col * Inches(3.2)
    y = Inches(1.4) + row * Inches(2.6)

    add_shape(slide, x, y, Inches(3.0), Inches(2.45),
              LIGHT_GRAY, border_color=color, border_width=Pt(1.5))
    # color top
    add_shape(slide, x, y, Inches(3.0), Inches(0.5),
              color, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x + Inches(0.1), y + Inches(0.05),
             Inches(0.7), Inches(0.4),
             limit, font_size=18, bold=True, color=EU_YELLOW)
    add_text(slide, x + Inches(0.8), y + Inches(0.05),
             Inches(2.1), Inches(0.4),
             name, font_size=12, bold=True, color=WHITE)
    add_text(slide, x + Inches(0.15), y + Inches(0.6),
             Inches(2.7), Inches(1.0),
             desc, font_size=11, color=DARK_GRAY)
    # Performance badge
    add_shape(slide, x + Inches(0.15), y + Inches(1.75),
              Inches(2.7), Inches(0.55),
              color)
    add_text(slide, x + Inches(0.2), y + Inches(1.82),
             Inches(2.6), Inches(0.45),
             perf, font_size=11, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)

# Note on right side (8th cell)
x_n = Inches(0.5) + 3 * Inches(3.2)
y_n = Inches(1.4) + 1 * Inches(2.6)
add_shape(slide, x_n, y_n, Inches(3.0), Inches(2.45),
          NAVY)
add_text(slide, x_n + Inches(0.2), y_n + Inches(0.2),
         Inches(2.6), Inches(0.4),
         "DRAIGON Impact",
         font_size=14, bold=True, color=EU_YELLOW)
note_items = [
    "QC before classification",
    "Novelty alerts for new pathogens",
    "Detect mosaic/chimeric plasmids",
    "Time-resolved phylodynamics",
    "Robust thresholds across sites",
]
add_paragraphs(slide, x_n + Inches(0.2), y_n + Inches(0.7),
               Inches(2.7), Inches(1.7),
               note_items, font_size=11, color=WHITE,
               bullet_char="•", bullet_color=EU_YELLOW,
               space_after=Pt(4))


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 11: Clinical Deployment / Software
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 11, TOTAL_SLIDES,
               "Clinical Deployment: Streamlit GUI + Cross-Platform Installers",
               "Ready for use at DRAIGON clinical partner sites")
add_footer(slide)

# Three columns: Software stack, Deployment, User workflow
col_w = Inches(4.0)
col_h = Inches(5.1)
gap = Inches(0.2)
start = Inches(0.5)

cols = [
    ("Software Stack", DARK_BLUE, [
        "Python 3.10+ Streamlit GUI",
        "BioPython for FASTA handling",
        "scikit-learn KNN classifier",
        "AMRFinderPlus v4.2.5 integration",
        "MOB-suite for mobility",
        "minimap2 for recombination",
        "MinCED for CRISPR analysis",
        "matplotlib + plotly visualisation",
    ]),
    ("Cross-Platform Deployment", TEAL, [
        "macOS one-click installer (.sh)",
        "Linux installer with conda env",
        "Windows .bat installer",
        "Single command: streamlit run plin_app.py",
        "Standard hardware (16 GB RAM)",
        "<30 min full pipeline runtime",
        "No specialist bioinformatics required",
        "Open-source GPL-3.0 on GitHub",
    ]),
    ("User Workflow at Clinical Sites", ORANGE, [
        "Upload FASTA / multi-FASTA",
        "Auto-detect plasmid contigs",
        "Inc/Rep group classification",
        "Permanent pLIN code assigned",
        "AMR + virulence + stress profiling",
        "Clinical risk tier (3-level)",
        "Outbreak detection (within batch)",
        "Downloadable PDF / CSV reports",
    ]),
]

for i, (title, color, items) in enumerate(cols):
    x = start + i * (col_w + gap)
    y = Inches(1.4)
    add_shape(slide, x, y, col_w, col_h, LIGHT_GRAY,
              border_color=color, border_width=Pt(1.5))
    add_shape(slide, x, y, col_w, Inches(0.55), color,
              shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, x + Inches(0.2), y + Inches(0.1),
             col_w - Inches(0.4), Inches(0.4),
             title, font_size=16, bold=True, color=WHITE)
    add_paragraphs(slide, x + Inches(0.25), y + Inches(0.7),
                   col_w - Inches(0.5), col_h - Inches(0.85),
                   items, font_size=12, color=DARK_GRAY,
                   bullet_color=color, space_after=Pt(6))

# Bottom: GitHub link
add_shape(slide, Inches(0.5), Inches(6.6), Inches(12.3), Inches(0.5),
          NAVY)
add_text(slide, Inches(0.7), Inches(6.65), Inches(12), Inches(0.4),
         "github.com/xavierbasilbritto-hub/pLIN-plasmid-classification  |  "
         "Open-source under GPL-3.0  |  Free for DRAIGON partner deployment",
         font_size=12, bold=True, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 12: Year 2 Deliverables
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 12, TOTAL_SLIDES,
               "Year 2 Deliverables: Status Report",
               "What was promised, what was delivered")
add_footer(slide)

deliverables = [
    ("D1", "Expanded reference database",
     "From 6,998 (3 Inc groups) to 79,305 plasmids (28 groups)",
     "DELIVERED", GREEN),
    ("D2", "Gram-positive coverage",
     "S. aureus + Enterococcus rep typing added (4 groups)",
     "DELIVERED", GREEN),
    ("D3", "ESKAPE expansion",
     "A. baumannii + P. aeruginosa rep types (4 groups)",
     "DELIVERED", GREEN),
    ("D4", "Manuscript submission",
     "Submitted to Lancet Microbe (also CMI, Microbial Genomics, BIB)",
     "SUBMITTED", GREEN),
    ("D5", "Seven analytical modules",
     "Completeness, coverage, recombination, novelty, evolution, stability, MGE",
     "DELIVERED", GREEN),
    ("D6", "Cross-platform installers",
     "macOS, Linux, Windows one-click setup",
     "DELIVERED", GREEN),
    ("D7", "Outbreak validation dataset",
     "27 published studies, 13 countries, 7 resistance mechanisms",
     "DELIVERED", GREEN),
    ("D8", "Clinical risk stratification",
     "3-tier system (Critical/High/Moderate) for AMR carriage",
     "DELIVERED", GREEN),
    ("D9", "FastANI validation Gram-pos/ESKAPE",
     "Currently running for 4 new group families",
     "IN PROGRESS", ORANGE),
    ("D10", "Long-read integration",
     "Nanopore-specific QC and contig handling tuning",
     "PLANNED Y3", MED_GRAY),
]

dx = Inches(0.5)
dy = Inches(1.4)
dw = Inches(12.3)
dh = Inches(0.5)

# Header
add_shape(slide, dx, dy, dw, dh, DARK_BLUE,
          shape_type=MSO_SHAPE.RECTANGLE)
hdr = ["ID", "Deliverable", "Description", "Status"]
hdr_w = [Inches(0.7), Inches(3.5), Inches(6.5), Inches(1.6)]
for j, h in enumerate(hdr):
    cum_x = dx + sum(hdr_w[:j], Inches(0))
    add_text(slide, cum_x + Inches(0.1), dy + Inches(0.1),
             hdr_w[j] - Inches(0.2), Inches(0.3),
             h, font_size=12, bold=True, color=WHITE)

# Rows
for i, (id_, name, desc, status, color) in enumerate(deliverables):
    y = dy + dh + i * Inches(0.45)
    bg = LIGHT_GRAY if i % 2 == 0 else WHITE

    # Row background
    add_shape(slide, dx, y, dw, Inches(0.45), bg,
              shape_type=MSO_SHAPE.RECTANGLE,
              border_color=RGBColor(0xCC, 0xCC, 0xCC), border_width=Pt(0.3))

    # Cells
    add_text(slide, dx + Inches(0.1), y + Inches(0.1),
             hdr_w[0] - Inches(0.2), Inches(0.3),
             id_, font_size=11, bold=True, color=DARK_BLUE)
    add_text(slide, dx + hdr_w[0] + Inches(0.1), y + Inches(0.1),
             hdr_w[1] - Inches(0.2), Inches(0.3),
             name, font_size=11, bold=True, color=DARK_GRAY)
    add_text(slide, dx + hdr_w[0] + hdr_w[1] + Inches(0.1), y + Inches(0.1),
             hdr_w[2] - Inches(0.2), Inches(0.3),
             desc, font_size=10, color=MED_GRAY)
    # Status badge
    sx = dx + hdr_w[0] + hdr_w[1] + hdr_w[2] + Inches(0.1)
    add_shape(slide, sx, y + Inches(0.07), Inches(1.4), Inches(0.31),
              color)
    add_text(slide, sx, y + Inches(0.1), Inches(1.4), Inches(0.25),
             status, font_size=10, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 13: Year 3 Roadmap
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 13, TOTAL_SLIDES,
               "Year 3 Roadmap: Next Steps",
               "Integration with DRAIGON partner workflows")
add_footer(slide)

roadmap = [
    ("Q1", "Long-Read Integration",
     "Nanopore-specific QC, polishing-aware contig handling, raw-signal compatibility",
     DARK_BLUE),
    ("Q1", "Clinical Site Pilots",
     "Deploy at UMCG (Groningen), Isala (Zwolle), OSS (Vienna) — BSI & PJI cohorts",
     TEAL),
    ("Q2", "FastANI Gram-Positive / ESKAPE",
     "Complete external validation for 8 new groups",
     ORANGE),
    ("Q2", "Mayo & Johns Hopkins integration",
     "US-side validation cohort, comparison with local typing pipelines",
     PURPLE),
    ("Q3", "MIC Prediction Bridge",
     "Connect pLIN lineage → expected resistance phenotype → MIC prediction",
     GREEN),
    ("Q3", "Albania Pilot (UHSN)",
     "Validation in medium-income setting with high AMR burden",
     RED),
    ("Q4", "Health Technology Assessment",
     "Health-Ecore HTA — economic + clinical impact analysis",
     NAVY),
    ("Q4", "v3.0 Release",
     "Public DRAIGON-validated release with all partner data",
     EU_BLUE),
]

for i, (q, name, desc, color) in enumerate(roadmap):
    col = i % 2
    row = i // 2
    x = Inches(0.5) + col * Inches(6.3)
    y = Inches(1.4) + row * Inches(1.35)

    # Quarter badge
    add_shape(slide, x, y, Inches(0.7), Inches(1.2), color)
    add_text(slide, x, y + Inches(0.4), Inches(0.7), Inches(0.4),
             q, font_size=22, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)

    # Content box
    add_shape(slide, x + Inches(0.7), y, Inches(5.3), Inches(1.2),
              LIGHT_GRAY, border_color=color, border_width=Pt(1))
    add_text(slide, x + Inches(0.85), y + Inches(0.1),
             Inches(5), Inches(0.4),
             name, font_size=14, bold=True, color=color)
    add_text(slide, x + Inches(0.85), y + Inches(0.5),
             Inches(5), Inches(0.65),
             desc, font_size=11, color=DARK_GRAY)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 14: Partner Asks / Discussion Points
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, 14, TOTAL_SLIDES,
               "Partner Asks: How You Can Help",
               "Concrete actions for DRAIGON consortium members")
add_footer(slide)

asks = [
    {
        "partner": "Clinical Sites\n(UMCG, Isala, OSS, UHSN, Mayo, JHU)",
        "color": DARK_BLUE,
        "asks": [
            "Share BSI / PJI plasmid sequences for prospective validation",
            "Pilot pLIN GUI in routine workflow — collect usability feedback",
            "Provide clinical metadata (resistance phenotype, patient outcome)",
            "Identify high-priority plasmid lineages of local concern",
        ]
    },
    {
        "partner": "Bioinformatics & Industry\n(Sandoz, Camtech, Health-Ecore)",
        "color": TEAL,
        "asks": [
            "Benchmark pLIN against in-house typing tools",
            "Integrate pLIN API into existing diagnostic platforms",
            "Co-develop MIC prediction layer using lineage features",
            "HTA modelling: cost-effectiveness of lineage surveillance",
        ]
    },
    {
        "partner": "Coordination & WP Leads\n(EVI, all WP leads)",
        "color": ORANGE,
        "asks": [
            "Endorse pLIN as DRAIGON's reference plasmid classifier",
            "Support manuscript dissemination at upcoming conferences",
            "Coordinate cross-site data-sharing agreements",
            "Plan joint workshop / training session at next meeting",
        ]
    },
]

for i, ask_block in enumerate(asks):
    x = Inches(0.5) + i * Inches(4.27)
    y = Inches(1.4)
    add_shape(slide, x, y, Inches(4.0), Inches(5.2),
              LIGHT_GRAY, border_color=ask_block["color"], border_width=Pt(1.5))
    # Header
    add_shape(slide, x, y, Inches(4.0), Inches(1.0),
              ask_block["color"])
    add_text(slide, x + Inches(0.2), y + Inches(0.15),
             Inches(3.7), Inches(0.8),
             ask_block["partner"], font_size=14, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)

    # Items
    add_paragraphs(slide, x + Inches(0.25), y + Inches(1.15),
                   Inches(3.5), Inches(3.9),
                   ask_block["asks"], font_size=12, color=DARK_GRAY,
                   bullet_color=ask_block["color"], bullet_char="→",
                   space_after=Pt(8))

# Bottom CTA
add_shape(slide, Inches(0.5), Inches(6.7), Inches(12.3), Inches(0.4),
          NAVY)
add_text(slide, Inches(0.5), Inches(6.75), Inches(12.3), Inches(0.35),
         "Contact: basilbritto.xavier@umcg.nl  |  Repository:  github.com/xavierbasilbritto-hub/pLIN-plasmid-classification",
         font_size=12, bold=True, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 15: Acknowledgements + Q&A
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, NAVY)

# EU bar
add_shape(slide, Inches(0), Inches(0), prs.slide_width, Inches(0.5),
          EU_BLUE, shape_type=MSO_SHAPE.RECTANGLE)
add_text(slide, Inches(0.5), Inches(0.1), Inches(12), Inches(0.35),
         "★  Funded by the European Union  |  Horizon Europe  |  Grant Agreement No. 101137383  ★",
         font_size=11, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)

# Big Thank You
add_text(slide, Inches(0.5), Inches(0.9), Inches(12.3), Inches(1.0),
         "Thank You",
         font_size=60, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(2.0), Inches(12.3), Inches(0.5),
         "Questions & Discussion",
         font_size=22, color=EU_YELLOW, alignment=PP_ALIGN.CENTER, italic=True)

# Acknowledgements
add_shape(slide, Inches(1.0), Inches(2.9), Inches(11.3), Inches(2.7),
          DARK_BLUE, border_color=MED_BLUE)

add_text(slide, Inches(1.2), Inches(3.0), Inches(11), Inches(0.4),
         "Acknowledgements",
         font_size=18, bold=True, color=EU_YELLOW)

ack_items = [
    "DRAIGON Consortium partners and Work Package leads",
    "European Vaccine Initiative (EVI) — Project Coordination",
    "AGE Research Group, UMCG — Department of Medical Microbiology and Infection Prevention",
    "PLSDB and NCBI RefSeq — sequence data submitters worldwide",
    "Authors of the 27 published outbreak studies that enabled validation",
]
add_paragraphs(slide, Inches(1.3), Inches(3.5), Inches(10.8), Inches(2.0),
               ack_items, font_size=13, color=WHITE,
               bullet_color=EU_YELLOW, bullet_char="•",
               space_after=Pt(6))

# Disclaimer
add_text(slide, Inches(1.0), Inches(5.8), Inches(11.3), Inches(0.6),
         "Views and opinions expressed are those of the authors only and do not necessarily reflect "
         "those of the European Union or HADEA. Neither the European Union nor the granting authority "
         "can be held responsible for them.",
         font_size=10, color=MED_GRAY, alignment=PP_ALIGN.CENTER, italic=True)

# Contact info
add_text(slide, Inches(0.5), Inches(6.6), Inches(12.3), Inches(0.4),
         "basilbritto.xavier@umcg.nl  |  github.com/xavierbasilbritto-hub/pLIN-plasmid-classification",
         font_size=14, bold=True, color=EU_YELLOW, alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(0.5), Inches(7.05), Inches(12.3), Inches(0.4),
         "UMCG | University of Groningen | The Netherlands",
         font_size=12, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════
out_path = os.path.join(OUT_DIR, "pLIN_DRAIGON_Annual_Meeting.pptx")
prs.save(out_path)
print(f"Slide deck saved to: {out_path}")
print(f"Total slides: {len(prs.slides)}")
