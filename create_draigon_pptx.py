#!/usr/bin/env python3
"""
Generate a short PPTX slide deck showing how pLIN fulfils DRAIGON objectives.
Output: output/pLIN_DRAIGON_Objectives.pptx
"""

import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE

# ── Colours ──────────────────────────────────────────────────────────────────
DARK_BLUE   = RGBColor(0x0D, 0x47, 0xA1)
MED_BLUE    = RGBColor(0x1E, 0x88, 0xE5)
LIGHT_BLUE  = RGBColor(0xBB, 0xDE, 0xFB)
ACCENT_BLUE = RGBColor(0xE3, 0xF2, 0xFD)
WHITE       = RGBColor(0xFF, 0xFF, 0xFF)
BLACK       = RGBColor(0x00, 0x00, 0x00)
DARK_GRAY   = RGBColor(0x33, 0x33, 0x33)
MED_GRAY    = RGBColor(0x75, 0x75, 0x75)
LIGHT_GRAY  = RGBColor(0xF5, 0xF5, 0xF5)
GREEN       = RGBColor(0x2E, 0x7D, 0x32)
ORANGE      = RGBColor(0xEF, 0x6C, 0x00)
RED         = RGBColor(0xC6, 0x28, 0x28)
TEAL        = RGBColor(0x00, 0x79, 0x6B)
PURPLE      = RGBColor(0x6A, 0x1B, 0x9A)
EU_BLUE     = RGBColor(0x00, 0x3F, 0x9A)  # EU flag blue

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR  = os.path.join(BASE_DIR, "output")
os.makedirs(OUT_DIR, exist_ok=True)

prs = Presentation()
prs.slide_width  = Inches(13.333)
prs.slide_height = Inches(7.5)

# ── Helper functions ─────────────────────────────────────────────────────────

def add_bg(slide, color=WHITE):
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_shape(slide, left, top, width, height, fill_color,
              shape_type=MSO_SHAPE.ROUNDED_RECTANGLE, border_color=None, border_width=Pt(1)):
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
             anchor=MSO_ANCHOR.TOP):
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    tf.auto_size = None
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    return txBox


def add_bullet_box(slide, left, top, width, height, title, bullets,
                   title_size=16, bullet_size=13, title_color=DARK_BLUE,
                   bullet_color=DARK_GRAY, bg_color=None, border_color=None,
                   icon=None):
    """Add a box with title and bullet points."""
    if bg_color:
        box = add_shape(slide, left, top, width, height, bg_color,
                        border_color=border_color)

    txBox = slide.shapes.add_textbox(
        left + Inches(0.2), top + Inches(0.15),
        width - Inches(0.4), height - Inches(0.3))
    tf = txBox.text_frame
    tf.word_wrap = True

    # Title
    p = tf.paragraphs[0]
    if icon:
        p.text = f"{icon}  {title}"
    else:
        p.text = title
    p.font.size = Pt(title_size)
    p.font.bold = True
    p.font.color.rgb = title_color
    p.font.name = "Calibri"
    p.space_after = Pt(6)

    # Bullets
    for bullet in bullets:
        p = tf.add_paragraph()
        p.text = bullet
        p.font.size = Pt(bullet_size)
        p.font.color.rgb = bullet_color
        p.font.name = "Calibri"
        p.space_before = Pt(2)
        p.space_after = Pt(2)
        p.level = 0
        # Add bullet character
        pPr = p._pPr
        if pPr is None:
            from pptx.oxml.ns import qn
            pPr = p._p.get_or_add_pPr()

    return txBox


def add_check_item(tf, text, font_size=13, color=DARK_GRAY, check_color=GREEN):
    """Add a checkmark bullet item to a text frame."""
    p = tf.add_paragraph()
    run1 = p.add_run()
    run1.text = "\u2713  "
    run1.font.size = Pt(font_size)
    run1.font.color.rgb = check_color
    run1.font.bold = True
    run1.font.name = "Calibri"
    run2 = p.add_run()
    run2.text = text
    run2.font.size = Pt(font_size)
    run2.font.color.rgb = color
    run2.font.name = "Calibri"
    p.space_before = Pt(3)
    p.space_after = Pt(3)
    return p


def add_header_bar(slide, text, subtitle=None):
    """Add a dark blue header bar across the top."""
    add_shape(slide, Inches(0), Inches(0), prs.slide_width, Inches(1.15),
              DARK_BLUE, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, Inches(0.6), Inches(0.15), Inches(12), Inches(0.55),
             text, font_size=28, bold=True, color=WHITE, font_name="Calibri")
    if subtitle:
        add_text(slide, Inches(0.6), Inches(0.7), Inches(12), Inches(0.4),
                 subtitle, font_size=15, color=LIGHT_BLUE, font_name="Calibri")


def add_footer(slide):
    """Add a thin footer bar."""
    add_shape(slide, Inches(0), Inches(7.1), prs.slide_width, Inches(0.4),
              LIGHT_GRAY, shape_type=MSO_SHAPE.RECTANGLE)
    add_text(slide, Inches(0.5), Inches(7.13), Inches(5), Inches(0.3),
             "DRAIGON | Horizon Europe Grant No. 101137383",
             font_size=9, color=MED_GRAY, font_name="Calibri")
    add_text(slide, Inches(8), Inches(7.13), Inches(5), Inches(0.3),
             "pLIN v2.1 | University Medical Center Groningen",
             font_size=9, color=MED_GRAY, font_name="Calibri",
             alignment=PP_ALIGN.RIGHT)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 1: Title Slide
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
add_bg(slide, DARK_BLUE)

# EU funding stripe at top
add_shape(slide, Inches(0), Inches(0), prs.slide_width, Inches(0.5),
          EU_BLUE, shape_type=MSO_SHAPE.RECTANGLE)
add_text(slide, Inches(0.6), Inches(0.08), Inches(12), Inches(0.35),
         "Funded by the European Union | Horizon Europe | Grant Agreement No. 101137383",
         font_size=11, color=WHITE, font_name="Calibri")

# Main title
add_text(slide, Inches(1), Inches(1.6), Inches(11.3), Inches(1.2),
         "pLIN: Fulfilling DRAIGON Objectives",
         font_size=40, bold=True, color=WHITE, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)

# Subtitle
add_text(slide, Inches(1), Inches(2.9), Inches(11.3), Inches(0.8),
         "Hierarchical Plasmid Classification and AMR Surveillance System",
         font_size=22, color=LIGHT_BLUE, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)

# DRAIGON full name
add_shape(slide, Inches(1.5), Inches(4.0), Inches(10.3), Inches(1.3),
          RGBColor(0x0A, 0x3A, 0x8A), border_color=MED_BLUE)
add_text(slide, Inches(1.8), Inches(4.15), Inches(9.7), Inches(1.0),
         "DRAIGON: Diagnosing Infections with Multi-Drug-Resistant Microorganisms\n"
         "using AI-Powered Genomic Antibiotic Susceptibility Prediction\n"
         "from Long-Read Sequencing Data",
         font_size=16, color=LIGHT_BLUE, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)

# Authors
add_text(slide, Inches(1), Inches(5.6), Inches(11.3), Inches(0.5),
         "Basil Britto Xavier, Anurag Kumar Bari, Bhanu Sinha, John W.A. Rossen",
         font_size=15, color=WHITE, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(1), Inches(6.05), Inches(11.3), Inches(0.4),
         "on behalf of the DRAIGON Consortium",
         font_size=13, color=LIGHT_BLUE, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)
add_text(slide, Inches(1), Inches(6.4), Inches(11.3), Inches(0.4),
         "Department of Medical Microbiology and Infection Prevention | UMCG, University of Groningen",
         font_size=12, color=MED_GRAY, font_name="Calibri",
         alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: DRAIGON Objectives Overview
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "DRAIGON Objectives", "What the consortium aims to achieve")
add_footer(slide)

# Left column: DRAIGON core objectives
obj_data = [
    ("O1: Rapid MDR Diagnostics",
     "End-to-end workflow combining long-read\nWGS with AI for rapid pathogen ID\nand antibiotic susceptibility prediction",
     DARK_BLUE),
    ("O2: AMR Surveillance",
     "Comprehensive antibiogram data and\nresistance gene detection from genomic\ndata in a single assay",
     TEAL),
    ("O3: Outbreak Detection",
     "Early detection system to prevent\ncross-border pathogen spread using\ngenomic cluster information",
     ORANGE),
    ("O4: Antibiotic Stewardship",
     "Right antibiotic, right dose, right\npatient, right time \u2014 guided by\ngenomic resistance profiles",
     GREEN),
    ("O5: Clinical Validation",
     "Validated across 5 hospitals in\nNetherlands, Austria, Albania, USA\nfor BSI and PJI",
     PURPLE),
]

box_w = Inches(3.7)
box_h = Inches(1.0)
start_x = Inches(0.5)
gap = Inches(0.15)

for i, (title, desc, color) in enumerate(obj_data):
    y = Inches(1.45) + i * (box_h + gap)

    # Colour accent bar
    add_shape(slide, start_x, y, Inches(0.08), box_h, color,
              shape_type=MSO_SHAPE.RECTANGLE)

    # Box
    add_shape(slide, start_x + Inches(0.08), y, box_w, box_h,
              LIGHT_GRAY, border_color=RGBColor(0xDD, 0xDD, 0xDD))

    # Title
    add_text(slide, start_x + Inches(0.25), y + Inches(0.05),
             box_w - Inches(0.3), Inches(0.3),
             title, font_size=13, bold=True, color=color)
    # Description
    add_text(slide, start_x + Inches(0.25), y + Inches(0.32),
             box_w - Inches(0.3), Inches(0.65),
             desc, font_size=10, color=MED_GRAY)

# Right column: key numbers
add_text(slide, Inches(4.8), Inches(1.45), Inches(4), Inches(0.4),
         "DRAIGON Consortium", font_size=18, bold=True, color=DARK_BLUE)

# Partners box
add_shape(slide, Inches(4.8), Inches(1.95), Inches(8), Inches(4.7),
          ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(0.5))

partners = [
    ("EVI (Coordinator)", "Heidelberg, Germany"),
    ("UMCG", "Groningen, Netherlands"),
    ("Isala Hospital", "Zwolle, Netherlands"),
    ("OSS", "Vienna, Austria"),
    ("Mayo Clinic", "Rochester, USA"),
    ("Johns Hopkins University", "Baltimore, USA"),
    ("UHSN", "Tirana, Albania"),
    ("Health-Ecore B.V.", "Netherlands"),
    ("Sandoz Pharma (Assoc.)", "Switzerland"),
]

add_text(slide, Inches(5.1), Inches(2.05), Inches(3.5), Inches(0.35),
         "9 Partners across 6 Countries", font_size=14, bold=True, color=DARK_BLUE)

for i, (name, loc) in enumerate(partners):
    y_p = Inches(2.45) + i * Inches(0.42)
    add_text(slide, Inches(5.3), y_p, Inches(3.0), Inches(0.22),
             name, font_size=11, bold=True, color=DARK_GRAY)
    add_text(slide, Inches(8.3), y_p, Inches(3.5), Inches(0.22),
             loc, font_size=11, color=MED_GRAY)

# Funding
add_text(slide, Inches(5.1), Inches(6.4), Inches(7), Inches(0.35),
         "Total EU Contribution: \u20ac6.35 million  |  Duration: 2024\u20132027",
         font_size=13, bold=True, color=DARK_BLUE)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: pLIN \u2192 DRAIGON Mapping (core slide)
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "How pLIN Fulfils DRAIGON Objectives",
               "Mapping pLIN capabilities to consortium goals")
add_footer(slide)

mappings = [
    {
        "objective": "O1: Rapid MDR Diagnostics",
        "obj_color": DARK_BLUE,
        "items": [
            "KNN classifier (91.1% accuracy) identifies Inc/Rep group in seconds",
            "Processes FASTA input through full pipeline in <30 min on standard hardware",
            "28 Inc/Rep groups: Gram-neg, Gram-pos, A. baumannii, P. aeruginosa",
            "Auto contig identification separates plasmid from chromosome",
        ]
    },
    {
        "objective": "O2: AMR Surveillance",
        "obj_color": TEAL,
        "items": [
            "AMRFinderPlus integration: 64,891 gene hits across 83.1% of plasmids",
            "1,635 carbapenemase + 1,804 ESBL + 204 mcr + 2,315 PMQR detections",
            "Lineage-specific resistance profiles link AMR genes to plasmid backbone",
            "79,305-plasmid reference database with 57,886 unique pLIN codes",
        ]
    },
    {
        "objective": "O3: Outbreak Detection",
        "obj_color": ORANGE,
        "items": [
            "6-level hierarchical codes: L6 (99.9% ANI) detects outbreak clusters",
            "Validated against 74 plasmids from 27 published outbreaks, 13 countries",
            "85.1% high-confidence classifications; 9 intra-study clusters detected",
            "Combined chromosomal + plasmid typing distinguishes clonal vs HGT spread",
        ]
    },
    {
        "objective": "O4: Antibiotic Stewardship",
        "obj_color": GREEN,
        "items": [
            "3-tier clinical risk stratification (Critical / High / Moderate)",
            "High-risk lineage alerts: pLIN 671 (100% KPC-2), pLIN 860 (44% mcr)",
            "Simpson's D = 0.985 \u2014 1.54\u00d7 improvement over Inc typing (D = 0.641)",
            "Actionable reports linking plasmid lineage to resistance profile",
        ]
    },
    {
        "objective": "O5: Clinical Deployability",
        "obj_color": PURPLE,
        "items": [
            "Interactive Streamlit GUI \u2014 no bioinformatics expertise required",
            "Cross-platform: macOS, Linux, Windows one-click installers",
            "Open-source (GPL-3.0) at GitHub; runs on standard laptop/workstation",
            "7 analytical modules: completeness, novelty, recombination, MGE, stability",
        ]
    },
]

col_w = Inches(6.0)
row_h = Inches(1.05)
left_x = Inches(0.4)
right_x = Inches(6.8)

for i, m in enumerate(mappings):
    y = Inches(1.35) + i * (row_h + Inches(0.08))

    # Objective label (left)
    add_shape(slide, left_x, y, Inches(0.06), row_h, m["obj_color"],
              shape_type=MSO_SHAPE.RECTANGLE)
    add_shape(slide, left_x + Inches(0.06), y, Inches(2.4), row_h,
              LIGHT_GRAY, border_color=RGBColor(0xDD, 0xDD, 0xDD))
    txBox = slide.shapes.add_textbox(
        left_x + Inches(0.2), y + Inches(0.05),
        Inches(2.2), row_h - Inches(0.1))
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = m["objective"]
    p.font.size = Pt(13)
    p.font.bold = True
    p.font.color.rgb = m["obj_color"]
    p.font.name = "Calibri"

    # Arrow
    add_shape(slide, Inches(2.95), y + Inches(0.3), Inches(0.5), Inches(0.35),
              m["obj_color"], shape_type=MSO_SHAPE.RIGHT_ARROW)

    # pLIN deliverables (right)
    add_shape(slide, Inches(3.55), y, Inches(9.3), row_h,
              ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(0.5))

    txBox2 = slide.shapes.add_textbox(
        Inches(3.7), y + Inches(0.02),
        Inches(9.0), row_h - Inches(0.04))
    tf2 = txBox2.text_frame
    tf2.word_wrap = True

    for j, item in enumerate(m["items"]):
        if j == 0:
            p2 = tf2.paragraphs[0]
        else:
            p2 = tf2.add_paragraph()
        run_check = p2.add_run()
        run_check.text = "\u2713  "
        run_check.font.size = Pt(11)
        run_check.font.color.rgb = GREEN
        run_check.font.bold = True
        run_check.font.name = "Calibri"
        run_txt = p2.add_run()
        run_txt.text = item
        run_txt.font.size = Pt(11)
        run_txt.font.color.rgb = DARK_GRAY
        run_txt.font.name = "Calibri"
        p2.space_before = Pt(1)
        p2.space_after = Pt(1)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: Key Results at a Glance
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Key Results at a Glance",
               "pLIN performance metrics aligned with DRAIGON deliverables")
add_footer(slide)

metrics = [
    ("79,305", "Plasmids\nClassified", DARK_BLUE, "Reference database\nwith 57,886 unique codes"),
    ("28", "Inc/Rep\nGroups", TEAL, "Gram-neg, Gram-pos,\nA. baumannii, P. aeruginosa"),
    ("91.1%", "Classifier\nAccuracy", GREEN, "KNN k=5, cosine metric\n5-fold cross-validated"),
    ("0.985", "Simpson's\nDiversity (D)", ORANGE, "1.54\u00d7 improvement\nover Inc typing alone"),
    ("64,891", "AMR Gene\nDetections", RED, "Including 1,635 carba-\npenemase genes"),
    ("27", "Outbreak\nStudies Validated", PURPLE, "74 plasmids, 13 countries\n85.1% high-confidence"),
]

box_w = Inches(1.85)
box_h = Inches(2.6)
start_x = Inches(0.55)
gap_x = Inches(0.25)

for i, (number, label, color, detail) in enumerate(metrics):
    x = start_x + i * (box_w + gap_x)
    y = Inches(1.5)

    # Card
    add_shape(slide, x, y, box_w, box_h, LIGHT_GRAY,
              border_color=color, border_width=Pt(1.5))
    # Color top accent
    add_shape(slide, x, y, box_w, Inches(0.06), color,
              shape_type=MSO_SHAPE.RECTANGLE)

    # Number
    add_text(slide, x, y + Inches(0.2), box_w, Inches(0.6),
             number, font_size=32, bold=True, color=color,
             alignment=PP_ALIGN.CENTER)
    # Label
    add_text(slide, x, y + Inches(0.8), box_w, Inches(0.6),
             label, font_size=13, bold=True, color=DARK_GRAY,
             alignment=PP_ALIGN.CENTER)
    # Detail
    add_text(slide, x + Inches(0.1), y + Inches(1.5), box_w - Inches(0.2), Inches(0.9),
             detail, font_size=10, color=MED_GRAY,
             alignment=PP_ALIGN.CENTER)

# Bottom summary bar
add_shape(slide, Inches(0.5), Inches(4.5), Inches(12.3), Inches(0.7),
          ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(0.5))
add_text(slide, Inches(0.8), Inches(4.55), Inches(11.7), Inches(0.6),
         "pLIN provides the first permanent, hierarchical plasmid classification system integrating "
         "multi-resolution typing with AMR surveillance and automated outbreak detection \u2014 "
         "directly enabling DRAIGON\u2019s mission for AI-powered genomic AMR diagnostics.",
         font_size=13, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

# Analytical modules bar
add_text(slide, Inches(0.5), Inches(5.45), Inches(12.3), Inches(0.4),
         "7 Analytical Modules Addressing Key Surveillance Gaps",
         font_size=16, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

modules = [
    ("Assembly\nCompleteness", "94.2% correct\ncategorisation"),
    ("Database\nCoverage", "Traffic-light\nnovelty alerts"),
    ("Recombination\nDetection", "87.3% sensitivity\n94.1% specificity"),
    ("Novel Group\nDiscovery", "18 putative groups\nfrom 2,109 sequences"),
    ("Evolutionary\nRate", "1.8\u20138.4 \u00d7 10\u207b\u2076\nsubs/site/year"),
    ("Cluster\nStability", "Mean ARI > 0.85\n20/28 groups"),
    ("MGE Boundary\nDetection", "94.7% IS element\nsensitivity"),
]

mod_w = Inches(1.65)
mod_h = Inches(1.1)
mod_start_x = Inches(0.55)
mod_gap = Inches(0.15)

for i, (name, perf) in enumerate(modules):
    x = mod_start_x + i * (mod_w + mod_gap)
    y = Inches(5.85)
    add_shape(slide, x, y, mod_w, mod_h, DARK_BLUE,
              border_color=MED_BLUE, border_width=Pt(0.5))
    add_text(slide, x, y + Inches(0.05), mod_w, Inches(0.45),
             name, font_size=10, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.55), mod_w, Inches(0.45),
             perf, font_size=9, color=LIGHT_BLUE,
             alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: High-Risk Lineage Discovery (clinical impact)
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "High-Risk Lineage Discovery",
               "Clinically actionable intelligence invisible to conventional typing")
add_footer(slide)

# Two example lineages
lineages = [
    {
        "name": "pLIN 671 (IncN)",
        "color": RED,
        "stats": [
            ("Members", "n = 90"),
            ("blaKPC-2 carriage", "100%"),
            ("Mean AMR genes", "13.2 per plasmid"),
            ("Risk tier", "CRITICAL"),
        ],
        "note": "A dominant KPC-2-carrying IncN lineage spanning multiple\n"
                "countries \u2014 invisible under conventional 'IncN' typing.\n"
                "All members carry identical carbapenemase backbone.",
    },
    {
        "name": "pLIN 860 (multi-Inc)",
        "color": ORANGE,
        "stats": [
            ("Members", "n = 142 (5 Inc groups)"),
            ("mcr carriage", "44.4%"),
            ("Mean AMR genes", "14.4 per plasmid"),
            ("Risk tier", "CRITICAL"),
        ],
        "note": "Cross-Inc lineage spanning IncN, IncHI2, IncFII, IncHI1, IncX1.\n"
                "Carries mobile colistin resistance alongside carbapenemases.\n"
                "Demonstrates multi-drug resistance stacking on a single backbone.",
    },
]

for idx, lin in enumerate(lineages):
    x_base = Inches(0.5) + idx * Inches(6.4)
    y_base = Inches(1.4)

    # Box
    add_shape(slide, x_base, y_base, Inches(6.0), Inches(3.6),
              LIGHT_GRAY, border_color=lin["color"], border_width=Pt(2))
    # Color accent top
    add_shape(slide, x_base, y_base, Inches(6.0), Inches(0.06),
              lin["color"], shape_type=MSO_SHAPE.RECTANGLE)

    # Name
    add_text(slide, x_base + Inches(0.3), y_base + Inches(0.15),
             Inches(5.4), Inches(0.4),
             lin["name"], font_size=20, bold=True, color=lin["color"])

    # Stats
    for j, (key, val) in enumerate(lin["stats"]):
        y_s = y_base + Inches(0.65) + j * Inches(0.38)
        add_text(slide, x_base + Inches(0.4), y_s,
                 Inches(2.5), Inches(0.3),
                 key + ":", font_size=12, bold=True, color=DARK_GRAY)
        add_text(slide, x_base + Inches(3.0), y_s,
                 Inches(2.8), Inches(0.3),
                 val, font_size=12, bold=True, color=lin["color"])

    # Note
    add_text(slide, x_base + Inches(0.3), y_base + Inches(2.35),
             Inches(5.4), Inches(1.1),
             lin["note"], font_size=11, color=MED_GRAY)

# Bottom: DRAIGON relevance
add_shape(slide, Inches(0.5), Inches(5.3), Inches(12.3), Inches(1.5),
          ACCENT_BLUE, border_color=MED_BLUE, border_width=Pt(0.5))

add_text(slide, Inches(0.8), Inches(5.4), Inches(11.7), Inches(0.35),
         "DRAIGON Relevance: From Genomic Data to Clinical Action",
         font_size=16, bold=True, color=DARK_BLUE)

relevance_items = [
    "Conventional typing reports 'IncN' or 'IncFII' \u2014 pLIN resolves these into 3,073 distinct lineages with specific resistance profiles",
    "Enables clinicians to distinguish single-plasmid outbreaks from independent resistance acquisition events",
    "Combined chromosomal (MLST) + plasmid (pLIN) typing classifies transmission as clonal spread vs. horizontal transfer",
    "Directly supports DRAIGON\u2019s goal: 'right antibiotic, right dose, right patient, right time' through lineage-linked AMR profiles",
]

txBox = slide.shapes.add_textbox(
    Inches(0.8), Inches(5.8), Inches(11.7), Inches(0.95))
tf = txBox.text_frame
tf.word_wrap = True

for k, item in enumerate(relevance_items):
    if k == 0:
        p = tf.paragraphs[0]
    else:
        p = tf.add_paragraph()
    run_a = p.add_run()
    run_a.text = "\u25b8  "
    run_a.font.size = Pt(11)
    run_a.font.color.rgb = DARK_BLUE
    run_a.font.name = "Calibri"
    run_b = p.add_run()
    run_b.text = item
    run_b.font.size = Pt(11)
    run_b.font.color.rgb = DARK_GRAY
    run_b.font.name = "Calibri"
    p.space_after = Pt(2)


# ═══════════════════════════════════════════════════════════════════════════════
# SLIDE 6: Summary / Conclusions
# ═══════════════════════════════════════════════════════════════════════════════
slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, DARK_BLUE)

# Title
add_text(slide, Inches(0.6), Inches(0.4), Inches(12), Inches(0.6),
         "Summary: pLIN Delivers on DRAIGON Objectives",
         font_size=30, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

# Separator line
add_shape(slide, Inches(4), Inches(1.1), Inches(5.3), Inches(0.03),
          MED_BLUE, shape_type=MSO_SHAPE.RECTANGLE)

conclusions = [
    {
        "icon": "\u2713",
        "title": "AI-Powered Genomic Diagnostics",
        "text": "KNN classifier (91.1%) + 7 analytical modules provide comprehensive, automated plasmid characterisation from raw sequence data",
    },
    {
        "icon": "\u2713",
        "title": "Integrated AMR Surveillance",
        "text": "64,891 resistance gene detections mapped to hierarchical lineage codes \u2014 tracking resistance dissemination through plasmid backbone, not species",
    },
    {
        "icon": "\u2713",
        "title": "Cross-Border Outbreak Detection",
        "text": "Validated across 27 studies, 13 countries, 7 resistance mechanisms \u2014 permanent codes enable inter-institutional comparison",
    },
    {
        "icon": "\u2713",
        "title": "Clinical Deployability",
        "text": "Open-source Streamlit GUI runs on standard hardware in <30 min, no bioinformatics expertise needed \u2014 ready for DRAIGON clinical sites",
    },
    {
        "icon": "\u2713",
        "title": "Scalable & Future-Proof",
        "text": "79,305 plasmids classified; permanent codes guarantee longitudinal surveillance compatibility as DRAIGON expands across partner sites",
    },
]

for i, c in enumerate(conclusions):
    y = Inches(1.5) + i * Inches(1.05)
    # Check icon
    add_shape(slide, Inches(0.8), y, Inches(0.5), Inches(0.5),
              GREEN, shape_type=MSO_SHAPE.OVAL)
    add_text(slide, Inches(0.8), y + Inches(0.05), Inches(0.5), Inches(0.45),
             c["icon"], font_size=20, bold=True, color=WHITE,
             alignment=PP_ALIGN.CENTER)
    # Title
    add_text(slide, Inches(1.6), y + Inches(0.0), Inches(10.5), Inches(0.35),
             c["title"], font_size=16, bold=True, color=WHITE)
    # Text
    add_text(slide, Inches(1.6), y + Inches(0.4), Inches(10.5), Inches(0.55),
             c["text"], font_size=13, color=LIGHT_BLUE)

# Funding acknowledgement
add_shape(slide, Inches(0.5), Inches(6.75), Inches(12.3), Inches(0.55),
          RGBColor(0x0A, 0x3A, 0x8A), shape_type=MSO_SHAPE.RECTANGLE)
add_text(slide, Inches(0.8), Inches(6.8), Inches(11.7), Inches(0.45),
         "Funded by the European Union\u2019s Horizon Europe programme | Grant No. 101137383 | "
         "Views expressed are those of the authors and do not necessarily reflect those of the EU or HADEA.",
         font_size=10, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════
out_path = os.path.join(OUT_DIR, "pLIN_DRAIGON_Objectives.pptx")
prs.save(out_path)
print(f"Slide deck saved to: {out_path}")
print(f"Total slides: {len(prs.slides)}")


if __name__ == "__main__":
    pass  # All code runs at module level
