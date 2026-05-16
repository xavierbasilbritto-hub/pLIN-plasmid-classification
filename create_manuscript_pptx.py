#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""Generate PowerPoint presentation for the pLIN manuscript."""

from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
import os

# Colors
DARK_BLUE = RGBColor(0x0D, 0x47, 0xA1)
MED_BLUE = RGBColor(0x1E, 0x88, 0xE5)
LIGHT_BLUE = RGBColor(0xBB, 0xDE, 0xFB)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
BLACK = RGBColor(0x00, 0x00, 0x00)
DARK_GRAY = RGBColor(0x33, 0x33, 0x33)
MED_GRAY = RGBColor(0x75, 0x75, 0x75)
LIGHT_GRAY = RGBColor(0xEC, 0xEF, 0xF1)
RED = RGBColor(0xE5, 0x39, 0x35)
GREEN = RGBColor(0x43, 0xA0, 0x47)
ORANGE = RGBColor(0xFB, 0x8C, 0x00)
PURPLE = RGBColor(0x8E, 0x24, 0xAA)
TEAL = RGBColor(0x00, 0x96, 0x88)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(BASE_DIR, "output", "figures")

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)


def add_bg(slide, color=WHITE):
    """Set slide background."""
    bg = slide.background
    fill = bg.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_box(slide, left, top, width, height, fill_color, border_color=None, border_width=Pt(1)):
    """Add a rounded rectangle."""
    shape = slide.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if border_color:
        shape.line.color.rgb = border_color
        shape.line.width = border_width
    else:
        shape.line.fill.background()
    return shape


def add_text(slide, left, top, width, height, text, font_size=14, bold=False,
             color=DARK_GRAY, alignment=PP_ALIGN.LEFT, font_name="Calibri"):
    """Add a text box."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = alignment
    return txBox


def add_multiline(slide, left, top, width, height, lines, font_size=12,
                  color=DARK_GRAY, bold_first=False, spacing=1.0, font_name="Calibri"):
    """Add multi-line text box."""
    txBox = slide.shapes.add_textbox(left, top, width, height)
    tf = txBox.text_frame
    tf.word_wrap = True
    for i, line in enumerate(lines):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = line
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.font.name = font_name
        p.space_after = Pt(spacing * 2)
        if bold_first and i == 0:
            p.font.bold = True
    return txBox


def add_arrow(slide, left, top, width, height, color=MED_BLUE):
    """Add a down arrow."""
    shape = slide.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_right_arrow(slide, left, top, width, height, color=MED_BLUE):
    """Add a right arrow."""
    shape = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = color
    shape.line.fill.background()
    return shape


def add_header_bar(slide, text, subtitle=""):
    """Add a dark blue header bar at top."""
    add_box(slide, Inches(0), Inches(0), Inches(13.333), Inches(1.2), DARK_BLUE)
    add_text(slide, Inches(0.6), Inches(0.15), Inches(12), Inches(0.6),
             text, font_size=32, bold=True, color=WHITE)
    if subtitle:
        add_text(slide, Inches(0.6), Inches(0.7), Inches(12), Inches(0.4),
                 subtitle, font_size=16, color=LIGHT_BLUE)


def add_figure(slide, filename, left, top, width, height=None):
    """Embed a PNG figure if it exists, otherwise add a placeholder box."""
    fig_path = os.path.join(FIG_DIR, filename)
    if os.path.isfile(fig_path):
        if height:
            slide.shapes.add_picture(fig_path, left, top, width, height)
        else:
            slide.shapes.add_picture(fig_path, left, top, width=width)
    else:
        # Placeholder box if figure not found
        h = height if height else Inches(4.0)
        box = add_box(slide, left, top, width, h, LIGHT_GRAY, border_color=MED_GRAY, border_width=Pt(2))
        add_text(slide, left, top + h // 2 - Inches(0.2), width, Inches(0.4),
                 f"[Figure: {filename}]", font_size=14, bold=True,
                 color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 1: Title Slide
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
add_bg(slide, DARK_BLUE)

add_text(slide, Inches(0.8), Inches(1.0), Inches(11.5), Inches(1),
         "pLIN: Plasmid Lineage Identification Number",
         font_size=40, bold=True, color=WHITE)

add_text(slide, Inches(0.8), Inches(2.1), Inches(11.5), Inches(0.8),
         "A hierarchical plasmid classification and surveillance system\n"
         "linking antimicrobial resistance to transmissible lineages",
         font_size=22, color=LIGHT_BLUE)

add_text(slide, Inches(0.8), Inches(3.2), Inches(11.5), Inches(0.5),
         "Basil Britto Xavier, Anurag Kumar Bari, Bhanu Sinha, John W A Rossen",
         font_size=18, color=RGBColor(0x90, 0xCA, 0xF9))

# Hero stat boxes
hero_stats = [
    ("79,305", "Plasmids"),
    ("28", "Inc/Rep Groups"),
    ("57,886", "Unique Codes"),
    ("91.1%", "Accuracy"),
    ("64,891", "AMR Detections"),
]
for i, (num, label) in enumerate(hero_stats):
    x = Inches(0.8 + i * 2.5)
    y = Inches(4.2)
    box = add_box(slide, x, y, Inches(2.2), Inches(1.8), RGBColor(0x0A, 0x2F, 0x6E),
                  border_color=MED_BLUE, border_width=Pt(2))
    add_text(slide, x, y + Inches(0.2), Inches(2.2), Inches(0.7),
             num, font_size=36, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.9), Inches(2.2), Inches(0.7),
             label, font_size=14, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)

add_text(slide, Inches(0.8), Inches(6.5), Inches(11.5), Inches(0.5),
         "Manuscript Presentation", font_size=14, color=MED_GRAY)


# ======================================================================
# SLIDE 2: The Problem — AMR & Plasmid Classification Gap
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "The Problem: No Unified Plasmid Nomenclature",
               "AMR dissemination via plasmid HGT demands a standardized classification system")

# Left side: 3 key stat boxes
stat_items = [
    ("4.95M", "Deaths associated\nwith AMR (2019)", RED),
    ("1.27M", "Deaths directly\ncaused by AMR", RGBColor(0xC6, 0x28, 0x28)),
    ("HGT", "Plasmid horizontal gene\ntransfer = primary AMR\ndissemination mechanism", ORANGE),
]
for i, (num, desc, color) in enumerate(stat_items):
    y = Inches(1.5 + i * 1.6)
    add_box(slide, Inches(0.5), y, Inches(5.5), Inches(1.4), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, Inches(0.7), y + Inches(0.15), Inches(1.5), Inches(0.5),
             num, font_size=28, bold=True, color=color)
    add_text(slide, Inches(2.3), y + Inches(0.2), Inches(3.5), Inches(1.0),
             desc, font_size=13, color=DARK_GRAY)

# Right side: Current tools are fragmented
add_text(slide, Inches(6.5), Inches(1.5), Inches(6.5), Inches(0.5),
         "Current Tools Are Fragmented", font_size=18, bold=True, color=DARK_BLUE)

tools_frag = [
    ("PlasmidFinder", "Flat, single-level replicon typing", MED_GRAY),
    ("pMLST", "Only 6 Inc groups supported", MED_GRAY),
    ("MOB-suite", "Codes change when DB updates", MED_GRAY),
    ("COPLA", "Only 41% of plasmids assignable", MED_GRAY),
    ("mge-cluster", "No permanent codes assigned", MED_GRAY),
]
for i, (tool, issue, color) in enumerate(tools_frag):
    y = Inches(2.1 + i * 0.75)
    add_box(slide, Inches(6.5), y, Inches(6.3), Inches(0.6), WHITE, border_color=color, border_width=Pt(1))
    add_text(slide, Inches(6.7), y + Inches(0.1), Inches(2.0), Inches(0.4),
             tool, font_size=13, bold=True, color=RED, font_name="Consolas")
    add_text(slide, Inches(8.8), y + Inches(0.1), Inches(3.8), Inches(0.4),
             issue, font_size=12, color=DARK_GRAY)

# Bottom callout
add_box(slide, Inches(0.5), Inches(6.2), Inches(12.4), Inches(0.9), RGBColor(0xFF, 0xEB, 0xEE),
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.35), Inches(12), Inches(0.6),
         "No tool combines: hierarchical classification + code permanence + AMR integration + outbreak detection",
         font_size=15, bold=True, color=RED, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 3: The Solution — pLIN Framework
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "The Solution: pLIN Framework",
               "Six-level hierarchical coding with permanent, stable plasmid identifiers")

# Hierarchy levels
levels = [
    ("L1", "d\u22640.150", "~85% ANI", "Family", RGBColor(0xE5, 0x39, 0x35)),
    ("L2", "d\u22640.100", "~90% ANI", "Subfamily", RGBColor(0xFB, 0x8C, 0x00)),
    ("L3", "d\u22640.050", "~95% ANI", "Cluster", RGBColor(0xFD, 0xD8, 0x35)),
    ("L4", "d\u22640.020", "~98% ANI", "Subcluster", RGBColor(0x43, 0xA0, 0x47)),
    ("L5", "d\u22640.010", "~99% ANI", "Clone complex", RGBColor(0x1E, 0x88, 0xE5)),
    ("L6", "d\u22640.001", "~99.9% ANI", "Strain", RGBColor(0x8E, 0x24, 0xAA)),
]

# Level flow
for i, (lvl, thresh, ani, rank, color) in enumerate(levels):
    x = Inches(0.5 + i * 2.05)
    y = Inches(1.6)
    add_box(slide, x, y, Inches(1.85), Inches(2.0), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, x, y + Inches(0.1), Inches(1.85), Inches(0.4),
             lvl, font_size=22, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.55), Inches(1.85), Inches(0.3),
             thresh, font_size=13, color=DARK_GRAY, alignment=PP_ALIGN.CENTER, font_name="Consolas")
    add_text(slide, x, y + Inches(0.85), Inches(1.85), Inches(0.3),
             ani, font_size=12, color=MED_GRAY, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(1.2), Inches(1.85), Inches(0.3),
             rank, font_size=12, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    if i < len(levels) - 1:
        add_right_arrow(slide, x + Inches(1.9), y + Inches(0.8), Inches(0.15), Inches(0.3), color)

# Example code
add_box(slide, Inches(0.5), Inches(4.0), Inches(12.4), Inches(0.9), RGBColor(0xE3, 0xF2, 0xFD),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(4.1), Inches(3), Inches(0.3),
         "Example pLIN Code:", font_size=14, bold=True, color=DARK_BLUE)
add_text(slide, Inches(4.0), Inches(4.1), Inches(8.5), Inches(0.7),
         "1 . 1 . 2 . 15 . 48 . 671", font_size=28, bold=True, color=DARK_BLUE,
         alignment=PP_ALIGN.CENTER, font_name="Consolas")

# 4 key innovations
innovations = [
    ("Permanent Codes", "Codes never change\nwhen new sequences\nare added", MED_BLUE),
    ("6-Level Hierarchy", "Family to Strain\nresolution in a\nsingle code", GREEN),
    ("AMR Integration", "64,891 AMR detections\nmapped to specific\nplasmid lineages", RED),
    ("Outbreak Detection", "Automated risk-tier\nassignment from\nshared pLIN + AMR", PURPLE),
]
for i, (title, desc, color) in enumerate(innovations):
    x = Inches(0.5 + i * 3.15)
    y = Inches(5.2)
    add_box(slide, x, y, Inches(2.9), Inches(0.5), color)
    add_text(slide, x, y + Inches(0.07), Inches(2.9), Inches(0.35),
             title, font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, y + Inches(0.5), Inches(2.9), Inches(1.5), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.1), y + Inches(0.6), Inches(2.7), Inches(1.3),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 4: Methodology — Pipeline Overview
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "pLIN Classification Pipeline",
               "From FASTA input to hierarchical plasmid codes")

add_figure(slide, "Figure8_pipeline_overview.png",
           Inches(0.5), Inches(1.5), Inches(12.3), Inches(5.0))

add_text(slide, Inches(0.5), Inches(6.7), Inches(12.3), Inches(0.5),
         "4-mer frequency vectors (256 features) \u2192 cosine distance \u2192 single-linkage hierarchical clustering \u2192 6-level pLIN code",
         font_size=13, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 5: Dataset Overview (Figure 1)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Training Dataset: 8,077 Plasmids Across 28 Inc/Rep Groups",
               "Curated from NCBI RefSeq with replicon-typed training sequences")

add_figure(slide, "Figure1_dataset_overview.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(5.5))

# Key stats box (right side)
add_box(slide, Inches(9.6), Inches(1.5), Inches(3.4), Inches(5.5), LIGHT_GRAY,
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(9.8), Inches(1.6), Inches(3.0), Inches(0.4),
         "Key Statistics", font_size=16, bold=True, color=DARK_BLUE)

stats_lines = [
    "IncFII: 66.1%",
    "(4,629 plasmids)",
    "",
    "IncN: 15.7%",
    "(1,064 plasmids)",
    "",
    "IncX1: 10.1%",
    "(701 plasmids)",
    "",
    "17 additional groups:",
    "8.1% combined",
    "",
    "Total: 8,077 plasmids",
    "28 Inc/Rep groups",
]
add_multiline(slide, Inches(9.8), Inches(2.1), Inches(3.0), Inches(4.5),
              stats_lines, font_size=12, color=DARK_GRAY)


# ======================================================================
# SLIDE 6: Classification Performance (Figure 12)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Classification Performance",
               "KNN classifier validated with 5-fold stratified cross-validation")

add_figure(slide, "Figure12_inc_training_performance.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(5.5))

# Stats panel (right)
add_box(slide, Inches(9.6), Inches(1.5), Inches(3.4), Inches(5.5), LIGHT_GRAY,
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(9.8), Inches(1.6), Inches(3.0), Inches(0.4),
         "Performance", font_size=16, bold=True, color=GREEN)

perf_items = [
    ("KNN Accuracy", "91.1%", MED_BLUE),
    ("XGBoost F1", "0.896", GREEN),
    ("Simpson's D", "0.985", PURPLE),
    ("vs Inc typing", "0.641", RED),
    ("Improvement", "1.54\u00d7", ORANGE),
]
for i, (label, value, color) in enumerate(perf_items):
    y = Inches(2.2 + i * 0.9)
    add_text(slide, Inches(9.8), y, Inches(3.0), Inches(0.3),
             label, font_size=12, color=MED_GRAY)
    add_text(slide, Inches(9.8), y + Inches(0.3), Inches(3.0), Inches(0.4),
             value, font_size=24, bold=True, color=color, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 7: Hierarchical Resolution (Figure 6)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Hierarchical Resolution \u2014 3,073 Unique Strain-Level Codes",
               "Six-level hierarchy resolves plasmid diversity at multiple scales")

add_figure(slide, "Figure6_pLIN_hierarchy.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(5.5))

# Stats panel (right)
add_box(slide, Inches(9.6), Inches(1.5), Inches(3.4), Inches(2.5), RGBColor(0xE3, 0xF2, 0xFD),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(9.8), Inches(1.6), Inches(3.0), Inches(0.4),
         "Resolution", font_size=16, bold=True, color=DARK_BLUE)
add_multiline(slide, Inches(9.8), Inches(2.1), Inches(3.0), Inches(1.8), [
    "97.8% Inc concordance",
    "76.0% singletons",
    "Simpson's D = 0.985",
    "3,073 unique L6 codes",
], font_size=13, color=DARK_GRAY)

add_box(slide, Inches(9.6), Inches(4.3), Inches(3.4), Inches(2.7), RGBColor(0xF3, 0xE5, 0xF5),
        border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(9.8), Inches(4.4), Inches(3.0), Inches(0.4),
         "Hierarchy Levels", font_size=14, bold=True, color=PURPLE)
add_multiline(slide, Inches(9.8), Inches(4.9), Inches(3.0), Inches(2.0), [
    "L1: 8 families",
    "L2: 16 subfamilies",
    "L3: 112 clusters",
    "L4: 571 subclusters",
    "L5: 1,315 clone groups",
    "L6: 3,073 strains",
], font_size=12, color=DARK_GRAY)


# ======================================================================
# SLIDE 8: Method Comparison (Figure 10)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "pLIN vs Existing Tools",
               "Comprehensive comparison across 7 key criteria")

add_figure(slide, "Figure10_method_comparison.png",
           Inches(0.3), Inches(1.5), Inches(8.5), Inches(5.5))

# Comparison summary (right side)
add_box(slide, Inches(9.1), Inches(1.5), Inches(3.9), Inches(5.5), LIGHT_GRAY,
        border_color=DARK_BLUE, border_width=Pt(2))
add_text(slide, Inches(9.3), Inches(1.6), Inches(3.5), Inches(0.4),
         "7 Key Criteria", font_size=16, bold=True, color=DARK_BLUE)

criteria = [
    "Hierarchical classification",
    "Permanent code stability",
    "Reference-free approach",
    "AMR gene integration",
    "Outbreak detection",
    "Multi-level resolution",
    "Cross-Inc lineage tracking",
]
for i, crit in enumerate(criteria):
    y = Inches(2.2 + i * 0.65)
    # Green check circle
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(9.3), y + Inches(0.05), Inches(0.25), Inches(0.25))
    circ.fill.solid()
    circ.fill.fore_color.rgb = GREEN
    circ.line.fill.background()
    add_text(slide, Inches(9.7), y, Inches(3.1), Inches(0.35),
             crit, font_size=12, color=DARK_GRAY)

add_box(slide, Inches(9.1), Inches(6.0), Inches(3.9), Inches(0.9), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(9.3), Inches(6.15), Inches(3.5), Inches(0.6),
         "pLIN is the ONLY tool\nmeeting all 7 criteria",
         font_size=14, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 9: AMR Landscape Overview (Figure 2)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "AMR Gene Surveillance: 64,891 Detections",
               "Comprehensive antimicrobial resistance profiling across all classified plasmids")

add_figure(slide, "Figure2_AMR_prevalence.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(5.5))

# Stats panel
add_box(slide, Inches(9.6), Inches(1.5), Inches(3.4), Inches(5.5), LIGHT_GRAY,
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(9.8), Inches(1.6), Inches(3.0), Inches(0.4),
         "Detection Summary", font_size=15, bold=True, color=RED)

amr_stats = [
    ("29,583", "AMR genes", RED),
    ("6,286", "Virulence factors", PURPLE),
    ("29,022", "Stress resistance", ORANGE),
    ("83.1%", "Plasmids with hits", MED_BLUE),
]
for i, (num, label, color) in enumerate(amr_stats):
    y = Inches(2.2 + i * 1.1)
    add_text(slide, Inches(9.8), y, Inches(3.0), Inches(0.4),
             num, font_size=24, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, Inches(9.8), y + Inches(0.4), Inches(3.0), Inches(0.3),
             label, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 10: Critical Resistance Determinants (Figure 3)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Clinically Critical Resistance: Carbapenemases, ESBLs, mcr, PMQR",
               "High-priority resistance determinants mapped across plasmid lineages")

add_figure(slide, "Figure3_critical_AMR.png",
           Inches(0.3), Inches(1.5), Inches(8.5), Inches(4.5))

# 4 stat boxes at bottom
crit_stats = [
    ("1,635", "Carbapenemases", RED),
    ("1,804", "ESBLs", ORANGE),
    ("204", "mcr (colistin)", PURPLE),
    ("2,315", "PMQR", MED_BLUE),
]
for i, (num, label, color) in enumerate(crit_stats):
    x = Inches(0.5 + i * 3.15)
    y = Inches(6.2)
    add_box(slide, x, y, Inches(2.9), Inches(1.0), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, x, y + Inches(0.05), Inches(2.9), Inches(0.4),
             num, font_size=24, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.5), Inches(2.9), Inches(0.3),
             label, font_size=13, bold=True, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 11: High-Risk Lineages (Figure 4)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "High-Risk Plasmid Lineages Revealed by pLIN",
               "AMR gene cargo mapped to specific pLIN codes")

add_figure(slide, "Figure4_pLIN_AMR_heatmap.png",
           Inches(0.3), Inches(1.5), Inches(8.5), Inches(4.3))

# Two highlighted lineages
add_box(slide, Inches(0.5), Inches(6.0), Inches(6.0), Inches(1.2), RGBColor(0xFF, 0xEB, 0xEE),
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.1), Inches(5.6), Inches(0.3),
         "pLIN 671 (IncN, n=90)", font_size=14, bold=True, color=RED)
add_text(slide, Inches(0.7), Inches(6.45), Inches(5.6), Inches(0.6),
         "100% blaKPC-2 carriage | 13.2 mean AMR genes per plasmid",
         font_size=12, color=DARK_GRAY)

add_box(slide, Inches(6.8), Inches(6.0), Inches(6.0), Inches(1.2), RGBColor(0xF3, 0xE5, 0xF5),
        border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(7.0), Inches(6.1), Inches(5.6), Inches(0.3),
         "pLIN 860 (5 Inc groups, n=142)", font_size=14, bold=True, color=PURPLE)
add_text(slide, Inches(7.0), Inches(6.45), Inches(5.6), Inches(0.6),
         "44.4% mcr carriage | 14.4 mean AMR genes | Cross-Inc lineage",
         font_size=12, color=DARK_GRAY)


# ======================================================================
# SLIDE 12: AMR Burden & Virulence (Figure 5)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "AMR Burden and Virulence Co-occurrence",
               "Quantifying resistance gene load and virulence factor co-carriage per Inc group")

add_figure(slide, "Figure5_AMR_burden_virulence.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(4.5))

# Highlights at bottom
add_box(slide, Inches(0.5), Inches(6.2), Inches(12.4), Inches(1.0), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.3), Inches(12), Inches(0.3),
         "Key Findings", font_size=14, bold=True, color=GREEN)
add_text(slide, Inches(0.7), Inches(6.65), Inches(12), Inches(0.4),
         "5 Inc groups with 100% AMR carriage  |  IncAC2: 12.6 mean AMR genes  |  IncFII: 23.2% virulence co-carriage",
         font_size=13, color=DARK_GRAY)


# ======================================================================
# SLIDE 13: FastANI Validation
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "FastANI Validation: 4,970 Pairs Across 20 Inc Groups",
               "Orthogonal validation of composition-based pLIN thresholds against true ANI")

# Key result boxes (top row)
ani_results = [
    ("Spearman \u03c1", "-0.348", "P < 10\u207b\u00b9\u2070\u2070", MED_BLUE),
    ("L6 (d\u22640.001)", "Median ANI", "99.9%", GREEN),
    ("Groups P<0.001", "15 of 20", "Inc groups significant", ORANGE),
    ("Strongest", "IncFIC", "\u03c1 = -0.883", PURPLE),
]
for i, (label, value, sub, color) in enumerate(ani_results):
    x = Inches(0.5 + i * 3.15)
    y = Inches(1.5)
    add_box(slide, x, y, Inches(2.9), Inches(1.6), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, x, y + Inches(0.1), Inches(2.9), Inches(0.3),
             label, font_size=12, color=MED_GRAY, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.45), Inches(2.9), Inches(0.5),
             value, font_size=22, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(1.0), Inches(2.9), Inches(0.3),
             sub, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Threshold vs observed ANI table
add_box(slide, Inches(0.5), Inches(3.5), Inches(12.4), Inches(0.5), DARK_BLUE)
table_headers = ["pLIN Level", "Distance Threshold", "Observed Median ANI", "N Pairs"]
for i, h in enumerate(table_headers):
    add_text(slide, Inches(0.6 + i * 3.1), Inches(3.53), Inches(3.0), Inches(0.4),
             h, font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

ani_rows = [
    ("L6 (Strain)", "d \u2264 0.001", "99.9%", "n = 726"),
    ("L5 (Clone complex)", "d \u2264 0.010", "98.2%", "n = 2,696"),
    ("L4 (Subcluster)", "d \u2264 0.020", "97.9%", "n = 3,566"),
    ("L3 (Cluster)", "d \u2264 0.050", "95.1%", "n = 4,208"),
    ("L2 (Subfamily)", "d \u2264 0.100", "91.3%", "n = 4,652"),
    ("L1 (Family)", "d \u2264 0.150", "86.7%", "n = 4,970"),
]
for j, (level, thresh, ani_val, n_pairs) in enumerate(ani_rows):
    y = Inches(4.0 + j * 0.5)
    bg = LIGHT_GRAY if j % 2 == 0 else WHITE
    add_box(slide, Inches(0.5), y, Inches(12.4), Inches(0.5), bg)
    vals = [level, thresh, ani_val, n_pairs]
    for i, v in enumerate(vals):
        add_text(slide, Inches(0.6 + i * 3.1), y + Inches(0.05), Inches(3.0), Inches(0.35),
                 v, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Bottom note
add_text(slide, Inches(0.5), Inches(7.0), Inches(12.4), Inches(0.3),
         "pLIN thresholds correspond closely to expected ANI ranges, validating composition-based distance as a proxy for genome similarity",
         font_size=11, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 14: Outbreak Cross-Validation (Figure 14)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Cross-Validation: 74 Plasmids from 27 Published Outbreak Studies",
               "Independent validation across 7 resistance mechanisms, 13 countries, 4 continents")

add_figure(slide, "Figure14_outbreak_validation.png",
           Inches(0.3), Inches(1.5), Inches(8.5), Inches(4.5))

# Key findings (right side)
add_box(slide, Inches(9.1), Inches(1.5), Inches(3.9), Inches(4.5), LIGHT_GRAY,
        border_color=DARK_BLUE, border_width=Pt(2))
add_text(slide, Inches(9.3), Inches(1.6), Inches(3.5), Inches(0.4),
         "Key Findings", font_size=16, bold=True, color=DARK_BLUE)
add_multiline(slide, Inches(9.3), Inches(2.1), Inches(3.5), Inches(3.5), [
    "3 global lineages matched:",
    "pLIN 671 (KPC-2, d=0.0000)",
    "pLIN 860 (MDR hub, d=0.0004)",
    "pLIN 1688 (OXA-48, 3 countries)",
    "",
    "9 intra-study clusters:",
    "12/13 HK NDM \u2192 pLIN 475",
    "5/6 OXA-48 \u2192 pLIN 1688",
    "4 mcr-1 \u2192 pLIN 87 (IncI2)",
    "",
    "85.1% high-confidence",
], font_size=12, color=DARK_GRAY)

# Bottom summary
add_box(slide, Inches(0.5), Inches(6.3), Inches(12.4), Inches(0.9), RGBColor(0xE3, 0xF2, 0xFD),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.45), Inches(12), Inches(0.5),
         "74 plasmids | 26 studies | 13 countries | 42 unique codes | 85.1% high-conf | 3 global lineage matches",
         font_size=14, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 15: Outbreak Validation Details
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Expanded Outbreak Validation \u2014 Key Results by Resistance Mechanism",
               "26 studies across 13 countries validate pLIN classification across 7 resistance mechanisms")

studies = [
    ("KPC (n=13, 92% high-conf)",
     "Germany (61 hospitals) + NIH USA + China + Italy + Brazil",
     "pLIN 671 hotspot confirmed across 2 continents; NIH KPC-3 \u2192 pLIN 672",
     RED),
    ("NDM (n=36, 83% high-conf)",
     "Hong Kong ICU + Germany polyclonal + Colombia + China",
     "12/13 HK plasmids share pLIN 475; 4 German NDM \u2192 pLIN 492",
     MED_BLUE),
    ("OXA-48 + mcr-1 (n=12, 92%)",
     "OXA-48: Turkey/Netherlands/France | mcr-1: China/Europe/USA",
     "5/6 OXA-48 \u2192 pLIN 1688 (3 countries); IncI2 vs IncX4 mcr separated",
     GREEN),
    ("CTX-M + VIM + IMP (n=13, 85%)",
     "CTX-M-15: USA/UK/India/France | VIM: Italy | IMP: Japan/Australia",
     "USA outbreak: 3 share pLIN 1482; pLIN 860 MDR hub match (Australia)",
     PURPLE),
]

for i, (title, context, result, color) in enumerate(studies):
    x = Inches(0.5) if i % 2 == 0 else Inches(6.8)
    y = Inches(1.5) if i < 2 else Inches(4.1)
    add_box(slide, x, y, Inches(6.0), Inches(2.3), WHITE, border_color=color, border_width=Pt(2))
    # Colored header band
    add_box(slide, x, y, Inches(6.0), Inches(0.5), color)
    add_text(slide, x + Inches(0.1), y + Inches(0.07), Inches(5.8), Inches(0.35),
             title, font_size=14, bold=True, color=WHITE)
    add_text(slide, x + Inches(0.2), y + Inches(0.6), Inches(5.6), Inches(0.4),
             context, font_size=11, color=DARK_GRAY)
    add_text(slide, x + Inches(0.2), y + Inches(1.1), Inches(5.6), Inches(0.4),
             "Result:", font_size=12, bold=True, color=color)
    add_text(slide, x + Inches(0.2), y + Inches(1.5), Inches(5.6), Inches(0.6),
             result, font_size=12, bold=True, color=DARK_GRAY)

# Summary box at bottom
add_box(slide, Inches(0.5), Inches(6.6), Inches(12.4), Inches(0.6), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.7), Inches(12), Inches(0.4),
         "74 plasmids | 26 studies | 13 countries | 42 unique codes | 9 intra-study clusters | 3 global lineages",
         font_size=14, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 16: Method Resolution Comparison (Figure 15)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Resolution Comparison: pLIN vs Existing Classification Methods",
               "Same 74 outbreak plasmids classified by 6 different methods")

add_figure(slide, "Figure15_method_resolution.png",
           Inches(0.15), Inches(1.4), Inches(9.5), Inches(5.0))

# Key message (right side)
add_box(slide, Inches(9.8), Inches(1.4), Inches(3.4), Inches(5.0), LIGHT_GRAY,
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(9.9), Inches(1.5), Inches(3.2), Inches(0.4),
         "Resolution Summary", font_size=14, bold=True, color=GREEN)
add_multiline(slide, Inches(9.9), Inches(2.0), Inches(3.2), Inches(4.0), [
    "PlasmidFinder: 11 groups",
    "(flat Inc labels only)",
    "",
    "pMLST: 10 types",
    "(39% lack scheme)",
    "",
    "MOB-suite: 8 clusters",
    "(39% unclassified)",
    "",
    "COPLA/PTU: 6 groups",
    "(43% unclassified)",
    "",
    "mge-cluster: 22 clusters",
    "(no permanent codes)",
    "",
    "pLIN: 42 unique codes",
    "(100% classified,",
    " hierarchical + permanent)",
], font_size=10, color=DARK_GRAY)

# Bottom highlight
add_box(slide, Inches(0.5), Inches(6.6), Inches(12.4), Inches(0.6), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.7), Inches(12), Inches(0.4),
         "pLIN provides 3.8\u00d7 more resolution than PlasmidFinder with 100% coverage and permanent hierarchical codes",
         font_size=13, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 17: Reference Database Expansion (Figure 13)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Scalability: 79,305 Plasmids, 57,886 Unique Codes",
               "Full-scale deployment on NCBI RefSeq plasmid database")

add_figure(slide, "Figure13_reference_database.png",
           Inches(0.3), Inches(1.5), Inches(9.0), Inches(4.5))

# Stats (right side)
scale_stats = [
    ("97.3%", "Classification\nrate", GREEN),
    ("<30 min", "Runtime on\nstandard laptop", MED_BLUE),
    ("2.3 GB", "Peak memory\nusage", ORANGE),
    ("18.8\u00d7", "Code expansion\n(3,073 \u2192 57,886)", PURPLE),
]
for i, (num, label, color) in enumerate(scale_stats):
    y = Inches(1.5 + i * 1.3)
    add_box(slide, Inches(9.6), y, Inches(3.4), Inches(1.1), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, Inches(9.8), y + Inches(0.05), Inches(3.0), Inches(0.4),
             num, font_size=22, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, Inches(9.8), y + Inches(0.5), Inches(3.0), Inches(0.5),
             label, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Bottom bar
add_box(slide, Inches(0.5), Inches(6.3), Inches(12.4), Inches(0.9), RGBColor(0xE3, 0xF2, 0xFD),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.45), Inches(12), Inches(0.5),
         "79,305 plasmids from NCBI RefSeq classified into 57,886 unique pLIN codes with 97.3% assignment rate",
         font_size=14, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 17: Clinical Application — Outbreak Detection
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Integrated Outbreak Detection Module",
               "Automated risk stratification for plasmid-mediated AMR dissemination")

# Two-tier system
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "Basic Detection", font_size=18, bold=True, color=MED_BLUE)

basic_items = [
    ("Same L6 code + Same AMR fingerprint", "Cluster identified"),
    ("\u22653 shared AMR genes", "HIGH risk"),
    ("<3 shared AMR genes", "MODERATE risk"),
]
for i, (cond, result) in enumerate(basic_items):
    y = Inches(2.0 + i * 0.55)
    add_box(slide, Inches(0.5), y, Inches(5.8), Inches(0.45), LIGHT_GRAY if i % 2 == 0 else WHITE)
    add_text(slide, Inches(0.7), y + Inches(0.05), Inches(3.8), Inches(0.35),
             cond, font_size=11, color=DARK_GRAY)
    add_text(slide, Inches(4.5), y + Inches(0.05), Inches(1.8), Inches(0.35),
             result, font_size=11, bold=True, color=RED if "HIGH" in result else MED_BLUE,
             alignment=PP_ALIGN.RIGHT)

add_text(slide, Inches(0.5), Inches(3.8), Inches(6), Inches(0.4),
         "Temporal Detection", font_size=18, bold=True, color=GREEN)

temporal_items = [
    ("\u22653 AMR genes + \u22647 days window", "CRITICAL risk"),
    ("\u22653 AMR genes + \u226430 days window", "HIGH risk"),
    ("<3 AMR genes + time window", "MODERATE risk"),
]
for i, (cond, result) in enumerate(temporal_items):
    y = Inches(4.3 + i * 0.55)
    add_box(slide, Inches(0.5), y, Inches(5.8), Inches(0.45), LIGHT_GRAY if i % 2 == 0 else WHITE)
    add_text(slide, Inches(0.7), y + Inches(0.05), Inches(3.8), Inches(0.35),
             cond, font_size=11, color=DARK_GRAY)
    color_r = RED if "CRITICAL" in result else (ORANGE if "HIGH" in result else MED_BLUE)
    add_text(slide, Inches(4.5), y + Inches(0.05), Inches(1.8), Inches(0.35),
             result, font_size=11, bold=True, color=color_r, alignment=PP_ALIGN.RIGHT)

# SNP sub-typing (right side)
add_text(slide, Inches(7.0), Inches(1.5), Inches(6), Inches(0.4),
         "SNP Sub-typing Within L6 Clusters", font_size=18, bold=True, color=PURPLE)

snp_tiers = [
    ("0 SNPs", "Clonal", "Identical backbone, recent transfer", RED),
    ("1\u20135 SNPs", "Recent divergence", "Likely shared origin within weeks", ORANGE),
    ("6\u201320 SNPs", "Indirect relationship", "Shared ancestor, diverging lineages", MED_BLUE),
    (">20 SNPs", "Independent", "Independently acquired, convergent", MED_GRAY),
]
for i, (snps, label, desc, color) in enumerate(snp_tiers):
    y = Inches(2.0 + i * 1.05)
    add_box(slide, Inches(7.0), y, Inches(5.8), Inches(0.9), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, Inches(7.2), y + Inches(0.05), Inches(1.3), Inches(0.35),
             snps, font_size=14, bold=True, color=color)
    add_text(slide, Inches(8.5), y + Inches(0.05), Inches(4.1), Inches(0.35),
             label, font_size=13, bold=True, color=DARK_GRAY)
    add_text(slide, Inches(8.5), y + Inches(0.4), Inches(4.1), Inches(0.35),
             desc, font_size=11, color=MED_GRAY)

# Bottom key message
add_box(slide, Inches(0.5), Inches(6.3), Inches(12.4), Inches(0.9), RGBColor(0xFF, 0xEB, 0xEE),
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.45), Inches(12), Inches(0.5),
         "Key question answered: Is resistance spreading on one plasmid or many?",
         font_size=16, bold=True, color=RED, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 18: Unique Advantages
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "What Makes pLIN Unique",
               "Six key innovations distinguishing pLIN from all existing tools")

advantages = [
    ("First LIN Application\nto Plasmids",
     "Life Identification Numbers (LIN) have been applied to bacteria and viruses, "
     "but pLIN is the first system adapting this framework for plasmid classification. "
     "Composition-based distances replace whole-genome alignment.",
     MED_BLUE),
    ("Permanent, Hierarchical\nCodes (6 Levels)",
     "pLIN codes are stable and never change when new sequences are added. "
     "The 6-level hierarchy (Family to Strain) provides resolution at multiple "
     "biological scales in a single unified code.",
     GREEN),
    ("Integrated AMR +\nVirulence + Stress",
     "Native integration with AMRFinderPlus enables direct mapping of 64,891 "
     "resistance gene detections, virulence factors, and stress genes onto "
     "specific plasmid lineages.",
     RED),
    ("Automated Outbreak\nDetection with Risk Tiers",
     "Two-tier outbreak detection (basic + temporal) with CRITICAL/HIGH/MODERATE "
     "risk stratification. SNP sub-typing provides sub-strain resolution for "
     "epidemiological investigations.",
     ORANGE),
    ("Cross-Inc Lineage\nDetection (pLIN 860)",
     "pLIN reveals lineages spanning multiple Inc groups, such as pLIN 860 "
     "(detected across 5 Inc groups with 44.4% mcr carriage). No other tool "
     "captures these cross-Inc relationships.",
     PURPLE),
    ("Cross-Validated Against\nReal-World Outbreaks",
     "74 plasmids from 27 published studies across 13 countries and "
     "7 resistance mechanisms validated pLIN. 3 globally disseminated "
     "lineages (671, 860, 1688) matched across continents.",
     TEAL),
]

for i, (title, desc, color) in enumerate(advantages):
    col = i % 3
    row = i // 3
    x = Inches(0.4 + col * 4.25)
    y = Inches(1.5 + row * 2.9)
    # Title bar
    add_box(slide, x, y, Inches(4.0), Inches(0.7), color)
    add_text(slide, x + Inches(0.1), y + Inches(0.05), Inches(3.8), Inches(0.6),
             title, font_size=13, bold=True, color=WHITE)
    # Description box
    add_box(slide, x, y + Inches(0.7), Inches(4.0), Inches(2.0), WHITE,
            border_color=color, border_width=Pt(1))
    add_text(slide, x + Inches(0.15), y + Inches(0.85), Inches(3.7), Inches(1.7),
             desc, font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 19: Phase 2 Features — Advanced Capabilities
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Phase 2 Features — Advanced Capabilities",
               "Additional modules extending pLIN beyond classification")

# 4 feature cards in 2x2 grid
phase2_features = [
    ("MOBsuite Mobility Typing", ORANGE,
     "3-tier priority cascade for mobility prediction",
     [
         "Priority 1: MOBsuite mob_typer (relaxase families: MOBF, MOBH, MOBP, MOBQ, MOBC, MOBV)",
         "Priority 2: AMRFinderPlus stress/virulence markers for transfer-associated genes",
         "Priority 3: Default non-mobilisable classification",
         "MPF type detection (MPFT, MPFF, MPFI, MPFG)",
         "Integrated conjugability + relaxase + MPF risk assessment",
     ]),
    ("CRISPR Spacer-Based Host Inference", TEAL,
     "Inferring plasmid-host relationships from CRISPR array matching",
     [
         "MinCED extracts CRISPR spacers from bacterial genomes",
         "BLASTn matches spacers against plasmid reference database",
         "Softmax probability ranking across candidate hosts",
         "Confidence categories: High (>=0.7), Moderate, Low",
         "Two modes: user-uploaded host genomes or reference DB",
     ]),
    ("Nucleotide Transformer LLM Integration", PURPLE,
     "Genomic foundation model for Inc group and AMR class prediction",
     [
         "InstaDeep Nucleotide Transformer (500M parameter model)",
         "6-mer tokenization of plasmid sequences (max 6 kb window)",
         "LLM-based Inc group prediction alongside traditional KNN",
         "AMR drug class prediction from sequence embeddings",
         "Cross-method validation: KNN vs LLM agreement scoring",
     ]),
    ("DRAGNOME Buddy — AI Assistant", RGBColor(0x00, 0x69, 0x5C),
     "Local LLM-powered chatbot for plasmid biology Q&A",
     [
         "Runs locally via Ollama (llama3.2, mistral) — no API keys needed",
         "Context-aware: feeds pLIN results into LLM for interpretation",
         "Answers questions about Inc groups, AMR genes, plasmid biology",
         "Generates clinical interpretation summaries",
         "All data stays local — no cloud upload required",
     ]),
]

for idx, (title, color, subtitle, items) in enumerate(phase2_features):
    col = idx % 2
    row = idx // 2
    x = Inches(0.4 + col * 6.5)
    y = Inches(1.5 + row * 2.9)

    # Title bar
    add_box(slide, x, y, Inches(6.1), Inches(0.55), color)
    add_text(slide, x + Inches(0.15), y + Inches(0.05), Inches(5.8), Inches(0.45),
             title, font_size=15, bold=True, color=WHITE)

    # Content box
    add_box(slide, x, y + Inches(0.55), Inches(6.1), Inches(2.2), WHITE,
            border_color=color, border_width=Pt(1))
    add_text(slide, x + Inches(0.15), y + Inches(0.6), Inches(5.8), Inches(0.3),
             subtitle, font_size=11, bold=True, color=color)
    add_multiline(slide, x + Inches(0.15), y + Inches(0.95), Inches(5.8), Inches(1.7),
                  items, font_size=10, color=DARK_GRAY)


# ======================================================================
# SLIDE 20: Limitations & Future Directions
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Limitations & Future Directions",
               "Current constraints and planned improvements")

# Left: Limitations
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.5),
         "Current Limitations", font_size=20, bold=True, color=RED)

limitations = [
    ("Uneven Inc Group Training Sizes",
     "IncFII has 4,629 samples while IncFIBK has only 11. Smaller classes may have lower accuracy."),
    ("Composition-Only Features",
     "4-mer captures global composition but ignores gene content, synteny, and structural rearrangements."),
    ("Single-Linkage Chaining",
     "Single-linkage clustering can produce chain-like clusters at coarser thresholds. Alternative linkage methods may require re-calibration."),
    ("Requires Complete Plasmid Sequences",
     "Short contigs from fragmented assemblies have noisy k-mer profiles. No multi-contig reconstruction."),
    ("Static Training Set",
     "KNN classifier uses a fixed reference. New NCBI submissions require periodic retraining."),
    ("Chromosomal Sequence Sensitivity",
     "Whole-genome uploads may include chromosomal contigs (>500 kb) that receive spurious classifications."),
]

for i, (title, desc) in enumerate(limitations):
    y = Inches(2.1 + i * 0.82)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.05), Inches(0.3), Inches(0.3))
    circ.fill.solid()
    circ.fill.fore_color.rgb = RED
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = str(i + 1)
    tf.paragraphs[0].font.size = Pt(11)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    add_text(slide, Inches(0.95), y, Inches(5.8), Inches(0.3),
             title, font_size=12, bold=True, color=RED)
    add_text(slide, Inches(0.95), y + Inches(0.3), Inches(5.8), Inches(0.45),
             desc, font_size=9, color=DARK_GRAY)

# Right: Future Directions
add_text(slide, Inches(7.0), Inches(1.5), Inches(6), Inches(0.5),
         "Future Directions", font_size=20, bold=True, color=GREEN)

futures = [
    ("Expand Inc Group Coverage",
     "Target 30+ Inc groups (IncP, IncW, IncL/M) with 20,000+ balanced training sequences.",
     GREEN),
    ("Hybrid Features",
     "Combine 4-mer composition with gene presence/absence and synteny for higher resolution.",
     MED_BLUE),
    ("Metagenomic Support",
     "Handle multi-contig plasmid bins from metagenomic assemblies (PlasFlow, MOB-recon).",
     PURPLE),
    ("Online Database & API",
     "Web database for pLIN lookups and REST API for LIMS integration.",
     TEAL),
    ("Multi-Hospital Surveillance",
     "Cross-institutional outbreak tracking with spatial-temporal epidemiological metadata.",
     RED),
    ("WHO GLASS Integration",
     "Standardised plasmid nomenclature for national and international AMR surveillance networks.",
     ORANGE),
]

for i, (title, desc, color) in enumerate(futures):
    y = Inches(2.1 + i * 0.82)
    arrow = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(7.0), y + Inches(0.1), Inches(0.35), Inches(0.25))
    arrow.fill.solid()
    arrow.fill.fore_color.rgb = color
    arrow.line.fill.background()
    add_text(slide, Inches(7.5), y, Inches(5.3), Inches(0.3),
             title, font_size=12, bold=True, color=color)
    add_text(slide, Inches(7.5), y + Inches(0.3), Inches(5.3), Inches(0.45),
             desc, font_size=9, color=DARK_GRAY)


# ======================================================================
# SLIDE 21: Roadmap — Phase 1 (Complete) vs Phase 2 (Planned)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Development Roadmap — v2.1 Feature Matrix",
               "Phase 1 features (shipped) and Phase 2 features (planned)")

# Phase 1: Complete
add_box(slide, Inches(0.4), Inches(1.5), Inches(6.2), Inches(0.5), GREEN)
add_text(slide, Inches(0.5), Inches(1.55), Inches(6.0), Inches(0.4),
         "Phase 1 — Shipped (v2.1)", font_size=16, bold=True, color=WHITE)

phase1_items = [
    ("1", "Adaptive Calibration Default", "Per-Inc-group thresholds enabled by default for all 20 groups", GREEN),
    ("2", "Sequence Length Warnings", "Flags <5 kb (noisy) and >500 kb (chromosomal) sequences", ORANGE),
    ("3", "Metadata CSV Upload", "Patient/sample collection dates for epidemiological analysis", MED_BLUE),
    ("4", "Mash/MinHash ANI", "Fast approximate ANI validation (k=21, s=10,000)", PURPLE),
    ("5", "Single-Plasmid Query Mode", "Classify individual plasmids against the 8,077 reference database", TEAL),
    ("6", "Docker Deployment", "One-command containerised deployment with all biotools included", RGBColor(0x00, 0x69, 0x5C)),
]

for i, (num, title, desc, color) in enumerate(phase1_items):
    y = Inches(2.15 + i * 0.52)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.03), Inches(0.28), Inches(0.28))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(10)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    add_text(slide, Inches(0.9), y, Inches(2.5), Inches(0.25),
             title, font_size=11, bold=True, color=color)
    add_text(slide, Inches(0.9), y + Inches(0.24), Inches(5.6), Inches(0.22),
             desc, font_size=9, color=MED_GRAY)

# Phase 2: Planned
add_box(slide, Inches(0.4), Inches(5.4), Inches(6.2), Inches(0.5), RED)
add_text(slide, Inches(0.5), Inches(5.45), Inches(6.0), Inches(0.4),
         "Phase 2 — Advanced Genomic Resolution (Planned)", font_size=16, bold=True, color=WHITE)

phase2_items = [
    ("7", "FastANI True ANI", "Ground-truth ANI validation with --fragLen 1000 (v1.34)", MED_BLUE),
    ("8", "SNP Sub-typing (L6)", "minimap2 -cx asm5 within L6 clusters for outbreak confirmation", GREEN),
    ("9", "Temporal Outbreak Clustering", "30-day sliding window + 3-tier risk (CRITICAL/HIGH/MODERATE)", RED),
]

for i, (num, title, desc, color) in enumerate(phase2_items):
    y = Inches(6.05 + i * 0.45)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.03), Inches(0.28), Inches(0.28))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(10)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    add_text(slide, Inches(0.9), y, Inches(2.5), Inches(0.25),
             title, font_size=11, bold=True, color=color)
    add_text(slide, Inches(0.9), y + Inches(0.24), Inches(5.6), Inches(0.22),
             desc, font_size=9, color=MED_GRAY)

# Right side: Feature summary matrix
add_text(slide, Inches(7.0), Inches(1.5), Inches(6), Inches(0.5),
         "Impact & Tool Requirements", font_size=18, bold=True, color=DARK_BLUE)

# Table header
add_box(slide, Inches(7.0), Inches(2.1), Inches(5.9), Inches(0.4), DARK_BLUE)
imp_headers = ["Feature", "Tool", "Status", "Impact"]
imp_widths = [1.8, 1.2, 0.9, 2.0]
x = Inches(7.1)
for h, w in zip(imp_headers, imp_widths):
    add_text(slide, x, Inches(2.15), Inches(w), Inches(0.3),
             h, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    x += Inches(w)

impact_rows = [
    ("Adaptive Calibration", "Built-in", "Done", "Per-Inc accuracy"),
    ("Length Warnings", "Built-in", "Done", "User safety"),
    ("Metadata Upload", "Built-in", "Done", "Epi context"),
    ("Mash ANI", "mash", "Done", "ANI validation"),
    ("Single-Query Mode", "Built-in", "Done", "Clinical use"),
    ("Docker Deploy", "Docker", "Done", "Easy setup"),
    ("FastANI", "fastANI", "Planned", "True ANI"),
    ("SNP Sub-typing", "minimap2", "Planned", "Outbreak resolution"),
    ("Temporal Clusters", "Built-in", "Planned", "Surveillance"),
    ("MOBsuite", "mob_typer", "Planned", "Mobility typing"),
    ("CRISPR Hosts", "MinCED", "Planned", "Host inference"),
    ("NT LLM", "PyTorch", "Planned", "LLM validation"),
    ("DRAGNOME Buddy", "Ollama", "Planned", "AI assistant"),
]

for i, (feat, tool, status, impact) in enumerate(impact_rows):
    y = Inches(2.55 + i * 0.36)
    bg = LIGHT_GRAY if i % 2 == 0 else WHITE
    add_box(slide, Inches(7.0), y, Inches(5.9), Inches(0.36), bg)
    x = Inches(7.1)
    for j, (val, w) in enumerate(zip([feat, tool, status, impact], imp_widths)):
        if val == "Built-in":
            c = GREEN
        elif val == "Done":
            c = GREEN
        elif val == "Planned":
            c = ORANGE
        else:
            c = DARK_GRAY
        add_text(slide, x, y + Inches(0.05), Inches(w), Inches(0.22),
                 val, font_size=9, color=c, alignment=PP_ALIGN.CENTER)
        x += Inches(w)

# Version badge
add_box(slide, Inches(7.0), Inches(7.3 - 0.6), Inches(5.9), Inches(0.5), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.2), Inches(7.3 - 0.55), Inches(5.5), Inches(0.4),
         "pLIN v2.1 — 6 shipped + 7 planned = 13 total features",
         font_size=13, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 22: Gram-Positive Plasmid Expansion
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Gram-Positive Plasmid Expansion",
               "Extending pLIN from 20 Gram-negative Inc groups to 28 groups with Gram-positive, Acinetobacter, and Pseudomonas rep types")

# Figure (left side)
add_figure(slide, "figure17_gram_positive_expansion.png",
           Inches(0.3), Inches(1.5), Inches(7.5), Inches(5.5))

# Stats panel (right side)
add_box(slide, Inches(8.2), Inches(1.5), Inches(4.8), Inches(2.8), LIGHT_GRAY,
        border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(1.6), Inches(4.4), Inches(0.4),
         "Expansion Summary", font_size=16, bold=True, color=PURPLE)

expansion_stats = [
    ("20 \u2192 28", "Inc/Rep Groups", MED_BLUE),
    ("8,077", "Training Sequences", GREEN),
    ("91.1%", "KNN Accuracy", ORANGE),
]
for i, (num, label, color) in enumerate(expansion_stats):
    y = Inches(2.2 + i * 0.75)
    add_text(slide, Inches(8.4), y, Inches(4.4), Inches(0.35),
             num, font_size=22, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, Inches(8.4), y + Inches(0.35), Inches(4.4), Inches(0.3),
             label, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# New rep type groups box
add_box(slide, Inches(8.2), Inches(4.6), Inches(4.8), Inches(2.4), WHITE,
        border_color=TEAL, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(4.7), Inches(4.4), Inches(0.4),
         "4 New Rep Type Groups", font_size=14, bold=True, color=TEAL)

new_groups = [
    ("repSA_large", "S. aureus large plasmids (conjugative)"),
    ("repSA_small", "S. aureus small plasmids (mobilisable)"),
    ("repEF_conj", "E. faecium/faecalis conjugative plasmids"),
    ("repEF_res", "E. faecium/faecalis resistance plasmids"),
]
for i, (grp, desc) in enumerate(new_groups):
    y = Inches(5.2 + i * 0.42)
    add_text(slide, Inches(8.5), y, Inches(1.8), Inches(0.35),
             grp, font_size=11, bold=True, color=TEAL, font_name="Consolas")
    add_text(slide, Inches(10.3), y, Inches(2.5), Inches(0.35),
             desc, font_size=10, color=DARK_GRAY)


# ======================================================================
# SLIDE 23: Assembly Completeness (L3)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L3 \u2014 Assembly Completeness Assessment",
               "Composite quality score (0\u2013100) for plasmid assembly completeness")

# 5 criteria boxes across the top
criteria_items = [
    ("Contig Count", "Single vs\nmulti-contig", MED_BLUE),
    ("N50 Ratio", "N50 / total\nassembly length", GREEN),
    ("Circular Signal", "Overlap detection\nfor circular closure", ORANGE),
    ("Coding Density", "CDS coverage\nrelative to length", PURPLE),
    ("N-Gaps", "Ambiguous base\ncount (N content)", RED),
]
for i, (title, desc, color) in enumerate(criteria_items):
    x = Inches(0.3 + i * 2.55)
    y = Inches(1.6)
    add_box(slide, x, y, Inches(2.35), Inches(0.5), color)
    add_text(slide, x + Inches(0.05), y + Inches(0.07), Inches(2.25), Inches(0.35),
             title, font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, y + Inches(0.5), Inches(2.35), Inches(1.0), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.1), y + Inches(0.6), Inches(2.15), Inches(0.8),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)

# Arrow down to composite score
add_arrow(slide, Inches(6.3), Inches(3.3), Inches(0.5), Inches(0.5), MED_BLUE)

# Composite score box
add_box(slide, Inches(3.5), Inches(4.0), Inches(6.3), Inches(0.7),
        RGBColor(0xE3, 0xF2, 0xFD), border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(3.7), Inches(4.1), Inches(5.9), Inches(0.5),
         "Composite Score = weighted sum of 5 criteria \u2192 0\u2013100 scale",
         font_size=15, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

# 4 tier boxes
tiers = [
    ("COMPLETE", "Score \u2265 90", "Closed, single-contig,\ncircular plasmid", GREEN),
    ("NEAR-COMPLETE", "Score 70\u201389", "Minor fragmentation,\nhigh coding density", MED_BLUE),
    ("FRAGMENTED", "Score 40\u201369", "Multiple contigs,\nmoderate N50 ratio", ORANGE),
    ("POOR", "Score < 40", "Highly fragmented,\nlow coding density", RED),
]
for i, (tier, score, desc, color) in enumerate(tiers):
    x = Inches(0.4 + i * 3.2)
    y = Inches(5.0)
    add_box(slide, x, y, Inches(3.0), Inches(0.45), color)
    add_text(slide, x + Inches(0.05), y + Inches(0.05), Inches(2.9), Inches(0.35),
             tier, font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, y + Inches(0.45), Inches(3.0), Inches(1.4), WHITE,
            border_color=color, border_width=Pt(1))
    add_text(slide, x + Inches(0.1), y + Inches(0.5), Inches(2.8), Inches(0.3),
             score, font_size=13, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_multiline(slide, x + Inches(0.1), y + Inches(0.85), Inches(2.8), Inches(0.9),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)

# Bottom note
add_text(slide, Inches(0.5), Inches(7.0), Inches(12.4), Inches(0.3),
         "Assembly completeness tiers inform confidence of downstream pLIN classification and AMR detection",
         font_size=12, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 24: Plasmid Contig Identification
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Plasmid Contig Identification",
               "Automatic classification of plasmid vs chromosomal contigs")

# --- Section title: 4 Scoring Signals ---
add_text(slide, Inches(0.5), Inches(1.4), Inches(12.4), Inches(0.4),
         "Multi-Signal Scoring System — 4 Complementary Signals",
         font_size=18, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

# --- Four signal boxes across ---
signal_items = [
    ("Sequence Length", [
        ">1 Mb:  \u221250 pts",
        ">500 kb: \u221230 pts",
        "<300 kb: +20 pts",
        "<20 kb:  +10 pts",
    ], MED_BLUE),
    ("4-mer Distance", [
        "Cosine distance to",
        "nearest training",
        "plasmid vector in",
        "256-dim k-mer space",
    ], GREEN),
    ("Header Keywords", [
        "Detect plasmid /",
        "chromosome keywords",
        "in FASTA headers",
        "(regex-based scan)",
    ], ORANGE),
    ("Inc Confidence", [
        "KNN classifier",
        "confidence score",
        "for Inc/Rep group",
        "(k=5, cosine dist.)",
    ], PURPLE),
]
for i, (title, desc_lines, color) in enumerate(signal_items):
    x = Inches(0.3 + i * 3.2)
    y = Inches(2.0)
    # Colored header strip
    add_box(slide, x, y, Inches(3.0), Inches(0.5), color)
    add_text(slide, x + Inches(0.05), y + Inches(0.07), Inches(2.9), Inches(0.35),
             title, font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    # Body box
    add_box(slide, x, y + Inches(0.5), Inches(3.0), Inches(1.4), WHITE,
            border_color=color, border_width=Pt(2))
    add_multiline(slide, x + Inches(0.15), y + Inches(0.6), Inches(2.7), Inches(1.2),
                  desc_lines, font_size=11, color=DARK_GRAY)

# Arrow from signals to classification
add_arrow(slide, Inches(6.3), Inches(4.1), Inches(0.5), Inches(0.5), MED_BLUE)

# --- Composite scoring box ---
add_box(slide, Inches(2.5), Inches(4.75), Inches(8.3), Inches(0.55),
        RGBColor(0xE3, 0xF2, 0xFD), border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(2.7), Inches(4.82), Inches(7.9), Inches(0.4),
         "Composite Score = weighted sum of 4 signals \u2192 classification decision",
         font_size=14, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)

# --- Three classification outcome boxes ---
class_items = [
    ("PLASMID", "\u2265 95% confidence", "pLIN code + all modules\n(AMR, mobility, CRISPR)",
     GREEN, RGBColor(0xE8, 0xF5, 0xE9)),
    ("INCOMPLETE PLASMID", "< 95% confidence", "AMR/mobility only\n(no pLIN assignment)",
     ORANGE, RGBColor(0xFF, 0xF3, 0xE0)),
    ("CHROMOSOME", "excluded", "Excluded entirely\nfrom plasmid analysis",
     RED, RGBColor(0xFF, 0xEB, 0xEE)),
]
for i, (label, conf, desc, color, bg_color) in enumerate(class_items):
    x = Inches(0.4 + i * 4.3)
    y = Inches(5.55)
    # Colored header
    add_box(slide, x, y, Inches(3.9), Inches(0.45), color)
    add_text(slide, x + Inches(0.05), y + Inches(0.05), Inches(3.8), Inches(0.35),
             label, font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    # Detail box with tinted background
    add_box(slide, x, y + Inches(0.45), Inches(3.9), Inches(1.2), bg_color,
            border_color=color, border_width=Pt(2))
    add_text(slide, x + Inches(0.15), y + Inches(0.55), Inches(3.6), Inches(0.3),
             conf, font_size=13, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_multiline(slide, x + Inches(0.15), y + Inches(0.9), Inches(3.6), Inches(0.7),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)

# Bottom note
add_text(slide, Inches(0.5), Inches(7.1), Inches(12.4), Inches(0.3),
         "Multi-signal scoring system combining 4 complementary signals to reliably separate plasmid from chromosomal contigs",
         font_size=12, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 25: Database Coverage & Novelty Detection (L4)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L4 \u2014 Database Coverage & Novelty Detection",
               "Nearest-neighbour distance percentile with traffic-light confidence system")

# Left: NN distance explanation
add_text(slide, Inches(0.5), Inches(1.5), Inches(6.0), Inches(0.4),
         "Nearest-Neighbour Distance Percentile", font_size=18, bold=True, color=DARK_BLUE)

add_box(slide, Inches(0.5), Inches(2.0), Inches(6.0), Inches(1.5), LIGHT_GRAY,
        border_color=MED_BLUE, border_width=Pt(2))
add_multiline(slide, Inches(0.7), Inches(2.1), Inches(5.6), Inches(1.3), [
    "For each query plasmid, compute cosine distance to its nearest",
    "neighbour in the reference database. Rank this distance against",
    "the distribution of all pairwise NN distances within the assigned",
    "Inc group to produce a percentile score (0\u2013100th).",
], font_size=12, color=DARK_GRAY)

# Traffic light system
add_text(slide, Inches(0.5), Inches(3.8), Inches(6.0), Inches(0.4),
         "Traffic-Light Confidence System", font_size=18, bold=True, color=DARK_BLUE)

traffic_lights = [
    ("GREEN", "\u2264 75th percentile", "Well-represented in reference DB;\nhigh-confidence classification", GREEN),
    ("YELLOW", "75th\u201395th percentile", "Moderate novelty; classification valid\nbut underrepresented region", ORANGE),
    ("RED", "> 95th percentile", "Potential novel lineage; distant from\nall known references \u2192 flag for review", RED),
]
for i, (light, thresh, desc, color) in enumerate(traffic_lights):
    y = Inches(4.3 + i * 1.0)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.1),
                                  Inches(0.4), Inches(0.4))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    add_text(slide, Inches(1.1), y, Inches(1.5), Inches(0.35),
             light, font_size=14, bold=True, color=color)
    add_text(slide, Inches(2.5), y, Inches(2.0), Inches(0.35),
             thresh, font_size=12, color=DARK_GRAY, font_name="Consolas")
    add_multiline(slide, Inches(2.5), y + Inches(0.35), Inches(4.0), Inches(0.6),
                  desc.split("\n"), font_size=11, color=MED_GRAY)

# Right: Novelty flagging
add_text(slide, Inches(7.0), Inches(1.5), Inches(6.0), Inches(0.4),
         "Novelty Flagging", font_size=18, bold=True, color=PURPLE)

add_box(slide, Inches(7.0), Inches(2.0), Inches(5.8), Inches(5.0), WHITE,
        border_color=PURPLE, border_width=Pt(2))

novelty_items = [
    "Query plasmids exceeding the 95th percentile",
    "NN distance are flagged as PUTATIVE NOVEL",
    "",
    "These may represent:",
    "\u2022 Genuinely novel plasmid lineages",
    "\u2022 Recombinant / chimeric plasmids",
    "\u2022 Under-sampled Inc group diversity",
    "\u2022 Misassembled sequences",
    "",
    "Action: Flagged plasmids are routed to",
    "the L6 novel group discovery module",
    "for hierarchical clustering analysis",
    "",
    "Integration: RED-flagged plasmids with",
    "low assembly completeness (L3 POOR)",
    "are deprioritised automatically",
]
add_multiline(slide, Inches(7.2), Inches(2.1), Inches(5.4), Inches(4.8),
              novelty_items, font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 25: Recombination Detection (L5)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L5 \u2014 Recombination Detection",
               "minimap2 PAF alignment analysis for identifying mosaic plasmid architectures")

# Figure (left)
add_figure(slide, "figure18_recombination_evolutionary_rate.png",
           Inches(0.3), Inches(1.5), Inches(7.5), Inches(5.5))

# Method panel (right)
add_box(slide, Inches(8.2), Inches(1.5), Inches(4.8), Inches(2.5), LIGHT_GRAY,
        border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(1.6), Inches(4.4), Inches(0.4),
         "Alignment Analysis", font_size=16, bold=True, color=ORANGE)
add_multiline(slide, Inches(8.4), Inches(2.1), Inches(4.4), Inches(1.8), [
    "minimap2 -cx asm5 aligns query",
    "vs nearest reference plasmid",
    "",
    "Three metrics computed from PAF:",
    "\u2022 Alignment coverage (%)",
    "\u2022 Number of alignment blocks",
    "\u2022 Inter-block gap lengths",
], font_size=11, color=DARK_GRAY)

# Recombination tiers
add_box(slide, Inches(8.2), Inches(4.3), Inches(4.8), Inches(2.8), WHITE,
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(4.4), Inches(4.4), Inches(0.4),
         "Recombination Flags", font_size=14, bold=True, color=RED)

recomb_tiers = [
    ("None", "\u22651 block, >95% cov", GREEN),
    ("Low", "2\u20133 blocks, >80% cov", MED_BLUE),
    ("Medium", "4\u20136 blocks or gaps >5 kb", ORANGE),
    ("High", ">6 blocks, <70% cov, large gaps", RED),
]
for i, (level, criteria, color) in enumerate(recomb_tiers):
    y = Inches(4.9 + i * 0.5)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(8.4), y + Inches(0.05),
                                  Inches(0.25), Inches(0.25))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    add_text(slide, Inches(8.8), y, Inches(1.2), Inches(0.35),
             level, font_size=12, bold=True, color=color)
    add_text(slide, Inches(10.0), y, Inches(2.8), Inches(0.35),
             criteria, font_size=10, color=DARK_GRAY)


# ======================================================================
# SLIDE 26: Novel Inc Group Discovery (L6)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L6 \u2014 Novel Inc Group Discovery",
               "Hierarchical clustering of low-confidence plasmids to identify putative novel groups")

# Left: Method description
add_text(slide, Inches(0.5), Inches(1.5), Inches(6.0), Inches(0.4),
         "Discovery Pipeline", font_size=18, bold=True, color=DARK_BLUE)

# Pipeline steps
discovery_steps = [
    ("1. Collect", "Gather all plasmids flagged RED by\nthe L4 novelty detection module", MED_BLUE),
    ("2. Cluster", "Hierarchical clustering (cosine distance,\nsingle-linkage) on 4-mer frequency vectors", GREEN),
    ("3. Evaluate", "Apply silhouette score and cluster\nsize thresholds to identify stable groups", ORANGE),
    ("4. Assign", "Putative novel groups receive provisional\nInc designations pending wet-lab validation", PURPLE),
]
for i, (step, desc, color) in enumerate(discovery_steps):
    y = Inches(2.0 + i * 1.2)
    add_box(slide, Inches(0.5), y, Inches(1.3), Inches(0.45), color)
    add_text(slide, Inches(0.55), y + Inches(0.05), Inches(1.2), Inches(0.35),
             step, font_size=12, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_multiline(slide, Inches(2.0), y + Inches(0.05), Inches(4.5), Inches(0.9),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)
    if i < len(discovery_steps) - 1:
        add_arrow(slide, Inches(1.0), y + Inches(0.55), Inches(0.3), Inches(0.45), color)

# Right: Putative novel groups panel
add_text(slide, Inches(7.0), Inches(1.5), Inches(6.0), Inches(0.4),
         "Putative Novel Groups", font_size=18, bold=True, color=PURPLE)

add_box(slide, Inches(7.0), Inches(2.0), Inches(5.8), Inches(4.8), LIGHT_GRAY,
        border_color=PURPLE, border_width=Pt(2))
add_multiline(slide, Inches(7.2), Inches(2.1), Inches(5.4), Inches(4.6), [
    "Plasmids that cannot be confidently assigned",
    "to any of the 24 known Inc/Rep groups are",
    "pooled and clustered independently.",
    "",
    "Criteria for novel group designation:",
    "\u2022 Cluster size \u2265 5 plasmids",
    "\u2022 Silhouette score \u2265 0.3",
    "\u2022 Intra-cluster distance < L1 threshold",
    "\u2022 No overlap with existing Inc groups",
    "",
    "Output:",
    "\u2022 Provisional group ID (e.g., Novel_A, Novel_B)",
    "\u2022 Representative sequence for each group",
    "\u2022 Recommended for replicon typing / PCR validation",
    "",
    "Goal: Systematic expansion of pLIN coverage",
    "to capture emerging plasmid diversity",
], font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 27: Evolutionary Rate Estimation (L7)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L7 \u2014 Evolutionary Rate Estimation",
               "SNP accumulation vs time regression for molecular clock calibration")

# Figure (left)
add_figure(slide, "figure18_recombination_evolutionary_rate.png",
           Inches(0.3), Inches(1.5), Inches(7.5), Inches(5.5))

# Method panel (right)
add_box(slide, Inches(8.2), Inches(1.5), Inches(4.8), Inches(2.5), LIGHT_GRAY,
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(1.6), Inches(4.4), Inches(0.4),
         "Method", font_size=16, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(8.4), Inches(2.1), Inches(4.4), Inches(1.8), [
    "Within each L5/L6 cluster with",
    "temporal metadata (collection dates):",
    "",
    "1. Compute pairwise SNP distances",
    "2. Compute pairwise time differences",
    "3. Linear regression: SNPs ~ time",
], font_size=11, color=DARK_GRAY)

# Output panel
add_box(slide, Inches(8.2), Inches(4.3), Inches(4.8), Inches(2.8), WHITE,
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(4.4), Inches(4.4), Inches(0.4),
         "Output Metrics", font_size=14, bold=True, color=GREEN)
add_multiline(slide, Inches(8.4), Inches(4.9), Inches(4.4), Inches(2.0), [
    "\u2022 Substitution rate (subs/site/year)",
    "\u2022 R\u00b2 for temporal signal strength",
    "\u2022 95% confidence interval",
    "",
    "Application:",
    "\u2022 Estimate TMRCA within clusters",
    "\u2022 Calibrate outbreak timelines",
    "\u2022 Compare rates across Inc groups",
    "\u2022 Detect rate heterogeneity",
], font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 28: Cluster Stability Assessment (L8)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L8 \u2014 Cluster Stability Assessment",
               "Bootstrap resampling and adaptive thresholds for robust clustering")

# Figure (left)
add_figure(slide, "figure19_cluster_stability.png",
           Inches(0.3), Inches(1.5), Inches(7.5), Inches(5.5))

# Bootstrap method (right top)
add_box(slide, Inches(8.2), Inches(1.5), Inches(4.8), Inches(2.3), LIGHT_GRAY,
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(1.6), Inches(4.4), Inches(0.4),
         "Bootstrap Resampling", font_size=16, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(8.4), Inches(2.1), Inches(4.4), Inches(1.5), [
    "Resample 4-mer feature vectors with",
    "replacement (N=100 iterations)",
    "",
    "For each bootstrap replicate:",
    "\u2022 Recompute cosine distance matrix",
    "\u2022 Re-cluster with single-linkage",
    "\u2022 Compare cluster assignments to",
    "  original using ARI",
], font_size=11, color=DARK_GRAY)

# ARI and adaptive thresholds (right bottom)
add_box(slide, Inches(8.2), Inches(4.1), Inches(4.8), Inches(1.4), WHITE,
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(4.2), Inches(4.4), Inches(0.35),
         "ARI Linkage Comparison", font_size=14, bold=True, color=GREEN)
add_multiline(slide, Inches(8.4), Inches(4.6), Inches(4.4), Inches(0.8), [
    "Adjusted Rand Index (ARI) measures",
    "agreement between original and",
    "bootstrap cluster assignments",
], font_size=11, color=DARK_GRAY)

add_box(slide, Inches(8.2), Inches(5.7), Inches(4.8), Inches(1.4), WHITE,
        border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(8.4), Inches(5.8), Inches(4.4), Inches(0.35),
         "Adaptive Thresholds", font_size=14, bold=True, color=PURPLE)
add_multiline(slide, Inches(8.4), Inches(6.2), Inches(4.4), Inches(0.8), [
    "Per-Inc-group threshold adjustment",
    "based on bootstrap stability analysis;",
    "groups with low ARI receive tighter",
    "distance thresholds automatically",
], font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 29: MGE Boundaries (L10)
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "L10 \u2014 Mobile Genetic Element Boundary Detection",
               "IS elements, composite transposons, and color-coded gene maps")

# Left: Method
add_text(slide, Inches(0.5), Inches(1.5), Inches(6.0), Inches(0.4),
         "MGE Detection Pipeline", font_size=18, bold=True, color=DARK_BLUE)

mge_steps = [
    ("IS Element Detection", "ISfinder-based HMM scan identifies insertion\nsequences (IS1, IS26, IS903, ISEcp1, etc.)", RED),
    ("Composite Transposon\nIdentification", "Paired IS elements flanking AMR gene\ncassettes define composite transposons", ORANGE),
    ("Gene Map Annotation", "Prokka/Bakta CDS annotation combined\nwith AMR and IS element coordinates", MED_BLUE),
    ("Boundary Visualisation", "Color-coded linear gene maps showing\nIS, AMR, backbone, and hypothetical genes", GREEN),
]
for i, (title, desc, color) in enumerate(mge_steps):
    y = Inches(2.0 + i * 1.2)
    add_box(slide, Inches(0.5), y, Inches(5.8), Inches(1.05), WHITE,
            border_color=color, border_width=Pt(2))
    add_box(slide, Inches(0.5), y, Inches(0.08), Inches(1.05), color)
    add_text(slide, Inches(0.8), y + Inches(0.05), Inches(2.3), Inches(0.45),
             title, font_size=12, bold=True, color=color)
    add_multiline(slide, Inches(3.1), y + Inches(0.05), Inches(3.0), Inches(0.9),
                  desc.split("\n"), font_size=11, color=DARK_GRAY)

# Right: Gene map legend and output
add_text(slide, Inches(7.0), Inches(1.5), Inches(6.0), Inches(0.4),
         "Color-Coded Gene Map Legend", font_size=18, bold=True, color=DARK_BLUE)

gene_colors = [
    ("AMR genes", "Red arrows", RED),
    ("IS elements", "Orange arrows", ORANGE),
    ("Virulence factors", "Purple arrows", PURPLE),
    ("Backbone / replication", "Blue arrows", MED_BLUE),
    ("Hypothetical proteins", "Grey arrows", MED_GRAY),
    ("Transfer / conjugation", "Green arrows", GREEN),
]
for i, (gene_type, arrow_desc, color) in enumerate(gene_colors):
    y = Inches(2.0 + i * 0.65)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(7.0), y + Inches(0.05),
                                  Inches(0.3), Inches(0.3))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    add_text(slide, Inches(7.5), y, Inches(2.5), Inches(0.35),
             gene_type, font_size=12, bold=True, color=color)
    add_text(slide, Inches(10.0), y, Inches(2.8), Inches(0.35),
             arrow_desc, font_size=11, color=DARK_GRAY)

# Output box
add_box(slide, Inches(7.0), Inches(6.0), Inches(5.8), Inches(1.0),
        RGBColor(0xE3, 0xF2, 0xFD), border_color=MED_BLUE, border_width=Pt(2))
add_multiline(slide, Inches(7.2), Inches(6.1), Inches(5.4), Inches(0.8), [
    "Output: SVG/PNG gene maps per plasmid with IS boundaries,",
    "composite transposon annotations, and AMR gene cassettes highlighted",
], font_size=11, bold_first=True, color=DARK_BLUE)


# ======================================================================
# SLIDE 30: Limitations Addressed Summary
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Limitations Addressed \u2014 Comprehensive Overview",
               "10 original limitations: 3 addressed by MLST, 7 by new modules; 3 remaining true limitations")

# Figure (top)
add_figure(slide, "figure16_limitations_addressed.png",
           Inches(0.3), Inches(1.4), Inches(12.7), Inches(3.0))

# Bottom panel: 3-column summary
# Column 1: MLST-addressed
add_box(slide, Inches(0.3), Inches(4.6), Inches(4.0), Inches(0.45), GREEN)
add_text(slide, Inches(0.4), Inches(4.65), Inches(3.8), Inches(0.35),
         "3 Addressed by MLST", font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_box(slide, Inches(0.3), Inches(5.05), Inches(4.0), Inches(2.0), WHITE,
        border_color=GREEN, border_width=Pt(1))
mlst_items = [
    "L1: Host range estimation",
    "L2: Phenotype-genotype validation",
    "L9: Transmission mode inference",
]
add_multiline(slide, Inches(0.5), Inches(5.15), Inches(3.6), Inches(1.8),
              mlst_items, font_size=11, color=DARK_GRAY, spacing=4.0)

# Column 2: New modules
add_box(slide, Inches(4.6), Inches(4.6), Inches(4.0), Inches(0.45), MED_BLUE)
add_text(slide, Inches(4.7), Inches(4.65), Inches(3.8), Inches(0.35),
         "7 Addressed by New Modules", font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_box(slide, Inches(4.6), Inches(5.05), Inches(4.0), Inches(2.0), WHITE,
        border_color=MED_BLUE, border_width=Pt(1))
new_module_items = [
    "L3: Assembly completeness",
    "L4: Database coverage & novelty",
    "L5: Recombination detection",
    "L6: Novel Inc group discovery",
    "L7: Evolutionary rate estimation",
    "L8: Cluster stability assessment",
    "L10: MGE boundary detection",
]
add_multiline(slide, Inches(4.8), Inches(5.15), Inches(3.6), Inches(1.8),
              new_module_items, font_size=10, color=DARK_GRAY, spacing=2.0)

# Column 3: Remaining true limitations
add_box(slide, Inches(8.9), Inches(4.6), Inches(4.1), Inches(0.45), RED)
add_text(slide, Inches(9.0), Inches(4.65), Inches(3.9), Inches(0.35),
         "3 Remaining True Limitations", font_size=13, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
add_box(slide, Inches(8.9), Inches(5.05), Inches(4.1), Inches(2.0), WHITE,
        border_color=RED, border_width=Pt(1))
remaining_items = [
    "Single-linkage chaining at coarse",
    "  thresholds (inherent to method)",
    "Static training set requires",
    "  periodic retraining",
    "Chromosomal sequence sensitivity",
    "  (>500 kb contigs)",
]
add_multiline(slide, Inches(9.1), Inches(5.15), Inches(3.7), Inches(1.8),
              remaining_items, font_size=10, color=DARK_GRAY, spacing=2.0)

# Bottom summary bar
add_box(slide, Inches(0.3), Inches(7.15), Inches(12.7), Inches(0.25), DARK_BLUE)
add_text(slide, Inches(0.5), Inches(7.15), Inches(12.3), Inches(0.25),
         "10/10 limitations addressed (3 MLST + 7 new modules) | 3 inherent constraints remain as true limitations",
         font_size=11, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 31: Conclusions
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Conclusions",
               "Five key takeaways from the pLIN classification system")

conclusions = [
    ("1", "pLIN resolves 3,073 unique codes (Simpson's D = 0.985) \u2014 1.54\u00d7 improvement "
          "over Inc typing alone, providing strain-level resolution across 28 Inc/Rep groups.",
     MED_BLUE),
    ("2", "64,891 AMR detections mapped to specific lineages, including 1,635 carbapenemases "
          "and 204 mcr (colistin resistance) genes linked to transmissible plasmid backbones.",
     RED),
    ("3", "High-risk lineages identified: pLIN 671 (100% KPC-2, IncN, n=90) and pLIN 860 "
          "(44.4% mcr, 5 Inc groups, n=142) represent priority targets for surveillance.",
     ORANGE),
    ("4", "Cross-validated against 74 outbreak plasmids from 26 studies across 13 countries. "
          "3 globally disseminated lineages and 9 intra-study clusters detected.",
     GREEN),
    ("5", "Scalable to 79,305 plasmids (57,886 codes) on standard hardware in <30 minutes. "
          "97.3% classification rate with 2.3 GB peak memory.",
     PURPLE),
]

for i, (num, text, color) in enumerate(conclusions):
    y = Inches(1.5 + i * 1.05)
    # Number circle
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.1), Inches(0.5), Inches(0.5))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(18)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

    add_box(slide, Inches(1.2), y, Inches(11.7), Inches(0.9), WHITE,
            border_color=color, border_width=Pt(2))
    # Colored left accent
    add_box(slide, Inches(1.2), y, Inches(0.1), Inches(0.9), color)
    add_text(slide, Inches(1.5), y + Inches(0.1), Inches(11.2), Inches(0.7),
             text, font_size=13, color=DARK_GRAY)

# Bottom: Open source
add_box(slide, Inches(0.5), Inches(6.8), Inches(12.4), Inches(0.5), DARK_BLUE)
add_text(slide, Inches(0.7), Inches(6.85), Inches(12), Inches(0.4),
         "Open source: github.com/xavierbasilbritto-hub/pLIN-plasmid-classification (GPL-3.0)",
         font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 23: Chromosomal-Plasmid Integration
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Combined Chromosomal-Plasmid Typing",
               "Integrating MLST sequence types with pLIN codes for transmission mode inference")

# Left side: Four transmission modes
add_text(slide, Inches(0.5), Inches(1.5), Inches(6.0), Inches(0.4),
         "Transmission Mode Inference", font_size=18, bold=True, color=DARK_BLUE)

transmission_modes = [
    ("Clonal spread",
     "Same MLST ST + same pLIN L6",
     "Vertical co-transmission of chromosome and plasmid",
     RED),
    ("Horizontal plasmid transfer",
     "Different STs + same pLIN L6",
     "Conjugation / HGT disseminating the plasmid across lineages",
     ORANGE),
    ("Same strain, different plasmids",
     "Same ST + different pLIN L6",
     "Plasmid replacement or acquisition within a clonal background",
     MED_BLUE),
    ("Independent",
     "Different STs + different pLIN L6",
     "No epidemiological link detected",
     MED_GRAY),
]

for i, (mode, criteria, interpretation, color) in enumerate(transmission_modes):
    y = Inches(2.0 + i * 1.15)
    # Colored indicator circle
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.1),
                                  Inches(0.35), Inches(0.35))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    # Mode name
    add_text(slide, Inches(1.0), y, Inches(5.5), Inches(0.35),
             mode, font_size=14, bold=True, color=color)
    # Criteria
    add_text(slide, Inches(1.0), y + Inches(0.35), Inches(5.5), Inches(0.3),
             criteria, font_size=11, color=DARK_GRAY, font_name="Consolas")
    # Interpretation
    add_text(slide, Inches(1.0), y + Inches(0.65), Inches(5.5), Inches(0.3),
             interpretation, font_size=10, color=MED_GRAY)

# Right side: Retrospective validation highlights
add_text(slide, Inches(7.0), Inches(1.5), Inches(6.0), Inches(0.4),
         "Retrospective Validation Highlights", font_size=18, bold=True, color=DARK_BLUE)

# Summary stats box
add_box(slide, Inches(7.0), Inches(2.0), Inches(5.8), Inches(1.2),
        RGBColor(0xE3, 0xF2, 0xFD), border_color=MED_BLUE, border_width=Pt(2))
add_multiline(slide, Inches(7.2), Inches(2.1), Inches(5.4), Inches(1.0), [
    "13 studies, 7 species, 20 STs",
    "92.3% overall concordance (12/13 studies)",
], font_size=13, bold_first=True, color=DARK_BLUE)

# Key validation cases
validation_cases = [
    ("Conlan 2014", "ST258 + pLIN 672", "Clonal spread", RED),
    ("Yao 2023", "ST11 / ST131", "Horizontal gene transfer", ORANGE),
    ("Weber 2019", "4 species, 9 STs", "HGT across species barrier", PURPLE),
    ("Jousset 2019", "OXA-48 cross-species", "HGT confirmed", TEAL),
]

for i, (study, detail, mode, color) in enumerate(validation_cases):
    y = Inches(3.5 + i * 0.85)
    add_box(slide, Inches(7.0), y, Inches(5.8), Inches(0.7), WHITE,
            border_color=color, border_width=Pt(2))
    # Colored left accent
    add_box(slide, Inches(7.0), y, Inches(0.08), Inches(0.7), color)
    add_text(slide, Inches(7.2), y + Inches(0.05), Inches(2.0), Inches(0.3),
             study, font_size=12, bold=True, color=color)
    add_text(slide, Inches(9.2), y + Inches(0.05), Inches(2.0), Inches(0.3),
             detail, font_size=11, color=DARK_GRAY, font_name="Consolas")
    add_text(slide, Inches(11.2), y + Inches(0.05), Inches(1.5), Inches(0.3),
             mode, font_size=10, bold=True, color=color, alignment=PP_ALIGN.RIGHT)

# Bottom summary
add_box(slide, Inches(0.5), Inches(6.3), Inches(12.4), Inches(0.9),
        RGBColor(0xE8, 0xF5, 0xE9), border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.45), Inches(12), Inches(0.5),
         "Combining MLST ST + pLIN L6 enables mechanistic inference of AMR transmission pathways "
         "without requiring WGS phylogenetics",
         font_size=13, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 24: Method Stability & Universality
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "pLIN: Method Stability and Universality",
               "Design principles ensuring long-term reliability and broad applicability")

# Three key properties
add_text(slide, Inches(0.5), Inches(1.5), Inches(12.4), Inches(0.4),
         "Three Key Design Properties", font_size=20, bold=True, color=DARK_BLUE)

properties = [
    ("Code Permanence",
     "Nearest-neighbour assignment ensures that once a pLIN code is issued, it never changes. "
     "New sequences are assigned to the closest existing cluster or receive a new code. "
     "Previously assigned codes remain stable regardless of database growth.",
     MED_BLUE,
     "Codes never change"),
    ("Metric Universality",
     "Tetranucleotide (4-mer) frequency is an intrinsic property of any DNA sequence. "
     "It requires no gene annotation, no reference alignment, and no species-specific databases. "
     "The cosine distance metric operates identically on any input sequence.",
     GREEN,
     "Intrinsic sequence property"),
    ("Taxonomic Generality",
     "The pLIN framework is applicable to any DNA replicon: plasmids, phages, chromosomes, "
     "or metagenomic contigs. The 6-level hierarchical structure adapts to any genomic entity "
     "where sequence composition reflects evolutionary relationships.",
     PURPLE,
     "Any DNA sequence"),
]

for i, (title, desc, color, tag) in enumerate(properties):
    x = Inches(0.4 + i * 4.25)
    y = Inches(2.0)
    # Title bar
    add_box(slide, x, y, Inches(4.0), Inches(0.55), color)
    add_text(slide, x + Inches(0.1), y + Inches(0.07), Inches(3.8), Inches(0.4),
             title, font_size=15, bold=True, color=WHITE)
    # Tag subtitle
    add_box(slide, x, y + Inches(0.55), Inches(4.0), Inches(0.35), RGBColor(0xF5, 0xF5, 0xF5))
    add_text(slide, x + Inches(0.1), y + Inches(0.58), Inches(3.8), Inches(0.3),
             tag, font_size=11, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    # Description box
    add_box(slide, x, y + Inches(0.9), Inches(4.0), Inches(2.2), WHITE,
            border_color=color, border_width=Pt(1))
    add_text(slide, x + Inches(0.15), y + Inches(1.0), Inches(3.7), Inches(2.0),
             desc, font_size=11, color=DARK_GRAY)

# Vulnerability & Mitigation section
add_text(slide, Inches(0.5), Inches(5.3), Inches(12.4), Inches(0.4),
         "Known Vulnerability & Mitigation", font_size=18, bold=True, color=ORANGE)

# Vulnerability box (left)
add_box(slide, Inches(0.5), Inches(5.8), Inches(5.8), Inches(1.3),
        RGBColor(0xFF, 0xF3, 0xE0), border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(5.9), Inches(5.4), Inches(0.3),
         "Vulnerability: Threshold Calibration", font_size=14, bold=True, color=ORANGE)
add_text(slide, Inches(0.7), Inches(6.25), Inches(5.4), Inches(0.7),
         "Current pLIN thresholds (L1-L6) were calibrated on Enterobacterales plasmids. "
         "Applying these thresholds to phylogenetically distant replicons (e.g., Gram-positive, "
         "Acinetobacter) may produce suboptimal clustering resolution.",
         font_size=11, color=DARK_GRAY)

# Mitigation box (right)
add_box(slide, Inches(6.8), Inches(5.8), Inches(5.8), Inches(1.3),
        RGBColor(0xE8, 0xF5, 0xE9), border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.0), Inches(5.9), Inches(5.4), Inches(0.3),
         "Mitigation: Modular Architecture", font_size=14, bold=True, color=GREEN)
add_text(slide, Inches(7.0), Inches(6.25), Inches(5.4), Inches(0.7),
         "pLIN's modular design enables recalibration of distance thresholds per taxonomic group "
         "without disrupting previously assigned codes. New threshold sets can be deployed as "
         "configuration profiles while preserving backward compatibility.",
         font_size=11, color=DARK_GRAY)


# ======================================================================
# SLIDE 25: Acknowledgments & Contact
# ======================================================================

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, DARK_BLUE)

add_text(slide, Inches(0.8), Inches(0.8), Inches(11.5), Inches(0.8),
         "Acknowledgments & Contact",
         font_size=36, bold=True, color=WHITE)

# Authors with affiliations
add_box(slide, Inches(0.8), Inches(1.8), Inches(11.5), Inches(2.2), RGBColor(0x0A, 0x2F, 0x6E),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(1.0), Inches(1.9), Inches(11.1), Inches(0.4),
         "Authors", font_size=18, bold=True, color=LIGHT_BLUE)

authors = [
    "Basil Britto Xavier \u2014 University Medical Center Groningen (UMCG), Department of Medical Microbiology and Infection Prevention",
    "Anurag Kumar Bari \u2014 University Medical Center Groningen (UMCG), Department of Medical Microbiology and Infection Prevention",
    "Bhanu Sinha \u2014 University Medical Center Groningen (UMCG), Department of Medical Microbiology and Infection Prevention",
    "John W A Rossen \u2014 University Medical Center Groningen (UMCG), Department of Medical Microbiology and Infection Prevention",
]
add_multiline(slide, Inches(1.0), Inches(2.4), Inches(11.1), Inches(1.5),
              authors, font_size=12, color=RGBColor(0x90, 0xCA, 0xF9), spacing=3.0)

# GitHub
add_box(slide, Inches(0.8), Inches(4.3), Inches(5.4), Inches(1.2), RGBColor(0x0A, 0x2F, 0x6E),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(1.0), Inches(4.4), Inches(5.0), Inches(0.4),
         "GitHub Repository", font_size=16, bold=True, color=LIGHT_BLUE)
add_text(slide, Inches(1.0), Inches(4.8), Inches(5.0), Inches(0.5),
         "github.com/xavierbasilbritto-hub/pLIN-plasmid-classification",
         font_size=13, color=WHITE, font_name="Consolas")

# License
add_box(slide, Inches(6.5), Inches(4.3), Inches(5.8), Inches(1.2), RGBColor(0x0A, 0x2F, 0x6E),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(6.7), Inches(4.4), Inches(5.4), Inches(0.4),
         "License", font_size=16, bold=True, color=LIGHT_BLUE)
add_text(slide, Inches(6.7), Inches(4.8), Inches(5.4), Inches(0.5),
         "GPL-3.0 + Citation clause (see CITATION.cff)",
         font_size=13, color=WHITE)

# Deployment
add_box(slide, Inches(0.8), Inches(5.8), Inches(11.5), Inches(1.0), RGBColor(0x0A, 0x2F, 0x6E),
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(1.0), Inches(5.9), Inches(11.1), Inches(0.4),
         "Availability", font_size=16, bold=True, color=LIGHT_BLUE)
add_text(slide, Inches(1.0), Inches(6.3), Inches(11.1), Inches(0.4),
         "Available as Docker container, standalone desktop app, or cloud deployment via Streamlit Community Cloud",
         font_size=14, color=WHITE)

add_text(slide, Inches(0.8), Inches(7.0), Inches(11.5), Inches(0.3),
         "pLIN \u2014 Plasmid Lineage Identification Number \u2014 Manuscript Presentation",
         font_size=11, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ======================================================================
# SLIDE 9b: Classification System (Figure 9) — BONUS embedded as supplemental
# This is inserted to embed Figure9_classification_system.png
# ======================================================================

# Note: The 20 slides above cover the full plan. Figures 9 is referenced
# contextually in the pipeline and classification slides already.


# ======================================================================
# Save
# ======================================================================

out_path = os.path.join(BASE_DIR, "output", "pLIN_Manuscript_Presentation.pptx")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
prs.save(out_path)
print(f"Saved: {out_path}")
print(f"Slides: {len(prs.slides)}")
