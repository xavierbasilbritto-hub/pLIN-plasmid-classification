#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""Generate PowerPoint presentation of pLIN tool architecture."""

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


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 1: Title Slide
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])  # blank
add_bg(slide, DARK_BLUE)

add_text(slide, Inches(0.8), Inches(1.5), Inches(11.5), Inches(1),
         "pLIN — Plasmid Life Identification Number",
         font_size=40, bold=True, color=WHITE)

add_text(slide, Inches(0.8), Inches(2.6), Inches(11.5), Inches(0.8),
         "A Hierarchical, Reference-Free Classification System\nfor Bacterial Plasmid Genomes with AMR Surveillance",
         font_size=22, color=LIGHT_BLUE)

# Feature boxes
features = [
    ("6,998", "Training\nPlasmids"),
    ("20", "Inc Groups\n(20 families)"),
    ("2,232", "Unique pLIN\nCodes"),
    ("92.2%", "Inc Detection\nAccuracy"),
    ("27,465", "AMR Gene\nDetections"),
]
for i, (num, label) in enumerate(features):
    x = Inches(0.8 + i * 2.5)
    y = Inches(4.0)
    box = add_box(slide, x, y, Inches(2.2), Inches(1.8), RGBColor(0x0A, 0x2F, 0x6E),
                  border_color=MED_BLUE, border_width=Pt(2))
    add_text(slide, x, y + Inches(0.2), Inches(2.2), Inches(0.7),
             num, font_size=36, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, y + Inches(0.9), Inches(2.2), Inches(0.7),
             label, font_size=14, color=LIGHT_BLUE, alignment=PP_ALIGN.CENTER)

add_text(slide, Inches(0.8), Inches(6.5), Inches(11.5), Inches(0.5),
         "Tool Architecture Overview", font_size=14, color=MED_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: System Overview / Architecture Diagram
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "System Architecture", "End-to-end pipeline from FASTA input to classified plasmids + AMR profile")

# Input
add_box(slide, Inches(0.4), Inches(1.6), Inches(2.4), Inches(1.4), LIGHT_BLUE, border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.5), Inches(1.7), Inches(2.2), Inches(0.4),
         "INPUT", font_size=14, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(0.5), Inches(2.1), Inches(2.2), Inches(0.8), [
    "Plasmid FASTA files",
    "(.fasta, .fa, .fna)",
], font_size=12, color=DARK_GRAY)

# Arrow 1
add_right_arrow(slide, Inches(2.95), Inches(2.0), Inches(0.6), Inches(0.4))

# Inc Detection
add_box(slide, Inches(3.7), Inches(1.6), Inches(2.6), Inches(1.4), RGBColor(0xFF, 0xF3, 0xE0), border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(3.8), Inches(1.7), Inches(2.4), Inches(0.4),
         "INC GROUP DETECTION", font_size=12, bold=True, color=ORANGE, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(3.8), Inches(2.1), Inches(2.4), Inches(0.8), [
    "KNN Classifier (k=5)",
    "Cosine distance on 4-mers",
    "Trained on 6,998 plasmids",
], font_size=11, color=DARK_GRAY)

# Arrow 2
add_right_arrow(slide, Inches(6.45), Inches(2.0), Inches(0.6), Inches(0.4))

# K-mer + Clustering
add_box(slide, Inches(7.2), Inches(1.6), Inches(2.8), Inches(1.4), RGBColor(0xE8, 0xF5, 0xE9), border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.3), Inches(1.7), Inches(2.6), Inches(0.4),
         "pLIN ASSIGNMENT", font_size=12, bold=True, color=GREEN, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(7.3), Inches(2.1), Inches(2.6), Inches(0.8), [
    "4-mer frequency vectors (256D)",
    "Cosine distance matrix",
    "Single-linkage clustering",
    "6 threshold cuts \u2192 A.B.C.D.E.F",
], font_size=11, color=DARK_GRAY)

# Arrow 3
add_right_arrow(slide, Inches(10.15), Inches(2.0), Inches(0.6), Inches(0.4))

# Output
add_box(slide, Inches(10.9), Inches(1.6), Inches(2.1), Inches(1.4), RGBColor(0xE3, 0xF2, 0xFD), border_color=DARK_BLUE, border_width=Pt(2))
add_text(slide, Inches(11.0), Inches(1.7), Inches(1.9), Inches(0.4),
         "pLIN CODES", font_size=12, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(11.0), Inches(2.1), Inches(1.9), Inches(0.8), [
    "Hierarchical codes",
    "e.g. 1.1.1.1.2.5",
    "Permanent & stable",
], font_size=11, color=DARK_GRAY)

# AMR Branch (parallel path)
add_box(slide, Inches(0.4), Inches(3.5), Inches(2.4), Inches(1.4), RGBColor(0xFC, 0xE4, 0xEC), border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.5), Inches(3.6), Inches(2.2), Inches(0.4),
         "AMRFinderPlus", font_size=14, bold=True, color=RED, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(0.5), Inches(4.0), Inches(2.2), Inches(0.8), [
    "NCBI AMR detection",
    "AMR + Stress + Virulence",
    "Per-contig gene calls",
], font_size=11, color=DARK_GRAY)

# Arrow from AMR
add_right_arrow(slide, Inches(2.95), Inches(3.9), Inches(0.6), Inches(0.4))

# Integration
add_box(slide, Inches(3.7), Inches(3.5), Inches(2.6), Inches(1.4), RGBColor(0xF3, 0xE5, 0xF5), border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(3.8), Inches(3.6), Inches(2.4), Inches(0.4),
         "INTEGRATION", font_size=12, bold=True, color=PURPLE, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(3.8), Inches(4.0), Inches(2.4), Inches(0.8), [
    "pLIN codes + AMR genes",
    "Per-plasmid AMR profile",
    "Lineage AMR summary",
], font_size=11, color=DARK_GRAY)

# Arrow to Visualization
add_right_arrow(slide, Inches(6.45), Inches(3.9), Inches(0.6), Inches(0.4))

# Visualization
add_box(slide, Inches(7.2), Inches(3.5), Inches(2.8), Inches(1.4), RGBColor(0xE0, 0xF7, 0xFA), border_color=TEAL, border_width=Pt(2))
add_text(slide, Inches(7.3), Inches(3.6), Inches(2.6), Inches(0.4),
         "VISUALIZATION", font_size=12, bold=True, color=TEAL, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(7.3), Inches(4.0), Inches(2.6), Inches(0.8), [
    "Cladograms (4 types)",
    "AMR heatmaps (Plotly)",
    "Drug class pie charts",
    "Critical gene alerts",
], font_size=11, color=DARK_GRAY)

# Arrow to Export
add_right_arrow(slide, Inches(10.15), Inches(3.9), Inches(0.6), Inches(0.4))

# Export
add_box(slide, Inches(10.9), Inches(3.5), Inches(2.1), Inches(1.4), LIGHT_GRAY, border_color=MED_GRAY, border_width=Pt(2))
add_text(slide, Inches(11.0), Inches(3.6), Inches(1.9), Inches(0.4),
         "EXPORT", font_size=12, bold=True, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)
add_multiline(slide, Inches(11.0), Inches(4.0), Inches(1.9), Inches(0.8), [
    "TSV tables",
    "PNG / PDF figures",
    "ZIP bundle",
], font_size=11, color=DARK_GRAY)

# GUI bar at bottom
add_box(slide, Inches(0.4), Inches(5.5), Inches(12.6), Inches(1.5), RGBColor(0xE8, 0xEA, 0xF6), border_color=RGBColor(0x53, 0x4B, 0xAE), border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(5.6), Inches(3), Inches(0.4),
         "Streamlit Web GUI (plin_app.py)", font_size=14, bold=True, color=RGBColor(0x53, 0x4B, 0xAE))
tabs = ["Upload & Config", "Overview Tab", "Results Tab", "Cladogram Tab", "AMR Analysis Tab", "Export Tab"]
for i, tab in enumerate(tabs):
    x = Inches(0.6 + i * 2.05)
    add_box(slide, x, Inches(6.1), Inches(1.9), Inches(0.6), WHITE, border_color=RGBColor(0x53, 0x4B, 0xAE))
    add_text(slide, x, Inches(6.15), Inches(1.9), Inches(0.5),
             tab, font_size=11, bold=True, color=RGBColor(0x53, 0x4B, 0xAE), alignment=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: pLIN Classification System
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "pLIN Classification System", "Six-level hierarchical coding from tetranucleotide composition")

# Threshold table
levels = [
    ("A", "L1 (Family)", "0.150", "~85%", "#E53935"),
    ("B", "L2 (Subfamily)", "0.100", "~90%", "#FB8C00"),
    ("C", "L3 (Cluster)", "0.050", "~95%", "#FDD835"),
    ("D", "L4 (Subcluster)", "0.020", "~98%", "#43A047"),
    ("E", "L5 (Clone)", "0.010", "~99%", "#1E88E5"),
    ("F", "L6 (Strain)", "0.001", "~99.9%", "#8E24AA"),
]

# Table header
add_box(slide, Inches(0.5), Inches(1.6), Inches(5.5), Inches(0.5), DARK_BLUE)
headers = ["Position", "Level", "Threshold (d\u2264)", "ANI Equiv."]
for i, h in enumerate(headers):
    add_text(slide, Inches(0.6 + i * 1.35), Inches(1.63), Inches(1.3), Inches(0.4),
             h, font_size=12, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

for j, (pos, level, thresh, ani, hex_col) in enumerate(levels):
    y = Inches(2.1 + j * 0.45)
    bg = LIGHT_GRAY if j % 2 == 0 else WHITE
    add_box(slide, Inches(0.5), y, Inches(5.5), Inches(0.45), bg)
    r, g, b = int(hex_col[1:3], 16), int(hex_col[3:5], 16), int(hex_col[5:7], 16)
    # Color dot
    dot = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.7), y + Inches(0.1), Inches(0.25), Inches(0.25))
    dot.fill.solid()
    dot.fill.fore_color.rgb = RGBColor(r, g, b)
    dot.line.fill.background()
    add_text(slide, Inches(1.0), y + Inches(0.05), Inches(0.9), Inches(0.35),
             pos, font_size=13, bold=True, color=RGBColor(r, g, b), alignment=PP_ALIGN.CENTER)
    add_text(slide, Inches(1.95), y + Inches(0.05), Inches(1.3), Inches(0.35),
             level, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)
    add_text(slide, Inches(3.3), y + Inches(0.05), Inches(1.3), Inches(0.35),
             f"d \u2264 {thresh}", font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER, font_name="Consolas")
    add_text(slide, Inches(4.65), y + Inches(0.05), Inches(1.3), Inches(0.35),
             ani, font_size=12, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

# Example pLIN code
add_box(slide, Inches(0.5), Inches(5.0), Inches(5.5), Inches(1.0), RGBColor(0xE3, 0xF2, 0xFD), border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(5.1), Inches(5), Inches(0.4),
         "Example pLIN Code:", font_size=13, bold=True, color=DARK_BLUE)
add_text(slide, Inches(0.7), Inches(5.5), Inches(5), Inches(0.4),
         "1 . 1 . 1 . 1 . 2 . 5", font_size=24, bold=True, color=DARK_BLUE, alignment=PP_ALIGN.CENTER, font_name="Consolas")

# Right side: Algorithm steps
add_text(slide, Inches(6.8), Inches(1.6), Inches(6), Inches(0.5),
         "Classification Algorithm", font_size=20, bold=True, color=DARK_BLUE)

steps = [
    ("1", "Parse FASTA", "Read plasmid sequences using BioPython", GREEN),
    ("2", "Compute 4-mer Vectors", "256-dimensional normalised frequency vectors", ORANGE),
    ("3", "Cosine Distance Matrix", "Pairwise distance between all sequences", MED_BLUE),
    ("4", "Single-Linkage Clustering", "Hierarchical agglomerative clustering (scipy)", PURPLE),
    ("5", "Threshold Cutting", "Cut dendrogram at 6 thresholds \u2192 6 cluster IDs", RED),
    ("6", "Assign pLIN Codes", "Concatenate cluster IDs: A.B.C.D.E.F", DARK_BLUE),
]

for i, (num, title, desc, color) in enumerate(steps):
    y = Inches(2.2 + i * 0.72)
    # Number circle
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(6.8), y, Inches(0.4), Inches(0.4))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(14)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

    add_text(slide, Inches(7.4), y - Inches(0.02), Inches(5.5), Inches(0.3),
             title, font_size=14, bold=True, color=color)
    add_text(slide, Inches(7.4), y + Inches(0.28), Inches(5.5), Inches(0.3),
             desc, font_size=11, color=MED_GRAY)

    if i < len(steps) - 1:
        # Connecting line
        line = slide.shapes.add_connector(1, Inches(7.0), y + Inches(0.4), Inches(7.0), y + Inches(0.72))
        line.line.color.rgb = LIGHT_BLUE
        line.line.width = Pt(2)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: Inc Group Auto-Detection
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Inc Group Auto-Detection", "KNN classifier trained on 6,998 RefSeq plasmids across 20 Inc groups")

# Training data
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.5),
         "Training Dataset (Top 6 of 20 Inc Groups)", font_size=18, bold=True, color=DARK_BLUE)

groups = [
    ("IncFII", "4,629", "66.1%", RED),
    ("IncN", "1,064", "15.2%", GREEN),
    ("IncX1", "701", "10.0%", MED_BLUE),
    ("IncF", "148", "2.1%", ORANGE),
    ("IncI1", "72", "1.0%", PURPLE),
    ("Others (14)", "384", "5.5%", MED_GRAY),
]

for i, (name, count, pct, color) in enumerate(groups):
    y = Inches(2.1 + i * 0.7)
    bar_width = max(0.5, float(pct.replace("%", "")) / 100 * 5.5)
    add_box(slide, Inches(0.5), y, Inches(bar_width), Inches(0.5), color)
    add_text(slide, Inches(0.7), y + Inches(0.07), Inches(4), Inches(0.35),
             f"{name}  —  {count} ({pct})", font_size=12, bold=True, color=WHITE)

# Classifier details
add_text(slide, Inches(0.5), Inches(6.4), Inches(6), Inches(0.5),
         "Classifier Details", font_size=16, bold=True, color=DARK_BLUE)
details = [
    "Algorithm: K-Nearest Neighbors (k=5, cosine, distance-weighted)",
    "5-fold stratified cross-validation accuracy: 92.2%",
    "Dynamic k: adapts to smallest class size",
    "Model file: data/inc_classifier.npz",
]
add_multiline(slide, Inches(0.5), Inches(6.8), Inches(6), Inches(0.8), details, font_size=11, color=DARK_GRAY)

# Right side: How it works
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.5),
         "How Auto-Detection Works", font_size=20, bold=True, color=DARK_BLUE)

steps_inc = [
    "Upload new plasmid FASTA",
    "Compute 4-mer frequency vector (256 features)",
    "Find 5 nearest neighbors in training set (cosine distance)",
    "Distance-weighted vote \u2192 Inc group prediction",
    "Return predicted group + confidence probability",
]
for i, step in enumerate(steps_inc):
    y = Inches(2.3 + i * 0.65)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(7.0), y, Inches(0.35), Inches(0.35))
    circ.fill.solid()
    circ.fill.fore_color.rgb = ORANGE
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = str(i + 1)
    tf.paragraphs[0].font.size = Pt(12)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

    add_text(slide, Inches(7.5), y, Inches(5.5), Inches(0.35),
             step, font_size=13, color=DARK_GRAY)

# Output example
add_box(slide, Inches(7.0), Inches(5.7), Inches(5.8), Inches(1.3), RGBColor(0xFf, 0xF8, 0xE1), border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(7.2), Inches(5.8), Inches(5.4), Inches(0.3),
         "Example Output:", font_size=12, bold=True, color=ORANGE)
add_multiline(slide, Inches(7.2), Inches(6.1), Inches(5.4), Inches(0.8), [
    "Plasmid: SP12_P2  \u2192  IncX1 (confidence: 0.79)",
    "  Top 5: IncX1: 0.79 | IncN: 0.12 | IncX3: 0.04 | IncFII: 0.03 | IncI1: 0.02",
], font_size=11, color=DARK_GRAY, font_name="Consolas")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: AMRFinderPlus Integration
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "AMRFinderPlus Integration", "NCBI antimicrobial resistance gene detection pipeline")

# Left: AMR pipeline
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.5),
         "AMR Detection Pipeline", font_size=20, bold=True, color=RED)

amr_steps = [
    ("Auto-Detect", "Search PATH, conda envs for amrfinder binary", RGBColor(0x78, 0x90, 0x9C)),
    ("Run per File", "amrfinder -n {fasta} -d {db} --plus -o {output}", RED),
    ("Parse Results", "Read TSV output: gene, class, type, coordinates", ORANGE),
    ("Integrate", "Merge AMR hits with pLIN assignments per plasmid", PURPLE),
    ("Visualize", "Heatmaps, bar charts, pie charts, alerts", TEAL),
]

for i, (title, desc, color) in enumerate(amr_steps):
    y = Inches(2.3 + i * 0.85)
    add_box(slide, Inches(0.5), y, Inches(5.8), Inches(0.7), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, Inches(0.7), y + Inches(0.05), Inches(1.8), Inches(0.3),
             title, font_size=13, bold=True, color=color)
    add_text(slide, Inches(0.7), y + Inches(0.35), Inches(5.4), Inches(0.3),
             desc, font_size=11, color=MED_GRAY)

# Right: Detection types
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.5),
         "Detection Categories", font_size=20, bold=True, color=RED)

categories = [
    ("AMR Genes", "Antibiotic resistance genes\nblaNDM, blaKPC, blaSHV, mcr, etc.", RED, "27,465"),
    ("Stress Genes", "Biocide/metal resistance\nbleomycin, mercury, arsenic, etc.", ORANGE, "5,834"),
    ("Virulence", "Virulence-associated factors\nType III secretion, toxins, etc.", PURPLE, "1,200+"),
]

for i, (title, desc, color, count) in enumerate(categories):
    y = Inches(2.3 + i * 1.4)
    add_box(slide, Inches(7.0), y, Inches(5.8), Inches(1.2), WHITE, border_color=color, border_width=Pt(2))
    dot = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(7.3), y + Inches(0.35), Inches(0.4), Inches(0.4))
    dot.fill.solid()
    dot.fill.fore_color.rgb = color
    dot.line.fill.background()
    add_text(slide, Inches(7.9), y + Inches(0.1), Inches(3), Inches(0.35),
             title, font_size=15, bold=True, color=color)
    add_text(slide, Inches(7.9), y + Inches(0.5), Inches(4.5), Inches(0.6),
             desc, font_size=11, color=MED_GRAY)
    add_text(slide, Inches(11.0), y + Inches(0.2), Inches(1.5), Inches(0.8),
             count, font_size=22, bold=True, color=color, alignment=PP_ALIGN.RIGHT)

# Critical gene alert box
add_box(slide, Inches(7.0), Inches(6.0), Inches(5.8), Inches(1.0), RGBColor(0xFF, 0xEB, 0xEE), border_color=RED, border_width=Pt(2))
add_text(slide, Inches(7.2), Inches(6.1), Inches(5.4), Inches(0.3),
         "\u26a0\ufe0f  Critical Gene Alert System", font_size=13, bold=True, color=RED)
add_text(slide, Inches(7.2), Inches(6.45), Inches(5.4), Inches(0.4),
         "Carbapenemases (KPC, NDM, OXA-48, VIM, IMP)  |  ESBLs (CTX-M, SHV-12, TEM)  |  Colistin (mcr)",
         font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 6: Visualization & GUI
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Streamlit Web GUI", "Interactive web application for plasmid analysis (plin_app.py)")

# GUI Layout diagram
add_text(slide, Inches(0.5), Inches(1.6), Inches(12), Inches(0.5),
         "Application Layout", font_size=20, bold=True, color=DARK_BLUE)

# Upload area
add_box(slide, Inches(0.5), Inches(2.2), Inches(8.5), Inches(1.2), LIGHT_BLUE, border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(2.3), Inches(4), Inches(0.4),
         "\u2191 Upload FASTA Files", font_size=14, bold=True, color=DARK_BLUE)
add_text(slide, Inches(5.5), Inches(2.3), Inches(3), Inches(0.4),
         "Inc Group: [Auto-detect \u25bc]", font_size=12, color=DARK_GRAY)
add_text(slide, Inches(5.5), Inches(2.7), Inches(3), Inches(0.4),
         "\u2611 Run AMRFinderPlus", font_size=12, color=DARK_GRAY)
add_box(slide, Inches(0.7), Inches(2.9), Inches(2), Inches(0.4), MED_BLUE)
add_text(slide, Inches(0.7), Inches(2.92), Inches(2), Inches(0.35),
         "\u25b6 Run Analysis", font_size=12, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

# Tabs
tabs_info = [
    ("Overview", "pLIN description\nThreshold table\nMetrics dashboard", MED_BLUE),
    ("Results", "Interactive data table\nSearch/filter\npLIN distribution", GREEN),
    ("Cladogram", "Rectangular\nCircular\nHeatmap\nAMR Annotated", ORANGE),
    ("AMR", "Gene prevalence\nDrug class pie\nCritical alerts", RED),
    ("Epidemiology", "Mobility prediction\nOutbreak detection\nRisk stratification", TEAL),
    ("CRISPR Host", "Spacer extraction\nHost-plasmid heatmap\nProbability ranking", RGBColor(0x00, 0x69, 0x5C)),
    ("Buddy", "AI chatbot\nContext-aware Q&A\nOllama LLM", RGBColor(0x6A, 0x1B, 0x9A)),
    ("Export", "TSV tables\nPNG/PDF figures\nZIP bundle", PURPLE),
]

for i, (name, desc, color) in enumerate(tabs_info):
    x = Inches(0.4 + i * 1.6)
    y = Inches(3.8)
    add_box(slide, x, y, Inches(1.5), Inches(0.4), color)
    add_text(slide, x, y + Inches(0.03), Inches(1.5), Inches(0.35),
             name, font_size=9, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, y + Inches(0.4), Inches(1.5), Inches(2.5), WHITE, border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.05), y + Inches(0.5), Inches(1.4), Inches(2.3),
                  desc.split("\n"), font_size=8, color=DARK_GRAY)

# Sidebar
add_box(slide, Inches(9.3), Inches(2.2), Inches(3.7), Inches(4.8), LIGHT_GRAY, border_color=MED_GRAY, border_width=Pt(1))
add_text(slide, Inches(9.5), Inches(2.3), Inches(3.3), Inches(0.4),
         "Sidebar", font_size=16, bold=True, color=DARK_GRAY)
sidebar_items = [
    "Analysis status & metrics",
    "Plasmid count",
    "Unique pLIN codes",
    "Detected Inc groups",
    "",
    "\U0001f504 Clear & Reset",
    "\U0001f4c2 New Analysis",
]
add_multiline(slide, Inches(9.5), Inches(2.8), Inches(3.3), Inches(4), sidebar_items, font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 7: File Structure & Project Organization
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Project Structure", "File organization and data flow")

# Left column: Scripts
add_text(slide, Inches(0.5), Inches(1.6), Inches(4), Inches(0.5),
         "Python Scripts", font_size=18, bold=True, color=DARK_BLUE)

scripts = [
    ("plin_app.py", "Streamlit GUI (main app)", MED_BLUE),
    ("assign_pLIN.py", "Batch pLIN assignment", GREEN),
    ("integrate_pLIN_AMR.py", "pLIN + AMR integration", PURPLE),
    ("generate_figures.py", "Publication figures", ORANGE),
    ("build_inc_centroids.py", "Train Inc classifier", RED),
    ("test_pLIN.py", "Test pLIN on 22 plasmids", MED_GRAY),
    ("test_cladogram.py", "Test cladogram generation", MED_GRAY),
    ("test_integrate_and_cladogram.py", "Test AMR + cladogram", MED_GRAY),
]

for i, (fname, desc, color) in enumerate(scripts):
    y = Inches(2.1 + i * 0.46)
    dot = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.6), y + Inches(0.08), Inches(0.2), Inches(0.2))
    dot.fill.solid()
    dot.fill.fore_color.rgb = color
    dot.line.fill.background()
    add_text(slide, Inches(0.9), y, Inches(2), Inches(0.35),
             fname, font_size=11, bold=True, color=DARK_GRAY, font_name="Consolas")
    add_text(slide, Inches(3.0), y, Inches(2), Inches(0.35),
             desc, font_size=10, color=MED_GRAY)

# Middle column: Data
add_text(slide, Inches(5.2), Inches(1.6), Inches(3.5), Inches(0.5),
         "Data & Training", font_size=18, bold=True, color=DARK_BLUE)

data_items = [
    "plasmid_sequences_for_training/",
    "  \u251c\u2500 20 Inc group folders/fastas/",
    "  \u2514\u2500 6,998 training files total",
    "",
    "reference/",
    "  \u2514\u2500 72,556 individual plasmid FASTAs",
    "",
    "data/",
    "  \u251c\u2500 inc_classifier.npz",
    "  \u2514\u2500 inc_centroids.npz",
    "",
    "test_plasmids/",
    "  \u2514\u2500 IncX/  (22 FASTA files)",
]

add_multiline(slide, Inches(5.2), Inches(2.1), Inches(3.8), Inches(4.5),
              data_items, font_size=11, color=DARK_GRAY, font_name="Consolas")

# Right column: Output
add_text(slide, Inches(9.5), Inches(1.6), Inches(3.5), Inches(0.5),
         "Output", font_size=18, bold=True, color=DARK_BLUE)

output_items = [
    "output/",
    "  \u251c\u2500 pLIN_assignments.tsv",
    "  \u251c\u2500 amrfinder/",
    "  \u2502   \u2514\u2500 amrfinder_all_plasmids.tsv",
    "  \u251c\u2500 integrated/",
    "  \u2502   \u251c\u2500 pLIN_AMR_integrated.tsv",
    "  \u2502   \u2514\u2500 pLIN_lineage_AMR_summary.tsv",
    "  \u251c\u2500 figures/ (7 figs, PNG+PDF)",
    "  \u2514\u2500 test/ (test results)",
]

add_multiline(slide, Inches(9.5), Inches(2.1), Inches(3.5), Inches(4.5),
              output_items, font_size=11, color=DARK_GRAY, font_name="Consolas")

# Bottom: Shell scripts
add_box(slide, Inches(0.5), Inches(6.0), Inches(12.4), Inches(1.0), LIGHT_GRAY, border_color=MED_GRAY)
add_text(slide, Inches(0.7), Inches(6.05), Inches(3), Inches(0.35),
         "Automation Scripts:", font_size=13, bold=True, color=DARK_GRAY)
shell_items = [
    "setup.sh / setup.bat  \u2014  Create venv & install deps",
    "run_all.sh / run_all.bat  \u2014  Full pipeline: pLIN \u2192 AMR \u2192 Integration \u2192 Figures",
    "run_amrfinder_all.sh  \u2014  Batch AMRFinderPlus on all training data",
]
add_multiline(slide, Inches(0.7), Inches(6.35), Inches(11.5), Inches(0.6),
              shell_items, font_size=11, color=DARK_GRAY, font_name="Consolas")


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 8: Technology Stack
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Technology Stack", "Core dependencies and tools")

categories_tech = [
    ("Core Scientific", [
        ("NumPy \u2265 2.0", "Array operations, k-mer vectors"),
        ("Pandas \u2265 2.0", "DataFrames, TSV I/O"),
        ("SciPy \u2265 1.12", "pdist, linkage, fcluster, dendrogram"),
        ("BioPython \u2265 1.80", "FASTA parsing (SeqIO)"),
    ], MED_BLUE),
    ("Machine Learning", [
        ("scikit-learn \u2265 1.4", "KNN classifier for Inc detection"),
        ("XGBoost \u2265 2.0", "Validation classifier (F1=0.903)"),
        ("Optuna \u2265 3.5", "Hyperparameter optimization"),
    ], GREEN),
    ("Visualization", [
        ("Matplotlib \u2265 3.8", "Cladograms, static figures"),
        ("Seaborn \u2265 0.13", "Statistical plots, heatmaps"),
        ("Plotly \u2265 5.18", "Interactive charts (AMR tab)"),
    ], ORANGE),
    ("Web & External", [
        ("Streamlit \u2265 1.31", "Web GUI framework"),
        ("AMRFinderPlus", "NCBI AMR gene detection"),
        ("python-pptx", "PowerPoint generation"),
    ], PURPLE),
]

for i, (cat_name, items, color) in enumerate(categories_tech):
    x = Inches(0.5 + i * 3.2)
    add_box(slide, x, Inches(1.6), Inches(3.0), Inches(0.5), color)
    add_text(slide, x, Inches(1.63), Inches(3.0), Inches(0.4),
             cat_name, font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

    for j, (lib, desc) in enumerate(items):
        y = Inches(2.3 + j * 0.75)
        add_box(slide, x, y, Inches(3.0), Inches(0.65), WHITE, border_color=color, border_width=Pt(1))
        add_text(slide, x + Inches(0.1), y + Inches(0.03), Inches(2.8), Inches(0.3),
                 lib, font_size=12, bold=True, color=color, font_name="Consolas")
        add_text(slide, x + Inches(0.1), y + Inches(0.33), Inches(2.8), Inches(0.3),
                 desc, font_size=10, color=MED_GRAY)

# Key metrics at bottom
add_box(slide, Inches(0.5), Inches(5.6), Inches(12.4), Inches(1.4), RGBColor(0xE8, 0xF5, 0xE9), border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(5.7), Inches(12), Inches(0.4),
         "Key Performance Metrics", font_size=16, bold=True, color=GREEN)

metrics = [
    ("Simpson's D", "0.979", "Discriminatory power of pLIN system"),
    ("CV Accuracy", "92.2%", "Inc group auto-detection (20 groups, 5-fold)"),
    ("XGBoost F1", "0.903", "ML validation of k-mer features"),
    ("Concordance", "99.5%", "Composition vs known Inc groups"),
]

for i, (name, value, desc) in enumerate(metrics):
    x = Inches(0.7 + i * 3.1)
    add_text(slide, x, Inches(6.1), Inches(2.8), Inches(0.35),
             f"{name}: {value}", font_size=14, bold=True, color=DARK_BLUE)
    add_text(slide, x, Inches(6.45), Inches(2.8), Inches(0.3),
             desc, font_size=11, color=MED_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 9: Advanced Features (NEW)
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Advanced Analytical Features", "Epidemiological intelligence integrated into the pLIN pipeline")

# Feature 1: Adaptive Thresholds
add_box(slide, Inches(0.4), Inches(1.5), Inches(6.2), Inches(2.5), WHITE, border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(1.6), Inches(5.8), Inches(0.4),
         "Per-Inc-Group Adaptive Thresholds", font_size=16, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(0.6), Inches(2.1), Inches(5.8), Inches(1.6), [
    "Calibrates pLIN distance thresholds per Inc group from training data",
    "Uses quantile-based calibration on within-group distance distributions",
    "Addresses the limitation that fixed thresholds may not suit all Inc groups",
    "Example: IncFII (diverse, 4,581 plasmids) needs wider thresholds than IncX1 (compact, 701)",
    "Automatically selects thresholds for the dominant Inc group in each analysis",
], font_size=11, color=DARK_GRAY)

# Feature 2: Linkage Method Selection
add_box(slide, Inches(6.9), Inches(1.5), Inches(6.0), Inches(2.5), WHITE, border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.1), Inches(1.6), Inches(5.6), Inches(0.4),
         "Multi-Linkage Clustering", font_size=16, bold=True, color=GREEN)
add_multiline(slide, Inches(7.1), Inches(2.1), Inches(5.6), Inches(1.6), [
    "User-selectable linkage method: Single, Complete, Average, Weighted",
    "Single (default): traditional chaining, sensitive to intermediates",
    "Complete: maximum inter-cluster distance, produces tighter clusters",
    "Average (UPGMA): balanced, widely used in phylogenetics",
    "Reduces chaining artifacts that can merge distinct lineages",
], font_size=11, color=DARK_GRAY)

# Feature 3: Mobility Prediction
add_box(slide, Inches(0.4), Inches(4.3), Inches(6.2), Inches(2.8), WHITE, border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(4.4), Inches(5.8), Inches(0.4),
         "Plasmid Mobility Prediction", font_size=16, bold=True, color=ORANGE)
add_multiline(slide, Inches(0.6), Inches(4.9), Inches(5.8), Inches(2.0), [
    "3-tier priority cascade: MOBsuite \u2192 AMRFinderPlus \u2192 Non-mobilizable",
    "MOBsuite mob_typer: relaxase families (MOBF/H/P/Q/C/V) + MPF types",
    "AMRFinderPlus scan: tra/trb/mob/virB1-11/trwA-N/pilX/taxC/nikC-E",
    "Conjugative + AMR plasmids flagged as HIGH RISK for dissemination",
    "Relaxase family and MPF type charts in Epidemiology tab",
    "Enables risk-stratified surveillance of AMR-carrying plasmids",
], font_size=11, color=DARK_GRAY)

# Feature 4: Outbreak Detection
add_box(slide, Inches(6.9), Inches(4.3), Inches(6.0), Inches(2.8), WHITE, border_color=RED, border_width=Pt(2))
add_text(slide, Inches(7.1), Inches(4.4), Inches(5.6), Inches(0.4),
         "Outbreak / Clone Detection", font_size=16, bold=True, color=RED)
add_multiline(slide, Inches(7.1), Inches(4.9), Inches(5.6), Inches(2.0), [
    "Flags suspected outbreak clusters automatically",
    "Detection criteria: same pLIN strain code (F-level) + identical AMR profile",
    "Risk levels: HIGH (3+ shared AMR genes) / MODERATE (1-2 shared genes)",
    "Supports real-time surveillance for plasmid-mediated AMR spread",
    "Integrated into new Epidemiology tab with dissemination risk matrix",
    "Exportable as JSON for downstream analysis pipelines",
], font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 10: MOBsuite Mobility Typing (NEW)
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "MOBsuite Mobility Typing",
               "Relaxase family classification and Mating Pair Formation typing")

# Left: 3-tier cascade
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.5),
         "3-Tier Mobility Classification Cascade", font_size=18, bold=True, color=ORANGE)

tiers = [
    ("Priority 1: MOBsuite", "mob_typer classifies relaxase family (MOBF, MOBH, MOBP, MOBQ, MOBC, MOBV)\n"
     "and MPF type (Type T, F, I, G). Most reliable when available.",
     ORANGE, "HIGHEST"),
    ("Priority 2: AMRFinderPlus", "Scans detected genes for transfer/mobilization markers:\n"
     "Conjugative: tra/trb + virB1-11, trwA-N, pilX, taxC | Mobilizable: mob, nikC-E, oriT",
     MED_BLUE, "MEDIUM"),
    ("Priority 3: Non-mobilizable", "If no transfer genes detected by either method,\n"
     "plasmid classified as Non-mobilizable (default).",
     MED_GRAY, "DEFAULT"),
]

for i, (title, desc, color, priority) in enumerate(tiers):
    y = Inches(2.2 + i * 1.5)
    add_box(slide, Inches(0.5), y, Inches(6.0), Inches(1.3), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, Inches(0.7), y + Inches(0.1), Inches(4.5), Inches(0.35),
             title, font_size=14, bold=True, color=color)
    add_text(slide, Inches(5.0), y + Inches(0.1), Inches(1.3), Inches(0.35),
             priority, font_size=10, bold=True, color=color, alignment=PP_ALIGN.RIGHT)
    add_text(slide, Inches(0.7), y + Inches(0.5), Inches(5.6), Inches(0.7),
             desc, font_size=11, color=DARK_GRAY)
    if i < len(tiers) - 1:
        add_arrow(slide, Inches(3.2), y + Inches(1.3), Inches(0.4), Inches(0.2), color=MED_GRAY)

# Right: Relaxase families + MPF types
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.5),
         "Relaxase Families & MPF Types", font_size=18, bold=True, color=ORANGE)

relaxases = [
    ("MOBF", "F-type relaxases (IncF family)"),
    ("MOBH", "H-type relaxases (IncHI family)"),
    ("MOBP", "P-type relaxases (broad host range)"),
    ("MOBQ", "Q-type relaxases (small plasmids)"),
    ("MOBC", "C-type relaxases (ColE-like)"),
    ("MOBV", "V-type relaxases (Vibrio-associated)"),
]

for i, (name, desc) in enumerate(relaxases):
    y = Inches(2.2 + i * 0.5)
    add_text(slide, Inches(7.2), y, Inches(1.2), Inches(0.35),
             name, font_size=12, bold=True, color=ORANGE, font_name="Consolas")
    add_text(slide, Inches(8.5), y, Inches(4.3), Inches(0.35),
             desc, font_size=11, color=DARK_GRAY)

add_text(slide, Inches(7.0), Inches(5.3), Inches(6), Inches(0.4),
         "MPF Types (Mating Pair Formation)", font_size=14, bold=True, color=ORANGE)

mpf_types = [("Type T", "T4SS-related"), ("Type F", "F-pilus"), ("Type I", "I-pilus"), ("Type G", "Gram-positive")]
for i, (mtype, desc) in enumerate(mpf_types):
    x = Inches(7.0 + i * 1.55)
    add_box(slide, x, Inches(5.7), Inches(1.4), Inches(0.7), WHITE, border_color=ORANGE, border_width=Pt(1))
    add_text(slide, x, Inches(5.75), Inches(1.4), Inches(0.3),
             mtype, font_size=11, bold=True, color=ORANGE, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, Inches(6.05), Inches(1.4), Inches(0.3),
             desc, font_size=9, color=MED_GRAY, alignment=PP_ALIGN.CENTER)

# Risk stratification box at bottom
add_box(slide, Inches(0.5), Inches(6.5), Inches(12.4), Inches(0.7), RGBColor(0xFF, 0xEB, 0xEE), border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.55), Inches(5), Inches(0.3),
         "Risk Stratification:", font_size=12, bold=True, color=RED)
add_text(slide, Inches(0.7), Inches(6.85), Inches(12), Inches(0.3),
         "Conjugative + AMR \u2192 HIGH RISK    |    Mobilizable + AMR \u2192 MODERATE RISK    |    Non-mobilizable \u2192 LOWER RISK",
         font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 11: CRISPR Host Inference (NEW)
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "CRISPR Spacer-Based Host Inference",
               "Inferring plasmid-host relationships from CRISPR array matching")

# Left: Pipeline
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.5),
         "Host Inference Pipeline", font_size=18, bold=True, color=TEAL)

crispr_steps = [
    ("1", "Upload Host Genomes", "Bacterial genome FASTAs (one or more hosts)", TEAL),
    ("2", "Extract CRISPR Spacers", "MinCED identifies CRISPR arrays and extracts spacer sequences", GREEN),
    ("3", "BLAST Spacers vs Plasmids", "BLASTN-short with -dust no -word_size 7 -evalue 1e-5", MED_BLUE),
    ("4", "Stringent Filtering", "\u226595% identity, \u226525bp alignment, \u22641 mismatch, 0 gaps", ORANGE),
    ("5", "Probability Ranking", "Softmax normalization per host-plasmid pair (temperature=1.0)", PURPLE),
]

for i, (num, title, desc, color) in enumerate(crispr_steps):
    y = Inches(2.2 + i * 0.9)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y, Inches(0.35), Inches(0.35))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(12)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE

    add_text(slide, Inches(1.0), y - Inches(0.02), Inches(5.5), Inches(0.3),
             title, font_size=13, bold=True, color=color)
    add_text(slide, Inches(1.0), y + Inches(0.28), Inches(5.5), Inches(0.3),
             desc, font_size=10, color=MED_GRAY)

    if i < len(crispr_steps) - 1:
        line = slide.shapes.add_connector(1, Inches(0.68), y + Inches(0.35), Inches(0.68), y + Inches(0.9))
        line.line.color.rgb = LIGHT_BLUE
        line.line.width = Pt(2)

# Right: Features
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.5),
         "Key Features", font_size=18, bold=True, color=TEAL)

feature_boxes = [
    ("Two Source Modes", "Match spacers against uploaded plasmids\nor built-in reference DB (72,556 plasmids)",
     RGBColor(0xE0, 0xF2, 0xF1), TEAL),
    ("Confidence Categories", "High (\u22650.7): Strong CRISPR evidence\n"
     "Medium (\u22650.4): Moderate support\nLow (<0.4): Weak association",
     RGBColor(0xE8, 0xF5, 0xE9), GREEN),
    ("Rich Visualizations", "Host-plasmid probability heatmap\n"
     "Confidence pie chart, spacers bar chart\nFiltered BLAST hits expander",
     RGBColor(0xE3, 0xF2, 0xFD), MED_BLUE),
    ("Full Export", "CRISPR host predictions TSV\n"
     "Extracted spacers TSV\nIncluded in ZIP bundle",
     RGBColor(0xF3, 0xE5, 0xF5), PURPLE),
]

for i, (title, desc, bg_color, text_color) in enumerate(feature_boxes):
    y = Inches(2.2 + i * 1.2)
    add_box(slide, Inches(7.0), y, Inches(5.8), Inches(1.05), bg_color, border_color=text_color, border_width=Pt(1))
    add_text(slide, Inches(7.2), y + Inches(0.05), Inches(5.4), Inches(0.3),
             title, font_size=13, bold=True, color=text_color)
    add_text(slide, Inches(7.2), y + Inches(0.35), Inches(5.4), Inches(0.65),
             desc, font_size=10, color=DARK_GRAY)

# Bottom info box
add_box(slide, Inches(0.5), Inches(6.6), Inches(12.4), Inches(0.6), RGBColor(0xE0, 0xF2, 0xF1), border_color=TEAL, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(6.65), Inches(12), Inches(0.5),
         "Tools: MinCED 0.4.2 (spacer extraction) + BLASTN (NCBI BLAST+) | "
         "Spacer IDs: {genome}__spacer_{N} for host tracing | "
         "Numerically stable softmax with max-subtraction",
         font_size=10, color=TEAL)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 12: Nucleotide Transformer LLM Integration
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Nucleotide Transformer LLM Integration",
               "Optional genomic language model for enhanced Inc group and AMR prediction")

# Architecture diagram area
add_box(slide, Inches(0.4), Inches(1.5), Inches(12.5), Inches(2.8), LIGHT_GRAY,
        border_color=MED_GRAY, border_width=Pt(1))
add_text(slide, Inches(0.6), Inches(1.55), Inches(12), Inches(0.35),
         "Dual Prediction Architecture", font_size=16, bold=True, color=DARK_BLUE)

# Left pipeline: Traditional KNN
add_box(slide, Inches(0.7), Inches(2.1), Inches(5.5), Inches(1.9), WHITE,
        border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.9), Inches(2.15), Inches(5.1), Inches(0.35),
         "Traditional Pipeline (KNN)", font_size=14, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(0.9), Inches(2.55), Inches(5.1), Inches(1.3), [
    "FASTA \u2192 4-mer frequency vectors (256D)",
    "KNN classifier (k=5, cosine distance)",
    "92.2% cross-validation accuracy (20 groups)",
    "Instant inference (\u003C100ms)",
], font_size=11, color=DARK_GRAY)

# Arrow between pipelines
arrow = slide.shapes.add_shape(MSO_SHAPE.LEFT_RIGHT_ARROW,
                                Inches(6.35), Inches(2.8), Inches(0.6), Inches(0.35))
arrow.fill.solid()
arrow.fill.fore_color.rgb = DARK_GRAY
arrow.line.fill.background()

# Right pipeline: NT LLM
add_box(slide, Inches(7.1), Inches(2.1), Inches(5.5), Inches(1.9), WHITE,
        border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(7.3), Inches(2.15), Inches(5.1), Inches(0.35),
         "Nucleotide Transformer (LLM)", font_size=14, bold=True, color=PURPLE)
add_multiline(slide, Inches(7.3), Inches(2.55), Inches(5.1), Inches(1.3), [
    "FASTA \u2192 6-mer tokenization \u2192 Transformer",
    "Chunked embedding (5kb chunks, mean-pooled)",
    "Linear probe classifiers on frozen embeddings",
    "Models: 50M / 100M / 250M / 500M params",
], font_size=11, color=DARK_GRAY)

# Bottom row: feature boxes
features_nt = [
    ("Sequence Understanding", "NT captures long-range sequence\ncontext and motif patterns\nbeyond simple k-mer counts", MED_BLUE),
    ("Chunked Embedding", "Plasmids split into 5,000 bp\noverlapping chunks (stride 2,500)\nMean-pooled across chunks", GREEN),
    ("Inc Group Prediction", "LogisticRegression probe on\nNT embeddings predicts Inc\ngroup with CV accuracy report", ORANGE),
    ("AMR Class Prediction", "Multi-label probe predicts\nAMR drug classes from\nsequence composition alone", RED),
    ("KNN vs NT Comparison", "Agreement rate shown in\nOverview tab. Disagreements\nflagged for manual review", PURPLE),
    ("Graceful Degradation", "Fully optional \u2014 works without\ntransformers/torch installed.\nFallback to KNN if NT fails", TEAL),
]

for i, (title, desc, color) in enumerate(features_nt):
    x = Inches(0.4 + i * 2.15)
    y = Inches(4.6)
    add_box(slide, x, y, Inches(2.0), Inches(0.4), color)
    add_text(slide, x, y + Inches(0.05), Inches(2.0), Inches(0.3),
             title, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, y + Inches(0.4), Inches(2.0), Inches(1.4), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.05), y + Inches(0.5), Inches(1.9), Inches(1.2),
                  desc.split("\n"), font_size=9, color=DARK_GRAY)

# Training workflow at bottom
add_box(slide, Inches(0.4), Inches(6.3), Inches(12.5), Inches(0.9), WHITE,
        border_color=DARK_BLUE, border_width=Pt(1))
add_text(slide, Inches(0.6), Inches(6.35), Inches(12), Inches(0.3),
         "Training Workflow", font_size=13, bold=True, color=DARK_BLUE)
add_text(slide, Inches(0.6), Inches(6.65), Inches(12), Inches(0.45),
         "pip install transformers torch  \u2192  python train_nt_classifier.py --model 50m  "
         "\u2192  Extracts embeddings from training data  \u2192  Trains Inc + AMR probes  "
         "\u2192  Saves data/nt_inc_probe.pkl & data/nt_amr_probe.pkl",
         font_size=10, color=MED_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 11: Deployment & Distribution
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Deployment & Distribution",
               "Multiple deployment options for different user environments")

# Row 1: Three deployment options
deploy_options = [
    ("One-Click Desktop Launch",
     "Platform-specific launchers that\nauto-install dependencies and\nlaunch the Streamlit web app",
     [
         "Windows: launch_pLIN.bat (double-click)",
         "macOS: launch_pLIN.command (double-click)",
         "Linux: launch_pLIN.sh (chmod +x, run)",
         "Auto-creates virtual environment",
         "Installs all dependencies on first run",
         "No command-line experience required",
     ], MED_BLUE),
    ("Docker Container",
     "Containerized deployment for\nreproducible, isolated execution\non any cloud or server",
     [
         "docker build -t plin .",
         "docker run -p 8501:8501 plin",
         "Docker Compose support included",
         "Health check endpoint built-in",
         "Volumes for data persistence",
         "Deploy to AWS, GCP, Azure, etc.",
     ], GREEN),
    ("Streamlit Community Cloud",
     "Free public web app deployment\ndirectly from GitHub repository\n\u2014 zero infrastructure needed",
     [
         "Deploy at share.streamlit.io",
         "Sign in with GitHub account",
         "Select repo + plin_app.py",
         "Public URL auto-generated",
         "Auto-redeploys on git push",
         "Free tier available",
     ], PURPLE),
]

for i, (title, subtitle, items, color) in enumerate(deploy_options):
    x = Inches(0.4 + i * 4.25)
    # Title box
    add_box(slide, x, Inches(1.5), Inches(3.95), Inches(0.5), color)
    add_text(slide, x, Inches(1.55), Inches(3.95), Inches(0.4),
             title, font_size=14, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    # Subtitle
    add_box(slide, x, Inches(2.0), Inches(3.95), Inches(1.0), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.1), Inches(2.1), Inches(3.75), Inches(0.85),
                  subtitle.split("\n"), font_size=10, color=MED_GRAY)
    # Items
    add_box(slide, x, Inches(3.0), Inches(3.95), Inches(2.3), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.1), Inches(3.1), Inches(3.75), Inches(2.1),
                  items, font_size=10, color=DARK_GRAY)

# Architecture flow at bottom
add_box(slide, Inches(0.4), Inches(5.6), Inches(12.5), Inches(1.6), LIGHT_GRAY,
        border_color=MED_GRAY, border_width=Pt(1))
add_text(slide, Inches(0.6), Inches(5.65), Inches(12), Inches(0.35),
         "Distribution Architecture", font_size=14, bold=True, color=DARK_BLUE)

# Flow boxes
flow_items = [
    ("GitHub Repo", "Source code\nLICENSE\nCITATION.cff", DARK_BLUE),
    ("pip install", "requirements.txt\nPython 3.10+\nvenv isolation", MED_BLUE),
    ("Streamlit App", "plin_app.py\n6 analysis tabs\nInteractive GUI", GREEN),
    ("Docker Image", "Dockerfile\nSelf-contained\nCloud-ready", ORANGE),
    ("Public Cloud", "Streamlit Cloud\nPublic URL\nAuto-deploy", PURPLE),
]

for i, (title, desc, color) in enumerate(flow_items):
    x = Inches(0.6 + i * 2.5)
    add_box(slide, x, Inches(6.1), Inches(2.1), Inches(0.35), color)
    add_text(slide, x, Inches(6.12), Inches(2.1), Inches(0.3),
             title, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    add_box(slide, x, Inches(6.45), Inches(2.1), Inches(0.65), WHITE,
            border_color=color, border_width=Pt(1))
    add_multiline(slide, x + Inches(0.05), Inches(6.5), Inches(2.0), Inches(0.55),
                  desc.split("\n"), font_size=8, color=DARK_GRAY)

    # Arrow between flow items
    if i < len(flow_items) - 1:
        arr = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW,
                                      x + Inches(2.15), Inches(6.55), Inches(0.3), Inches(0.2))
        arr.fill.solid()
        arr.fill.fore_color.rgb = DARK_GRAY
        arr.line.fill.background()


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 12: Novelty
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Novelty of pLIN System", "What distinguishes pLIN from existing plasmid classification approaches")

# Left column: Novelty points
novelty_items = [
    ("Reference-Free Classification",
     "Unlike MOB-suite, PlasmidFinder, or replicon typing that rely on curated reference databases, "
     "pLIN uses intrinsic sequence composition (4-mer frequencies). New or divergent plasmids that lack "
     "known replicons are still classifiable — no database gaps.",
     MED_BLUE),
    ("Permanent, Hierarchical Codes",
     "pLIN codes are stable and never change when new sequences are added. Existing plasmid typing schemes "
     "(e.g., pMLST) can reassign types as databases grow. The six-level hierarchy (L1 \u2192 L6) provides "
     "resolution at multiple scales in a single code — no existing system offers this.",
     GREEN),
    ("Composition-Based Inc Group Detection",
     "Most Inc group detection tools (PlasmidFinder, COPLA) require gene-level BLAST searches. pLIN's KNN "
     "classifier predicts Inc group from global k-mer composition alone (92% accuracy, 20 groups), enabling classification "
     "even when replicon genes are fragmented, truncated, or absent from assemblies.",
     ORANGE),
    ("Integrated AMR Surveillance Pipeline",
     "pLIN is the first system to natively integrate plasmid hierarchical classification with AMRFinderPlus "
     "gene detection in a single automated pipeline. This enables direct mapping of resistance gene cargo "
     "onto the plasmid phylogeny — linking plasmid backbone evolution with AMR gene acquisition.",
     RED),
    ("Interactive GUI with Epidemiological Intelligence",
     "Unlike command-line-only tools (MOB-suite, plasmidfinder CLI), pLIN ships a Streamlit web app with "
     "upload, classification, AMR screening, cladogram visualization, mobility prediction, outbreak detection, "
     "adaptive thresholds, and multi-linkage clustering — accessible to non-bioinformaticians.",
     PURPLE),
    ("Genomic LLM Integration (Nucleotide Transformer)",
     "First plasmid classification tool to integrate a genomic foundation model (InstaDeep Nucleotide Transformer). "
     "Provides LLM-based Inc group and AMR class predictions alongside traditional KNN, enabling cross-method "
     "validation. Chunked embedding strategy handles full-length plasmids (up to 300+ kb).",
     TEAL),
]

for i, (title, desc, color) in enumerate(novelty_items):
    y = Inches(1.5 + i * 0.95)
    # Colored left bar
    add_box(slide, Inches(0.4), y, Inches(0.12), Inches(0.95), color)
    add_text(slide, Inches(0.7), y + Inches(0.02), Inches(12.2), Inches(0.35),
             title, font_size=15, bold=True, color=color)
    add_text(slide, Inches(0.7), y + Inches(0.38), Inches(12.2), Inches(0.55),
             desc, font_size=11, color=DARK_GRAY)

# Comparison table at bottom
add_box(slide, Inches(0.4), Inches(7.0), Inches(12.5), Inches(0.35), DARK_BLUE)
comp_headers = ["Feature", "pLIN", "MOB-suite", "PlasmidFinder", "pMLST"]
col_widths = [2.8, 2.2, 2.5, 2.5, 2.5]
x_pos = 0.4
for j, (h, w) in enumerate(zip(comp_headers, col_widths)):
    add_text(slide, Inches(x_pos + 0.05), Inches(7.02), Inches(w), Inches(0.3),
             h, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    x_pos += w


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 13: Limitations & Future Directions
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Limitations & Future Directions", "Current constraints and planned improvements")

# Limitations (left side)
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.5),
         "Current Limitations", font_size=20, bold=True, color=RED)

limitations = [
    ("Uneven Inc Group Class Sizes",
     "Classifier covers 20 Inc groups but with uneven training data: IncFII has 4,629 samples "
     "while IncFIBK has only 11. Smaller classes may have lower per-class accuracy. "
     "Balanced augmentation is needed for underrepresented groups."),
    ("Composition-Only Features",
     "4-mer frequency captures global sequence composition but ignores gene content, synteny, "
     "and structural rearrangements. Two plasmids with similar base composition but different "
     "gene cargo may receive similar pLIN codes at coarse levels."),
    ("Linkage Method Sensitivity",
     "Single-linkage clustering (default) can produce chain-like clusters. While the tool supports "
     "complete, average, and weighted linkage as alternatives, the pLIN thresholds were originally "
     "calibrated for single-linkage \u2014 using other methods may require re-calibration."),
    ("No Fragmented Assembly Handling",
     "pLIN expects complete or near-complete plasmid sequences. Short contigs from fragmented "
     "assemblies will have noisy k-mer profiles, reducing classification accuracy. No scaffolding "
     "or multi-contig plasmid reconstruction is performed."),
    ("CRISPR Host Inference Depends on DB",
     "CRISPR-based host prediction quality depends on completeness of the host genome "
     "database and CRISPR array presence. Hosts lacking CRISPR systems will not be detected. "
     "Reference DB mode requires the 72,556-plasmid sequences.fasta file."),
    ("No Real-Time Database Updates",
     "The KNN classifier uses a static training set. New plasmid submissions to GenBank/RefSeq "
     "are not automatically incorporated. Periodic retraining is needed to maintain accuracy "
     "as novel plasmid lineages emerge."),
]

for i, (title, desc) in enumerate(limitations):
    y = Inches(2.15 + i * 0.82)
    # Red number circle
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
    add_text(slide, Inches(0.95), y + Inches(0.3), Inches(5.8), Inches(0.5),
             desc, font_size=9, color=DARK_GRAY)

# Future Directions (right side)
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.5),
         "Future Directions", font_size=20, bold=True, color=GREEN)

futures = [
    ("Balance Inc Group Training",
     "Expand underrepresented groups (IncFIBK, ColE, IncI2) with targeted "
     "data collection. Target: 30+ Inc groups, 20,000+ balanced training sequences.",
     GREEN),
    ("Hybrid Features",
     "Combine k-mer composition with gene presence/absence and synteny features "
     "for higher-resolution classification at fine levels (E, F).",
     MED_BLUE),
    ("Threshold Cross-Validation",
     "Systematically validate adaptive thresholds across all 20 Inc groups. "
     "Benchmark linkage methods against known plasmid phylogenies.",
     ORANGE),
    ("Metagenomic Support",
     "Handle multi-contig plasmid bins from metagenomic assemblies. "
     "Integrate with plasmid prediction tools (PlasFlow, MOB-recon).",
     PURPLE),
    ("Online Database & API",
     "Web database for pLIN code lookups and new sequence submission. "
     "REST API for programmatic access and LIMS integration.",
     TEAL),
    ("Multi-Hospital Outbreak Tracking",
     "Extend outbreak detection to cross-institutional datasets. "
     "Integrate with epidemiological metadata for spatial-temporal tracking.",
     RED),
]

for i, (title, desc, color) in enumerate(futures):
    y = Inches(2.15 + i * 0.82)
    # Arrow shape
    arrow = slide.shapes.add_shape(MSO_SHAPE.RIGHT_ARROW, Inches(7.0), y + Inches(0.1), Inches(0.35), Inches(0.25))
    arrow.fill.solid()
    arrow.fill.fore_color.rgb = color
    arrow.line.fill.background()

    add_text(slide, Inches(7.5), y, Inches(5.3), Inches(0.3),
             title, font_size=12, bold=True, color=color)
    add_text(slide, Inches(7.5), y + Inches(0.3), Inches(5.3), Inches(0.5),
             desc, font_size=9, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 14: DRAGNOME Buddy — AI Chatbot
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "DRAGNOME Buddy — AI Assistant",
               "Local LLM-powered chatbot for plasmid biology Q&A")

# Left side: Architecture diagram
add_text(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(0.4),
         "How It Works", font_size=18, bold=True, color=DARK_BLUE)

# Chat flow boxes
flow_items = [
    ("User Question", "Natural language query about\nplasmids or analysis results", LIGHT_BLUE, DARK_BLUE),
    ("Context Builder", "Extracts analysis data:\npLIN codes, Inc groups, AMR,\nmobility, outbreaks", RGBColor(0xE8, 0xF5, 0xE9), GREEN),
    ("Ollama LLM", "Local model (llama3.2, mistral)\nNo API keys, data stays local", RGBColor(0xFF, 0xF3, 0xE0), ORANGE),
    ("Streamed Response", "Real-time answer with\ncontext-aware explanations", RGBColor(0xE3, 0xF2, 0xFD), MED_BLUE),
]

for i, (title, desc, fill, border) in enumerate(flow_items):
    y = Inches(2.1 + i * 1.2)
    box = add_box(slide, Inches(0.5), y, Inches(5.5), Inches(1.0), fill, border_color=border, border_width=Pt(2))
    add_text(slide, Inches(0.7), y + Inches(0.15), Inches(5.1), Inches(0.3),
             title, font_size=14, bold=True, color=border)
    add_text(slide, Inches(0.7), y + Inches(0.45), Inches(5.1), Inches(0.5),
             desc, font_size=11, color=DARK_GRAY)
    if i < len(flow_items) - 1:
        add_arrow(slide, Inches(3.0), y + Inches(1.0), Inches(0.4), Inches(0.2), border)

# Right side: Features
add_text(slide, Inches(6.5), Inches(1.6), Inches(6.5), Inches(0.4),
         "Key Features", font_size=18, bold=True, color=DARK_BLUE)

features = [
    ("Privacy-First", "Runs 100% locally via Ollama — no data leaves your machine", GREEN),
    ("Context-Aware", "Understands your analysis: Inc groups, AMR genes, outbreaks", MED_BLUE),
    ("Multiple Models", "Choose: llama3.2, mistral, mixtral, phi3, gemma2", ORANGE),
    ("Streaming Chat", "Real-time response streaming for smooth UX", PURPLE),
    ("Suggested Prompts", "Pre-written questions to get started quickly", TEAL),
]

for i, (title, desc, color) in enumerate(features):
    y = Inches(2.1 + i * 0.95)
    # Colored bullet
    bullet = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(6.5), y + Inches(0.1), Inches(0.2), Inches(0.2))
    bullet.fill.solid()
    bullet.fill.fore_color.rgb = color
    bullet.line.fill.background()
    add_text(slide, Inches(6.85), y, Inches(5.8), Inches(0.3),
             title, font_size=13, bold=True, color=color)
    add_text(slide, Inches(6.85), y + Inches(0.35), Inches(5.8), Inches(0.5),
             desc, font_size=11, color=DARK_GRAY)

# Example questions box
add_box(slide, Inches(6.5), Inches(6.0), Inches(6.3), Inches(1.2), RGBColor(0xF3, 0xE5, 0xF5),
        border_color=PURPLE, border_width=Pt(1))
add_text(slide, Inches(6.7), Inches(6.1), Inches(6), Inches(0.3),
         "Example Questions:", font_size=12, bold=True, color=PURPLE)
examples = '"Summarize my results" • "Which plasmids are high-risk?" • "Explain the AMR genes found"'
add_text(slide, Inches(6.7), Inches(6.4), Inches(5.9), Inches(0.7),
         examples, font_size=10, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 15: Unknown/Novel Inc Type Detection
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Unknown/Novel Inc Type Detection",
               "Confidence-based flagging prevents misclassification")

# Problem statement
add_box(slide, Inches(0.5), Inches(1.6), Inches(6), Inches(1.4), RGBColor(0xFF, 0xEB, 0xEE),
        border_color=RED, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(1.7), Inches(5.6), Inches(0.3),
         "The Problem", font_size=14, bold=True, color=RED)
add_text(slide, Inches(0.7), Inches(2.05), Inches(5.6), Inches(0.9),
         "KNN always assigns to the nearest class, even if the match is poor.\n"
         "Plasmids with Inc types not in training data (IncL, IncU, novel)\n"
         "get misclassified with no warning.",
         font_size=11, color=DARK_GRAY)

# Solution
add_box(slide, Inches(0.5), Inches(3.2), Inches(6), Inches(1.4), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(0.7), Inches(3.3), Inches(5.6), Inches(0.3),
         "The Solution", font_size=14, bold=True, color=GREEN)
add_text(slide, Inches(0.7), Inches(3.65), Inches(5.6), Inches(0.9),
         "Confidence threshold (40%): predictions below this are flagged\n"
         'as "Unknown/Novel" instead of being forced into a wrong class.\n'
         "Top 5 nearest candidates shown for manual verification.",
         font_size=11, color=DARK_GRAY)

# Flow diagram
add_text(slide, Inches(0.5), Inches(4.9), Inches(6), Inches(0.3),
         "Classification Flow", font_size=14, bold=True, color=DARK_BLUE)

flow_boxes = [
    ("Query Sequence", LIGHT_BLUE),
    ("KNN Classifier", MED_BLUE),
    ("Confidence Check", ORANGE),
]
for i, (label, color) in enumerate(flow_boxes):
    x = Inches(0.5 + i * 2.0)
    box = add_box(slide, x, Inches(5.3), Inches(1.8), Inches(0.6), color)
    add_text(slide, x, Inches(5.45), Inches(1.8), Inches(0.3),
             label, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    if i < len(flow_boxes) - 1:
        add_right_arrow(slide, x + Inches(1.8), Inches(5.45), Inches(0.2), Inches(0.3), color)

# Branch arrows and outcomes
add_text(slide, Inches(5.3), Inches(5.15), Inches(1.5), Inches(0.3),
         "≥40%", font_size=11, bold=True, color=GREEN)
add_box(slide, Inches(5.3), Inches(5.45), Inches(1.5), Inches(0.5), GREEN)
add_text(slide, Inches(5.3), Inches(5.55), Inches(1.5), Inches(0.3),
         "Assigned Inc", font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

add_text(slide, Inches(5.3), Inches(6.05), Inches(1.5), Inches(0.3),
         "<40%", font_size=11, bold=True, color=RED)
add_box(slide, Inches(5.3), Inches(6.35), Inches(1.5), Inches(0.5), RED)
add_text(slide, Inches(5.3), Inches(6.45), Inches(1.5), Inches(0.3),
         "Unknown/Novel", font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)

# Right side: Output example
add_text(slide, Inches(7.0), Inches(1.6), Inches(6), Inches(0.4),
         "Output Example", font_size=18, bold=True, color=DARK_BLUE)

# Example table
add_box(slide, Inches(7.0), Inches(2.1), Inches(6), Inches(2.8), LIGHT_GRAY,
        border_color=DARK_GRAY, border_width=Pt(1))

# Table header
add_box(slide, Inches(7.0), Inches(2.1), Inches(6), Inches(0.4), DARK_BLUE)
headers = ["Plasmid", "Inc Type", "Conf.", "Top 5 Candidates"]
widths = [1.2, 1.2, 0.7, 2.9]
x = Inches(7.1)
for h, w in zip(headers, widths):
    add_text(slide, x, Inches(2.15), Inches(w), Inches(0.3),
             h, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    x += Inches(w)

# Table rows
rows = [
    ("pBR322", "IncFII", "92%", "IncFII (92%), IncN (5%), IncX1 (2%)...", GREEN),
    ("pXYZ_novel", "Unknown", "35%", "IncFII (35%), IncN (28%), IncX1 (18%)...", RED),
    ("pABC_low", "Unknown", "41%", "IncN (41%), IncFII (30%), IncX (15%)...", ORANGE),
]
for i, (pid, inc, conf, top5, color) in enumerate(rows):
    y = Inches(2.55 + i * 0.5)
    row_color = RGBColor(0xFF, 0xEB, 0xEE) if "Unknown" in inc else WHITE
    add_box(slide, Inches(7.0), y, Inches(6), Inches(0.45), row_color)
    x = Inches(7.1)
    for val, w in zip([pid, inc, conf, top5], widths):
        c = RED if "Unknown" in inc and val == inc else (color if val == conf else DARK_GRAY)
        add_text(slide, x, y + Inches(0.1), Inches(w), Inches(0.25),
                 val, font_size=9, color=c, alignment=PP_ALIGN.CENTER)
        x += Inches(w)

# Benefits
add_text(slide, Inches(7.0), Inches(5.0), Inches(6), Inches(0.4),
         "Benefits", font_size=14, bold=True, color=DARK_BLUE)

benefits = [
    ("Honest uncertainty", "Tool admits when it doesn't know"),
    ("Prevents misclassification", "No more wrong Inc assignments"),
    ("Top 5 candidates", "Guides manual verification"),
    ("Actionable warnings", "Suggests PlasmidFinder for verification"),
]
for i, (title, desc) in enumerate(benefits):
    y = Inches(5.4 + i * 0.45)
    add_text(slide, Inches(7.2), y, Inches(2.5), Inches(0.25),
             f"✓ {title}", font_size=10, bold=True, color=GREEN)
    add_text(slide, Inches(9.7), y, Inches(3.2), Inches(0.25),
             desc, font_size=10, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 16: Complete Feature Summary
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Complete Feature Summary",
               "All capabilities of the pLIN classification system")

# Feature grid - 4 columns x 4 rows
features_grid = [
    # Row 1: Core
    [("pLIN Classification", "6-level hierarchical codes\n(L1\u2192L6)", DARK_BLUE),
     ("Inc Auto-Detection", "KNN classifier, 92.2%\naccuracy, 20 Inc groups", MED_BLUE),
     ("AMRFinderPlus", "AMR/stress/virulence\ngene detection", RED),
     ("Prodigal Annotation", "Full ORF prediction,\ncoding density stats", GREEN)],
    # Row 2: Analysis
    [("Mobility Prediction", "MOBsuite + AMRFinderPlus\n3-tier cascade", ORANGE),
     ("Outbreak Detection", "Strain + AMR profile\ncluster identification", PURPLE),
     ("CRISPR Host Inference", "Spacer-based plasmid-host\nsoftmax probability", RGBColor(0x00, 0x69, 0x5C)),
     ("Adaptive Thresholds", "Per-Inc calibrated\ndistance thresholds", TEAL)],
    # Row 3: AI/ML
    [("Nucleotide Transformer", "LLM-based Inc/AMR\nprediction (optional)", RGBColor(0x6A, 0x1B, 0x9A)),
     ("DRAGNOME Buddy", "AI chatbot for Q&A\nvia local Ollama", RGBColor(0x00, 0x69, 0x5C)),
     ("Unknown Detection", "Low-confidence flagging\n+ top 5 candidates", RGBColor(0xBF, 0x36, 0x0C)),
     ("Context-Aware AI", "Analysis results feed\ninto LLM responses", RGBColor(0x1A, 0x23, 0x7E))],
    # Row 4: Deployment
    [("Streamlit GUI", "8-tab web interface\nno CLI required", RGBColor(0xFF, 0x4B, 0x4B)),
     ("Docker Support", "Containerized deployment\nfor servers", RGBColor(0x00, 0x97, 0xA7)),
     ("One-Click Launch", "Windows/macOS/Linux\nauto-setup launchers", RGBColor(0x7B, 0x1F, 0xA2)),
     ("Export Options", "TSV, PNG, PDF, ZIP\nbundle downloads", RGBColor(0x2E, 0x7D, 0x32))],
]

for row_idx, row in enumerate(features_grid):
    for col_idx, (title, desc, color) in enumerate(row):
        x = Inches(0.4 + col_idx * 3.2)
        y = Inches(1.6 + row_idx * 1.4)
        box = add_box(slide, x, y, Inches(3.0), Inches(1.25), color)
        add_text(slide, x, y + Inches(0.15), Inches(3.0), Inches(0.35),
                 title, font_size=12, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
        add_text(slide, x, y + Inches(0.55), Inches(3.0), Inches(0.6),
                 desc, font_size=10, color=RGBColor(0xE0, 0xE0, 0xE0), alignment=PP_ALIGN.CENTER)

# Footer
add_text(slide, Inches(0.5), Inches(7.0), Inches(12.3), Inches(0.3),
         "pLIN v2.1 — Plasmid Life Identification Number System • GPL-3.0 + Citation Clause • github.com/xavierbasilbritto-hub/pLIN-plasmid-classification",
         font_size=10, color=MED_GRAY, alignment=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 17: Phase 1 Upgrades — ANI Validation & Metadata Integration
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Phase 1 Upgrades — ANI Validation & Metadata",
               "Mash/MinHash integration, metadata upload, sequence length warnings")

# Feature 1: Mash/MinHash ANI Estimation
add_box(slide, Inches(0.4), Inches(1.5), Inches(6.2), Inches(2.5), WHITE, border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(1.6), Inches(5.8), Inches(0.4),
         "Mash/MinHash ANI Estimation", font_size=16, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(0.6), Inches(2.1), Inches(5.8), Inches(1.6), [
    "Fast pairwise ANI estimation via MinHash sketching",
    "Parameters: k=21, sketch size=10,000",
    "ANI estimate = (1 - Mash distance) x 100",
    "Validates composition-based pLIN thresholds against ANI",
    "Auto-detects mash binary from PATH / conda environments",
    "Mean/min/max ANI metrics displayed in Epidemiology tab",
], font_size=11, color=DARK_GRAY)

# Feature 2: Metadata CSV Upload
add_box(slide, Inches(6.9), Inches(1.5), Inches(6.0), Inches(2.5), WHITE, border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.1), Inches(1.6), Inches(5.6), Inches(0.4),
         "Metadata CSV/TSV Upload", font_size=16, bold=True, color=GREEN)
add_multiline(slide, Inches(7.1), Inches(2.1), Inches(5.6), Inches(1.6), [
    "Upload patient/sample metadata alongside plasmid sequences",
    "Auto-detects join column: plasmid_id, filename, sample_id",
    "Auto-parses date columns for temporal outbreak clustering",
    "Merges metadata into integrated results table",
    "Enables location + date based epidemiological analysis",
    "Supports both CSV and TSV formats",
], font_size=11, color=DARK_GRAY)

# Feature 3: Sequence Length Warning
add_box(slide, Inches(0.4), Inches(4.3), Inches(6.2), Inches(2.8), WHITE, border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(4.4), Inches(5.8), Inches(0.4),
         "Sequence Length Warning System", font_size=16, bold=True, color=ORANGE)
add_multiline(slide, Inches(0.6), Inches(4.9), Inches(5.8), Inches(2.0), [
    "Threshold: 5,000 bp (SHORT_PLASMID_THRESHOLD)",
    "Short plasmids have noisy 4-mer frequency profiles",
    "Warning displayed in Overview tab with expandable detail",
    "Shows: plasmid ID, Inc type, length, pLIN code",
    "Recommends caution for classification of short sequences",
    "Does not block analysis \u2014 informational warning only",
], font_size=11, color=DARK_GRAY)

# Feature 4: Adaptive Calibration Default
add_box(slide, Inches(6.9), Inches(4.3), Inches(6.0), Inches(2.8), WHITE, border_color=PURPLE, border_width=Pt(2))
add_text(slide, Inches(7.1), Inches(4.4), Inches(5.6), Inches(0.4),
         "Adaptive Calibration (Now Default)", font_size=16, bold=True, color=PURPLE)
add_multiline(slide, Inches(7.1), Inches(4.9), Inches(5.6), Inches(2.0), [
    "Adaptive thresholds now enabled by default (was opt-in)",
    "Calibrates pLIN thresholds per Inc group automatically",
    "Uses quantile-based calibration on training distances",
    "Better handles diverse Inc groups (e.g. IncFII vs IncX1)",
    "Recommended for most analyses \u2014 fixed thresholds available",
    "Checkbox still available for universal threshold users",
], font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 18: Phase 2 Upgrades — FastANI, SNP Sub-typing & Temporal Outbreaks
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Phase 2 Upgrades — Genomic Resolution",
               "FastANI true ANI, minimap2 SNP sub-typing, temporal outbreak clustering")

# Feature 1: FastANI Integration
add_box(slide, Inches(0.4), Inches(1.5), Inches(4.0), Inches(2.7), WHITE, border_color=MED_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(1.6), Inches(3.6), Inches(0.4),
         "FastANI True ANI", font_size=16, bold=True, color=MED_BLUE)
add_multiline(slide, Inches(0.6), Inches(2.1), Inches(3.6), Inches(1.8), [
    "All-vs-all true ANI computation",
    "Optimized for plasmids: --fragLen 1000",
    "4 threads for parallel execution",
    "Pairwise ANI + fragment statistics",
    "Auto-detects fastANI binary",
    "Ground-truth validation of pLIN",
], font_size=11, color=DARK_GRAY)

# Feature 2: SNP Sub-typing within L6
add_box(slide, Inches(4.7), Inches(1.5), Inches(4.0), Inches(2.7), WHITE, border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(4.9), Inches(1.6), Inches(3.6), Inches(0.4),
         "SNP Sub-typing (L6)", font_size=16, bold=True, color=GREEN)
add_multiline(slide, Inches(4.9), Inches(2.1), Inches(3.6), Inches(1.8), [
    "minimap2 alignment within L6 clusters",
    "Preset: -cx asm5 (closely related)",
    "Counts mismatches from CS tags",
    "0-SNP pairs flagged as clonal",
    "Resolves within-strain diversity",
    "Critical for outbreak investigation",
], font_size=11, color=DARK_GRAY)

# Feature 3: Temporal Outbreak Clustering
add_box(slide, Inches(9.0), Inches(1.5), Inches(4.0), Inches(2.7), WHITE, border_color=RED, border_width=Pt(2))
add_text(slide, Inches(9.2), Inches(1.6), Inches(3.6), Inches(0.4),
         "Temporal Outbreaks", font_size=16, bold=True, color=RED)
add_multiline(slide, Inches(9.2), Inches(2.1), Inches(3.6), Inches(1.8), [
    "30-day sliding window detection",
    "Requires: L6 code + AMR + dates",
    "3 risk levels: CRITICAL/HIGH/MOD",
    "CRITICAL: \u22653 AMR + \u22647 days",
    "Leverages metadata CSV upload",
    "Real-time surveillance ready",
], font_size=11, color=DARK_GRAY)

# Bottom: Analysis Pipeline Flow
add_box(slide, Inches(0.4), Inches(4.5), Inches(12.5), Inches(2.7), LIGHT_GRAY,
        border_color=MED_GRAY, border_width=Pt(1))
add_text(slide, Inches(0.6), Inches(4.6), Inches(12), Inches(0.4),
         "Enhanced Analysis Pipeline (Steps 8\u201311)", font_size=16, bold=True, color=DARK_BLUE)

pipeline_steps = [
    ("Step 8b", "Temporal\nOutbreak\nClustering", RED),
    ("Step 9", "Mash ANI\nEstimation\n(MinHash)", MED_BLUE),
    ("Step 10", "FastANI\nTrue ANI\nComputation", GREEN),
    ("Step 11", "minimap2\nSNP Sub-\ntyping", ORANGE),
]

for i, (step, desc, color) in enumerate(pipeline_steps):
    x = Inches(0.8 + i * 3.1)
    add_box(slide, x, Inches(5.2), Inches(2.5), Inches(1.6), WHITE, border_color=color, border_width=Pt(2))
    add_text(slide, x, Inches(5.25), Inches(2.5), Inches(0.35),
             step, font_size=13, bold=True, color=color, alignment=PP_ALIGN.CENTER)
    add_text(slide, x, Inches(5.65), Inches(2.5), Inches(0.9),
             desc, font_size=11, color=DARK_GRAY, alignment=PP_ALIGN.CENTER)

    if i < len(pipeline_steps) - 1:
        add_right_arrow(slide, x + Inches(2.55), Inches(5.8), Inches(0.5), Inches(0.3), color=MED_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 19: Complete Upgrade Summary
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Upgrade Summary — v2.1 Feature Matrix",
               "7 new features across Phase 1 (immediate) and Phase 2 (advanced)")

# Phase 1 header
add_box(slide, Inches(0.4), Inches(1.5), Inches(6.2), Inches(0.5), MED_BLUE)
add_text(slide, Inches(0.5), Inches(1.55), Inches(6.0), Inches(0.4),
         "Phase 1 \u2014 Immediate, High-Value", font_size=16, bold=True, color=WHITE)

phase1_items = [
    ("1", "Adaptive Calibration Default", "Per-Inc thresholds now enabled by default", MED_BLUE),
    ("2", "Sequence Length Warning", "Flags plasmids < 5 kb with noisy profiles", ORANGE),
    ("3", "Metadata CSV Upload", "Patient/sample data for epi analysis", GREEN),
    ("4", "Mash/MinHash ANI", "Fast ANI validation (k=21, s=10000)", PURPLE),
]

for i, (num, title, desc, color) in enumerate(phase1_items):
    y = Inches(2.15 + i * 0.6)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.05), Inches(0.3), Inches(0.3))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(11)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    add_text(slide, Inches(0.95), y, Inches(2.5), Inches(0.3),
             title, font_size=12, bold=True, color=color)
    add_text(slide, Inches(0.95), y + Inches(0.28), Inches(5.5), Inches(0.25),
             desc, font_size=10, color=MED_GRAY)

# Phase 2 header
add_box(slide, Inches(0.4), Inches(4.6), Inches(6.2), Inches(0.5), RED)
add_text(slide, Inches(0.5), Inches(4.65), Inches(6.0), Inches(0.4),
         "Phase 2 \u2014 Advanced Genomic Resolution", font_size=16, bold=True, color=WHITE)

phase2_items = [
    ("5", "FastANI True ANI", "Ground-truth ANI with --fragLen 1000", MED_BLUE),
    ("6", "SNP Sub-typing (L6)", "minimap2 -cx asm5 within L6 clusters", GREEN),
    ("7", "Temporal Outbreak Clustering", "30-day window + AMR fingerprint matching", RED),
]

for i, (num, title, desc, color) in enumerate(phase2_items):
    y = Inches(5.25 + i * 0.6)
    circ = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(0.5), y + Inches(0.05), Inches(0.3), Inches(0.3))
    circ.fill.solid()
    circ.fill.fore_color.rgb = color
    circ.line.fill.background()
    tf = circ.text_frame
    tf.paragraphs[0].text = num
    tf.paragraphs[0].font.size = Pt(11)
    tf.paragraphs[0].font.bold = True
    tf.paragraphs[0].font.color.rgb = WHITE
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    add_text(slide, Inches(0.95), y, Inches(2.5), Inches(0.3),
             title, font_size=12, bold=True, color=color)
    add_text(slide, Inches(0.95), y + Inches(0.28), Inches(5.5), Inches(0.25),
             desc, font_size=10, color=MED_GRAY)

# Right side: Impact matrix
add_text(slide, Inches(7.0), Inches(1.5), Inches(6), Inches(0.5),
         "Impact & Tool Requirements", font_size=18, bold=True, color=DARK_BLUE)

# Impact table header
add_box(slide, Inches(7.0), Inches(2.1), Inches(5.8), Inches(0.4), DARK_BLUE)
imp_headers = ["Feature", "Tool", "Impact"]
imp_widths = [2.0, 1.5, 2.3]
x = Inches(7.1)
for h, w in zip(imp_headers, imp_widths):
    add_text(slide, x, Inches(2.15), Inches(w), Inches(0.3),
             h, font_size=10, bold=True, color=WHITE, alignment=PP_ALIGN.CENTER)
    x += Inches(w)

impact_rows = [
    ("Adaptive Default", "Built-in", "Better per-Inc accuracy"),
    ("Length Warning", "Built-in", "User confidence"),
    ("Metadata Upload", "Built-in", "Epi context"),
    ("Mash ANI", "mash", "ANI validation"),
    ("FastANI", "fastANI", "True ANI"),
    ("SNP Sub-typing", "minimap2", "Outbreak resolution"),
    ("Temporal Clusters", "Built-in", "Surveillance"),
]

for i, (feat, tool, impact) in enumerate(impact_rows):
    y = Inches(2.55 + i * 0.42)
    bg = LIGHT_GRAY if i % 2 == 0 else WHITE
    add_box(slide, Inches(7.0), y, Inches(5.8), Inches(0.42), bg)
    x = Inches(7.1)
    for val, w in zip([feat, tool, impact], imp_widths):
        c = GREEN if val == "Built-in" else DARK_GRAY
        add_text(slide, x, y + Inches(0.08), Inches(w), Inches(0.25),
                 val, font_size=10, color=c, alignment=PP_ALIGN.CENTER)
        x += Inches(w)

# Version badge
add_box(slide, Inches(7.0), Inches(5.7), Inches(5.8), Inches(1.5), RGBColor(0xE8, 0xF5, 0xE9),
        border_color=GREEN, border_width=Pt(2))
add_text(slide, Inches(7.2), Inches(5.8), Inches(5.4), Inches(0.4),
         "pLIN v2.1 \u2014 7 New Features", font_size=16, bold=True, color=GREEN)
add_multiline(slide, Inches(7.2), Inches(6.3), Inches(5.4), Inches(0.8), [
    "All Phase 2 tools are optional \u2014 graceful degradation when not installed",
    "Phase 1 features require no additional software",
    "Total: 19 slides, 28+ features, full GUI integration",
], font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# Save
# ══════════════════════════════════════════════════════════════════════════════

out_path = os.path.join(os.path.dirname(__file__), "output", "pLIN_Tool_Architecture.pptx")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
prs.save(out_path)
print(f"Saved: {out_path}")
print(f"Slides: {len(prs.slides)}")
