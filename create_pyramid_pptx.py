#!/usr/bin/env python3
"""
Create a Pyramid Principle PowerPoint presentation for pLIN.

Structure:
1. SITUATION: The AMR crisis and plasmid-mediated resistance
2. COMPLICATION: Current plasmid classification methods are inadequate
3. QUESTION: How can we improve plasmid surveillance?
4. ANSWER: pLIN - A hierarchical, permanent classification system

Author: Basil Xavier
"""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import nsmap
from pptx.oxml import parse_xml

# ══════════════════════════════════════════════════════════════════════════════
# COLOR PALETTE
# ══════════════════════════════════════════════════════════════════════════════

WHITE = RGBColor(255, 255, 255)
BLACK = RGBColor(0, 0, 0)
DARK_BLUE = RGBColor(0x0D, 0x47, 0xA1)
MED_BLUE = RGBColor(0x19, 0x76, 0xD2)
LIGHT_BLUE = RGBColor(0xE3, 0xF2, 0xFD)
DARK_GRAY = RGBColor(0x42, 0x42, 0x42)
MED_GRAY = RGBColor(0x75, 0x75, 0x75)
LIGHT_GRAY = RGBColor(0xF5, 0xF5, 0xF5)
GREEN = RGBColor(0x2E, 0x7D, 0x32)
LIGHT_GREEN = RGBColor(0xE8, 0xF5, 0xE9)
RED = RGBColor(0xC6, 0x28, 0x28)
LIGHT_RED = RGBColor(0xFF, 0xEB, 0xEE)
ORANGE = RGBColor(0xE6, 0x51, 0x00)
LIGHT_ORANGE = RGBColor(0xFF, 0xF3, 0xE0)
PURPLE = RGBColor(0x6A, 0x1B, 0x9A)
LIGHT_PURPLE = RGBColor(0xF3, 0xE5, 0xF5)
TEAL = RGBColor(0x00, 0x69, 0x5C)
GOLD = RGBColor(0xFF, 0xA0, 0x00)

# ══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

def add_bg(slide, color):
    """Set slide background color."""
    background = slide.background
    fill = background.fill
    fill.solid()
    fill.fore_color.rgb = color


def add_gradient_bg(slide, color1, color2):
    """Add a subtle gradient background."""
    background = slide.background
    fill = background.fill
    fill.gradient()
    fill.gradient_angle = 90
    fill.gradient_stops[0].color.rgb = color1
    fill.gradient_stops[1].color.rgb = color2


def add_text(slide, left, top, width, height, text, font_size=14,
             bold=False, italic=False, color=BLACK, align=PP_ALIGN.LEFT,
             font_name="Calibri", valign=MSO_ANCHOR.TOP):
    """Add a text box."""
    shape = slide.shapes.add_textbox(left, top, width, height)
    tf = shape.text_frame
    tf.word_wrap = True
    tf.auto_size = None
    p = tf.paragraphs[0]
    p.text = text
    p.font.size = Pt(font_size)
    p.font.bold = bold
    p.font.italic = italic
    p.font.color.rgb = color
    p.font.name = font_name
    p.alignment = align
    tf.vertical_anchor = valign
    return shape


def add_box(slide, left, top, width, height, fill_color, border_color=None, border_width=Pt(1)):
    """Add a rounded rectangle box."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE, left, top, width, height
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    if border_color:
        shape.line.color.rgb = border_color
        shape.line.width = border_width
    else:
        shape.line.fill.background()
    # Adjust corner radius
    shape.adjustments[0] = 0.1
    return shape


def add_header_bar(slide, title, subtitle=None):
    """Add a colored header bar at top of slide."""
    # Header background
    header = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(0), Inches(0), Inches(13.33), Inches(1.2)
    )
    header.fill.solid()
    header.fill.fore_color.rgb = DARK_BLUE
    header.line.fill.background()

    # Title
    add_text(slide, Inches(0.5), Inches(0.25), Inches(12), Inches(0.5),
             title, font_size=32, bold=True, color=WHITE)

    if subtitle:
        add_text(slide, Inches(0.5), Inches(0.75), Inches(12), Inches(0.35),
                 subtitle, font_size=16, color=RGBColor(0xBB, 0xDE, 0xFB))


def add_pyramid_shape(slide, left, top, width, height, fill_color, text, font_size=14):
    """Add a pyramid/triangle shape with text."""
    shape = slide.shapes.add_shape(
        MSO_SHAPE.ISOSCELES_TRIANGLE, left, top, width, height
    )
    shape.fill.solid()
    shape.fill.fore_color.rgb = fill_color
    shape.line.color.rgb = RGBColor(0x00, 0x00, 0x00)
    shape.line.width = Pt(1)

    # Add text on top
    text_top = top + height * 0.4
    add_text(slide, left, text_top, width, Inches(0.5),
             text, font_size=font_size, bold=True, color=WHITE, align=PP_ALIGN.CENTER)


def add_arrow(slide, start_left, start_top, end_left, end_top, color=DARK_GRAY):
    """Add an arrow connector."""
    # Use a simple line with arrow
    connector = slide.shapes.add_connector(
        1,  # straight connector
        start_left, start_top, end_left, end_top
    )
    connector.line.color.rgb = color
    connector.line.width = Pt(2)


def add_bullet_text(slide, left, top, width, height, items, font_size=14, color=DARK_GRAY):
    """Add bulleted text box."""
    shape = slide.shapes.add_textbox(left, top, width, height)
    tf = shape.text_frame
    tf.word_wrap = True

    for i, item in enumerate(items):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()
        p.text = f"• {item}"
        p.font.size = Pt(font_size)
        p.font.color.rgb = color
        p.font.name = "Calibri"
        p.space_after = Pt(6)

    return shape


# ══════════════════════════════════════════════════════════════════════════════
# CREATE PRESENTATION
# ══════════════════════════════════════════════════════════════════════════════

prs = Presentation()
prs.slide_width = Inches(13.33)
prs.slide_height = Inches(7.5)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 1: TITLE SLIDE
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_gradient_bg(slide, LIGHT_BLUE, WHITE)

# DNA helix decoration (simplified as colored bars)
for i in range(8):
    bar = slide.shapes.add_shape(
        MSO_SHAPE.RECTANGLE, Inches(0.3 + i * 0.15), Inches(2.5 + (i % 2) * 0.3),
        Inches(0.1), Inches(0.8)
    )
    bar.fill.solid()
    bar.fill.fore_color.rgb = MED_BLUE if i % 2 == 0 else TEAL
    bar.line.fill.background()
    bar.rotation = 15

# Main title
add_text(slide, Inches(1.5), Inches(2.2), Inches(10), Inches(1.2),
         "pLIN: Plasmid Life Identification Number",
         font_size=44, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)

add_text(slide, Inches(1.5), Inches(3.5), Inches(10), Inches(0.8),
         "A Hierarchical, Permanent Classification System for\nBacterial Plasmids and AMR Surveillance",
         font_size=24, color=MED_GRAY, align=PP_ALIGN.CENTER)

# Author info
add_text(slide, Inches(1.5), Inches(5.5), Inches(10), Inches(0.5),
         "Basil Xavier", font_size=20, bold=True, color=DARK_GRAY, align=PP_ALIGN.CENTER)

add_text(slide, Inches(1.5), Inches(6.0), Inches(10), Inches(0.5),
         "2025", font_size=16, color=MED_GRAY, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 2: THE PYRAMID PRINCIPLE OVERVIEW
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Executive Summary", "The Pyramid Principle")

# Main answer box at top (the pyramid tip)
answer_box = add_box(slide, Inches(3.5), Inches(1.5), Inches(6.33), Inches(1.2), DARK_BLUE)
add_text(slide, Inches(3.7), Inches(1.65), Inches(5.93), Inches(1.0),
         "pLIN provides a permanent, hierarchical classification\nsystem that solves plasmid tracking for AMR surveillance",
         font_size=18, bold=True, color=WHITE, align=PP_ALIGN.CENTER, valign=MSO_ANCHOR.MIDDLE)

# Three supporting arguments
support_data = [
    ("SITUATION", "AMR is a global crisis;\nplasmids drive 80% of\nresistance spread", LIGHT_RED, RED),
    ("COMPLICATION", "Current methods lack\npermanent IDs and\nmiss relationships", LIGHT_ORANGE, ORANGE),
    ("SOLUTION", "pLIN: 6-level hierarchy\nwith ML-powered\nclassification", LIGHT_GREEN, GREEN),
]

for i, (title, desc, fill, border) in enumerate(support_data):
    x = Inches(0.8 + i * 4.2)
    box = add_box(slide, x, Inches(3.2), Inches(3.8), Inches(1.8), fill, border, Pt(3))
    add_text(slide, x + Inches(0.15), Inches(3.35), Inches(3.5), Inches(0.4),
             title, font_size=16, bold=True, color=border, align=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.15), Inches(3.75), Inches(3.5), Inches(1.2),
             desc, font_size=14, color=DARK_GRAY, align=PP_ALIGN.CENTER)

# Key benefits row
add_text(slide, Inches(0.5), Inches(5.3), Inches(12), Inches(0.4),
         "Key Benefits", font_size=18, bold=True, color=DARK_BLUE)

benefits = [
    ("Permanent IDs", "Codes never change\nas database grows"),
    ("92% Accuracy", "KNN classifier (20 Inc groups)\nauto Inc detection"),
    ("Outbreak Detection", "Strain-level clustering\nfinds transmission"),
    ("AMR Integration", "Links plasmids to\nresistance genes"),
]

for i, (title, desc) in enumerate(benefits):
    x = Inches(0.5 + i * 3.2)
    box = add_box(slide, x, Inches(5.7), Inches(3.0), Inches(1.3), LIGHT_BLUE, MED_BLUE, Pt(2))
    add_text(slide, x + Inches(0.1), Inches(5.85), Inches(2.8), Inches(0.35),
             title, font_size=14, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.1), Inches(6.2), Inches(2.8), Inches(0.7),
             desc, font_size=12, color=DARK_GRAY, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 3: SITUATION - The AMR Crisis
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "SITUATION: The AMR Crisis", "Antimicrobial resistance threatens modern medicine")

# Key statistics
stats = [
    ("1.27M", "Deaths directly caused\nby AMR in 2019", RED),
    ("4.95M", "Deaths associated\nwith AMR in 2019", ORANGE),
    ("$100B+", "Annual global\neconomic burden", PURPLE),
    ("10M", "Projected annual deaths\nby 2050 if unchecked", RGBColor(0x4A, 0x14, 0x8C)),
]

add_text(slide, Inches(0.5), Inches(1.5), Inches(12), Inches(0.4),
         "The Scale of the Problem", font_size=20, bold=True, color=DARK_BLUE)

for i, (number, desc, color) in enumerate(stats):
    x = Inches(0.5 + i * 3.2)
    box = add_box(slide, x, Inches(2.0), Inches(3.0), Inches(1.6), color)
    add_text(slide, x, Inches(2.15), Inches(3.0), Inches(0.7),
             number, font_size=36, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.1), Inches(2.85), Inches(2.8), Inches(0.7),
             desc, font_size=12, color=WHITE, align=PP_ALIGN.CENTER)

# Why it matters section
add_text(slide, Inches(0.5), Inches(3.9), Inches(6), Inches(0.4),
         "Why Plasmids Matter", font_size=18, bold=True, color=DARK_BLUE)

plasmid_facts = [
    "Plasmids carry 60-80% of clinically relevant resistance genes",
    "Horizontal gene transfer spreads resistance across species barriers",
    "A single plasmid can carry multiple resistance genes (MDR)",
    "Conjugative plasmids spread resistance without selective pressure",
    "Plasmids persist in environments long after antibiotic use stops",
]
add_bullet_text(slide, Inches(0.5), Inches(4.3), Inches(6), Inches(2.5),
                plasmid_facts, font_size=14, color=DARK_GRAY)

# Visual: Plasmid transmission diagram
add_text(slide, Inches(7), Inches(3.9), Inches(6), Inches(0.4),
         "Plasmid-Mediated AMR Spread", font_size=18, bold=True, color=DARK_BLUE)

# Bacteria icons (simplified as circles)
bacteria_data = [
    (Inches(8), Inches(4.8), "Donor\nBacterium", GREEN),
    (Inches(11), Inches(4.8), "Recipient\nBacterium", MED_GRAY),
    (Inches(11), Inches(6.0), "Now\nResistant!", RED),
]

for x, y, label, color in bacteria_data:
    circle = slide.shapes.add_shape(MSO_SHAPE.OVAL, x, y, Inches(1.2), Inches(0.8))
    circle.fill.solid()
    circle.fill.fore_color.rgb = color
    circle.line.color.rgb = RGBColor(0, 0, 0)
    add_text(slide, x - Inches(0.1), y + Inches(0.15), Inches(1.4), Inches(0.6),
             label, font_size=10, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

# Plasmid (small circle) inside donor
plasmid = slide.shapes.add_shape(MSO_SHAPE.OVAL, Inches(8.4), Inches(5.0), Inches(0.4), Inches(0.4))
plasmid.fill.solid()
plasmid.fill.fore_color.rgb = RED
plasmid.line.color.rgb = RGBColor(0, 0, 0)

# Arrow indicating transfer
add_text(slide, Inches(9.3), Inches(4.95), Inches(1.5), Inches(0.4),
         "→ Conjugation →", font_size=12, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)

add_text(slide, Inches(10.7), Inches(5.6), Inches(1), Inches(0.3),
         "↓", font_size=20, bold=True, color=RED, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 4: SITUATION - Plasmid Biology Background
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "SITUATION: Plasmid Biology", "Understanding the vectors of resistance")

# Left column: What are plasmids?
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "What Are Plasmids?", font_size=20, bold=True, color=DARK_BLUE)

plasmid_info = [
    "Extrachromosomal circular DNA molecules (1 kb - 2 Mb)",
    "Self-replicating genetic elements in bacteria",
    "Carry accessory genes: resistance, virulence, metabolism",
    "Can transfer between bacteria via conjugation",
    "Multiple plasmids can coexist if compatible",
]
add_bullet_text(slide, Inches(0.5), Inches(1.9), Inches(6), Inches(2.2),
                plasmid_info, font_size=14, color=DARK_GRAY)

# Incompatibility groups
add_text(slide, Inches(0.5), Inches(4.2), Inches(6), Inches(0.4),
         "Incompatibility (Inc) Groups", font_size=18, bold=True, color=DARK_BLUE)

inc_info = [
    "Plasmids with same replication machinery cannot coexist",
    "Inc group = classification by replication incompatibility",
    "20+ Inc groups classified in pLIN (Enterobacteriaceae)",
    "Clinical importance: IncF, IncN, IncI, IncA/C, IncX, IncHI",
    "Same Inc group often indicates related plasmid lineages",
]
add_bullet_text(slide, Inches(0.5), Inches(4.6), Inches(6), Inches(2.2),
                inc_info, font_size=14, color=DARK_GRAY)

# Right column: Key plasmid types table
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Clinically Important Plasmid Types", font_size=18, bold=True, color=DARK_BLUE)

# Table header
table_header_box = add_box(slide, Inches(7), Inches(1.95), Inches(5.8), Inches(0.45), DARK_BLUE)
headers = ["Inc Group", "Size", "Key Resistance", "Mobility"]
for j, header in enumerate(headers):
    add_text(slide, Inches(7.1 + j * 1.45), Inches(2.0), Inches(1.4), Inches(0.35),
             header, font_size=11, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

# Table data
table_data = [
    ("IncF", "50-200 kb", "ESBLs, Carbapenemases", "Conjugative"),
    ("IncN", "40-60 kb", "ESBLs (CTX-M)", "Conjugative"),
    ("IncA/C", "100-200 kb", "MDR, Carbapenemases", "Conjugative"),
    ("IncX", "30-50 kb", "Carbapenemases (NDM)", "Conjugative"),
    ("IncHI", "200-300 kb", "Colistin (mcr)", "Conjugative"),
    ("IncI", "80-100 kb", "ESBLs, AmpC", "Conjugative"),
]

for i, (inc, size, resistance, mobility) in enumerate(table_data):
    y = Inches(2.45 + i * 0.5)
    fill = LIGHT_BLUE if i % 2 == 0 else WHITE
    row_box = add_box(slide, Inches(7), y, Inches(5.8), Inches(0.45), fill, MED_BLUE, Pt(0.5))
    row_data = [inc, size, resistance, mobility]
    for j, cell in enumerate(row_data):
        add_text(slide, Inches(7.1 + j * 1.45), y + Inches(0.05), Inches(1.4), Inches(0.35),
                 cell, font_size=10, color=DARK_GRAY, align=PP_ALIGN.CENTER)

# Key resistance genes
add_text(slide, Inches(7), Inches(5.6), Inches(6), Inches(0.35),
         "Plasmid-Borne Resistance Genes of Concern", font_size=14, bold=True, color=DARK_BLUE)

genes = [
    ("blaCTX-M", "ESBL"),
    ("blaKPC", "Carbapenemase"),
    ("blaNDM", "Carbapenemase"),
    ("mcr-1", "Colistin resistance"),
    ("qnr", "Quinolone resistance"),
]

for i, (gene, desc) in enumerate(genes):
    x = Inches(7 + (i % 3) * 2)
    y = Inches(6.0 + (i // 3) * 0.5)
    add_text(slide, x, y, Inches(1.8), Inches(0.4),
             f"{gene}: {desc}", font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 5: COMPLICATION - Current Methods Are Inadequate
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "COMPLICATION: Current Methods Fall Short",
               "Existing plasmid classification approaches have critical limitations")

# Current methods comparison
add_text(slide, Inches(0.5), Inches(1.5), Inches(12), Inches(0.4),
         "Existing Plasmid Classification Methods", font_size=20, bold=True, color=DARK_BLUE)

methods = [
    ("Inc Typing\n(PCR-based)",
     ["Well-established protocol", "Targets replicon genes", "Widely used clinically"],
     ["Only detects known replicons", "Misses mosaic plasmids", "No resolution within Inc groups", "Cannot track evolution"],
     LIGHT_BLUE, MED_BLUE),
    ("pMLST\n(Plasmid MLST)",
     ["Sequence-based typing", "Higher resolution than Inc", "Standardized schemes"],
     ["Only available for few Inc groups", "Requires known alleles", "No universal coverage", "Cannot compare across Inc groups"],
     LIGHT_GREEN, GREEN),
    ("MOB Typing\n(Mobility genes)",
     ["Classifies by transfer genes", "Predicts conjugation ability", "Useful for HGT studies"],
     ["Only for mobile plasmids", "Ignores non-conjugative", "Does not capture content", "Limited evolutionary signal"],
     LIGHT_PURPLE, PURPLE),
    ("Whole Plasmid\nAlignment",
     ["Uses complete sequence", "Captures all variation", "Gold standard accuracy"],
     ["Computationally expensive", "Sensitive to rearrangements", "No standardized IDs", "Hard to compare distantly related"],
     LIGHT_ORANGE, ORANGE),
]

for i, (name, pros, cons, fill, border) in enumerate(methods):
    x = Inches(0.4 + i * 3.2)
    # Method box
    box = add_box(slide, x, Inches(2.0), Inches(3.0), Inches(4.8), fill, border, Pt(2))
    add_text(slide, x + Inches(0.1), Inches(2.1), Inches(2.8), Inches(0.7),
             name, font_size=14, bold=True, color=border, align=PP_ALIGN.CENTER)

    # Pros
    add_text(slide, x + Inches(0.1), Inches(2.75), Inches(2.8), Inches(0.3),
             "✓ Strengths", font_size=11, bold=True, color=GREEN)
    for j, pro in enumerate(pros):
        add_text(slide, x + Inches(0.15), Inches(3.0 + j * 0.35), Inches(2.7), Inches(0.35),
                 f"• {pro}", font_size=10, color=DARK_GRAY)

    # Cons
    add_text(slide, x + Inches(0.1), Inches(4.1), Inches(2.8), Inches(0.3),
             "✗ Limitations", font_size=11, bold=True, color=RED)
    for j, con in enumerate(cons):
        add_text(slide, x + Inches(0.15), Inches(4.35 + j * 0.35), Inches(2.7), Inches(0.35),
                 f"• {con}", font_size=10, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 6: COMPLICATION - The Core Problems
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "COMPLICATION: Three Core Problems",
               "Why we need a new approach to plasmid classification")

problems = [
    ("No Permanent Identifiers",
     "Plasmid 'names' change between studies. The same plasmid may be called "
     "pKP-123, pCTX-M-15_NYC, or Plasmid_A depending on who sequences it. "
     "This makes cross-study comparison and outbreak tracking nearly impossible.",
     ["GenBank accessions are sample-specific, not plasmid-type-specific",
      "No universal naming convention exists",
      "Literature searches for 'the same plasmid' are unreliable",
      "Epidemiological links are missed due to naming inconsistency"],
     RED),
    ("Limited Hierarchical Resolution",
     "Inc groups are too coarse for epidemiology. Thousands of distinct plasmids "
     "share the same Inc type. We cannot distinguish outbreak clones from "
     "unrelated plasmids within an Inc group.",
     ["IncF contains >50% of Enterobacteriaceae plasmids",
      "pMLST only exists for ~10 Inc groups",
      "No way to relate plasmids across Inc groups",
      "Evolutionary relationships are invisible"],
     ORANGE),
    ("No AMR-Integrated View",
     "Plasmid classification and AMR gene detection are separate workflows. "
     "We know a patient has blaNDM, but not which plasmid carries it or how "
     "that plasmid relates to regional or global outbreaks.",
     ["AMRFinderPlus gives genes, not plasmid context",
      "ResFinder doesn't cluster related plasmids",
      "No tool links plasmid lineage to AMR cargo",
      "Transmission risk assessment is ad hoc"],
     PURPLE),
]

for i, (title, desc, points, color) in enumerate(problems):
    y = Inches(1.5 + i * 2.0)

    # Problem number
    num_box = add_box(slide, Inches(0.5), y, Inches(0.6), Inches(0.6), color)
    add_text(slide, Inches(0.5), y + Inches(0.1), Inches(0.6), Inches(0.5),
             str(i + 1), font_size=24, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Title
    add_text(slide, Inches(1.3), y, Inches(4), Inches(0.5),
             title, font_size=18, bold=True, color=color)

    # Description
    add_text(slide, Inches(1.3), y + Inches(0.5), Inches(5), Inches(1.4),
             desc, font_size=12, color=DARK_GRAY)

    # Bullet points on right
    for j, point in enumerate(points):
        add_text(slide, Inches(6.8), y + Inches(0.1 + j * 0.4), Inches(6), Inches(0.4),
                 f"• {point}", font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 7: QUESTION - What Do We Need?
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, LIGHT_BLUE)
add_header_bar(slide, "QUESTION: What Do We Need?",
               "Defining the requirements for effective plasmid surveillance")

# Central question
question_box = add_box(slide, Inches(2), Inches(1.6), Inches(9.33), Inches(1.2), WHITE, DARK_BLUE, Pt(3))
add_text(slide, Inches(2.2), Inches(1.75), Inches(8.93), Inches(1.0),
         "How can we create a universal, permanent, hierarchical classification\n"
         "system that enables global plasmid tracking and AMR surveillance?",
         font_size=20, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)

# Requirements
add_text(slide, Inches(0.5), Inches(3.1), Inches(12), Inches(0.4),
         "Essential Requirements", font_size=20, bold=True, color=DARK_BLUE)

requirements = [
    ("Permanent IDs", "Codes must never change\nwhen database expands", "🔒"),
    ("Hierarchical", "Multi-level resolution from\nfamily to strain", "📊"),
    ("Universal", "Works for ALL plasmids,\nnot just known Inc groups", "🌍"),
    ("Automated", "Minimal manual curation,\nscalable to millions", "⚡"),
    ("Integrated", "Links classification to\nAMR and mobility", "🔗"),
    ("Open", "Free, reproducible,\ncommunity-driven", "📖"),
]

for i, (title, desc, icon) in enumerate(requirements):
    x = Inches(0.5 + (i % 3) * 4.3)
    y = Inches(3.6 + (i // 3) * 1.8)

    box = add_box(slide, x, y, Inches(4.0), Inches(1.5), WHITE, MED_BLUE, Pt(2))
    add_text(slide, x + Inches(0.15), y + Inches(0.15), Inches(0.5), Inches(0.5),
             icon, font_size=24, align=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.7), y + Inches(0.2), Inches(3.1), Inches(0.4),
             title, font_size=16, bold=True, color=DARK_BLUE)
    add_text(slide, x + Inches(0.7), y + Inches(0.6), Inches(3.1), Inches(0.8),
             desc, font_size=12, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 8: ANSWER - Introducing pLIN
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "ANSWER: Introducing pLIN",
               "Plasmid Life Identification Number — A permanent, hierarchical classification")

# The pLIN code explanation
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "The pLIN Code Structure", font_size=20, bold=True, color=DARK_BLUE)

# Example pLIN code visual
code_box = add_box(slide, Inches(0.5), Inches(2.0), Inches(5.5), Inches(1.0), LIGHT_BLUE, MED_BLUE, Pt(2))
add_text(slide, Inches(0.7), Inches(2.15), Inches(5.1), Inches(0.7),
         "1 . 3 . 7 . 12 . 45 . 128",
         font_size=32, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)

# Level labels
levels = [
    ("A", "Family", "d ≤ 0.30"),
    ("B", "Subfamily", "d ≤ 0.20"),
    ("C", "Clade", "d ≤ 0.10"),
    ("D", "Subclade", "d ≤ 0.05"),
    ("E", "Lineage", "d ≤ 0.02"),
    ("F", "Strain", "d ≤ 0.01"),
]

for i, (pos, level, thresh) in enumerate(levels):
    x = Inches(0.65 + i * 0.9)
    add_text(slide, x, Inches(3.1), Inches(0.8), Inches(0.3),
             pos, font_size=14, bold=True, color=MED_BLUE, align=PP_ALIGN.CENTER)
    add_text(slide, x, Inches(3.4), Inches(0.8), Inches(0.3),
             level, font_size=11, color=DARK_GRAY, align=PP_ALIGN.CENTER)
    add_text(slide, x, Inches(3.65), Inches(0.8), Inches(0.3),
             thresh, font_size=9, color=MED_GRAY, align=PP_ALIGN.CENTER)

# How it works
add_text(slide, Inches(0.5), Inches(4.2), Inches(6), Inches(0.4),
         "How pLIN Works", font_size=18, bold=True, color=DARK_BLUE)

steps = [
    "1. Compute tetranucleotide (4-mer) frequency vector for each plasmid",
    "2. Calculate pairwise cosine distances between all plasmids",
    "3. Perform hierarchical clustering (single-linkage by default)",
    "4. Cut dendrogram at 6 thresholds to assign hierarchical codes",
    "5. Codes are permanent: new plasmids get new numbers, never reassign",
]
add_bullet_text(slide, Inches(0.5), Inches(4.55), Inches(6), Inches(2.5),
                [s[3:] for s in steps], font_size=13, color=DARK_GRAY)

# Right side: Key innovations
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Key Innovations", font_size=20, bold=True, color=DARK_BLUE)

innovations = [
    ("Alignment-Free", "4-mer composition captures plasmid 'signature'\nwithout computationally expensive alignment",
     GREEN),
    ("Permanent Codes", "Hierarchical clustering ensures codes never change\nwhen new plasmids are added to database",
     MED_BLUE),
    ("Inc-Aware", "KNN classifier auto-detects Inc group (92% accuracy, 20 groups)\nbefore pLIN assignment",
     PURPLE),
    ("AMR-Integrated", "Links plasmid lineage to AMRFinderPlus results\nfor comprehensive resistance profiling",
     ORANGE),
]

for i, (title, desc, color) in enumerate(innovations):
    y = Inches(1.95 + i * 1.35)
    box = add_box(slide, Inches(7), y, Inches(5.8), Inches(1.2), WHITE, color, Pt(2))
    add_text(slide, Inches(7.15), y + Inches(0.15), Inches(5.5), Inches(0.35),
             title, font_size=15, bold=True, color=color)
    add_text(slide, Inches(7.15), y + Inches(0.5), Inches(5.5), Inches(0.65),
             desc, font_size=12, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 9: ANSWER - The Complete Pipeline
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "ANSWER: The pLIN Pipeline",
               "From FASTA upload to epidemiological insights")

# Pipeline steps
pipeline = [
    ("1. Upload", "FASTA files", LIGHT_BLUE, MED_BLUE),
    ("2. Classify", "Inc group (KNN)", LIGHT_GREEN, GREEN),
    ("3. Vectorize", "4-mer frequencies", LIGHT_PURPLE, PURPLE),
    ("4. Cluster", "Hierarchical", LIGHT_ORANGE, ORANGE),
    ("5. Assign", "pLIN codes", LIGHT_BLUE, DARK_BLUE),
    ("6. Annotate", "AMR + Mobility", LIGHT_RED, RED),
    ("7. Analyze", "Outbreaks", LIGHT_GREEN, TEAL),
]

for i, (step, desc, fill, border) in enumerate(pipeline):
    x = Inches(0.3 + i * 1.85)
    box = add_box(slide, x, Inches(1.6), Inches(1.7), Inches(1.3), fill, border, Pt(2))
    add_text(slide, x + Inches(0.05), Inches(1.7), Inches(1.6), Inches(0.5),
             step, font_size=13, bold=True, color=border, align=PP_ALIGN.CENTER)
    add_text(slide, x + Inches(0.05), Inches(2.15), Inches(1.6), Inches(0.6),
             desc, font_size=11, color=DARK_GRAY, align=PP_ALIGN.CENTER)

    # Arrow to next
    if i < len(pipeline) - 1:
        add_text(slide, x + Inches(1.7), Inches(2.0), Inches(0.3), Inches(0.5),
                 "→", font_size=18, bold=True, color=DARK_GRAY)

# Features grid
add_text(slide, Inches(0.5), Inches(3.2), Inches(12), Inches(0.4),
         "Integrated Capabilities", font_size=18, bold=True, color=DARK_BLUE)

features = [
    ("Inc Group Detection", "KNN classifier trained on 6,998 plasmids\n92.2% accuracy across 20 Inc groups"),
    ("MOBsuite Typing", "3-tier mobility cascade\nRelaxase families + MPF types"),
    ("CRISPR Host Inference", "Spacer-based plasmid-host matching\nSoftmax probability ranking"),
    ("AMR Gene Detection", "AMRFinderPlus integration\nLinks resistance to plasmid lineages"),
    ("Outbreak Detection", "Strain-level clustering (F threshold)\nIdentifies transmission events"),
    ("DRAGNOME Buddy", "Local LLM-powered chatbot\nContext-aware Q&A about results"),
]

for i, (title, desc) in enumerate(features):
    x = Inches(0.5 + (i % 3) * 4.2)
    y = Inches(3.7 + (i // 3) * 1.6)
    box = add_box(slide, x, y, Inches(4.0), Inches(1.4), LIGHT_GRAY, MED_GRAY, Pt(1))
    add_text(slide, x + Inches(0.15), y + Inches(0.1), Inches(3.7), Inches(0.4),
             title, font_size=13, bold=True, color=DARK_BLUE)
    add_text(slide, x + Inches(0.15), y + Inches(0.5), Inches(3.7), Inches(0.85),
             desc, font_size=11, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 10: COSINE PAIRWISE DISTANCE EXPLAINED
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "The Mathematics: Cosine Pairwise Distance",
               "How pLIN measures plasmid similarity using k-mer composition")

# Left side: What is cosine distance?
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "What is Cosine Distance?", font_size=20, bold=True, color=DARK_BLUE)

add_text(slide, Inches(0.5), Inches(1.95), Inches(6), Inches(1.2),
         "Cosine distance measures the angle between two vectors in high-dimensional space. "
         "It captures how similar two plasmids are based on their k-mer composition, "
         "regardless of sequence length.",
         font_size=13, color=DARK_GRAY)

# Formula box
formula_box = add_box(slide, Inches(0.5), Inches(3.0), Inches(5.5), Inches(1.4), LIGHT_BLUE, MED_BLUE, Pt(2))
add_text(slide, Inches(0.7), Inches(3.1), Inches(5.1), Inches(0.4),
         "Cosine Similarity:", font_size=14, bold=True, color=DARK_BLUE)
add_text(slide, Inches(0.7), Inches(3.45), Inches(5.1), Inches(0.5),
         "cos(θ) = (A · B) / (||A|| × ||B||)", font_size=20, bold=True, color=DARK_BLUE, align=PP_ALIGN.CENTER)
add_text(slide, Inches(0.7), Inches(3.95), Inches(5.1), Inches(0.4),
         "Cosine Distance = 1 - cos(θ)    [Range: 0 to 1]", font_size=12, color=DARK_GRAY, align=PP_ALIGN.CENTER)

# Interpretation
add_text(slide, Inches(0.5), Inches(4.6), Inches(6), Inches(0.4),
         "Interpretation", font_size=16, bold=True, color=DARK_BLUE)

interpretations = [
    ("d = 0.00", "Identical k-mer profiles → Same plasmid", GREEN),
    ("d ≤ 0.01", "Nearly identical → Same strain (F level)", RGBColor(0x2E, 0x7D, 0x32)),
    ("d ≤ 0.05", "Very similar → Same subclade (D level)", MED_BLUE),
    ("d ≤ 0.10", "Similar → Same clade (C level)", PURPLE),
    ("d ≤ 0.30", "Related → Same family (A level)", ORANGE),
    ("d > 0.30", "Distant → Different plasmid families", RED),
]

for i, (dist, meaning, color) in enumerate(interpretations):
    y = Inches(5.0 + i * 0.38)
    add_text(slide, Inches(0.6), y, Inches(1.0), Inches(0.35),
             dist, font_size=11, bold=True, color=color)
    add_text(slide, Inches(1.7), y, Inches(4.5), Inches(0.35),
             meaning, font_size=11, color=DARK_GRAY)

# Right side: Visual explanation
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Step-by-Step Process", font_size=20, bold=True, color=DARK_BLUE)

# Step 1: K-mer vectors
step1_box = add_box(slide, Inches(7), Inches(1.95), Inches(5.8), Inches(1.3), LIGHT_GREEN, GREEN, Pt(2))
add_text(slide, Inches(7.15), Inches(2.05), Inches(5.5), Inches(0.35),
         "1. Convert Sequences to 4-mer Vectors", font_size=13, bold=True, color=GREEN)
add_text(slide, Inches(7.15), Inches(2.4), Inches(5.5), Inches(0.8),
         "Each plasmid → 256-dimensional vector (4⁴ possible 4-mers)\n"
         "Vector[i] = frequency of k-mer i in the sequence\n"
         "Example: ACGT appears 1,234 times → normalize by sequence length",
         font_size=11, color=DARK_GRAY)

# Step 2: Pairwise distances
step2_box = add_box(slide, Inches(7), Inches(3.4), Inches(5.8), Inches(1.3), LIGHT_PURPLE, PURPLE, Pt(2))
add_text(slide, Inches(7.15), Inches(3.5), Inches(5.5), Inches(0.35),
         "2. Calculate Pairwise Distances", font_size=13, bold=True, color=PURPLE)
add_text(slide, Inches(7.15), Inches(3.85), Inches(5.5), Inches(0.8),
         "Compare every plasmid to every other plasmid\n"
         "N plasmids → N×(N-1)/2 pairwise comparisons\n"
         "Result: Distance matrix (condensed form for clustering)",
         font_size=11, color=DARK_GRAY)

# Step 3: Hierarchical clustering
step3_box = add_box(slide, Inches(7), Inches(4.85), Inches(5.8), Inches(1.3), LIGHT_ORANGE, ORANGE, Pt(2))
add_text(slide, Inches(7.15), Inches(4.95), Inches(5.5), Inches(0.35),
         "3. Hierarchical Clustering (Single-Linkage)", font_size=13, bold=True, color=ORANGE)
add_text(slide, Inches(7.15), Inches(5.3), Inches(5.5), Inches(0.8),
         "Build dendrogram from distance matrix\n"
         "Single-linkage: cluster distance = minimum pairwise distance\n"
         "Cut at 6 thresholds → 6-level pLIN code (A.B.C.D.E.F)",
         font_size=11, color=DARK_GRAY)

# Why cosine? box at bottom
why_box = add_box(slide, Inches(7), Inches(6.3), Inches(5.8), Inches(0.9), RGBColor(0xE8, 0xF5, 0xE9), GREEN, Pt(1))
add_text(slide, Inches(7.15), Inches(6.4), Inches(5.5), Inches(0.35),
         "Why Cosine Distance?", font_size=12, bold=True, color=GREEN)
add_text(slide, Inches(7.15), Inches(6.7), Inches(5.5), Inches(0.45),
         "✓ Length-independent  ✓ Fast computation  ✓ Captures composition, not order",
         font_size=10, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 11: SUPPORTING ARGUMENT 1 - Permanent Hierarchical Codes
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Supporting Argument 1: Permanent Hierarchical Codes",
               "Codes never change — enabling reproducible, cross-study comparisons")

# Comparison: Before vs After
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "Before pLIN: Naming Chaos", font_size=18, bold=True, color=RED)

before_examples = [
    "• Same plasmid, different names across studies:",
    "   - pKPC-NYC-2023, pKPC_outbreak_3, Plasmid_A, unnamed",
    "• No way to know if they're related",
    "• Literature review is a nightmare",
    "• Outbreaks are missed or duplicated in reports",
]
for i, line in enumerate(before_examples):
    add_text(slide, Inches(0.5), Inches(1.9 + i * 0.4), Inches(6), Inches(0.4),
             line, font_size=13, color=DARK_GRAY)

add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "After pLIN: Permanent Identity", font_size=18, bold=True, color=GREEN)

after_examples = [
    "• Every plasmid gets a permanent pLIN:",
    "   - 1.3.7.12.45.128 (always this code, everywhere)",
    "• Hierarchy reveals relationships instantly",
    "• Cross-study comparison is trivial",
    "• Global surveillance becomes possible",
]
for i, line in enumerate(after_examples):
    add_text(slide, Inches(7), Inches(1.9 + i * 0.4), Inches(6), Inches(0.4),
             line, font_size=13, color=DARK_GRAY)

# Visual: Hierarchy interpretation
add_text(slide, Inches(0.5), Inches(4.0), Inches(12), Inches(0.4),
         "Interpreting pLIN Codes", font_size=18, bold=True, color=DARK_BLUE)

# Example comparison
examples = [
    ("1.3.7.12.45.128", "1.3.7.12.45.129", "Same lineage (E), different strains — likely recent divergence"),
    ("1.3.7.12.45.128", "1.3.7.12.46.200", "Same subclade (D), different lineages — common ancestor"),
    ("1.3.7.12.45.128", "1.3.8.15.50.130", "Same subfamily (B), different clades — distant relatives"),
    ("1.3.7.12.45.128", "2.5.10.20.60.150", "Different families — unrelated plasmids"),
]

for i, (code1, code2, interpretation) in enumerate(examples):
    y = Inches(4.5 + i * 0.7)
    add_text(slide, Inches(0.5), y, Inches(2.5), Inches(0.4),
             code1, font_size=12, bold=True, color=MED_BLUE)
    add_text(slide, Inches(3.2), y, Inches(0.5), Inches(0.4),
             "vs", font_size=12, color=DARK_GRAY, align=PP_ALIGN.CENTER)
    add_text(slide, Inches(3.8), y, Inches(2.5), Inches(0.4),
             code2, font_size=12, bold=True, color=MED_BLUE)
    add_text(slide, Inches(6.5), y, Inches(6.5), Inches(0.4),
             interpretation, font_size=12, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 11: SUPPORTING ARGUMENT 2 - ML-Powered Classification
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Supporting Argument 2: ML-Powered Classification",
               "92% accurate Inc group detection across 20 groups with confidence scoring")

# KNN Classifier
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "KNN Classifier for Inc Groups", font_size=18, bold=True, color=DARK_BLUE)

knn_features = [
    "Trained on 6,998 characterized plasmids across 20 Inc groups",
    "Uses 4-mer frequency vectors (256 features)",
    "k=5 neighbors with distance-weighted voting (cosine metric)",
    "Achieves 92.2% stratified cross-validated accuracy",
    "Provides confidence scores for each prediction",
]
add_bullet_text(slide, Inches(0.5), Inches(1.9), Inches(5.5), Inches(2.0),
                knn_features, font_size=13, color=DARK_GRAY)

# Unknown detection
add_text(slide, Inches(0.5), Inches(4.0), Inches(6), Inches(0.4),
         "Unknown/Novel Detection", font_size=18, bold=True, color=DARK_BLUE)

unknown_features = [
    "Confidence threshold: 40%",
    "Low-confidence predictions flagged as 'Unknown/Novel'",
    "Shows top 5 candidate Inc groups with probabilities",
    "Prevents forced misclassification of novel plasmids",
    "Researchers can investigate further with full candidate list",
]
add_bullet_text(slide, Inches(0.5), Inches(4.4), Inches(5.5), Inches(2.0),
                unknown_features, font_size=13, color=DARK_GRAY)

# Right side: Confidence example
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Example: Confidence Scoring", font_size=18, bold=True, color=DARK_BLUE)

# High confidence example
high_conf_box = add_box(slide, Inches(7), Inches(1.95), Inches(5.5), Inches(1.8), LIGHT_GREEN, GREEN, Pt(2))
add_text(slide, Inches(7.15), Inches(2.05), Inches(5.2), Inches(0.35),
         "High Confidence (>40%)", font_size=14, bold=True, color=GREEN)
add_text(slide, Inches(7.15), Inches(2.4), Inches(5.2), Inches(1.3),
         "Predicted: IncFII (78%)\n"
         "Top 5: IncFII (78%), IncFIA (12%), IncFIB (5%), IncN (3%), IncI (2%)\n\n"
         "→ Confident classification as IncFII",
         font_size=12, color=DARK_GRAY)

# Low confidence example
low_conf_box = add_box(slide, Inches(7), Inches(3.9), Inches(5.5), Inches(1.8), LIGHT_ORANGE, ORANGE, Pt(2))
add_text(slide, Inches(7.15), Inches(4.0), Inches(5.2), Inches(0.35),
         "Low Confidence (<40%)", font_size=14, bold=True, color=ORANGE)
add_text(slide, Inches(7.15), Inches(4.35), Inches(5.2), Inches(1.3),
         "Predicted: Unknown/Novel (32%)\n"
         "Top 5: IncX3 (32%), IncX4 (25%), IncL (18%), IncN (15%), IncP (10%)\n\n"
         "→ Flagged for manual review, likely novel plasmid type",
         font_size=12, color=DARK_GRAY)

# Optional: Nucleotide Transformer
add_text(slide, Inches(7), Inches(5.9), Inches(5.5), Inches(0.35),
         "Optional: Nucleotide Transformer LLM", font_size=14, bold=True, color=PURPLE)
add_text(slide, Inches(7), Inches(6.25), Inches(5.5), Inches(0.8),
         "Deep learning embeddings for enhanced prediction.\n"
         "Captures long-range sequence patterns beyond k-mers.",
         font_size=12, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 12: SUPPORTING ARGUMENT 3 - AMR Integration
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Supporting Argument 3: AMR Integration",
               "Linking plasmid lineages to resistance gene cargo")

# The problem
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "The Missing Link in AMR Surveillance", font_size=18, bold=True, color=RED)

add_text(slide, Inches(0.5), Inches(1.9), Inches(6), Inches(1.0),
         "Traditional AMR analysis tells us WHAT genes are present, but not:\n"
         "• Which plasmid carries them\n"
         "• How that plasmid relates to others\n"
         "• Whether the same plasmid is spreading regionally/globally",
         font_size=13, color=DARK_GRAY)

# The solution
add_text(slide, Inches(0.5), Inches(3.1), Inches(6), Inches(0.4),
         "pLIN Closes the Loop", font_size=18, bold=True, color=GREEN)

solution_points = [
    "AMRFinderPlus integration detects resistance genes",
    "Prodigal annotation identifies all ORFs",
    "pLIN codes link AMR cargo to plasmid lineage",
    "Epidemiology tab identifies high-risk plasmids",
    "Mobility prediction assesses transmission potential",
]
add_bullet_text(slide, Inches(0.5), Inches(3.5), Inches(6), Inches(2.0),
                solution_points, font_size=13, color=DARK_GRAY)

# Right side: Example output
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Integrated Analysis Example", font_size=18, bold=True, color=DARK_BLUE)

# Example table
table_box = add_box(slide, Inches(7), Inches(1.95), Inches(5.8), Inches(2.5), WHITE, MED_BLUE, Pt(1))

example_data = [
    ("pLIN", "Inc", "AMR Genes", "Mobility", "Risk"),
    ("1.3.7.12.45.128", "IncFII", "blaNDM-1, aac(6')", "Conjugative", "HIGH"),
    ("1.3.7.12.45.129", "IncFII", "blaNDM-1, aac(6')", "Conjugative", "HIGH"),
    ("2.5.10.20.60.150", "IncN", "blaCTX-M-15", "Conjugative", "MED"),
    ("3.8.15.30.90.300", "IncX4", "mcr-1", "Mobilizable", "HIGH"),
]

for i, row in enumerate(example_data):
    y = Inches(2.0 + i * 0.45)
    fill = DARK_BLUE if i == 0 else (LIGHT_BLUE if i % 2 == 1 else WHITE)
    text_color = WHITE if i == 0 else DARK_GRAY

    widths = [1.8, 0.8, 1.5, 1.0, 0.6]
    x_offset = 7.05
    for j, (cell, w) in enumerate(zip(row, widths)):
        add_text(slide, Inches(x_offset), y, Inches(w), Inches(0.4),
                 cell, font_size=10 if i > 0 else 11,
                 bold=(i == 0), color=text_color)
        x_offset += w

# Insight callout
insight_box = add_box(slide, Inches(7), Inches(4.6), Inches(5.8), Inches(1.2), LIGHT_GREEN, GREEN, Pt(2))
add_text(slide, Inches(7.15), Inches(4.7), Inches(5.5), Inches(0.35),
         "Insight", font_size=14, bold=True, color=GREEN)
add_text(slide, Inches(7.15), Inches(5.05), Inches(5.5), Inches(0.7),
         "Plasmids 1.3.7.12.45.128 and 1.3.7.12.45.129 share the same lineage (E)\n"
         "and carry identical AMR genes → Likely outbreak transmission!",
         font_size=12, color=DARK_GRAY)

# Bottom: Epidemiological value
add_text(slide, Inches(0.5), Inches(5.8), Inches(12), Inches(0.4),
         "Epidemiological Value", font_size=16, bold=True, color=DARK_BLUE)
add_text(slide, Inches(0.5), Inches(6.15), Inches(12), Inches(0.8),
         "pLIN enables plasmid-centric surveillance: track specific resistance-carrying plasmid lineages across institutions, "
         "regions, and countries. Identify emerging threats before they become widespread outbreaks.",
         font_size=13, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 13: IMPLEMENTATION - GUI & Deployment
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Implementation: Easy-to-Use GUI",
               "No command line required — web-based interface for all users")

# Streamlit GUI features
add_text(slide, Inches(0.5), Inches(1.5), Inches(6), Inches(0.4),
         "8-Tab Streamlit Interface", font_size=18, bold=True, color=DARK_BLUE)

tabs = [
    ("📋 Overview", "pLIN hierarchy explanation and method description"),
    ("📊 Results", "Interactive table with filtering, sorting, search"),
    ("🌳 Cladogram", "Hierarchical tree visualization with download"),
    ("💊 AMR Analysis", "Resistance gene details and heatmaps"),
    ("🔬 Epidemiology", "Outbreak detection and risk assessment"),
    ("🧫 CRISPR Host", "Spacer-based plasmid-host inference"),
    ("🧬 DRAGNOME Buddy", "AI chatbot for Q&A about results"),
    ("📥 Export", "Download TSV, PNG, PDF, ZIP bundle"),
]

for i, (tab, desc) in enumerate(tabs):
    y = Inches(1.9 + i * 0.48)
    add_text(slide, Inches(0.5), y, Inches(1.8), Inches(0.4),
             tab, font_size=13, bold=True, color=MED_BLUE)
    add_text(slide, Inches(2.3), y, Inches(4), Inches(0.4),
             desc, font_size=12, color=DARK_GRAY)

# Deployment options
add_text(slide, Inches(7), Inches(1.5), Inches(6), Inches(0.4),
         "Deployment Options", font_size=18, bold=True, color=DARK_BLUE)

deployments = [
    ("Local", "pip install + streamlit run\nRuns on your machine", LIGHT_BLUE, MED_BLUE),
    ("Docker", "docker-compose up\nContainerized deployment", LIGHT_GREEN, GREEN),
    ("Cloud", "Streamlit Cloud (free)\nOne-click deployment", LIGHT_PURPLE, PURPLE),
]

for i, (title, desc, fill, border) in enumerate(deployments):
    y = Inches(1.95 + i * 1.4)
    box = add_box(slide, Inches(7), y, Inches(5.5), Inches(1.2), fill, border, Pt(2))
    add_text(slide, Inches(7.15), y + Inches(0.15), Inches(5.2), Inches(0.35),
             title, font_size=14, bold=True, color=border)
    add_text(slide, Inches(7.15), y + Inches(0.5), Inches(5.2), Inches(0.6),
             desc, font_size=12, color=DARK_GRAY)

# Quick start
add_text(slide, Inches(0.5), Inches(5.8), Inches(12), Inches(0.4),
         "Quick Start", font_size=16, bold=True, color=DARK_BLUE)

code_box = add_box(slide, Inches(0.5), Inches(6.15), Inches(12.3), Inches(0.9), RGBColor(0x26, 0x32, 0x38))
add_text(slide, Inches(0.7), Inches(6.25), Inches(12), Inches(0.7),
         "git clone https://github.com/your-repo/pLIN.git && cd pLIN\n"
         "pip install -r requirements.txt\n"
         "streamlit run plin_app.py",
         font_size=12, color=RGBColor(0x4C, 0xAF, 0x50))


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 14: New Capabilities — Host Inference & Enhanced Mobility (NEW)
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "New Capabilities: Host Inference & Enhanced Mobility",
               "Transforming pLIN from classification tool to comprehensive plasmid epidemiology platform")

# Left: CRISPR
add_box(slide, Inches(0.4), Inches(1.5), Inches(6.2), Inches(4.8), WHITE,
        border_color=TEAL, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(1.6), Inches(5.8), Inches(0.4),
         "CRISPR Spacer-Based Host Inference", font_size=16, bold=True, color=TEAL)

crispr_items = [
    "MinCED extracts CRISPR spacers from host bacterial genomes",
    "BLASTN-short matches spacers against plasmid sequences",
    "Stringent filtering: \u226595% identity, \u226525bp, \u22641 mismatch",
    "Softmax probability ranking per host-plasmid pair",
    "Confidence categories: High (\u22650.7), Medium (\u22650.4), Low",
    "Two source modes: uploaded plasmids or reference DB (72,556)",
    "Heatmap, pie chart, and bar chart visualizations",
    "Full export: predictions, spacers, and summary TSVs",
]
add_bullet_text(slide, Inches(0.6), Inches(2.1), Inches(5.8), Inches(4.0),
                crispr_items, font_size=12, color=DARK_GRAY)

# Right: MOBsuite
add_box(slide, Inches(6.9), Inches(1.5), Inches(6.0), Inches(4.8), WHITE,
        border_color=ORANGE, border_width=Pt(2))
add_text(slide, Inches(7.1), Inches(1.6), Inches(5.6), Inches(0.4),
         "MOBsuite 3-Tier Mobility Cascade", font_size=16, bold=True, color=ORANGE)

mob_items = [
    "Priority 1: MOBsuite mob_typer (when available)",
    "  \u2192 Relaxase families: MOBF, MOBH, MOBP, MOBQ, MOBC, MOBV",
    "  \u2192 MPF types: Type T, F, I, G",
    "Priority 2: AMRFinderPlus gene scan",
    "  \u2192 Conjugative: tra/trb + virB1-11, trwA-N, pilX, taxC",
    "  \u2192 Mobilizable: mob, nikC-E, MOBF-V, oriT",
    "Priority 3: Non-mobilizable (default)",
    "Risk stratification integrated into Epidemiology tab",
]
add_bullet_text(slide, Inches(7.1), Inches(2.1), Inches(5.6), Inches(4.0),
                mob_items, font_size=12, color=DARK_GRAY)

# Bottom banner
add_box(slide, Inches(0.4), Inches(6.5), Inches(12.5), Inches(0.7),
        RGBColor(0xE3, 0xF2, 0xFD), border_color=DARK_BLUE, border_width=Pt(2))
add_text(slide, Inches(0.6), Inches(6.55), Inches(12.1), Inches(0.6),
         "Together, these additions enable researchers to answer: Which host carries this plasmid? "
         "Can this plasmid transfer? What resistance does it carry? \u2014 All in a single analysis run.",
         font_size=12, bold=True, color=DARK_BLUE)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 15: CONCLUSION - Call to Action
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "Conclusion: A New Era in Plasmid Surveillance",
               "pLIN enables global, standardized plasmid tracking")

# Summary pyramid
add_text(slide, Inches(0.5), Inches(1.5), Inches(5.5), Inches(0.4),
         "Summary", font_size=20, bold=True, color=DARK_BLUE)

summary_box = add_box(slide, Inches(0.5), Inches(1.95), Inches(5.5), Inches(3.0), LIGHT_BLUE, MED_BLUE, Pt(2))

summary_points = [
    ("SITUATION", "AMR is a global crisis; plasmids are the\nprimary vectors of resistance spread"),
    ("COMPLICATION", "Current methods lack permanent IDs,\nhierarchical resolution, and AMR integration"),
    ("SOLUTION", "pLIN provides permanent, hierarchical codes\nwith ML classification and AMR integration"),
]

for i, (label, text) in enumerate(summary_points):
    y = Inches(2.1 + i * 0.95)
    add_text(slide, Inches(0.65), y, Inches(1.3), Inches(0.35),
             label, font_size=12, bold=True, color=DARK_BLUE)
    add_text(slide, Inches(2.0), y, Inches(3.9), Inches(0.85),
             text, font_size=12, color=DARK_GRAY)

# Key takeaways
add_text(slide, Inches(0.5), Inches(5.1), Inches(5.5), Inches(0.4),
         "Key Takeaways", font_size=16, bold=True, color=DARK_BLUE)

takeaways = [
    "Permanent, universal plasmid identification",
    "92% accurate Inc group classification (20 groups)",
    "Integrated AMR and MOBsuite-enhanced mobility analysis",
    "CRISPR-based plasmid-host inference",
    "Outbreak detection at strain level",
    "Easy-to-use 8-tab GUI, no bioinformatics expertise needed",
]
add_bullet_text(slide, Inches(0.5), Inches(5.45), Inches(5.5), Inches(1.8),
                takeaways, font_size=12, color=DARK_GRAY)

# Call to action
add_text(slide, Inches(7), Inches(1.5), Inches(5.8), Inches(0.4),
         "Next Steps", font_size=20, bold=True, color=DARK_BLUE)

cta_items = [
    ("Try pLIN", "Upload your plasmid sequences today\nhttps://github.com/your-repo/pLIN", GREEN),
    ("Contribute", "Report issues, suggest features,\nsubmit training data", MED_BLUE),
    ("Cite", "If you use pLIN, please cite our work\n(see documentation)", PURPLE),
    ("Collaborate", "Contact us for integration into\nyour surveillance pipeline", ORANGE),
]

for i, (title, desc, color) in enumerate(cta_items):
    y = Inches(1.95 + i * 1.25)
    box = add_box(slide, Inches(7), y, Inches(5.5), Inches(1.1), WHITE, color, Pt(2))
    add_text(slide, Inches(7.15), y + Inches(0.1), Inches(5.2), Inches(0.35),
             title, font_size=14, bold=True, color=color)
    add_text(slide, Inches(7.15), y + Inches(0.45), Inches(5.2), Inches(0.6),
             desc, font_size=11, color=DARK_GRAY)

# Contact
add_text(slide, Inches(7), Inches(6.3), Inches(5.5), Inches(0.4),
         "Contact: basil.xavier@example.com", font_size=12, color=MED_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SLIDE 15: REFERENCES
# ══════════════════════════════════════════════════════════════════════════════

slide = prs.slides.add_slide(prs.slide_layouts[6])
add_bg(slide, WHITE)
add_header_bar(slide, "References & Resources", "Key literature and tools")

references = [
    "Murray CJ et al. (2022). Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. Lancet.",
    "Carattoli A et al. (2014). In silico detection and typing of plasmids using PlasmidFinder. Antimicrob Agents Chemother.",
    "Jolley KA & Maiden MC (2010). BIGSdb: Scalable analysis of bacterial genome variation. BMC Bioinformatics.",
    "Robertson J & Nash JHE (2018). MOB-suite: software tools for clustering, reconstruction and typing of plasmids. Microb Genom.",
    "Feldgarden M et al. (2021). AMRFinderPlus and the Reference Gene Catalog. Sci Rep.",
    "Dalla-Costa LM et al. (2023). LIN codes for genomic epidemiology of bacterial pathogens. Nat Microbiol.",
    "Hyatt D et al. (2010). Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics.",
]

add_text(slide, Inches(0.5), Inches(1.5), Inches(12), Inches(0.4),
         "Key References", font_size=18, bold=True, color=DARK_BLUE)

for i, ref in enumerate(references):
    add_text(slide, Inches(0.5), Inches(1.95 + i * 0.5), Inches(12), Inches(0.5),
             f"{i+1}. {ref}", font_size=11, color=DARK_GRAY)

# Tools used
add_text(slide, Inches(0.5), Inches(5.6), Inches(12), Inches(0.4),
         "Tools & Technologies", font_size=16, bold=True, color=DARK_BLUE)

tools = [
    ("Streamlit", "Web interface"),
    ("scikit-learn", "KNN classifier"),
    ("SciPy", "Hierarchical clustering"),
    ("AMRFinderPlus", "AMR detection"),
    ("Prodigal", "Gene prediction"),
    ("Ollama", "Local LLM"),
]

for i, (tool, desc) in enumerate(tools):
    x = Inches(0.5 + (i % 3) * 4.2)
    y = Inches(6.0 + (i // 3) * 0.5)
    add_text(slide, x, y, Inches(4), Inches(0.4),
             f"• {tool}: {desc}", font_size=12, color=DARK_GRAY)


# ══════════════════════════════════════════════════════════════════════════════
# SAVE PRESENTATION
# ══════════════════════════════════════════════════════════════════════════════

output_path = "pLIN_Pyramid_Presentation.pptx"
prs.save(output_path)
print(f"Presentation saved to: {output_path}")
print(f"Total slides: {len(prs.slides)}")
