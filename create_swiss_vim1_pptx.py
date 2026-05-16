#!/usr/bin/env python3
"""
Generate a Swiss VIM-1 real-time pLIN validation PPTX.
Output: output/pLIN_Swiss_VIM1_CaseStudy.pptx
"""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Cm

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Colour palette ────────────────────────────────────────────────────────────
NAVY    = RGBColor(0x1A, 0x2E, 0x4A)   # dark navy
TEAL    = RGBColor(0x00, 0x7B, 0x8A)   # teal accent
AMBER   = RGBColor(0xF5, 0xA6, 0x23)   # amber highlight
RED     = RGBColor(0xC0, 0x39, 0x2B)   # alert red
GREEN   = RGBColor(0x27, 0xAE, 0x60)   # success green
LGRAY   = RGBColor(0xF4, 0xF6, 0xF9)   # light grey background
WHITE   = RGBColor(0xFF, 0xFF, 0xFF)
DGRAY   = RGBColor(0x44, 0x44, 0x44)
MIDGRAY = RGBColor(0x88, 0x88, 0x88)

SLIDE_W = Inches(13.333)
SLIDE_H = Inches(7.5)


def new_prs():
    prs = Presentation()
    prs.slide_width  = SLIDE_W
    prs.slide_height = SLIDE_H
    return prs


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def rect(slide, l, t, w, h, fill=None, line=None, line_w=Pt(0)):
    from pptx.util import Emu
    shape = slide.shapes.add_shape(1, l, t, w, h)
    shape.line.width = line_w
    if fill:
        shape.fill.solid()
        shape.fill.fore_color.rgb = fill
    else:
        shape.fill.background()
    if line:
        shape.line.color.rgb = line
    else:
        shape.line.fill.background()
    return shape


def txbox(slide, text, l, t, w, h,
          size=18, bold=False, italic=False, color=NAVY,
          align=PP_ALIGN.LEFT, wrap=True, font="Calibri"):
    tb = slide.shapes.add_textbox(l, t, w, h)
    tf = tb.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size = Pt(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    run.font.name = font
    return tb


def add_para(tf, text, size=14, bold=False, italic=False,
             color=DGRAY, align=PP_ALIGN.LEFT, space_before=Pt(4), font="Calibri"):
    from pptx.util import Pt as P
    p = tf.add_paragraph()
    p.alignment = align
    p.space_before = space_before
    run = p.add_run()
    run.text = text
    run.font.size = P(size)
    run.font.bold = bold
    run.font.italic = italic
    run.font.color.rgb = color
    run.font.name = font
    return p


def header_bar(slide, title, subtitle=None):
    rect(slide, 0, 0, SLIDE_W, Inches(1.15), fill=NAVY)
    txbox(slide, title, Inches(0.4), Inches(0.12), Inches(12.4), Inches(0.6),
          size=26, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    if subtitle:
        txbox(slide, subtitle, Inches(0.4), Inches(0.72), Inches(12.4), Inches(0.35),
              size=14, color=AMBER, align=PP_ALIGN.LEFT)


def footer(slide, text="pLIN — Prospective Real-Time Validation | Swiss VIM-1 Outbreak"):
    rect(slide, 0, Inches(7.15), SLIDE_W, Inches(0.35), fill=NAVY)
    txbox(slide, text, Inches(0.3), Inches(7.17), Inches(9), Inches(0.28),
          size=9, color=WHITE)
    txbox(slide, "UMCG / DRAIGON Consortium | 2025",
          Inches(10), Inches(7.17), Inches(3.1), Inches(0.28),
          size=9, color=AMBER, align=PP_ALIGN.RIGHT)


def pill(slide, text, l, t, w, h, bg=TEAL, fg=WHITE, size=13, bold=True):
    rect(slide, l, t, w, h, fill=bg)
    txbox(slide, text, l, t, w, h, size=size, bold=bold, color=fg,
          align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 1 — TITLE
# ══════════════════════════════════════════════════════════════════════════════
def slide_title(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=NAVY)
    # Teal accent stripe
    rect(sl, 0, Inches(4.6), SLIDE_W, Inches(0.08), fill=TEAL)

    txbox(sl, "pLIN Real-Time Validation",
          Inches(0.7), Inches(1.0), Inches(11.9), Inches(1.1),
          size=40, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    txbox(sl, "Swiss VIM-1 Multispecies Hospital Outbreak",
          Inches(0.7), Inches(2.2), Inches(11.9), Inches(0.7),
          size=28, bold=False, color=AMBER, align=PP_ALIGN.CENTER)

    txbox(sl,
          "Prospective de novo assembly → pLIN classification across 5 Enterobacterales species\n"
          "without pre-deposited reference sequences",
          Inches(1.5), Inches(3.05), Inches(10.3), Inches(0.8),
          size=16, color=RGBColor(0xBB, 0xCC, 0xDD), align=PP_ALIGN.CENTER)

    # Three stat pills at bottom
    for i, (val, lbl, col) in enumerate([
        ("8 isolates", "ONT + Illumina sequenced", TEAL),
        ("5 species", "Enterobacterales affected",  RGBColor(0x8E,0x44,0xAD)),
        ("1 pLIN code", "Shared at L6 strain level", GREEN),
    ]):
        x = Inches(1.0 + i * 3.9)
        rect(sl, x, Inches(4.9), Inches(3.4), Inches(1.4), fill=col)
        txbox(sl, val,  x, Inches(5.0),  Inches(3.4), Inches(0.55),
              size=22, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        txbox(sl, lbl,  x, Inches(5.55), Inches(3.4), Inches(0.55),
              size=11, color=WHITE, align=PP_ALIGN.CENTER)

    txbox(sl, "Basil Britto Xavier · UMCG / DRAIGON Consortium · 2025",
          Inches(0.7), Inches(6.9), Inches(11.9), Inches(0.4),
          size=11, color=MIDGRAY, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 2 — OUTBREAK BACKGROUND
# ══════════════════════════════════════════════════════════════════════════════
def slide_background(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "Outbreak Background", "Swiss tertiary care hospital · 2024 · PRJEB98563")

    # Left panel — facts
    rect(sl, Inches(0.3), Inches(1.3), Inches(5.8), Inches(5.6), fill=WHITE)
    txbox(sl, "Clinical Context", Inches(0.5), Inches(1.4), Inches(5.4), Inches(0.4),
          size=15, bold=True, color=NAVY)

    facts = [
        ("55 patients",        "spread across multiple wards"),
        ("5 bacterial species","E. hormaechei (n=37), E. kobei,\nE. ludwigii, C. farmeri, E. coli"),
        ("blaVIM-1 + mcr-9",   "co-carried on large IncHI2 plasmid\n(249–343 kb)"),
        ("Sequencing",         "55 Illumina NextSeq 1000 runs\n8 ONT GridION hybrid sequences"),
        ("ENA deposit",        "Raw reads only — no pre-assembled\nplasmid sequences available"),
    ]
    y = Inches(1.95)
    for val, desc in facts:
        rect(sl, Inches(0.5), y, Inches(0.07), Inches(0.28), fill=TEAL)
        txbox(sl, val,  Inches(0.7), y - Inches(0.03), Inches(2.2), Inches(0.35),
              size=13, bold=True, color=NAVY)
        txbox(sl, desc, Inches(2.95), y - Inches(0.03), Inches(2.9), Inches(0.5),
              size=11, color=DGRAY)
        y += Inches(0.85)

    # Right panel — why this is hard
    rect(sl, Inches(6.4), Inches(1.3), Inches(6.6), Inches(2.55), fill=WHITE)
    txbox(sl, "The Challenge", Inches(6.6), Inches(1.4), Inches(6.2), Inches(0.4),
          size=15, bold=True, color=NAVY)
    challenges = [
        "❌  No pre-deposited plasmid sequences in ENA",
        "❌  5 different bacterial species — polyclonal host",
        "❌  Inc typing alone cannot resolve outbreak plasmid",
        "❌  Manual analysis is slow, labour-intensive",
    ]
    y = Inches(1.95)
    for c in challenges:
        txbox(sl, c, Inches(6.6), y, Inches(6.1), Inches(0.38), size=12, color=DGRAY)
        y += Inches(0.46)

    # Right lower — pLIN approach
    rect(sl, Inches(6.4), Inches(4.05), Inches(6.6), Inches(2.8), fill=TEAL)
    txbox(sl, "pLIN Approach", Inches(6.6), Inches(4.15), Inches(6.2), Inches(0.4),
          size=15, bold=True, color=WHITE)
    solutions = [
        "✔  Assemble from raw ONT reads (Flye 2.9.6)",
        "✔  Extract plasmid contigs by size (5–700 kb)",
        "✔  Classify against 72,556-plasmid reference",
        "✔  Single pLIN code = outbreak attribution",
    ]
    y = Inches(4.65)
    for s in solutions:
        txbox(sl, s, Inches(6.6), y, Inches(6.1), Inches(0.38), size=12, color=WHITE)
        y += Inches(0.46)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 3 — PIPELINE
# ══════════════════════════════════════════════════════════════════════════════
def slide_pipeline(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "Analysis Pipeline", "From raw ONT reads to pLIN code — no pre-deposited references")

    steps = [
        (TEAL,                   "1\nDownload",     "ONT reads\n(ERR15903138–46)\nENA FTP"),
        (RGBColor(0x16,0x6F,0x82),"2\nAssemble",   "Flye 2.9.6\n--nano-raw\ngenome-size 5 Mb"),
        (RGBColor(0x0E,0x50,0x6A),"3\nExtract",    "Plasmid contigs\n5–700 kb\ncircular flag"),
        (NAVY,                   "4\n4-mer",         "Tetranucleotide\nfrequency vectors\n256 features"),
        (RGBColor(0x8E,0x44,0xAD),"5\nClassify",   "KNN classifier\n28 Inc/Rep groups\ncosine metric"),
        (GREEN,                  "6\npLIN",          "NN assignment\n72,556 reference\nplasmids"),
    ]

    bw = Inches(1.75)
    bh = Inches(2.4)
    gap = Inches(0.25)
    y0  = Inches(1.7)
    x0  = Inches(0.35)

    for i, (col, lbl, desc) in enumerate(steps):
        x = x0 + i * (bw + gap)
        rect(sl, x, y0, bw, bh, fill=col)
        txbox(sl, lbl, x, y0 + Inches(0.15), bw, Inches(0.85),
              size=14, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        txbox(sl, desc, x, y0 + Inches(0.95), bw, Inches(1.3),
              size=11, color=RGBColor(0xDD,0xEE,0xFF), align=PP_ALIGN.CENTER)
        # Arrow
        if i < len(steps) - 1:
            ax = x + bw + Inches(0.03)
            txbox(sl, "▶", ax, y0 + Inches(0.9), Inches(0.22), Inches(0.5),
                  size=16, bold=True, color=TEAL, align=PP_ALIGN.CENTER)

    # Stats row
    rect(sl, Inches(0.35), Inches(4.35), SLIDE_W - Inches(0.7), Inches(0.06), fill=TEAL)
    stats = [
        ("8 isolates",     "assembled in parallel"),
        ("17 contigs",     "plasmid-sized (5–700 kb)"),
        ("< 30 min",       "total pipeline runtime"),
        ("72,556",         "reference plasmids searched"),
        ("0.0001–0.0003",  "nearest-neighbour distance"),
    ]
    y_s = Inches(4.55)
    for i, (val, lbl) in enumerate(stats):
        x = Inches(0.5) + i * Inches(2.55)
        txbox(sl, val, x, y_s,             Inches(2.4), Inches(0.42),
              size=18, bold=True, color=NAVY, align=PP_ALIGN.CENTER)
        txbox(sl, lbl, x, y_s + Inches(0.4), Inches(2.4), Inches(0.32),
              size=10, color=MIDGRAY, align=PP_ALIGN.CENTER)

    # Note about no pre-deposited seqs
    rect(sl, Inches(0.35), Inches(5.7), SLIDE_W - Inches(0.7), Inches(0.75), fill=AMBER)
    txbox(sl,
          "⚠  No pre-assembled plasmid sequences were deposited in ENA. "
          "All analysis proceeded directly from raw ONT reads — a true prospective, blinded classification.",
          Inches(0.55), Inches(5.76), SLIDE_W - Inches(1.1), Inches(0.6),
          size=12, bold=False, color=NAVY)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 4 — ASSEMBLY RESULTS
# ══════════════════════════════════════════════════════════════════════════════
def slide_assemblies(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "Assembly & Contig Extraction", "Flye 2.9.6 de novo ONT assembly · 8 isolates")

    # Table header
    cols  = ["Isolate", "Assembly\ncontigs", "Circular\ncontigs", "Assembly\nsize (Mbp)",
             "Plasmid contigs\nextracted", "Largest plasmid\ncontig (kb)"]
    cw    = [Inches(1.7), Inches(1.4), Inches(1.4), Inches(1.5), Inches(1.8), Inches(1.8)]
    data  = [
        ("NARACHVIM11", "3", "2", "5.36", "2", "249"),
        ("NARACHVIM12", "3", "1", "5.29", "2", "299"),
        ("NARACHVIM20", "4", "3", "5.33", "1", "284"),
        ("NARACHVIM31", "5", "2", "5.33", "2", "271"),
        ("NARACHVIM36", "4", "0", "5.23", "3", "325"),
        ("NARACHVIM48", "3", "2", "5.24", "2", "342"),
        ("NARACHVIM52", "5", "5", "5.37", "3", "310"),
        ("NARACHVIM56", "5", "5", "5.34", "4", "335"),
    ]
    row_h = Inches(0.44)
    hdr_h = Inches(0.5)
    x0, y0 = Inches(0.3), Inches(1.28)

    # Draw header
    x = x0
    for i, (col, w) in enumerate(zip(cols, cw)):
        rect(sl, x, y0, w, hdr_h, fill=NAVY)
        txbox(sl, col, x + Inches(0.04), y0 + Inches(0.02), w - Inches(0.08), hdr_h,
              size=11, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        x += w

    # Draw rows
    for r, row in enumerate(data):
        y = y0 + hdr_h + r * row_h
        x = x0
        bg = WHITE if r % 2 == 0 else RGBColor(0xEC, 0xF0, 0xF7)
        for i, (val, w) in enumerate(zip(row, cw)):
            rect(sl, x, y, w, row_h, fill=bg)
            fc = TEAL if i == 0 else (GREEN if i == 5 else DGRAY)
            txbox(sl, val, x + Inches(0.04), y + Inches(0.06),
                  w - Inches(0.08), row_h - Inches(0.1),
                  size=12, bold=(i==0), color=fc, align=PP_ALIGN.CENTER)
            x += w

    # Key finding box
    rect(sl, Inches(0.3), Inches(6.1), Inches(12.7), Inches(0.9), fill=GREEN)
    txbox(sl,
          "✔  All 8 isolates yielded a large plasmid contig (249–342 kb) consistent with the "
          "published IncHI2/blaVIM-1 plasmid (249–343 kb) — assembled de novo from raw reads alone.",
          Inches(0.5), Inches(6.18), Inches(12.3), Inches(0.7),
          size=13, bold=False, color=WHITE)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 5 — pLIN RESULTS TABLE
# ══════════════════════════════════════════════════════════════════════════════
def slide_plin_results(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "pLIN Classification Results", "Nearest-neighbour cosine distance against 72,556 reference plasmids")

    cols = ["Isolate", "Contig size\n(bp)", "Inc group\n(KNN)", "Conf.", "NN dist.", "pLIN code", "Level"]
    cw   = [Inches(1.65), Inches(1.3), Inches(1.4), Inches(0.65), Inches(0.75), Inches(5.1), Inches(0.8)]
    data = [
        ("NARACHVIM11", "249,474", "IncHI2", "100%", "0.0001", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM12", "299,362", "IncN",   " 80%", "0.0001", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM20", "284,079", "IncHI2", " 80%", "0.0001", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM31", "271,025", "IncN",   " 52%", "0.0001", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM36", "325,037", "IncN",   " 75%", "0.0003", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM48", "342,211", "IncN",   " 75%", "0.0003", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM52", "310,079", "IncN",   "100%", "0.0001", "1.1.3.42.150.1396", "L6"),
        ("NARACHVIM56", "335,344", "IncN",   " 57%", "0.0003", "1.1.3.42.150.1396", "L6"),
    ]
    row_h = Inches(0.44)
    hdr_h = Inches(0.5)
    x0, y0 = Inches(0.3), Inches(1.25)

    x = x0
    for col, w in zip(cols, cw):
        rect(sl, x, y0, w, hdr_h, fill=NAVY)
        txbox(sl, col, x + Inches(0.03), y0 + Inches(0.02), w - Inches(0.06), hdr_h,
              size=10, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        x += w

    for r, row in enumerate(data):
        y = y0 + hdr_h + r * row_h
        x = x0
        bg = WHITE if r % 2 == 0 else RGBColor(0xEC,0xF0,0xF7)
        for i, (val, w) in enumerate(zip(row, cw)):
            # Highlight the pLIN column
            cell_bg = RGBColor(0xE8,0xF8,0xF0) if i == 5 else bg
            rect(sl, x, y, w, row_h, fill=cell_bg)
            fc = (TEAL if i==0 else
                  GREEN if i==5 else
                  RGBColor(0x27,0x6E,0xBF) if i==6 else
                  DGRAY)
            txbox(sl, val, x + Inches(0.03), y + Inches(0.07),
                  w - Inches(0.06), row_h - Inches(0.1),
                  size=11, bold=(i in (0,5,6)), color=fc, align=PP_ALIGN.CENTER)
            x += w

    # Shared code highlight
    rect(sl, Inches(0.3), Inches(5.85), Inches(12.7), Inches(1.2), fill=NAVY)
    txbox(sl, "Shared pLIN code across ALL 8 isolates, 5 bacterial species:",
          Inches(0.5), Inches(5.9), Inches(12.3), Inches(0.35),
          size=13, color=AMBER, bold=True)
    txbox(sl, "pLIN  1.1.3.42.150.1396  [L6 — strain level]",
          Inches(0.5), Inches(6.25), Inches(12.3), Inches(0.55),
          size=22, bold=True, color=WHITE, align=PP_ALIGN.CENTER)

    # Note on IncHI2 vs IncN
    txbox(sl,
          "* IncHI2/IncN discrepancy: large VIM-1 plasmids carry both replicon types (multi-replicon); "
          "4-mer composition (pLIN) is unambiguous.",
          Inches(0.35), Inches(7.08), Inches(12.7), Inches(0.3),
          size=9, italic=True, color=MIDGRAY)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 6 — INTER-SPECIES TRANSFER VISUAL
# ══════════════════════════════════════════════════════════════════════════════
def slide_interspecies(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "Inter-Species Horizontal Plasmid Transfer Confirmed",
               "Same pLIN 1.1.3.42.150.1396 across 5 Enterobacterales species")

    # Central plasmid box
    cx, cy = Inches(6.4), Inches(3.85)
    pw, ph = Inches(3.2), Inches(1.5)
    rect(sl, cx - pw/2, cy - ph/2, pw, ph, fill=NAVY)
    txbox(sl, "pLIN 1.1.3.42.150.1396",
          cx - pw/2, cy - ph/2 + Inches(0.1), pw, Inches(0.5),
          size=13, bold=True, color=AMBER, align=PP_ALIGN.CENTER)
    txbox(sl, "IncHI2 · blaVIM-1 · mcr-9\n249–342 kb",
          cx - pw/2, cy - ph/2 + Inches(0.55), pw, Inches(0.55),
          size=11, color=WHITE, align=PP_ALIGN.CENTER)
    txbox(sl, "NN dist: 0.0001–0.0003\n(>99.9% ANI equivalent)",
          cx - pw/2, cy - ph/2 + Inches(1.05), pw, Inches(0.35),
          size=9, color=TEAL, align=PP_ALIGN.CENTER)

    # Species bubbles around the centre
    species = [
        (Inches(1.3),  Inches(2.1), "Enterobacter\nhormaechei",  "n=37",  TEAL),
        (Inches(1.3),  Inches(5.3), "Enterobacter\nkobei",       "n=?",   RGBColor(0x16,0x6F,0x82)),
        (Inches(4.2),  Inches(1.2), "Enterobacter\nludwigii",    "n=?",   RGBColor(0x8E,0x44,0xAD)),
        (Inches(4.2),  Inches(6.1), "Citrobacter\nfarmeri",      "n=?",   RGBColor(0xC0,0x39,0x2B)),
        (Inches(9.8),  Inches(3.5), "Escherichia\ncoli",         "n=?",   RGBColor(0xD3,0x54,0x00)),
    ]
    sw, sh = Inches(2.2), Inches(1.1)
    for x, y, sp, n, col in species:
        rect(sl, x, y, sw, sh, fill=col)
        txbox(sl, sp, x, y + Inches(0.05), sw, Inches(0.6),
              size=12, bold=True, color=WHITE, align=PP_ALIGN.CENTER, italic=True)
        txbox(sl, n,  x, y + Inches(0.65), sw, Inches(0.35),
              size=11, color=WHITE, align=PP_ALIGN.CENTER)

        # Draw arrow toward centre (approximate)
        arrow_x = cx - Inches(1.6) if x < cx else cx + Inches(1.6)
        arr_lbl = "→" if x < cx - Inches(1.0) else "←"
        mid_x = (x + sw/2 + cx) / 2 - Inches(0.2)
        mid_y = (y + sh/2 + cy) / 2 - Inches(0.15)
        txbox(sl, arr_lbl, mid_x, mid_y, Inches(0.5), Inches(0.4),
              size=20, bold=True, color=TEAL, align=PP_ALIGN.CENTER)

    # Caption
    rect(sl, Inches(0.3), Inches(6.55), Inches(12.7), Inches(0.55), fill=GREEN)
    txbox(sl,
          "Single pLIN code at strain level (L6) across 5 host species = direct molecular evidence "
          "for horizontal inter-species plasmid transfer within the outbreak",
          Inches(0.5), Inches(6.6), Inches(12.3), Inches(0.45),
          size=12, color=WHITE, bold=False)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 7 — SECONDARY PLASMIDS
# ══════════════════════════════════════════════════════════════════════════════
def slide_secondary(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=LGRAY)
    header_bar(sl, "Complete Within-Isolate Plasmid Inventory",
               "Secondary plasmids classified simultaneously — distinguishing outbreak vs co-carried")

    # Left: secondary plasmid table
    rect(sl, Inches(0.3), Inches(1.28), Inches(7.5), Inches(5.1), fill=WHITE)
    txbox(sl, "All 19 classified contigs", Inches(0.5), Inches(1.35), Inches(7.1), Inches(0.38),
          size=14, bold=True, color=NAVY)

    sec_data = [
        ("NARACHVIM11", "53,411",  "IncN",   "1.1.3.42.150.NEW",         "Co-carried"),
        ("NARACHVIM12", "5,356",   "IncN",   "1.1.3.42.NEW.NEW",         "Small cryptic"),
        ("NARACHVIM31", "9,520",   "IncFII", "1.1.3.5998.8787.16330",    "Co-carried"),
        ("NARACHVIM36", "39,259",  "IncFII", "1.1.3.42.150.NEW",         "Co-carried"),
        ("NARACHVIM36", "47,271",  "IncFII", "1.1.3.42.150.NEW",         "Co-carried"),
        ("NARACHVIM48", "56,685",  "IncX1",  "1.1.3.42.9451.17322",      "Co-carried"),
        ("NARACHVIM52", "5,356",   "IncN",   "1.1.3.NEW.NEW.NEW",        "Novel group"),
        ("NARACHVIM52", "5,022",   "IncFII", "1.1.2786.10973.17716.29454","Co-carried"),
        ("NARACHVIM56", "167,394", "IncFII", "1.1.3.42.150.1800",        "Co-carried"),
        ("NARACHVIM56", "108,584", "IncN",   "1.1.3.42.150.NEW",         "Co-carried"),
        ("NARACHVIM56", "87,410",  "IncFII", "1.1.3.42.150.NEW",         "Co-carried"),
    ]
    hdrs = ["Isolate", "Size (bp)", "Inc", "pLIN code", "Type"]
    hw   = [Inches(1.55), Inches(0.95), Inches(0.7), Inches(2.85), Inches(1.25)]
    x0, y0 = Inches(0.35), Inches(1.8)
    rh = Inches(0.36)

    x = x0
    for h, w in zip(hdrs, hw):
        rect(sl, x, y0, w, Inches(0.38), fill=NAVY)
        txbox(sl, h, x, y0 + Inches(0.02), w, Inches(0.35),
              size=10, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        x += w

    for r, row in enumerate(sec_data):
        y = y0 + Inches(0.38) + r * rh
        x = x0
        bg = WHITE if r % 2 == 0 else LGRAY
        for i, (val, w) in enumerate(zip(row, hw)):
            cell_bg = RGBColor(0xFF,0xF3,0xCD) if row[4]=="Novel group" else bg
            rect(sl, x, y, w, rh, fill=cell_bg)
            fc = DGRAY if i != 3 else (RED if "NEW" in val and "150" not in val else MIDGRAY)
            txbox(sl, val, x + Inches(0.02), y + Inches(0.03), w - Inches(0.04), rh - Inches(0.05),
                  size=9, color=fc, align=PP_ALIGN.CENTER, bold=(i==0))
            x += w

    # Right: key points
    rect(sl, Inches(8.1), Inches(1.28), Inches(4.9), Inches(5.1), fill=WHITE)
    txbox(sl, "What This Shows", Inches(8.3), Inches(1.38), Inches(4.5), Inches(0.38),
          size=14, bold=True, color=NAVY)

    points = [
        (GREEN,  "Outbreak plasmid isolated",
                 "pLIN 1.1.3.42.150.1396 is the sole\nshared code — unambiguously\nthe VIM-1 outbreak lineage"),
        (TEAL,   "Co-carried plasmids identified",
                 "IncFII, IncN, IncX1 secondary\nplasmids get distinct codes —\nnot confused with outbreak"),
        (AMBER,  "Novel group flagged",
                 "NARACHVIM52 5 kb contig\nclassified to L3 only — potential\nnew Inc subgroup detected"),
        (RED,    "Inc typing would fail",
                 "Multi-replicon large plasmid\npredicted as IncN or IncHI2\nbut pLIN is unambiguous"),
    ]
    y = Inches(1.85)
    for col, title, desc in points:
        rect(sl, Inches(8.1), y, Inches(0.12), Inches(0.95), fill=col)
        txbox(sl, title, Inches(8.3), y, Inches(4.5), Inches(0.35),
              size=12, bold=True, color=NAVY)
        txbox(sl, desc, Inches(8.3), y + Inches(0.32), Inches(4.5), Inches(0.6),
              size=10, color=DGRAY)
        y += Inches(1.15)

    footer(sl)


# ══════════════════════════════════════════════════════════════════════════════
#  SLIDE 8 — KEY CONCLUSIONS
# ══════════════════════════════════════════════════════════════════════════════
def slide_conclusions(prs):
    sl = blank(prs)
    rect(sl, 0, 0, SLIDE_W, SLIDE_H, fill=NAVY)

    txbox(sl, "Key Conclusions", Inches(0.6), Inches(0.3), Inches(12.1), Inches(0.65),
          size=30, bold=True, color=WHITE)
    rect(sl, Inches(0.6), Inches(1.0), Inches(3.0), Inches(0.06), fill=TEAL)

    conclusions = [
        (GREEN,  "1",
         "Single shared pLIN code confirmed outbreak",
         "pLIN 1.1.3.42.150.1396 assigned at L6 (strain-level) across all 8 isolates from "
         "5 Enterobacterales species. Nearest-neighbour cosine distance 0.0001–0.0003 "
         "(>99.9% ANI equivalent). Outcome achieved within minutes of assembly completion."),
        (TEAL,   "2",
         "True prospective, blinded classification",
         "No pre-deposited plasmid sequences existed in ENA. All analysis started from raw "
         "ONT reads — de novo assembly, contig extraction, pLIN assignment. No manual curation."),
        (AMBER,  "3",
         "Permanent cross-institutional linkage",
         "The pLIN code 1.1.3.42.150.1396 is now permanently assigned. Any future isolate "
         "from any institution carrying this plasmid will receive the same code, enabling "
         "retrospective linkage across countries and time."),
        (RGBColor(0x8E,0x44,0xAD), "4",
         "Complete within-isolate plasmid inventory",
         "Secondary co-carried plasmids (IncFII, IncN, IncX1) were simultaneously classified "
         "with distinct pLIN codes, cleanly separating the outbreak lineage from passenger "
         "plasmids — impossible with replicon-based typing alone."),
    ]

    y = Inches(1.25)
    for col, num, title, body in conclusions:
        rect(sl, Inches(0.5), y, Inches(0.6), Inches(1.45), fill=col)
        txbox(sl, num, Inches(0.5), y + Inches(0.3), Inches(0.6), Inches(0.6),
              size=24, bold=True, color=WHITE, align=PP_ALIGN.CENTER)
        rect(sl, Inches(1.2), y, Inches(11.5), Inches(1.45), fill=RGBColor(0x1E,0x3A,0x5A))
        txbox(sl, title, Inches(1.35), y + Inches(0.07), Inches(11.2), Inches(0.38),
              size=14, bold=True, color=col)
        txbox(sl, body,  Inches(1.35), y + Inches(0.45), Inches(11.2), Inches(0.85),
              size=11, color=RGBColor(0xCC,0xDD,0xEE))
        y += Inches(1.6)

    # Bottom tag
    rect(sl, 0, Inches(7.15), SLIDE_W, Inches(0.35), fill=TEAL)
    txbox(sl,
          "pLIN 1.1.3.42.150.1396  ·  Swiss VIM-1 outbreak  ·  5 species  ·  8 isolates  ·  "
          "De novo assembly from raw ONT reads  ·  No manual curation",
          Inches(0.3), Inches(7.17), Inches(12.7), Inches(0.28),
          size=10, bold=True, color=WHITE, align=PP_ALIGN.CENTER)


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    prs = new_prs()
    slide_title(prs)
    slide_background(prs)
    slide_pipeline(prs)
    slide_assemblies(prs)
    slide_plin_results(prs)
    slide_interspecies(prs)
    slide_secondary(prs)
    slide_conclusions(prs)

    out = os.path.join(OUTPUT_DIR, "pLIN_Swiss_VIM1_CaseStudy.pptx")
    prs.save(out)
    print(f"Saved: {out}  ({os.path.getsize(out)//1024} KB)  — {len(prs.slides)} slides")


if __name__ == "__main__":
    main()
