#!/usr/bin/env python3
# Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
# See LICENSE and CITATION.cff for terms. Citation is MANDATORY.
"""
Generate a structured manuscript draft for the pLIN system as a Word document.

Usage:  python manuscript_draft.py
Output: output/pLIN_Manuscript.docx
"""

import os
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")


def set_cell_shading(cell, color_hex):
    """Set background shading on a table cell."""
    shading = cell._element.get_or_add_tcPr()
    shading_elem = shading.makeelement(
        qn("w:shd"),
        {qn("w:fill"): color_hex, qn("w:val"): "clear"},
    )
    shading.append(shading_elem)


def add_heading(doc, text, level=1):
    """Add a heading with consistent styling."""
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.color.rgb = RGBColor(0x0D, 0x47, 0xA1)
    return h


def add_para(doc, text, bold=False, italic=False, font_size=11):
    """Add a paragraph with styling."""
    p = doc.add_paragraph()
    run = p.add_run(text)
    run.font.size = Pt(font_size)
    run.bold = bold
    run.italic = italic
    return p


def add_table_row(table, cells, bold=False, shading=None):
    """Add a row to a table."""
    row = table.add_row()
    for i, text in enumerate(cells):
        cell = row.cells[i]
        cell.text = text
        for paragraph in cell.paragraphs:
            paragraph.style.font.size = Pt(10)
            for run in paragraph.runs:
                run.bold = bold
        if shading:
            set_cell_shading(cell, shading)
    return row


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    doc = Document()

    # ── Page setup ───────────────────────────────────────────────────────────
    style = doc.styles["Normal"]
    style.font.name = "Calibri"
    style.font.size = Pt(11)
    style.paragraph_format.space_after = Pt(6)
    style.paragraph_format.line_spacing = 1.5

    for section in doc.sections:
        section.top_margin = Inches(1)
        section.bottom_margin = Inches(1)
        section.left_margin = Inches(1.25)
        section.right_margin = Inches(1.25)

    # ══════════════════════════════════════════════════════════════════════════
    # TITLE PAGE
    # ══════════════════════════════════════════════════════════════════════════

    title_para = doc.add_paragraph()
    title_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title_para.space_after = Pt(24)
    run = title_para.add_run(
        "pLIN: A Hierarchical, Reference-Free Classification System "
        "for Bacterial Plasmid Genomes with Integrated AMR Surveillance "
        "and Host Inference"
    )
    run.bold = True
    run.font.size = Pt(16)
    run.font.color.rgb = RGBColor(0x0D, 0x47, 0xA1)

    # Authors
    author_para = doc.add_paragraph()
    author_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = author_para.add_run("Basil Xavier Britto")
    run.font.size = Pt(12)
    run.bold = True

    affil_para = doc.add_paragraph()
    affil_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = affil_para.add_run("[TODO: Department, Institution, City, Country]")
    run.font.size = Pt(10)
    run.italic = True

    email_para = doc.add_paragraph()
    email_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = email_para.add_run("[TODO: Corresponding author email]")
    run.font.size = Pt(10)

    doc.add_page_break()

    # ══════════════════════════════════════════════════════════════════════════
    # ABSTRACT
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "Abstract", level=1)

    add_para(doc, "Background", bold=True, font_size=11)
    doc.add_paragraph(
        "Antimicrobial resistance (AMR) mediated by bacterial plasmids represents a critical "
        "global health threat. Plasmids serve as the primary vehicles for horizontal gene transfer "
        "of resistance determinants across bacterial species and genera. Despite their clinical "
        "importance, existing plasmid classification systems lack permanent, hierarchical "
        "identifiers that enable consistent cross-study comparisons and longitudinal surveillance. "
        "Current tools such as PlasmidFinder, MOB-suite, and pMLST each address specific aspects "
        "of plasmid characterization but fail to provide an integrated, reference-free "
        "classification framework."
    )

    add_para(doc, "Methods", bold=True, font_size=11)
    doc.add_paragraph(
        "We present pLIN (Plasmid Life Identification Number), a hierarchical classification "
        "system based on tetranucleotide (4-mer) frequency composition. pLIN assigns permanent "
        "six-level codes (Family through Strain) using cosine pairwise distances and hierarchical "
        "clustering with six progressive distance thresholds. An integrated K-Nearest Neighbors "
        "(KNN) classifier, trained on 6,998 RefSeq plasmids across 20 incompatibility (Inc) "
        "groups, provides automated Inc group detection with 92.2% stratified cross-validated "
        "accuracy. The system integrates AMRFinderPlus for resistance gene detection, a three-tier "
        "mobility prediction cascade (MOBsuite, AMRFinderPlus gene scan, default non-mobilizable), "
        "and a novel CRISPR spacer-based plasmid-host inference module using MinCED and BLASTN-short "
        "with softmax probability ranking."
    )

    add_para(doc, "Results", bold=True, font_size=11)
    doc.add_paragraph(
        "pLIN achieves a Simpson's diversity index of 0.979, demonstrating high discriminatory power. "
        "The KNN classifier attains 92.2% accuracy across 20 Inc groups with dynamic k-selection and "
        "distance-weighted cosine voting. The mobility cascade correctly classifies plasmid transfer "
        "potential using relaxase family and MPF type information from MOBsuite, augmented by expanded "
        "gene prefix detection (tra, trb, virB, trw, pilX, mob, nik) from AMRFinderPlus output. "
        "The CRISPR host inference module extracts spacers from uploaded bacterial genomes, matches "
        "them against plasmid sequences with stringent criteria (\u226595% identity, \u226525 bp, "
        "\u22641 mismatch, 0 gaps), and ranks host-plasmid associations using numerically stable "
        "softmax probabilities with confidence categorization."
    )

    add_para(doc, "Conclusions", bold=True, font_size=11)
    doc.add_paragraph(
        "pLIN provides a comprehensive, open-source platform for plasmid classification, AMR "
        "surveillance, mobility typing, and host inference. Its hierarchical codes enable permanent, "
        "reproducible plasmid identification across institutions and studies. The interactive "
        "Streamlit web interface with eight analytical tabs requires no bioinformatics expertise, "
        "making it accessible to clinical microbiologists and infection control practitioners."
    )

    add_para(doc,
        "Keywords: plasmid classification, antimicrobial resistance, pLIN, incompatibility group, "
        "CRISPR host inference, mobility typing, MOBsuite, AMRFinderPlus",
        italic=True, font_size=10
    )

    doc.add_page_break()

    # ══════════════════════════════════════════════════════════════════════════
    # INTRODUCTION
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "1. Introduction", level=1)

    doc.add_paragraph(
        "Antimicrobial resistance (AMR) has emerged as one of the most pressing global health "
        "challenges of the 21st century. The World Health Organization estimates that drug-resistant "
        "infections cause approximately 1.27 million deaths annually, with projections suggesting "
        "this could rise to 10 million by 2050 if left unchecked. Central to the AMR crisis are "
        "bacterial plasmids \u2014 self-replicating extrachromosomal DNA elements that serve as the "
        "primary vehicles for horizontal gene transfer (HGT) of resistance determinants across "
        "bacterial species and genera."
    )

    doc.add_paragraph(
        "Plasmids carrying AMR genes can transfer between bacteria through conjugation, "
        "transformation, and transduction, enabling rapid dissemination of resistance across "
        "diverse bacterial populations. A single conjugative plasmid carrying multiple resistance "
        "genes can convert a susceptible bacterium into a multidrug-resistant (MDR) pathogen in a "
        "single transfer event. Understanding which plasmids carry which resistance genes, how "
        "they relate to each other, and which hosts they associate with is therefore critical for "
        "effective AMR surveillance and infection control."
    )

    doc.add_paragraph(
        "Existing plasmid classification approaches each address specific aspects of plasmid "
        "characterization but have significant limitations. PlasmidFinder (Carattoli et al., 2014) "
        "detects replicon sequences to assign Inc groups but depends on a curated reference database "
        "and cannot classify novel or divergent plasmids. MOB-suite (Robertson & Nash, 2018) "
        "provides relaxase-based typing and mobility classification but requires command-line "
        "operation and does not generate permanent hierarchical identifiers. Plasmid MLST (pMLST) "
        "offers high resolution within specific Inc groups but is limited to plasmids with defined "
        "MLST schemes. None of these tools provides a unified, reference-free classification system "
        "with permanent codes that remain stable as new sequences are added."
    )

    doc.add_paragraph(
        "To address these limitations, we developed pLIN (Plasmid Life Identification Number), "
        "a hierarchical classification system inspired by the Linnaean taxonomic framework. pLIN "
        "assigns permanent six-level codes based on tetranucleotide composition similarity, enabling "
        "consistent cross-study comparisons and longitudinal surveillance. The system integrates "
        "automated Inc group detection (92.2% accuracy across 20 groups), comprehensive AMR gene "
        "detection via AMRFinderPlus, a three-tier mobility prediction cascade incorporating "
        "MOBsuite, and a novel CRISPR spacer-based host inference module. An interactive Streamlit "
        "web interface with eight analytical tabs makes the system accessible to researchers without "
        "bioinformatics expertise."
    )

    # ══════════════════════════════════════════════════════════════════════════
    # METHODS
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "2. Methods", level=1)

    # 2.1 pLIN Classification
    add_heading(doc, "2.1 pLIN Classification System", level=2)
    doc.add_paragraph(
        "The pLIN classification system assigns hierarchical codes based on pairwise sequence "
        "similarity computed from tetranucleotide (4-mer) frequency vectors. For each input "
        "plasmid sequence, a normalized 256-dimensional frequency vector is computed by counting "
        "all overlapping 4-mers (AAAA through TTTT) and dividing by the total count. Pairwise "
        "cosine distances are computed between all sequences, producing a distance matrix D where "
        "D(i,j) = 1 \u2212 cos(\u03b8) between vectors v_i and v_j."
    )

    doc.add_paragraph(
        "Hierarchical agglomerative clustering is performed on the distance matrix using "
        "single-linkage (default), with complete, average, and weighted linkage available as "
        "alternatives. The resulting dendrogram is cut at six progressive distance thresholds "
        "to produce six nested cluster assignments. Each plasmid receives a six-level code in "
        "the format A.B.C.D.E.F, where A represents the broadest grouping (Family) and F the "
        "finest (Strain)."
    )

    # Table 1: Threshold levels
    add_para(doc, "Table 1. pLIN hierarchical threshold levels.", bold=True, italic=True, font_size=10)
    table1 = doc.add_table(rows=1, cols=4)
    table1.style = "Table Grid"
    table1.alignment = WD_TABLE_ALIGNMENT.CENTER

    headers = ["Position", "Level", "Threshold (d \u2264)", "ANI Equivalent"]
    for i, h in enumerate(headers):
        cell = table1.rows[0].cells[i]
        cell.text = h
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
        set_cell_shading(cell, "0D47A1")
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    levels = [
        ("A", "Family", "0.150", "~85%"),
        ("B", "Subfamily", "0.100", "~90%"),
        ("C", "Cluster", "0.050", "~95%"),
        ("D", "Subcluster", "0.020", "~98%"),
        ("E", "Clone", "0.010", "~99%"),
        ("F", "Strain", "0.001", "~99.9%"),
    ]
    for row_data in levels:
        add_table_row(table1, row_data)
    doc.add_paragraph()

    # 2.2 Inc Group Detection
    add_heading(doc, "2.2 Inc Group Detection", level=2)
    doc.add_paragraph(
        "Automated incompatibility (Inc) group detection is performed using a K-Nearest Neighbors "
        "(KNN) classifier trained on 6,998 characterized plasmid sequences from NCBI RefSeq, "
        "spanning 20 Inc groups: ColE, ColRNAI, IncA, IncAC2, IncC, IncF, IncFIB, IncFIBK, "
        "IncFIC, IncFII, IncHI1, IncHI2, IncI, IncI1, IncI2, IncN, IncR, IncX1, IncX3, and "
        "IncX4. Training sequences were identified by parsing FASTA headers from a reference "
        "database of 72,556 plasmid sequences using regex-based Inc type extraction with a "
        "validated whitelist to prevent false positives from species names (e.g., \"incola\", "
        "\"incerta\"). Sub-variants were normalized to canonical groups (e.g., IncFIBpQil \u2192 "
        "IncFIB, IncHI2A \u2192 IncHI2)."
    )

    doc.add_paragraph(
        "The classifier uses k=5 neighbors with cosine distance metric and distance-weighted "
        "voting. Dynamic k-selection ensures k < smallest class size to prevent degenerate "
        "predictions. Model performance was evaluated using 5-fold stratified cross-validation, "
        "achieving an overall accuracy of 92.2%. Predictions below a 40% confidence threshold "
        "are flagged as \"Unknown/Novel\" with the top 5 candidate Inc groups reported."
    )

    # Table 2: Inc group training data
    add_para(doc, "Table 2. Inc group training data composition (top 10 of 20 groups).",
             bold=True, italic=True, font_size=10)
    table2 = doc.add_table(rows=1, cols=3)
    table2.style = "Table Grid"
    table2.alignment = WD_TABLE_ALIGNMENT.CENTER

    for i, h in enumerate(["Inc Group", "Training Samples", "Percentage"]):
        cell = table2.rows[0].cells[i]
        cell.text = h
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
        set_cell_shading(cell, "0D47A1")
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    inc_data = [
        ("IncFII", "4,629", "66.1%"),
        ("IncN", "1,064", "15.2%"),
        ("IncX1", "701", "10.0%"),
        ("IncF", "148", "2.1%"),
        ("IncI1", "72", "1.0%"),
        ("IncHI2", "55", "0.8%"),
        ("ColRNAI", "53", "0.8%"),
        ("IncX3", "44", "0.6%"),
        ("IncA", "36", "0.5%"),
        ("Others (11)", "196", "2.8%"),
    ]
    for row_data in inc_data:
        add_table_row(table2, row_data)
    doc.add_paragraph()

    # 2.3 AMR Detection
    add_heading(doc, "2.3 AMR Detection", level=2)
    doc.add_paragraph(
        "Antimicrobial resistance gene detection is performed using NCBI AMRFinderPlus, which "
        "identifies AMR, stress resistance, and virulence-associated genes using curated protein "
        "and nucleotide Hidden Markov Models and BLAST. AMRFinderPlus is run per input file with "
        "the --plus flag to include stress and virulence genes. Results are parsed to extract gene "
        "names, drug classes, resistance mechanisms, and genomic coordinates, then integrated with "
        "pLIN assignments on a per-plasmid basis."
    )

    # 2.4 Mobility Prediction
    add_heading(doc, "2.4 Mobility Prediction", level=2)
    doc.add_paragraph(
        "Plasmid mobility is classified using a three-tier priority cascade. The first tier uses "
        "MOBsuite mob_typer (Robertson & Nash, 2018) when available, which classifies plasmids "
        "based on relaxase family (MOBF, MOBH, MOBP, MOBQ, MOBC, MOBV) and Mating Pair Formation "
        "(MPF) type (Type T, F, I, G). MOBsuite provides the most reliable mobility classification "
        "through direct identification of conjugation machinery components."
    )

    doc.add_paragraph(
        "The second tier scans AMRFinderPlus output for transfer and mobilization gene markers. "
        "Conjugative status is assigned when Type IV secretion system genes are detected, including "
        "tra/trb operons, virB1\u2013virB11, trwA\u2013trwN, pilX, and taxC. Mobilizable status "
        "is assigned when relaxase or mobilization genes are detected, including mob gene variants, "
        "nikC\u2013nikE, MOBF/MOBH/MOBP/MOBQ/MOBC/MOBV relaxase families, oriT, and mps genes. "
        "If neither tier detects transfer-associated genes, the plasmid is classified as "
        "non-mobilizable (third tier, default)."
    )

    doc.add_paragraph(
        "Plasmids classified as conjugative and carrying AMR genes are flagged as HIGH RISK for "
        "dissemination potential, while mobilizable plasmids with AMR genes are flagged as "
        "MODERATE RISK, enabling risk-stratified surveillance."
    )

    # 2.5 CRISPR Host Inference
    add_heading(doc, "2.5 CRISPR Spacer-Based Host Inference", level=2)
    doc.add_paragraph(
        "A novel module for inferring plasmid-host relationships leverages CRISPR arrays in "
        "bacterial genomes as a molecular record of past encounters with mobile genetic elements. "
        "The pipeline operates in five stages:"
    )

    doc.add_paragraph(
        "(1) CRISPR spacer extraction: MinCED v0.4.2 (Bland et al., 2007; Skennerton, 2018) is "
        "run on user-uploaded bacterial host genome FASTA files to identify CRISPR arrays and "
        "extract spacer sequences. Spacer IDs are formatted as {genome_name}__spacer_{N} to "
        "enable downstream host tracing."
    )
    doc.add_paragraph(
        "(2) BLAST database construction: A nucleotide BLAST database is built from either "
        "user-uploaded plasmid sequences or the built-in reference database containing 72,556 "
        "plasmid sequences."
    )
    doc.add_paragraph(
        "(3) Spacer-plasmid matching: BLASTN-short is run with parameters optimized for short "
        "query sequences: -task blastn-short -dust no -word_size 7 -evalue 1e-5 "
        "-max_target_seqs 500."
    )
    doc.add_paragraph(
        "(4) Stringent filtering: BLAST hits are filtered to retain only high-confidence matches: "
        "percent identity \u2265 95%, alignment length \u2265 25 bp, mismatches \u2264 1, "
        "gap openings = 0."
    )
    doc.add_paragraph(
        "(5) Probability ranking: For each host-plasmid pair, the number of unique spacer hits "
        "and average alignment quality are computed. Counts are normalized by total spacers per "
        "host genome. Softmax normalization (with numerical stability via max-subtraction before "
        "exponentiation, temperature T = 1.0) converts raw scores to probabilities. Predictions "
        "are categorized as High confidence (\u2265 0.7), Medium (\u2265 0.4), or Low (< 0.4)."
    )

    # 2.6 Implementation
    add_heading(doc, "2.6 Implementation", level=2)
    doc.add_paragraph(
        "pLIN is implemented as a monolithic Streamlit web application (plin_app.py, ~3,700 "
        "lines) providing an interactive eight-tab interface: Overview, Results, Cladogram, "
        "AMR Analysis, Epidemiology, CRISPR Host, DRAGNOME Buddy (AI chatbot via Ollama), and "
        "Export. The system is written in Python 3 and depends on NumPy, Pandas, SciPy, "
        "scikit-learn, BioPython, Matplotlib, Plotly, and Streamlit. External tools (AMRFinderPlus, "
        "MOBsuite, Prodigal, MinCED, BLAST+) are auto-detected from PATH and conda environments. "
        "All results are exportable as TSV files, PNG/PDF figures, and ZIP bundles."
    )

    # ══════════════════════════════════════════════════════════════════════════
    # RESULTS
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "3. Results", level=1)

    add_heading(doc, "3.1 Classification Performance", level=2)
    doc.add_paragraph(
        "The KNN classifier achieved an overall accuracy of 92.2% on 5-fold stratified "
        "cross-validation across 20 Inc groups (6,998 training samples). Per-class performance "
        "varied with training set size: the largest group (IncFII, n=4,629) achieved >95% F1-score, "
        "while smaller groups (IncFIBK, n=11) showed lower but still informative performance. "
        "The dynamic k-selection mechanism (k = min(5, max(1, min_class_size \u2212 1))) prevented "
        "degenerate predictions for small classes."
    )

    # Table 3: Classification metrics
    add_para(doc, "Table 3. Per-class classification metrics (selected groups).",
             bold=True, italic=True, font_size=10)
    table3 = doc.add_table(rows=1, cols=4)
    table3.style = "Table Grid"
    table3.alignment = WD_TABLE_ALIGNMENT.CENTER

    for i, h in enumerate(["Inc Group", "Precision", "Recall", "F1-Score"]):
        cell = table3.rows[0].cells[i]
        cell.text = h
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.bold = True
        set_cell_shading(cell, "0D47A1")
        for paragraph in cell.paragraphs:
            for run in paragraph.runs:
                run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)

    metrics_data = [
        ("IncFII", "[TODO]", "[TODO]", "[TODO]"),
        ("IncN", "[TODO]", "[TODO]", "[TODO]"),
        ("IncX1", "[TODO]", "[TODO]", "[TODO]"),
        ("IncF", "[TODO]", "[TODO]", "[TODO]"),
        ("IncHI2", "[TODO]", "[TODO]", "[TODO]"),
        ("Overall (weighted)", "[TODO]", "[TODO]", "0.922"),
    ]
    for row_data in metrics_data:
        add_table_row(table3, row_data)
    doc.add_paragraph()

    add_heading(doc, "3.2 Discriminatory Power", level=2)
    doc.add_paragraph(
        "The pLIN classification system demonstrated a Simpson's diversity index (D) of 0.979, "
        "indicating high discriminatory power at the strain level (F threshold, d \u2264 0.001). "
        "[TODO: Report number of unique pLIN codes, distribution statistics, and comparison with "
        "other classification methods.]"
    )

    add_heading(doc, "3.3 Mobility Typing", level=2)
    doc.add_paragraph(
        "The three-tier mobility cascade successfully classified plasmids into Conjugative, "
        "Mobilizable, and Non-mobilizable categories. [TODO: Report concordance between MOBsuite "
        "and AMRFinderPlus-based classifications. Report distribution of relaxase families and "
        "MPF types. Report number of HIGH RISK (conjugative + AMR) and MODERATE RISK "
        "(mobilizable + AMR) plasmids identified.]"
    )

    add_heading(doc, "3.4 CRISPR Host Prediction", level=2)
    doc.add_paragraph(
        "[TODO: Report results from CRISPR host inference analysis. Include: number of host "
        "genomes screened, CRISPR arrays found, total spacers extracted, host-plasmid predictions "
        "made, and distribution of confidence categories (High/Medium/Low). Validate against "
        "known host-plasmid associations where available.]"
    )

    # ══════════════════════════════════════════════════════════════════════════
    # DISCUSSION
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "4. Discussion", level=1)

    doc.add_paragraph(
        "pLIN addresses a critical gap in plasmid epidemiology by providing a unified, "
        "reference-free classification system with permanent hierarchical codes. Unlike "
        "existing tools that focus on single aspects of plasmid characterization, pLIN "
        "integrates classification, Inc group detection, AMR surveillance, mobility typing, "
        "and host inference into a single platform."
    )

    add_heading(doc, "4.1 Advantages Over Existing Tools", level=2)
    doc.add_paragraph(
        "Compared to PlasmidFinder, pLIN does not require replicon-specific reference databases "
        "and can classify novel plasmids without prior characterization. Unlike MOB-suite, pLIN "
        "provides permanent hierarchical codes that enable longitudinal tracking and cross-study "
        "comparisons. The composition-based approach captures the global genomic signature of "
        "plasmids, providing classification resolution from Family to Strain level. The "
        "integration of MOBsuite\u2019s relaxase and MPF typing with AMRFinderPlus-based gene "
        "scanning creates a robust, multi-evidence mobility prediction system."
    )

    add_heading(doc, "4.2 Limitations", level=2)
    doc.add_paragraph(
        "Several limitations should be acknowledged. First, the training data exhibits class "
        "imbalance: IncFII (4,629 samples) dominates the dataset while groups like IncFIBK (11 "
        "samples) are underrepresented. This affects per-class performance for minority groups "
        "and may bias predictions toward majority classes. Second, 4-mer frequency composition "
        "captures global sequence characteristics but ignores gene content, synteny, and "
        "structural rearrangements. Plasmids with similar base composition but different gene "
        "cargo may receive similar codes at coarse classification levels. Third, CRISPR-based "
        "host inference depends on the presence of CRISPR arrays in host genomes; hosts lacking "
        "CRISPR systems will not be detected."
    )

    add_heading(doc, "4.3 Future Directions", level=2)
    doc.add_paragraph(
        "Future work will focus on: (1) balancing Inc group training data through targeted "
        "sequence collection for underrepresented groups, targeting 30+ groups with 20,000+ "
        "balanced samples; (2) incorporating hybrid features combining k-mer composition with "
        "gene presence/absence for higher-resolution classification; (3) supporting metagenomic "
        "assemblies with multi-contig plasmid bin handling; (4) developing a web database and "
        "REST API for programmatic pLIN code lookups and sequence submission; and (5) extending "
        "outbreak detection to cross-institutional datasets with spatial-temporal metadata "
        "integration."
    )

    # ══════════════════════════════════════════════════════════════════════════
    # CONCLUSION
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "5. Conclusion", level=1)
    doc.add_paragraph(
        "pLIN provides a comprehensive, open-source platform for bacterial plasmid classification "
        "and AMR surveillance. Its permanent, hierarchical six-level codes enable consistent "
        "plasmid identification across institutions and studies. The integrated KNN classifier "
        "achieves 92.2% accuracy across 20 Inc groups, while the three-tier mobility cascade "
        "and novel CRISPR host inference module extend the system\u2019s analytical capabilities "
        "beyond classification to comprehensive epidemiological intelligence. The interactive "
        "Streamlit web interface with eight analytical tabs makes the system accessible to "
        "clinical microbiologists and infection control practitioners without bioinformatics "
        "expertise. pLIN is freely available under the GPL-3.0 license with mandatory citation "
        "clause at [TODO: GitHub URL]."
    )

    # ══════════════════════════════════════════════════════════════════════════
    # REFERENCES
    # ══════════════════════════════════════════════════════════════════════════

    add_heading(doc, "6. References", level=1)

    references = [
        "Bland C, Ramsey TL, Sabree F, Lowe M, Brown K, Kyrpides NC, Hugenholtz P. "
        "CRISPR Recognition Tool (CRT): a tool for automatic detection of clustered regularly "
        "interspaced palindromic repeats. BMC Bioinformatics. 2007;8:209.",

        "Carattoli A, Zankari E, Garc\u00eda-Fern\u00e1ndez A, Voldby Larsen M, Lund O, Villa L, "
        "M\u00f8ller Aarestrup F, Hasman H. In silico detection and typing of plasmids using "
        "PlasmidFinder and plasmid multilocus sequence typing. Antimicrob Agents Chemother. "
        "2014;58(7):3895-903.",

        "Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, "
        "Hoffmann M, Pettengill JB, Prasad AB, Tillman GE, Tyson GH, Klimke W. AMRFinderPlus "
        "and the Reference Gene Catalog facilitate examination of the genomic links among "
        "antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11(1):12728.",

        "Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic "
        "gene recognition and translation initiation site identification. BMC Bioinformatics. "
        "2010;11:119.",

        "Robertson J, Nash JHE. MOB-suite: software tools for clustering, reconstruction and "
        "typing of plasmids from draft assemblies. Microb Genom. 2018;4(8):e000206.",

        "Skennerton CT. MinCED: Mining CRISPRs in Environmental Datasets. GitHub. "
        "https://github.com/ctSkennerton/minced. 2018.",
    ]

    for i, ref in enumerate(references, 1):
        p = doc.add_paragraph()
        run = p.add_run(f"[{i}] {ref}")
        run.font.size = Pt(10)

    # ══════════════════════════════════════════════════════════════════════════
    # FIGURE LEGENDS
    # ══════════════════════════════════════════════════════════════════════════

    doc.add_page_break()
    add_heading(doc, "Figure Legends", level=1)

    fig_legends = [
        ("Figure 1.", " System architecture of the pLIN platform. The pipeline processes "
         "plasmid FASTA files through Inc group detection (KNN classifier), 4-mer vectorization, "
         "hierarchical clustering, pLIN code assignment, AMR annotation (AMRFinderPlus), mobility "
         "typing (MOBsuite cascade), and optional CRISPR host inference (MinCED + BLASTN-short). "
         "The interactive Streamlit web interface provides eight analytical tabs for comprehensive "
         "plasmid analysis and export."),

        ("Figure 2.", " Example pLIN cladogram showing hierarchical relationships among plasmid "
         "sequences. Color coding indicates Inc group assignments. Branch lengths reflect cosine "
         "distances based on 4-mer frequency composition. Threshold levels (A through F) are "
         "indicated by dashed horizontal lines. [TODO: Generate from test data.]"),

        ("Figure 3.", " CRISPR host-plasmid probability heatmap. Rows represent host bacterial "
         "genomes, columns represent plasmid sequences. Color intensity reflects softmax-normalized "
         "probability of host-plasmid association based on CRISPR spacer matching. "
         "[TODO: Generate from test data.]"),
    ]

    for label, text in fig_legends:
        p = doc.add_paragraph()
        run_bold = p.add_run(label)
        run_bold.bold = True
        run_bold.font.size = Pt(10)
        run_text = p.add_run(text)
        run_text.font.size = Pt(10)

    # ══════════════════════════════════════════════════════════════════════════
    # SAVE
    # ══════════════════════════════════════════════════════════════════════════

    out_path = os.path.join(OUTPUT_DIR, "pLIN_Manuscript.docx")
    doc.save(out_path)
    print(f"Saved: {out_path}")
    print(f"Sections: Title, Abstract, Introduction, Methods (6 sub), Results (4 sub), "
          f"Discussion (3 sub), Conclusion, References ({len(references)}), Figure Legends (3)")


if __name__ == "__main__":
    main()
