# Main Novel Findings of the pLIN Classification System

## 1. First Integrated Plasmid Classification + AMR Surveillance + Outbreak Detection System

pLIN is the first tool that simultaneously provides:
- **Hierarchical multi-resolution classification** (6 levels, L1-L6, from family to strain)
- **Integrated AMR gene profiling** (AMRFinderPlus v4.2.5)
- **Automated outbreak detection** with clinical risk stratification (CRITICAL/HIGH/MODERATE)
- **Permanent nomenclature** that never changes with database updates

No existing tool (PlasmidFinder, MOB-suite, COPLA, mge-cluster, pMLST) combines all four capabilities. PlasmidFinder provides replicon typing only. MOB-suite clusters by relaxase but reassigns codes with every database update. Neither integrates AMR data or provides outbreak detection.

---

## 2. Discovery of High-Risk Plasmid Lineages Invisible to Conventional Typing

### The "KPC-IncN" Lineage (pLIN 1.1.2.15.48.671)
- **90 plasmids** carrying blaKPC-2 on **100%** of members
- Tightly conserved MDR cassette: blaTEM-1 (100%), aph(3'')-Ib (98%), aac(3)-IId (97%)
- **13.2 mean AMR genes** per plasmid
- Confers simultaneous resistance to carbapenems, penicillins, and aminoglycosides
- **Clinical impact:** Conventional typing reports this as "IncN" -- a label shared by >1,000 plasmids with vastly different resistance profiles. pLIN uniquely identifies this 100% carbapenemase-positive subpopulation as a concrete surveillance target

### Cross-Inc Resistance Hub (pLIN 1.1.2.15.48.860)
- **142 plasmids** spanning **5 incompatibility groups** (IncFII, IncHI1, IncHI2, IncN, IncX1)
- Top AMR genes: sul1 (61.3%), floR (56.3%), mph(A)/mrx(A) (54.2%), blaTEM-1 (38.0%), mcr-1.1 (36.6%)
- **14.4 mean AMR genes** per plasmid, **44.4% mcr colistin resistance** carriage (63/142 any mcr variant)
- **Novel concept:** Represents a convergence zone where resistance cassettes are exchanged between structurally different plasmid families -- entirely invisible to replicon-based or relaxase-based typing

---

## 3. Comprehensive AMR Landscape Across 6,998 Clinical Plasmids

- **83.1% of plasmids** carry resistance determinants (5,816 of 6,998)
- **64,891 total gene detections:** 29,583 AMR + 6,286 virulence + 29,022 stress response
- **1,635 carbapenemase** detections (blaKPC-2 n=824, blaNDM-1 n=228)
- **1,804 ESBL** detections (blaCTX-M-15 n=505, blaCTX-M-65 n=319, blaSHV-12 n=277)
- **204 colistin resistance (mcr)** detections -- threatening the last-resort polymyxin class
- **2,315 PMQR** detections (qnrS1 n=732, aac(6')-Ib-cr5 n=737)
- **412 plasmids** co-carry carbapenemase + ESBL genes -- pan-beta-lactam resistance
- **64.7% of AMR+ plasmids** carry resistance to 3+ drug classes; **39.7%** to 5+ drug classes

---

## 4. 1.54-Fold Improvement in Discriminatory Power Over Inc/Rep Typing

- pLIN resolves **3,073 unique strain-level codes** from 8,077 training samples (28 groups)
- Simpson's Diversity Index: **D = 0.985** (pLIN) vs **D = 0.641** (Inc/rep typing alone)
- **97.8%** of L6 codes contain plasmids from a single group, confirming pLIN captures known taxonomy as an emergent property
- The hierarchical structure enables multi-resolution investigation unavailable from any flat typing scheme:
  - L6: Same plasmid? (outbreak confirmation)
  - L5: Same clone complex? (ward-level spread)
  - L4-L3: Same sublineage? (regional epidemic context)

---

## 5. Largest Plasmid Classification Resource (79,305 Plasmids)

- Scaled from 8,056 unique training to **79,305 plasmids** (+ 71,249 PLSDB/NCBI RefSeq references)
- Resolves **57,886 unique pLIN codes**
- KNN classifier achieves **97.3% classification rate** (91.1% accuracy across 28 groups, k=5, cosine distance)
- Per-Inc-group clustering reduces peak memory from ~24 GB to **2.3 GB**
- Completes in **<30 minutes** on a standard laptop -- no HPC/cloud required
- Provides foundation for prospective global plasmid surveillance

---

## 6. Two-Tier Outbreak Detection With Clinical Risk Stratification

**First automated plasmid outbreak detection system.** No competing tool offers this.

### Tier 1 -- Basic Module
- Identifies groups of >=2 plasmids sharing both identical L6 pLIN code AND identical AMR fingerprint
- Dual criteria minimise false positives: same backbone but different AMR cargo = independent acquisitions, not outbreak

### Tier 2 -- Temporal Module
- Adds collection date windowing (default 30 days) and location metadata
- Three-tier risk classification:

| Risk | Criteria | Action |
|------|----------|--------|
| **CRITICAL** | >=3 AMR genes AND <=7 days apart | Immediate enhanced isolation + contact tracing |
| **HIGH** | >=3 AMR genes OR <=7 days apart | Infection control review + enhanced screening |
| **MODERATE** | Same L6 + same AMR within 30 days | Surveillance monitoring + antibiogram review |

### SNP Sub-typing for Definitive Confirmation
- 0 SNPs: potentially clonal (strongest evidence for direct transmission)
- 1-5 SNPs: highly related (recent divergence during transmission chain)

---

## 7. Permanent Code Stability for Longitudinal Surveillance

- pLIN codes are **permanent** -- once assigned, they never change regardless of future database additions
- Inherited from the LIN nearest-neighbour assignment rule
- Solves a critical barrier to multi-centre surveillance: hospitals in different regions can directly compare lineage codes without reconciling database versions
- MOB-suite and COPLA reassign codes with every update, making longitudinal tracking impossible
- Essential for WHO GLASS and national reference laboratory networks

---

## 8. Co-Selection of Virulence and AMR on the Same Plasmid

- IncFII plasmids (66.1% of GN dataset) carry **both** the highest virulence factor rate (23.2%) and substantial AMR (62.3%)
- Virulence factors include Salmonella spv genes, iron acquisition systems (iucA, iutA), and serum resistance (traT)
- Antibiotic pressure co-selects for virulence: prescribing any antibiotic to which the plasmid carries resistance simultaneously selects for virulence factors on the same replicon
- IncX1 plasmids carry alpha-haemolysin (hlyA 48.2%) alongside AMR genes -- tissue damage + resistance

---

## 9. Colistin Resistance Lineage Tracking

- **204 mcr-positive plasmids** mapped to specific pLIN lineages
- Predominance of mcr-1.1 (n=83) on specific lineages suggests **clonal expansion** rather than independent emergence
- **18 plasmids** co-carry mcr + carbapenemase genes -- the most therapeutically dire scenario (pan-drug resistance with no reliable treatment)
- Lineage-level tracking enables targeted surveillance of colistin resistance dissemination

---

## 10. Accessible Open-Source Clinical Tool

- Freely available as an interactive **Streamlit web application**
- Requires **no bioinformatics expertise** to operate
- Runs on **standard hardware** (laptop) in **<30 minutes**
- Complete pipeline: FASTA input --> Inc classification --> pLIN assignment --> AMR profiling --> outbreak detection --> clinical reporting
- Open source under GPL-3.0 licence
- Immediately deployable in clinical microbiology laboratories worldwide

---

## Summary of What pLIN Uniquely Reveals

| What pLIN shows | What existing tools show |
|----------------|------------------------|
| pLIN 1.1.2.15.48.671 (n=90): 100% KPC-positive, 13.2 AMR genes, immediate isolation trigger | "IncN" (shared by >1,000 plasmids with vastly different profiles) |
| Cross-Inc resistance hub exchanging carbapenemase + pan-aminoglycoside resistance | Individual replicon types with no cross-Inc linkage |
| CRITICAL outbreak alert: 5 ICU patients, same plasmid, same AMR, within 7 days | No outbreak detection capability |
| Same pLIN code in two countries, 2 years apart -- international lineage spread confirmed | Codes changed with database updates -- comparison impossible |
| 39.7% of AMR+ plasmids carry resistance to 5+ drug classes in a single transferable unit | Individual gene detections without plasmid-level integration |
