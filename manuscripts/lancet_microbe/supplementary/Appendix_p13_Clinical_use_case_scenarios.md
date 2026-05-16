# Appendix p 13: Clinical use case scenarios

## Scenario 1: ICU carbapenemase outbreak investigation

A hospital clinical microbiology laboratory sequences plasmids from five CRE isolates detected in an ICU over two weeks. Using pLIN:

1. Upload FASTA sequences → pLIN assigns codes within minutes
2. All five plasmids receive pLIN code 1.1.2.15.48.671 (L6 identity confirmed)
3. AMRFinderPlus detects blaKPC-2 on all five (identical AMR fingerprint)
4. Temporal module: 5 isolates within 14 days → CRITICAL risk
5. SNP sub-typing: 0–2 SNPs between pairs → clonal spread confirmed
6. Action: Enhanced contact precautions, environmental screening, antibiotic stewardship review

**Without pLIN:** Laboratory reports "IncN plasmid with KPC" for each isolate. Whether these represent the same plasmid (clonal spread requiring isolation) or five independent IncN-KPC acquisitions (requiring stewardship review) cannot be determined without specialist bioinformatics analysis.

## Scenario 2: Regional colistin resistance surveillance

A public health reference laboratory receives plasmid sequences from MCR-positive isolates across 12 hospitals over six months. Using pLIN:

1. Assign pLIN codes to all MCR-positive plasmids
2. Three hospitals share plasmids with identical L6 codes and mcr-1.1 → same lineage circulating regionally
3. Two other hospitals have plasmids with different L6 codes but same L4 prefix → related but distinct sublineages
4. Remaining hospitals have unrelated MCR plasmids → independent acquisitions
5. Hierarchical zoom-out: L4 analysis reveals all related plasmids belong to a single sublineage expanding since month 2

**Without pLIN:** All isolates reported as "MCR-positive". The regional spread pattern and distinction between clonal expansion vs independent emergence is invisible.

## Scenario 3: Longitudinal surveillance across database updates

- **Year 1:** Hospital A assigns pLIN code 4.11.83.485.1175.3594 to an IncFII plasmid.
- **Year 3:** Hospital B (different country) detects a plasmid with identical L6 code.
- **Interpretation:** Same strain-level lineage detected in two countries — possible international spread.

**Key property:** The pLIN code assigned in Year 1 has NOT changed despite two years of database updates (71,249 new reference sequences added). This longitudinal stability is impossible with MOB-suite (cluster codes reassigned with each database release) or COPLA (PTU codes recomputed with each update).
