# Appendix p 11: Outbreak detection algorithm — detailed specification

## Basic outbreak detection module

```
INPUT:  pLIN results table (pLIN code, AMR genes per plasmid)
OUTPUT: Outbreak clusters with risk levels

FOR each unique L6 pLIN code:
    GROUP plasmids sharing this L6 code
    FOR each unique AMR fingerprint within the group:
        IF count >= 2:
            CREATE outbreak cluster
            IF number of shared AMR genes >= 3:
                ASSIGN risk = HIGH
            ELSE:
                ASSIGN risk = MODERATE
```

## Temporal outbreak detection module

```
INPUT:  pLIN results table + metadata with collection dates
OUTPUT: Time-windowed outbreak clusters with 3-tier risk

FOR each unique L6 pLIN code:
    GROUP plasmids sharing this L6 code
    FOR each unique AMR fingerprint within the group:
        IF count >= 2:
            SORT by collection date
            FOR each pair within time_window (default 30 days):
                CREATE temporal cluster
                days_apart = |date_i - date_j|
                n_amr = count of shared AMR genes
                IF n_amr >= 3 AND days_apart <= 7:
                    ASSIGN risk = CRITICAL
                ELIF n_amr >= 3 OR days_apart <= 7:
                    ASSIGN risk = HIGH
                ELSE:
                    ASSIGN risk = MODERATE
```

## Design rationale for dual criteria (L6 identity + AMR fingerprint)

Plasmids sharing the same L6 code (cosine distance <= 0.001, ~99.9% ANI) but carrying different AMR gene complements likely represent independent acquisitions of resistance elements onto a common plasmid backbone — rather than clonal spread of a single AMR-carrying plasmid. By requiring both compositional identity AND functional (AMR) identity, false positive outbreak alerts are minimised.

## SNP sub-typing for outbreak confirmation

Within flagged L6 clusters, minimap2 alignment (-cx asm5, --cs tag) provides nucleotide-level resolution:
- **0 SNPs:** Potentially clonal; strongest evidence for direct transmission
- **1–5 SNPs:** Highly related; consistent with recent divergence during transmission chain
- **6–20 SNPs:** Related but divergent; possible indirect transmission or shared source
- **>20 SNPs:** Distinct within the strain-level cluster; independent events
