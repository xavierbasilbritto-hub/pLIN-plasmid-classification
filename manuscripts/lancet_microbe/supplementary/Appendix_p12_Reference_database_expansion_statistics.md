# Appendix p 12: Reference database expansion — detailed statistics

## KNN Inc type classification performance

| Metric | Value |
|--------|-------|
| Total reference sequences | 71,249 |
| Classified (confidence >=40%) | 69,140 (97.0%) |
| Unknown/Novel (confidence <40%) | 2,109 (3.0%) |
| Multi-replicon flagged (2+ Inc types >25%) | 10,703 |
| Mean classification confidence | 81.3% |
| High confidence (>=95%) | 31,122 (43.3%) |
| Moderate confidence (80–95%) | 9,855 (13.7%) |
| Low-moderate confidence (40–80%) | 28,588 (39.8%) |

## Per-Inc-group clustering results (full 79,305 plasmids)

| Inc group | Training | Reference | Total | L6 pLIN codes | Distance matrix size |
|-----------|----------|-----------|-------|---------------|---------------------|
| IncFII | 4,629 | 29,407 | 34,036 | 14,403 | 2.3 GB (chunked) |
| IncX1 | 705 | 25,449 | 26,154 | 14,060 | 1.4 GB (chunked) |
| IncN | 1,097 | 4,206 | 5,303 | 1,912 | 56 MB |
| ColRNAI | 91 | 1,836 | 1,927 | 618 | 7.4 MB |
| IncX3 | 56 | 1,403 | 1,459 | 323 | 4.3 MB |
| IncFIB | 97 | 1,022 | 1,119 | 279 | 2.5 MB |
| IncI1 | 27 | 920 | 947 | 82 | 1.8 MB |
| IncI2 | 25 | 851 | 876 | 229 | 1.5 MB |
| IncHI1 | 16 | 697 | 713 | 115 | 1.0 MB |
| IncAC2 | 14 | 693 | 707 | 262 | 1.0 MB |
| IncX4 | 24 | 669 | 693 | 195 | 0.96 MB |
| IncHI2 | 36 | 607 | 643 | 152 | 0.82 MB |
| IncF | 75 | 412 | 487 | 120 | 0.47 MB |
| IncA | 14 | 331 | 345 | 188 | 0.24 MB |
| IncC | 16 | 324 | 340 | 107 | 0.23 MB |
| ColE | 19 | 307 | 326 | 74 | 0.21 MB |
| IncI | 11 | 305 | 316 | 60 | 0.20 MB |
| IncR | 21 | 70 | 91 | 67 | 0.02 MB |
| IncFIC | 14 | 28 | 42 | 22 | <0.01 MB |
| IncFIBK | 11 | 14 | 25 | 8 | <0.01 MB |

## Runtime (Apple M-series laptop, single-threaded)

- Phase 1 (4-mer computation, 71,249 sequences): 23.5 min (~51 sequences/sec)
- Phase 2 (KNN classification): 7 sec
- Phase 3 (per-group clustering, 20 Gram-negative groups): 4.2 min
- Phase 4 (output): <1 sec
- **Total: 28 minutes**
