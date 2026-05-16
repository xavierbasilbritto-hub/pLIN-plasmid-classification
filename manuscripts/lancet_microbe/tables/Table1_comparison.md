# Table 1: Comparative evaluation of plasmid classification systems

| Feature | PlasmidFinder | pMLST | MOB-suite | COPLA | mge-cluster | pLIN |
|---------|:---:|:---:|:---:|:---:|:---:|:---:|
| **Classification approach** | Replicon typing | Allelic profiling | Relaxase clustering | Host-range + mobility | Reference-free k-mer | Hierarchical 4-mer |
| **Resolution levels** | 1 (Inc group) | 1 (sequence type) | 1 (cluster) | 1 (PTU) | 1 (cluster) | 6 (L1--L6) |
| **Multi-resolution hierarchy** | No | No | No | Partial | No | Yes |
| **Stable nomenclature** | Yes | Yes | No* | No* | No | Yes |
| **Discriminatory power (Simpson's D)** | 0.641 | NA | NA | NA | NA | 0.985 |
| **Integrated AMR profiling** | No | No | No | No | No | Yes |
| **Automated outbreak detection** | No | No | No | No | No | Yes |
| **Clinical risk stratification** | No | No | No | No | No | Yes |
| **Reference-free operation** | No | No | No | No | Yes | Yes** |
| **Open-source availability** | Yes | Yes | Yes | Yes | Yes | Yes |
| **No. Inc/Rep groups supported** | >30 | ~10 | NA | NA | NA | 28 |
| **Scalability (>50,000 plasmids)** | Yes | Yes | Yes | Yes | Yes | Yes |

\* MOB-suite and COPLA reassign cluster identifiers with each database update, precluding longitudinal surveillance.

\*\* pLIN operates reference-free for 4-mer computation and hierarchical clustering; the Inc group classifier requires the reference database.

NA = not applicable or not reported.
