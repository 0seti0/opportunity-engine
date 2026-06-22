# Covalent best-in-class opportunity engine (Layer 1)

Mines public data for **covalently-druggable, validated, covalent-open** drug targets and ranks them
by structural tractability + selectivity + lane status + momentum.

## Pipeline

| Stage | Script | What it does | Sources | LLM? |
|---|---|---|---|---|
| Covalent detection | *(DuckDB+RDKit)* | warhead-from-chemistry SMARTS on patent compounds → `covalent_warhead_feed_clean.json` | SureChEMBL bulk | no |
| Target attribution | `attribute.py` | patents → biomedical_entities → official HGNC symbols (validated, de-noised) | SureChEMBL, HGNC | no |
| Patent radar / crowding | `window.py` | distinct-assignee crowding + clinical maturity → OPEN/CONTESTED/CLOSED | SureChEMBL, ClinicalTrials.gov | no |
| Structural opportunity | `cysdb.py` | CysDB ligandable cysteine × FDA-validated × covalent-open × non-essential | CysDB, DepMap, CEG2 | no |
| Pocket cysteines | `structural_cys.py` | KLIFS kinase hinge/P-loop/front-pocket cysteines (recovers CysDB scout-fragment blind spots) | KLIFS | no |
| Emerging momentum | `momentum.py` | literature acceleration (2025-26 vs 2021-23) + preprints + new trials | Europe PMC, ClinicalTrials.gov | no |
| **Fact-check** | `factcheck.py` | **deterministic** residue + paralog-conservation check vs UniProt (catches LLM hallucinations) | UniProt, Biopython | no |
| Report | `report.py` | per-target evidence + provenance → `opportunity_report.json` | (assembles above) | no |
| Verification / ranking | *(Claude Workflows)* | adversarial web-grounded verification, deep analysis, ranking | web search | **yes** |

The `.py` stages are **deterministic and run standalone**. The verification/deep-analysis/ranking steps
are LLM-orchestrated (Claude agents with web search) and are not standalone.

## Run

```bash
uv run --with duckdb --with rdkit --python 3.12 python attribute.py
uv run --with duckdb --python 3.12 python window.py
uv run --with openpyxl --python 3.12 python cysdb.py
uv run --with openpyxl --python 3.12 python structural_cys.py
uv run --with openpyxl --python 3.12 python momentum.py
uv run --with biopython --python 3.12 python factcheck.py
uv run --with openpyxl --python 3.12 python report.py
```

## Data (not in repo — download separately)
- SureChEMBL bulk parquet: https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/
- CysDB supplement (NIHMS1893018-supplement-2.xlsx): from the CysDB paper (Cell Chem Biol 2023)
- HGNC complete set, DepMap common essentials, Hart CEG2 — see script headers

## Note
`.epo_creds.json` (EPO OPS keys) is gitignored — never commit it.
