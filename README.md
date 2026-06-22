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

## Weekly run — no API key (rides your Claude Max subscription)

The LLM judgment runs through Claude Code headless (`claude -p`), authenticated by your Max
subscription, so no Anthropic API key is needed for personal use:

- `llm_claude_code.py` — `judge(evidence, schema)` shells out to `claude -p --output-format json`.
  Swap this one file for an API-key client when productizing for customers (same signature).
- `assemble_evidence.py` — deterministic per-target evidence bundle (handles, fact-check, trials, patents, momentum).
- `run_weekly.py` — discovery → evidence → `claude -p` judgment → fact-check override → rank → snapshot + weekly delta.
- `report_gen.py` — renders the **persuasion report** (`covalent_report_<date>.md`): per target **who's moving** (patents+programs)
  · **why now** (momentum+biology) · **sources/IDs** (linked). Two tiers — open opportunities + a fresh-credible-filings
  fast-follow tier (degrader-filtered, title-validated "Company X just filed on target X").

```bash
uv run --python 3.12 python run_weekly.py
```

Schedule it (Mondays 8am) via `crontab -e`:
```
0 8 * * 1 cd /path/to/layer1 && ~/.local/bin/uv run --python 3.12 python run_weekly.py >> runs/cron.log 2>&1
```

The deterministic fact-check (`factcheck.py`) overrides the model on any residue/selectivity fact.
For a multi-customer product, swap `claude -p` for an Anthropic/OpenAI API key (console.anthropic.com) —
the subscription path is for individual use and is rate-limited.

## Note
`.epo_creds.json` (EPO OPS keys) is gitignored — never commit it.
