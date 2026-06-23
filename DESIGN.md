# Covalent Early-Warning Radar — Design

## Goal

**Primary — covalent-move early warning.** The moment a credible player makes a *specifically covalent*
move — a **patent**, a **Phase 1**, or a **pre-Phase-1** signal — alert us, with the chemistry where we
have it, so we can **fast-follow**.

**Secondary — validated-target watchlist.** Well-documented targets with many programs across phases
(RAF1-like) where a covalent **best-in-class** is worth a look.

## Principles

1. **Two legs, opposite timing, cross-linked.** Patents lag ~18–34 months but carry the *chemistry*;
   clinical/PR is the *earliest* signal but carries *no structure*. Run both; link by company + target.
2. **Deterministic core, bounded LLM.** Everything verifiable is code. The LLM only *judges* ambiguous
   cases, and the UniProt/structural fact-check *overrides* it. (~85% of the engine is deterministic.)
3. **Target-centric, not company-centric.** Monitor *targets* + warhead keywords and let the chemistry
   surface the company — a 50-company list misses stealth Series-A biotechs, which is where disruptive
   covalent chemistry often originates.

## Architecture

```
LEG 1 — PATENT (late, has chemistry)        LEG 2 — CLINICAL (earliest, no chemistry)
 USPTO/EPO/WIPO bulk XML                      ClinicalTrials.gov diff (new Ph1 / Early-Ph1)
   ├ claims-text tripwire (named warheads)      covalency by:
   └ Examples resolver → SMILES                   (a) patent cross-link
       (InChI→RDKit | IUPAC→OPSIN | img→DECIMER)  (b) live web search
   → warhead SMARTS on EXEMPLIFIED cpds            (c) clinical-design fingerprint
        │                                              │
        └───────────────┬──────────────────────────────┘
                        ▼
         UNIFY → one covalent-move log (cross-linked, earliest-dated) → ALERT digest

  DISCOVERY (target-centric): "<target>" × (covalent|irreversible|warhead|TPD) → surfaces the company
  WATCHLIST (structural): PDB ensemble → SASA → PROPKA → fpocket(cryptic) → Covalent-BIC priority
```

## Leg 1 — Patent (Day-1, claims-focused, structure-aware)

- **Sources:** USPTO / EPO / WIPO **bulk XML**, at publication. (SureChEMBL becomes a *later* cross-check,
  not the source — its structure extraction has its own lag.)
- **Tripwire (text, instant):** parse XML → isolate **independent claims** (no dependency back-reference) →
  regex a warhead lexicon (common **and IUPAC** forms, e.g. `acrylamide`/`prop-2-enamide`,
  `propiolamide`/`but-2-ynamide`, `chloroacetamide`, `vinyl sulfonamide`, `epoxide`), **bound to chemical
  syntax** (preceded by locants/brackets/hyphens, e.g. `\]\-acryl`, `\-[0-9a-zA-Z]*\-amide`), with an
  **exclusion list** (`gel|polymer|excipient|column|purification`). → catches *named* warheads.
- **Resolver (structure):** locate the **Examples / Biological-Data / Embodiments** section → extract each
  exemplified compound in whichever of three formats it appears, preferring the most reliable:
  1. embedded **SMILES/InChI** → RDKit (exact)
  2. **IUPAC name** → **OPSIN** (deterministic — the default; born-digital text, ~no error)
  3. **2D image** → **DECIMER** OCSR (ML — fallback, only when image-only)
  → unified SMILES list of the *exemplified* compounds. → catches *drawn* warheads (which claims text misses).
- **Confirm:** run warhead **SMARTS** on those SMILES; require the warhead on a **meaningful fraction** of
  the examples (not a fringe R-group in a broad Markush); tie to **IC50/Ki** from the potency table where present.
- **Output:** Day-1 covalent-patent event, with the lead structures.

## Leg 2 — Clinical (covalency without the chemistry)

- **Trigger:** ClinicalTrials.gov diff → new **Phase 1 / Early-Phase-1** interventional trials.
- **Covalency classification (stacked, cheapest first):**
  1. **Patent cross-link** — does the sponsor have a covalent patent on the same target? (most reliable)
  2. **Live web** — search the code name (NOT an LLM-from-memory lookup; new code names aren't in training data)
  3. **Clinical-design fingerprint** — the design betrays covalency with zero chemistry:
     - target-**occupancy time-course** ("% occupancy at 24/48 h") — covalent occupancy outlasts drug PK
     - PD **sustained past plasma clearance**
     - **PBMC mass-spec** / covalent-adduct assays
     - dosing tied to **protein resynthesis rate**, not drug half-life
     → weighted classifier; the *combination* is highly specific (occupancy alone is not).

## Discovery — target-centric (surfaces stealth companies)

Per target of interest: search `"<target>" AND (covalent OR irreversible OR warhead OR "targeted protein
degrader")` across news/PR, bioRxiv, conference abstracts (AACR/ASCO), CT.gov. **The company is the output**;
credibility is a filter on what surfaces, not the search key.

## Unify + Alert

One **covalent-move event log** — deduped, **cross-linked by company + target**, dated by the *earliest*
signal, carrying the chemistry once the patent lands. Prioritize by stage-earliness × credibility ×
Axiom-relevance. → digest (email/Slack).

## Watchlist — validated *and* structurally covalent-tractable

For each heavily-pursued target (high clinical/literature count):
1. **RCSB PDB API** → all experimental structures, esp. **holo** (ligand-bound) = a conformational ensemble.
   (No PDB → AlphaFold single-model fallback, weaker.)
2. **SASA** (Biopython/MDAnalysis) per cysteine across the ensemble — buried in *all* states → drop;
   **spikes in some holo state → conformationally-dependent opportunity** (the RAF1 pattern).
3. **PROPKA** → flag cysteines with predicted **pKa ≪ 8.3** (depressed = reactive thiolate).
4. **fpocket** on apo vs holo → cysteine on a druggable pocket, especially a **cryptic** one (present in
   holo, absent in apo). (CryptoSite = heavier ML alternative.)
- **Score:** heavily-pursued × [exposed in some state] + [low pKa] + [on a (cryptic) pocket] → top Covalent-BIC.

## Determinism & reproducibility

- **Deterministic (no LLM):** patent parsing/regex/resolver/SMARTS, CT.gov diff, the clinical-design
  fingerprint, the structural stack, the UniProt fact-check.
- **LLM (bounded, cached, fact-check-overridden):** covalency judgment on ambiguous cases; watchlist prose.
- **Pin:** data-source snapshot dates, tool versions; write a per-run `manifest.json`; pull knobs (date
  cutoff, warhead lexicon, SMARTS, fraction threshold, SASA/pKa cutoffs) into one `config.yaml`.

## Build order

| # | Build | Why first |
|---|---|---|
| 0+1 | Patent **claims-text tripwire** + **Examples resolver** (OPSIN→DECIMER) + SMARTS | 80% exists; gives a Day-1 patent radar |
| 2 | **Clinical-design fingerprint** (CT.gov text only) | highest leverage, needs no chemistry |
| 3 | **Unify + alert** | turns signals into a product |
| 4 | **Target-centric discovery** | catches stealth biotechs |
| 5 | **Structural watchlist** (PDB ensemble · SASA · PROPKA · fpocket) | the secondary BIC list |

## Honest limits

- Claims text often **draws** the warhead (not names it) → OCSR is the necessary complement; OCSR is
  error-prone but **tolerable for binary warhead detection** (you only need the reactive group recognized).
- **OPSIN** fails on malformed names; a single patent **mixes formats** → all three converters are needed.
- **Occupancy assays aren't covalent-exclusive** → rely on the *combination* of fingerprints.
- The structural stack is **heuristic triage** (SASA/PROPKA/fpocket), **necessary-not-sufficient** → a
  flagged cysteine still needs chemoproteomic/experimental confirmation (the RAF1 lesson).
- **Patents can't beat their own ~18-month publication lag** — the clinical/PR leg is the only sub-18-month signal.
