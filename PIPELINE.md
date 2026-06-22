# Layer 1 — Target Radar: Detailed Pipeline Description

The end-to-end data flow for Layer 1, stage by stage, with data shapes, rules, cross-cutting
guarantees, and failure handling. Companion to the [Layer 1 course](../layer1-course/README.md).

---

## 0. Two kinds of source (the distinction that organizes everything)

- **Reference sources** — HGNC, UniProt, Ensembl/MANE, VRS/SPDI tooling. They build the
  **identity index**; they do *not* produce claims. Refreshed on database release.
- **Evidence sources** — ClinicalTrials.gov, openFDA, patents (SureChEMBL), financings, genetics.
  They produce **claims** about targets. Polled incrementally.

The identity index is built first and consulted by every evidence connector. Everything below
assumes that split.

```
REFERENCE ──▶ identity index ──┐
                               ▼ (consulted during normalization)
EVIDENCE ──▶ extract ──▶ normalize+resolve ──▶ observed claims ──▶ resolver ──▶ canonical graph ──▶ query ──▶ actionable targets
                │                                     │                              │            │
             raw store                          (append-only log)               (projection)  corroboration+timing
```

---

## 1. Extraction (connectors)

One connector per source; the only code that knows a source's wire format.

**Responsibilities:** auth · pagination · rate-limit + exponential backoff · incremental cursor ·
fetch raw · write raw to the **raw store** · pass raw record downstream.

**Per-source cursor (incremental refresh):** CT.gov → `lastUpdatePostDate`; openFDA → label
effective date; patents → publication date. The connector pulls only records newer than the last
successful cursor, never the whole corpus.

**Output:** raw records keyed `(source, source_record_id, retrieved_at, content_hash)`. The raw
store is immutable — kept for replay and audit, never mutated.

**Idempotency:** if `content_hash` is already seen for that `source_record_id`, the record is a
no-op (re-running a connector never duplicates).

---

## 2. Normalization + identity resolution (the transform core)

Two sub-steps per raw record.

**2a. Field extraction.** A per-source field map pulls the structured facts out of the source's
idiosyncratic shape:

- candidate **target mention(s)** (e.g., CT.gov `conditionsModule.conditions` + intervention/
  keyword text; openFDA `openfda.*` + indication text)
- **actor** (sponsor / patent assignee)
- **dates** (first-post, filing, approval)
- **evidence payload** (NCT id, document URL, label section)

**2b. Identity resolution** (enforces "never join on a name"):

1. Normalize the mention string (uppercase, strip).
2. Look up in the identity index across `{symbol, prev_symbol, alias_symbol, xref}`.
3. Resolve through the ID **lifecycle** (withdrawn/merged/split → current).
4. If the mention is a variant (e.g., "HER2 A775_G776insYVMA"), parse to an `allele_exact` via
   HGVS→VRS/SPDI; otherwise resolve to the `gene` tier.
5. Return a canonical `EntityRef`.

**Outcomes are explicit, never guessed:** `resolved` | `ambiguous` (mention maps to >1 node →
flagged, not auto-picked) | `unresolved` (logged; a spike in the unresolved rate is the
schema-drift alarm).

---

## 3. Claim construction

Each resolved mention becomes one **observed claim** — an immutable, source-faithful assertion:

```
observed_claim {
  subject        : EntityRef            # what it's about (canonical, typed)
  predicate      : has_trial | is_approved_for | has_patent | has_genetic_assoc | ...
  object         : EntityRef?           # optional (e.g., disease/indication)
  evidence_type  : trial | approval | patent | financing | genetics   # the sensor / latent variable
  source         : "clinicaltrials.gov"
  primary_knowledge_source : "...",  evidence : {doc_id, url}
  actor          : sponsor/assignee
  valid_time     : 2009-05-12          # when true in the world (first-post / filing / approval)
  observed_at    : 2026-06-11          # when we recorded it (transaction time)
  source_release : "ctgov-2026-06-11"
  coordinate_context : {assembly/transcript}?   # only for allele claims
  confidence     : 0.0–1.0
  claim_key      : hash(source, source_record_id, predicate, subject)   # idempotency
}
```

Append-only. The same trial seen again with no change re-hashes to the same `claim_key` → no-op;
a changed trial appends a *new* claim (history preserved).

---

## 4. Persistence (four stores)

| Store | Contents | Mutability |
|---|---|---|
| **Raw store** | verbatim API payloads | immutable (replay/audit) |
| **Identity index** | reference-derived name/xref → canonical node + lifecycle | rebuilt per reference release |
| **Observed-claim log** | all observed claims | append-only (source of truth) |
| **Canonical store** | the reconciled graph | derived, recomputable |

Tech: DuckDB/SQLite (Postgres at scale) for claims; canonical graph as a typed nodes+edges table
(graph DB optional later).

---

## 5. Resolution (observed → canonical)

The **resolver** projects the observed-claim log into the canonical graph:

- **Node canonicalization** — group claims whose subjects resolve to the same `primary_id` into
  one canonical node (synonyms already collapsed by Stage 2).
- **Edge construction** — claims become typed edges/attributes on nodes; identity edges
  (`encodes`, `variant_of`) come from reference data.
- **Dedup + conflict** — collapse duplicate assertions; record disagreements.
- **Confidence aggregation** — v1 trivial: a fact is canonical if backed by claims above a
  confidence threshold.

**Recomputable:** canonical = pure function of (observed log, resolver version). Drop it and
replay the log to rebuild; bump the resolver and recompute with full history intact. The resolver
version is stamped on every canonical record.

---

## 6. Canonical graph

Typed nodes (`gene / protein / protein_isoform / allele_exact` in v1) + typed edges + claim-derived
evidence attached to each node. Each node carries the set of evidence claims pointing at it
directly or reachable by roll-up edges, each with its `evidence_type` and `valid_time`.

---

## 7. Query: corroboration + timing + comparison tier

Not precomputed — evaluated per query so the comparison tier is tunable.

**Input:** `comparison_tier` (default `gene`), time `window`, independence requirement.

**Process per candidate node:**

1. **Roll up** claims to the comparison tier along typed identity edges only (never across
   fusion/complex/allele boundaries).
2. **Corroboration:** require ≥2 **independent evidence types** within `window`
   (trial ∧ approval = yes; patent ∧ patent = no).
3. **Recency:** earliest `valid_time` among supporting claims (rewards upstream capture).
4. **Stage/crowding:** derive from the claim mix (approved vs phase-1; sponsor count).

**Output:** ranked actionable targets, each with its corroborating claim set, earliest-signal
date, and full provenance.

---

## 8. Handoff

Each actionable target → Layers 2–6 (validation, differentiation, FTO, covalent feasibility,
scoring), carrying its corroborating claims, earliest-signal date, and provenance, so downstream
layers never re-fetch.

---

## Cross-cutting guarantees

- **Idempotent + incremental** ingestion (content hash + per-source cursor).
- **Provenance / audit** — every canonical fact → its claims; every resolved ID → its lifecycle
  record.
- **Recomputability** — canonical is disposable; the observed log is the truth.
- **Bitemporality** — `valid_time` vs `observed_at` on every claim → point-in-time reconstruction.
- **Observability** — per-run counts (pulled / normalized / resolved / **unresolved** / ambiguous);
  a rising unresolved rate = schema-drift alarm.
- **Isolation** — one connector failing doesn't corrupt others or the graph; partial runs are safe.
- **Versioning** — connector output schema and resolver both versioned and stamped.

---

## Failure modes → handling

| Failure | Handling |
|---|---|
| Name mismatch (ERBB2 vs HER2) | resolved to same `HGNC:3430` via alias index (Stage 2) |
| Source schema drift | unresolved-rate alarm; only that connector affected (Stage 1/2) |
| Ambiguous alias | flagged `ambiguous`, never auto-picked (Stage 2) |
| Allele vs class conflation | typed nodes + hard boundaries; no cross-tier equate (Stage 6–7) |
| Duplicate ingestion | `claim_key` / `content_hash` no-op (Stage 1/3) |
| Source outage | last good cursor retained; resume on next run |
| Bad source discovered later | down-weight + recompute canonical from the log (Stage 5) |

---

## End-to-end trace (one signal)

BI's 2021 HER2 patent and Iambic's later trial both **normalize → resolve to `HGNC:3430`** (names
never used as the key) → two observed claims (`evidence_type` = patent, trial) → the resolver
attaches both to the HER2 gene node → a query at the gene tier finds **two independent evidence
types within the window → corroborated**, recency scored at **2021** → HER2 is emitted as an
actionable target to Layers 2–6.

---

## v1 scope (what this pipeline builds first)

- **Nodes:** `gene`, `protein`, `protein_isoform`, `allele_exact`.
- **Reference sources:** HGNC (+ UniProt, MANE for xrefs).
- **Evidence sources:** ClinicalTrials.gov, openFDA. (Patents/SureChEMBL next; they are the real
  upstream signal but need a heavier ingestion path.)
- **Resolver:** trivial (union above a confidence threshold).
- **Query:** corroboration at the `gene` tier, recency on earliest `valid_time`.

Deferred but modelled: `proteoform`, `allele_class`, `complex`, `fusion`, `pathway`,
`stateful_entity`, `ncRNA`; sophisticated resolver; full comparison-tier engine; patent OCSR/Markush.
