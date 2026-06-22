# Target Radar (Layer 1) — Technical Approach & Architecture

**Audience:** a senior engineer evaluating the design cold. No prior context assumed.
**Status:** greenfield. A throwaway probe (Phase 0) empirically retired the single biggest
risk — that public data can't yield targets — but none of its code is load-bearing; this
document designs the durable system from scratch.

---

## 1. What it does, and why it's hard

**Goal.** Continuously ingest public biomedical sources — clinical-trial registries
(ClinicalTrials.gov), drug approvals/labels (openFDA), later patents — and maintain a
queryable knowledge graph of **corroborated, time-stamped, target-level claims**. On query,
return a ranked list of drug targets *gaining momentum*, e.g. *"HER2 — corroborated by an
active trial and an approval, earliest signal 2021,"* each with full provenance.

**Why it's not a CRUD app.** It's three well-understood hard problems composed, plus a domain wrinkle:

1. **Entity resolution / master data management.** The same target appears as `HER2`,
   `ERBB2`, `NEU`, `c-erbB2`; the same drug as `doxorubicin`, `DOXORUBICIN HYDROCHLORIDE`,
   `Adriamycin`. Signals only corroborate if these collapse to one canonical identity. The
   governing rule is **never join on a name.**
2. **Event sourcing.** Sources revise records, get re-scored, and get superseded; we must
   keep history and reprocess it. So the source of truth is an **append-only log of immutable
   claims**, and the queryable graph is a **derived, recomputable projection**.
3. **Bitemporal modeling.** Ranking on "earliest signal" and answering "what did we know, and
   when" requires tracking two timelines per claim: when it was *true in the world* and when we
   *recorded it*.

**The domain wrinkle.** The target is **latent**: these sources are keyed on *disease* and
*drug*, never on the molecular target. So extraction must *recover* the target two ways — from
biomarker text, and by mapping the drug to its target via a reference database. The Phase 0
probe measured this: the drug→target path was ~100% precise on an independently checkable set;
the binding constraint is **recall** (drug-name resolution), not correctness. That result shapes
the whole design.

---

## 2. Design goals (the invariants) and non-goals

These are the properties that justify the architecture. Each is expensive to retrofit, so they're
adopted from day one.

- **Idempotent + incremental ingestion** — re-running a connector never duplicates; each run
  pulls only what changed (content hashing + per-source cursor).
- **Provenance / audit** — every canonical fact resolves to the claims behind it; every resolved
  ID resolves to its lifecycle record. No orphan assertions.
- **Recomputability** — `canonical = f(observed_log, resolver_version)`. The graph is disposable;
  drop it and replay the log to rebuild, bit-for-bit.
- **Bitemporality** — `valid_time` vs `observed_at` on every claim → point-in-time reconstruction.
- **Never join on a name** — every cross-record join goes through a canonical ID.
- **Isolation** — one connector failing degrades only its own feed; partial runs are safe.
- **Observability** — per-run counts (pulled / normalized / resolved / **unresolved** / ambiguous);
  a rising unresolved rate is the schema-drift alarm.
- **Versioning** — connector output schema and resolver are both versioned and stamped on outputs.

**Non-goals (v1).** Not a general biomedical KG (Open Targets exists for that); not real-time
(hourly/daily batch is fine); not ML-driven ranking (deterministic, explainable scoring first).
Deferred: patents, sophisticated resolver tiers, the full comparison-tier engine, a graph DB.

---

## 3. Architecture overview

The organizing distinction is **two kinds of source**:

- **Reference sources** (HGNC, UniProt, MANE; ChEMBL/GtoPdb for drugs) build the **identity
  index**. They do not produce claims. Refreshed per database *release*.
- **Evidence sources** (ClinicalTrials.gov, openFDA, patents) produce **claims**. Polled
  *incrementally*.

```
REFERENCE ──▶ identity index ──┐  (consulted during resolution)
                               ▼
EVIDENCE ──▶ connector ──▶ extract ──▶ resolve ──▶ observed claims ──▶ resolver ──▶ canonical ──▶ query
              │              (2a)        (2b)      (append-only log)   (projection)   graph     (corroboration
           raw store                                  = source of truth   (pure fn)             + timing)
                                                                                                    │
                                                                                          ranked targets + provenance
```

Pattern-wise this is **event sourcing + CQRS** (the claim log is the event store/write model;
the canonical graph is a rebuildable read model) layered over a **master-data-management** hub
(the identity index + resolver) with a **bitemporal** record on every event.

---

## 4. Data model (the core)

Three types carry the system.

**`EntityRef`** — a resolved, typed identity. Produced by the resolver; the only thing other
components are allowed to join on.
```
EntityRef { primary_id   # e.g. "HGNC:3430"
            type         # gene | protein | protein_isoform | allele_exact   (hard boundaries)
            lifecycle    # current | merged_from(..) | withdrawn
            ref_release  # which reference release resolved it (reproducibility)
            tier }       # the comparison tier this id sits at
```
Typing is load-bearing: it prevents conflating an allele-specific signal with the whole gene,
and roll-up happens **only along typed identity edges**, never across allele/fusion/complex
boundaries.

**`ObservedClaim`** — the immutable unit of evidence. Append-only; never updated.
```
ObservedClaim {
  subject        : EntityRef                      # canonical, typed
  predicate      : has_trial | is_approved_for | has_patent | ...
  object         : EntityRef?                      # e.g. disease/indication
  evidence_type  : trial | approval | patent | genetics | financing   # the "sensor"
  source, primary_knowledge_source, evidence{doc_id,url}
  actor          : sponsor | assignee              # for actor-independence in corroboration
  extraction_path: text_biomarker | drug_to_target # how we recovered it (≠ confidence)
  valid_time     : date     # when true in the world (first-post / filing / approval)
  observed_at    : date     # transaction time (when we recorded it)
  source_release : str
  confidence     : 0..1
  claim_key      : hash(source, source_record_id, predicate, subject, content_hash) }
```
Note `claim_key` **includes `content_hash`** (deliberate fix to the original spec): an unchanged
re-observation re-hashes identically → no-op; a *changed* source record produces a new hash →
appends a new claim, preserving history. Without the content hash, history is silently lost.

**Canonical node / edge** — the derived projection. Typed nodes (`gene/protein/isoform/allele`)
+ typed edges (`encodes`, `variant_of` from reference data; evidence edges from claims), each
node carrying the claim set pointing at it. Every canonical record is stamped with
`resolver_version` so it's recomputable and diffable across resolver upgrades.

---

## 5. Components — responsibility, interface, key decisions

### 5.1 Connectors (one per evidence source)
The only code that knows a source's wire format. Responsibilities: auth · pagination ·
rate-limit + exponential backoff · incremental cursor · fetch raw · write raw verbatim · emit
downstream.
- **Incremental cursor:** a per-source high-water mark (`CT.gov.lastUpdatePostDate`,
  openFDA label effective date, patent publication date). At-least-once fetch + idempotent
  writes = effectively-once.
- **Idempotency:** raw records keyed `(source, source_record_id, content_hash)`; a seen hash is
  a no-op.
- **Error discipline (first-class, learned the hard way):** throttle *every* call; retry
  transient failures (429/5xx/network) with backoff honoring `Retry-After`; **never silently
  swallow** — distinguish a *deterministic empty* (legitimately no data) from a *transient
  error* (retryable), because conflating them silently deflates yield and reads as a real gap.

### 5.2 Identity layer (the spine)
Two sub-parts, both reusing the **same machinery** for two entity types:
- **Reference ingestion → identity index.** HGNC (`symbol/prev_symbol/alias_symbol/withdrawn`)
  + UniProt + MANE for gene↔protein xrefs → `{normalized name/xref → canonical node + lifecycle}`,
  rebuilt per release.
- **The resolver** — a *pure function* with explicit outcomes:
  ```
  resolve(mention, context) -> Resolution {
      status: resolved | ambiguous | unresolved
      entity_ref?            # on resolved
      candidates?            # on ambiguous (>1 node) -- flagged, NEVER auto-picked
  }
  ```
  Steps: normalize → look up across `{symbol, prev_symbol, alias_symbol, xref}` → walk the
  lifecycle (withdrawn/merged/split → current) → return. **Precision controls** (from the probe):
  case-sensitive candidate detection (`MET` the gene vs "met" the word), a risk-tier/stoplist
  for collision-prone short symbols, and context-cue requirements for the risky ones.
- **Gold-set gate.** The resolver ships behind a hand-labeled `mention → expected outcome` set
  (synonyms, lifecycle, genuinely-ambiguous-must-flag, and negative controls). This is the
  definition of done; it stops the resolver from being graded on its own output.
- **Two instances:** gene-identity (HGNC) and drug-identity (drug name → canonical drug, the #1
  recall lever — salts/code-names/combos). Same contract, same gate.

### 5.3 Extraction / normalization
Per-source field map (2a) pulls candidate mentions, actor, dates, evidence payload; then
resolution (2b) turns mentions into `EntityRef`s. Two recovery paths:
- **Path A — text biomarker:** scan the record's own text (conditions, **title, eligibility**)
  for a named target. Must distinguish *inclusion/stratification* biomarkers (the trial's target)
  from *exclusion / concomitant-med* mentions ("no prior anti-TNF" is **not** a TNF trial).
- **Path B — drug→target:** resolve the drug (5.2) then map to its target via a **dual
  ChEMBL + GtoPdb** bridge (ChEMBL is precise but recall-limited; GtoPdb, independently curated,
  recovers code names/salts/biologics ChEMBL misses), collapsing **family/complex explosions**
  (a taxane is "tubulin," not 15 `TUBB*` genes; an ADC's target is the antibody's, not the
  payload's).

The two paths capture *different notions of target* — Path A the selection biomarker, Path B the
drug mechanism — so both are kept as **typed, distinct claims**, not merged.

### 5.4 Claim construction
Each resolved mention becomes one `ObservedClaim`, appended to the log. The same record seen
unchanged is a no-op (claim_key); a changed record appends a new claim.

### 5.5 Persistence — four stores
| Store | Contents | Mutability | Tech |
|---|---|---|---|
| Raw store | verbatim API payloads | immutable (replay/audit) | object store / files, content-addressed |
| Identity index | name/xref → node + lifecycle | rebuilt per reference release | DuckDB/SQLite table |
| Observed-claim log | all claims | **append-only (source of truth)** | DuckDB/Parquet (Postgres at scale) |
| Canonical store | reconciled graph | derived, recomputable | typed nodes+edges table (graph DB later) |

Rationale: the graph is small and the dominant access pattern is *analytical replay over the
log*, which a columnar embedded store (DuckDB/Parquet) does extremely well with near-zero ops.
Reach for Postgres when concurrency/scale demands it, and a graph DB only when query traversals
genuinely need it — not by default.

### 5.6 Resolver (observed → canonical)
Projects the log into the graph: node canonicalization (group claims by `primary_id`), edge
construction, dedup + conflict recording, confidence aggregation (v1: a fact is canonical above a
confidence threshold). **Pure and versioned:** `canonical = f(log, resolver_version)`. Bump the
version → recompute with full history intact; drop the store → replay to rebuild.

### 5.7 Query engine (corroboration + timing)
Evaluated **per query, not precomputed**, so the comparison tier and window stay tunable. Per
candidate node: roll up claims to the comparison tier *along typed identity edges only* →
**corroboration** (≥2 *independent* evidence types within the window) → **recency** (earliest
`valid_time`) → **stage/crowding** (from the claim mix). Output: ranked targets + corroborating
claim set + earliest-signal date + full provenance. *Open design question: "independent" should
probably mean **actor**-independent (3 companies) not merely evidence-type-independent (a trial +
its own resulting approval is one program) — a validity decision, surfaced not buried.*

---

## 6. Why these choices (trade-offs)

- **Event sourcing over a mutable store.** Sources revise and get superseded; resolvers improve.
  An append-only log + recomputable projection gives provenance, point-in-time reconstruction,
  and "bad source discovered later → down-weight and recompute" for free. Cost: more storage and
  an eventually-consistent read model. Worth it — auditability and reprocessing are core, not
  optional, for a system whose outputs inform investment decisions.
- **Resolve, never join on a name.** Names are ambiguous and drift; without an identity hub,
  `HER2` and `ERBB2` never corroborate and `MET`/"met" silently corrupts. Cost: a resolver + gold
  set, and resolution errors propagate — mitigated by explicit `ambiguous`/`unresolved` outcomes
  and the unresolved-rate alarm.
- **Bitemporality recorded now, reconstruction deferred.** The two timestamps are nearly free to
  store and impossible to backfill; the point-in-time *query* logic can come later.
- **Per-query corroboration.** Precomputing freezes the tier/window into the data. Evaluating per
  query keeps them tunable; materialize later if query cost bites.
- **Dual drug→target source.** Phase 0 showed a single source (ChEMBL) is precise but blind to
  exactly the high-value assets (emerging, code-named). GtoPdb's independent curation closes much
  of that gap; both feed one drug-identity resolver.
- **Python.** I/O-bound (APIs), best-in-class biomedical/data ecosystem (HGVS→VRS/SPDI tooling,
  pandas/DuckDB), typed models via pydantic/dataclasses. Throughput is gated by source rate
  limits, not CPU, so the language's speed is irrelevant.

---

## 7. Testing & validation strategy

- **Gold sets** for resolution (the gate in 5.2). The hardest, most valuable artifact.
- **Fixtures:** saved raw payloads → deterministic, network-free extraction tests that double as
  schema-drift regression.
- **Invariant tests:** *idempotency* (run twice → identical state) and **recomputability** (drop
  canonical, replay log, assert byte-identical) — the single highest-leverage test in the system.
- **Independent extraction validation:** grade Path B against a *non-circular* source (GtoPdb —
  curated separately from ChEMBL) plus an LLM-judge over the record's own primary text; never
  against ChEMBL-derived sources (Open Targets), which would be self-confirming.
- **Observability as a live test:** the unresolved-rate alarm catches source drift in production.

---

## 8. Failure modes → handling
| Failure | Handling |
|---|---|
| Name mismatch (ERBB2 vs HER2) | resolved to one `HGNC:3430` via the alias index |
| Drug-name miss (salt/code-name/combo) | dual ChEMBL+GtoPdb; salt/parent walk; flagged `unresolved` if neither |
| Source schema drift | unresolved-rate alarm; only that connector affected |
| Ambiguous alias | flagged `ambiguous`, never auto-picked |
| Allele vs class conflation | typed nodes + hard boundaries; no cross-tier roll-up |
| Family/complex explosion | collapse to the specifying target; ADC → antibody not payload |
| Transient API error | retry/backoff; **never cached as a terminal negative** |
| Deterministic empty (e.g. no curated target) | recorded as a clean negative, distinct from error |
| Duplicate ingestion | `content_hash` / `claim_key` no-op |
| Source outage | last good cursor retained; resume next run |
| Bad source discovered later | down-weight + recompute canonical from the log |

---

## 9. Build sequence (gated)

Each phase ships behind an acceptance gate; never build N+1 on an unvalidated N.

0. **Feasibility probe** — *done*; GO decision (extraction is accurate, recall-limited).
1. **Identity spine** — resolver + index + gold-set gate. *Gate:* passes the gold set; `HER2`
   and `ERBB2` resolve to the same node.
2. **Walking skeleton, one source (CT.gov)** — connector → raw store → extract → resolve →
   append-only log → trivial canonical → list query. *Gate:* one record flows end-to-end with
   provenance; idempotent re-run adds nothing.
3. **Second source + corroboration** — add openFDA; the corroboration/timing query. *Gate:* the
   headline trace ("HER2 corroborated by trial ∧ approval, earliest-signal date") runs as a test.
4. **Harden** — idempotency, cursors, bitemporal stamping, alarms, versioning. *Gate:*
   drop-and-replay is byte-identical; injected drift trips the alarm.
5. **Patents** — the earliest competitive signal (leads trials ~2 yrs) and the code-named-asset
   fix. Heavy ingestion; OCSR/Markush modelled but deferred.
6. **Depth & handoff** — richer node types, comparison-tier engine, graph DB if needed; hand
   ranked targets to Layers 2–6 (validation, FTO, covalent feasibility) with full provenance.

---

## 10. Risks & open questions (honest)

- **Entity-resolution precision at scale** — the central technical risk. Gold sets gate it, but
  coverage grows faster than the gold set; the unresolved/ambiguous rates are the early warning.
- **Drug-name recall** — the main yield lever; dual-source helps, but brand-new code-named assets
  live in *no* reference DB and require patents/pipeline disclosures (Phase 5).
- **Family/complex collapse** — "when is a complex a target?" is a genuine modeling decision, not
  a bug.
- **Biomarker vs mechanism** — Path A and Path B mean different things by "target"; the product
  must choose, or carry both typed (current plan).
- **Corroboration independence** — type- vs actor-independence materially changes signal validity.
- **Scale thresholds** — when per-query corroboration or the relational graph needs materialization
  or a graph DB. Deliberately deferred until measured.

---

## 11. Stack (concrete)
Python 3.11+ · httpx (connectors) · pydantic/dataclasses (typed models) · DuckDB + Parquet
(claim log + canonical, columnar replay) · SQLite/Postgres at scale · HGVS→VRS/SPDI tooling for
allele claims · pytest (fixtures + invariant tests). No framework, no orchestration engine until
the DAG warrants one (then Dagster/Prefect). Deploy as scheduled batch jobs; the system is
embarrassingly resumable by construction.
```
targetradar/
  core/        # EntityRef, ObservedClaim, versioning, types
  sources/     # one connector per evidence source
  identity/    # reference ingest, index, resolver, gold sets
  extract/     # field maps, path A/B, claim construction
  store/        # raw, claim log, identity index, canonical
  resolve/     # observed -> canonical projection (pure, versioned)
  query/       # corroboration + timing
```
Boundaries mirror the data flow; `resolve/` and `query/` are pure functions over `store/`, which
is what makes the whole thing testable and recomputable.
```
