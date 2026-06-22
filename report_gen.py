"""Generate the engine's persuasion report (clean Markdown) for the top targets. Three pillars per
target — WHO'S MOVING (patents + programs), WHY NOW (momentum + biology), SOURCES/IDs (every claim
linked). Assembles from existing artifacts; deterministic except the program/biology text, which comes
from the verification analyses already on disk.
"""
import datetime
import json
import re
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).parent


def _load(name, default=None):
    p = HERE / name
    return json.load(open(p)) if p.exists() else default


LEGAL = {"INC", "LLC", "LTD", "LIMITED", "CO", "CORP", "CORPORATION", "COMPANY", "GMBH", "AG", "SA",
         "SAS", "BV", "NV", "KG", "KGAA", "PLC", "LP", "LLP", "PTE", "PTY", "SRL", "SPA", "OY", "AB", "AS", "ULC"}
# degraders/glues: warhead = linker/E3 chemistry, not covalent target engagement -> not a covalent-inhibitor race
DEGRADER = re.compile(r"DEGRADER|PROTAC|MOLECULAR GLUE|BIFUNCTIONAL|CEREBLON|\bCRBN\b|VHL LIGAND|GLUTARIMIDE|"
                       r"\bIMID\b|E3 (UBIQUITIN )?LIGASE|ISOINDOLIN|PIPERIDINE-2,6-DIONE")


def _gpatent(pn):                                  # SureChEMBL id -> Google Patents url
    return f"https://patents.google.com/patent/{pn.replace('-', '')}"


def _firm(a):                                      # normalize assignee: one company = one player
    a = re.sub(r"\(.*?\)", " ", (a or "").upper())
    a = re.sub(r"[.,/]", " ", a)
    toks = a.split()
    while toks and toks[-1] in LEGAL:
        toks.pop()
    return " ".join(toks).strip()


def run():
    ranking = _load("covalent_ranking_final.json", [])
    deep = {d["sym"]: d for d in (_load("leveled_ranking.json", {}) or {}).get("deep", [])}
    fc = _load("factcheck.json", {})
    scan = {r["sym"]: r for r in (_load("emerging_open_ranked.json", {}) or {}).get("all", [])}
    hot = {r["sym"]: r for r in (_load("emerging_open_ranked.json", {}) or {}).get("hot_ranked", [])}
    tg = _load("patent_targets.json", {})
    feed = {e[0]: e for e in (_load("covalent_warhead_feed_clean.json", []) or [])}

    # patents (title/claims) per target: company + id + date
    pats = defaultdict(list)
    for pn, d in tg.items():
        if d.get("sym") and d.get("conf") in ("title", "claims") and pn in feed:
            e = feed[pn]
            pats[d["sym"]].append({"pn": pn, "assignee": e[2].split(';')[0].strip(), "date": e[1], "title": e[3]})

    today = datetime.date.today().isoformat()
    out = [f"# Covalent target radar — {today}",
           f"\n_{len(ranking)} top targets. Each: **who's moving** · **why now** · **sources**. "
           "Selectivity is fact-checked against UniProt; patents from SureChEMBL; momentum from Europe PMC._\n",
           "## At a glance\n",
           "| # | Target | Selectivity | Lane | Patents | Momentum |",
           "|---|--------|-------------|------|---------|----------|"]

    for r in ranking:
        sym = r["target"]
        dd = deep.get(sym, {})
        npat = len(pats.get(sym, []))
        accel = scan.get(sym, {}).get("accel")
        mom = f"×{accel}" if accel and accel >= 1.15 else "—"
        out.append(f"| {r['rank']} | **{sym}** | {r.get('selectivity_factchecked','?')} | "
                   f"{dd.get('lane','?')} | {npat or '—'} | {mom} |")

    out.append("\n---\n")

    for r in ranking:
        sym = r["target"]
        dd = deep.get(sym, {})
        acc = fc.get(sym, {}).get("uniprot")
        out.append(f"\n## #{r['rank']} · {sym}  —  `{dd.get('verdict', r.get('verdict','?'))}` · lane: **{dd.get('lane','?')}**")
        out.append(f"> {dd.get('covalent_handle', r.get('note',''))[:200]}")

        # ── WHO'S MOVING ──────────────────────────────────────────────
        out.append("\n**🏢 Who's moving**")
        ps = sorted(pats.get(sym, []), key=lambda x: x["date"], reverse=True)
        if ps:
            for p in ps[:5]:
                out.append(f"- Patent **{p['assignee']}** — [{p['pn']}]({_gpatent(p['pn'])}) ({p['date']}) · _{p['title'][:60]}_")
        else:
            out.append("- No covalent patent in the SureChEMBL feed → **open IP lane (first-mover)**")
        ec = dd.get("existing_covalent")
        if ec:
            out.append(f"- Programs: {ec}")

        # ── WHY NOW ───────────────────────────────────────────────────
        out.append("\n**📈 Why now**")
        sc = scan.get(sym, {}); h = hot.get(sym, {})
        if sc.get("accel"):
            line = f"- Literature **×{sc['accel']}** ({sc.get('recent')} papers '25-26 vs {sc.get('base')} '21-23)"
            if h.get("preprints_25_26"):
                line += f" · {h['preprints_25_26']} preprints"
            if h.get("new_trials_24_26"):
                line += f" · {h['new_trials_24_26']} new trials"
            out.append(line)
        out.append(f"- Biology / validation: {dd.get('validation', '—')}")

        # ── SOURCES ───────────────────────────────────────────────────
        pdbs = sorted(set(re.findall(r'PDB[:\s]+([0-9][A-Za-z0-9]{3})\b', dd.get("covalent_handle", "") + " " + dd.get("validation", ""))))
        src = []
        if acc:
            src.append(f"[UniProt {acc}](https://www.uniprot.org/uniprotkb/{acc})")
        for pid in pdbs[:3]:
            src.append(f"[PDB {pid.upper()}](https://www.rcsb.org/structure/{pid.upper()})")
        src.append(f"[trials](https://clinicaltrials.gov/search?term={sym})")
        src.append(f"[literature](https://europepmc.org/search?query={sym})")
        if ps:
            src.append("SureChEMBL patents (linked above)")
        src.append("CysDB/KLIFS handle · fact-check vs UniProt")
        out.append("\n**🔗 Sources** · " + " · ".join(src))
        out.append("")

    # ── FRESH CREDIBLE FILINGS (fast-follow tier) ────────────────────────────────
    # credible company JUST filed covalent IP (2025+), few players, lane not yet crowded.
    top_syms = {r["target"] for r in ranking}
    prog = defaultdict(lambda: {"firms": set(), "pats": []})
    for pn, d in tg.items():
        if d.get("sym") and d.get("conf") in ("title", "claims") and pn in feed:
            e = feed[pn]
            title = (e[3] or "").upper()
            if (e[1] or "")[:4] < "2024" or DEGRADER.search(title):
                continue                                        # skip pre-2024 + degrader/glue patents
            common = (d.get("common") or "").upper()
            if d["sym"].upper() not in title and (not common or common not in title):
                continue                                        # title must NAME the target = "filed on X", not just mentions it
            firm = _firm(e[2].split(';')[0])
            if not firm:
                continue
            prog[d["sym"]]["firms"].add(firm)
            prog[d["sym"]]["pats"].append({"pn": pn, "firm": firm, "date": e[1], "title": e[3]})
    fresh = []
    for sym, p in prog.items():
        if sym in top_syms:
            continue
        players, latest = len(p["firms"]), max(x["date"] for x in p["pats"])
        if 1 <= players <= 3 and latest >= "2025":             # credible (feed-filtered) · uncrowded · fresh
            fresh.append((sym, players, latest, p))
    fresh.sort(key=lambda x: (x[1], x[2]), reverse=False)       # fewest players first
    fresh.sort(key=lambda x: x[2], reverse=True)                # then most recent

    out += ["\n---\n", "## ⚡ Fresh credible filings — fast-follow tier",
            "\n_A credible company just filed covalent IP (2025+), only 1–3 players, lane not yet crowded "
            "— someone moved; the race is open._\n"]
    for sym, players, latest, p in fresh[:12]:
        ps = sorted(p["pats"], key=lambda x: x["date"], reverse=True)
        accel = scan.get(sym, {}).get("accel")
        out.append(f"\n### ⚡ {sym}  —  {players} player(s) · latest filing **{latest}**")
        out.append("\n**🏢 Who's moving**")
        for x in ps[:4]:
            out.append(f"- **{x['firm']}** — [{x['pn']}]({_gpatent(x['pn'])}) ({x['date']}) · _{x['title'][:62]}_")
        out.append(f"\n**📈 Why now** · filed **{latest}** (fresh) · only **{players}** player(s) — race just starting"
                   + (f" · literature ×{accel}" if accel and accel >= 1.15 else ""))
        out.append(f"**🔗 Sources** · [trials](https://clinicaltrials.gov/search?term={sym}) · "
                   f"[literature](https://europepmc.org/search?query={sym}) · SureChEMBL patents (linked above)")

    md = "\n".join(out)
    (HERE / f"covalent_report_{today}.md").write_text(md)
    print(md[:2000])
    print(f"\n... saved covalent_report_{today}.md ({len(md)} chars, {len(ranking)} targets)")


if __name__ == "__main__":
    run()
