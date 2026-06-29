"""Refresh the static covalent-patent census to the latest SureChEMBL snapshot.

The census (covalent_warhead_feed_clean.json -> attribute.py -> patent_targets.json) is frozen to ONE
bi-weekly SureChEMBL snapshot. This DISCOVERS the latest *complete* snapshot (SureChEMBL sometimes lists a
date whose parquets aren't uploaded yet — e.g. 2026-06-16), (incrementally) rebuilds the covalent-warhead
feed, re-attributes + grades it, and writes snapshot.txt so attribute.py picks it up. Cron it bi-weekly.

  python refresh.py           # incremental: SMARTS only over patents newer than the current feed, then re-attribute
  python refresh.py --full    # full rebuild (SMARTS over ALL compounds — slow, scans the multi-GB compounds parquet)

Heavy steps (a ~GB biomedical_locations download + RDKit over the new compounds) — built for an offline/cron
run, not interactive. The warhead definitions are reused from patent_radar, so the feed and the live radar
agree on what 'covalent' means.
"""
import datetime
import json
import re
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import duckdb
from rdkit import Chem

from patent_radar import _SMARTS                # the same compiled warhead SMARTS the live radar uses

HERE = Path(__file__).parent
BULK = "https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data"
FEED = HERE / "covalent_warhead_feed_clean.json"
MIN_CPDS = 3                                    # a patent is "covalent" if >= this many of its compounds carry a warhead


def latest_snapshot(today=None):
    """Newest snapshot dir (date <= today) whose parquets are actually uploaded (DESCRIBE succeeds)."""
    today = today or datetime.date.today().isoformat()
    html = urllib.request.urlopen(BULK + "/", timeout=60).read().decode("utf-8", "ignore")
    con = duckdb.connect(); con.sql("INSTALL httpfs; LOAD httpfs;")
    for d in sorted({m for m in re.findall(r"20\d\d-\d\d-\d\d", html) if m <= today}, reverse=True):
        try:
            con.sql(f"DESCRIBE SELECT * FROM read_parquet('{BULK}/{d}/patents.parquet') LIMIT 1")
            return d
        except Exception:
            continue                            # listed but not uploaded yet -> skip to the previous one
    raise RuntimeError("no complete SureChEMBL snapshot found")


def download_bio(snap):
    """attribute.py reads the gene-NER parquets locally — fetch them for this snapshot (locations is ~GB)."""
    for f in ("biomedical_entities", "biomedical_locations"):
        print(f"  downloading {f}.parquet ...", flush=True)
        urllib.request.urlretrieve(f"{BULK}/{snap}/{f}.parquet", f"/tmp/{f}.parquet")


def _covalent_patents(rows):
    """[(pn, date, assignee, title, smiles), ...] -> feed rows for patents with >= MIN_CPDS covalent compounds.
    Split out from the remote query so the warhead logic is testable offline (the query is the slow part)."""
    per = defaultdict(lambda: {"meta": None, "wh": set(), "n": 0})
    for pn, dt, asg, title, smi in rows:
        m = Chem.MolFromSmiles(smi) if smi else None
        hits = [k for k, pat in _SMARTS.items() if m and m.HasSubstructMatch(pat)] if m else []
        if hits:
            e = per[pn]; e["meta"] = (dt, asg, title); e["wh"].update(hits); e["n"] += 1
    return [[pn, str(e["meta"][0]), e["meta"][1] or "", e["meta"][2] or "", e["n"], sorted(e["wh"])]
            for pn, e in per.items() if e["n"] >= MIN_CPDS]


def build_feed(snap, since=None, country=None, limit=None):
    """SMARTS over the snapshot's compounds -> [pn, date, assignee, title, n_covalent_cpds, warheads] for
    patents with >= MIN_CPDS covalent compounds. `since` (publication_date) makes it INCREMENTAL. The join
    streams the multi-GB compounds parquet, so this is a heavy offline/cron step, not interactive."""
    con = duckdb.connect(); con.sql("INSTALL httpfs; LOAD httpfs;")
    B = f"{BULK}/{snap}"
    conds = [c for c in (f"p.publication_date > '{since}'" if since else None,
                         f"p.country = '{country}'" if country else None) if c]
    where = ("WHERE " + " AND ".join(conds)) if conds else ""
    rows = con.sql(f"""
        SELECT p.patent_number, p.publication_date, p.assignee, p.title, c.smiles
        FROM read_parquet('{B}/patents.parquet') p
        JOIN read_parquet('{B}/patent_compound_map.parquet') m ON m.patent_id = p.id
        JOIN read_parquet('{B}/compounds.parquet') c ON c.id = m.compound_id
        {where} {f'LIMIT {limit}' if limit else ''}""").fetchall()
    return _covalent_patents(rows)


def run(full=False):
    snap = latest_snapshot()
    (HERE / "snapshot.txt").write_text(snap + "\n")
    print(f"latest complete SureChEMBL snapshot: {snap}  (-> snapshot.txt; attribute.py now uses it)")

    feed = json.loads(FEED.read_text()) if FEED.exists() else []
    have = {e[0] for e in feed}
    since = None if full else max((e[1] for e in feed), default=None)
    print(f"{'FULL rebuild (all compounds)' if full else f'incremental: patents after {since}'} ...")

    download_bio(snap)
    found = build_feed(snap, since=since)
    added = [e for e in found if e[0] not in have]
    feed = found if full else feed + added
    FEED.write_text(json.dumps(feed))
    print(f"feed: +{len(added)} new covalent patents -> {len(feed)} total")

    import attribute                              # re-attribute + grade the whole feed against the new snapshot
    attribute.run()
    print("census refreshed: patent_targets.json regraded; consumers read it via trusted().")


if __name__ == "__main__":
    if "--selftest" in sys.argv:                # offline check of the feed-build core (no network)
        demo = [("WO-1", "2026-06-12", "ACME", "BTK inhibitors", "C=CC(=O)Nc1ccccc1"),     # acrylamide
                ("WO-1", "2026-06-12", "ACME", "BTK inhibitors", "C=CC(=O)NCc1ccncc1"),     # acrylamide
                ("WO-1", "2026-06-12", "ACME", "BTK inhibitors", "ClCC(=O)Nc1ccccc1"),      # chloroacetamide
                ("WO-2", "2026-06-12", "X", "aspirin", "CC(=O)Oc1ccccc1C(=O)O")]            # not covalent
        out = _covalent_patents(demo)
        assert [e[0] for e in out] == ["WO-1"] and out[0][4] == 3, out   # WO-1: 3 covalent cpds; WO-2 dropped
        assert set(out[0][5]) == {"acrylamide", "chloroacetamide"}, out
        print("refresh self-check OK:", out)
    else:
        run(full="--full" in sys.argv)
