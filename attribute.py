"""Attribute each covalent feed-patent to its REAL target gene using SureChEMBL's text-mined
GeneOrProtein annotations (HGNC-resolved), validated + relabeled against the authoritative HGNC set.

Chain: feed.patent_number -> patents.id -> biomedical_locations.patent_id -> biomedical_entities
(entity_type_id=1, resolved_form LIKE 'HGNC:%').

Two denoising steps over the raw annotations:
  1. VALIDATE: keep a mention only if its surface text is a recognized HGNC form (approved/alias/prev
     symbol or name) of the resolved gene -> kills NER false-positives like the French word "leurs"
     mis-resolved to HGNC:17095, "dddd", etc.
  2. RANK title-first: primary target = gene named in the TITLE (a patent titled "X inhibitors" IS X),
     tie-break claims -> abstract -> body. Beats per-mention weighting (a pathway gene cited 50x in the
     description used to outscore the real target named once in the title).
Labels use the official HGNC symbol.  Writes patent_targets.json: {pn: {hgnc, sym, conf, others:[...]}}.
"""
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import duckdb

HERE = Path(__file__).parent
B = "https://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/bulk_data/2026-06-15"
ENT, LOC, HGNC = "/tmp/biomedical_entities.parquet", "/tmp/biomedical_locations.parquet", "/tmp/hgnc.txt"


def run():
    feed = json.load(open(HERE / "covalent_warhead_feed_clean.json"))
    pns = sorted({e[0] for e in feed})
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs;")

    # HGNC reference: official symbol per id, and every recognized surface form (symbol/alias/prev/name)
    con.sql(f"""CREATE TEMP TABLE hgnc_sym AS
        SELECT hgnc_id AS hgnc, symbol FROM read_csv('{HGNC}', delim='\t', header=true, quote='', all_varchar=true)""")
    # collision-prone forms to drop: <3 chars (DF/DR/"1") + EN/FR/DE function words that are also HGNC
    # aliases (the French "leurs" maps to LARS2's alias "LEURS"; "set"/"son"/"was"/"est"/"des"...).
    STOP = ("THE AND FOR USE USES USED NEW ALL ONE TWO CAN MAY NOT ARE WAS WERE OUR HER HIS ITS WHO HOW WHY "
            "OFF OUT VIA PER ANY HAS HAD OWN SET SON LEURS DES AUX SUR PAR POUR UNE SES EST PAS QUE QUI AVEC "
            "SANS DANS LES CES NOS VOS PLUS FAIT DER DIE DAS UND FUR MIT EIN VON AUF IST SIE DEN DEM EINE").split()
    stopsql = ",".join(f"'{w}'" for w in STOP)
    con.sql(f"""CREATE TEMP TABLE hgnc_form AS
        SELECT DISTINCT hgnc, upper(trim(form)) AS form FROM (
          SELECT hgnc_id AS hgnc, unnest(list_concat([symbol, name],
                 string_split(coalesce(alias_symbol,''), '|'), string_split(coalesce(alias_name,''), '|'),
                 string_split(coalesce(prev_symbol,''), '|'), string_split(coalesce(prev_name,''), '|'))) AS form
          FROM read_csv('{HGNC}', delim='\t', header=true, quote='', all_varchar=true))
        WHERE length(trim(form)) >= 3 AND upper(trim(form)) NOT IN ({stopsql})""")

    # feed patent_numbers -> internal patents.id (remote projection, filtered to the feed only)
    con.sql("CREATE TEMP TABLE feed(pn VARCHAR)")
    con.executemany("INSERT INTO feed VALUES (?)", [(p,) for p in pns])
    con.sql(f"""CREATE TEMP TABLE fp AS
        SELECT p.id, p.patent_number AS pn
        FROM read_parquet('{B}/patents.parquet') p JOIN feed f ON f.pn = p.patent_number""")
    print("feed patents matched to ids:", con.sql("SELECT count(*) FROM fp").fetchone()[0], "/", len(pns))

    # VALIDATED gene/protein entities: surface text is a recognized HGNC form of the resolved gene
    con.sql(f"""CREATE TEMP TABLE ent AS
        SELECT e.id, e.resolved_form AS hgnc, e.corrected_text AS sfc
        FROM read_parquet('{ENT}') e
        WHERE e.entity_type_id = 1 AND e.resolved_form LIKE 'HGNC:%'
          AND EXISTS (SELECT 1 FROM hgnc_form h WHERE h.hgnc = e.resolved_form
                      AND h.form IN (upper(e.corrected_text), upper(e.original_text)))""")
    print("validated gene entities:", con.sql("SELECT count(*) FROM ent").fetchone()[0])

    # all validated gene mentions in feed patents (one local scan of the 2.2GB link table) -> tiny table
    con.sql(f"""CREATE TEMP TABLE gm AS
        SELECT fp.pn, ent.hgnc, ent.sfc, l.field_id AS f, l.count AS c
        FROM read_parquet('{LOC}') l JOIN fp ON fp.id = l.patent_id JOIN ent ON ent.id = l.entity_id""")
    print("gene-mention rows:", con.sql("SELECT count(*) FROM gm").fetchone()[0])

    def clean(s):                                    # gene symbols are uppercase; leave full-name phrases alone
        return s.upper() if (s and " " not in s and len(s) <= 8) else s
    # common name per gene = its dominant surface form (what trials/people actually call it; better CT.gov term)
    common = {h: clean(s) for h, s in con.sql("""
        SELECT hgnc, sfc FROM (SELECT hgnc, sfc, row_number() OVER (PARTITION BY hgnc ORDER BY sum(c) DESC) rn
        FROM gm GROUP BY hgnc, sfc) WHERE rn = 1""").fetchall()}

    # primary target per patent: title-first (1=desc 2=clms 3=abst 4=ttl), labeled by official HGNC symbol
    rows = con.sql("""
        SELECT g.pn, g.hgnc, s.symbol, g.tw, g.cw
        FROM (SELECT pn, hgnc, sum(c) FILTER (f=4) tw, sum(c) FILTER (f=2) cw,
                     sum(c) FILTER (f=3) aw, sum(c) allc
              FROM gm GROUP BY pn, hgnc) g
        JOIN hgnc_sym s ON s.hgnc = g.hgnc
        ORDER BY g.pn, g.tw DESC, g.cw DESC, g.aw DESC, g.allc DESC""").fetchall()

    # title-substring QC: a TITLE-ranked gene must actually appear as a token in the (English) title,
    # else it is a phantom tag (multilingual boilerplate / alias collision, e.g. F2R from "IRAK degraders")
    # -> skip it so the real claims/body target wins. Claims/body genes are trusted (no claims text on hand).
    titles = {e[0]: (e[3] or "").upper() for e in json.load(open(HERE / "covalent_warhead_feed_clean.json"))}
    forms = defaultdict(set)
    for h, fm in con.sql("SELECT hgnc, form FROM hgnc_form WHERE hgnc IN (SELECT DISTINCT hgnc FROM gm)").fetchall():
        forms[h].add(fm)

    def in_title(hgnc, title):
        return any(re.search(r"(?<![A-Z0-9])" + re.escape(f) + r"(?![A-Z0-9])", title) for f in forms.get(hgnc, ()))

    ranked = defaultdict(list)
    for pn, hgnc, sym, tw, cw in rows:
        ranked[pn].append((hgnc, sym, tw, cw))
    by_pn, skipped = {}, 0
    for pn, genes in ranked.items():
        title = titles.get(pn, "")
        primary, others = None, []
        for hgnc, sym, tw, cw in genes:
            if tw and not in_title(hgnc, title):
                skipped += 1
                continue                                   # phantom title tag
            if primary is None:
                primary = (hgnc, sym, "title" if tw else "claims" if cw else "body")
            elif len(others) < 4:
                others.append(sym)
        if primary:
            h, sym, conf = primary
            by_pn[pn] = {"sym": sym, "common": common.get(h, sym), "conf": conf, "others": others}
    print(f"title-QC skipped {skipped} phantom title-gene tags")

    json.dump(by_pn, open(HERE / "patent_targets.json", "w"), indent=1)
    print(f"attributed {len(by_pn)}/{len(pns)} feed patents to a validated HGNC target "
          f"({len(pns)-len(by_pn)} have no usable gene annotation)")
    print("confidence:", dict(Counter(d["conf"] for d in by_pn.values())))
    print("top targets:", Counter(d["sym"] for d in by_pn.values()).most_common(18))


if __name__ == "__main__":
    run()
