"""Minimal EPO OPS client — OAuth2 + biblio search + full-text. Verified against the live API.

Creds live in .epo_creds.json (gitignored): {"key": ..., "secret": ...}; used only to authenticate.
OPS full-text coverage is PARTIAL (many recent EP applications 404) — callers treat None as "no full-text".
"""
import base64
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

_BASE = "https://ops.epo.org/3.2/rest-services"
_AUTH = "https://ops.epo.org/3.2/auth/accesstoken"
_CREDS = Path(__file__).parent / ".epo_creds.json"
_tok = None


def _token(force=False):
    global _tok
    if _tok and not force:
        return _tok
    c = json.load(open(_CREDS))
    a = base64.b64encode(f"{c['key']}:{c['secret']}".encode()).decode()
    req = urllib.request.Request(_AUTH, data=b"grant_type=client_credentials",
                                headers={"Authorization": f"Basic {a}",
                                         "Content-Type": "application/x-www-form-urlencoded"})
    _tok = json.loads(urllib.request.urlopen(req, timeout=30).read())["access_token"]
    return _tok


def _get(path, accept="application/xml"):
    """GET an OPS path -> (http_status, text). Re-auths once on 401/403 (token ~20 min TTL)."""
    for attempt in (1, 2):
        req = urllib.request.Request(_BASE + path,
                                     headers={"Authorization": f"Bearer {_token(force=attempt == 2)}",
                                              "Accept": accept})
        try:
            r = urllib.request.urlopen(req, timeout=45)
            return r.getcode(), r.read().decode("utf-8", "ignore")
        except urllib.error.HTTPError as e:
            if e.code in (401, 403) and attempt == 1:
                continue                                    # token expired -> force-refresh and retry once
            return e.code, ""
    return None, ""


def search(cql, n=25):
    """OPS biblio search -> list of epodoc publication numbers (e.g. 'EP4722284'), de-duped, in order."""
    _, xml = _get(f"/published-data/search?q={urllib.parse.quote(cql)}&Range=1-{n}")
    refs = re.findall(r'<country>(\w+)</country>\s*<doc-number>(\d+)</doc-number>\s*<kind>\w+</kind>', xml)
    return list(dict.fromkeys(f"{c}{num}" for c, num in refs))


def fulltext(pn, part="claims"):
    """part: 'claims' | 'description'. Returns XML text, or None when OPS has no full-text (404/empty)."""
    code, body = _get(f"/published-data/publication/epodoc/{pn}/{part}")
    return body if code == 200 and len(body) > 200 else None


def _get_bytes(path):
    """GET a binary OPS endpoint (image/PDF) -> (status, bytes). Re-auths once on 401/403."""
    for attempt in (1, 2):
        req = urllib.request.Request(_BASE + path,
                                     headers={"Authorization": f"Bearer {_token(force=attempt == 2)}"})
        try:
            r = urllib.request.urlopen(req, timeout=60)
            return r.getcode(), r.read()
        except urllib.error.HTTPError as e:
            if e.code in (401, 403) and attempt == 1:
                continue
            return e.code, b""
    return None, b""


def images_inventory(pn):
    """Image inventory -> {'link', 'pages', 'sections': {NAME: start_page}} or None. (Images have wider
    coverage than full-text — present even for many patents that 404 on /description.)"""
    code, xml = _get(f"/published-data/publication/epodoc/{pn}/images")
    m = re.search(r'<ops:document-instance[^>]*desc="FullDocument"[^>]*>', xml) if code == 200 else None
    if not m:
        return None
    link = re.search(r'link="([^"]+)"', m.group(0))
    pages = re.search(r'number-of-pages="(\d+)"', m.group(0))
    return {"link": link.group(1) if link else None,
            "pages": int(pages.group(1)) if pages else None,
            "sections": {s: int(p) for s, p in re.findall(r'name="(\w+)"\s+start-page="(\d+)"', xml)}}


def page_pdf(link, page):
    """One page of a patent as PDF bytes (OPS serves exactly one page per call). `link` from images_inventory."""
    code, data = _get_bytes(f"/{link}.pdf?Range={page}")
    return data if code == 200 and data[:4] == b"%PDF" else None


def biblio(pn):
    """-> (publication_date 'YYYY-MM-DD', applicant, english title) from OPS biblio, best-effort (None on miss)."""
    code, xml = _get(f"/published-data/publication/epodoc/{pn}/biblio")
    if code != 200:
        return (None, None, None)
    pr = re.search(r"<publication-reference\b.*?</publication-reference>", xml, re.S)
    dm = re.search(r"<date>\s*(\d{8})\s*</date>", pr.group(0)) if pr else None
    appl = re.search(r"<applicant-name>\s*<name[^>]*>([^<]+)</name>", xml)
    ttl = (re.search(r'<invention-title[^>]*lang="en"[^>]*>([^<]+)</invention-title>', xml)
           or re.search(r"<invention-title[^>]*>([^<]+)</invention-title>", xml))
    d = dm.group(1) if dm else None
    return (f"{d[:4]}-{d[4:6]}-{d[6:8]}" if d else None,
            appl.group(1).strip() if appl else None,
            ttl.group(1).strip() if ttl else None)


if __name__ == "__main__":
    # live smoke: auth + search + full-text availability (exercises the re-auth retry + the search regex)
    pubs = search('txt="acrylamide" and ic=C07D and pn=EP', n=10)
    print(f"search -> {len(pubs)} pubs:", pubs[:5])
    assert pubs, "search returned nothing — check creds / OPS quota"
    cov = [(p, fulltext(p, "claims") is not None) for p in pubs[:5]]
    print("full-text claims available:", cov)
    assert any(has for _, has in cov), "no full-text among first 5 (coverage gap or auth issue)"
    print("OK")
