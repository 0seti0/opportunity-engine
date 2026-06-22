"""Gene-identity resolver: a name mention -> canonical HGNC id, built from the HGNC complete set."""
import csv
from collections import namedtuple
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional

HGNC_TSV = "/tmp/hgnc.tsv"  # ponytail: reads the file we already downloaded; add download-if-missing when run unattended

Index = namedtuple("Index", "symbols aliases")  # symbols: name->id (unique); aliases: name->set of ids


class Outcome(Enum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    UNRESOLVED = "unresolved"


@dataclass
class Resolution:
    outcome: Outcome
    hgnc_id: Optional[str] = None
    candidates: list = field(default_factory=list)


def build_index(rows):
    """Official symbols (unique) kept apart from aliases/prev-symbols (which collide). All uppercased."""
    symbols, aliases = {}, {}
    for r in rows:
        gid = r["hgnc_id"]
        symbols[r["symbol"].upper()] = gid
        for name in (r.get("alias_symbol", "") + "|" + r.get("prev_symbol", "")).split("|"):
            if name:
                aliases.setdefault(name.upper(), set()).add(gid)
    return Index(symbols, aliases)


def resolve(mention, index):
    key = mention.strip().upper()
    if key in index.symbols:                                       # official symbol wins (unique -> never ambiguous)
        return Resolution(Outcome.RESOLVED, hgnc_id=index.symbols[key])
    ids = index.aliases.get(key)
    if not ids:
        return Resolution(Outcome.UNRESOLVED)
    if len(ids) == 1:
        return Resolution(Outcome.RESOLVED, hgnc_id=next(iter(ids)))
    return Resolution(Outcome.AMBIGUOUS, candidates=sorted(ids))


def load_rows(path=HGNC_TSV):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


GOLD = [  # (mention, expected outcome, expected id) -- the acceptance gate; curate freely
    ("ERBB2", Outcome.RESOLVED,   "HGNC:3430"),   # official symbol
    ("HER2",  Outcome.RESOLVED,   "HGNC:3430"),   # alias
    ("her2",  Outcome.RESOLVED,   "HGNC:3430"),   # case-insensitive (decision #2: case-folded)
    ("NGL",   Outcome.RESOLVED,   "HGNC:3430"),   # prev-symbol / lifecycle
    ("CCR10", Outcome.RESOLVED,   "HGNC:4474"),   # symbol beats alias (decision #1 = B)
    ("ASP",   Outcome.AMBIGUOUS,  None),          # true alias-vs-alias ambiguity
    ("xyzzy", Outcome.UNRESOLVED, None),          # negative control
]


if __name__ == "__main__":
    index = build_index(load_rows())
    for mention, want, want_id in GOLD:
        got = resolve(mention, index)
        assert got.outcome is want and got.hgnc_id == want_id, \
            f"{mention}: got {got.outcome.value}/{got.hgnc_id}, want {want.value}/{want_id}"
    print(f"gold set {len(GOLD)}/{len(GOLD)} green | "
          f"{len(index.symbols):,} symbols + {len(index.aliases):,} alias-names indexed")
