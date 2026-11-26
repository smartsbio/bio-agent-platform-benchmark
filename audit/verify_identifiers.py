#!/usr/bin/env python3
"""Resolve every DOI, PMID and PDB id cited in the transcripts.

Rationale: the blinded LLM judge raised four hallucination flags across batches
A and B. All four were false positives, three of them on identifiers that were
simply unfamiliar to it — most notably the bioRxiv/medRxiv DOI prefix 10.64898,
which belongs to openRxiv and which the judge called "suspect". A model cannot
adjudicate whether an accession exists; a registry lookup can.

So identifier existence is now a DETERMINISTIC check and the judge's
hallucination flag is only a pointer to something worth resolving here.

    python3 -I audit/verify_identifiers.py                  # every transcript
    python3 -I audit/verify_identifiers.py --task lit-001
    python3 -I audit/verify_identifiers.py --arm biomni

Writes results/identifier_audit.csv: one row per (run, identifier, verdict).
Responses are cached in audit/.id_cache.json so re-runs are free and the
registries are not hammered. Needs network on first run.

What this does and does not establish: it proves an identifier is registered
and reports its title. It does NOT prove the citation supports the claim it was
attached to — a real paper cited for something it does not say is still wrong,
and that remains a job for human review.
"""

import argparse
import csv
import json
import pathlib
import re
import sys
import time
import urllib.error
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
CACHE = ROOT / "audit" / ".id_cache.json"
OUT = ROOT / "results" / "identifier_audit.csv"

# Trailing sentence punctuation is stripped; DOIs legitimately contain '.' so
# only a trailing run of it is removed, not internal dots.
DOI_RE = re.compile(r'\b10\.\d{4,6}/[^\s)\]\},;"\'<>]+')
PMID_RE = re.compile(r'\bPMID[:\s]*(\d{6,8})\b', re.I)
PDB_RE = re.compile(r'\b(?:PDB[\s:]+(?:ID[\s:]+)?)([0-9][A-Za-z0-9]{3})\b')
# Added after the first run reported "every identifier resolves" across eight
# molecule runs while in fact matching none of them: those answers cite PubChem
# CIDs and ChEMBL ids, which the three patterns above do not cover.
CHEMBL_RE = re.compile(r'\b(CHEMBL\d+)\b', re.I)
# CID 0 is not a compound; it shows up as a placeholder in table cells.
CID_RE = re.compile(r'\b(?:PubChem\s+)?CID[:\s]*([1-9]\d{0,8})\b', re.I)
# A bare-accession pattern also matches molecular formulas — it flagged
# aspirin's C9H8O4 as a fabricated UniProt accession — so an explicit
# "UniProt"/"accession" cue is required. This under-matches by design:
# a missed identifier is a gap, a false fabrication charge is a retraction.
UNIPROT_RE = re.compile(
    r'(?:UniProt(?:KB)?|accession)\s*(?:ID|accession)?[:\s#]*'
    r'([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b', re.I)

UA = {"User-Agent": "smartsbio-platform-benchmark/1.0 (identifier audit)"}


def load_cache():
    if CACHE.exists():
        try:
            return json.loads(CACHE.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def get_json(url, timeout=30):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=timeout) as fh:
        return json.load(fh)


def resolve(kind, ident, cache):
    """-> (verdict, title). Verdict is 'registered', 'absent' or 'error:…'."""
    key = f"{kind}:{ident}"
    if key in cache:
        return tuple(cache[key])
    try:
        if kind == "doi":
            if ident.lower().startswith("10.5281/zenodo"):
                # Zenodo registers with DataCite, not Crossref.
                d = get_json(f"https://api.datacite.org/dois/{ident}")
                titles = d["data"]["attributes"].get("titles") or [{}]
                out = ("registered", (titles[0].get("title") or "")[:90])
            else:
                m = get_json(f"https://api.crossref.org/works/{ident}")["message"]
                out = ("registered", (m.get("title") or [""])[0][:90])
        elif kind == "pmid":
            d = get_json("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
                         f"esummary.fcgi?db=pubmed&id={ident}&retmode=json")
            r = d["result"].get(ident, {})
            if not r or r.get("error"):
                out = ("absent", "")
            else:
                out = ("registered", f"{r.get('title','')} [{r.get('source','')}]"[:90])
        elif kind == "pdb":
            d = get_json(f"https://data.rcsb.org/rest/v1/core/entry/{ident.upper()}")
            out = ("registered",
                   (d.get("struct", {}).get("title") or "")[:90])
        elif kind == "chembl":
            # A CHEMBL id can be a molecule, a target, an assay or a document.
            # Querying only /molecule reported six real targets (COX-1, COX-2,
            # CXCR1/2, 5-LOX) as fabrications, so all four are tried.
            out = ("absent", "")
            for ep in ("molecule", "target", "assay", "document"):
                try:
                    d = get_json("https://www.ebi.ac.uk/chembl/api/data/"
                                 f"{ep}/{ident.upper()}?format=json")
                except urllib.error.HTTPError:
                    continue
                name = (d.get("pref_name") or d.get("description")
                        or d.get("title") or ident)
                out = ("registered", f"{name} [{ep}]"[:90])
                break
        elif kind == "cid":
            d = get_json("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
                         f"cid/{ident}/property/Title,MolecularFormula/JSON")
            p = d["PropertyTable"]["Properties"][0]
            out = ("registered",
                   f"{p.get('Title','')} [{p.get('MolecularFormula','')}]"[:90])
        elif kind == "uniprot":
            d = get_json(f"https://rest.uniprot.org/uniprotkb/{ident.upper()}.json")
            nm = (d.get("proteinDescription", {}).get("recommendedName", {})
                  .get("fullName", {}).get("value") or "")
            out = ("registered", f"{nm} [{d.get('primaryAccession','')}]"[:90])
        else:
            out = ("error:unknown-kind", "")
    except urllib.error.HTTPError as e:
        out = ("absent" if e.code in (404, 400) else f"error:http{e.code}", "")
    except Exception as e:                                  # noqa: BLE001
        out = (f"error:{type(e).__name__}", "")
    cache[key] = list(out)
    time.sleep(0.15)                                        # be polite
    return out


def final_answer(path):
    text = None
    with open(path) as fh:
        for line in fh:
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            if ev.get("type") == "result":
                text = ev.get("result")
    return text or ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", default="")
    ap.add_argument("--arm", default="")
    args = ap.parse_args()

    runs = ROOT / "results" / "runs.csv"
    rows = list(csv.DictReader(runs.open())) if runs.exists() else []
    if args.task:
        rows = [r for r in rows if r["task_id"] in set(args.task.split(","))]
    if args.arm:
        rows = [r for r in rows if r["arm"] in set(args.arm.split(","))]

    cache = load_cache()
    out, absent = [], []
    for r in rows:
        tp = ROOT / (r.get("transcript_path") or
                     f"transcripts/{r['task_id']}__{r['arm']}__r{r['rep']}.jsonl")
        if not tp.exists():
            continue
        text = final_answer(tp)
        found = ([("doi", d.rstrip('.')) for d in DOI_RE.findall(text)] +
                 [("pmid", p) for p in PMID_RE.findall(text)] +
                 [("pdb", p) for p in PDB_RE.findall(text)] +
                 [("chembl", c.upper()) for c in CHEMBL_RE.findall(text)] +
                 [("cid", c) for c in CID_RE.findall(text)] +
                 [("uniprot", u.upper()) for u in UNIPROT_RE.findall(text)])
        for kind, ident in sorted(set(found)):
            verdict, title = resolve(kind, ident, cache)
            out.append({"task_id": r["task_id"], "arm": r["arm"], "rep": r["rep"],
                        "kind": kind, "identifier": ident,
                        "verdict": verdict, "title": title,
                        "judge_hallucination": r.get("hallucination", "")})
            if verdict == "absent":
                absent.append((r["task_id"], r["arm"], kind, ident))
            print(f"{r['task_id']:<10}{r['arm']:<15}{kind:<5}{ident:<34}"
                  f"{verdict:<12}{title[:40]}")

    CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["task_id", "arm", "rep", "kind",
                                           "identifier", "verdict", "title",
                                           "judge_hallucination"])
        w.writeheader()
        w.writerows(out)

    print(f"\n{len(out)} identifier(s) checked across {len(rows)} run(s)"
          f" -> {OUT.relative_to(ROOT)}")
    if absent:
        print(f"\n{len(absent)} UNRESOLVABLE — these are real fabrications:")
        for t, a, k, i in absent:
            print(f"  {t}/{a}: {k} {i}")
    elif not out:
        # Saying "everything resolves" after checking nothing is how the first
        # version of this script reported eight molecule runs as clean while
        # matching none of their identifiers. An empty result is not a pass.
        print("\nNO identifiers matched in these runs. That is NOT a clean "
              "result — either the answers cite nothing checkable, or the "
              "patterns in this script do not cover the identifier types they "
              "use. Inspect a transcript before reading this as a pass.")
    else:
        print(f"\nAll {len(out)} identifier(s) cited in these runs resolve to "
              f"real records.")
    errs = [o for o in out if o["verdict"].startswith("error")]
    if errs:
        print(f"\n{len(errs)} lookup(s) errored (network/rate limit), not a "
              f"verdict — re-run to retry:")
        for o in errs[:10]:
            print(f"  {o['task_id']}/{o['arm']}: {o['kind']} {o['identifier']} "
                  f"{o['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
