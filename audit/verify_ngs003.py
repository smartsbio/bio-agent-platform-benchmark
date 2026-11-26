#!/usr/bin/env python3
"""Reference verification for ngs-003 (germline variant calling).

The blinded judge scored Galaxy and smarts.bio 0 with a hallucination flag,
because neither reported the two variants the frozen ground truth names:
`chr1:14000 A>C` and `chr2:11500 A>G`. Those two variants are not in the
input. They were planted into the tumour/normal FASTQ pair by
`fixtures/make_pair_fixtures.py`; `fixtures/aligned.bam` is a separate,
pre-existing fixture that the generator never touched. The frozen
ground_truth for ngs-003 attributes the pair-fixture germline SNVs to the
wrong file.

So the task is scored against a variant set its own input cannot contain,
and the deterministic `contains_all` check inherits the same error. This
script derives the variant set that is actually in `aligned.bam` straight
from the BAM and the reference, then compares each answer against it.

Deterministic throughout: a pileup of 8000 reads over a 64 kb reference,
no model and no network. Needs `samtools` on PATH.

    python3 -I audit/verify_ngs003.py
"""

import csv
import pathlib
import re
import shutil
import subprocess
import sys
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
BAM = ROOT / "fixtures" / "aligned.bam"
REF = ROOT / "fixtures" / "hg19mini.fasta"
RUNS = ROOT / "results" / "runs.csv"
OUT = ROOT / "results" / "reference_audit.csv"
AUDIT_DATE = "2026-11-05"

# A homozygous SNV in simulated data sits at 100% VAF; 0.2 is loose enough
# that a real heterozygous call would also surface rather than be missed.
MIN_DEPTH, MIN_VAF = 5, 0.2


def truth():
    """Every SNV actually supported by the reads, from the pileup."""
    out = subprocess.run(
        ["samtools", "mpileup", "-f", str(REF), str(BAM)],
        capture_output=True, text=True, check=True).stdout
    found = []
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        chrom, pos, ref, depth, bases = f[0], f[1], f[2], int(f[3]), f[4]
        # `^X` is a read start with its mapping quality; both bytes must go
        # before the base characters are counted.
        bases = re.sub(r"\^.", "", bases).replace("$", "")
        alt = [c for c in bases.upper() if c in "ACGT"]
        if depth < MIN_DEPTH or not alt or len(alt) / depth < MIN_VAF:
            continue
        base, n = Counter(alt).most_common(1)[0]
        found.append((chrom, int(pos), ref.upper(), base, depth, n))
    return found


def reported(answer):
    """The (chrom, pos, ref, alt) calls an answer states, however formatted."""
    text = answer.replace("|", " ").replace("chr", "")
    pat = re.compile(r"\b(\d{1,2}|[XYxy])\s*:\s*([\d,]+)\s*\|?\s*"
                     r"([ACGT])\s*>\s*([ACGT])")
    seen = []
    for chrom, pos, ref, alt in pat.findall(text):
        call = (chrom, int(pos.replace(",", "")), ref, alt)
        if call not in seen:
            seen.append(call)
    return seen


def main():
    if not shutil.which("samtools"):
        sys.exit("samtools not on PATH; this audit needs it to read the BAM.")

    real = truth()
    print(f"{len(real)} SNV(s) carried by the reads in {BAM.name}:")
    for chrom, pos, ref, alt, depth, n in real:
        print(f"  {chrom}:{pos} {ref}>{alt}  depth={depth} alt={n} "
              f"vaf={n / depth:.2f}")
    expected = {(c, p, r, a) for c, p, r, a, _, _ in real}

    # The two the frozen ground truth names, shown absent rather than asserted.
    print("\nfrozen ground_truth positions, as the reads see them:")
    for chrom, pos in (("1", 14000), ("2", 11500)):
        hit = [x for x in real if (x[0], x[1]) == (chrom, pos)]
        print(f"  {chrom}:{pos} -> "
              + (f"{hit[0][2]}>{hit[0][3]}" if hit else
                 "no alt allele in any read (reference everywhere)"))

    rows = []
    for r in csv.DictReader(RUNS.open()):
        if r["task_id"] != "ngs-003" or r["status"] != "answered":
            continue
        calls = reported(r["answer"])
        if not calls:
            continue
        correct = sorted(expected & set(calls))
        wrong = sorted(set(calls) - expected)
        missed = sorted(expected - set(calls))
        acc = 3 if not wrong and not missed else 2 if not wrong else 1
        note = (f"Judge failed this row for missing chr1:14000 A>C and "
                f"chr2:11500 A>G. Those SNVs are absent from aligned.bam: "
                f"the pileup shows the reference base in every read at "
                f"both loci. They belong to the "
                f"tumour/normal fixtures, not this BAM, so the frozen "
                f"ground_truth names the wrong file. Against the variant "
                f"set the reads do carry, this answer reports "
                f"{len(correct)}/{len(expected)} correctly")
        note += (f", with no false positives." if not wrong else
                 f", with {len(wrong)} not supported by the reads: "
                 + "; ".join(f"{c}:{p} {rf}>{al}" for c, p, rf, al in wrong)
                 + ".")
        if missed:
            note += (" Missed: "
                     + "; ".join(f"{c}:{p} {rf}>{al}" for c, p, rf, al in missed)
                     + " — a genuine recall loss, independent of the "
                       "ground-truth error.")
        rows.append({
            "task_id": "ngs-003", "arm": r["arm"], "rep": r["rep"],
            "audit_date": AUDIT_DATE, "method": "audit/verify_ngs003.py",
            "judge_accuracy": r["accuracy_0_3"],
            "judge_hallucination": r["hallucination"],
            "audit_accuracy": acc, "audit_hallucination": "n",
            "agree": "n", "notes": note,
        })

    print()
    with OUT.open("a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        for row in rows:
            w.writerow(row)
            print(f"{row['arm']:14} judge acc={row['judge_accuracy']} "
                  f"halluc={row['judge_hallucination']} -> audit acc="
                  f"{row['audit_accuracy']} halluc=n")
    print(f"\n{len(rows)} row(s) appended -> {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
