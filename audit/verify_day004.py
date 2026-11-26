#!/usr/bin/env python3
"""Reference verification for day-004 (CRISPR guide design).

The blinded judge scored both smarts.bio and Biomni 0 with a hallucination
flag on this task, on the grounds that they "fabricated guide sequences,
coordinates and off-target scores". That is a strong claim about a published
result, so every checkable assertion in both answers is re-derived here from
the GRCh38 reference via NCBI E-utilities.

This is not a human audit and does not replace one. It is a deterministic
check of the claims that *can* be checked mechanically: does the protospacer
exist, at the stated coordinates, on the stated strand, with the stated PAM,
and do the named off-target loci carry the stated mismatches. Scores that are
model outputs by construction (CFD, MIT specificity, Doench) are NOT checked
and are reported as unverified.

    python3 -I audit/verify_day004.py

Needs network. Fetches ~1 kb from NCBI; no API key required.
"""

import sys
import urllib.request

EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={acc}&rettype=fasta&retmode=text"
          "&seq_start={start}&seq_stop={stop}")

CHR7, CHR1, CHRX = "NC_000007.14", "NC_000001.11", "NC_000023.11"


def fetch(acc, start, stop):
    """One region of the reference, uppercased, with its 1-based offset."""
    with urllib.request.urlopen(EFETCH.format(acc=acc, start=start, stop=stop),
                                timeout=60) as fh:
        body = fh.read().decode()
    seq = "".join(l.strip() for l in body.splitlines() if not l.startswith(">"))
    return seq.upper(), start


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def mismatches(a, b):
    return [i + 1 for i, (x, y) in enumerate(zip(a, b)) if x != y]


def locate(protospacer, regions):
    """Find a protospacer on either strand. Returns (strand, start, end, seq, off)."""
    for seq, off in regions:
        for probe, strand in ((protospacer, "+"), (rc(protospacer), "-")):
            i = seq.find(probe)
            if i != -1:
                return strand, off + i, off + i + len(probe) - 1, seq, off
    return None


def pam_of(strand, start, end, seq, off):
    """The three bases 3' of the protospacer, read on the protospacer's strand."""
    if strand == "+":
        return seq[end - off + 1:end - off + 4]
    return rc(seq[start - off - 3:start - off])


# (arm, id, protospacer, claimed PAM, claimed strand, claimed interval)
# smarts.bio reports the protospacer interval; Biomni reports protospacer+PAM.
GUIDES = [
    ("smartsbio", "g1", "AATGGTGCCAGGCATAATCC", "AGG", "-", 117559558, 117559577, "spacer"),
    ("smartsbio", "g2", "CAAAGCATGCCAACTAGAAG", "AGG", "+", 117559634, 117559653, "spacer"),
    ("smartsbio", "g3", "GGGAGAACTGGAGCCTTCAG", "AGG", "+", 117559487, 117559506, "spacer"),
    ("smartsbio", "g4", "GAGGGTAAAATTAAGCACAG", "TGG", "+", 117559506, 117559525, "spacer"),
    ("smartsbio", "g5", "ATTAAAGAAAATATCATCTT", "TGG", "+", 117559575, 117559594, "spacer"),
    ("smartsbio", "g6", "TCTGTATCTATATTCATCAT", "AGG", "-", 117559606, 117559625, "spacer"),
    ("smartsbio", "g7", "AGTTTCTTACCTCTTCTAGT", "TGG", "-", 117559646, 117559665, "spacer"),
    ("biomni", "g1", "CCAGAAAAAACATCGCCGAA", "GGG", "-", 117509084, 117509106, "spacer+pam"),
    ("biomni", "g2", "CTATGACCCGGATAACAAGG", "AGG", "+", 117530949, 117530971, "spacer+pam"),
    ("biomni", "g3", "ATCGCGATTTATCTAGGCAT", "AGG", "+", 117530980, 117531002, "spacer+pam"),
    ("biomni", "g4", "TTCCTATGACCCGGATAACA", "AGG", "+", 117530946, 117530968, "spacer+pam"),
    ("biomni", "g5", "TCCAGAAAAAACATCGCCGA", "AGG", "-", 117509085, 117509107, "spacer+pam"),
    ("biomni", "g6", "AATGGTGCCAGGCATAATCC", "AGG", "-", 117559555, 117559577, "spacer+pam"),
]


def main():
    print("Fetching GRCh38 reference regions from NCBI ...")
    regions = [fetch(CHR7, 117559450, 117559700),
               fetch(CHR7, 117508950, 117509250),
               fetch(CHR7, 117530820, 117531130)]

    print(f"\n{'arm':<11}{'id':<4}{'exists':<8}{'coords':<8}{'strand':<8}{'PAM':<8}")
    print("-" * 47)
    failures = []
    for arm, gid, ps, pam, strand, lo, hi, kind in GUIDES:
        hit = locate(ps, regions)
        if not hit:
            print(f"{arm:<11}{gid:<4}{'NO':<8}{'-':<8}{'-':<8}{'-':<8}")
            failures.append(f"{arm}/{gid}: protospacer absent from reference")
            continue
        st, s, e, seq, off = hit
        actual_pam = pam_of(st, s, e, seq, off)
        iv = (s, e + 3) if st == "+" else (s - 3, e)
        claimed_iv = (lo, hi) if kind == "spacer+pam" else (lo, hi)
        got_iv = iv if kind == "spacer+pam" else (s, e)
        ok_c, ok_s, ok_p = got_iv == claimed_iv, st == strand, actual_pam == pam
        for ok, what in ((ok_c, "coordinates"), (ok_s, "strand"), (ok_p, "PAM")):
            if not ok:
                failures.append(f"{arm}/{gid}: {what} wrong "
                                f"(claimed {claimed_iv if what=='coordinates' else (strand, pam)}, "
                                f"actual {got_iv if what=='coordinates' else (st, actual_pam)})")
        y = lambda b: "OK" if b else "WRONG"
        print(f"{arm:<11}{gid:<4}{'YES':<8}{y(ok_c):<8}{y(ok_s):<8}{y(ok_p):<8}")

    # The F508 codon, where the two arms disagree with each other.
    seq, off = regions[0]
    at = lambda a, b: seq[a - off:b - off + 1]
    print("\nF508 codon (the two arms give different coordinates):")
    print(f"  smarts.bio  chr7:117,559,593-595 = {at(117559593, 117559595)}"
          f"   {'correct (TTT = Phe)' if at(117559593,117559595)=='TTT' else 'WRONG'}")
    print(f"  biomni      chr7:117,559,594-596 = {at(117559594, 117559596)}"
          f"   {'correct (TTT = Phe)' if at(117559594,117559596)=='TTT' else 'WRONG — off by one'}")
    if at(117559594, 117559596) != "TTT":
        failures.append("biomni: F508 codon coordinate off by one, in a claim it "
                        "stated it had confirmed against the reference")

    # One off-target claim from each arm, the least verifiable part of both answers.
    print("\nOff-target spot checks:")
    xseq, xoff = fetch(CHRX, 15042920, 15042960)
    site = xseq[15042930 - xoff:15042950 - xoff]
    g2 = "CAAAGCATGCCAACTAGAAG"
    opam = xseq[15042950 - xoff:15042953 - xoff]
    mm = mismatches(g2, site)
    print(f"  smarts.bio g2 -> chrX:15,042,930-949: {len(mm)} mismatch at {mm} "
          f"(claimed 1 at pos 11), PAM {opam} (claimed GCA), "
          f"NGG/NAG: {opam[1:] in ('GG','AG')} -> not cleavable as claimed")
    if mm != [11] or opam != "GCA":
        failures.append("smarts.bio: g2 off-target claim does not match reference")

    c1seq, c1off = fetch(CHR1, 42457238, 42457278)
    g1 = "CCAGAAAAAACATCGCCGAA"
    best = None
    for start in range(c1off, c1off + len(c1seq) - 20):
        window = c1seq[start - c1off:start - c1off + 20]
        if len(window) < 20:
            continue
        for cand, st in ((window, "+"), (rc(window), "-")):
            d = mismatches(g1, cand)
            if best is None or len(d) < len(best[0]):
                best = (d, st, start, cand)
    print(f"  biomni g1     -> chr1:42,457,248: {len(best[0])} mismatches at {best[0]} "
          f"strand {best[1]} (claimed 3 at pos 5, 6, 17)")
    if best[0] != [5, 6, 17]:
        failures.append("biomni: g1 off-target claim does not match reference")

    print("\nNOT verified here, because they are model predictions with no "
          "reference answer: CFD, MIT specificity, Doench/Moreno-Mateos/Bae "
          "scores, and the gene names assigned to off-target loci.")

    print("\n" + "=" * 60)
    if failures:
        print(f"{len(failures)} discrepancy(ies):")
        for f in failures:
            print(f"  - {f}")
    else:
        print("Every mechanically checkable claim in both answers is correct.")
    print("=" * 60)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
