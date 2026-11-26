#!/usr/bin/env python3
"""Derive the exact per-gene count matrix for ngs-005 from the fixtures.

The frozen ground_truth for ngs-005 asks only that "a real count matrix per
gene is produced", which makes it a rubric row wearing a numeric row's clothes:
an arm that invents six plausible integers scores the same as one that runs
featureCounts. The counts are in fact exactly determined by the fixtures, so
this script computes them and the row can be checked on the numbers.

No aligner is installed and none is needed. `reads_R{1,2}.fastq.gz` are
error-free 100 bp cuts (every quality is `I`), so each read occurs verbatim at
one position on one strand. Locating it by seed-and-extend IS the alignment for
this data, and it is exact rather than heuristic.

One wrinkle, worth recording because it confirms a separate finding: 270 of the
8000 reads are NOT substrings of `hg19mini.fasta`. They are substrings of the
reference with the six SNVs of `aligned.bam` applied, so these reads and that
BAM are a matched pair and those six SNVs are the fixture's planted variant
set -- independent corroboration of the ngs-003 correction. Placement therefore
allows one mismatch, which covers a variant-carrying read without needing the
variant set to be known in advance.

COUNTING RULE, which the check must state because tools differ on it:
  - a mate is assigned to a gene when it overlaps >= 1 bp of any exon of that
    gene (featureCounts default `-t exon -g gene_id`, min overlap 1);
  - the six genes in `annotation.gtf` do not overlap each other, so no mate can
    be ambiguous between genes -- verified below rather than assumed;
  - fragment counts (`-p --countReadPairs`) assign a pair once if either mate
    lands in the gene; read counts assign each mate separately.
Both are reported, and the check accepts either. This script originally
computed only the read-level column and `checks.json` was frozen on it; all
three tool-running arms reported the fragment column, so three correct count
matrices were failed. The difference is a constant ~396 per gene -- the
4 x (readlen - 1) positions where a mate straddles an exon edge -- which is
what a constant offset across every gene should have suggested. A tool run
with intron-aware or union-intersection semantics may differ again, so the
check accepts either column rather than one exact integer.

    python3 -I audit/derive_ngs005_truth.py

Reads only the fixtures. No model, no network, no platform.
"""

import collections
import gzip
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIX = ROOT / "fixtures"
REF = FIX / "hg19mini.fasta"
GTF = FIX / "annotation.gtf"
R1, R2 = FIX / "reads_R1.fastq.gz", FIX / "reads_R2.fastq.gz"

SEED = 25   # unique in a 64 kb genome; verified by the full-length recheck
MAX_MM = 1  # a read crossing a planted SNV differs from the reference by one


def mismatches(a, b):
    return sum(x != y for x, y in zip(a, b))


def read_fasta(path):
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def fastq(path):
    """(name, sequence) per record, mate suffix stripped from the name."""
    out = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].strip().split()[0].rsplit("/", 1)[0]
            elif i % 4 == 1:
                out.append((name, line.strip().upper()))
    return out


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def exons():
    """{gene: [(chrom, start, end), ...]} 1-based inclusive, as the GTF gives."""
    by_gene = collections.defaultdict(list)
    for line in GTF.read_text().splitlines():
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        gene = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        by_gene[gene].append((f[0], int(f[3]), int(f[4])))
    return dict(by_gene)


def check_disjoint(by_gene):
    """No two genes may share a base, or 'ambiguous' would need defining."""
    spans = []
    for gene, ex in by_gene.items():
        for chrom, s, e in ex:
            spans.append((chrom, s, e, gene))
    clashes = []
    for i, (c1, s1, e1, g1) in enumerate(spans):
        for c2, s2, e2, g2 in spans[i + 1:]:
            if c1 == c2 and g1 != g2 and s1 <= e2 and s2 <= e1:
                clashes.append((g1, g2, c1, max(s1, s2), min(e1, e2)))
    return clashes


def build_index(genome):
    """seed -> [(chrom, pos)], pos 0-based. Seeds seen twice are dropped."""
    idx, dup = {}, set()
    for chrom, seq in genome.items():
        for i in range(len(seq) - SEED + 1):
            k = seq[i:i + SEED]
            if "N" in k:
                continue
            if k in idx:
                dup.add(k)
            else:
                idx[k] = (chrom, i)
    for k in dup:
        idx.pop(k, None)
    return idx, len(dup)


def candidates(seq, genome, idx):
    """Every position this read could have come from, best-scoring first.

    Returns [(mismatches, chrom, start, end)] 1-based inclusive. Both strands
    and every seed offset are tried, and nothing is discarded, because the
    caller needs to know whether a read places in one place or several: in
    hg19mini's telomeric repeat tracts many reads place equally well in
    several, and which one an aligner picks is a property of the aligner.
    Hiding that behind a single choice is what made earlier versions of this
    script produce counts that moved when the tie-break changed.
    """
    out = {}
    for strand in (seq, rc(seq)):
        for off in range(0, len(strand) - SEED + 1, 5):
            hit = idx.get(strand[off:off + SEED])
            if not hit:
                continue
            chrom, pos = hit
            begin = pos - off
            if begin < 0:
                continue
            window = genome[chrom][begin:begin + len(strand)]
            if len(window) != len(strand):
                continue
            mm = mismatches(window, strand)
            if mm <= MAX_MM:
                key = (chrom, begin)
                out[key] = min(mm, out.get(key, mm))
    best = min(out.values(), default=None)
    return [(mm, c, b + 1, b + len(seq))
            for (c, b), mm in sorted(out.items()) if mm == best]


def assign(loc, by_gene):
    """Genes whose exons this interval overlaps by >= 1 bp."""
    chrom, start, end = loc
    return {g for g, ex in by_gene.items()
            for c, s, e in ex if c == chrom and start <= e and s <= end}


SWEEP = ((25, 5, 1), (25, 1, 1), (20, 5, 1), (30, 5, 1), (25, 5, 0))


def fragment_of(name):
    """The pair a mate belongs to, however the simulator spelled the suffix."""
    return name.split()[0].split("/")[0]


def tally(genome, by_gene, mates, seed, step, maxmm):
    """One pass of placement-and-counting at the given parameters.

    Returns read-level counts (each mate counted separately, featureCounts
    default) and fragment-level counts (`-p --countReadPairs`, a pair counted
    once if either mate lands in the gene). Which one a tool reports is a
    flag, not a result, so both are carried and the check accepts either.
    """
    global SEED, MAX_MM
    SEED, MAX_MM = seed, maxmm
    idx, _ = build_index(genome)
    counts, unplaced = collections.Counter(), 0
    frags = collections.defaultdict(list)
    for name, seq in mates:
        found = {}
        for strand in (seq, rc(seq)):
            for off in range(0, len(strand) - seed + 1, step):
                hit = idx.get(strand[off:off + seed])
                if not hit:
                    continue
                chrom, pos = hit
                begin = pos - off
                if begin < 0:
                    continue
                window = genome[chrom][begin:begin + len(strand)]
                if len(window) != len(strand):
                    continue
                mm = mismatches(window, strand)
                if mm <= maxmm:
                    found[(chrom, begin)] = min(mm, found.get((chrom, begin), mm))
        if not found:
            unplaced += 1
            continue
        best = min(found.values())
        at_best = [k for k, v in found.items() if v == best]
        if len(at_best) != 1:
            continue
        chrom, begin = at_best[0]
        loc = (chrom, begin + 1, begin + len(seq))
        genes = assign(loc, by_gene)
        if len(genes) == 1:
            counts[genes.pop()] += 1
        frags[fragment_of(name)].append(loc)

    # A fragment spans the leftmost to the rightmost base its mates cover, so
    # a pair straddling an exon edge counts once for that gene.
    fcounts = collections.Counter()
    for locs in frags.values():
        if len({l[0] for l in locs}) != 1:
            continue
        span = (locs[0][0], min(l[1] for l in locs), max(l[2] for l in locs))
        genes = assign(span, by_gene)
        if len(genes) == 1:
            fcounts[genes.pop()] += 1
    return counts, fcounts, unplaced


def main():
    genome = read_fasta(REF)
    by_gene = exons()
    clashes = check_disjoint(by_gene)
    print(f"reference: {', '.join(f'{k}={len(v)}' for k, v in genome.items())}")
    print(f"annotation: {len(by_gene)} genes, "
          f"{sum(len(v) for v in by_gene.values())} exons")
    if clashes:
        print("  OVERLAPPING GENES — an ambiguity policy is needed:")
        for g1, g2, c, a, b in clashes:
            print(f"    {g1}/{g2} at {c}:{a}-{b}")
    else:
        print("  no two genes share a base, so no mate is ambiguous between "
              "genes")

    mates = fastq(R1) + fastq(R2)
    print(f"reads: {len(mates)}")

    # A count is only ground truth if it survives the arbitrary choices in
    # this script. Each sweep setting is a defensible way to place a read; a
    # gene whose count moves between them is not exactly determined, and
    # saying otherwise would hand the benchmark a number two correct aligners
    # would disagree on.
    print(f"\nplacement sweep over {len(SWEEP)} parameter settings "
          "(seed, step, max mismatch):")
    runs, fruns = [], []
    for seed, step, maxmm in SWEEP:
        counts, fcounts, unplaced = tally(genome, by_gene, mates,
                                          seed, step, maxmm)
        runs.append(counts)
        fruns.append(fcounts)
        print(f"  seed={seed:<3} step={step:<2} mm={maxmm}  "
              f"unplaced={unplaced:<4} "
              + " ".join(f"{g}={counts[g]}/{fcounts[g]}"
                         for g in sorted(by_gene)))
    print("  (each cell is reads/fragments)")

    print("\nGROUND TRUTH for ngs-005")
    stable, band = [], []
    for gene in sorted(by_gene):
        vals = [r[gene] for r in runs]
        fvals = [r[gene] for r in fruns]
        lo, hi = min(vals), max(vals)
        flo, fhi = min(fvals), max(fvals)
        if lo == hi and flo == fhi:
            stable.append((gene, lo, flo))
            print(f"  {gene:8} reads={lo:<6} frags={flo:<6} exact — both "
                  f"conventions unchanged across every setting")
        else:
            band.append((gene, lo, hi, flo, fhi))
            print(f"  {gene:8} reads={lo}-{hi:<5} frags={flo}-{fhi:<5} "
                  f"depends on repeat-tract placement")

    print(f"\n  {len(stable)} of {len(by_gene)} genes are exactly determined "
          f"by the fixtures.")
    if band:
        names = ", ".join(r[0] for r in band)
        print(f"  {names} sit against hg19mini's telomeric repeat tracts, "
              f"where which\n  position a read is given is a property of the "
              f"aligner, not of the data.")
    print("\n  SO THE CHECK SHOULD: for each exact gene, accept EITHER "
          "convention's\n  value, because which one a tool reports is a flag "
          "and not a result;\n  accept the band for the rest; and above all "
          "require that a real matrix\n  was produced at all — six invented "
          "integers must not pass. An arm\n  reporting these counts in "
          "either column has demonstrably run the tool.")
    if stable:
        print("\n  pairs spec for checks.json (gene, reads, fragments):")
        for gene, reads, frags in stable:
            print(f'    ["{gene}", {reads}, {frags}],')
    print("\nnothing was written; this is a read-only derivation.")


if __name__ == "__main__":
    main()
