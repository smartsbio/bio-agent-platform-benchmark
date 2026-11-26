#!/usr/bin/env python3
"""Re-derive the ngs-004/005/006 truth sets from the fixtures themselves.

`ngs-003` was scored against germline variants that were not in its input: the
frozen ground_truth attributed the tumour/normal fixtures' SNVs to
`aligned.bam`, which never carried them. See the Known limitations entry in
the README. The remaining NGS rows draw their ground truth from those same
pair fixtures, so each one is checked here BEFORE it runs rather than audited
after.

The question is not "what did the generator intend to plant" -- that is in
`pair_fixtures_truth.json` -- but "what can a caller actually recover from the
reads as delivered". A variant inside an N block, or at a depth no caller would
call, is not a testable expectation however it was planted.

METHOD. `make_pair_fixtures.py` applies the variants to the genome and then
cuts error-free 100 bp reads from it (every quality is `I`), so a read covering
a variant contains the alt allele exactly. No aligner is needed and none is
installed: a 41 bp window centred on each locus is taken from the reference,
the alt version is built, and both are counted as exact substrings of the
reads in each FASTQ, on both strands. That is a depth and VAF measurement with
no alignment step to go wrong.

    python3 -I audit/derive_ngs004_truth.py

Reads only the fixtures. No model, no network.
"""

import gzip
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
FIX = ROOT / "fixtures"
TRUTH = FIX / "pair_fixtures_truth.json"
REF = FIX / "hg19mini.fasta"
PEAKS = FIX / "peaks_hg19mini.bed"
GTF = FIX / "annotation.gtf"

FLANK = 20          # 41 bp window: long enough to be unique in a 64 kb genome
MIN_CALLABLE = 10   # below this depth no caller is expected to make a call


def read_fasta(path):
    seqs, name = {}, None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def reads(path):
    with gzip.open(path, "rt") as fh:
        return [l.strip().upper() for i, l in enumerate(fh) if i % 4 == 1]


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def count(windows, pool):
    """How many reads contain each window, either strand."""
    hits = {k: 0 for k in windows}
    for r in pool:
        for k, w in windows.items():
            if w in r or rc(w) in r:
                hits[k] += 1
    return hits


def observe(genome, pool, variants):
    """Measured depth and alt fraction at each locus, from the reads alone."""
    out = []
    for v in variants:
        chrom, pos, ref, alt = v["chrom"], v["pos"], v["ref"], v["alt"]
        seq = genome[chrom]
        left, right = seq[pos - 1 - FLANK:pos - 1], seq[pos:pos + FLANK]
        actual = seq[pos - 1]
        win = {"ref": left + ref + right, "alt": left + alt + right}
        hits = count(win, pool)
        depth = hits["ref"] + hits["alt"]
        out.append({
            "chrom": chrom, "pos": pos, "ref": ref, "alt": alt,
            "ref_base_matches_reference": actual == ref,
            "reference_base_here": actual,
            "n_flank": seq[pos - 1 - FLANK:pos + FLANK].count("N"),
            "ref_reads": hits["ref"], "alt_reads": hits["alt"],
            "depth": depth,
            "vaf": round(hits["alt"] / depth, 3) if depth else None,
        })
    return out


def show(label, rows):
    print(f"\n{label}")
    for r in rows:
        flag = []
        if not r["ref_base_matches_reference"]:
            flag.append(f"REF MISMATCH (reference has "
                        f"{r['reference_base_here']})")
        if r["n_flank"]:
            flag.append(f"{r['n_flank']} N in flank")
        if r["depth"] < MIN_CALLABLE:
            flag.append("NOT CALLABLE: depth below "
                        f"{MIN_CALLABLE}")
        print(f"  {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}  "
              f"depth={r['depth']:<4} alt={r['alt_reads']:<4} "
              f"vaf={r['vaf']}  {'; '.join(flag) or 'ok'}")


def main():
    truth = json.loads(TRUTH.read_text())
    genome = read_fasta(REF)
    print(f"reference: {', '.join(f'{k}={len(v)}bp' for k, v in genome.items())}")

    tumour = reads(FIX / "tumor_R1.fastq.gz") + reads(FIX / "tumor_R2.fastq.gz")
    normal = reads(FIX / "normal_R1.fastq.gz") + reads(FIX / "normal_R2.fastq.gz")
    print(f"tumour reads: {len(tumour)}   normal reads: {len(normal)}")

    som = truth["somatic_variants"]
    germ = truth["germline_variants"]

    print("\n" + "=" * 68)
    print("ngs-004 — somatic calling (tumour/normal)")
    print("=" * 68)
    show("SOMATIC, in TUMOUR (must be called):", observe(genome, tumour, som))
    show("SOMATIC, in NORMAL (must be absent):", observe(genome, normal, som))
    show("GERMLINE, in TUMOUR (present, must NOT be called somatic):",
         observe(genome, tumour, germ))
    show("GERMLINE, in NORMAL (present — this is what subtracts them):",
         observe(genome, normal, germ))

    print("\n" + "=" * 68)
    print("ngs-005 — RNA-seq quantification")
    print("=" * 68)
    feats = [l.split("\t") for l in GTF.read_text().splitlines()
             if l.strip() and not l.startswith("#")]
    genes = sorted({re.search(r'gene_id "([^"]+)"', f[8]).group(1)
                    for f in feats if len(f) > 8
                    and re.search(r'gene_id "([^"]+)"', f[8])})
    exons = [f for f in feats if f[2] == "exon"]
    print(f"  annotation.gtf: {len(feats)} features, {len(exons)} exons, "
          f"{len(genes)} gene(s): {', '.join(genes)}")
    print("  reads for this task are the uniform ngs-001/002 pair, simulated "
          "with\n  no transcript model, so per-gene counts are whatever falls "
          "in each exon\n  by chance. There is no planted expected count: the "
          "frozen ground_truth\n  asks only that a real count matrix is "
          "produced. That is checkable; an\n  exact count is not.")

    print("\n" + "=" * 68)
    print("ngs-006 — ChIP peak calling")
    print("=" * 68)
    bed = [l.split("\t") for l in PEAKS.read_text().splitlines() if l.strip()]
    print(f"  peaks_hg19mini.bed: {len(bed)} planted region(s)")
    chip = reads(FIX / "chip_R1.fastq.gz") + reads(FIX / "chip_R2.fastq.gz")
    ctrl = reads(FIX / "input_R1.fastq.gz") + reads(FIX / "input_R2.fastq.gz")
    print(f"  chip reads: {len(chip)}   input reads: {len(ctrl)}")
    for chrom, start, end, *rest in bed:
        start, end = int(start), int(end)
        seq = genome[chrom]
        if seq[start:end].count("N"):
            print(f"  {chrom}:{start}-{end}  N BLOCK — not callable")
            continue
        # Enrichment is measured by read starts inside the region, which for
        # error-free reads is an exact substring count of the first 30 bp.
        probes = {f"p{i}": seq[p:p + 30]
                  for i, p in enumerate(range(start, end - 30, 50))}
        c = sum(count(probes, chip).values())
        n = sum(count(probes, ctrl).values())
        ratio = round(c / n, 2) if n else "inf"
        print(f"  {chrom}:{start}-{end}  chip={c:<5} input={n:<5} "
              f"fold={ratio}")

    print("\nnothing was written; this is a read-only derivation.")


if __name__ == "__main__":
    main()
