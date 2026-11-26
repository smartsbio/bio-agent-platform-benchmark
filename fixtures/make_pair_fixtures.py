#!/usr/bin/env python3
"""Generate the paired fixtures that somatic-variant-calling and chip-seq need.

Both pipelines compare two samples. Until now the suite passed the SAME reads as
both sides, so the only correct answer was "nothing found" -- the expectations
literally said "~0 somatic calls (tumor = normal)" and "few/no peaks is fine".
A pipeline that silently produced nothing therefore passed, which is exactly the
failure mode that let RoseTTAFold3 discard every prediction for weeks.

These fixtures have known ground truth instead:

  somatic   3 somatic SNVs present only in the tumour, plus 2 germline SNVs
            present in BOTH. A caller that reports the germline ones as somatic
            is wrong in a way "tumour = normal" could never reveal.

  chip      ChIP reads enriched ~10x inside the 6 regions of peaks_hg19mini.bed;
            input reads uniform. Peaks must be called there and only there.

Everything is seeded, so re-running reproduces the files byte for byte. Reads are
error-free with uniform high quality: the point is to test the pipelines' logic,
not a basecaller.

Usage:  python3 -I make_pair_fixtures.py          (writes into this directory)
"""

import gzip
import json
import pathlib
import random

HERE = pathlib.Path(__file__).resolve().parent
REFERENCE = HERE / "hg19mini.fasta"
PEAKS = HERE / "peaks_hg19mini.bed"

READ_LEN = 100
FRAGMENT = 300
SEED = 20261101

# hg19mini is 4 x 16 kb but only 10001-16000 of chr1, chr2 and chr4 is real
# sequence -- the rest is N, and chr3 is N end to end. Variants and reads must
# stay inside the real blocks: nothing aligns to an N run, so a variant planted
# there is invisible to any caller. (The same caveat applies to the SNPs the
# skill documents for aligned.bam, most of which are in N regions.)
REAL_BLOCK = (10001, 16000)
REAL_CHROMS = ("1", "2", "4")

# Somatic SNVs: tumour only. Spread across the usable contigs so a per-contig
# bug is visible, and clear of the ChIP peak regions.
# Alt bases are transversions away from the reference base at each locus
# (1:11000 C, 2:13000 A, 4:12000 G, 1:14000 A, 2:11500 A).
SOMATIC = [
    ("1", 11000, "T"),
    ("2", 13000, "C"),
    ("4", 12000, "T"),
]
# Germline SNVs: in BOTH samples. A somatic caller must subtract these.
GERMLINE = [
    ("1", 14000, "C"),
    ("2", 11500, "G"),
]

COVERAGE = 30          # per sample, per mate
PEAK_ENRICHMENT = 0.60  # fraction of ChIP fragments drawn from peak regions


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def apply_variants(seqs, variants):
    """Return a copy of the genome with `variants` (1-based) substituted in."""
    out = {k: list(v) for k, v in seqs.items()}
    for chrom, pos, alt in variants:
        ref = out[chrom][pos - 1]
        if ref.upper() == "N":
            raise SystemExit(
                f"variant {chrom}:{pos} sits in an N run; no read can support it"
            )
        if ref == alt:
            raise SystemExit(
                f"variant {chrom}:{pos} {ref}>{alt} is not a change; pick another base"
            )
        out[chrom][pos - 1] = alt
    return {k: "".join(v) for k, v in out.items()}


def revcomp(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def load_peaks():
    regions = []
    for line in PEAKS.read_text().splitlines():
        if not line.strip():
            continue
        chrom, start, end = line.split("\t")[:3]
        regions.append((chrom, int(start), int(end)))
    return regions


def write_pair(prefix, fragments, genome):
    """Write <prefix>_R1.fastq.gz / _R2.fastq.gz for a list of (chrom, start)."""
    qual = "I" * READ_LEN
    r1_path = HERE / f"{prefix}_R1.fastq.gz"
    r2_path = HERE / f"{prefix}_R2.fastq.gz"
    # mtime=0 so a re-run is byte-identical rather than differing by timestamp.
    with gzip.GzipFile(r1_path, "wb", mtime=0) as r1, \
         gzip.GzipFile(r2_path, "wb", mtime=0) as r2:
        for i, (chrom, start) in enumerate(fragments, 1):
            frag = genome[chrom][start:start + FRAGMENT]
            if len(frag) < FRAGMENT:
                continue
            name = f"{prefix}_read{i:06d}"
            fwd = frag[:READ_LEN]
            rev = revcomp(frag[-READ_LEN:])
            r1.write(f"@{name}/1\n{fwd}\n+\n{qual}\n".encode())
            r2.write(f"@{name}/2\n{rev}\n+\n{qual}\n".encode())
    return r1_path, r2_path


def uniform_fragments(rng, genome, n):
    """Fragments drawn only from the real (non-N) block of the usable contigs."""
    lo, hi = REAL_BLOCK
    out = []
    for _ in range(n):
        chrom = rng.choice(REAL_CHROMS)
        # 0-based start, keeping the whole fragment inside the real block.
        out.append((chrom, rng.randrange(lo - 1, hi - FRAGMENT)))
    return out


def enriched_fragments(rng, genome, n, regions):
    lo, hi = REAL_BLOCK
    out = []
    n_peak = int(n * PEAK_ENRICHMENT)
    for _ in range(n_peak):
        chrom, s, e = rng.choice(regions)
        # Centre the fragment on the peak so the pileup has a clear summit,
        # clamped to the real block so no read comes from an N run.
        centre = rng.randrange(s, e)
        start = max(lo - 1, min(centre - FRAGMENT // 2, hi - FRAGMENT))
        out.append((chrom, start))
    out += uniform_fragments(rng, genome, n - n_peak)
    rng.shuffle(out)
    return out


def main():
    genome = read_fasta(REFERENCE)
    lo, hi = REAL_BLOCK
    # Coverage is over the sequence reads can actually align to, not the 64 kb
    # of mostly-N contigs, or the real depth would be under half what it says.
    total_bp = (hi - lo + 1) * len(REAL_CHROMS)
    n_frags = total_bp * COVERAGE // READ_LEN
    regions = load_peaks()
    off_target = [r for r in regions if r[0] not in REAL_CHROMS]
    if off_target:
        raise SystemExit(f"peak regions on unusable contigs: {off_target}")

    normal_genome = apply_variants(genome, GERMLINE)
    tumour_genome = apply_variants(genome, GERMLINE + SOMATIC)

    rng = random.Random(SEED)
    written = []

    # --- somatic: tumour vs normal, same depth, independent draws
    written += write_pair("normal", uniform_fragments(rng, normal_genome, n_frags), normal_genome)
    written += write_pair("tumor", uniform_fragments(rng, tumour_genome, n_frags), tumour_genome)

    # --- chip: enriched vs uniform input
    written += write_pair("chip", enriched_fragments(rng, genome, n_frags, regions), genome)
    written += write_pair("input", uniform_fragments(rng, genome, n_frags), genome)

    truth = {
        "reference": REFERENCE.name,
        "seed": SEED,
        "coverage_per_sample": COVERAGE,
        "read_length": READ_LEN,
        "fragment_length": FRAGMENT,
        "somatic_variants": [
            {"chrom": c, "pos": p, "ref": genome[c][p - 1], "alt": a,
             "expected": "called as somatic"}
            for c, p, a in SOMATIC
        ],
        "germline_variants": [
            {"chrom": c, "pos": p, "ref": genome[c][p - 1], "alt": a,
             "expected": "present in both samples; must NOT be called somatic"}
            for c, p, a in GERMLINE
        ],
        "chip_peak_regions": [
            {"chrom": c, "start": s, "end": e} for c, s, e in regions
        ],
        "chip_enrichment_fraction": PEAK_ENRICHMENT,
    }
    truth_path = HERE / "pair_fixtures_truth.json"
    truth_path.write_text(json.dumps(truth, indent=2) + "\n")

    for p in written + [truth_path]:
        print(f"  {p.name:28} {p.stat().st_size:>9,} bytes")
    print(f"\n{n_frags:,} fragments per sample "
          f"(~{COVERAGE}x over {total_bp:,} bp across {len(genome)} contigs)")


if __name__ == "__main__":
    main()
